# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: YUCCA BREVIFOLIA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage 
# in Yucca brevifolia (data files from Royer et al. 2016). Genetic data sourced from 
# CSV reformatted to STRUCTURE input file; geographic coordinates are converted and
# cleaned before being utilized for correlations.

pacman::p_load(adegenet, terra, parallel, RColorBrewer, viridis, scales, vcfR, usedist)

# Read in relevant functions
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/Code/'
setwd(GeoGenCorr_wd)
source('Scripts/functions_GeoGenCoverage.R')

# ---- VARIABLES ----
# Specify number of resampling replicates
num_reps <- 5
# ---- BUFFER SIZES
# Specify geographic buffer size in meters
geo_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500))
# Specify ecological buffer size in meters
eco_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500))

# ---- PARALLELIZATION
# Set up relevant cores
num_cores <- detectCores() - 4
cl <- makeCluster(num_cores)
# Make sure libraries (adegenet, terra, etc.) are on cluster (but avoid printing output)
invisible(clusterEvalQ(cl, library('adegenet')))
invisible(clusterEvalQ(cl, library('terra')))
invisible(clusterEvalQ(cl, library('parallel')))
invisible(clusterEvalQ(cl, library('usedist')))

# %%% CONDUCT RESAMPLING %%% ----
# ---- READ IN DATA ----
# Specify filepath for QULO geographic and genetic data
YUBR_filePath <- paste0(GeoGenCorr_wd, 'Datasets/YUBR/')

# ---- GENETIC MATRIX
# Read in the genetic matrix file provided in the Dryad repository for 
# Royer et al. 2016. This is a STRUCTURE input file of 319 individuals and 9,516
# SNP loci that has been filtered for missing data (only SNPs present in at least 80%
# of individuals, only individuals missing no more than 30% of all SNPs). However,
# the original CSV needs to be cleaned before it can be processed using the 
# adegenet::read.structure function.

# Read in the original CSV
YUBR_struTable <- read.csv2(paste0(YUBR_filePath, 'Genetic/GoodLoci319Trees.csv'),
                          header = TRUE, sep = ",")
# Remove 2nd column of the genetic matrix, which stores the "frequency" values
YUBR_struTable <- YUBR_struTable[,-2]
# Write a ".stru" file based on the genetic matrix, converting the commas to spaces
# and removing the marker names. Then, delete the YUBR_struTable object (for space)
write.table(YUBR_struTable, file=paste0(YUBR_filePath, 'Genetic/GoodLoci319Trees.stru'),
            sep =" ", quote = FALSE, row.names = FALSE, col.names = FALSE)
rm(YUBR_struTable)
# Read in the reformatted STRUCTURE file 
YUBR_genind <- read.structure(file=paste0(YUBR_filePath, 'Genetic/GoodLoci319Trees.stru'), 
                              n.ind = 319, n.loc = 9516, onerowperind = FALSE, col.lab = 1, 
                              col.pop = 0, row.marknames = 0, sep = ' ', ask = FALSE)
# Delete the ".stru" file just written to the directory
file.remove(paste0(YUBR_filePath, 'Genetic/GoodLoci319Trees.stru'))
# Reprocess the names of the individuals (in order to allow for matching with
# coordinate values, later on): remove "cat_" prefix, where present
indNames(YUBR_genind) <- gsub("cat_", "", indNames(YUBR_genind))

# ---- GEOGRAPHIC/ECOLOGICAL DATA FILES
# The original coordinates file for this Yucca dataset needs to be processed such that decimal
# minutes (in two different columns) is reorganized into decimal degrees (in a single column).
# Check if the processed file (called YUBR_coordinates.csv) already exists; if not, then 
# run the necessary processing steps.
if(file.exists(paste0(YUBR_filePath, 'Geographic/YUBR_coordinates.csv'))){
  # Read in the CSV of processed coordinates. The first column contains row numbers
  YUBR_coordinates <- read.csv(
    paste0(YUBR_filePath, 'Geographic/YUBR_coordinates.csv'), header=TRUE)
} else {
  # Read in and process the wild occurrence points. This CSV has 5 columns: sample names, 
  # latitude (decimal minutes, 2 columns), and longitude (decimal minutes, 2 columns)  
  YUBR_points <- read.csv(
    paste0(YUBR_filePath, 'Geographic/YUBR_coordinates_Original.csv'), header=TRUE)
  # Convert decimal minutes to decimal degrees by dividing the 3rd and 5th columns by
  # 60, then adding them to the degrees columns (2nd and 4th)
  lats <- round(YUBR_points[,3]/60 + YUBR_points[,2], 4)
  # For longitudes, multiply by -1 to properly locate coordinates
  longs <- -1*(round(YUBR_points[,5]/60 + YUBR_points[,4], 4))
  # Reformat coordinate values based on conversions, and rename
  YUBR_points <- cbind.data.frame(YUBR_points[,1], lats, longs)
  colnames(YUBR_points) <- c('Sample', 'decimalLatitude', 'decimalLongitude')
  # Subset the coordinates data.frame to strictly the samples included in the genetic matrix
  YUBR_coordinates <- YUBR_points[which(YUBR_points$Sample %in% indNames(YUBR_genind)),]
  # Reorder the coordinate values to match the order of samples in the genind file
  YUBR_coordinates <- YUBR_coordinates[match(indNames(YUBR_genind), YUBR_coordinates$Sample),]
  # Rename the rows numerically, to match the order of samples in the genind file
  rownames(YUBR_coordinates) <- as.character(seq(1, nrow(YUBR_coordinates)))
  # Write resulting coordinates data.frame as a CSV to disk, for future runs
  write.csv(YUBR_coordinates, file=paste0(YUBR_filePath,'Geographic/YUBR_coordinates.csv'), 
            row.names = FALSE)
}
# Read in raster data, for SDM
YUBR_sdm <- terra::rast(paste0(YUBR_filePath,'Geographic/YUBR_319inds_rast_Carver.tif'))

# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <-
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.gpkg')))
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = YUBR_coordinates)
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage 
# (solely in the U.S.)
ecoregion_poly <-
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level4/us_eco_l4.shp')))
# Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
# exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
world_poly_clip_W <- wrap(world_poly_clip)
ecoregion_poly_W <- wrap(ecoregion_poly)
YUBR_sdm_W <- wrap(YUBR_sdm)

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) 
# to the cluster
clusterExport(cl, varlist = c('YUBR_genind','YUBR_coordinates','num_reps','geo_buffSize',
                              'eco_buffSize','world_poly_clip_W', 'ecoregion_poly_W','YUBR_sdm_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres', 
                              'eco.intersectBuff','eco.compareBuff','gen.getAlleleCategories', 
                              'gen.buildDistMat', 'gen.calcGenDistCov', 'eco.totalEcoregionCount',
                              'calculateCoverage','exSituResample.Par', 'geo.gen.Resample.Par'))
# Specify file path, for saving resampling array
arrayDir <- paste0(YUBR_filePath, 'resamplingData/YUBR_SMBO3_G2G2E_resampArr.Rdata')
# Run resampling (in parallel)
YUBR_demoArray_Par <- 
  geo.gen.Resample.Par(genObj = YUBR_genind, genDistFlag=TRUE, geoFlag=TRUE, coordPts = YUBR_coordinates, 
                       geoBuff = geo_buffSize,SDMrast = YUBR_sdm_W, boundary=world_poly_clip_W, 
                       ecoFlag = TRUE, ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, 
                       ecoLayer = 'US', reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
# Close cores
stopCluster(cl)

# Run resampling not in parallel (for function testing purposes)
# YUBR_demoArray_IND <-
#   geo.gen.Resample(gen_obj = YUBR_genind, geoFlag = TRUE, coordPts = YUBR_coordinates, geoBuff = geo_buffSize, 
#                    boundary = world_poly_clip, ecoFlag = TRUE, ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly,
#                    ecoLayer = "US", reps = 1)

# %%% ANALYZE DATA %%% ----
# Specify filepath for QULO geographic and genetic data, including resampling data
YUBR_filePath <- paste0(GeoGenCorr_wd, 'Datasets/YUBR/')
arrayDir <- paste0(YUBR_filePath, 'resamplingData/YUBR_1km_GE_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
YUBR_demoArray_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values
YUBR_DF <- resample.array2dataframe(YUBR_demoArray_Par)
# Calculate normalized root mean square value
YUBR_nrmse_geo <- nrmse_func(obs=YUBR_DF$Geo, pred=YUBR_DF$Total) ; YUBR_nrmse_geo
YUBR_nrmse_eco <- nrmse_func(obs=YUBR_DF$Eco, pred=YUBR_DF$Total) ; YUBR_nrmse_eco

# ---- PLOTTING ----
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(YUBR_demoArray_Par, allValues = TRUE)
# Subset matrix of all average values to just Total allelic, geographic, and ecological coverage
averageValueMat_TEG <- averageValueMat[,c(1,6,7)]
# Calculate the absolute difference between genetic and geographic/ecological, and add to data.frame
averageValueMat_TEG <- cbind(averageValueMat_TEG, abs(averageValueMat_TEG$Total-averageValueMat_TEG$Geo))
averageValueMat_TEG <- cbind(averageValueMat_TEG, abs(averageValueMat_TEG$Total-averageValueMat_TEG$Eco))
names(averageValueMat_TEG) <- c(names(averageValueMat_TEG)[1:3], 'Geo_Difference', 'Eco_Difference')

# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','darkblue', 'purple')
plotColors_Fade <- alpha(plotColors, 0.45)
plotColors_Sub <- plotColors_Fade[-(2:4)]
# Two plots in a single window
par(mfrow=c(2,1))
# ---- CORRELATION PLOTS
plot(averageValueMat_TEG$Geo, averageValueMat_TEG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     main='Y. brevifolia: Geographic by genetic coverage',xlab='', ylab='',
     col=plotColors_Fade[[5]])
mtext(text='319 Individuals; 1 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Geographic/Ecological coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add points for ecological coverage
points(x=averageValueMat$Eco, y=averageValueMat$Total, pch=20, col=plotColors_Fade[[6]])
# Add NRMSE values for each comparison
text(x = 76, y = 22, labels = paste0('NRMSE: ', YUBR_nrmse_geo), col='darkblue', cex=0.9)
text(x = 76, y = 8, labels = paste0('NRMSE: ', YUBR_nrmse_eco), col='purple', cex=0.9)
# Add legend
legend(x=57, y=41, inset = 0.05, xpd=TRUE,
       legend = c('Geographic', 'Ecological'),
       col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 0.8)
# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TEG[,1:3], ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='Y. brevifolia: Geo-Eco-Gen Coverage', line=1.5)
mtext(text='319 Individuals; 1 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=200, y=60, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (1 km buffer)',
                  'Ecological coverage (1 km buffer, EPA Level IV)'),
       col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 0.8)
# ---- DIFFERENCE PLOTS
# Plot difference between geographic and genetic coverage
matplot(averageValueMat_TEG[4:5], col=plotColors[5:6], pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='Y. brevifolia: Genetic-Geographic-Ecological Coverage Difference', line=1.5)
mtext(text='319 Individuals; 1 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.2)
mtext(text='Difference in Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=200, y=35, inset = 0.05,
       legend = c('Genographic coverage difference', 'Ecological coverage difference'),
       col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 1)

# %%%% SDM AND TOTAL BUFFER COMPARISON ----
# Specify filepath for MIGU geographic and genetic data, including resampling array
YUBR_filePath <- paste0(GeoGenCorr_wd, 'Datasets/YUBR/')
arrayDir <- paste0(YUBR_filePath, 'resamplingData/YUBR_1km_G2E_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
YUBR_geoComp_1km_array <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values
YUBR_geoComp_1km_DF <- resample.array2dataframe(YUBR_geoComp_1km_array)
# Calculate normalized root mean square value
YUBR_nrmse_geo_totalBuff <-
  nrmse_func(obs=YUBR_geoComp_1km_DF$Geo_Buff, pred=YUBR_geoComp_1km_DF$Total) ; YUBR_nrmse_geo_totalBuff
YUBR_nrmse_geo_SDM <-
  nrmse_func(obs=YUBR_geoComp_1km_DF$Geo_SDM, pred=YUBR_geoComp_1km_DF$Total) ; YUBR_nrmse_geo_SDM

# ---- PLOTTING
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic and geographic proportions
YUBR_geoComp_1km_averageValueMat <- meanArrayValues(YUBR_geoComp_1km_array)
# Calculate the absolute difference between genetic and geographic approaches, and add to data.frame
YUBR_geoComp_1km_averageValueMat <-
  cbind(YUBR_geoComp_1km_averageValueMat, abs(YUBR_geoComp_1km_averageValueMat$Total-YUBR_geoComp_1km_averageValueMat$Geo_Buff))
YUBR_geoComp_1km_averageValueMat <-
  cbind(YUBR_geoComp_1km_averageValueMat, abs(YUBR_geoComp_1km_averageValueMat$Total-YUBR_geoComp_1km_averageValueMat$Geo_SDM))
names(YUBR_geoComp_1km_averageValueMat) <-
  c(names(YUBR_geoComp_1km_averageValueMat)[1:4],'Geo_Buff_Difference', 'Geo_SDM_Difference')
# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue')
plotColors_fade <- alpha(c('red','red4','darkorange3','coral','purple', 'darkblue'), 0.45)
# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1))

# ---- CORRELATION PLOT
# Plot genetic coverage against geographic coverage using total buffer approach
plot(YUBR_geoComp_1km_averageValueMat$Geo_Buff, YUBR_geoComp_1km_averageValueMat$Total, pch=20,
     main='Y. brevifolia: Geographic by genetic coverage',
     xlab='Geographic coverage (%)', ylab='Genetic coverage (%)',
     col=plotColors_fade[[2]])
# Add points for SDM approach
points(x=YUBR_geoComp_1km_averageValueMat$Geo_SDM, y=YUBR_geoComp_1km_averageValueMat$Total,
       pch=20, col=plotColors_fade[[3]])
# Subtitle
mtext(text='319 Individuals; 1 km buffer; 5 replicates', side=3, line=0.3)
# Add NRMSE values for each comparison
text(x = 86, y = 59, labels = paste0('NRMSE: ', YUBR_nrmse_geo_totalBuff), col='red4', cex=0.9)
text(x = 86, y = 50, labels = paste0('NRMSE: ', YUBR_nrmse_geo_SDM), col='darkorange3', cex=0.9)
# Add legend
legend(x=65, y=69, inset = 0.05, xpd=TRUE,
       legend = c('Total buffer approach', 'SDM approach'),
       col=plotColors[2:3], pch = c(20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.8)
# ---- COVERAGE PLOT
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(YUBR_geoComp_1km_averageValueMat[,1:3], ylim=c(0,100), col=plotColors_fade,
        pch=16, ylab='Coverage (%)')
# Add title and x-axis labels to the graph
title(main='Y. brevifolia: Coverage Values by Sample Size', line=1.5)
mtext(text='319 Individuals; 1 km buffer; 5 replicates', side=3, line=0.3)
mtext(text='Number of individuals', side=1, line=2.4)
# Add legend
legend(x=200, y=85, inset = 0.05, xpd=TRUE,
       legend = c('Genetic coverage', 'Geographic, Total buffer (1 km)', 'Geographic, SDM (1 km)'),
       col=plotColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.8)

# %%%% SMBO: MULTIPLE BUFFER SIZES ----
# Specify filepath for YUBR geographic and genetic data, including resampling array
YUBR_filePath <- paste0(GeoGenCorr_wd, 'Datasets/YUBR/')
arrayDir <- paste0(YUBR_filePath, 'resamplingData/YUBR_SMBO3_G2G2E_resampArr.Rdata')
# Read in array and build a data.frame of values
YUBR_SMBO3_array <- readRDS(arrayDir)

# ---- CALCULATIONS ----
# Build a data.frame from array values
YUBR_SMBO3_DF <- resample.array2dataframe(YUBR_SMBO3_array)
# Build tables of NRSMSE values, calculated based on data.frame
YUBR_NRMSE_Mat_CV <- buildNRMSEmatrix(resampDF=YUBR_SMBO3_DF, genCovType='CV', sdmFlag=TRUE)
YUBR_NRMSE_Mat_GD <- buildNRMSEmatrix(resampDF=YUBR_SMBO3_DF, genCovType='GD', sdmFlag=TRUE)
# Combine the results of the NRMSE values calculated using allelic coverages and using
# genetic distances, and then rename the columns accordingly
YUBR_NRMSE_Mat <- cbind(YUBR_NRMSE_Mat_CV, YUBR_NRMSE_Mat_GD)
# Store the matrix as a CSV to disk
write.table(YUBR_NRMSE_Mat,
            file=paste0(YUBR_filePath, 'resamplingData/YUBR_SMBO3_NRMSE.csv'), sep=',')

# ---- PLOTTING ----
# Specify plot colors
plotColors <- colorRampPalette(c("darkred","azure4","lightgray"))(129)
# Build a matrix of mean values (based on array; custom function returns a data.frame)
YUBR_SMBO2_meanValues <- as.matrix(meanArrayValues(YUBR_MultBuff_array))
# Subset the matrix of mean values according to each coverage type
YUBR_SMBO2_GeoBuffMeans <- YUBR_SMBO2_meanValues[,grep('Geo_Buff',colnames(YUBR_SMBO2_meanValues))]
YUBR_SMBO2_GeoSDMMeans <- YUBR_SMBO2_meanValues[,grep('Geo_SDM',colnames(YUBR_SMBO2_meanValues))]
YUBR_SMBO2_EcoBuffMeans <- YUBR_SMBO2_meanValues[,grep('Eco_Buff',colnames(YUBR_SMBO2_meanValues))]
# Create colors based on the NRMSE values in matrix. Make all points transparent (alpha)
GeoBuffCols <-
  alpha(plotColors[as.numeric(cut(YUBR_NRMSE_Mat[,1], breaks = length(plotColors)))], 0.15)
GeoSDMCols <-
  alpha(plotColors[as.numeric(cut(YUBR_NRMSE_Mat[,2], breaks = length(plotColors)))], 0.15)
EcoBuffCols <-
  alpha(plotColors[as.numeric(cut(YUBR_NRMSE_Mat[,3], breaks = length(plotColors)))], 0.15)
# For each color vector, decrease transparency of points corresponding to the lowest NRMSE
GeoBuffCols[[which.min(YUBR_NRMSE_Mat[,1])]] <-
  alpha(GeoBuffCols[[which.min(YUBR_NRMSE_Mat[,1])]], 0.65)
GeoSDMCols[[which.min(YUBR_NRMSE_Mat[,2])]] <-
  alpha(GeoSDMCols[[which.min(YUBR_NRMSE_Mat[,2])]], 0.65)
EcoBuffCols[[which.min(YUBR_NRMSE_Mat[,3])]] <-
  alpha(EcoBuffCols[[which.min(YUBR_NRMSE_Mat[,3])]], 0.65)

# Use matplot to plot values for different coverages
# GeoBuff
matplot(YUBR_SMBO2_GeoBuffMeans, ylim=c(0,100), col=GeoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Y. brevifolia: Geographic Coverages (Total buffer)')
# Add points for genetic values, subtitle, optimal buffer size, and legend
points(YUBR_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
mtext(text='319 Individuals; 41 buffer sizes (0.5km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal geographic buffer size: 4 km', side=1, line=-2, at=60, cex=1.1)
legend(x=550, y=55, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Genetic values'),
       y.intersp = 0.75)
# GeoSDM
matplot(YUBR_SMBO2_GeoSDMMeans, ylim=c(0,100), col=GeoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Y. brevifolia: Geographic Coverages (SDM)')
# Add points for genetic values, subtitle, optimal buffer size, and legend
points(YUBR_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
mtext(text='319 Individuals; 41 buffer sizes (0.5km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal geographic buffer size: 20 km', side=1, line=-1.7, at=60, cex=1.1)
legend(x=550, y=55, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Genetic values'),
       y.intersp = 0.75)
# EcoBuff
matplot(YUBR_SMBO2_EcoBuffMeans, ylim=c(0,100), col=EcoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Y. brevifolia: Ecological Coverages')
# Add points for genetic values, subtitle, optimal buffer size, and legend
points(YUBR_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
mtext(text='319 Individuals; 41 buffer sizes (0.5km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal ecological buffer size: 25 km', side=1, line=-1.7, at=60, cex=1.1)
legend(x=550, y=55, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Genetic values'),
       y.intersp = 0.75)
