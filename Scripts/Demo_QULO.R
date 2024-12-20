# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: QUERCUS LOBATA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script demonstrating approach for calculating the correlation 
# between genetic, geographic, and ecological coverage. 
# Uses data files from Gugger et al. 2021 to pull in genetic data 
# (as a table) and geographic coordinates (in a CSV) to conduct correlation analyses.

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

# ---- READ IN DATA ----
# Specify filepath for QULO geographic and genetic data
QULO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QULO/')

# ---- GENETIC MATRIX
# Read in the SNPs genetic matrix file
QULO_tab <- read.table(paste0(QULO_filePath, 'Genetic/SNPs80.forR'), header = TRUE)
# Make sample names the row names
rownames(QULO_tab) <- QULO_tab[,1] ; QULO_tab <- QULO_tab[,-1]
# Convert to genind. ncode based off similar practice for other geninds 
# (values of 1, 2, and 3 generate identical results)
QULO_genind <- df2genind(QULO_tab, ncode = 3, ploidy = 2)

# ---- GEOGRAPHIC/ECOLOGICAL DATA FILES
# Read in wild occurrence points. This CSV has 3 columns: sample name, latitude, and longitude. 
# The sample names (and order) have to match the sample names/order of the genind object 
# (rownams of the genetic matrix) read in below.
QULO_points <- read.csv(paste0(QULO_filePath, 'Geographic/QULO_coordinates.csv'), header=TRUE)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- grabWorldAdmin(GeoGenCorr_wd = GeoGenCorr_wd, fileExtentsion = ".gpkg", overwrite = FALSE)
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = QULO_points)
# Read in raster data, for SDM
QULO_sdm <- terra::rast(paste0(QULO_filePath,'Geographic/QULO_436inds_rast_Carver.tif'))
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage 
# (solely in the U.S.)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level4/us_eco_l4.shp')))

# ---- PARALLELIZATION
# Flag for running resampling steps in parallel
parFlag <- TRUE

# If running in parallel, set up cores and export required libraries
if(parFlag==TRUE){
  # Set up relevant cores 
  num_cores <- detectCores() - 4 
  cl <- makeCluster(num_cores)
  # Make sure libraries (adegenet, terra, etc.) are on cluster (but avoid printing output)
  invisible(clusterEvalQ(cl, library('adegenet')))
  invisible(clusterEvalQ(cl, library('terra')))
  invisible(clusterEvalQ(cl, library('parallel')))
  invisible(clusterEvalQ(cl, library('usedist')))
  invisible(clusterEvalQ(cl, library('ape')))
  # Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
  # exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
  world_poly_clip_W <- wrap(world_poly_clip)
  ecoregion_poly_W <- wrap(ecoregion_poly)
  QULO_sdm_W <- wrap(QULO_sdm)
}

# ---- RESAMPLING ----
if(parFlag==TRUE){
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('QULO_points','QULO_genind','QULO_sdm_W', 'num_reps','geo_buffSize', 
                              'eco_buffSize', 'world_poly_clip_W', 'ecoregion_poly_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres', 
                              'eco.intersectBuff','eco.compareBuff','gen.getAlleleCategories', 
                              'gen.buildDistMat', 'gen.calcGenDistCov', 'eco.totalEcoregionCount',
                              'calculateCoverage','exSituResample.Par', 'geo.gen.Resample.Par'))
# Specify file path, for saving resampling array
arrayDir <- paste0(QULO_filePath, 'resamplingData/QULO_SMBO3_G2G2E_5r_resampArr.Rdata')
# Run resampling (in parallel)
QULO_demoArray_Par <- 
  geo.gen.Resample.Par(genObj = QULO_genind, genDistFlag=TRUE, geoFlag = TRUE, coordPts = QULO_points, 
                       geoBuff = geo_buffSize, SDMrast=QULO_sdm_W, boundary=world_poly_clip_W, 
                       ecoFlag = TRUE, ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, 
                       ecoLayer = 'US', reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
# Close cores
stopCluster(cl)
} else {
  # Reduce the number of samples, to allow for more efficient testing
  randomSamp <- sample(indNames(QULO_genind), size=430)
  # Subset geographic coordinate dataframe
  QULO_points_small <- QULO_points[which(QULO_points$sampleID %in% randomSamp),]
  # Subset genind object, based on smaller geo coordinates dataframe
  QULO_genind_small <- QULO_genind[QULO_points_small[,1], drop=TRUE]
  # Run resampling not in parallel (for function testing purposes)
  QULO_demoArray_IND <-
    geo.gen.Resample(genObj=QULO_genind_small, genDistFlag=TRUE, SDMrast=NA, geoFlag=TRUE, 
                     coordPts=QULO_points_small, geoBuff=geo_buffSize, boundary=world_poly_clip, 
                     ecoFlag=FALSE, ecoBuff=eco_buffSize, ecoRegions=ecoregion_poly, ecoLayer='US', reps=1)
}

# %%% ANALYZE DATA %%% ----
# Specify filepath for QULO geographic and genetic data, including resampling array
QULO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QULO/')
arrayDir <- paste0(QULO_filePath, 'resamplingData/QULO_50km_GE_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
QULO_demoArray_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values
QULO_DF <- resample.array2dataframe(QULO_demoArray_Par)
# Calculate normalized root mean square value
QULO_nrmse_geo <- nrmse_func(obs=QULO_DF$Geo, pred=QULO_DF$Total) ; QULO_nrmse_geo
QULO_nrmse_eco <- nrmse_func(obs=QULO_DF$Eco, pred=QULO_DF$Total) ; QULO_nrmse_eco

# ---- PLOTTING ----
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(QULO_demoArray_Par, allValues = TRUE)
# Subset matrix of all average values to just Total allelic, geographic, and ecological coverage
averageValueMat_TEG <- averageValueMat[,c(1,6,7)]
# Calculate the absolute difference between genetic and geographic/ecological, and add to data.frame
averageValueMat_TEG <- cbind(averageValueMat_TEG, abs(averageValueMat_TEG$Total-averageValueMat_TEG$Geo))
averageValueMat_TEG <- cbind(averageValueMat_TEG, abs(averageValueMat_TEG$Total-averageValueMat_TEG$Eco))
names(averageValueMat_TEG) <- c(names(averageValueMat_TEG)[1:3], 'Geo_Difference', 'Eco_Difference')
# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','darkblue', 'purple')
plotColors <- alpha(plotColors, 0.45)
plotColors_Sub <- plotColors[-(2:4)]
# Two plots in a single window
par(mfrow=c(2,1))

# ---- CORRELATION PLOTS
plot(averageValueMat_TEG$Geo, averageValueMat_TEG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     main='Q. lobata: Geographic/Ecological by genetic coverage',xlab='', ylab='', col=plotColors[[5]])
mtext(text='436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Geographic/Ecological coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add points for ecological coverage
points(x=averageValueMat$Eco, y=averageValueMat$Total, pch=20, col=plotColors[[6]])
# Add NRMSE values for each comparison
text(x = 76, y = 35, labels = paste0('NRMSE: ', QULO_nrmse_geo), col='darkblue', cex=0.9)
text(x = 76, y = 20, labels = paste0('NRMSE: ', QULO_nrmse_eco), col='purple', cex=0.9)
# Add legend
legend(x=57, y=53, inset = 0.05, xpd=TRUE,
       legend = c('Geographic', 'Ecological'),
       col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.8)
# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TEG[,1:3], ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='Q. lobata: Geo-Eco-Gen Coverage', line=1.5)
mtext(text='436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=275, y=60, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (50 km buffer)',
                  'Ecological coverage (50 km buffer, EPA Level IV)'),
       col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 0.8)

# %%%% SDM AND TOTAL BUFFER COMPARISON ----
# Specify filepath for QULO geographic and genetic data, including resampling array
QULO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QULO/')
arrayDir <- paste0(QULO_filePath, 'resamplingData/QULO_50km_G2E_LozaMX_5r_resampArr.Rdata')
# Read in array and build a data.frame of values
QULO_geoComp_50km_array <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values
QULO_geoComp_50km_DF <- resample.array2dataframe(QULO_geoComp_50km_array)
# Calculate normalized root mean square value
QULO_nrmse_geo_totalBuff <-
  nrmse_func(obs=QULO_geoComp_50km_DF$Geo_Buff, pred=QULO_geoComp_50km_DF$Total) ; QULO_nrmse_geo_totalBuff
QULO_nrmse_geo_SDM <-
  nrmse_func(obs=QULO_geoComp_50km_DF$Geo_SDM, pred=QULO_geoComp_50km_DF$Total) ; QULO_nrmse_geo_SDM

# ---- PLOTTING
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic and geographic proportions
QULO_geoComp_50km_averageValueMat <- meanArrayValues(QULO_geoComp_50km_array)
# Calculate the absolute difference between genetic and geographic approaches, and add to data.frame
QULO_geoComp_50km_averageValueMat <-
  cbind(QULO_geoComp_50km_averageValueMat, abs(QULO_geoComp_50km_averageValueMat$Total-QULO_geoComp_50km_averageValueMat$Geo_Buff))
QULO_geoComp_50km_averageValueMat <-
  cbind(QULO_geoComp_50km_averageValueMat, abs(QULO_geoComp_50km_averageValueMat$Total-QULO_geoComp_50km_averageValueMat$Geo_SDM))
names(QULO_geoComp_50km_averageValueMat) <-
  c(names(QULO_geoComp_50km_averageValueMat)[1:4],'Geo_Buff_Difference', 'Geo_SDM_Difference')
# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue')
plotColors_Fade <- alpha(plotColors, 0.45)

# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1))

# ---- CORRELATION PLOT
# Plot genetic coverage against geographic coverage using total buffer approach
plot(QULO_geoComp_50km_averageValueMat$Geo_Buff, QULO_geoComp_50km_averageValueMat$Total, pch=20,
     main='Q. lobata: Geographic by genetic coverage',
     xlab='Geographic coverage (%)', ylab='Genetic coverage (%)',
     col=plotColors_Fade[2])
# Add points for SDM approach
points(x=QULO_geoComp_50km_averageValueMat$Geo_SDM, y=QULO_geoComp_50km_averageValueMat$Total,
       pch=20, col=plotColors_Fade[3])
# Subtitle
mtext(text='SDM: Loza (MaxEnt); 436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.1)
# Add NRMSE values for each comparison
text(x = 80, y = 65, labels = paste0('NRMSE: ', QULO_nrmse_geo_totalBuff), col='red4', cex=0.9)
text(x = 80, y = 58, labels = paste0('NRMSE: ', QULO_nrmse_geo_SDM), col='darkorange3', cex=0.9)
# Add legend
legend(x=58, y=74, inset = 0.05, xpd=TRUE,
       legend = c('Total buffer approach', 'SDM approach'),
       col=plotColors[2:3], pch = c(20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 0.8)
# ---- COVERAGE PLOT
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(QULO_geoComp_50km_averageValueMat[,1:3], ylim=c(0,100), col=plotColors_Fade,
        pch=16, ylab='Coverage (%)')
# Add title and x-axis labels to the graph
title(main='Q. lobata: Coverage Values by Sample Size', line=1.5)
mtext(text='SDM: Loza (MaxEnt); 436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3)
mtext(text='Number of individuals', side=1, line=2.4)
# Add legend
legend(x=300, y=63, inset = 0.05, xpd=TRUE,
       legend = c('Genetic coverage', 'Geographic, Total buffer (50 km)', 'Geographic, SDM (50 km)'),
       col=plotColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.8)
# # ---- MAP OF SDM AND SAMPLED POINTS
# makeAMap(QULO_points, raster = QULO_sdm, buffer = geo_buffSize)

# %%%% SMBO: MULTIPLE BUFFER SIZES ----
# Specify filepath for QULO geographic and genetic data, including resampling array
QULO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QULO/')
arrayDir <- paste0(QULO_filePath, 'resamplingData/QULO_SMBO3_G2G2E_5r_resampArr.Rdata')
# Read in array
QULO_SMBO3_array <- readRDS(arrayDir)

# ---- CALCULATIONS ----
# Build a data.frame from array values
QULO_SMBO3_DF <- resample.array2dataframe(QULO_SMBO3_array)
# Build tables of NRSMSE values, calculated based on data.frame
QULO_NRMSE_Mat_CV <- buildNRMSEmatrix(resampDF=QULO_SMBO3_DF, genCovType='CV', sdmFlag=TRUE)
QULO_NRMSE_Mat_GD <- buildNRMSEmatrix(resampDF=QULO_SMBO3_DF, genCovType='GD', sdmFlag=TRUE)
# Combine the results of the NRMSE values calculated using allelic coverages and using
# genetic distances, and then rename the columns accordingly
QULO_NRMSE_Mat <- cbind(QULO_NRMSE_Mat_CV, QULO_NRMSE_Mat_GD)
# Store the matrix as a CSV to disk
write.table(QULO_NRMSE_Mat,
            file=paste0(QULO_filePath, 'resamplingData/QULO_SMBO3_NRMSE.csv'), sep=',')

# ---- PLOTTING ----
# SMBO2 ----
# Specify plot colors
plotColors <- colorRampPalette(c("darkred","azure4","lightgray"))(129)
# Build a matrix of mean values (based on array; custom function returns a data.frame)
QULO_SMBO2_meanValues <- as.matrix(meanArrayValues(QULO_SMBO2_array))
# Subset the matrix of mean values according to each coverage type
QULO_SMBO2_GeoBuffMeans <- QULO_SMBO2_meanValues[,grep('Geo_Buff',colnames(QULO_SMBO2_meanValues))]
QULO_SMBO2_GeoSDMMeans <- QULO_SMBO2_meanValues[,grep('Geo_SDM',colnames(QULO_SMBO2_meanValues))]
QULO_SMBO2_EcoBuffMeans <- QULO_SMBO2_meanValues[,grep('Eco_Buff',colnames(QULO_SMBO2_meanValues))]
# Create colors based on the NRMSE values in matrix. Make all points transparent (alpha)
GeoBuffCols <-
  alpha(plotColors[as.numeric(cut(QULO_NRMSE_Mat[,1], breaks = length(plotColors)))], 0.15)
GeoSDMCols <-
  alpha(plotColors[as.numeric(cut(QULO_NRMSE_Mat[,2], breaks = length(plotColors)))], 0.15)
EcoBuffCols <-
  alpha(plotColors[as.numeric(cut(QULO_NRMSE_Mat[,3], breaks = length(plotColors)))], 0.15)
# For each color vector, decrease transparency of points corresponding to the lowest NRMSE
GeoBuffCols[[which.min(QULO_NRMSE_Mat[,1])]] <-
  alpha(GeoBuffCols[[which.min(QULO_NRMSE_Mat[,1])]], 0.65)
GeoSDMCols[[which.min(QULO_NRMSE_Mat[,2])]] <-
  alpha(GeoSDMCols[[which.min(QULO_NRMSE_Mat[,2])]], 0.65)
EcoBuffCols[[which.min(QULO_NRMSE_Mat[,3])]] <-
  alpha(EcoBuffCols[[which.min(QULO_NRMSE_Mat[,3])]], 0.65)

# Use matplot to plot values for different coverages
# GeoBuff
matplot(QULO_SMBO2_GeoBuffMeans, ylim=c(0,100), col=GeoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Q. lobata: Geographic Coverages (Total buffer)')
# Add points for genetic values, subtitle, optimal buffer size, and legend
points(QULO_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
mtext(text='436 Individuals; 41 buffer sizes (0.5km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal geographic buffer size: 250 km', side=1, line=-2, at=200, cex=1.1)
legend(x=300, y=55, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Genetic values'),
       y.intersp = 0.75)
# GeoSDM
matplot(QULO_SMBO2_GeoSDMMeans, ylim=c(0,100), col=GeoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Q. lobata: Geographic Coverages (SDM)')
# Add points for genetic values, subtitle, optimal buffer size, and legend
points(QULO_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
mtext(text='436 Individuals; 41 buffer sizes (0.5km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal geographic buffer size: 120 km', side=1, line=-1.7, at=200, cex=1.1)
legend(x=300, y=55, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Genetic values'),
       y.intersp = 0.75)
# EcoBuff
matplot(QULO_SMBO2_EcoBuffMeans, ylim=c(0,100), col=EcoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Q. lobata: Ecological Coverages')
# Add points for genetic values, subtitle, optimal buffer size, and legend
points(QULO_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
mtext(text='436 Individuals; 41 buffer sizes (0.5km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal ecological buffer size: 130 km', side=1, line=-1.7, at=200, cex=1.1)
legend(x=300, y=55, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Genetic values'),
       y.intersp = 0.75)
# Add arrows
arrows(x0=75, y0=40, x1=25, y1=80)

# Update SDM plots: remove data for first few buffer sizes (<10km), to make plots clearer
QULO_SMBO2_GeoSDMMeans <- QULO_SMBO2_GeoSDMMeans[,-(1:6)]
# Create new color vector
GeoSDMCols <-
  alpha(plotColors[as.numeric(cut(QULO_NRMSE_Mat[,2], breaks = length(plotColors)))], 0.15)
GeoSDMCols[[which.min(QULO_NRMSE_Mat[,2])]] <-
  alpha(GeoSDMCols[[which.min(QULO_NRMSE_Mat[,2])]], 0.65)

# Call new plot
matplot(QULO_SMBO2_GeoSDMMeans, ylim=c(0,100), col=GeoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Q. lobata: Geographic Coverages (SDM)')
# Add points for genetic values, subtitle, optimal buffer size, and legend
points(QULO_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
mtext(text='436 Individuals; 35 buffer sizes (10km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal geographic buffer size: 120 km', side=1, line=-1.7, at=200, cex=1.1)
legend(x=300, y=55, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Genetic values'),
       y.intersp = 0.75)

# SMBO2: OPTIMAL BUFFER SIZES ----
# Read in QULO SMBO2 resampling array amd convert to data.frame
QULO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QULO/')
QULO_arrayDir <- paste0(QULO_filePath, 'resamplingData/SMBO2/QULO_SMBO2_G2E_5r_resampArr.Rdata')
# From QULO resampling array, return a matrix of average coverage values for optimal buffer sizes
QULO_optCovMat <- extractOptCovs(QULO_arrayDir)
# Calculate MSSEs: minimum number of samples for 95% of each coverage type
QULO_Gen_MSSE <- min(which(QULO_optCovMat[,1] > 95)) ; QULO_Gen_MSSE
QULO_GeoBuff_MSSE <- min(which(QULO_optCovMat[,2] > 95)) ; QULO_GeoBuff_MSSE
QULO_GeoSDM_MSSE <- min(which(QULO_optCovMat[,3] > 95)) ; QULO_GeoSDM_MSSE
QULO_Eco_MSSE <- min(which(QULO_optCovMat[,4] > 95)) ; QULO_Eco_MSSE

# PLOTTING
# Specify plot colors
plotColors <- c('red', 'darkblue','darkorange3', 'purple')
plotColors_Fade <- alpha(plotColors, c(0.45, rep(0.85, length(plotColors)-1)))
# Set plotting window to stack 3 graphs vertically
par(mfcol=c(3,1), oma=rep(0.2,4), mar=c(2,4,3,1)+0.1)
# Geo Buff
matplot(QULO_optCovMat[,c(1,2)], ylim=c(0,110), col=plotColors_Fade[c(1, 2)], pch=16, ylab='')
abline(h=95, col="black", lty=3)
abline(v=QULO_Gen_MSSE, col="red") ; abline(v=QULO_GeoBuff_MSSE, col="darkblue")
mtext(text=paste0('MSSE: ', QULO_Gen_MSSE), side=1, line=-1, at=QULO_Gen_MSSE+15, cex=0.7, col='red')
mtext(text=paste0(' MSSE: ', QULO_GeoBuff_MSSE), line=-1.5, side=1, at=QULO_GeoBuff_MSSE-15, cex=0.7, col='darkblue')
title('Quercus lobata: Coverages at Optimal Buffer Sizes', cex.sub=1.2, line = 2)
mtext(text='Geographic (Total Buffer): 250 km', side=3, at=80, cex=0.8)
# Geo SDM
par(mar=c(2,4,2,1)+0.1)
matplot(QULO_optCovMat[,c(1,3)], ylim=c(0,110), col=plotColors_Fade[c(1, 3)], pch=16, ylab='')
abline(h=95, col="black", lty=3)
abline(v=QULO_Gen_MSSE, col="red") ; abline(v=QULO_GeoSDM_MSSE, col="darkorange3")
mtext(text=paste0('MSSE: ', QULO_Gen_MSSE), side=1, line=-1, at=QULO_Gen_MSSE+15, cex=0.7, col='red')
mtext(text=paste0(' MSSE: ', QULO_GeoSDM_MSSE), line=-1.5, side=1, at=QULO_GeoSDM_MSSE-15, cex=0.7, col='darkorange3')
mtext(text='Geographic (SDM): 120 km', side=3, at=c(80), cex=0.8)
mtext(text="Coverage (%)", side=2, line=2.6, cex=1.2, srt=90)
# Legend
legend(x=300, y=125, xpd=TRUE, cex=1.2, pch=rep(19,4),
       col=c('red','darkblue','darkorange3', 'purple'),
       legend=c('Genetic', 'Geographic (Total Buffer)','Geographic (SDM)', 'Ecological'),
       y.intersp = 0.3, bty='n')
# Eco Buff
par(mar=c(3,4,2,1)+0.1)
matplot(QULO_optCovMat[,c(1,4)], ylim=c(0,110), col=plotColors_Fade[c(1, 4)], pch=16, ylab='')
abline(h=95, col="black", lty=3)
abline(v=QULO_Gen_MSSE, col="red") ; abline(v=QULO_Eco_MSSE, col="purple")
mtext(text=paste0('MSSE: ', QULO_Gen_MSSE), side=1, line=-1, at=QULO_Gen_MSSE+15, cex=0.7, col='red')
mtext(text=paste0(' MSSE: ', QULO_Eco_MSSE), line=-1.5, side=1, at=QULO_Eco_MSSE-15, cex=0.7, col='purple')
mtext(text='Ecological: 130 km', side=3, at=c(80), cex=0.8)

# SMBO3 ----
# ALLELIC AND GENETIC DISTANCE COVERAGE
# Vector of colors for plotting
ExCols <- c('red', 'darkred', 'darkorchid4', 'darkorchid1')
ExCols_Fade <- c(alpha('red', 0.35), alpha('darkred', 0.45), alpha('darkorchid4', 0.45), alpha('darkorchid1', 0.35))
# Build a matrix of mean values (based on array; custom function returns a data.frame)
QULO_SMBO3_meanValues <- as.matrix(meanArrayValues(QULO_SMBO3_array))
# Subset the total mean Values to just Total, GenDist, GeoBuff_0.5km (optimal for Gen Dist),
# and GeoBuff 250 km (optimal for allelic coverage)
QULO_SMBO3_Gen2GeoValues <- QULO_SMBO3_meanValues[,c(1:3,42)]
# Plot allelic coverage, genetic distance coverage, and geographic coverage at optimal buffer size
matplot(QULO_SMBO3_Gen2GeoValues, ylim=c(0,100), col=ExCols_Fade2,
        pch=16, ylab='Coverage (%)', xlab='Number of individuals',
        main='Q. lobata: Allelic and Genetic Distance Coverage')
mtext(text='436 Individuals; 41 buffer sizes (0.5km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
abline(h=95, lty=3)
legend(x=c(320,415), y=c(0,43), inset = 0.05, title='Coverages',
       legend = c('Allelic', 'Genetic distance',
                  'Geographic (Total, 0.5 km, GD*)', 'Geographic (Total, 250 km, CV*)'),
       col=ExCols, pch = c(20,20,20,20), cex=0.9, pt.cex = 2, y.intersp = 0.8)

# ALL SMBO3 RESULTS
# Specify plot colors
plotColors <- colorRampPalette(c("darkred","azure4","lightgray"))(129)
# Build a matrix of mean values (based on array; custom function returns a data.frame)
QULO_SMBO3_meanValues <- as.matrix(meanArrayValues(QULO_SMBO3_array))
# Subset the matrix of mean values according to each coverage type
QULO_SMBO3_GeoBuffMeans <- QULO_SMBO3_meanValues[,grep('Geo_Buff',colnames(QULO_SMBO3_meanValues))]
QULO_SMBO3_GeoSDMMeans <- QULO_SMBO3_meanValues[,grep('Geo_SDM',colnames(QULO_SMBO3_meanValues))]
QULO_SMBO3_EcoBuffMeans <- QULO_SMBO3_meanValues[,grep('Eco_Buff',colnames(QULO_SMBO3_meanValues))]
# Create colors based on the NRMSE values in matrix (allelic coverages). Make all points transparent (alpha)
GeoBuffCols <-
  alpha(plotColors[as.numeric(cut(QULO_NRMSE_Mat[,1], breaks = length(plotColors)))], 0.15)
GeoSDMCols <-
  alpha(plotColors[as.numeric(cut(QULO_NRMSE_Mat[,2], breaks = length(plotColors)))], 0.15)
EcoBuffCols <-
  alpha(plotColors[as.numeric(cut(QULO_NRMSE_Mat[,3], breaks = length(plotColors)))], 0.15)
# For each color vector, decrease transparency of points corresponding to the lowest NRMSE
GeoBuffCols[[which.min(QULO_NRMSE_Mat[,1])]] <-
  alpha(GeoBuffCols[[which.min(QULO_NRMSE_Mat[,1])]], 0.65)
GeoSDMCols[[which.min(QULO_NRMSE_Mat[,2])]] <-
  alpha(GeoSDMCols[[which.min(QULO_NRMSE_Mat[,2])]], 0.65)
EcoBuffCols[[which.min(QULO_NRMSE_Mat[,3])]] <-
  alpha(EcoBuffCols[[which.min(QULO_NRMSE_Mat[,3])]], 0.65)

# Use matplot to plot values for different coverages
# GeoBuff
matplot(QULO_SMBO3_GeoBuffMeans, ylim=c(0,100), col=GeoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Q. lobata: Geographic Coverages (Total buffer)')
# Add points for allelic coverage values, genetic distance values, subtitle, optimal buffer size, and legend
points(QULO_SMBO3_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
points(QULO_SMBO3_meanValues[,2], col=alpha('blue', 0.55), pch=20)
mtext(text='436 Individuals; 41 buffer sizes (0.5km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal geo. buffer: 0.5 km', side=1, line=-10.0, at=77, cex=1.0)
legend(x=59, y=35, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4','blue'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Allelic coverage', 'Genetic distance coverage'),
       y.intersp = 0.75)
# GeoSDM
matplot(QULO_SMBO3_GeoSDMMeans, ylim=c(0,100), col=GeoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Q. lobata: Geographic Coverages (SDM)')
# Add points for genetic values, subtitle, optimal buffer size, and legend
points(QULO_SMBO3_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
points(QULO_SMBO3_meanValues[,2], col=alpha('blue', 0.55), pch=20)
mtext(text='436 Individuals; 41 buffer sizes (0.5km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal geo. buffer: 10 km', side=1, line=-10.0, at=77, cex=1.0)
legend(x=59, y=35, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4','blue'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Allelic coverage', 'Genetic distance coverage'),
       y.intersp = 0.75)
# EcoBuff
matplot(QULO_SMBO3_EcoBuffMeans, ylim=c(0,100), col=EcoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Q. lobata: Ecological Coverages')
# Add points for genetic values, subtitle, optimal buffer size, and legend
points(QULO_SMBO3_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
points(QULO_SMBO3_meanValues[,2], col=alpha('blue', 0.55), pch=20)
mtext(text='436 Individuals; 41 buffer sizes (0.5km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal ecological buffer size: 35 km', side=1, line=-10.0, at=77, cex=1.0)
legend(x=59, y=35, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4','blue'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Allelic coverage', 'Genetic distance coverage'),
       y.intersp = 0.75)
# Add arrows
arrows(x0=75, y0=40, x1=25, y1=80)

# Update SDM plots: remove data for first few buffer sizes (<10km), to make plots clearer
QULO_SMBO3_GeoSDMMeans <- QULO_SMBO3_GeoSDMMeans[,-(1:6)]
# Create new color vector
GeoSDMCols <-
  alpha(plotColors[as.numeric(cut(QULO_NRMSE_Mat[,2], breaks = length(plotColors)))], 0.15)
GeoSDMCols[[which.min(QULO_NRMSE_Mat[,2])]] <-
  alpha(GeoSDMCols[[which.min(QULO_NRMSE_Mat[,2])]], 0.65)

# Call new plot
matplot(QULO_SMBO3_GeoSDMMeans, ylim=c(0,100), col=GeoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Q. lobata: Geographic Coverages (SDM)')
# Add points for genetic values, subtitle, optimal buffer size, and legend
points(QULO_SMBO3_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
mtext(text='436 Individuals; 35 buffer sizes (10km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal geographic buffer size: 35 km', side=1, line=-1.7, at=200, cex=1.1)
legend(x=300, y=55, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Genetic values'),
       y.intersp = 0.75)
