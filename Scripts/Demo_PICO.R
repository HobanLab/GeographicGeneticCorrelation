# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: PINUS CONTORTA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage
# for Pinus contorta. Uses data files from MacLachlan et al. 2021 to pull in genetic data 
# (as a STRUCTURE input file) and geographic coordinates (included in a CSV) to conduct 
# correlation analyses. We also process the geographic coordinates to match the order of 
# the genetic samples, prior to passing both along to the correlation analyses.

library(adegenet)
library(terra)
library(parallel)
library(RColorBrewer)
library(scales)

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
# Make sure libraries (adegenet + terra) are on cluster
clusterEvalQ(cl, library('adegenet'))
clusterEvalQ(cl, library('terra'))
clusterEvalQ(cl, library('parallel'))

# %%% CONDUCT RESAMPLING %%% ----
# ---- READ IN DATA ----
# Specify filepath for PICO geographic and genetic data
PICO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/PICO/')

# ---- GENETIC MATRIX
# The 1st script in the MacLachlan et al. 2021 supplement 
# (1_MacLachlan_etal_Pine_GPA_ped&mapfile_formatting_Jan10th2021.R) generates a .ped and .map file. 
# This is for 929 individuals (control individuals are not included), and it includes 32,449 loci, 
# after filtering for minor alleles and missing data. The .ped file was then are passed into PLINK 
# in order to generate a STRUCTURE input file, which is converted into a genind object in the call below.

# (Because the STRUCTURE input file is greater than 50 MB in size, it is not included on the GitHub repository)
PICO_genind <- read.structure(file=paste0(PICO_filePath, 'Genetic/Pine_NaturalComp_85SNPfilter.stru'), 
                              n.ind = 929, n.loc = 32449, onerowperind = TRUE, col.lab = 1, 
                              col.pop = 0, row.marknames = 0, sep = ' ', ask = FALSE)
# Capture the names of samples as they're ordered in the genind object. This is for the steps below.
PICO_sampleNames <- indNames(PICO_genind)

# ---- GEOGRAPHIC/ECOLOGICAL DATA FILES
# The supplement of MacLachlan et al. 2021 includes a csv that contains climate data for all of the 
# seedlings in the study, as well as coordinate information. This dataset needs to first be subset 
# just to the 929 samples included in the genind object (above), and then ordered to match the order 
# of samples in that genind object. These steps are taken below.
PICO_coordinates <- 
  read.csv(file=paste0(PICO_filePath, 'Geographic/Pine_NaturalOrchard_ClimateData_Sept7th2015.csv'), header = TRUE)
# Start by subsetting the CSV to just the variables we need: sample names, latitude, and longitude
PICO_coordinates <- 
  PICO_coordinates[which(PICO_coordinates$Internal_ID %in% PICO_sampleNames),c('Internal_ID','Latitude','Longitude')]
# Reorder the coordinate values to match the order of samples in the genind file
PICO_coordinates <- PICO_coordinates[order(match(PICO_coordinates$Internal_ID, PICO_sampleNames)),]
# Rename the columns of the geographic coordinates data.frame (because geo.compareBuff function 
# expects certain strings)
colnames(PICO_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')
# Read in raster data, for SDM
PICO_sdm <- terra::rast(paste0(PICO_filePath,'Geographic/PICO_929inds_rast_Carver.tif'))

# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.gpkg')))
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = PICO_coordinates)
# Read in the EPA Level III ecoregion shapefile, which is used for calculating ecological coverage 
# (in North America)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level3/NA_CEC_Eco_Level3.shp')))
# Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
# exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
world_poly_clip_W <- wrap(world_poly_clip)
ecoregion_poly_W <- wrap(ecoregion_poly)
PICO_sdm_W <- wrap(PICO_sdm)

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('PICO_coordinates','PICO_genind','num_reps','geo_buffSize', 'eco_buffSize',
                              'world_poly_clip_W', 'ecoregion_poly_W', 'PICO_sdm_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres', 
                              'eco.intersectBuff','eco.compareBuff','gen.getAlleleCategories',
                              'eco.totalEcoregionCount','calculateCoverage','exSituResample.Par',
                              'geo.gen.Resample.Par'))
# Specify file path, for saving resampling array
arrayDir <- paste0(PICO_filePath, 'resamplingData/PICO_SMBO2_GE_5r_resampArr.Rdata')

# Run resampling (in parallel)
PICO_demoArray_Par <- 
  geo.gen.Resample.Par(genObj = PICO_genind, geoFlag = TRUE, coordPts = PICO_coordinates, 
                       geoBuff = geo_buffSize, SDMrast = NA, boundary=world_poly_clip_W, 
                       ecoFlag = TRUE, ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, 
                       ecoLayer = 'NA', reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
# Close cores
stopCluster(cl)

# Run resampling not in parallel (for function testing purposes)
# PICO_demoArray_IND <-
#   geo.gen.Resample(gen_obj = PICO_genind, SDMrast = PICO_sdm, geoFlag = TRUE, coordPts = PICO_coordinates,
#                    geoBuff = geo_buffSize, boundary = world_poly_clip, ecoFlag = TRUE, ecoBuff=eco_buffSize,
#                    ecoRegions = ecoregion_poly, ecoLayer = 'NA', reps = 1)

# %%% ANALYZE DATA %%% ----
# Specify filepath for PICO geographic and genetic data, including resampling array
PICO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/PICO/')
arrayDir <- paste0(PICO_filePath, 'resamplingData/PICO_50km_GE_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
PICO_demoArray_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values
PICO_DF <- resample.array2dataframe(PICO_demoArray_Par)
# Calculate normalized root mean square value
PICO_nrmse_geo <- nrmse_func(obs=PICO_DF$Geo, pred=PICO_DF$Total) ; PICO_nrmse_geo
PICO_nrmse_eco <- nrmse_func(obs=PICO_DF$Eco, pred=PICO_DF$Total) ; PICO_nrmse_eco

# ---- PLOTTING ----
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(PICO_demoArray_Par, allValues = TRUE)
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
     main='P. contorta: Geographic by genetic coverage',xlab='', ylab='',
     col=plotColors_Fade[[5]])
mtext(text='929 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Geographic/Ecological coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add points for ecological coverage
points(x=averageValueMat$Eco, y=averageValueMat$Total, pch=20, col=plotColors[[6]])
# Add Spearman's r values for each comparison
text(x = 76, y = 35, labels = paste0('NRMSE: ', PICO_nrmse_geo), col='darkblue', cex=0.9)
text(x = 76, y = 20, labels = paste0('NRMSE: ', PICO_nrmse_eco), col='purple', cex=0.9)
# Add legend
legend(x=57, y=53, inset = 0.05, xpd=TRUE,
       legend = c('Geographic', 'Ecological'),
       col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 0.8)
# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TEG[,1:3], ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='P. contorta: Geo-Eco-Gen Coverage', line=1.5)
mtext(text='929 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=580, y=70, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (50 km buffer)',
                  'Ecological coverage (50 km buffer, EPA Level III)'),
       col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 0.8)

# %%%% SDM AND TOTAL BUFFER COMPARISON ----
# Specify filepath for PICO geographic and genetic data, including resampling array
PICO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/PICO/')
arrayDir <- paste0(PICO_filePath, 'resamplingData/PICO_50km_G2E_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
PICO_geoComp_50km_array <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values
PICO_geoComp_50km_DF <- resample.array2dataframe(PICO_geoComp_50km_array)
# Calculate normalized root mean square value
PICO_nrmse_geo_totalBuff <-
  nrmse_func(obs=PICO_geoComp_50km_DF$Geo_Buff, pred=PICO_geoComp_50km_DF$Total) ; PICO_nrmse_geo_totalBuff
PICO_nrmse_geo_SDM <-
  nrmse_func(obs=PICO_geoComp_50km_DF$Geo_SDM, pred=PICO_geoComp_50km_DF$Total) ; PICO_nrmse_geo_SDM

# ---- PLOTTING
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic and geographic proportions
PICO_geoComp_50km_averageValueMat <- meanArrayValues(PICO_geoComp_50km_array)
# Calculate the absolute difference between genetic and geographic approaches, and add to data.frame
PICO_geoComp_50km_averageValueMat <-
  cbind(PICO_geoComp_50km_averageValueMat, abs(PICO_geoComp_50km_averageValueMat$Total-PICO_geoComp_50km_averageValueMat$Geo_Buff))
PICO_geoComp_50km_averageValueMat <-
  cbind(PICO_geoComp_50km_averageValueMat, abs(PICO_geoComp_50km_averageValueMat$Total-PICO_geoComp_50km_averageValueMat$Geo_SDM))
names(PICO_geoComp_50km_averageValueMat) <-
  c(names(PICO_geoComp_50km_averageValueMat)[1:4],'Geo_Buff_Difference', 'Geo_SDM_Difference')
# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue')
plotColors_fade <- alpha(c('red','red4','darkorange3','coral','purple', 'darkblue'), 0.45)
# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1))

# ---- CORRELATION PLOT
# Plot genetic coverage against geographic coverage using total buffer approach
plot(PICO_geoComp_50km_averageValueMat$Geo_Buff, PICO_geoComp_50km_averageValueMat$Total, pch=20,
     main='P. contorta: Geographic by genetic coverage',
     xlab='Geographic coverage (%)', ylab='Genetic coverage (%)',
     col=plotColors_fade[[2]])
# Add points for SDM approach
points(x=PICO_geoComp_50km_averageValueMat$Geo_SDM, y=PICO_geoComp_50km_averageValueMat$Total,
       pch=20, col=plotColors_fade[[3]])
# Subtitle
mtext(text='929 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3)
# Add NRMSE values for each comparison
text(x = 82, y = 62, labels = paste0('NRMSE: ', PICO_nrmse_geo_totalBuff), col='red4', cex=0.9)
text(x = 82, y = 56, labels = paste0('NRMSE: ', PICO_nrmse_geo_SDM), col='darkorange3', cex=0.9)
# Add legend
legend(x=57, y=71, inset = 0.05, xpd=TRUE,
       legend = c('Total buffer approach', 'SDM approach'),
       col=plotColors[2:3], pch = c(20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.8)
# ---- COVERAGE PLOT
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(PICO_geoComp_50km_averageValueMat[,1:3], ylim=c(0,100), col=plotColors_fade,
        pch=16, ylab='Coverage (%)')
# Add title and x-axis labels to the graph
title(main='P. contorta: Coverage Values by Sample Size', line=1.5)
mtext(text='929 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3)
mtext(text='Number of individuals', side=1, line=2.4)
# Add legend
legend(x=500, y=54, inset = 0.05, xpd=TRUE,
       legend = c('Genetic coverage', 'Geographic, Total buffer (50 km)', 'Geographic, SDM (50 km)'),
       col=plotColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.8)

# %%%% SMBO: MULTIPLE BUFFER SIZES ----
# Specify filepath for PICO geographic and genetic data, including resampling array
PICO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/PICO/')
arrayDir <- paste0(PICO_filePath, 'resamplingData/PICO_SMBO2_GE_5r_resampArr.Rdata')
# Read in array and build a data.frame of values
PICO_SMBO2_array <- readRDS(arrayDir)
# Specify geographic buffer size in meters (used above)
geo_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500))

# ---- CALCULATIONS ----
# Build a data.frame from array values
PICO_SMBO2_DF <- resample.array2dataframe(PICO_SMBO2_array)
# Build a matrix to capture NRMSE values
PICO_NRMSE_Mat <- matrix(NA, nrow=length(geo_buffSize), ncol=2)
# The names of this matrix match the different parts of the dataframe names
colnames(PICO_NRMSE_Mat) <- c('Geo_Buff','Eco_Buff')
rownames(PICO_NRMSE_Mat) <- paste0(geo_buffSize/1000, 'km')
# Loop through the dataframe columns. The first two columns are skipped, as they're sampleNumber and the
# predictve variable (genetic coverages)
for(i in 3:ncol(PICO_SMBO2_DF)){
  # Calculate NRMSE for the current column in the dataframe
  PICO_NRMSEvalue <- nrmse.func(PICO_SMBO2_DF[,i], pred = PICO_SMBO2_DF$Total)
  # Get the name of the current dataframe column
  dataName <- unlist(strsplit(names(PICO_SMBO2_DF)[[i]],'_'))
  # Match the data name to the relevant rows/columns of the receiving matrix
  matRow <- which(rownames(PICO_NRMSE_Mat) == dataName[[3]])
  matCol <- which(colnames(PICO_NRMSE_Mat) == paste0(dataName[[1]],'_',dataName[[2]]))
  # Locate the NRMSE value accordingly
  PICO_NRMSE_Mat[matRow,matCol] <- PICO_NRMSEvalue
}
print(PICO_NRMSE_Mat)
# Store the matrix as a CSV to disk
write.table(PICO_NRMSE_Mat,
            file=paste0(PICO_filePath, 'resamplingData/PICO_SMBO2_NRMSE.csv'), sep=',')
