# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: QUERCUS LOBATA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script demonstrating approach for calculating the correlation 
# between genetic, geographic, and ecological coverage. 
# Uses data files from Gugger et al. 2021 to pull in genetic data 
# (as a table) and geographic coordinates (in a CSV) to conduct correlation analyses.

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
num_reps <- 2
# ---- BUFFER SIZES
# Specify geographic buffer size in meters 
geo_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500))
# Specify ecological buffer size in meters 
eco_buffSize <- 1000*(c(5,40))
# eco_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500))

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
QULO_points <- read.csv(paste0(QULO_filePath, 'Geographic/Quercus_lobata.csv'), header=TRUE)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- grabWorldAdmin(GeoGenCorr_wd = GeoGenCorr_wd, fileExtentsion = ".gpkg", overwrite = FALSE)
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = QULO_points)
# Read in raster data, for SDM
QULO_sdm <- terra::rast(paste0(QULO_filePath,'Geographic/QULO_436inds_rast_Loza_MX.tif'))
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage 
# (solely in the U.S.)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level4/us_eco_l4.shp')))

# ---- PARALLELIZATION
# Flag for running resampling steps in parallel
parFlag <- FALSE

# If running in parallel, set up cores and export required libraries
if(parFlag==TRUE){
  # Set up relevant cores 
  num_cores <- detectCores() - 4 
  cl <- makeCluster(num_cores)
  # Make sure libraries (adegenet + terra) are on cluster
  clusterEvalQ(cl, library('adegenet'))
  clusterEvalQ(cl, library('terra'))
  clusterEvalQ(cl, library('parallel'))
  # Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
  # exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
  world_poly_clip_W <- wrap(world_poly_clip)
  ecoregion_poly_W <- wrap(ecoregion_poly)
  QULO_sdm_W <- wrap(QULO_sdm)
}

# ---- RESAMPLING ----
if(parFlag==TRUE){
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('QULO_points','QULO_genind','QULO_sdm', 'num_reps','geo_buffSize', 
                              'eco_buffSize', 'world_poly_clip_W', 'ecoregion_poly_W', 'QULO_sdm_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers', 'geo.compareBuff', 'geo.compareBuffSDM', 'geo.checkSDMres', 
                              'eco.intersectBuff', 'eco.compareBuff', 'gen.getAlleleCategories',
                              'calculateCoverage', 'exSituResample.Par', 'geo.gen.Resample.Par'))
# Specify file path, for saving resampling array
arrayDir <- paste0(QULO_filePath, 'resamplingData/QULO_SMBO2_GE_5r_resampArr.Rdata')
# Run resampling (in parallel)
QULO_demoArray_Par <- 
  geo.gen.Resample.Par(genObj = QULO_genind, geoFlag = TRUE, coordPts = QULO_points, 
                       geoBuff = geo_buffSize, SDMrast=NA, boundary=world_poly_clip_W, 
                       ecoFlag = TRUE, ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, 
                       ecoLayer = 'US', reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
# Close cores
stopCluster(cl)
} else {
  # Reduce the number of samples, to allow for more efficient testing
  randomSamp <- sample(indNames(QULO_genind), size=10)
  # Subset geographic coordinate dataframe
  QULO_points_small <- QULO_points[which(QULO_points$sampleID %in% randomSamp),]
  # Subset genind object, based on smaller geo coordinates dataframe
  QULO_genind_small <- QULO_genind[QULO_points_small[,1], drop=TRUE]
  # Run resampling not in parallel (for function testing purposes)
  QULO_demoArray_IND <-
    geo.gen.Resample(genObj = QULO_genind_small, SDMrast = NA, geoFlag = FALSE, coordPts = QULO_points_small,
                     geoBuff = geo_buffSize, boundary = world_poly_clip, ecoFlag = TRUE, 
                     ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly, ecoLayer = 'US', reps = 1)
}

# # %%% ANALYZE DATA %%% ----
# # Specify filepath for QULO geographic and genetic data, including resampling array
# QULO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QULO/')
# arrayDir <- paste0(QULO_filePath, 'resamplingData/QULO_50km_GE_5r_resampArr.Rdata')
# # Read in the resampling array .Rdata object, saved to disk
# QULO_demoArray_Par <- readRDS(arrayDir)
# 
# # ---- CORRELATION ----
# # Build a data.frame from array values, to pass to linear models
# QULO_DF <- resample.array2dataframe(QULO_demoArray_Par)
# # Calculate normalized root mean square value
# QULO_nrmse_geo <- nrmse_func(obs=QULO_DF$Geo, pred=QULO_DF$Total) ; QULO_nrmse_geo
# QULO_nrmse_eco <- nrmse_func(obs=QULO_DF$Eco, pred=QULO_DF$Total) ; QULO_nrmse_eco
# 
# # ---- PLOTTING ----
# # Generate the average values (across replicates) for all proportions
# # This function has default arguments for returning just Total allelic geographic proportions
# averageValueMat <- meanArrayValues(QULO_demoArray_Par, allValues = TRUE)
# # Subset matrix of all average values to just Total allelic, geographic, and ecological coverage
# averageValueMat_TEG <- averageValueMat[,c(1,6,7)]
# # Calculate the absolute difference between genetic and geographic/ecological, and add to data.frame
# averageValueMat_TEG <- cbind(averageValueMat_TEG, abs(averageValueMat_TEG$Total-averageValueMat_TEG$Geo))
# averageValueMat_TEG <- cbind(averageValueMat_TEG, abs(averageValueMat_TEG$Total-averageValueMat_TEG$Eco))
# names(averageValueMat_TEG) <- c(names(averageValueMat_TEG)[1:3], 'Geo_Difference', 'Eco_Difference')
# # Specify plot colors
# plotColors <- c('red','red4','darkorange3','coral','darkblue', 'purple')
# plotColors <- alpha(plotColors, 0.45)
# plotColors_Sub <- plotColors[-(2:4)]
# # Two plots in a single window
# par(mfrow=c(2,1))
# 
# # ---- CORRELATION PLOTS
# plot(averageValueMat_TEG$Geo, averageValueMat_TEG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
#      main='Q. lobata: Geographic/Ecological by genetic coverage',xlab='', ylab='', col=plotColors[[5]])
# mtext(text='436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
# mtext(text='Geographic/Ecological coverage (%)', side=1, line=3, cex=1.6)
# mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# # Add points for ecological coverage
# points(x=averageValueMat$Eco, y=averageValueMat$Total, pch=20, col=plotColors[[6]])
# # Add NRMSE values for each comparison
# text(x = 76, y = 35, labels = paste0('NRMSE: ', QULO_nrmse_geo), col='darkblue', cex=0.9)
# text(x = 76, y = 20, labels = paste0('NRMSE: ', QULO_nrmse_eco), col='purple', cex=0.9)
# # Add legend
# legend(x=57, y=53, inset = 0.05, xpd=TRUE,
#        legend = c('Geographic', 'Ecological'),
#        col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.8)
# # ---- COVERAGE PLOTS
# # Use the matplot function to plot the matrix of average values, with specified settings
# matplot(averageValueMat_TEG[,1:3], ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# # Add title and x-axis labels to the graph
# title(main='Q. lobata: Geo-Eco-Gen Coverage', line=1.5)
# mtext(text='436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
# mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
# mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# # Add legend
# legend(x=275, y=60, inset = 0.05,
#        legend = c('Genetic coverage', 'Geographic coverage (50 km buffer)',
#                   'Ecological coverage (50 km buffer, EPA Level IV)'),
#        col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n',
#        y.intersp = 0.8)
# 
# # %%%% SDM AND TOTAL BUFFER COMPARISON ----
# # Specify filepath for QULO geographic and genetic data, including resampling array
# QULO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QULO/')
# arrayDir <- paste0(QULO_filePath, 'resamplingData/QULO_50km_G2E_LozaMX_5r_resampArr.Rdata')
# # Read in array and build a data.frame of values
# QULO_geoComp_50km_array <- readRDS(arrayDir)
# 
# # ---- CORRELATION ----
# # Build a data.frame from array values, to pass to linear models
# QULO_geoComp_50km_DF <- resample.array2dataframe(QULO_geoComp_50km_array)
# # Calculate normalized root mean square value
# QULO_nrmse_geo_totalBuff <-
#   nrmse_func(obs=QULO_geoComp_50km_DF$Geo_Buff, pred=QULO_geoComp_50km_DF$Total) ; QULO_nrmse_geo_totalBuff
# QULO_nrmse_geo_SDM <-
#   nrmse_func(obs=QULO_geoComp_50km_DF$Geo_SDM, pred=QULO_geoComp_50km_DF$Total) ; QULO_nrmse_geo_SDM
# 
# # ---- PLOTTING
# # Generate the average values (across replicates) for all proportions
# # This function has default arguments for returning just Total allelic and geographic proportions
# QULO_geoComp_50km_averageValueMat <- meanArrayValues(QULO_geoComp_50km_array)
# # Calculate the absolute difference between genetic and geographic approaches, and add to data.frame
# QULO_geoComp_50km_averageValueMat <-
#   cbind(QULO_geoComp_50km_averageValueMat, abs(QULO_geoComp_50km_averageValueMat$Total-QULO_geoComp_50km_averageValueMat$Geo_Buff))
# QULO_geoComp_50km_averageValueMat <-
#   cbind(QULO_geoComp_50km_averageValueMat, abs(QULO_geoComp_50km_averageValueMat$Total-QULO_geoComp_50km_averageValueMat$Geo_SDM))
# names(QULO_geoComp_50km_averageValueMat) <-
#   c(names(QULO_geoComp_50km_averageValueMat)[1:4],'Geo_Buff_Difference', 'Geo_SDM_Difference')
# # Specify plot colors
# plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue')
# plotColors_Fade <- alpha(plotColors, 0.45)
# 
# # Set plotting window to stack 2 graphs vertically
# par(mfcol=c(2,1))
# 
# # ---- CORRELATION PLOT
# # Plot genetic coverage against geographic coverage using total buffer approach
# plot(QULO_geoComp_50km_averageValueMat$Geo_Buff, QULO_geoComp_50km_averageValueMat$Total, pch=20,
#      main='Q. lobata: Geographic by genetic coverage',
#      xlab='Geographic coverage (%)', ylab='Genetic coverage (%)',
#      col=plotColors_Fade[2])
# # Add points for SDM approach
# points(x=QULO_geoComp_50km_averageValueMat$Geo_SDM, y=QULO_geoComp_50km_averageValueMat$Total,
#        pch=20, col=plotColors_Fade[3])
# # Subtitle
# mtext(text='SDM: Loza (MaxEnt); 436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.1)
# # Add NRMSE values for each comparison
# text(x = 80, y = 65, labels = paste0('NRMSE: ', QULO_nrmse_geo_totalBuff), col='red4', cex=0.9)
# text(x = 80, y = 58, labels = paste0('NRMSE: ', QULO_nrmse_geo_SDM), col='darkorange3', cex=0.9)
# # Add legend
# legend(x=58, y=74, inset = 0.05, xpd=TRUE,
#        legend = c('Total buffer approach', 'SDM approach'),
#        col=plotColors[2:3], pch = c(20,20), cex=0.9, pt.cex = 2, bty='n',
#        y.intersp = 0.8)
# # ---- COVERAGE PLOT
# # Use the matplot function to plot the matrix of average values, with specified settings
# matplot(QULO_geoComp_50km_averageValueMat[,1:3], ylim=c(0,100), col=plotColors_Fade,
#         pch=16, ylab='Coverage (%)')
# # Add title and x-axis labels to the graph
# title(main='Q. lobata: Coverage Values by Sample Size', line=1.5)
# mtext(text='SDM: Loza (MaxEnt); 436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3)
# mtext(text='Number of individuals', side=1, line=2.4)
# # Add legend
# legend(x=300, y=63, inset = 0.05, xpd=TRUE,
#        legend = c('Genetic coverage', 'Geographic, Total buffer (50 km)', 'Geographic, SDM (50 km)'),
#        col=plotColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.8)
# # # ---- MAP OF SDM AND SAMPLED POINTS
# # makeAMap(QULO_points, raster = QULO_sdm, buffer = geo_buffSize)
# 
# # %%%% SMBO: MULTIPLE BUFFER SIZES ----
# # Specify filepath for QULO geographic and genetic data, including resampling array
# QULO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QULO/')
# arrayDir <- paste0(QULO_filePath, 'resamplingData/QULO_SMBO2_G2E_5r_resampArr.Rdata')
# # Read in array and build a data.frame of values
# QULO_MultBuff_array <- readRDS(arrayDir)
# 
# # ---- CALCULATIONS ----
# # Build a data.frame from array values, to pass to linear models
# QULO_MultBuff_DF <- resample.array2dataframe(QULO_MultBuff_array)
# # Loop through the data.frame columns. The first two columns are skipped, as they're sampleNumber and the
# # predictve variable (genetic coverages)
# for(i in 3:ncol(QULO_MultBuff_DF)){
#   # Calculate NRMSE for the current column in the data.frame
#   QULO_NRMSEvalue <- nrmse.func(QULO_MultBuff_DF[,i], pred = QULO_MultBuff_DF$Total)
#   # Print result, for each explanatory variable in data.frame
#   print(paste0(names(QULO_MultBuff_DF)[[i]], ': ', QULO_NRMSEvalue))
# }
