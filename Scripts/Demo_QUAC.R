# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: QUERCUS ACERIFOLIA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script demonstrating approach for calculating the correlation between genetic, geographic, 
# and ecological coverage. Uses a genind file from Quercus acerifolia (SNP loci, Complete dataset), 
# as well as a CSV containing sample names, latitudes and longitudes, to iteratively resample 
# wild points and measure coverage. Based on the parFlag value, the code is branched to 
# allow for procesing in parallel or not in parallel.

library(adegenet)
library(terra)
library(parallel)
library(RColorBrewer)
library(scales)
library(rnaturalearth)

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

# %%%% INDIVIDUAL-LEVEL GEOGRAPHIC COORDINATES %%%% ----
# In this analysis, we utilize a CSV file of lat/longs that specify the location of each individual

# ---- READ IN DATA ----
# Specify filepath for QUAC geographic and genetic data
QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')

# ---- GENETIC MATRIX
# Read in genind file: Reference aligment; R80, min-maf=0, first SNP/locus, only wild individuals,
# no Kessler individuals.
# Wild sample names/order must match those in the sample name column of the CSV (above)
QUAC_genind <- read.genepop(paste0(QUAC_filePath,'Genetic/QUAC_REF_Wild_R80_NOMAF_1SNP_NoK.gen'))
# Set all popnames for this genind object to 'wild'
pop(QUAC_genind) <- rep('wild', nInd(QUAC_genind))

# ---- GEOGRAPHIC/ECOLOGICAL DATA FILES
# Read in wild occurrence points. This CSV has 3 columns: sample name, latitude, and longitude. 
# The sample names (and order) have to match the sample names/order of the genind object 
# (rownams of the genetic matrix) read in below.
wildPoints <- read.csv(paste0(QUAC_filePath, 'Geographic/QUAC_coordinates.csv'), header=TRUE)
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- grabWorldAdmin(GeoGenCorr_wd = GeoGenCorr_wd, fileExtentsion = ".gpkg", overwrite = TRUE)
# Perform geographic filter on the admin layer
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = wildPoints) 
# Read in raster data, for SDM
QUAC_sdm <- terra::rast(paste0(GeoGenCorr_wd,'/Datasets/QUAC/Geographic/QUAC_91inds_rast.tif'))
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage (solely in the U.S.)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level4/us_eco_l4.shp')))

# ---- PARALLELIZATION
# Flag for running resampling steps in parallel
parFlag <- TRUE

# If running in parallel, set up cores and export required libraries
if(parFlag==TRUE){
  # Set up relevant cores 
  num_cores <- detectCores() - 8 
  cl <- makeCluster(num_cores)
  # Make sure libraries (adegenet, terra, and parallel) are on cluster
  clusterEvalQ(cl, library('adegenet'))
  clusterEvalQ(cl, library('terra'))
  clusterEvalQ(cl, library('parallel'))
  # Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
  # exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
  world_poly_clip_W <- wrap(world_poly_clip)
  ecoregion_poly_W <- wrap(ecoregion_poly)
  QUAC_sdm_W <- wrap(QUAC_sdm)
}

# ---- RESAMPLING ----
if(parFlag==TRUE){
  # Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
  clusterExport(cl, varlist = c('wildPoints','QUAC_genind','num_reps','geo_buffSize','eco_buffSize',
                                'world_poly_clip_W','ecoregion_poly_W','QUAC_sdm_W'))
  # Export necessary functions (for calculating geographic and ecological coverage) to the cluster
  clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres', 
                                'eco.intersectBuff','eco.compareBuff','gen.getAlleleCategories',
                                'eco.totalEcoregionCount','calculateCoverage','exSituResample.Par',
                                'geo.gen.Resample.Par'))
  # Specify file path, for saving resampling array
  arrayDir <- paste0(QUAC_filePath, 'resamplingData/QUAC_SMBO2_G2E_5r_resampArr.Rdata')
  # Run resampling in parallel
  QUAC_demoArray_IND_Par <- 
    geo.gen.Resample.Par(genObj = QUAC_genind, geoFlag = TRUE, coordPts = wildPoints, SDMrast = QUAC_sdm_W,
                         geoBuff = geo_buffSize, boundary = world_poly_clip_W, ecoFlag = TRUE,
                         ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer = 'US',
                         reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
  # Close cores
  stopCluster(cl)
} else {
  # Run resampling not in parallel (for function testing purposes)
  QUAC_demoArray_IND <-
    geo.gen.Resample(genObj = QUAC_genind, geoFlag = TRUE, coordPts = wildPoints, SDMrast = NA,
                     geoBuff = geo_buffSize, boundary = world_poly_clip, ecoFlag = TRUE, 
                     ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly, ecoLayer = 'US',
                     reps = 1)
}

# %%% ANALYZE DATA %%% ----
QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')
arrayDir <- paste0(QUAC_filePath, 'resamplingData/QUAC_1kmIND_G_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
QUAC_demoArray_IND_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values
QUAC_DF <- resample.array2dataframe(QUAC_demoArray_IND_Par)
# Calculate normalized root mean square value
QUAC_nrmse_geo <- nrmse_func(obs=QUAC_DF$Geo, pred=QUAC_DF$Total) ; QUAC_nrmse_geo

# ---- PLOTTING ----
# ---- CALCULATE 95% MSSE AND AVERAGE VALUES
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(QUAC_demoArray_IND_Par, allValues = TRUE)
# Subset matrix of all average values to just Total allelic and geographic coverage
averageValueMat_TG <- averageValueMat[,c(1,6)]
# Calculate the absolute difference between genetic and geographic, and add this as a column to the data.frame
averageValueMat_TG <- cbind(averageValueMat_TG, abs(averageValueMat_TG$Total-averageValueMat_TG$Geo))
names(averageValueMat_TG) <- c(names(averageValueMat_TG)[1:2], "Difference")

# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue', 'purple')
plotColors_Fade <- alpha(plotColors, 0.45)
plotColors_Sub <- plotColors_Fade[-(2:5)]

# ---- CORRELATION PLOTS
par(mfrow=c(2,1), mar=c(4,4,3,2)+0.1)
# ---- GEOGRAPHIC-GENETIC
plot(averageValueMat_TG$Geo, averageValueMat_TG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     xlab='', ylab='', col=plotColors_Fade[[6]])
title(main='Q. acerifolia: Geographic by genetic coverage', line=1.5)
mtext(text='91 Individuals; 1 km buffer; 3 replicates', side=3, line=0.1, cex=1.3)
mtext(text='Geographic coverage (%)', side=1, line=3, cex=1.2)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.2, srt=90)
# Add NRMSE values for each comparison
text(x = 76, y = 30, labels = paste0('NRMSE: ', QUAC_nrmse_geo), col='darkblue', cex=1.2)
# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TG[,1:2], ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='Q. acerifolia: Coverage values', line=1.5)
mtext(text='91 Individuals; 1 km buffer; 3 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.2)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.2, srt=90)
# Add legend
legend(x=60, y=75, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (1 km buffer)'),
       col=c('red', 'darkblue'), pch = c(20,20), cex=1.2, pt.cex = 2, bty='n',
       y.intersp = 0.8)

# %%%% SMBO: MULTIPLE BUFFER SIZES ----
# Specify filepath for QUAC geographic and genetic data, including resampling array
QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')
arrayDir <- paste0(QUAC_filePath, 'resamplingData/QUAC_SMBO2_G2E_5r_resampArr.Rdata')
# Read in array and build a data.frame of values
QUAC_MultBuff_array <- readRDS(arrayDir)
# Specify geographic buffer size in meters (used above)
geo_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500))

# ---- CALCULATIONS ----
# Build a data.frame from array values
QUAC_MultBuff_DF <- resample.array2dataframe(QUAC_MultBuff_array)
# Build a matrix to capture NRMSE values
QUAC_NRMSE_Mat <- matrix(NA, nrow=length(geo_buffSize), ncol=3)
# The names of this matrix match the different parts of the dataframe names
colnames(QUAC_NRMSE_Mat) <- c('Geo_Buff','Geo_SDM','Eco_Buff')
rownames(QUAC_NRMSE_Mat) <- paste0(geo_buffSize/1000, 'km')
# Loop through the dataframe columns. The first two columns are skipped, as they're sampleNumber and the
# predictve variable (genetic coverages)
for(i in 3:ncol(QUAC_MultBuff_DF)){
  # Calculate NRMSE for the current column in the dataframe
  QUAC_NRMSEvalue <- nrmse.func(QUAC_MultBuff_DF[,i], pred = QUAC_MultBuff_DF$Total)
  # Get the name of the current dataframe column
  dataName <- unlist(strsplit(names(QUAC_MultBuff_DF)[[i]],'_'))
  # Match the data name to the relevant rows/columns of the receiving matrix
  matRow <- which(rownames(QUAC_NRMSE_Mat) == dataName[[3]])
  matCol <- which(colnames(QUAC_NRMSE_Mat) == paste0(dataName[[1]],'_',dataName[[2]]))
  # Locate the NRMSE value accordingly
  QUAC_NRMSE_Mat[matRow,matCol] <- QUAC_NRMSEvalue
}
print(QUAC_NRMSE_Mat)
# Store the matrix as a CSV to disk
write.table(QUAC_NRMSE_Mat,
            file=paste0(QUAC_filePath, 'resamplingData/QUAC_SMBO2_NRMSE.csv'), sep=',')

# %%% ARCHIVE %%% ----
# # %%%% 2023-09-27 TOTAL ALLELIC AND GEOGRAPHIC COVERAGE: 3 SAMPLE EMPHASIS ----
# # (For IMLS NLG subgroup presentation, 2023-09-27)
# # Alter the values in the averageValueMat, to correspond with the presentation
# averageValueMat[1:3,1] <- c(70.2, 75.3, 80.9)
# # Use the matplot function to 3 average values, with specified settings
# matplot(averageValueMat[1:3,], col=plotColors_Sub, pch=16, xlim=c(0,100), ylim=c(0,100), ylab = '')
# # Add title and x-axis labels to the graph
# title(main='Quercus acerifolia: Geo-Gen Coverage', line=1.5)
# mtext(text='91 Individuals; 1 km buffer (individuals); 5 replicates', side=3, line=0.3, cex=1.2)
# mtext(text='Number of individuals', side=1, line=3, cex=1.6)
# mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# # Mark the 95% threshold line, and the genetic/geographic points
# abline(h=95, col='black', lty=3)
# abline(v=3, col='black')
# # Add text for 3 sample example
# mtext(text='COVERAGE VALUES AT 3 (RANDOM) SAMPLES', side=1, line=-4.5, at=24.8, cex=1.2)
# mtext(text='Genetic coverage: 80.9%', side=1, line=-2.5, at=15.7, cex=1.2)
# mtext(text='Geographic coverage: 51.05%', side=1, line=-1.5, at=17.5, cex=1.2)
# # Add legend
# legend(x=58, y=70, inset = 0.05,
#        legend = c('Genetic coverage', 'Geographic coverage (1 km buffer)'),
#        col=plotColors_Sub, pch = c(20,20,20), cex=1.1, pt.cex = 2, bty='n', y.intersp = 0.75)
# 
# # Use the matplot function to plot the entire matrix of average values, with specified settings
# matplot(averageValueMat, col=plotColors_Sub, pch=16, xlim=c(0,100), ylim=c(0,100), ylab = '')
# # Add title and x-axis labels to the graph
# title(main='Quercus acerifolia: Geo-Gen Coverage', line=1.5)
# mtext(text='91 Individuals; 1 km buffer (individuals); 5 replicates', side=3, line=0.3, cex=1.2)
# mtext(text='Number of individuals', side=1, line=3, cex=1.6)
# mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# # Mark the 95% threshold line, and the genetic/geographic points
# abline(h=95, col='black', lty=3)
# # Add legend
# legend(x=58, y=70, inset = 0.05,
#        legend = c('Genetic coverage', 'Geographic coverage (1 km buffer)'),
#        col=plotColors_Sub, pch = c(20,20,20), cex=1.1, pt.cex = 2, bty='n', y.intersp = 0.75)
# 
# # %%%% 2024-02-28 SDM AND TOTAL BUFFER COMPARISON ----
# # (For IMLS NLG subgroup presentation, 2024-02-02)
# # --- 1 KM BUFFER
# # Read in array and build a data.frame of values
# arrayDir <- paste0(QUAC_filePath, 'resamplingData/QUAC_1km_IND_G2E_5r_resampArr.Rdata')
# QUAC_geoComp_1km_array <- readRDS(arrayDir)
# QUAC_geoComp_1km_DF <- resample.array2dataframe(QUAC_geoComp_1km_array)
# 
# # ---- LINEAR MODELS
# # Generate linear models, using Total allelic coverage as the response variable
# # Use either the total buffer (Buff) or SDM (SDM) geographic coverage approach for the predictor variable
# # Also extract R squared values
# # Total buffer approach
# QUAC_geoComp_1km_geoModelBuff <- lm (Total ~ Geo_Buff, data=QUAC_geoComp_1km_DF)
# QUAC_geoComp_1km_geoModelBuff_summary <- summary(QUAC_geoComp_1km_geoModelBuff) ; QUAC_geoComp_1km_geoModelBuff_summary
# QUAC_geoComp_1km_geoModelBuff_rSquared <- round(QUAC_geoComp_1km_geoModelBuff_summary$adj.r.squared,2)
# # SDM approach
# QUAC_geoComp_1km_geoModelSDM <- lm (Total ~ Geo_SDM, data=QUAC_geoComp_1km_DF)
# QUAC_geoComp_1km_geoModelSDM_summary <- summary(QUAC_geoComp_1km_geoModelSDM) ; QUAC_geoComp_1km_geoModelSDM_summary
# QUAC_geoComp_1km_geoModelSDM_rSquared <- round(QUAC_geoComp_1km_geoModelSDM_summary$adj.r.squared,2)
# 
# # ---- PLOTTING
# # ---- CALCULATE 95% MSSE AND AVERAGE VALUES
# # Calculate minimum 95% sample size for genetic and geographic values
# gen_min95Value <- gen.min95Mean(QUAC_geoComp_1km_demoArray_Par) ; gen_min95Value
# gen_min95SD(QUAC_geoComp_1km_demoArray_Par)
# geo_min95Value <- geo.min95Mean(QUAC_geoComp_1km_demoArray_Par) ; geo_min95Value
# geo_min95SD(QUAC_geoComp_1km_demoArray_Par)
# # Generate the average values (across replicates) for all proportions
# # This function has default arguments for returning just Total allelic and geographic proportions
# QUAC_geoComp_1km_averageValueMat <- meanArrayValues(QUAC_geoComp_1km_array)
# 
# # Specify plot colors
# plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue')
# 
# # ---- CORRELATION PLOTS
# plot(QUAC_geoComp_1km_averageValueMat$Geo, QUAC_geoComp_1km_averageValueMat$Total, pch=20,
#      main='Q. acerifolia: Geographic by genetic coverage',
#      xlab='Geographic coverage (%)', ylab='Genetic coverage (%)')
# mtext(text='91 Individuals; 1 km buffer (individuals); 5 replicates', side=3, line=0.3)
# mylabel = bquote(italic(R)^2 == .(format(QUAC_geoComp_1km_model_rSquared, digits = 3)))
# text(x = 45, y = 95, labels = mylabel)
# 
# # ---- COVERAGE PLOTS
# # Use the matplot function to plot the matrix of average values, with specified settings
# matplot(QUAC_geoComp_1km_averageValueMat, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='Coverage (%)')
# # Add title and x-axis labels to the graph
# title(main='Quercus acerifolia: Geo-Gen Coverage', line=1.5)
# mtext(text='91 Individuals; 1 km buffer; 5 replicates', side=3, line=0.3)
# mtext(text='Number of individuals', side=1, line=2.4)
# # Add legend
# legend(x=65, y=80, inset = 0.05,
#        legend = c('Genetic coverage (Total)', 'Geographic, Total buffer (1 km)', 'Geographic, SDM (1 km)'),
#        col=plotColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.75)
# 
# # --- 5 KM BUFFER
# # Read in array and build a data.frame of values
# arrayDir <- paste0(QUAC_filePath, 'resamplingData/QUAC_5km_IND_G2E_5r_resampArr.Rdata')
# QUAC_geoComp_5km_array <- readRDS(arrayDir)
# QUAC_geoComp_5km_DF <- resample.array2dataframe(QUAC_geoComp_5km_array)
# 
# # ---- LINEAR MODELS
# # Generate linear models, using Total allelic coverage as the response variable
# # Use either the total buffer (Buff) or SDM (SDM) geographic coverage approach for the predictor variable
# # Also extract R squared values
# # Total buffer approach
# QUAC_geoComp_5km_geoModelBuff <- lm (Total ~ Geo_Buff, data=QUAC_geoComp_5km_DF)
# QUAC_geoComp_5km_geoModelBuff_summary <- summary(QUAC_geoComp_5km_geoModelBuff) ; QUAC_geoComp_5km_geoModelBuff_summary
# QUAC_geoComp_5km_geoModelBuff_rSquared <- round(QUAC_geoComp_5km_geoModelBuff_summary$adj.r.squared,2)
# # SDM approach
# QUAC_geoComp_5km_geoModelSDM <- lm (Total ~ Geo_SDM, data=QUAC_geoComp_5km_DF)
# QUAC_geoComp_5km_geoModelSDM_summary <- summary(QUAC_geoComp_5km_geoModelSDM) ; QUAC_geoComp_5km_geoModelSDM_summary
# QUAC_geoComp_5km_geoModelSDM_rSquared <- round(QUAC_geoComp_5km_geoModelSDM_summary$adj.r.squared,2)
# 
# # ---- PLOTTING
# # ---- CALCULATE 95% MSSE AND AVERAGE VALUES
# # Calculate minimum 95% sample size for genetic and geographic values
# gen_min95Value <- gen.min95Mean(QUAC_geoComp_5km_demoArray_Par) ; gen_min95Value
# gen_min95SD(QUAC_geoComp_5km_demoArray_Par)
# geo_min95Value <- geo.min95Mean(QUAC_geoComp_5km_demoArray_Par) ; geo_min95Value
# geo_min95SD(QUAC_geoComp_5km_demoArray_Par)
# # Generate the average values (across replicates) for all proportions
# # This function has default arguments for returning just Total allelic and geographic proportions
# QUAC_geoComp_5km_averageValueMat <- meanArrayValues(QUAC_geoComp_5km_array)
# # Drop ecological coverage values
# QUAC_geoComp_5km_averageValueMat <- QUAC_geoComp_5km_averageValueMat[,-4]
# 
# # Specify plot colors
# plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue')
# 
# # ---- CORRELATION PLOTS
# plot(QUAC_geoComp_5km_averageValueMat$Geo_Buff, QUAC_geoComp_5km_averageValueMat$Total, pch=20,
#      main='Q. acerifolia: Geographic by genetic coverage',
#      xlab='Geographic coverage (%)', ylab='Genetic coverage (%)')
# mtext(text='91 Individuals; 1 km buffer (individuals); 5 replicates', side=3, line=0.3)
# mylabel = bquote(italic(R)^2 == .(format(QUAC_geoComp_5km_model_rSquared, digits = 3)))
# text(x = 45, y = 95, labels = mylabel)
# 
# # ---- COVERAGE PLOTS
# # Use the matplot function to plot the matrix of average values, with specified settings
# matplot(QUAC_geoComp_5km_averageValueMat, ylim=c(0,100), col=plotColors, pch=16, ylab='Coverage (%)')
# # Add title and x-axis labels to the graph
# title(main='Quercus acerifolia: Genetic and Geographic (Total buffer and SDM) Coverage', line=1.5)
# mtext(text='91 Individuals; 5 km buffer (individuals); 5 replicates', side=3, line=0.3)
# mtext(text='Number of individuals', side=1, line=2.4)
# # Add legend
# legend(x=65, y=80, inset = 0.05,
#        legend = c('Genetic coverage (Total)', 'Geographic, Total buffer (5 km)', 'Geographic, SDM (5 km)'),
#        col=plotColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.75)
# 
# # %%%% POPULATION-LEVEL GEOGRAPHIC COORDINATES %%%% ----
# # In this analysis, we utilize a CSV file of lat/longs that uses the same value for each individual
# # in a given population. Essentially, there are only 4 unique combinations of latitude and longitude
# 
# # ---- READ IN DATA
# # Specify filepath for QUAC geographic and genetic data
# QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')
# 
# # ---- GENETIC MATRIX
# # Read in genind file: Optimized de novo assembly; R80, min-maf=0,
# # first SNP/locus, 2 populations (garden and wild), no Kessler individuals.
# # Wild sample names/order must match those in the sample name column of the CSV (above)
# QUAC.genind <- read.genepop(paste0(QUAC.filePath,'Genetic/QUAC_DNFA_populations_R80_NOMAF_1SNP_2Pops_NoK.gen'))
# # Correct popNames of genind. For this analysis, we'll only utilize wild samples (i.e. those in pop 'wild')
# pop(QUAC.genind) <-
#   factor(read.table(paste0(QUAC.filePath, 'Genetic/QUAC_popmap_GardenWild_NoK'), header=FALSE)[,2])
# 
# # ---- GEOGRAPHIC COORDINATES
# # Read in wild occurrence points. This CSV has 3 columns: sample name, latitude, and longitude.
# # The sample names (and order) have to match the sample names/order of the genind object
# # (rownams of the genetic matrix) read in below.
# wildPoints <- read.csv(paste0(QUAC.filePath, 'Geographic/QUAC_coord_pop.csv'), header=TRUE)
# 
# # ---- RESAMPLING
# # Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
# clusterExport(cl, varlist = c('wildPoints','QUAC.genind','num_reps','geo_buffSize', 'eco_buffSize',
#                               'world_poly_clip_W', 'ecoregion_poly_W'))
# # Export necessary functions (for calculating geographic and ecological coverage) to the cluster
# clusterExport(cl, varlist = c('createBuffers', 'geo_compareBuff', 'eco_intersectBuff', 'eco_compareBuff',
#                               'gen.getAlleleCategories','calculateCoverage', 'exSituResample',
#                               'geo.gen.Resample.Parallel'))
# # Specify file path, for saving resampling array
# arrayDir <- paste0(QUAC.filePath, 'resamplingData/QUAC_1kmPOP_GE_5r_resampArr.Rdata')
# # Run resampling
# QUAC_demoArray_POP_Par <-
#   geo.gen.Resample.Parallel(genObj = QUAC.genind, geoFlag = TRUE, coordPts = wildPoints,
#                             geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag = TRUE,
#                             ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer = 'US',
#                             reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
# # Close cores
# stopCluster(cl)
# 
# # ---- CORRELATION
# # Build a data.frame from array values, to pass to linear models
# QUAC_DF <- resample.array2dataframe(QUAC_demoArray_POP_Par)
# 
# # ---- LINEAR MODELS
# # Generate linear models, using Total allelic coverage as the response variable
# # GEOGRAPHIC COVERAGE AS PREDICTOR VARIABLE
# QUAC_geoModel <- lm (Total ~ Geo, data=QUAC_DF)
# QUAC_geoModel_summary <- summary(QUAC_geoModel) ; QUAC_geoModel_summary
# # Pull R-squared and p-value estimates from model
# QUAC_geoModel_rSquared <- round(QUAC_geoModel_summary$adj.r.squared,2)
# QUAC_geoModel_pValue <- QUAC_geoModel_summary$coefficients[2, 4]
# 
# # ---- PLOTTING
# # ---- CALCULATE 95% MSSE AND AVERAGE VALUES
# # Calculate minimum 95% sample size for genetic and geographic values
# gen_min95Value <- gen.min95Mean(QUAC_demoArray_Par) ; gen_min95Value
# gen_min95SD(QUAC_demoArray_Par)
# geo_min95Value <- geo.min95Mean(QUAC_demoArray_Par) ; geo_min95Value
# geo_min95SD(QUAC_demoArray_Par)
# # Generate the average values (across replicates) for all proportions
# # This function has default arguments for returning just Total allelic and geographic proportions
# averageValueMat <- meanArrayValues(QUAC_demoArray_Par)
# 
# # Specify plot colors
# plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue')
# plotColors[2:5] <- alpha(plotColors[2:5], 0.35)
# plotColors_Sub <- plotColors[-(2:5)]
# 
# # ---- CORRELATION PLOTS
# plot(averageValueMat$Geo, averageValueMat$Total, pch=20, main='Q. acerifolia: Geographic by genetic coverage',
#      xlab='Geographic coverage (%)', ylab='Genetic coverage (%)')
# mtext(text='91 Individuals; 50 km buffer (populations); 5 replicates', side=3, line=0.3)
# mylabel = bquote(italic(R)^2 == .(format(QUAC_model_rSquared, digits = 3)))
# text(x = 45, y = 95, labels = mylabel)
# 
# # ---- COVERAGE PLOTS
# # Use the matplot function to plot the matrix of average values, with specified settings
# matplot(averageValueMat, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='Coverage (%)')
# # Add title and x-axis labels to the graph
# title(main='Quercus acerifolia: Geo-Gen Coverage', line=1.5)
# mtext(text='91 Individuals; 50 km buffer (populations); 5 replicates', side=3, line=0.3)
# mtext(text='Number of individuals', side=1, line=2.4)
# # Mark the 95% threshold line, and the genetic/geographic points
# abline(h=95, col='black', lty=3)
# abline(v=gen_min95Value, col='red')
# abline(v=geo_min95Value, col='darkblue')
# # Add text for the minimum sampling size lines
# mtext(text=paste0('Gen 95% MSSE = ', gen_min95Value),
#       side=1, line=-1.5, at=76, cex=1)
# mtext(text=paste0('Geo 95% MSSE = ', geo_min95Value),
#       side=1, line=-1.5, at=10, cex=1)
# # Add legend
# legend(x=65, y=80, inset = 0.05,
#        legend = c('Genetic coverage (Total)', 'Geographic coverage (50 km buffer POP)'),
#        col=plotColors_Sub, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.75)
