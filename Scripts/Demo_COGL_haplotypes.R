# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: CONRADINA GLABRA USING HAPLOTYPES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage
# for Conradina glabra, a dataset originally provided by Lauren Eserman. This script explores
# the impact of using haplotypes of different lengths on genetic resampling curves. To do this,
# resampling is performed multiple times using different input files, with the results of each
# haplotype length stored to a different array object.

library(adegenet)
library(terra)
library(parallel)
library(RColorBrewer)
library(scales)
library(vcfR) # for vcfR2genind function

# Read in relevant functions
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/Code/'
setwd(GeoGenCorr_wd)
source('Scripts/functions_GeoGenCoverage.R')

# ---- VARIABLES ----
# Specify number of resampling replicates
num_reps <- 5
# ---- BUFFER SIZES
# Specify geographic buffer size in meters 
geo_buffSize <- 1000
# Specify ecological buffer size in meters 
eco_buffSize <- 1000
# ---- SHAPEFILES
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.shp')))
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage 
# (solely in the U.S.)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level4/us_eco_l4.shp')))
# Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
# exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
world_poly_clip_W <- wrap(world_poly_clip)
ecoregion_poly_W <- wrap(ecoregion_poly)

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
# Specify filepath for COGL geographic and genetic data
COGL_filePath <- paste0(GeoGenCorr_wd, 'Datasets/COGL/')

# ---- GENETIC MATRICES
# For each haplotype length (ranging from a single SNP to 10 concatenated SNPs), the genetic
# data is stored as a STRUCTURE input file. Each STRUCTURE file is read in and converted into a
# genind object. 
COGL_haplos_filePath <- paste0(COGL_filePath, 'Genetic/haplotypicData')
COGL_haplos_input <- list.files(COGL_haplos_filePath, full.names = TRUE)
# Reposition the haplotype file of 10 concatenated SNPs to be the last file (instead of the 2nd)
COGL_haplos_input <- c(COGL_haplos_input[[1]],COGL_haplos_input[3:10],COGL_haplos_input[[2]])
# Declare a list of integers, storing the number of loci for each STRUCTURE input file
COGL_lociLevels <- c(16291,7045,4001,2509,1627,1083,727,473,330,212)
# Predeclare an empty list to store genind objects, and a list of resampling arrays to write to
COGL_genList <- list()
# For each STRUCTURE input file, 
for(i in 1:length(COGL_haplos_input)){
  # Read the file in and convert it into a genind object
  COGL_genList[[i]] <- read.structure(file=COGL_haplos_input[[i]], n.ind=564, 
                                    n.loc=COGL_lociLevels[[i]],
                                    onerowperind=TRUE, col.lab=3, col.pop=1,
                                    col.others=2, row.marknames=0, NA.char=-9,
                                    sep=" ", ask=TRUE, quiet=TRUE)
  # Rename the individuals in the file. This is for sample names to match the coordinate file
  indNames(COGL_genList[[i]]) <- gsub("_sorted", "", indNames(COGL_genList[[i]]))
  indNames(COGL_genList[[i]]) <- paste0("Congla", indNames(COGL_genList[[i]]))
  # Remove two duplicate samples (with "_rep" in sample name), leaving 562 samples
  COGL_genList[[i]] <- COGL_genList[[i]][-grep('_rep', indNames(COGL_genList[[i]])), drop=TRUE]
  # Specify the population of every sample to be 'wild'
  pop(COGL_genList[[i]]) <- rep('wild', nInd(COGL_genList[[i]]))
}

# ---- COORDINATE POINTS
# CSV was provided by Lauren Eserman. Because individuals were sampled very close 
# together, coordinates were determined by obtaining coordinates for the center of a
# "clump" of individuals, and then individual coordinates were calculated based on the
# separation and direction of each individual from the center. Some samples have identical
# coordinate values.
# NOTE: these coordinates cannot be shared externally, due to the rare status of this species!
COGL_coordinates <- 
  read.csv(file=paste0(COGL_filePath, 'Geographic/Conradina_coord.csv'), header = TRUE)
# Rename the columns of the geographic coordinates data.frame (because geo.compareBuff function expects certain strings)
colnames(COGL_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')
# Remove two duplicate samples (with "_rep" in sample name), leaving 562 samples
COGL_coordinates<- COGL_coordinates[-grep('_rep', COGL_coordinates$Sample.Name),]

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('COGL_coordinates','COGL_genList','num_reps','geo_buffSize', 'eco_buffSize',
                              'world_poly_clip_W', 'ecoregion_poly_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers', 'geo.compareBuff', 'eco.intersectBuff', 'eco.compareBuff',
                              'gen.getAlleleCategories','calculateCoverage', 'exSituResample.Par', 
                              'geo.gen.Resample.Par'))
# Specify file paths, for saving resampling arrays
COGL_haplos_output <- dirname(gsub("Genetic", "resamplingData", COGL_haplos_input))
for(i in 1:length(COGL_haplos_output)){
  COGL_haplos_output[[i]] <- paste0(COGL_haplos_output[[i]],'/COGL_1km_GE_5r_Hap-', i, '_resampArr.Rdata')
}

# Run resampling (in parallel). Use a loop to loop through the different haplotype lengths
print('%%% BEGINNING RESAMPLING %%%')
for(i in 1:length(COGL_haplos_output)){
  print(paste0('%%% HAPLOTYPE LENGTH: ', i))
  # Run resampling
  COGL_demoArray_Par <- 
    geo.gen.Resample.Par(gen_obj = unlist(COGL_genList[[i]]), geoFlag = TRUE, coordPts = COGL_coordinates, 
                         geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag = TRUE, 
                         ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer = 'US', 
                         reps = num_reps, arrayFilepath = COGL_haplos_output[[i]], cluster = cl)
}
# Close cores
stopCluster(cl)

# Run resampling not in parallel (for function testing purposes)
# COGL_demoArray_IND <-
#   geo.gen.Resample(gen_obj = unlist(COGL_genList[[9]]), geoFlag = TRUE, coordPts = COGL_coordinates,
#                    geoBuff = geo_buffSize, boundary = world_poly_clip, ecoFlag = FALSE, reps = 1)

# # %%% ANALYZE DATA %%% ----
# # Specify filepath for COGL geographic and genetic data, including resampling data
# COGL_filePath <- paste0(GeoGenCorr_wd, 'Datasets/COGL/')
# arrayDir <- paste0(COGL_filePath, 'resamplingData/COGL_50km_GE_5r_resampArr.Rdata')
# # Read in the resampling array .Rdata object, saved to disk
# COGL_array <- readRDS(arrayDir)
# 
# # ---- CORRELATION ----
# # Build a data.frame from array values, to pass to linear models
# COGL_DF <- resample.array2dataframe(COGL_array)
# # Calculate Spearman's r for geographic coverage
# COGL_spearR_geo <- round(cor(COGL_DF$Geo, COGL_DF$Total, method = 'spearman'),3) ; COGL_spearR_geo
# 
# # # ---- LINEAR MODELS
# # # Generate linear models, using Total allelic coverage as the response variable
# # # GEOGRAPHIC COVERAGE AS PREDICTOR VARIABLE
# # COGL_geoModel <- lm (Total ~ Geo, data=COGL_DF)
# # COGL_geoModel_summary <- summary(COGL_geoModel) ; COGL_geoModel_summary
# # # Pull R-squared estimate from model
# # COGL_geoModel_rSquared <- round(COGL_geoModel_summary$adj.r.squared,2) ; COGL_geoModel_rSquared
# # # (Ecological coverage is 100% for a single sample, so not included in  these analyses)
# 
# # ---- PLOTTING ----
# # ---- CALCULATE 95% MSSE AND AVERAGE VALUES
# # Generate the average values (across replicates) for all proportions
# # This function has default arguments for returning just Total allelic geographic proportions
# averageValueMat <- meanArrayValues(COGL_array, allValues = TRUE)
# # Subset matrix of all average values to just Total allelic and geographic coverage
# averageValueMat_TG <- averageValueMat[,c(1,6)]
# # Calculate the absolute difference between genetic and geographic, and add this as a column to the data.frame
# averageValueMat_TG <- cbind(averageValueMat_TG, abs(averageValueMat_TG$Total-averageValueMat_TG$Geo))
# names(averageValueMat_TG) <- c(names(averageValueMat_TG)[1:2], "Difference")
# 
# # Specify plot colors
# plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue', 'purple')
# plotColors <- alpha(plotColors, 0.45)
# plotColors_Sub <- plotColors[-(2:5)]
# 
# # ---- CORRELATION PLOTS
# par(mfrow=c(2,1), mar=c(4,4,3,2)+0.1)
# # ---- GEOGRAPHIC-GENETIC
# plot(averageValueMat_TG$Geo, averageValueMat_TG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
#      xlab='', ylab='', col='darkblue')
# title(main='C. glabra: Geographic by genetic coverage', line=1.5)
# mtext(text='562 Individuals; 1 km buffer; 5 replicates', side=3, line=0.1, cex=1.3)
# mtext(text='Geographic coverage (%)', side=1, line=3, cex=1.2)
# mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.2, srt=90)
# # Add Spearman's r values for each comparison
# text(x = 76, y = 43, labels = paste0('Spearman r: ', COGL_spearR_geo), col='darkblue', cex=1.2)
# # ---- COVERAGE PLOTS
# # Use the matplot function to plot the matrix of average values, with specified settings
# matplot(averageValueMat_TG, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# # Add title and x-axis labels to the graph
# title(main='C. glabra: Coverage values', line=1.5)
# mtext(text='562 Individuals; 1 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
# mtext(text='Number of individuals', side=1, line=2.4, cex=1.2)
# mtext(text='Coverage (%)', side=2, line=2.3, cex=1.2, srt=90)
# # Add legend
# legend(x=350, y=80, inset = 0.05,
#        legend = c('Genetic coverage', 'Geographic coverage (1 km buffer)'),
#        col=c('red', 'darkblue'), pch = c(20,20), cex=1.2, pt.cex = 2, bty='n',
#        y.intersp = 0.4)
# # ---- DIFFERENCE PLOTS
# # Plot difference between geographic and genetic coverage
# matplot(averageValueMat_TG[[3]], col=plotColors[[6]], pch=16, ylab='Gen-Geo Difference')
# # Add title and x-axis labels to the graph
# title(main='C. glabra: Genetic Geographic Coverage Difference', line=1.5)
# mtext(text='562 Individuals; 1 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
# mtext(text='Number of individuals', side=1, line=2.4, cex=1.2)
