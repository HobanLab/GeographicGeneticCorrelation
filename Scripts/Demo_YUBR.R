# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: YUCCA BREVIFOLIA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage 
# in Yucca brevifolia (data files from Royer et al. 2016). Genetic data sourced from 
# CSV reformatted to STRUCTURE input file; geographic coordinates are converted and
# cleaned before being utilized for correlations.

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
geo_buffSize <- 50000
# Specify ecological buffer size in meters
eco_buffSize <- 50000
# ---- SHAPEFILES
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <-
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.shp')))
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage (solely in the U.S.)
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

# ---- GEOGRAPHIC COORDINATES
# Read in wild occurrence points. This CSV has 5 columns: sample names, 
# latitude (decimal minutes, 2 columns), and longitude (decimal minutes, 2 columns)  
YUBR_points <- read.csv(
  paste0(YUBR_filePath, 'Geographic/TreeCoordinates.csv'), header=TRUE)
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

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) 
# to the cluster
clusterExport(cl, varlist = c('YUBR_genind','YUBR_coordinates','num_reps',
                              'geo_buffSize', 'eco_buffSize','world_poly_clip_W', 
                              'ecoregion_poly_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers', 'geo.compareBuff', 'eco.intersectBuff', 'eco.compareBuff',
                              'gen.getAlleleCategories','calculateCoverage', 'exSituResample.Par', 
                              'geo.gen.Resample.Par'))
# Specify file path, for saving resampling array
arrayDir <- paste0(YUBR_filePath, 'resamplingData/YUBR_50km_GE_5r_resampArr.RData')
# Run resampling (in parallel)
YUBR_demoArray_Par <- 
  geo.gen.Resample.Par(gen_obj = YUBR_genind, geoFlag = TRUE, coordPts = YUBR_coordinates, 
                       geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag = TRUE, 
                       ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer = 'US', 
                       reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
# Close cores
stopCluster(cl)

# Run resampling not in parallel (for function testing purposes)
# YUBR_demoArray_IND <-
#   geo.gen.Resample(gen_obj = YUBR_genind, geoFlag = TRUE, coordPts = YUBR_coordinates, geoBuff = geo_buffSize, 
#                    boundary = world_poly_clip, ecoFlag = TRUE, ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly,
#                    ecoLayer = "US", reps = 1)

# %%% ANALYZE DATA %%% ----
# Read in the resampling array .Rdata object, saved to disk
YUBR_demoArray_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values, to pass to linear models
YUBR_DF <- resample.array2dataframe(YUBR_demoArray_Par)

# ---- LINEAR MODELS
# Generate linear models, using Total allelic coverage as the response variable
# GEOGRAPHIC COVERAGE AS PREDICTOR VARIABLE
YUBR_geoModel <- lm (Total ~ Geo, data=YUBR_DF)
YUBR_geoModel_summary <- summary(YUBR_geoModel) ; YUBR_geoModel_summary
# Pull R-squared estimate from model
YUBR_geoModel_rSquared <- round(YUBR_geoModel_summary$adj.r.squared,2)
# ECOLOGICAL COVERAGE AS PREDICTOR VARIABLE
YUBR_ecoModel <- lm (Total ~ Eco, data=YUBR_DF)
YUBR_ecoModel_summary <- summary(YUBR_ecoModel) ; YUBR_ecoModel_summary
# Pull R-squared estimate from model
YUBR_ecoModel_rSquared <- round(YUBR_ecoModel_summary$adj.r.squared, 2)

# ---- PLOTTING ----
# ---- CALCULATE 95% MSSE AND AVERAGE VALUES
# Calculate minimum 95% sample size for genetic and geographic values
gen_min95Value <- gen.min95Mean(YUBR_demoArray_Par) ; gen_min95Value
geo_min95Value <- geo.min95Mean(YUBR_demoArray_Par) ; geo_min95Value
eco_min95Value <- eco.min95Mean(YUBR_demoArray_Par) ; eco_min95Value
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(YUBR_demoArray_Par, allValues = TRUE)
# Subset matrix of all average values to just Total allelic, geographic, and ecological coverage
averageValueMat_TEG <- averageValueMat[,c(1,6,7)]

# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue', 'purple')
plotColors <- alpha(plotColors, 0.45)
plotColors_Sub <- plotColors[-(2:5)]

# ---- CORRELATION PLOTS
par(mfrow=c(2,1))
# ---- GEOGRAPHIC-GENETIC
plot(averageValueMat_TEG$Geo, averageValueMat_TEG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     main='Y. brevifolia: Geographic by genetic coverage',xlab='', ylab='')
mtext(text='319 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Geographic coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
mylabel = bquote(italic(R)^2 == .(format(YUBR_geoModel_rSquared, digits = 3)))
text(x = 2, y = 10, labels = mylabel, cex=0.8)
# ---- ECOLOGICAL-GENETIC
plot(averageValueMat_TEG$Eco, averageValueMat_TEG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     main='Y. brevifolia: Ecological by genetic coverage',xlab='', ylab='')
mtext(text='319 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Ecological coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
mylabel = bquote(italic(R)^2 == .(format(YUBR_ecoModel_rSquared, digits = 3)))
text(x = 2, y = 10, labels = mylabel, cex=0.8)

# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TEG, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='Yucca brevifolia: Geo-Eco-Gen Coverage', line=1.5)
mtext(text='319 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=50, y=60, inset = 0.05,
       legend = c('Genetic coverage (Total)', 'Geographic coverage (50 km buffer)', 'Ecological coverage (EPA Level IV)'),
       col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=1.2, pt.cex = 2, bty='n',
       y.intersp = 0.8)
