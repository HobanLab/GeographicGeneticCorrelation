# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: PINUS CONTORTA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage
# for Pinus contorta. Uses data files from MacLachlan et al. 2021 to pull in genetic data (as a STRUCTURE input file) 
# and geographic coordinates (included in a CSV) to conduct correlation analyses. We also process
# the geographic coordinates to match the order of the genetic samples, prior to passing both along
# to the correlation analyses.

library(adegenet)
library(terra)
library(parallel)
library(RColorBrewer)
library(scales)

# Read in relevant functions
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/Code/'
# option2 
GeoGenCorr_wd <- '~/Documents/GeographicGeneticCorrelation/' 
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
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.gpkg')))
# Read in the EPA Level III ecoregion shapefile, which is used for calculating ecological coverage (in North America)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level3/NA_CEC_Eco_Level3.shp')))
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
# Specify filepath for PICO geographic and genetic data. Because the STRUCTURe input file is greater than 50 MB 
# in size, it is not included on the GitHub repository
PICO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/PICO/')

# ---- GENETIC MATRIX
# The 1st script in the MacLachlan et al. 2021 supplement (1_MacLachlan_etal_Pine_GPA_ped&mapfile_formatting_Jan10th2021.R)
# generates a .ped and .map file. This is for 929 individuals (control individuals are not included), and it includes
# 32,449 loci, after filtering for minor alleles and missing data. The .ped file was then are passed into PLINK 
# in order to generate a STRUCTURE input file, which is converted into a genind object in the call below.
PICO_genind <- read.structure(file=paste0(PICO_filePath, 'Genetic/Pine_NaturalComp_85SNPfilter.stru'), 
                              n.ind = 929, n.loc = 32449, onerowperind = TRUE, col.lab = 1, 
                              col.pop = 0, row.marknames = 0, sep = ' ', ask = FALSE)
# Capture the names of samples as they're ordered in the genind object. This is for the steps below.
PICO_sampleNames <- indNames(PICO_genind)

# ---- GEOGRAPHIC COORDINATES
# The supplement of MacLachlan et al. 2021 includes a csv that contains climate data for all of the seedlings in the study,
# as well as coordinate information. This dataset needs to first be subset just to the 929 samples included in the
# genind object (above), and then ordered to match the order of samples in that genind object. These steps are taken below.
PICO_coordinates <- 
  read.csv(file=paste0(PICO_filePath, 'Geographic/Pine_NaturalOrchard_ClimateData_Sept7th2015.csv'), header = TRUE)
# Start by subsetting the CSV to just the variables we need: sample names, latitude, and longitude
PICO_coordinates <- 
  PICO_coordinates[which(PICO_coordinates$Internal_ID %in% PICO_sampleNames),c('Internal_ID','Latitude','Longitude')]
# Reorder the coordinate values to match the order of samples in the genind file
PICO_coordinates <- PICO_coordinates[order(match(PICO_coordinates$Internal_ID, PICO_sampleNames)),]
# Rename the columns of the geographic coordinates data.frame (because geo.compareBuff function expects certain strings)
colnames(PICO_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('PICO_coordinates','PICO_genind','num_reps','geo_buffSize', 'eco_buffSize',
                              'world_poly_clip_W', 'ecoregion_poly_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers', 'geo.compareBuff', 'eco.intersectBuff', 'eco.compareBuff',
                              'gen.getAlleleCategories','calculateCoverage', 'exSituResample.Par', 
                              'geo.gen.Resample.Par'))
# Specify file path, for saving resampling array
arrayDir <- paste0(PICO_filePath, 'resamplingData/PICO_50km_GE_5r_resampArr.Rdata')

# Run resampling (in parallel)
print("%%% BEGINNING RESAMPLING %%%")
PICO_demoArray_Par <- 
  geo.gen.Resample.Parallel(gen_obj = PICO_genind, geoFlag = TRUE, coordPts = PICO_coordinates, 
                            geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag = TRUE, 
                            ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer = 'NA', 
                            reps = num_reps, arrayFilepath = arrayDir, cluster = cl)

# Close cores
stopCluster(cl)

# %%% ANALYZE DATA %%% ----
# Read in the resampling array .Rdata object, saved to disk
PICO_demoArray_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values, to pass to linear models
PICO_DF <- resample.array2dataframe(PICO_demoArray_Par)

# ---- LINEAR MODELS
# Generate linear models, using Total allelic coverage as the response variable
# GEOGRAPHIC COVERAGE AS PREDICTOR VARIABLE
PICO_geoModel <- lm (Total ~ Geo, data=PICO_DF)
PICO_geoModel_summary <- summary(PICO_geoModel) ; PICO_geoModel_summary
# Pull R-squared estimate from model
PICO_geoModel_rSquared <- round(PICO_geoModel_summary$adj.r.squared,2)
# ECOLOGICAL COVERAGE AS PREDICTOR VARIABLE
PICO_ecoModel <- lm (Total ~ Eco, data=PICO_DF)
PICO_ecoModel_summary <- summary(PICO_ecoModel) ; PICO_ecoModel_summary
# Pull R-squared estimate from model
PICO_ecoModel_rSquared <- round(PICO_ecoModel_summary$adj.r.squared, 2)

# ---- PLOTTING ----
# ---- CALCULATE 95% MSSE AND AVERAGE VALUES
# Calculate minimum 95% sample size for genetic and geographic values
gen_min95Value <- gen.min95Mean(PICO_demoArray_Par) ; gen_min95Value
geo_min95Value <- geo.min95Mean(PICO_demoArray_Par) ; geo_min95Value
eco_min95Value <- eco.min95Mean(PICO_demoArray_Par) ; eco_min95Value
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(PICO_demoArray_Par, allValues = TRUE)
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
     main='P. contorta: Geographic by genetic coverage',xlab='', ylab='')
mtext(text='929 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Geographic coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
mylabel = bquote(italic(R)^2 == .(format(PICO_geoModel_rSquared, digits = 3)))
text(x = 2, y = 10, labels = mylabel, cex=0.8)
# ---- ECOLOGICAL-GENETIC
plot(averageValueMat_TEG$Eco, averageValueMat_TEG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     main='P. contorta: Ecological by genetic coverage',xlab='', ylab='')
mtext(text='929 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Ecological coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
mylabel = bquote(italic(R)^2 == .(format(PICO_ecoModel_rSquared, digits = 3)))
text(x = 2, y = 10, labels = mylabel, cex=0.8)

# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TEG, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='Pinus contorta: Gen-Geo-Eco Coverage', line=1.5)
mtext(text='929 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=405, y=60, inset = 0.05,
       legend = c('Genetic coverage (Total)', 'Geographic coverage (50 km buffer)', 'Ecological coverage (EPA Level III)'),
       col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=1.2, pt.cex = 2, bty='n', 
       y.intersp = 0.5)
