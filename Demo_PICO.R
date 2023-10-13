# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEO-ECO-GEN CORRELATION DEMO: PINUS CONTORTA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script demonstrating calculating the correlation between genetic, geographic, and ecological coverage
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
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/'
setwd(GeoGenCorr_wd)
source('Code/functions_GeoGenCoverage.R')

# ---- VARIABLES ----
# Specify number of resampling replicates. 
num_reps <- 5
# ---- BUFFER SIZES
# Specify geographic buffer size in meters 
# geo_buffSize <- 1000
geo_buffSize <- 50000
# Specify ecological buffer size in meters 
# eco_buffSize <- 1000
eco_buffSize <- 50000
# ---- SHAPEFILES
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.shp')))
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

# ---- READ IN DATA ----
# Specify filepath for PICO geographic and genetic data
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

# ---- COORDINATE POINTS
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
                              'getAlleleCategories','calculateCoverage', 'exSituResample', 
                              'geo.gen.Resample.Parallel'))
# Specify file path, for saving resampling array
arrayDir <- paste0(PICO_filePath, 'resamplingData/PICO_50km_GE_5r_resampArr.Rdata')

# Run resampling (in parallel)
print("%%% BEGINNING RESAMPLING %%%")
PICO_demoArray_Par <- 
  geo.gen.Resample.Parallel(gen_obj = PICO_genind, geoFlag = TRUE, coordPts = wildPoints, 
                            geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag = TRUE, 
                            ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer = 'NA', 
                            reps = num_reps, arrayFilepath = arrayDir, cluster = cl)

# Close cores
stopCluster(cl)
