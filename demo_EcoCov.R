# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% TESTING ECOLOGICAL COVERAGE %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script for testing out first draft of functions calculating ecological coverage

library(adegenet)
library(terra)
library(parallel)
library(RColorBrewer)
library(scales)

# Read in relevant functions
GeoGenCorr.wd <- "/home/akoontz/Documents/GeoGenCorr/"
setwd(GeoGenCorr.wd)
source("Code/functions_GeoGenCoverage.R")

# ---- VARIABLES ----
# Declare a directory within which to store .Rdata objects of resampling arrays
resamplingDataDir <- paste0(GeoGenCorr.wd, "Code/resamplingData/")
# Specify number of resampling replicates. 
# num_reps <- 5
num_reps <- 3
# ---- GEOGRAPHIC VARIABLES
# Specify geographic buffer size in meters 
geo_buffSize <- 1000
# geo_buffSize <- 50000
# Specify ecological buffer size in meters 
eco_buffSize <- 1000
# eco_buffSize <- 50000
# SHAPEFILES
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- 
  vect(file.path(paste0(GeoGenCorr.wd, "GIS_shpFiles/world_countries_10m/world_countries_10m.shp")))
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage (solely in the U.S.)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr.wd, "GIS_shpFiles/ecoregions_EPA_level4/us_eco_l4.shp")))
# Shapefiles are by default a "non-exportable" object, which means the must be processed before being
# exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
world_poly_clip_W <- wrap(world_poly_clip)
ecoregion_poly_W <- wrap(ecoregion_poly)

# ---- PARALLELIZATION
# Set up relevant cores 
num_cores <- detectCores() - 8 
cl <- makeCluster(num_cores)
# Make sure libraries (adegenet, terra, and parallel) are on cluster
clusterEvalQ(cl, library("adegenet"))
clusterEvalQ(cl, library("terra"))
clusterEvalQ(cl, library("parallel"))

# %%%% INDIVIDUAL-LEVEL GEOGRAPHIC COORDINATES %%%% ----
# In this analysis, we utilize a CSV file of lat/longs that specify the location of each individual
# ---- READ IN GEOGRAPHIC AND GENETIC DATA ----
# Specify filepath for QUAC geographic and genetic data
QUAC.filePath <- paste0(GeoGenCorr.wd, "Datasets/QUAC/")

# ---- GEOGRAPHIC/ECOLOGICAL
# Read in wild occurrence points. This CSV has 3 columns: sample name, latitude, and longitude. 
# The sample names (and order) have to match the sample names/order of the genind object 
# (rownams of the genetic matrix) read in below.
wildPoints <- read.csv(paste0(QUAC.filePath, "QUAC_coord_ind.csv"), header=TRUE)

# ---- GENETIC
# Read in genind file: Optimized de novo assembly; R80, min-maf=0, 
# first SNP/locus, 2 populations (garden and wild), no Kessler individuals.
# Wild sample names/order must match those in the sample name column of the CSV (above)
QUAC.genind <- read.genepop(paste0(QUAC.filePath,"QUAC_DNFA_populations_R80_NOMAF_1SNP_2Pops_NoK.gen"))
# Correct popNames of genind. For this analysis, we'll only utilize wild samples (i.e. those in pop "wild")
pop(QUAC.genind) <- 
  factor(read.table(paste0(QUAC.filePath, "QUAC_popmap_GardenWild_NoK"), header=FALSE)[,2])

# TEST: calculateCoverage
QUAC.genMat <- QUAC.genind@tab[which(pop(QUAC.genind) == "wild"),]
calculateCoverage(gen_mat = QUAC.genMat, coordPts = wildPoints, geoBuff = geo_buffSize,
                  boundary=world_poly_clip, ecoBuff = eco_buffSize, 
                  ecoRegions = ecoregion_poly, numSamples = 2)

# TEST: exSituResample
# Start time
start <- Sys.time() ; print(start)
# Function test
exSituResample(gen_obj = QUAC.genind, coordPts = wildPoints, geoBuff = geo_buffSize,
               boundary = world_poly_clip, ecoBuff = eco_buffSize, 
               ecoRegions = ecoregion_poly, parFlag = FALSE)
# End time
end <- Sys.time() ; print(end)
print(paste0("Run time: ", end-start))

# TEST: geo.gen.Resample
# Start time
start <- Sys.time() ; print(start)
# Function test
test <- geo.gen.Resample(gen_obj = QUAC.genind, geoFlag = TRUE, coordPts = wildPoints, geoBuff = 1000,
                 boundary=world_poly_clip, ecoFlag = TRUE, ecoBuff = 1000, ecoRegions = ecoregion_poly, 
                 ecoLayer = "US", reps = num_reps)
# End time
end <- Sys.time() ; print(end)
print(paste0("Run time: ", end-start))

# TEST: geo.gen.Resample.Parallel
arrayDir <- paste0(QUAC.filePath, "resamplingData/QUAC_1kmIND_3r_EcoGeoTest_resampArr.Rdata")
# Start time
start <- Sys.time() ; print(start)
# Function test
geo.gen.Resample.Parallel(gen_obj = QUAC.genind, geoFlag = TRUE, coordPts = wildPoints, geoBuff = 1000,
                          boundary=world_poly_clip_W, ecoFlag = TRUE, ecoBuff = 1000, ecoRegions = ecoregion_poly_W, 
                          ecoLayer = "US", reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
# End time
end <- Sys.time() ; print(end)
print(paste0("Run time: ", end-start))
