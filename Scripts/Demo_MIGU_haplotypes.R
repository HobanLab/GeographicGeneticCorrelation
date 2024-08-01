# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: MIMULUS GUTTATUS USING HAPLOTYPES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage
# for Mimulus guttatus. This script explores the impact of using haplotypes of different 
# lengths on genetic resampling curves. To do this, resampling is performed multiple times 
# using different input files, with the results of each haplotype length stored to a 
# different array object.

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
geo_buffSize <- 50000
# Specify ecological buffer size in meters 
eco_buffSize <- 50000

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
# Specify filepath for MIGU geographic and genetic data
MIGU_filePath <- paste0(GeoGenCorr_wd, 'Datasets/MIGU/')

# ---- GEOGRAPHIC/ECOLOGICAL DATA FILES
# The supplement of Vallejo-Martin et al. 2021 includes a PDF of the latitudes and longitudes for each population 
# (SupplementaryData1). This was first converted to a Google Sheet. In order to convert this to a CSV that had latitudes 
# and longitudes for each individual, individual names were added to the spreadsheet and matched manually 
# to the lat/long values for each population. Finally, this CSV was subset to include only the 255 individuals
# present in the native range of western North America.
MIGU_coordinates <- 
  read.csv(file=paste0(MIGU_filePath, 'Geographic/MIGU_LatLongs.csv'), header = TRUE)
# Rename the columns of the geographic coordinates data.frame (because geo.compareBuff function expects certain strings)
colnames(MIGU_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')
# Read in raster data, for SDM
MIGU_sdm <- terra::rast(paste0(MIGU_filePath,'Geographic/MIGU_474inds_rast_Carver.tif'))
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.shp')))
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = MIGU_coordinates)
# Read in the EPA Level III ecoregion shapefile, which is used for calculating ecological coverage 
# (in North America)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level3/NA_CEC_Eco_Level3.shp')))
# Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
# exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
world_poly_clip_W <- wrap(world_poly_clip)
ecoregion_poly_W <- wrap(ecoregion_poly)
MIGU_sdm_W <- wrap(MIGU_sdm)

# ---- GENETIC MATRICES
# For each haplotype length (ranging from a single SNP to 5 concatenated SNPs), the genetic
# data is stored as a STRUCTURE input file. Each STRUCTURE file is read in and converted into a
# genind object.
MIGU_haplos_filePath <- paste0(MIGU_filePath, 'Genetic/haplotypicData')
MIGU_haplos_input <- list.files(MIGU_haplos_filePath, full.names = TRUE)
# Declare a list of integers, storing the number of loci for each STRUCTURE input file
MIGU_lociLevels <- c(1498,608,326,190,115)

# Predeclare an empty list to store genind objects
MIGU_genList <- list()
# For each STRUCTURE input file, 
for(i in 1:length(MIGU_haplos_input)){
  # Read the file in and convert it into a genind object
  MIGU_genind_global <- read.structure(file=MIGU_haplos_input[[i]], n.ind=474, 
                                       n.loc=MIGU_lociLevels[[i]],onerowperind=TRUE, 
                                       col.lab=3, col.pop=1, col.others=2, row.marknames=0, 
                                       NA.char=-9, sep=" ", ask=TRUE, quiet=TRUE)
  
  # REMOVE INTRODUCED POPULATIONS: Subset global genind object to only contain individuals from native range. 
  # The 'drop' argument removes alleles no longer present in the dataset.
  MIGU_genList[[i]] <- MIGU_genind_global[MIGU_coordinates[,1], drop=TRUE]
  # Specify the population of every sample to be 'wild'
  pop(MIGU_genList[[i]]) <- rep('wild', nInd(MIGU_genList[[i]]))
}

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('MIGU_coordinates','MIGU_genList','num_reps','geo_buffSize', 'eco_buffSize',
                              'world_poly_clip_W', 'ecoregion_poly_W', 'MIGU_sdm_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers', 'geo.compareBuff', 'geo.compareBuffSDM', 
                              'eco.intersectBuff', 'eco.compareBuff', 'gen.getAlleleCategories',
                              'calculateCoverage', 'exSituResample.Par', 'geo.gen.Resample.Par'))
# Specify file paths, for saving resampling arrays (if/else statement is for leading zeros)
MIGU_haplos_output <- dirname(gsub("Genetic", "resamplingData", MIGU_haplos_input))
for(i in 1:length(MIGU_haplos_output)){
    MIGU_haplos_output[[i]] <- paste0(MIGU_haplos_output[[i]],'/MIGU_50km_GE_5r_Hap-', i, '_resampArr.Rdata')
}

# Run resampling (in parallel). Use a loop to iterate through the different haplotype lengths
print('%%% BEGINNING RESAMPLING %%%')
for(i in 1:length(MIGU_haplos_output)){
  print(paste0('%%% HAPLOTYPE LENGTH: ', i))
  # Run resampling
  MIGU_demoArray_Par <- 
    geo.gen.Resample.Par(gen_obj = MIGU_genList[[i]], geoFlag = TRUE, coordPts = MIGU_coordinates, 
                         geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag = TRUE, 
                         ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer = 'NA', 
                         reps = num_reps, arrayFilepath = MIGU_haplos_output[[i]], cluster = cl)
}
# Close cores
stopCluster(cl)

# %%% ANALYZE DATA %%% ----
# Specify filepath for MIGU geographic and genetic data, including resampling data
MIGU_filePath <- paste0(GeoGenCorr_wd, 'Datasets/MIGU/')
arrayDir <- paste0(MIGU_filePath, 'resamplingData/haplotypicData/')
# Read in the resampling array .Rdata objects, saved to disk. Then, calculate the average 
# value matrix (across resampling replicates) for each resampling array, and add that 
# matrix as an item to a list.
MIGU_haplos_output <- list.files(arrayDir, full.names = TRUE); MIGU_haplos_averageValueMats <- list()
for(i in 1:length(MIGU_haplos_output)){
  MIGU_haplos_averageValueMats[[i]] <- meanArrayValues(readRDS(MIGU_haplos_output[[i]]))
}
# Build a summary matrix that consists of the Total allelic representation values for each 
# haplotype dataset (5 columns), average Geographic representation values (1 column), and average
# Ecological representation values (1 column)
MIGU_haplos_summaryMat <- array(data=NA, dim=c(nrow(MIGU_haplos_averageValueMats[[1]]),7))
colnames(MIGU_haplos_summaryMat) <- c(paste0(rep('Total_Hap', 5),seq(1:5)),'Geo', 'Eco')
# Total allelic representation values: pull the "Total" columns from each average value matrix,
# and cbind these together
MIGU_haplos_summaryMat[,1:5] <- do.call(cbind,lapply(MIGU_haplos_averageValueMats, function(x) x$Total))
# For the geographic/ecological coverages: build a matrix of the coverage values
# from each average value matrix (using sapply). Then, calculate the means across the 
# rows of this matrix (using rowMeans), and pass this to the last column of the summary matrix
MIGU_haplos_summaryMat[,6] <- rowMeans(sapply(MIGU_haplos_averageValueMats, function(df) df[, 2]))
MIGU_haplos_summaryMat[,7] <- rowMeans(sapply(MIGU_haplos_averageValueMats, function(df) df[, 3]))
# Convert the summary data into a data.frame
MIGU_haplos_summaryMat <- as.data.frame(MIGU_haplos_summaryMat)
# %%% NORMALIZED ROOT MEAN SQUARE ERROR: calculate the NRMSE for each haplotype length, for both
# geographic and ecological coverage
MIGU_NRMSE_Values <- matrix(nrow=length(MIGU_haplos_output), ncol=2)

for(i in 1:nrow(MIGU_NRMSE_Values)){
  MIGU_NRMSE_Values[i,1] <- nrmse_func(MIGU_haplos_summaryMat$Geo, pred=MIGU_haplos_summaryMat[,i])
  MIGU_NRMSE_Values[i,2] <- nrmse_func(MIGU_haplos_summaryMat$Eco, pred=MIGU_haplos_summaryMat[,i]) 
}

# ---- PLOTTING ----
hapColors <- c('gold2','orange','salmon','darkorange2','tomato3','darkblue','purple')
hapColors_Fade <- alpha(hapColors, 0.5)
legText <- c(paste0(rep('Hap. length: ', 5), seq(1:5)), 'Geographic coverage (50 km buffer)', 
             'Ecological coverage (50 km buffer, EPA Level III)')
# ---- COVERAGE PLOTS
matplot(MIGU_haplos_summaryMat, ylim=c(0,110), col=hapColors_Fade, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='M. guttatus: Haplotypic Coverages', line=1.5)
mtext(text='255 Individuals; 50 km Geographic buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.2)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.2, srt=90)
# Add legend
legend(x=120, y=60, inset = 0.05, legend = legText, col=hapColors, pch = c(19,19), 
       cex=1.2, pt.cex = 2, bty='n', y.intersp = 0.5)
