# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% MIMULUS GUTTATUS : COVERAGE OF OPTIMIZED SUBSETS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script utilizes the Mimulus guttatus dataset from Vallejo-Martin et al. 2021 and calculates
# the coverages of "optimized" subsets of individuals. It first calculates the geographic coverage
# of genetic "core sets" (as determined by Chris Richards and Pat Reeves, using the M+ algorithm),
# and also maps the indviduals found in these core subsets. Next, it calculates the genetic coverage
# of geographically optimized subsets, as identified by Dan Carver.

pacman::p_load(adegenet, terra, parallel, RColorBrewer, viridis, scales, vcfR, usedist)

# Read in relevant functions
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/Code/'
setwd(GeoGenCorr_wd)
source('Scripts/functions_GeoGenCoverage.R')
# Specify filepath for MIGU geographic and genetic data
MIGU_filePath <- paste0(GeoGenCorr_wd, 'Datasets/MIGU/')

# %%% GEOGRAPHIC COVERAGE OF GENETICALLY OPTIMIZED SUBSETS %%% ----
# %%% READ IN RELEVANT DATA ----
# %%% GEOGRAPHIC/ECOLOGICAL DATA LAYERS
# Read in coordinate files for MIGU. This CSV was subset to include only the 255 individuals
# present in the native range of western North America.
MIGU_coordinates <- 
  read.csv(file=paste0(MIGU_filePath, 'Geographic/MIGU_coordinates.csv'), header = TRUE)
colnames(MIGU_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.shp')))
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = MIGU_coordinates)
# Declare buffer and point projections
ptProj='+proj=longlat +datum=WGS84'
buffProj='+proj=eqearth +datum=WGS84'
# Read in raster data, for SDM
MIGU_sdm <- terra::rast(paste0(MIGU_filePath,'Geographic/MIGU_255inds_rast_Carver.tif'))
# Read in the EPA Level III ecoregion shapefile, which is used for calculating ecological coverage 
# (in North America)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level3/NA_CEC_Eco_Level3.shp')))

# %%% GENETIC DATA 
# Read in the VCF file provided in Vallejo-Martin et al. 2021 using vcfR::read.vcfR
MIGU_vcf <- 
  read.vcfR(file=paste0(MIGU_filePath, 'Genetic/mgut_all_20180305_gut_filter_75.i50.recode.pruned.plink_20180326.vcf'))
# Convert the vcf to a genind; the return.alleles TRUE value is suggested in the function's help file
# This genind file is made up of 474 individuals and 1,498 loci
MIGU_genind_global <- vcfR2genind(MIGU_vcf, sep = "/", return.alleles = TRUE)
# REMOVE INTRODUCED POPULATIONS: Subset global genind object to only contain individuals from native range. 
# The 'drop' argument removes alleles no longer present in the dataset.
MIGU_genind <- MIGU_genind_global[MIGU_coordinates[,1], drop=TRUE]
# Read in lists of MIGU genetic core sets. These are CSVs adapted from results provided by Chris and Pat
MIGU_90genCoreSets <- read.csv(file=paste0(GeoGenCorr_wd,'../Datasets/Mimulus_guttatus/Genetic/MIGU_genCoreSets_0.9-subsets.csv'))

# %%% CALCULATE GEOGRAPHIC (TOTAL BUFFER) COVERAGES ----
# Declare vector of geographic buffer sizes
geoBuff <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500,1000,1500,2000))
# Build a matrix to store geographic coverage rates (columns different core sets, rows different buffer sizes)
geoBuffRateMat <- matrix(NA, nrow=length(geoBuff), ncol=length(MIGU_90genCoreSets))
colnames(geoBuffRateMat) <- colnames(MIGU_90genCoreSets)
rownames(geoBuffRateMat) <- paste0(geoBuff/1000, 'km')
# Loop through the genetic core sets, calculating the geographic coverages of each one
for(i in 1:length(MIGU_90genCoreSets)){
  coreSet <- MIGU_90genCoreSets[,i]
  geoRates <- 
    lapply(geoBuff, function(x) geo.compareBuff(totalWildPoints=MIGU_coordinates, sampVect=coreSet,
                                                buffSize=x, ptProj=ptProj, buffProj=buffProj, 
                                                boundary=world_poly_clip, parFlag=FALSE))
  # Store the results
  geoBuffRateMat[,i] <- unlist(geoRates)
}
# Save matrix of geographic (total buffer) coverages to disk
write.table(geoBuffRateMat, file=paste0(MIGU_filePath,'MIGU_CoreSets_GeoBuffRates.csv'), 
            row.names=TRUE, col.names=NA, sep=',')

# %%% CALCULATE GEOGRAPHIC (SDM) COVERAGES ----
# Declare vector of geographic buffer sizes
geoSDMBuff <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500,1000,1500,2000))
# The SDM needs to be resampled according to the buffer size, to ensure the buffer isn't smaller than the 
# SDM resolution. This function compares the buffer size to the SDM resolution, resamples the SDM if necessary,
# and then returns a list of SDMs.
SDMrast <- lapply(geoSDMBuff, function(x) geo.checkSDMres(buffSize=x, raster=MIGU_sdm))
# Build a matrix to store geographic coverage rates (columns different core sets, rows different buffer sizes)
geoSDMRateMat <- matrix(NA, nrow=length(geoSDMBuff), ncol=length(MIGU_90genCoreSets))
colnames(geoSDMRateMat) <- colnames(MIGU_90genCoreSets)
rownames(geoSDMRateMat) <- paste0(geoSDMBuff/1000, 'km')
# Loop through the genetic core sets, calculating the geographic coverages of each one
for(i in 1:length(MIGU_90genCoreSets)){
  coreSet <- MIGU_90genCoreSets[,i]
  
  geoSDMRates <- 
    mapply(function(b,r) geo.compareBuffSDM(totalWildPoints=MIGU_coordinates, sampVect=coreSet,
                                            buffSize=b, model=r, ptProj=ptProj, 
                                            buffProj=buffProj, boundary=world_poly_clip), 
           b=geoSDMBuff, r=SDMrast)
  
  # Store the results
  geoSDMRateMat[,i] <- unlist(geoSDMRates)
}
# Save matrix of geographic (SDM) coverages to disk
write.table(geoSDMRateMat, file=paste0(MIGU_filePath,'MIGU_CoreSets_GeoSDMRates.csv'), 
            row.names=TRUE, col.names=NA, sep=',')

# # %%% CALCULATE ECOLOGICAL COVERAGES ----
# # Declare vector of geographic buffer sizes
# ecoBuff <- 1000*(c(5,10))
# ecoBuff <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500,1000,1500,2000))
# # The eco.compareBuff function is designed such that the denominator is provided to it. So, the total count
# # of ecoregions needs to be calculated first
# ecoTotalCount <- lapply(ecoBuff, 
#                         function(x) eco.totalEcoregionCount(totalWildPoints=MIGU_coordinates, buffSize=x,
#                                                             ptProj=ptProj, buffProj=buffProj, 
#                                                             ecoRegion=ecoregion_poly, layerType='NA',
#                                                             boundary=world_poly_clip))

# # !!! This isn't working, and I'm not sure why. 0 ecoregions are being recorded. !!!

# # Build a matrix to store geographic coverage rates (columns different core sets, rows different buffer sizes)
# ecoRateMat <- matrix(NA, nrow=length(ecoBuff), ncol=length(MIGU_90genCoreSets))
# colnames(ecoRateMat) <- colnames(MIGU_90genCoreSets)
# rownames(ecoRateMat) <- paste0(ecoBuff/1000, 'km')
# # Loop through the genetic core sets, calculating the geographic coverages of each one
# for(i in 1:length(MIGU_90genCoreSets)){
#   coreSet <- MIGU_90genCoreSets[,i]
#   ecoRates <- 
#     lapply(ecoBuff, function(x) eco.compareBuff(totalWildPoints=MIGU_coordinates, sampVect=coreSet,
#                                                 buffSize=x, ecoTotalCount, ptProj=ptProj, buffProj=buffProj, 
#                                                 ecoRegion=ecoregion_poly, layerType='NA',parFlag=FALSE))
#   # Store the results
#   ecoRateMat[,i] <- unlist(geoecoRates)
# }
# # Save matrix of ecological coverages to disk
# write.table(ecoRateMat, file=paste0(MIGU_filePath,'MIGU_CoreSets_EcoRates.csv'), 
#             row.names=TRUE, col.names=NA, sep=',')

# # %%% BUILDING VCF OF STRICTLY NATIVE SAMPLES %%% ----
# # This is code which only had to be run once, in order to provide the necessary VCF for determining the genetic "core sets"
# # It is kept here for documentation purposes.
# # Read in coordinate file of strictly individuals from North America (255)
# MIGU_coordinates <- 
#   read.csv(file=paste0(MIGU_filePath, 'Geographic/MIGU_coordinates.csv'), header = TRUE)
# colnames(MIGU_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')
# # Read in original VCF, containing native and international MIGU individuals (474)
# MIGU_vcf <- 
#   read.vcfR(file=paste0(MIGU_filePath, 'Genetic/mgut_all_20180305_gut_filter_75.i50.recode.pruned.plink_20180326.vcf'))
# # Build a list of the names of international samples
# MIGU_intNames <- colnames(MIGU_vcf@gt)[which(!(colnames(MIGU_vcf@gt) %in% MIGU_coordinates$Sample.Name))]
# # Remove the "FORMAT" argument from the list
# MIGU_intNames <- MIGU_intNames[-1]
# # Checking that the length is correct
# 474 - length(MIGU_intNames) == 255
# # Append the vcftools "--remove-indv" flag before each sample name
# BCF_text <- paste0(rep('^', length(MIGU_intNames)), MIGU_intNames)
# # Collapse text into a single, long string
# BCF_text <- paste(BCF_text, collapse=" ")

# %%% GENETIC COVERAGE OF GEOGRAPHICALLY OPTIMIZED SUBSETS %%% ----
# %%% READ IN RELEVANT DATA ----
# %%% GEOGRAPHIC/ECOLOGICAL DATA LAYERS
# Read in coordinate files for MIGU. This CSV was subset to include only the 255 individuals
# present in the native range of western North America.
MIGU_coordinates <- 
  read.csv(file=paste0(MIGU_filePath, 'Geographic/MIGU_coordinates.csv'), header = TRUE)
colnames(MIGU_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.shp')))
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = MIGU_coordinates)
# Declare buffer and point projections
ptProj='+proj=longlat +datum=WGS84'
buffProj='+proj=eqearth +datum=WGS84'
# Read in raster data, for SDM
MIGU_sdm <- terra::rast(paste0(MIGU_filePath,'Geographic/MIGU_255inds_rast_Carver.tif'))
# Read in the EPA Level III ecoregion shapefile, which is used for calculating ecological coverage 
# (in North America)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level3/NA_CEC_Eco_Level3.shp')))

# %%% GENETIC DATA 
# Read in the VCF file provided in Vallejo-Martin et al. 2021 using vcfR::read.vcfR
MIGU_vcf <- 
  read.vcfR(file=paste0(MIGU_filePath, 'Genetic/mgut_all_20180305_gut_filter_75.i50.recode.pruned.plink_20180326.vcf'))
# Convert the vcf to a genind; the return.alleles TRUE value is suggested in the function's help file
# This genind file is made up of 474 individuals and 1,498 loci
MIGU_genind_global <- vcfR2genind(MIGU_vcf, sep = "/", return.alleles = TRUE)
# REMOVE INTRODUCED POPULATIONS: Subset global genind object to only contain individuals from native range. 
# The 'drop' argument removes alleles no longer present in the dataset.
MIGU_genind <- MIGU_genind_global[MIGU_coordinates[,1], drop=TRUE]

# Read in lists of MIGU geographic core sets. These are CSVs adapted from results provided by Dan Carver,
# with a unique group of core sets for each buffer size (5km, 10km, 25km, 50km, 250km, 500km)
MIGU_5km_geoCoreSets <- 
  read.csv(file='/home/akoontz/Documents/GeoGenCorr/Datasets/Mimulus_guttatus/Geographic/MIGU_geoCoreSets_005km.csv')
MIGU_10km_geoCoreSets <- 
  read.csv(file='/home/akoontz/Documents/GeoGenCorr/Datasets/Mimulus_guttatus/Geographic/MIGU_geoCoreSets_010km.csv')
MIGU_25km_geoCoreSets <- 
  read.csv(file='/home/akoontz/Documents/GeoGenCorr/Datasets/Mimulus_guttatus/Geographic/MIGU_geoCoreSets_025km.csv')
MIGU_50km_geoCoreSets <- 
  read.csv(file='/home/akoontz/Documents/GeoGenCorr/Datasets/Mimulus_guttatus/Geographic/MIGU_geoCoreSets_050km.csv')
MIGU_250km_geoCoreSets <- 
  read.csv(file='/home/akoontz/Documents/GeoGenCorr/Datasets/Mimulus_guttatus/Geographic/MIGU_geoCoreSets_250km.csv')
MIGU_500km_geoCoreSets <- 
  read.csv(file='/home/akoontz/Documents/GeoGenCorr/Datasets/Mimulus_guttatus/Geographic/MIGU_geoCoreSets_500km.csv')
# Make a list of both geograhpic core sets, to loop through later
MIGU_geoCoreSets  <- list(MIGU_5km_geoCoreSets, MIGU_10km_geoCoreSets, MIGU_25km_geoCoreSets, 
                          MIGU_50km_geoCoreSets, MIGU_250km_geoCoreSets, MIGU_500km_geoCoreSets)

# %%% CALCULATE GENETIC COVERAGES ----
# Build a matrix to store genetic coverage rates (columns different core sets, rows different buffer sizes)
genRateMat <- matrix(NA, nrow=length(MIGU_geoCoreSets), ncol=length(MIGU_geoCoreSets[[1]]))
colnames(genRateMat) <- colnames(MIGU_5km_geoCoreSets)
rownames(genRateMat) <- c('5km','10km','25km','50km','250km','500km')
# Loop through the geographic core sets, calculating the genetic coverages 
# for each buffer size
for(i in 1:length(MIGU_geoCoreSets)){
  # Specify which buffer size to analyze
  coreSet <- MIGU_geoCoreSets[[i]]
  # Loop through the subsets, of the specified core set
  for(j in 1:length(coreSet)){
    # Extract the same names in the relevant subset, and build the corresponding
    # genind object with the sample names (dropping absent loci)
    coreSubset <- coreSet[,j]
    MIGU_genSubset <- MIGU_genind[coreSubset,,drop=TRUE]
    # With the given genind object, calculate the geographic coverages, compared
    # to the complete genind object
    genRateMat[i,j] <- 
      gen.getAlleleCategories(genMat=MIGU_genind@tab, samp=MIGU_genSubset@tab)[1,3]
  }
}

# Save matrix of genetic coverages to disk
write.table(genRateMat, file=paste0(MIGU_filePath,'MIGU_CoreSets_GenRates.csv'), 
            row.names=TRUE, sep=',')

# %%% RESAMPLING FOR RANDOMIZED GENETIC COVERAGES (10 REPLICATES) %%% ----
# Specify geographic buffer size in meters 
geo_buffSize <- 1000*(c(2,3,4,5,seq(10,100,5),seq(110,250,10),500,1000,1500,2000))

# ---- PARALLELIZATION
# Set up relevant cores 
num_cores <- detectCores() - 4 
cl <- makeCluster(num_cores)
# Make sure libraries (adegenet, terra, etc.) are on cluster (but avoid printing output)
invisible(clusterEvalQ(cl, library('adegenet')))
invisible(clusterEvalQ(cl, library('terra')))
invisible(clusterEvalQ(cl, library('parallel')))
invisible(clusterEvalQ(cl, library('usedist')))
invisible(clusterEvalQ(cl, library('ape')))

# %%% CONDUCT RESAMPLING ----
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
  read.csv(file=paste0(MIGU_filePath, 'Geographic/MIGU_coordinates.csv'), header = TRUE)
# Rename the columns of the geographic coordinates data.frame (because geo.compareBuff function expects certain strings)
colnames(MIGU_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')
# Read in raster data, for SDM
MIGU_sdm <- terra::rast(paste0(MIGU_filePath,'Geographic/MIGU_255inds_rast_Carver.tif'))

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

# ---- GENETIC MATRIX
# Read in the VCF file provided in Vallejo-Martin et al. 2021 using vcfR::read.vcfR
MIGU_vcf <- 
  read.vcfR(file=paste0(MIGU_filePath, 'Genetic/mgut_all_20180305_gut_filter_75.i50.recode.pruned.plink_20180326.vcf'))
# Convert the vcf to a genind; the return.alleles TRUE value is suggested in the function's help file
# This genind file is made up of 474 individuals and 1,498 loci
MIGU_genind_global <- vcfR2genind(MIGU_vcf, sep = "/", return.alleles = TRUE)
# REMOVE INTRODUCED POPULATIONS: Subset global genind object to only contain individuals from native range. 
# The 'drop' argument removes alleles no longer present in the dataset.
MIGU_genind <- MIGU_genind_global[MIGU_coordinates[,1], drop=TRUE]

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('MIGU_coordinates','MIGU_genind','num_reps','geo_buffSize','eco_buffSize',
                              'world_poly_clip_W','ecoregion_poly_W','MIGU_sdm_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres', 
                              'eco.intersectBuff','eco.compareBuff','gen.getAlleleCategories', 
                              'gen.buildDistMat', 'gen.calcGenDistCov', 'eco.totalEcoregionCount',
                              'calculateCoverage','exSituResample.Par', 'geo.gen.Resample.Par'))
# Specify file path, for saving resampling array
GeoOptarrayDir <- paste0(MIGU_filePath, 'resamplingData/MIGU_GeoOpt_G2_10r_resampArr.Rdata')

# Run resampling for geographic optimization comparison
MIGU_demoArray_Par <- 
  geo.gen.Resample.Par(genObj = MIGU_genind, geoFlag = TRUE, coordPts = MIGU_coordinates, 
                       geoBuff = geo_buffSize, SDMrast = MIGU_sdm_W, 
                       boundary=world_poly_clip_W, ecoFlag = FALSE, reps = num_reps, 
                       arrayFilepath = GeoOptarrayDir, cluster = cl)
# Close cores
stopCluster(cl)

# %%% PROCESS THE RESAMPLING ARRAY ----
# Read in the resampling array .Rdata object, saved to disk
MIGU_geoOptArray <- readRDS(GeoOptarrayDir)

# Extract only the rows necessary: 7, 11, 49, 59, 63, and 67 individuals
MIGU_geoRandSubsets <- MIGU_geoOptArray[c(6,10,48,68,62,66),1,]
colnames(MIGU_geoRandSubsets) <- paste0("Subset",1:10)
rownames(MIGU_geoRandSubsets) <- c('07 individuals','11 individuals','49 individuals',
                                  '59 individuals','63 individuals','67 individuals')

# Save matrix of randomized genetic coverages to disk
write.table(MIGU_geoRandSubsets, file=paste0(MIGU_filePath,'MIGU_RandomSets_GenRates.csv'), 
            row.names=TRUE, sep=',')
