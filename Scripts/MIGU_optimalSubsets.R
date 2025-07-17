# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% MIMULUS GUTTATUS : COVERAGE OF OPTIMIZED SUBSETS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script utilizes the Mimulus guttatus dataset from Vallejo-Martin et al. 2021 and calculates
# the coverages of "optimized" subsets of individuals. It first calculates the genetic coverage
# of geographically-optimized sample sets, as identified by Dan Carver. Then, it compares those
# "geographically optimized" genetic coverages twice: once to "truly randomized" genetic coverages,
# which uses the random sampling built into the GeoGenCor workflow (where two individuals from the
# same population might be selected). Second, to a "semi randomized" set of genetic coverages, where
# individuals are randomly selected but each population is represented only once.

# Finally, it also calculates the geographic coverage of genetic "core sets" 
# (as determined by Chris Richards and Pat Reeves, using the M+ algorithm),
# and also maps the indviduals found in these core subsets. This code has been archived and commented
# out at the bottom of the script, as the primary focus of this script is the geographic optimization.

pacman::p_load(adegenet, terra, parallel, RColorBrewer, viridis, scales, vcfR, usedist)

# Read in relevant functions
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/Code/'
setwd(GeoGenCorr_wd)
source('Scripts/functions_GeoGenCoverage.R')
# Specify filepath for MIGU geographic and genetic data
MIGU_filePath <- paste0(GeoGenCorr_wd, 'Datasets/MIGU/')

# %%% GENETIC COVERAGE OF GEOGRAPHICALLY OPTIMIZED SUBSETS %%% ----
# %%% READ IN RELEVANT DATA ----
# %%% GEOGRAPHIC/ECOLOGICAL DATA LAYERS
# Read in coordinate files for MIGU. This CSV was subset to include only the 255 individuals
# present in the native range of western North America.
MIGU_coordinates <-
  read.csv(file=paste0(MIGU_filePath, 'Geographic/MIGU_coordinates.csv'), header = TRUE)
colnames(MIGU_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')

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
MIGU_geoCoreSets  <- list(MIGU_500km_geoCoreSets, MIGU_250km_geoCoreSets, MIGU_50km_geoCoreSets,
                          MIGU_25km_geoCoreSets, MIGU_10km_geoCoreSets, MIGU_5km_geoCoreSets)

# %%% CALCULATE GENETIC COVERAGES ----
# Build a matrix to store genetic coverage rates (columns different core sets, rows different buffer sizes)
MIGU_geoOptSubsets <- matrix(NA, nrow=length(MIGU_geoCoreSets), ncol=length(MIGU_geoCoreSets[[1]]))
colnames(MIGU_geoOptSubsets) <- paste0("Subset",1:10)
rownames(MIGU_geoOptSubsets) <- c('7n_500km','11n_250km','49n_50km', '59n_25km','63n_10km','67n_5km')
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
    MIGU_geoOptSubsets[i,j] <-
      gen.getAlleleCategories(genMat=MIGU_genind@tab, samp=MIGU_genSubset@tab)[1,3]
  }
}

# Save matrix of genetic coverages to disk
write.table(MIGU_geoOptSubsets, file=paste0(MIGU_filePath,'MIGU_CoreSets_GenRates.csv'),
            row.names=TRUE, sep=',')

# %%% RESAMPLING FOR RANDOMIZED GENETIC COVERAGES (10 REPLICATES) %%% ----
# Read in relevant functions
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/Code/'
setwd(GeoGenCorr_wd)
source('Scripts/functions_GeoGenCoverage.R')
# Specify filepath for MIGU geographic and genetic data
MIGU_filePath <- paste0(GeoGenCorr_wd, 'Datasets/MIGU/')
# Declare variable for whether or not resampling needs to be run
resamplingFlag <- FALSE
# If resampling has already occurred, read in the array; otherwise, run resampling workflow
if(resamplingFlag==TRUE){
  # Specify resampling replicates
  num_reps <- 10

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
  # Read in coordinates dataframe. This is necessary for properly subsetting the genetic data
  MIGU_coordinates <-
    read.csv(file=paste0(MIGU_filePath, 'Geographic/MIGU_coordinates.csv'), header = TRUE)
  # Rename the columns of the geographic coordinates data.frame (because geo.compareBuff function expects certain strings)
  colnames(MIGU_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')

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
  # Export necessary objects to cluster
  clusterExport(cl, varlist = c('MIGU_coordinates','MIGU_genind','num_reps'))
  # Export necessary functions (for calculating geographic and ecological coverage) to the cluster
  clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres',
                                'eco.intersectBuff','eco.compareBuff','gen.getAlleleCategories',
                                'gen.buildDistMat', 'gen.calcGenDistCov', 'eco.totalEcoregionCount',
                                'calculateCoverage','exSituResample.Par', 'geo.gen.Resample.Par'))
  # Specify file path, for saving resampling array
  GeoOpt_Rand_arrayDir <- paste0(MIGU_filePath, 'resamplingData/MIGU_GeoOpt_RandGen_10r_resampArr.Rdata')

  # Run resampling for geographic optimization comparison
  MIGU_GeoOpt_RandArr <-
    geo.gen.Resample.Par(genObj = MIGU_genind, geoFlag = FALSE,  reps = num_reps,
                         arrayFilepath = GeoOptarrayDir, cluster = cl)
  # Close cores
  stopCluster(cl)
} else {
  # Specify file path, for saving resampling array
  GeoOpt_Rand_arrayDir <- paste0(MIGU_filePath, 'resamplingData/MIGU_GeoOpt_RandGen_10r_resampArr.Rdata')
  # Read in the existin resampling array .Rdata object, saved to disk
  MIGU_GeoOpt_RandArr <- readRDS(GeoOpt_Rand_arrayDir)
}

# %%% PROCESS THE RESAMPLING ARRAY ----
# Extract only the rows necessary: 7, 11, 49, 59, 63, and 67 individuals
MIGU_geoRandSubsets <- MIGU_GeoOpt_RandArr[c(6,10,48,58,62,66),1,]
colnames(MIGU_geoRandSubsets) <- paste0("Subset",1:10)
rownames(MIGU_geoRandSubsets) <- c('7n_500km','11n_250km','49n_50km', '59n_25km','63n_10km','67n_5km')
# Save matrix of randomized genetic coverages to disk
write.table(MIGU_geoRandSubsets, file=paste0(MIGU_filePath,'MIGU_RandomSets_GenRates.csv'),
            row.names=TRUE, sep=',')

# %%% COMPARE OPTIMIZED AND RANDOMIZED COVERAGES ----
# Generate a matrix for storing results of Wilcox comparisons
WilcoxResultsMat <- matrix(nrow=nrow(MIGU_geoOptSubsets), ncol=2)
rownames(WilcoxResultsMat) <- rownames(MIGU_geoOptSubsets)
colnames(WilcoxResultsMat) <- c('Statistic', 'P_value')
# Loop through the rows (sample sizes) of the optimal coverage values
for(i in 1:nrow(MIGU_geoOptSubsets)){
  # Conduct the Wilcox test. Here, we're explicitly testing whether the
  # geographically optimized coverage values (X) are greater than the
  # randomized coverage values (Y).
  result <-
    wilcox.test(MIGU_geoOptSubsets[i,], MIGU_geoRandSubsets[i,], alternative = "greater")
  # The warning here "cannot compute exact p-value with ties" is because
  # there are some values which match between the two datasets, which leads
  # to a normal approximation; this should be safe for our data.
  WilcoxResultsMat[i,1] <- result$statistic
  WilcoxResultsMat[i,2] <- result$p.value
}
print(WilcoxResultsMat)

# %%% RESAMPLING FOR SEMI-RANDOMIZED GENETIC COVERAGES (10 REPLICATES) %%% ----
# This approach begins by first randomly selecting a single individual from each population, and subsetting
# the dataset down to these randomly selected individuals. Then, using these subsets of individuals, resampling
# is run (using the standard functions written). Because we want to randomly select individuals from each population
# independently ten times, a single resampling replicate has to be generated each time.

# Read in relevant functions
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/Code/'
setwd(GeoGenCorr_wd)
source('Scripts/functions_GeoGenCoverage.R')
# Specify filepath for MIGU geographic and genetic data
MIGU_filePath <- paste0(GeoGenCorr_wd, 'Datasets/MIGU/')
# Declare variable for whether or not resampling needs to be run
resamplingFlag <- FALSE
# If resampling has already occurred, read in the array; otherwise, run resampling workflow (not in parallel)
if(resamplingFlag==TRUE){
  # Specify resampling replicates; this will determine the number of loop iterations
  num_reps <- 10
  # Predeclare an empty array, to receive the results of each resampling run in each loop iteration
  MIGU_GeoOpt_semiRandArr <- array(NA, dim=c(99,7,10))
  
    # ---- READ IN DATA ----
  # Read in coordinates dataframe
  MIGU_coordinates <- 
    read.csv(file=paste0(MIGU_filePath, 'Geographic/MIGU_coordinates.csv'), header = TRUE)
  # Rename geographic dataframe columns, and sort
  colnames(MIGU_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')
  MIGU_coordinates <- MIGU_coordinates[order(MIGU_coordinates$Sample.Name),]
  # Read in the VCF file provided in Vallejo-Martin et al. 2021 using vcfR::read.vcfR
  MIGU_vcf <- 
    read.vcfR(file=paste0(MIGU_filePath, 'Genetic/mgut_all_20180305_gut_filter_75.i50.recode.pruned.plink_20180326.vcf'))
  # Convert the vcf to a genind; the return.alleles TRUE value is suggested in the function's help file
  # This genind file is made up of 474 individuals and 1,498 loci
  MIGU_genind_global <- vcfR2genind(MIGU_vcf, sep = "/", return.alleles = TRUE)
  # REMOVE INTRODUCED POPULATIONS: Subset global genind object to only contain individuals from native range. 
  # The 'drop' argument removes alleles no longer present in the dataset.
  MIGU_genind <- MIGU_genind_global[MIGU_coordinates[,1], drop=TRUE]
  
  for(i in 1:num_reps){
    # ---- SUBSETTING TO SINGLE INDIVIDUAL PER POPULATION
    # Read in table from Supplementary Data 1, which lists each population separately
    MIGUpops <- read.csv('/home/akoontz/Documents/GeoGenCorr/Datasets/Mimulus_guttatus/Geographic/SuppData1.csv')
    # Subset to only Alaskan and North American populations
    MIGUpops <- MIGUpops[which(MIGUpops$Region == 'ak' | MIGUpops$Region == 'nam'),]
    # There are 4 populations which are present ni this SuppData1 file, but not within the VCF of samples:
    # CAB, CMD, LIN, and 15_NAU. So, strike these population names from the MIGUpops vector
    MIGUpops <- MIGUpops[-which(MIGUpops$Population == 'LIN' | MIGUpops$Population == 'CAB' | 
                                  MIGUpops$Population == 'CMD' | MIGUpops$Population == '15_NAU'),]
    # Generate a vector listing each population once
    MIGUpops <- sort(MIGUpops$Population)
    # Predeclare an empty vector, which will be used to store the name of samples we want to retain
    sampsFromEachPop <- vector(length = length(MIGUpops))
    # Loop through the names of unique populations
    for(j in 1:length(MIGUpops)){
      # Extract all the individuals belonging to a certain population
      popInds <- MIGU_coordinates[grep(MIGUpops[[j]], MIGU_coordinates$Sample.Name),]
      # Randomly sample one of the individuals from that population
      sampsFromEachPop[j] <- popInds[sample(nrow(popInds), size=1),]$Sample.Name
    }
    # Now, subset the geographic coordinate dataframe and the genind object by this sample set
    MIGU_coordinatesSub <- MIGU_coordinates[MIGU_coordinates$Sample.Name %in% sampsFromEachPop,]
    MIGU_genindSub <- MIGU_genind[MIGU_coordinatesSub[,1], drop=TRUE]
    
    # Run semi-random resampling, to obtain genetic values for geographic comparison
    MIGU_GeoOpt_semiRandArr[,,i] <- 
      geo.gen.Resample(genObj = MIGU_genindSub, geoFlag = FALSE, reps = 1)
  }
  # Specify file path, for saving resampling array
  MIGU_GeoOpt_semiRand_arrayDir <- paste0(MIGU_filePath, 'resamplingData/MIGU_GeoOpt_SemiRandGen_10r_resampArr.Rdata')
  # Save resampling array to disk, dropping last two columns (which aren't used)
  saveRDS(MIGU_GeoOpt_semiRandArr[,-(6:7),], MIGU_GeoOpt_semiRand_arrayDir)
} else {
  # Specify file path, for saving resampling array
  MIGU_GeoOpt_semiRand_arrayDir <- paste0(MIGU_filePath, 'resamplingData/MIGU_GeoOpt_SemiRandGen_10r_resampArr.Rdata')
  # Read in the existin resampling array .Rdata object, saved to disk
  MIGU_GeoOpt_semiRandArr <- readRDS(MIGU_GeoOpt_semiRand_arrayDir)
}

# %%% PROCESS THE RESAMPLING ARRAY ----
# Extract only the rows necessary: 7, 11, 49, 59, 63, and 67 individuals
MIGU_geoSemiRandSubsets <- MIGU_GeoOpt_semiRandArr[c(6,10,48,58,62,66),1,]
colnames(MIGU_geoSemiRandSubsets) <- paste0("Subset",1:10)
rownames(MIGU_geoSemiRandSubsets) <- c('7n_500km','11n_250km','49n_50km', '59n_25km','63n_10km','67n_5km')
# Save matrix of randomized genetic coverages to disk
write.table(MIGU_geoSemiRandSubsets, file=paste0(MIGU_filePath,'MIGU_SemiRandomSets_GenRates.csv'),
            row.names=TRUE, sep=',')

# %%% COMPARE OPTIMIZED AND RANDOMIZED COVERAGES ----
# Generate a matrix for storing results of Wilcox comparisons
WilcoxResultsMat <- matrix(nrow=nrow(MIGU_geoOptSubsets), ncol=2)
rownames(WilcoxResultsMat) <- rownames(MIGU_geoOptSubsets)
colnames(WilcoxResultsMat) <- c('Statistic', 'P_value')
# Loop through the rows (sample sizes) of the optimal coverage values
for(i in 1:nrow(MIGU_geoOptSubsets)){
  # Conduct the Wilcox test. Here, we're explicitly testing whether the
  # geographically optimized coverage values (X) are greater than the
  # randomized coverage values (Y).
  result <-
    wilcox.test(MIGU_geoOptSubsets[i,], MIGU_geoSemiRandSubsets[i,], alternative = "greater")
  # The warning here "cannot compute exact p-value with ties" is because
  # there are some values which match between the two datasets, which leads
  # to a normal approximation; this should be safe for our data.
  WilcoxResultsMat[i,1] <- result$statistic
  WilcoxResultsMat[i,2] <- result$p.value
}
print(WilcoxResultsMat)

# %%% Test: Examining how many unique populations are represented in the MIGU_coordinates dataframe ----
unqPops <- gsub("[0-9]+$", "", MIGU_coordinates$Sample.Name)
unqPops <- gsub("_\\d+_N", "", unqPops)
unqPops <- gsub("_[A-Za-z]$", "", unqPops)
unqPops <- gsub("_[0-9]+$", "", unqPops)
unqPops <- sort(unique(unqPops))
length(unqPops) # 102 unique populations
length(MIGUpops) # 100 unique populations
unqPops[which(!(unqPops %in% MIGUpops))] # Difference is 16_AMS1 and 16_BHS1, which are actually samples
# from the populations 16_AMS and 16_DHS. So, there's alignment between the 2 vectors

# # %%% GEOGRAPHIC COVERAGE OF GENETICALLY OPTIMIZED SUBSETS %%% ----
# # %%% READ IN RELEVANT DATA 
# # %%% GEOGRAPHIC/ECOLOGICAL DATA LAYERS
# # Read in coordinate files for MIGU. This CSV was subset to include only the 255 individuals
# # present in the native range of western North America.
# MIGU_coordinates <-
#   read.csv(file=paste0(MIGU_filePath, 'Geographic/MIGU_coordinates.csv'), header = TRUE)
# colnames(MIGU_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')
# # Read in world countries layer (created as part of the gap analysis workflow)
# # This layer is used to clip buffers, to make sure they're not in the water
# world_poly_clip <-
#   vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.shp')))
# # Perform geographic filter on the admin layer.
# world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = MIGU_coordinates)
# # Declare buffer and point projections
# ptProj='+proj=longlat +datum=WGS84'
# buffProj='+proj=eqearth +datum=WGS84'
# # Read in raster data, for SDM
# MIGU_sdm <- terra::rast(paste0(MIGU_filePath,'Geographic/MIGU_255inds_rast_Carver.tif'))
# # Read in the EPA Level III ecoregion shapefile, which is used for calculating ecological coverage
# # (in North America)
# ecoregion_poly <-
#   vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level3/NA_CEC_Eco_Level3.shp')))
# 
# # %%% GENETIC DATA
# # Read in the VCF file provided in Vallejo-Martin et al. 2021 using vcfR::read.vcfR
# MIGU_vcf <-
#   read.vcfR(file=paste0(MIGU_filePath, 'Genetic/mgut_all_20180305_gut_filter_75.i50.recode.pruned.plink_20180326.vcf'))
# # Convert the vcf to a genind; the return.alleles TRUE value is suggested in the function's help file
# # This genind file is made up of 474 individuals and 1,498 loci
# MIGU_genind_global <- vcfR2genind(MIGU_vcf, sep = "/", return.alleles = TRUE)
# # REMOVE INTRODUCED POPULATIONS: Subset global genind object to only contain individuals from native range.
# # The 'drop' argument removes alleles no longer present in the dataset.
# MIGU_genind <- MIGU_genind_global[MIGU_coordinates[,1], drop=TRUE]
# # Read in lists of MIGU genetic core sets. These are CSVs adapted from results provided by Chris and Pat
# MIGU_90genCoreSets <- read.csv(file=paste0(GeoGenCorr_wd,'../Datasets/Mimulus_guttatus/Genetic/MIGU_genCoreSets_0.9-subsets.csv'))
# 
# # %%% CALCULATE GEOGRAPHIC (TOTAL BUFFER) COVERAGES 
# # Declare vector of geographic buffer sizes
# geoBuff <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500,1000,1500,2000))
# # Build a matrix to store geographic coverage rates (columns different core sets, rows different buffer sizes)
# geoBuffRateMat <- matrix(NA, nrow=length(geoBuff), ncol=length(MIGU_90genCoreSets))
# colnames(geoBuffRateMat) <- colnames(MIGU_90genCoreSets)
# rownames(geoBuffRateMat) <- paste0(geoBuff/1000, 'km')
# # Loop through the genetic core sets, calculating the geographic coverages of each one
# for(i in 1:length(MIGU_90genCoreSets)){
#   coreSet <- MIGU_90genCoreSets[,i]
#   geoRates <-
#     lapply(geoBuff, function(x) geo.compareBuff(totalWildPoints=MIGU_coordinates, sampVect=coreSet,
#                                                 buffSize=x, ptProj=ptProj, buffProj=buffProj,
#                                                 boundary=world_poly_clip, parFlag=FALSE))
#   # Store the results
#   geoBuffRateMat[,i] <- unlist(geoRates)
# }
# # Save matrix of geographic (total buffer) coverages to disk
# write.table(geoBuffRateMat, file=paste0(MIGU_filePath,'MIGU_CoreSets_GeoBuffRates.csv'),
#             row.names=TRUE, col.names=NA, sep=',')
# 
# # %%% CALCULATE GEOGRAPHIC (SDM) COVERAGES
# # Declare vector of geographic buffer sizes
# geoSDMBuff <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500,1000,1500,2000))
# # The SDM needs to be resampled according to the buffer size, to ensure the buffer isn't smaller than the
# # SDM resolution. This function compares the buffer size to the SDM resolution, resamples the SDM if necessary,
# # and then returns a list of SDMs.
# SDMrast <- lapply(geoSDMBuff, function(x) geo.checkSDMres(buffSize=x, raster=MIGU_sdm))
# # Build a matrix to store geographic coverage rates (columns different core sets, rows different buffer sizes)
# geoSDMRateMat <- matrix(NA, nrow=length(geoSDMBuff), ncol=length(MIGU_90genCoreSets))
# colnames(geoSDMRateMat) <- colnames(MIGU_90genCoreSets)
# rownames(geoSDMRateMat) <- paste0(geoSDMBuff/1000, 'km')
# # Loop through the genetic core sets, calculating the geographic coverages of each one
# for(i in 1:length(MIGU_90genCoreSets)){
#   coreSet <- MIGU_90genCoreSets[,i]
#   
#   geoSDMRates <-
#     mapply(function(b,r) geo.compareBuffSDM(totalWildPoints=MIGU_coordinates, sampVect=coreSet,
#                                             buffSize=b, model=r, ptProj=ptProj,
#                                             buffProj=buffProj, boundary=world_poly_clip),
#            b=geoSDMBuff, r=SDMrast)
#   
#   # Store the results
#   geoSDMRateMat[,i] <- unlist(geoSDMRates)
# }
# # Save matrix of geographic (SDM) coverages to disk
# write.table(geoSDMRateMat, file=paste0(MIGU_filePath,'MIGU_CoreSets_GeoSDMRates.csv'),
#             row.names=TRUE, col.names=NA, sep=',')
# 
# # # %%% CALCULATE ECOLOGICAL COVERAGES 
# # # Declare vector of geographic buffer sizes
# # ecoBuff <- 1000*(c(5,10))
# # ecoBuff <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500,1000,1500,2000))
# # # The eco.compareBuff function is designed such that the denominator is provided to it. So, the total count
# # # of ecoregions needs to be calculated first
# # ecoTotalCount <- lapply(ecoBuff,
# #                         function(x) eco.totalEcoregionCount(totalWildPoints=MIGU_coordinates, buffSize=x,
# #                                                             ptProj=ptProj, buffProj=buffProj,
# #                                                             ecoRegion=ecoregion_poly, layerType='NA',
# #                                                             boundary=world_poly_clip))
# 
# # # !!! This isn't working, and I'm not sure why. 0 ecoregions are being recorded. !!!
# 
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
# 
# # # %%% BUILDING VCF OF STRICTLY NATIVE SAMPLES %%% ----
# # # This is code which only had to be run once, in order to provide the necessary VCF for determining the genetic "core sets"
# # # It is kept here for documentation purposes.
# # # Read in coordinate file of strictly individuals from North America (255)
# # MIGU_coordinates <-
# #   read.csv(file=paste0(MIGU_filePath, 'Geographic/MIGU_coordinates.csv'), header = TRUE)
# # colnames(MIGU_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')
# # # Read in original VCF, containing native and international MIGU individuals (474)
# # MIGU_vcf <-
# #   read.vcfR(file=paste0(MIGU_filePath, 'Genetic/mgut_all_20180305_gut_filter_75.i50.recode.pruned.plink_20180326.vcf'))
# # # Build a list of the names of international samples
# # MIGU_intNames <- colnames(MIGU_vcf@gt)[which(!(colnames(MIGU_vcf@gt) %in% MIGU_coordinates$Sample.Name))]
# # # Remove the "FORMAT" argument from the list
# # MIGU_intNames <- MIGU_intNames[-1]
# # # Checking that the length is correct
# # 474 - length(MIGU_intNames) == 255
# # # Append the vcftools "--remove-indv" flag before each sample name
# # BCF_text <- paste0(rep('^', length(MIGU_intNames)), MIGU_intNames)
# # # Collapse text into a single, long string
# # BCF_text <- paste(BCF_text, collapse=" ")
