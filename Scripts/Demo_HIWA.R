# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: HIBISCUS WAIMEAE %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script demonstrating approach for calculating the correlation between genetic, geographic, and ecological coverage. 
# Uses data files obtained from Jeremie Fant, Seana Walsh, Susan Deans, and others, of Hibiscus waimeae individuals.

pacman::p_load(adegenet, terra, parallel, RColorBrewer, scales, vcfR, usedist)

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

# %%% CONDUCT RESAMPLING %%% ----
# ---- READ IN DATA ----
# Specify filepath for HIWA geographic and genetic data
HIWA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/HIWA/')

# ---- GEOGRAPHIC/ECOLOGICAL DATA FILES
# The original coordinates file for this Hibiscus dataset needs to be processed such that 
# ex situ individuals are removed.
# Check if the processed file (called HIWA_coordinates.csv) already exists; if not, then 
# run the necessary processing steps.
if(file.exists(paste0(HIWA_filePath, 'Geographic/HIWA_coordinates.csv'))){
  # Read in the CSV of processed coordinates. The first column contains row numbers
  HIWA_coordinates <- read.csv(
    paste0(HIWA_filePath, 'Geographic/HIWA_coordinates.csv'), header=TRUE)
} else {
  # Read in wild occurrence points. This CSV has 5 columns: sample name, sample code, locality, 
  # latitude, and longitude
  HIWA_points <- read.csv(paste0(HIWA_filePath, 'Geographic/HIWA_coordinates_Original.csv'), header=TRUE)
  # Remove the ex situ samples from the dataframe, and then drop sample code and locality columns
  HIWA_points <- HIWA_points[-which(HIWA_points$Locality=='Ex Situ'),]
  HIWA_coordinates <- HIWA_points[,-(2:3)]
  # Write resulting coordinates data.frame as a CSV to disk, for future runs
  write.csv(HIWA_coordinates, file=paste0(HIWA_filePath,'Geographic/HIWA_coordinates.csv'), 
            row.names = FALSE)
}
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- grabWorldAdmin(GeoGenCorr_wd = GeoGenCorr_wd, fileExtentsion = ".gpkg", overwrite = FALSE)
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = HIWA_coordinates)
# Read in the TNC global ecoregion shapefile, which is used for calculating ecological coverage 
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_globalTNC/Terrestrial_Ecoregions.shp')))
# Shapefiles are by default a 'non-exportable' object, which means they must be processed before being
# exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
world_poly_clip_W <- wrap(world_poly_clip)
ecoregion_poly_W <- wrap(ecoregion_poly)

# ---- GENETIC MATRIX
# Read in the VCF file containing Hibiscus individuals
HIWA_vcf <- read.vcfR(paste0(HIWA_filePath,'Genetic/hawaii.populations.snps.vcf'))
# Convert the vcf to a genind; the return.alleles FALSE value allows for downstream genetic distance calculations
HIWA_all_genind <- vcfR2genind(HIWA_vcf, return.alleles = FALSE)
# Using sample names, create a new genind object that doesn't include the ex situ samples
HIWA_genind <- HIWA_all_genind[HIWA_coordinates[,1], drop=TRUE]

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('HIWA_coordinates','HIWA_genind', 'num_reps','geo_buffSize', 
                              'eco_buffSize', 'world_poly_clip_W', 'ecoregion_poly_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres', 
                              'eco.intersectBuff','eco.compareBuff','gen.getAlleleCategories', 
                              'gen.buildDistMat', 'gen.calcGenDistCov', 'eco.totalEcoregionCount',
                              'calculateCoverage','exSituResample.Par', 'geo.gen.Resample.Par'))
# Specify file path, for saving resampling array
arrayDir <- paste0(HIWA_filePath, 'resamplingData/HIWA_SMBO3_G2GE_5r_resampArr.Rdata')
# Run resampling (in parallel)
HIWA_demoArray_Par <- 
  geo.gen.Resample.Par(genObj=HIWA_genind,  genDistFlag=TRUE, geoFlag=TRUE, coordPts=HIWA_coordinates, 
                       geoBuff=geo_buffSize, SDMrast=NA, boundary=world_poly_clip_W, 
                       ecoFlag=TRUE, ecoBuff=eco_buffSize, ecoRegions=ecoregion_poly_W, 
                       ecoLayer='GL', reps=num_reps, arrayFilepath=arrayDir, cluster=cl)
# Close cores
stopCluster(cl)

# Run resampling not in parallel (for function testing purposes)
# HIWA_demoArray_IND <-
#   geo.gen.Resample(genObj=HIWA_genind, geoFlag=TRUE, genDistFlag=TRUE, coordPts=HIWA_coordinates,
#                    geoBuff=geo_buffSize, boundary=world_poly_clip, ecoFlag=FALSE, reps = 1)

# # %%% ANALYZE DATA %%% ----
# # Specify filepath for HIWA geographic and genetic data, including resampling array
# HIWA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/HIWA/')
# arrayDir <- paste0(HIWA_filePath, 'resamplingData/HIWA_50km_GE_5r_resampArr.Rdata')
# # Read in the resampling array .Rdata object, saved to disk
# HIWA_demoArray_Par <- readRDS(arrayDir)
# 
# # %%%% SMBO: MULTIPLE BUFFER SIZES ----
# # Specify filepath for HIWA geographic and genetic data, including resampling array
# HIWA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/HIWA/')
# arrayDir <- paste0(HIWA_filePath, 'resamplingData/HIWA_SMBO3_G2GE_5r_resampArr.Rdata')
# # Read in array and build a data.frame of values
# HIWA_SMBO3_array <- readRDS(arrayDir)
# 
# # ---- CALCULATIONS ----
# # Build a data.frame from array values
# HIWA_SMBO3_DF <- resample.array2dataframe(HIWA_SMBO3_array)
# # Build tables of NRSMSE values, calculated based on data.frame
# HIWA_NRMSE_Mat_CV <- buildNRMSEmatrix(resampDF=HIWA_SMBO3_DF, genCovType='CV', sdmFlag=FALSE)
# HIWA_NRMSE_Mat_GD <- buildNRMSEmatrix(resampDF=HIWA_SMBO3_DF, genCovType='GD', sdmFlag=FALSE)
# # Combine the results of the NRMSE values calculated using allelic coverages and using
# # genetic distances, and then rename the columns accordingly
# HIWA_NRMSE_Mat <- cbind(HIWA_NRMSE_Mat_CV, HIWA_NRMSE_Mat_GD)
# # Store the matrix as a CSV to disk
# write.table(HIWA_NRMSE_Mat,
#             file=paste0(HIWA_filePath, 'resamplingData/HIWA_SMBO3_NRMSE.csv'), sep=',')
# 
# # SMBO2: OPTIMAL BUFFER SIZES ----
# # Read in HIWA SMBO2 resampling array amd convert to data.frame
# HIWA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/HIWA/')
# HIWA_arrayDir <- paste0(HIWA_filePath, 'resamplingData/HIWA_SMBO2_GE_5r_resampArr.Rdata')
# # From HIWA resampling array, return a matrix of average coverage values for optimal buffer sizes
# HIWA_optCovMat <- extractOptCovs(HIWA_arrayDir)
# # Calculate MSSEs: minimum number of samples for 95% of each coverage type
# HIWA_Gen_MSSE <- min(which(HIWA_optCovMat[,1] > 95)) ; HIWA_Gen_MSSE
# HIWA_GeoBuff_MSSE <- min(which(HIWA_optCovMat[,2] > 95)) ; HIWA_GeoBuff_MSSE
# HIWA_Eco_MSSE <- min(which(HIWA_optCovMat[,3] > 95)) ; HIWA_Eco_MSSE
