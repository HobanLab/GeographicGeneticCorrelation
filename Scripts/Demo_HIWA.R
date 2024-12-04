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
# Make sure libraries (adegenet + terra) are on cluster
clusterEvalQ(cl, library('adegenet'))
clusterEvalQ(cl, library('terra'))
clusterEvalQ(cl, library('parallel'))
clusterEvalQ(cl, library('usedist'))

# %%% CONDUCT RESAMPLING %%% ----
# ---- READ IN DATA ----
# Specify filepath for HIWA geographic and genetic data
HIWA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/HIWA/')

# ---- GEOGRAPHIC/ECOLOGICAL DATA FILES
# Read in wild occurrence points. This CSV has 5 columns: sample name, sample code, locality, 
# latitude, and longitude
HIWA_points <- read.csv(paste0(HIWA_filePath, 'Geographic/HIWA_coordinates.csv'), header=TRUE)
# Remove the ex situ samples from the dataframe, and then drop sample code and locality columns
HIWA_points <- HIWA_points[-which(HIWA_points$Locality=='Ex Situ'),]
HIWA_points <- HIWA_points[,-(2:3)]
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- grabWorldAdmin(GeoGenCorr_wd = GeoGenCorr_wd, fileExtentsion = ".gpkg", overwrite = FALSE)
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = HIWA_points)
# Read in the TNC global ecoregion shapefile, which is used for calculating ecological coverage 
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_globalTNC/Terrestrial_Ecoregions.shp')))
# Shapefiles are by default a 'non-exportable' object, which means they must be processed before being
# exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
world_poly_clip_W <- wrap(world_poly_clip)
ecoregion_poly_W <- wrap(ecoregion_poly)

# ---- GENETIC MATRIX
# # Read in the genepop file provided for this dataset
# HIWA_all_genind <- read.genepop(paste0(HIWA_filePath,'Genetic/final.recode.p.snps.gen'))
# # Using sample names, create a new genind object that doesn't include the ex situ samples
# HIWA_genind <- HIWA_all_genind[HIWA_points[,1], drop=TRUE]
# 
# Read in the VCF file containing Hibiscus individuals
HIWA_vcf <- read.vcfR(paste0(HIWA_filePath,'Genetic/hawaii.populations.snps.vcf'))
# Convert the vcf to a genind; the return.alleles FALSE value allows for downstream genetic distance calculations
HIWA_all_genind <- vcfR2genind(HIWA_vcf, return.alleles = FALSE)
# Using sample names, create a new genind object that doesn't include the ex situ samples
HIWA_genind <- HIWA_all_genind[HIWA_points[,1], drop=TRUE]

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('HIWA_points','HIWA_genind', 'num_reps','geo_buffSize', 
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
  geo.gen.Resample.Par(genObj=HIWA_genind,  genDistFlag=TRUE, geoFlag=TRUE, coordPts=HIWA_points, 
                       geoBuff=geo_buffSize, SDMrast=NA, boundary=world_poly_clip_W, 
                       ecoFlag=TRUE, ecoBuff=eco_buffSize, ecoRegions=ecoregion_poly_W, 
                       ecoLayer='GL', reps=num_reps, arrayFilepath=arrayDir, cluster=cl)
# Close cores
stopCluster(cl)

# Run resampling not in parallel (for function testing purposes)
# HIWA_demoArray_IND <-
#   geo.gen.Resample(genObj=HIWA_genind, geoFlag=TRUE, genDistFlag=TRUE, coordPts=HIWA_points,
#                    geoBuff=geo_buffSize, boundary=world_poly_clip, ecoFlag=FALSE, reps = 1)

# %%% ANALYZE DATA %%% ----
# Specify filepath for HIWA geographic and genetic data, including resampling array
HIWA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/HIWA/')
arrayDir <- paste0(HIWA_filePath, 'resamplingData/HIWA_50km_GE_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
HIWA_demoArray_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values, to pass to linear models
HIWA_DF <- resample.array2dataframe(HIWA_demoArray_Par)
# Calculate normalized root mean square value
HIWA_nrmse_geo <- nrmse_func(obs=HIWA_DF$Geo, pred=HIWA_DF$Total) ; HIWA_nrmse_geo
HIWA_nrmse_eco <- nrmse_func(obs=HIWA_DF$Eco, pred=HIWA_DF$Total) ; HIWA_nrmse_eco

# ---- PLOTTING ----
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(HIWA_demoArray_Par, allValues = TRUE)
# Subset matrix of all average values to just Total allelic, geographic, and ecological coverage
averageValueMat_TEG <- averageValueMat[,c(1,6,7)]
# Calculate the absolute difference between genetic and geographic/ecological, and add to data.frame
averageValueMat_TEG <- cbind(averageValueMat_TEG, abs(averageValueMat_TEG$Total-averageValueMat_TEG$Geo))
averageValueMat_TEG <- cbind(averageValueMat_TEG, abs(averageValueMat_TEG$Total-averageValueMat_TEG$Eco))
names(averageValueMat_TEG) <- c(names(averageValueMat_TEG)[1:3], 'Geo_Difference', 'Eco_Difference')
# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','darkblue', 'purple')
plotColors <- alpha(plotColors, 0.45)
plotColors_Sub <- plotColors[-(2:4)]
# Two plots in a single window
par(mfrow=c(2,1))

# ---- CORRELATION PLOTS
plot(averageValueMat_TEG$Geo, averageValueMat_TEG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     main='Q. lobata: Geographic/Ecological by genetic coverage',xlab='', ylab='', col=plotColors[[5]])
mtext(text='436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Geographic/Ecological coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add points for ecological coverage
points(x=averageValueMat$Eco, y=averageValueMat$Total, pch=20, col=plotColors[[6]])
# Add NRMSE values for each comparison
text(x = 76, y = 35, labels = paste0('NRMSE: ', HIWA_nrmse_geo), col='darkblue', cex=0.9)
text(x = 76, y = 20, labels = paste0('NRMSE: ', HIWA_nrmse_eco), col='purple', cex=0.9)
# Add legend
legend(x=57, y=53, inset = 0.05, xpd=TRUE,
       legend = c('Geographic', 'Ecological'),
       col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.8)
# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TEG[,1:3], ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='Q. lobata: Geo-Eco-Gen Coverage', line=1.5)
mtext(text='436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=275, y=60, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (50 km buffer)',
                  'Ecological coverage (50 km buffer, EPA Level IV)'),
       col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 0.8)

# %%%% SMBO: MULTIPLE BUFFER SIZES ----
# Specify filepath for HIWA geographic and genetic data, including resampling array
HIWA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/HIWA/')
arrayDir <- paste0(HIWA_filePath, 'resamplingData/HIWA_SMBO3_G2GE_5r_resampArr.Rdata')
# Read in array and build a data.frame of values
HIWA_SMBO3_array <- readRDS(arrayDir)

# ---- CALCULATIONS ----
# Build a data.frame from array values
HIWA_SMBO3_DF <- resample.array2dataframe(HIWA_SMBO3_array)
# Build tables of NRSMSE values, calculated based on data.frame
HIWA_NRMSE_Mat_CV <- buildNRMSEmatrix(resampDF=HIWA_SMBO3_DF, genCovType='CV', sdmFlag=FALSE)
HIWA_NRMSE_Mat_GD <- buildNRMSEmatrix(resampDF=HIWA_SMBO3_DF, genCovType='GD', sdmFlag=FALSE)
# Combine the results of the NRMSE values calculated using allelic coverages and using
# genetic distances, and then rename the columns accordingly
HIWA_NRMSE_Mat <- cbind(HIWA_NRMSE_Mat_CV, HIWA_NRMSE_Mat_GD)
# Store the matrix as a CSV to disk
write.table(HIWA_NRMSE_Mat,
            file=paste0(HIWA_filePath, 'resamplingData/HIWA_SMBO3_NRMSE.csv'), sep=',')
