# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: MIMULUS GUTTATUS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage
# for Mimulus guttatus. Uses data files from Vallejo-Martin et al. 2021--a VCF, for genetic data, 
# and a spreadsheet that was derived from supplemental data from the manuscript. This dataset is 
# unique in that it's global, and coordinates are only provided at the population-level, 
# not at the individual-level. This script filters out introduced (non-native) populations, and
# conducts resampling analyses only using native populations (from western North America).

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
clusterExport(cl, varlist = c('MIGU_coordinates','MIGU_genind','num_reps','geo_buffSize', 'eco_buffSize',
                              'world_poly_clip_W', 'ecoregion_poly_W', 'MIGU_sdm_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers', 'geo.compareBuff', 'geo.compareBuffSDM', 
                              'eco.intersectBuff', 'eco.compareBuff', 'gen.getAlleleCategories',
                              'calculateCoverage', 'exSituResample.Par', 'geo.gen.Resample.Par'))
# Specify file path, for saving resampling array
arrayDir <- paste0(MIGU_filePath, 'resamplingData/MIGU_50km_GE_5r_resampArr.Rdata')

# Run resampling (in parallel)
print("%%% BEGINNING RESAMPLING %%%")
MIGU_demoArray_Par <- 
  geo.gen.Resample.Par(gen_obj = MIGU_genind, geoFlag = TRUE, coordPts = MIGU_coordinates, 
                       geoBuff = geo_buffSize, SDMrast = MIGU_sdm_W, 
                       boundary=world_poly_clip_W, ecoFlag = TRUE, ecoBuff = eco_buffSize, 
                       ecoRegions = ecoregion_poly_W, ecoLayer = 'NA', reps = num_reps, 
                       arrayFilepath = arrayDir, cluster = cl)
# Close cores
stopCluster(cl)

# %%% ANALYZE DATA %%% ----
# Specify filepath for MIGU geographic and genetic data, including resampling array
MIGU_filePath <- paste0(GeoGenCorr_wd, 'Datasets/MIGU/')
arrayDir <- paste0(MIGU_filePath, 'resamplingData/MIGU_50km_GE_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
MIGU_demoArray_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values, to pass to linear models
MIGU_DF <- resample.array2dataframe(MIGU_demoArray_Par)
# # Calculate Spearman's r for geographic/ecological coverage
# MIGU_spearR_geo <- round(cor(MIGU_DF$Geo, MIGU_DF$Total, method = 'spearman'),3) ; MIGU_spearR_geo
# MIGU_spearR_eco <- round(cor(MIGU_DF$Eco, MIGU_DF$Total, method = 'spearman'),3) ; MIGU_spearR_eco
# Calculate normalized root mean square value
MIGU_nrmse_geo <- nrmse_func(obs=MIGU_DF$Geo, pred=MIGU_DF$Total) ; MIGU_nrmse_geo
MIGU_nrmse_eco <- nrmse_func(obs=MIGU_DF$Eco, pred=MIGU_DF$Total) ; MIGU_nrmse_eco

# ---- PLOTTING ----
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(MIGU_demoArray_Par, allValues = TRUE)
# Subset matrix of all average values to just Total allelic, geographic, and ecological coverage
averageValueMat_TEG <- averageValueMat[,c(1,6,7)]
# Calculate the absolute difference between genetic and geographic/ecological, and add to data.frame
averageValueMat_TEG <- cbind(averageValueMat_TEG, abs(averageValueMat_TEG$Total-averageValueMat_TEG$Geo))
averageValueMat_TEG <- cbind(averageValueMat_TEG, abs(averageValueMat_TEG$Total-averageValueMat_TEG$Eco))
names(averageValueMat_TEG) <- c(names(averageValueMat_TEG)[1:3], 'Geo_Difference', 'Eco_Difference')

# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','darkblue', 'purple')
plotColors_Fade <- alpha(plotColors, 0.45)
plotColors_Sub <- plotColors_Fade[-(2:4)]
# Two plots in a single window
par(mfrow=c(2,1))
# ---- CORRELATION PLOTS
plot(averageValueMat_TEG$Geo, averageValueMat_TEG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     main='M. guttatus: Geographic by genetic coverage',xlab='', ylab='', col=plotColors_Fade[[5]])
mtext(text='255 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Geographic/Ecological coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add points for ecological coverage
points(x=averageValueMat$Eco, y=averageValueMat$Total, pch=20, col=plotColors_Fade[[6]])
# Add NRMSE values for each comparison
text(x = 76, y = 35, labels = paste0('NRMSE: ', MIGU_nrmse_geo), col='darkblue', cex=0.9)
text(x = 76, y = 20, labels = paste0('NRMSE: ', MIGU_nrmse_eco), col='purple', cex=0.9)
# Add legend
legend(x=58, y=54, inset = 0.05, xpd=TRUE,
       legend = c('Geographic', 'Ecological'),
       col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n', 
       y.intersp = 0.8)
# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TEG[,1:3], ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='M. guttatus: Geo-Eco-Gen Coverage', line=1.5)
mtext(text='255 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=160, y=68, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (50 km buffer)', 
                  'Ecological coverage (50 km buffer, EPA Level III)'),
       col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 0.8)
# ---- DIFFERENCE PLOTS
# Plot difference between geographic and genetic coverage
matplot(averageValueMat_TEG[4:5], col=plotColors_Fade[5:6], pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='M. guttatus: Genetic-Geographic-Ecological Coverage Difference', line=1.5)
mtext(text='255 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.2)
mtext(text='Difference in Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=160, y=50, inset = 0.05,
       legend = c('Genographic coverage difference', 'Ecological coverage difference'), 
       col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 1)

# %%%% SDM AND TOTAL BUFFER COMPARISON ----
# Specify filepath for MIGU geographic and genetic data, including resampling array
MIGU_filePath <- paste0(GeoGenCorr_wd, 'Datasets/MIGU/')
arrayDir <- paste0(MIGU_filePath, 'resamplingData/MIGU_50km_G2E_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
MIGU_geoComp_50km_array <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values, to pass to linear models
MIGU_geoComp_50km_DF <- resample.array2dataframe(MIGU_geoComp_50km_array)
# # Calculate Spearman's R for total buffer/SDM coverage
# MIGU_spearR_geo_totalBuff <- 
#   round(cor(MIGU_geoComp_50km_DF$Geo_Buff, MIGU_geoComp_50km_DF$Total, method = 'spearman'),3)
# MIGU_spearR_geo_SDM <- 
#   round(cor(MIGU_geoComp_50km_DF$Geo_SDM, MIGU_geoComp_50km_DF$Total, method = 'spearman'),3)
# Calculate normalized root mean square value
MIGU_nrmse_geo_totalBuff <- 
  nrmse_func(obs=MIGU_geoComp_50km_DF$Geo_Buff, pred=MIGU_geoComp_50km_DF$Total) ; MIGU_nrmse_geo_totalBuff
MIGU_nrmse_geo_SDM <- 
  nrmse_func(obs=MIGU_geoComp_50km_DF$Geo_SDM, pred=MIGU_geoComp_50km_DF$Total) ; MIGU_nrmse_geo_SDM

# ---- PLOTTING
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic and geographic proportions
MIGU_geoComp_50km_averageValueMat <- meanArrayValues(MIGU_geoComp_50km_array)
# Calculate the absolute difference between genetic and geographic approaches, and add to data.frame
MIGU_geoComp_50km_averageValueMat <- 
  cbind(MIGU_geoComp_50km_averageValueMat, abs(MIGU_geoComp_50km_averageValueMat$Total-MIGU_geoComp_50km_averageValueMat$Geo_Buff))
MIGU_geoComp_50km_averageValueMat <- 
  cbind(MIGU_geoComp_50km_averageValueMat, abs(MIGU_geoComp_50km_averageValueMat$Total-MIGU_geoComp_50km_averageValueMat$Geo_SDM))
names(MIGU_geoComp_50km_averageValueMat) <- 
  c(names(MIGU_geoComp_50km_averageValueMat)[1:4],'Geo_Buff_Difference', 'Geo_SDM_Difference')
# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue')
plotColors_fade <- alpha(c('red','red4','darkorange3','coral','purple', 'darkblue'), 0.45)
# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1))

# ---- CORRELATION PLOT
# Plot genetic coverage against geographic coverage using total buffer approach
plot(MIGU_geoComp_50km_averageValueMat$Geo_Buff, MIGU_geoComp_50km_averageValueMat$Total, pch=20,
     main='M. guttatus: Geographic by genetic coverage',
     xlab='Geographic coverage (%)', ylab='Genetic coverage (%)',
     col=plotColors_fade[[2]])
# Add points for SDM approach
points(x=MIGU_geoComp_50km_averageValueMat$Geo_SDM, y=MIGU_geoComp_50km_averageValueMat$Total,
       pch=20, col=plotColors_fade[[3]])
# Subtitle
mtext(text='255 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3)
# Add NRMSE values for each comparison
text(x = 79, y = 71, labels = paste0('NRMSE: ', MIGU_nrmse_geo_totalBuff), col='red4', cex=0.9)
text(x = 79, y = 66, labels = paste0('NRMSE: ', MIGU_nrmse_geo_SDM), col='darkorange3', cex=0.9)
# Add legend
legend(x=57, y=78, inset = 0.05, xpd=TRUE,
       legend = c('Total buffer approach', 'SDM approach'),
       col=plotColors[2:3], pch = c(20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.8)
# ---- COVERAGE PLOT
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(MIGU_geoComp_50km_averageValueMat[,1:3], ylim=c(0,100), col=plotColors_fade, 
        pch=16, ylab='Coverage (%)')
# Add title and x-axis labels to the graph
title(main='M. guttatus: Coverage Values by Sample Size', line=1.5)
mtext(text='255 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3)
mtext(text='Number of individuals', side=1, line=2.4)
# Add legend
legend(x=175, y=50, inset = 0.05, xpd=TRUE,
       legend = c('Genetic coverage', 'Geographic, Total buffer (50 km)', 'Geographic, SDM (50 km)'),
       col=plotColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.8)
# ---- DIFFERENCE PLOTS
# Plot difference between geographic and genetic coverage
matplot(MIGU_geoComp_50km_averageValueMat[5:6], col=plotColors_fade[2:3], pch=16, ylab='')
# Add title, subtitle, and x-axis labels to the graph
title(main='M. guttatus: Genetic-Geographic Coverage Difference', line=1.5)
mtext(text='255 Individuals; 50 km buffer; 5 replicates', side=3, line=0.1)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.2)
mtext(text='Difference in Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=160, y=45, inset = 0.05,
       legend = c('Total buffer approach', 'SDM approach'), 
       col=c('red4', 'darkorange3'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 1)
