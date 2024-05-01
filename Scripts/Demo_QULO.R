# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: QUERCUS LOBATA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script demonstrating approach for calculating the correlation 
# between genetic, geographic, and ecological coverage. 
# Uses data files from Gugger et al. 2021 to pull in genetic data 
# (as a table) and geographic coordinates (in a CSV) to conduct correlation analyses.

library(adegenet)
library(terra)
library(parallel)
library(RColorBrewer)
library(scales)

# Read in relevant functions
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/Code/'
source('Scripts/functions_GeoGenCoverage.R')
source('Scripts/worldAdmin.R')

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
# Specify filepath for QULO geographic and genetic data
QULO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QULO/')

# ---- GENETIC MATRIX
# Read in the SNPs genetic matrix file
QULO_tab <- read.table(paste0(QULO_filePath, 'Genetic/SNPs80.forR'), header = TRUE)
# Make sample names the row names
rownames(QULO_tab) <- QULO_tab[,1] ; QULO_tab <- QULO_tab[,-1]
# Convert to genind. ncode based off similar practice for other geninds 
# (values of 1, 2, and 3 generate identical results)
QULO_genind <- df2genind(QULO_tab, ncode = 3, ploidy = 2)

# ---- GEOGRAPHIC/ECOLOGICAL DATA FILES
# Read in wild occurrence points. This CSV has 3 columns: sample name, latitude, and longitude. 
# The sample names (and order) have to match the sample names/order of the genind object 
# (rownams of the genetic matrix) read in below.
QULO_points <- read.csv(paste0(QULO_filePath, 'Geographic/Quercus_lobata.csv'), header=TRUE)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- grabWorldAdmin(GeoGenCorr_wd = GeoGenCorr_wd, fileExtentsion = ".gpkg", overwrite = FALSE)
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = QULO_points)
# Read in raster data, for SDM
QULO_sdm <- terra::rast(paste0(QULO_filePath,'Geographic/QULO_436inds_rast_Loza_MX.tif'))
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage 
# (solely in the U.S.)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level4/us_eco_l4.shp')))
# Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
# exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
world_poly_clip_W <- wrap(world_poly_clip)
ecoregion_poly_W <- wrap(ecoregion_poly)
QULO_sdm_W <- wrap(QULO_sdm)

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('QULO_points','QULO_genind','QULO_sdm', 'num_reps','geo_buffSize', 
                              'eco_buffSize', 'world_poly_clip_W', 'ecoregion_poly_W', 'QULO_sdm_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers', 'geo.compareBuff', 'geo.compareBuffSDM', 
                              'eco.intersectBuff', 'eco.compareBuff', 'gen.getAlleleCategories',
                              'calculateCoverage', 'exSituResample.Par', 'geo.gen.Resample.Par'))
# Specify file path, for saving resampling array
arrayDir <- paste0(QULO_filePath, 'resamplingData/QULO_50km_GE_5r_resampArr.Rdata')
# Run resampling (in parallel)
QULO_demoArray_Par <- 
  geo.gen.Resample.Par(gen_obj = QULO_genind, geoFlag = TRUE, coordPts = QULO_points, 
                       geoBuff = geo_buffSize, SDMrast=QULO_sdm_W, boundary=world_poly_clip_W, 
                       ecoFlag = TRUE, ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, 
                       ecoLayer = 'US', reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
# Close cores
stopCluster(cl)

# Run resampling not in parallel (for function testing purposes)
# QULO_demoArray_IND <-
#   geo.gen.Resample(gen_obj = QULO_genind, SDMrast = QULO_sdm, geoFlag = TRUE, coordPts = QULO_points,
#                    geoBuff = geo_buffSize, boundary = world_poly_clip, ecoFlag = FALSE, reps = 1)

# %%% ANALYZE DATA %%% ----
# Specify filepath for QULO geographic and genetic data, including resampling array
QULO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QULO/')
arrayDir <- paste0(QULO_filePath, 'resamplingData/QULO_50km_GE_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
QULO_demoArray_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values, to pass to linear models
QULO_DF <- resample.array2dataframe(QULO_demoArray_Par)
# Calculate Spearman's r for geographic/ecological coverage
QULO_spearR_geo <- round(cor(QULO_DF$Geo, QULO_DF$Total, method = 'spearman'),3) ; QULO_spearR_geo
QULO_spearR_eco <- round(cor(QULO_DF$Eco, QULO_DF$Total, method = 'spearman'),3) ; QULO_spearR_eco

# ---- LINEAR MODELS
# Generate linear models, using Total allelic coverage as the response variable
# GEOGRAPHIC COVERAGE AS PREDICTOR VARIABLE
QULO_geoModel <- lm (Total ~ Geo, data=QULO_DF)
QULO_geoModel_summary <- summary(QULO_geoModel) ; QULO_geoModel_summary
# Pull R-squared estimate from model
QULO_geoModel_rSquared <- round(QULO_geoModel_summary$adj.r.squared,2) ; QULO_geoModel_rSquared
# ECOLOGICAL COVERAGE AS PREDICTOR VARIABLE
QULO_ecoModel <- lm (Total ~ Eco, data=QULO_DF)
QULO_ecoModel_summary <- summary(QULO_ecoModel) ; QULO_ecoModel_summary
# Pull R-squared estimate from model
QULO_ecoModel_rSquared <- round(QULO_ecoModel_summary$adj.r.squared, 2) ; QULO_ecoModel_rSquared

# ---- PLOTTING ----
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(QULO_demoArray_Par, allValues = TRUE)
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
# Add Spearman's r values for each comparison
text(x = 76, y = 35, labels = paste0('Spearman r: ', QULO_spearR_geo), col='darkblue', cex=0.9)
text(x = 76, y = 20, labels = paste0('Spearman r: ', QULO_spearR_eco), col='purple', cex=0.9)
# Add legend
legend(x=60, y=235, inset = 0.05, xpd=TRUE,
       legend = c('Geographic', 'Ecological'),
       col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.04)
# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TEG, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='Q. lobata: Geo-Eco-Gen Coverage', line=1.5)
mtext(text='436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=275, y=220, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (50 km buffer)', 
                  'Ecological coverage (50 km buffer, EPA Level IV)'),
       col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 0.03)
# ---- DIFFERENCE PLOTS
# Plot difference between geographic and genetic coverage
matplot(averageValueMat_TEG[4:5], col=plotColors[5:6], pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='Q. lobata: Genetic-Geographic-Ecological Coverage Difference', line=1.5)
mtext(text='436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.2)
mtext(text='Difference in Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=275, y=45, inset = 0.05,
       legend = c('Genographic coverage difference', 'Ecological coverage difference'), 
       col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 1)

# %%%% 2024-03-11 SDM AND TOTAL BUFFER COMPARISON ----
# Specify filepath for QULO geographic and genetic data, including resampling array
QULO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QULO/')
arrayDir <- paste0(QULO_filePath, 'resamplingData/QULO_50km_G2E_Carver_5r_resampArr.Rdata')
# Read in array and build a data.frame of values
QULO_geoComp_50km_array <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values, to pass to linear models
QULO_geoComp_50km_DF <- resample.array2dataframe(QULO_geoComp_50km_array)
# Calculate Spearman's r for total buffer/SDM coverage
QULO_spearR_geo_totalBuff <- 
  round(cor(QULO_geoComp_50km_DF$Geo_Buff, QULO_geoComp_50km_DF$Total, method = 'spearman'),3)
QULO_spearR_geo_SDM <- 
  round(cor(QULO_geoComp_50km_DF$Geo_SDM, QULO_geoComp_50km_DF$Total, method = 'spearman'),3)

# # ---- LINEAR MODELS
# # Generate linear models, using Total allelic coverage as the response variable
# # Use either the total buffer (Buff) or SDM (SDM) geographic coverage approach for the predictor variable
# # Also extract R squared values
# # Total buffer approach
# QULO_geoComp_50km_geoModelBuff <- lm (Total ~ Geo_Buff, data=QULO_geoComp_50km_DF)
# QULO_geoComp_50km_geoModelBuff_summary <- summary(QULO_geoComp_50km_geoModelBuff) ; QULO_geoComp_50km_geoModelBuff_summary
# QULO_geoComp_50km_geoModelBuff_rSquared <- round(QULO_geoComp_50km_geoModelBuff_summary$adj.r.squared,2)
# # SDM approach
# QULO_geoComp_50km_geoModelSDM <- lm (Total ~ Geo_SDM, data=QULO_geoComp_50km_DF)
# QULO_geoComp_50km_geoModelSDM_summary <- summary(QULO_geoComp_50km_geoModelSDM) ; QULO_geoComp_50km_geoModelSDM_summary
# QULO_geoComp_50km_geoModelSDM_rSquared <- round(QULO_geoComp_50km_geoModelSDM_summary$adj.r.squared,2)

# ---- PLOTTING
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic and geographic proportions
QULO_geoComp_50km_averageValueMat <- meanArrayValues(QULO_geoComp_50km_array)
# Calculate the absolute difference between genetic and geographic approaches, and add to data.frame
QULO_geoComp_50km_averageValueMat <- 
  cbind(QULO_geoComp_50km_averageValueMat, abs(QULO_geoComp_50km_averageValueMat$Total-QULO_geoComp_50km_averageValueMat$Geo_Buff))
QULO_geoComp_50km_averageValueMat <- 
  cbind(QULO_geoComp_50km_averageValueMat, abs(QULO_geoComp_50km_averageValueMat$Total-QULO_geoComp_50km_averageValueMat$Geo_SDM))
names(QULO_geoComp_50km_averageValueMat) <- 
  c(names(QULO_geoComp_50km_averageValueMat)[1:4],'Geo_Buff_Difference', 'Geo_SDM_Difference')
# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue')
plotColors <- alpha(plotColors, 0.45)

# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1))

# ---- CORRELATION PLOT
# Plot genetic coverage against geographic coverage using total buffer approach
plot(QULO_geoComp_50km_averageValueMat$Geo_Buff, QULO_geoComp_50km_averageValueMat$Total, pch=20,
     main='Q. lobata: Geographic by genetic coverage',
     xlab='Geographic coverage (%)', ylab='Genetic coverage (%)',
     col='red4')
# Add points for SDM approach
points(x=QULO_geoComp_50km_averageValueMat$Geo_SDM, y=QULO_geoComp_50km_averageValueMat$Total,
       pch=20, col='darkorange3')
# Subtitle
mtext(text='SDM: Carver; 436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.1)
# Add Spearman's r values for each comparison
text(x = 80, y = 62, labels = paste0('Spearman r: ', QULO_spearR_geo_totalBuff), col='red4', cex=0.9)
text(x = 80, y = 57, labels = paste0('Spearman r: ', QULO_spearR_geo_SDM), col='darkorange3', cex=0.9)
# Add legend
legend(x=60, y=155, inset = 0.05, xpd=TRUE,
       legend = c('Total buffer approach', 'SDM approach'),
       col=plotColors[2:3], pch = c(20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.04)
# ---- COVERAGE PLOT
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(QULO_geoComp_50km_averageValueMat[,1:3], ylim=c(0,100), col=plotColors, pch=16, ylab='Coverage (%)')
# Add title and x-axis labels to the graph
title(main='Q. lobata: Coverage Values by Sample Size', line=1.5)
mtext(text='SDM: Carver; 436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3)
mtext(text='Number of individuals', side=1, line=2.4)
# Add legend
legend(x=300, y=225, inset = 0.05, xpd=TRUE,
       legend = c('Genetic coverage', 'Geographic, Total buffer (50 km)', 'Geographic, SDM (50 km)'),
       col=plotColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.04)
# ---- DIFFERENCE PLOTS
# Plot difference between geographic and genetic coverage
matplot(QULO_geoComp_50km_averageValueMat[5:6], col=plotColors[2:3], pch=16, ylab='')
# Add title, subtitle, and x-axis labels to the graph
title(main='Q. lobata: Genetic-Geographic Coverage Difference', line=1.5)
mtext(text='SDM: Carver; 436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.1)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.2)
mtext(text='Difference in Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=275, y=30, inset = 0.05,
       legend = c('Total buffer approach', 'SDM approach'), 
       col=c('red4', 'darkorange3'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 1)
# # ---- MAP OF SDM AND SAMPLED POINTS
# makeAMap(QULO_points, raster = QULO_sdm, buffer = geo_buffSize)
