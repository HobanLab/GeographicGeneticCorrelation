# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: QUERCUS LOBATA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script demonstrating first draft approach for calculating the correlation 
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
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.shp')))
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage (solely in the U.S.)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level4/us_eco_l4.shp')))
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
# Specify filepath for QULO geographic and genetic data
QULO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QULO/')

# ---- GENETIC MATRIX
# Processing input file
QULO_tab <- read.table(paste0(QULO_filePath, 'Genetic/SNPs80.forR'), header = TRUE)
# Make sample names the row names
rownames(QULO_tab) <- QULO_tab[,1] ; QULO_tab <- QULO_tab[,-1]
# Convert to genind. ncode based off similar practice for other geninds (values of 1, 2, and 3 generate identical results)
QULO_genind <- df2genind(QULO_tab, ncode = 3, ploidy = 2)

# ---- GEOGRAPHIC COORDINATES
# Read in wild occurrence points. This CSV has 3 columns: sample name, latitude, and longitude. 
# The sample names (and order) have to match the sample names/order of the genind object 
# (rownams of the genetic matrix) read in below.
wildPoints <- read.csv(paste0(QULO_filePath, 'Geographic/Quercus_lobata.csv'), header=TRUE)

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('wildPoints','QULO_genind','num_reps','geo_buffSize', 'eco_buffSize',
                              'world_poly_clip_W', 'ecoregion_poly_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers', 'geo.compareBuff', 'eco.intersectBuff', 'eco.compareBuff',
                              'gen.getAlleleCategories','calculateCoverage', 'exSituResample', 
                              'geo.gen.Resample.Parallel'))
# Specify file path, for saving resampling array
# arrayDir <- paste0(QUAC.filePath, 'resamplingData/QULO_1km_GE_5r_resampArr.Rdata')
arrayDir <- paste0(QULO_filePath, 'resamplingData/QULO_50km_GE_5r_resampArr.Rdata')
# Run resampling (in parallel)
QULO_demoArray_Par <- 
  geo.gen.Resample.Parallel(gen_obj = QULO_genind, geoFlag = TRUE, coordPts = wildPoints, 
                            geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag = TRUE, 
                            ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer = 'US', 
                            reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
# Close cores
stopCluster(cl)

# %%% ANALYZE DATA %%% ----
# Read in the resampling array .Rdata object, saved to disk
QULO_demoArray_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values, to pass to linear models
QULO_DF <- resample.array2dataframe(QULO_demoArray_Par)

# ---- LINEAR MODELS
# Generate linear models, using Total allelic coverage as the response variable
# GEOGRAPHIC COVERAGE AS PREDICTOR VARIABLE
QULO_geoModel <- lm (Total ~ Geo, data=QULO_DF)
QULO_geoModel_summary <- summary(QULO_geoModel) ; QULO_geoModel_summary
# Pull R-squared and p-value estimates from model
QULO_geoModel_rSquared <- round(QULO_geoModel_summary$adj.r.squared,2)
QULO_geoModel_pValue <- QULO_geoModel_summary$coefficients[2, 4]
# ECOLOGICAL COVERAGE AS PREDICTOR VARIABLE
QULO_ecoModel <- lm (Total ~ Eco, data=QULO_DF)
QULO_ecoModel_summary <- summary(QULO_ecoModel) ; QULO_ecoModel_summary
# Pull R-squared and p-value estimates from model
QULO_ecoModel_rSquared <- round(QULO_ecoModel_summary$adj.r.squared, 2)
QULO_ecoModel_pValue <- QULO_ecoModel_summary$coefficients[2, 4]

# ---- PLOTTING ----
# ---- CALCULATE 95% MSSE AND AVERAGE VALUES
# Calculate minimum 95% sample size for genetic and geographic values
gen_min95Value <- gen.min95Mean(QULO_demoArray_Par) ; gen_min95Value
geo_min95Value <- geo.min95Mean(QULO_demoArray_Par) ; geo_min95Value
eco_min95Value <- eco.min95Mean(QULO_demoArray_Par) ; eco_min95Value
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(QULO_demoArray_Par, allValues = TRUE)
# Subset matrix of all average values to just Total allelic, geographic, and ecological coverage
averageValueMat_TEG <- averageValueMat[,c(1,6,7)]

# Specify plot colors ()
plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue', 'purple')
plotColors <- alpha(plotColors, 0.45)
plotColors_Sub <- plotColors[-(2:5)]

# ---- CORRELATION PLOTS
par(mfrow=c(2,1))
# ---- GEOGRAPHIC-GENETIC
plot(averageValueMat_TEG$Geo, averageValueMat_TEG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     main='Q. lobata: Geographic by genetic coverage',xlab='', ylab='')
mtext(text='436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Geographic coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
mylabel = bquote(italic(R)^2 == .(format(QULO_geoModel_rSquared, digits = 3)))
text(x = 2, y = 10, labels = mylabel, cex=0.8)
# ---- ECOLOGICAL-GENETIC
plot(averageValueMat_TEG$Eco, averageValueMat_TEG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     main='Q. lobata: Ecological by genetic coverage',xlab='', ylab='')
mtext(text='436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Ecological coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
mylabel = bquote(italic(R)^2 == .(format(QULO_ecoModel_rSquared, digits = 3)))
text(x = 2, y = 10, labels = mylabel, cex=0.8)

# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TEG, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='Quercus lobata: Geo-Eco-Gen Coverage', line=1.5)
mtext(text='436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=205, y=60, inset = 0.05,
       legend = c('Genetic coverage (Total)', 'Geographic coverage (50 km buffer)', 'Ecological coverage (EPA Level IV)'),
       col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=1.2, pt.cex = 2, bty='n', 
       y.intersp = 0.8)

# ---- BOTH PLOTS (For IMLS NLG subgroup presentation, 2023-08-17) ----
par(mfrow=c(2,1))

# ---- TOTAL ALLELIC AND GEOGRAPHIC COVERAGE
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='Quercus lobata: Geo-Gen Coverage', line=1.5)
mtext(text='436 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Mark the 95% threshold line, and the genetic/geographic points
abline(h=95, col='black', lty=3)
# Add legend
legend(x=215, y=65, inset = 0.05,
       legend = c('Genetic coverage (Total)', 'Geographic coverage (50 km buffer)'),
       col=plotColors_Sub, pch = c(20,20,20), cex=1.2, pt.cex = 2, bty='n', y.intersp = 0.6)

# ---- GEOGRAPHIC-GENETIC CORRELATION
plot(averageValueMat$Geo, averageValueMat$Total, pch=20, main='',xlab='', ylab='')
mtext(text='Geographic coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
mylabel = bquote(italic(R)^2 == .(format(QULO_model_rSquared, digits = 3)))
text(x = 20, y = 85, labels = mylabel, cex=1.4)
