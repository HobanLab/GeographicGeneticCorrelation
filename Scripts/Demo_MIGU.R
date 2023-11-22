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

# %%% CONDUCT RESAMPLING %%% ----
# ---- READ IN DATA ----
# Specify filepath for MIGU geographic and genetic data
MIGU_filePath <- paste0(GeoGenCorr_wd, 'Datasets/MIGU/')

# ---- GENETIC MATRIX
# Read in the VCF file provided in Vallejo-Martin et al. 2021 using vcfR::read.vcfR
MIGU_vcf <- 
  read.vcfR(file=paste0(MIGU_filePath, 'Genetic/mgut_all_20180305_gut_filter_75.i50.recode.pruned.plink_20180326.vcf'))
# Convert the vcf to a genind; the return.alleles TRUE value is suggested in the function's help file
# This genind file is made up of 474 individuals and 1,498 loci
MIGU_genind_global <- vcfR2genind(MIGU_vcf, sep = "/", return.alleles = TRUE)

# ---- COORDINATE POINTS
# The supplement of Vallejo-Martin et al. 2021 includes a PDF of the latitudes and longitudes for each population 
# (SupplementaryData1). This was first converted to a Google Sheet. In order to convert this to a CSV that had latitudes 
# and longitudes for each individual, individual names were added to the spreadsheet and matched manually 
# to the lat/long values for each population. Finally, this CSV was subset to include only the 255 individuals
# present in the native range of western North America.
MIGU_coordinates <- 
  read.csv(file=paste0(MIGU_filePath, 'Geographic/MIGU_LatLongs.csv'), header = TRUE)
# Rename the columns of the geographic coordinates data.frame (because geo.compareBuff function expects certain strings)
colnames(MIGU_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')

# ---- REMOVE INTRODUCED POPULATIONS ----
# Subset global genind object to only contain individuals from native range. The 'drop' argument removes alleles
# no longer present in the dataset.
MIGU_genind <- MIGU_genind_global[MIGU_coordinates[,1], drop=TRUE]

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('MIGU_coordinates','MIGU_genind','num_reps','geo_buffSize', 'eco_buffSize',
                              'world_poly_clip_W', 'ecoregion_poly_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers', 'geo.compareBuff', 'eco.intersectBuff', 'eco.compareBuff',
                              'gen.getAlleleCategories','calculateCoverage', 'exSituResample', 
                              'geo.gen.Resample.Parallel'))
# Specify file path, for saving resampling array
arrayDir <- paste0(MIGU_filePath, 'resamplingData/MIGU_50km_GE_5r_resampArr.Rdata')

# Run resampling (in parallel)
print("%%% BEGINNING RESAMPLING %%%")
MIGU_demoArray_Par <- 
  geo.gen.Resample.Parallel(gen_obj = MIGU_genind, geoFlag = TRUE, coordPts = MIGU_coordinates, 
                            geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag = TRUE, 
                            ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer = 'NA', 
                            reps = num_reps, arrayFilepath = arrayDir, cluster = cl)

# Close cores
stopCluster(cl)

# %%% ANALYZE DATA %%% ----
# Read in the resampling array .Rdata object, saved to disk
MIGU_demoArray_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values, to pass to linear models
MIGU_DF <- resample.array2dataframe(MIGU_demoArray_Par)

# ---- LINEAR MODELS
# Generate linear models, using Total allelic coverage as the response variable
# GEOGRAPHIC COVERAGE AS PREDICTOR VARIABLE
MIGU_geoModel <- lm (Total ~ Geo, data=MIGU_DF)
MIGU_geoModel_summary <- summary(MIGU_geoModel) ; MIGU_geoModel_summary
# Pull R-squared and p-value estimates from model
MIGU_geoModel_rSquared <- round(MIGU_geoModel_summary$adj.r.squared,2)
MIGU_geoModel_pValue <- MIGU_geoModel_summary$coefficients[2, 4]
# ECOLOGICAL COVERAGE AS PREDICTOR VARIABLE
MIGU_ecoModel <- lm (Total ~ Eco, data=MIGU_DF)
MIGU_ecoModel_summary <- summary(MIGU_ecoModel) ; MIGU_ecoModel_summary
# Pull R-squared and p-value estimates from model
MIGU_ecoModel_rSquared <- round(MIGU_ecoModel_summary$adj.r.squared, 2)
MIGU_ecoModel_pValue <- MIGU_ecoModel_summary$coefficients[2, 4]

# ---- PLOTTING ----
# ---- CALCULATE 95% MSSE AND AVERAGE VALUES
# Calculate minimum 95% sample size for genetic and geographic values
gen_min95Value <- gen.min95Mean(MIGU_demoArray_Par) ; gen_min95Value
geo_min95Value <- geo.min95Mean(MIGU_demoArray_Par) ; geo_min95Value
eco_min95Value <- eco.min95Mean(MIGU_demoArray_Par) ; eco_min95Value
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(MIGU_demoArray_Par, allValues = TRUE)
# Subset matrix of all average values to just Total allelic, geographic, and ecological coverage
averageValueMat_TEG <- averageValueMat[,c(1,6,7)]

# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue', 'purple')
plotColors <- alpha(plotColors, 0.45)
plotColors_Sub <- plotColors[-(2:5)]

# ---- CORRELATION PLOTS
par(mfrow=c(2,1))
# ---- GEOGRAPHIC-GENETIC
plot(averageValueMat_TEG$Geo, averageValueMat_TEG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     main='M. guttatus: Geographic by genetic coverage',xlab='', ylab='')
mtext(text='255 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Geographic coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
mylabel = bquote(italic(R)^2 == .(format(MIGU_geoModel_rSquared, digits = 3)))
text(x = 2, y = 10, labels = mylabel, cex=0.8)
# ---- ECOLOGICAL-GENETIC
plot(averageValueMat_TEG$Eco, averageValueMat_TEG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     main='M. guttatus: Ecological by genetic coverage',xlab='', ylab='')
mtext(text='255 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Ecological coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
mylabel = bquote(italic(R)^2 == .(format(MIGU_ecoModel_rSquared, digits = 3)))
text(x = 2, y = 10, labels = mylabel, cex=0.8)

# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TEG, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='Mimulus guttatus: Gen-Geo-Eco Coverage', line=1.5)
mtext(text='255 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=105, y=35, inset = 0.05,
       legend = c('Genetic coverage (Total)', 'Geographic coverage (50 km buffer)', 'Ecological coverage (EPA Level III)'),
       col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=1.2, pt.cex = 2, bty='n',
       y.intersp = 0.5)
