# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: ARABIDOPSIS THALIANA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage
# for Arabidopsis thaliana. Uses data files from 1001 Genome Consortium et al. 2016--a VCF, for genetic data, 
# and a spreadsheet that was derived from supplemental data from the manuscript. 

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
# Specify filepath for ARTH geographic and genetic data
ARTH_filePath <- paste0(GeoGenCorr_wd, 'Datasets/ARTH/')

# ---- GEOGRAPHIC/ECOLOGICAL DATA FILES
# The metadata for the 1,135 accessions included in the analysis, including latitude and longitude values,
# was accessed from a CSV uploaded to the website here: https://1001genomes.org/accessions.html
ARTH_coordinates <- 
  read.csv(file=paste0(ARTH_filePath, 'Geographic/ARTH_accessionMetadata.csv'), header = FALSE)
# Remove unnecessary columns (CS number, collector, sequencer, etc.), and rename columns. Retain 
# country values, in order to filter out samples outside of the native range of Eurasia
ARTH_coordinates <- ARTH_coordinates[,-c(2:3,5,8:13)]
# Rename CSV columns, and drop unnecessary columns
colnames(ARTH_coordinates) <- c('Acc_ID','Country','decimalLatitude','decimalLongitude')
# Remove samples that come from the U.S. (USA) or Japan (JPN) (leaving 1,010 samples)
ARTH_coordinates <- ARTH_coordinates[-which(ARTH_coordinates$Country == 'USA'),]
ARTH_coordinates <- ARTH_coordinates[-which(ARTH_coordinates$Country == 'JPN'),]
# Remove the country column, and reformat sample names as characters (rather than numeric)
ARTH_coordinates <- ARTH_coordinates[,-2]
ARTH_coordinates[,1] <- as.character(ARTH_coordinates[,1])

# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.shp')))
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = ARTH_coordinates)
# Read in the TNC global ecoregion shapefile, which is used for calculating ecological coverage 
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_globalTNC/Terrestrial_Ecoregions.shp')))
# Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
# exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
world_poly_clip_W <- wrap(world_poly_clip)
ecoregion_poly_W <- wrap(ecoregion_poly)

# ---- GENETIC MATRIX
# Read in the VCF file provided via the 1001 Genomes Consortium 
# (website here: https://1001genomes.org/data/GMI-MPI/releases/v3.1/). Note that the VCF on the website
# contains the full genomic data for each individual; this was randomly subset to 10,000 loci using a simple
# BASH script
ARTH_vcf <- read.vcfR(file=paste0(ARTH_filePath, 'Genetic/ARTH_10k.vcf'))
# Convert the vcf to a genind; the return.alleles TRUE value is suggested in the function's help file
# This genind file is made up of 1,135 individuals and 10,000 loci
ARTH_genind_global <- vcfR2genind(ARTH_vcf, sep = "/", return.alleles = TRUE)
# REMOVE INTRODUCED POPULATIONS: Subset global genind object to only contain individuals from native range. 
# The 'drop' argument removes alleles no longer present in the dataset.
ARTH_genind <- ARTH_genind_global[ARTH_coordinates[,1], drop=TRUE]

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('ARTH_coordinates','ARTH_genind','num_reps','geo_buffSize', 'eco_buffSize',
                              'world_poly_clip_W', 'ecoregion_poly_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers', 'geo.compareBuff', 'geo.compareBuffSDM', 
                              'eco.intersectBuff', 'eco.compareBuff', 'gen.getAlleleCategories',
                              'calculateCoverage', 'exSituResample.Par', 'geo.gen.Resample.Par'))
# Specify file path, for saving resampling array
arrayDir <- paste0(ARTH_filePath, 'resamplingData/ARTH_50km_GE_5r_resampArr.Rdata')

# Run resampling (in parallel)
print("%%% BEGINNING RESAMPLING %%%")
ARTH_demoArray_Par <- 
  geo.gen.Resample.Par(gen_obj = ARTH_genind, geoFlag = TRUE, coordPts = ARTH_coordinates, 
                       geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag = TRUE, 
                       ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer = 'GL', 
                       reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
# Close cores
stopCluster(cl)

# %%% ANALYZE DATA %%% ----
# Specify filepath for ARTH geographic and genetic data, including resampling array
ARTH_filePath <- paste0(GeoGenCorr_wd, 'Datasets/ARTH/')
arrayDir <- paste0(ARTH_filePath, 'resamplingData/ARTH_50km_GE_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
ARTH_demoArray_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values, to pass to linear models
ARTH_DF <- resample.array2dataframe(ARTH_demoArray_Par)
# Calculate Spearman's r for geographic/ecological coverage
ARTH_spearR_geo <- round(cor(ARTH_DF$Geo, ARTH_DF$Total, method = 'spearman'),3) ; ARTH_spearR_geo
ARTH_spearR_eco <- round(cor(ARTH_DF$Eco, ARTH_DF$Total, method = 'spearman'),3) ; ARTH_spearR_eco

# ---- PLOTTING ----
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(ARTH_demoArray_Par, allValues = TRUE)
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
     main='A. thaliana: Geographic by genetic coverage',xlab='', ylab='', col=plotColors[[5]])
mtext(text='1,010 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Geographic/Ecological coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add points for ecological coverage
points(x=averageValueMat$Eco, y=averageValueMat$Total, pch=20, col=plotColors[[6]])
# Add Spearman's r values for each comparison
text(x = 76, y = 35, labels = paste0('Spearman r: ', ARTH_spearR_geo), col='darkblue', cex=0.9)
text(x = 76, y = 20, labels = paste0('Spearman r: ', ARTH_spearR_eco), col='purple', cex=0.9)
# Add legend
legend(x=55, y=70, inset = 0.05, xpd=TRUE,
       legend = c('Geographic', 'Ecological'),
       col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.35)
# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TEG[1:3], ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='A. thaliana: Geo-Eco-Gen Coverage', line=1.5)
mtext(text='1,010 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=600, y=70, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (50 km buffer)',
                  'Ecological coverage (TNC global ecoregions)'),
       col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 0.35)
# ---- DIFFERENCE PLOTS
# Plot difference between geographic and genetic coverage
matplot(averageValueMat_TEG[4:5], col=plotColors[5:6], pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='A. thaliana: Genetic-Geographic-Ecological Coverage Difference', line=1.5)
mtext(text='1,010 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.2)
mtext(text='Difference in Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=650, y=45, inset = 0.05,
       legend = c('Genographic coverage difference', 'Ecological coverage difference'),
       col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 1)
