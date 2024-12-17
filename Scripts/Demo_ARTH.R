# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: ARABIDOPSIS THALIANA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage
# for Arabidopsis thaliana. Uses data files from 1001 Genome Consortium et al. 2016--a VCF, for genetic data, 
# and a spreadsheet that was derived from supplemental data from the manuscript. 

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
geo_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500,1000,1500,2000))
# Specify ecological buffer size in meters 
eco_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500,1000,1500,2000))

# ---- READ IN DATA ----
# Specify filepath for ARTH geographic and genetic data
ARTH_filePath <- paste0(GeoGenCorr_wd, 'Datasets/ARTH/')

# ---- GEOGRAPHIC/ECOLOGICAL DATA FILES
# The original coordinates file for this Arabidopsis dataset needs to be processed such that 
# samples outside of the native range of Eurasia (in the U.S. and Japan) are removed
# Check if the processed file (called ARTH_coordinates.csv) already exists; if not, then 
# run the necessary processing steps.
if(file.exists(paste0(YUBR_filePath, 'Geographic/ARTH_coordinates.csv'))){
  # Read in the CSV of processed coordinates. The first column contains row numbers
  ARTH_coordinates <- read.csv(
    paste0(ARTH_filePath, 'Geographic/ARTH_coordinates.csv'), header=TRUE)
} else {
  # The metadata for the 1,135 accessions included in the analysis, including lat/long values,
  # was accessed from a CSV uploaded to the website here: https://1001genomes.org/accessions.html
  ARTH_coordinates <- 
    read.csv(file=paste0(ARTH_filePath, 'Geographic/ARTH_coordinates_Original.csv'), header = FALSE)
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
  # Write resulting coordinates data.frame as a CSV to disk, for future runs
  write.csv(ARTH_coordinates, file=paste0(ARTH_filePath,'Geographic/ARTH_coordinates.csv'), 
            row.names = FALSE)
}
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- grabWorldAdmin(GeoGenCorr_wd = GeoGenCorr_wd, fileExtentsion = ".gpkg", overwrite = TRUE)
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = ARTH_coordinates)
# Read in the TNC global ecoregion shapefile, which is used for calculating ecological coverage 
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_globalTNC/Terrestrial_Ecoregions.shp')))

# ---- GENETIC MATRIX
# Read in the VCF file provided via the 1001 Genomes Consortium 
# (website here: https://1001genomes.org/data/GMI-MPI/releases/v3.1/). Note that the VCF on the website
# contains the full genomic data for each individual; this was randomly subset to 10,000 loci using a simple
# BASH script
ARTH_vcf <- read.vcfR(file=paste0(ARTH_filePath, 'Genetic/ARTH_10k.vcf'))
# Convert the vcf to a genind; the return.alleles FALSE value allows for downstream genetic distance calculations
# This genind file is made up of 1,135 individuals and 10,000 loci
ARTH_genind_global <- vcfR2genind(ARTH_vcf, sep = "/", return.alleles = FALSE)
# REMOVE INTRODUCED POPULATIONS: Subset global genind object to only contain individuals from native range. 
# The 'drop' argument removes alleles no longer present in the dataset.
ARTH_genind <- ARTH_genind_global[ARTH_coordinates[,1], drop=TRUE]

# ---- PARALLELIZATION
# Flag for running resampling steps in parallel
parFlag <- TRUE

# If running in parallel, set up cores and export required libraries
if(parFlag==TRUE){
  # Set up relevant cores 
  num_cores <- detectCores() - 8 
  cl <- makeCluster(num_cores)
  # Make sure libraries (adegenet, terra, etc.) are on cluster (but avoid printing output)
  invisible(clusterEvalQ(cl, library('adegenet')))
  invisible(clusterEvalQ(cl, library('terra')))
  invisible(clusterEvalQ(cl, library('parallel')))
  invisible(clusterEvalQ(cl, library('usedist')))
  # Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
  # exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
  world_poly_clip_W <- wrap(world_poly_clip)
  ecoregion_poly_W <- wrap(ecoregion_poly)
}

# ---- RESAMPLING ----
if(parFlag==TRUE){
  # Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
  clusterExport(cl, varlist = c('ARTH_coordinates','ARTH_genind','num_reps','geo_buffSize', 'eco_buffSize',
                                'world_poly_clip_W', 'ecoregion_poly_W'))
  # Export necessary functions (for calculating geographic and ecological coverage) to the cluster
  clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres', 
                                'eco.intersectBuff','eco.compareBuff','gen.getAlleleCategories', 
                                'gen.buildDistMat', 'gen.calcGenDistCov', 'eco.totalEcoregionCount',
                                'calculateCoverage','exSituResample.Par', 'geo.gen.Resample.Par'))
  # Specify file path, for saving resampling array
  arrayDir <- paste0(ARTH_filePath, 'resamplingData/ARTH_SMBO3_G2GE_5r_resampArr.Rdata')
  
  # Run resampling (in parallel)
  ARTH_demoArray_Par <- 
    geo.gen.Resample.Par(genObj = ARTH_genind, genDistFlag = TRUE, geoFlag = TRUE, coordPts = ARTH_coordinates, 
                         geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag = TRUE, 
                         ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer = 'GL', 
                         reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
  # Close cores
  stopCluster(cl)
} else {
  # Specify file path, for saving resampling array
  arrayDir <- paste0(ARTH_filePath, 'resamplingData/ARTH_SMBO2_GE_5r_resampArr.Rdata')
  # Run resampling not in parallel (for function testing purposes)
  ARTH_demoArray <-
    geo.gen.Resample(genObj = ARTH_genind,  genDistFlag=TRUE, geoFlag = TRUE, coordPts = ARTH_coordinates,
                     geoBuff = geo_buffSize, boundary = world_poly_clip, ecoFlag = FALSE, 
                     ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly, ecoLayer = 'GL',
                     reps = 1)
  # Save the resampling array object to disk, for later usage
  saveRDS(ARTH_demoArray, file = arrayDir)
}

# %%% ANALYZE DATA %%% ----
# Specify filepath for ARTH geographic and genetic data, including resampling array
ARTH_filePath <- paste0(GeoGenCorr_wd, 'Datasets/ARTH/')
arrayDir <- paste0(ARTH_filePath, 'resamplingData/ARTH_50km_GE_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
ARTH_demoArray_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values
ARTH_DF <- resample.array2dataframe(ARTH_demoArray_Par)
# Calculate normalized root mean square value
ARTH_nrmse_geo <- nrmse_func(obs=ARTH_DF$Geo, pred=ARTH_DF$Total) ; ARTH_nrmse_geo
ARTH_nrmse_eco <- nrmse_func(obs=ARTH_DF$Eco, pred=ARTH_DF$Total) ; ARTH_nrmse_eco

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
# Add NRMSE values for each comparison
text(x = 76, y = 35, labels = paste0('NRMSE: ', ARTH_nrmse_geo), col='darkblue', cex=0.9)
text(x = 76, y = 20, labels = paste0('NRMSE: ', ARTH_nrmse_eco), col='purple', cex=0.9)
# Add legend
legend(x=55, y=55, inset = 0.05, xpd=TRUE,
       legend = c('Geographic', 'Ecological'),
       col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.8)
# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TEG[1:3], ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='A. thaliana: Geo-Eco-Gen Coverage', line=1.5)
mtext(text='1,010 Individuals; 50 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=600, y=65, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (50 km buffer)',
                  'Ecological coverage (TNC global ecoregions)'),
       col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 0.8)

# %%%% SMBO: MULTIPLE BUFFER SIZES ----
# Specify filepath for ARTH geographic and genetic data, including resampling array
ARTH_filePath <- paste0(GeoGenCorr_wd, 'Datasets/ARTH/')
arrayDir <- paste0(ARTH_filePath, 'resamplingData/ARTH_SMBO3_G2GE_5r_resampArr.Rdata')
# Read in array and build a data.frame of values
ARTH_SMBO3_array <- readRDS(arrayDir)

# ---- CALCULATIONS ----
# Build a data.frame from array values
ARTH_SMBO3_DF <- resample.array2dataframe(ARTH_SMBO3_array)
# Build tables of NRSMSE values, calculated based on data.frame
ARTH_NRMSE_Mat_CV <- buildNRMSEmatrix(resampDF=ARTH_SMBO3_DF, genCovType='CV', sdmFlag=FALSE)
ARTH_NRMSE_Mat_GD <- buildNRMSEmatrix(resampDF=ARTH_SMBO3_DF, genCovType='GD', sdmFlag=FALSE)
# Combine the results of the NRMSE values calculated using allelic coverages and using
# genetic distances, and then rename the columns accordingly
ARTH_NRMSE_Mat <- cbind(ARTH_NRMSE_Mat_CV, ARTH_NRMSE_Mat_GD)
# Store the matrix as a CSV to disk
write.table(ARTH_NRMSE_Mat,
            file=paste0(ARTH_filePath, 'resamplingData/ARTH_SMBO3_NRMSE.csv'), sep=',')

# SMBO2: OPTIMAL BUFFER SIZES ----
# Read in ARTH SMBO2 resampling array amd convert to data.frame
ARTH_filePath <- paste0(GeoGenCorr_wd, 'Datasets/ARTH/')
ARTH_arrayDir <- paste0(ARTH_filePath, 'resamplingData/ARTH_SMBO2_GE_5r_resampArr.Rdata')
# From ARTH resampling array, return a matrix of average coverage values for optimal buffer sizes
ARTH_optCovMat <- extractOptCovs(ARTH_arrayDir)
# Calculate MSSEs: minimum number of samples for 95% of each coverage type
ARTH_Gen_MSSE <- min(which(ARTH_optCovMat[,1] > 95)) ; ARTH_Gen_MSSE
ARTH_GeoBuff_MSSE <- min(which(ARTH_optCovMat[,2] > 95)) ; ARTH_GeoBuff_MSSE
ARTH_Eco_MSSE <- min(which(ARTH_optCovMat[,3] > 95)) ; ARTH_Eco_MSSE

# PLOTTING
# Specify plot colors
plotColors <- c('red', 'darkblue','darkorange3', 'purple')
plotColors_Fade <- alpha(plotColors, c(0.45, rep(0.85, length(plotColors)-1)))
# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1), oma=rep(0.2,4), mar=c(2,4,3,1)+0.1)
# Geo Buff
matplot(ARTH_optCovMat[,c(1,2)], ylim=c(0,110), col=plotColors_Fade[c(1, 2)], pch=16, ylab='')
abline(h=95, col="black", lty=3) 
abline(v=ARTH_Gen_MSSE, col="red") ; abline(v=ARTH_GeoBuff_MSSE, col="darkblue")
mtext(text=paste0('MSSE: ', ARTH_Gen_MSSE), side=1, line=-1, at=ARTH_Gen_MSSE+45, cex=0.7, col='red')
mtext(text=paste0(' MSSE: ', ARTH_GeoBuff_MSSE), line=-1.5, side=1, at=ARTH_GeoBuff_MSSE-45, cex=0.7, col='darkblue')
title('Arabidopsis thaliana: Coverages at Optimal Buffer Sizes', cex.sub=1.2, line = 2)
mtext(text='Geographic (Total Buffer): 500 km', side=3, at=80, cex=0.8)
# Eco Buff
par(mar=c(3,4,2,1)+0.1)
matplot(ARTH_optCovMat[,c(1,3)], ylim=c(0,110), col=plotColors_Fade[c(1, 4)], pch=16, ylab='')
abline(h=95, col="black", lty=3)
abline(v=ARTH_Gen_MSSE, col="red") ; abline(v=ARTH_Eco_MSSE, col="purple")
mtext(text=paste0('MSSE: ', ARTH_Gen_MSSE), side=1, line=-1, at=ARTH_Gen_MSSE+45, cex=0.7, col='red')
mtext(text=paste0(' MSSE: ', ARTH_Eco_MSSE), line=-1.5, side=1, at=ARTH_Eco_MSSE-45, cex=0.7, col='purple')
mtext(text='Ecological: 500 km', side=3, at=c(80), cex=0.8)
