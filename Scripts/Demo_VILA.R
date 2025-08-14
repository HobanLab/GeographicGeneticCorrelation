# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: VITIS LABRUSCA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage
# for Vitis labrusca. Uses data files provided by Dr. Zoe Migicovsky

pacman::p_load(adegenet, terra, parallel, RColorBrewer, scales, vcfR, usedist, dartR)

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
# Specify filepath for VILA geographic and genetic data
VILA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/VILA/')

# ---- GEOGRAPHIC/ECOLOGICAL DATA FILES
# The original coordinates file for this Vitis dataset needs to be processed such that 
# ex situ individuals are removed.
# Check if the processed file (called VILA_coordinates.csv) already exists; if not, then 
# run the necessary processing steps.
if(file.exists(paste0(VILA_filePath, 'Geographic/VILA_coordinates.csv'))){
  # Read in the CSV of processed coordinates. The first column contains row numbers
  VILA_coordinates <- read.csv(
    paste0(VILA_filePath, 'Geographic/VILA_coordinates.csv'), header=TRUE)
} else {
  # A text file was provided which specifies the sample name, latitude, longitude, 
  # source, and species for each individual. Read this in (using read.delim for TSV)
  VILA_coordinates <- 
    read.table(paste0(VILA_filePath, 'Geographic/samples_keyfile_austin.txt'), header = TRUE)
  # Subset coordinate dataframe to strictly wild samples by removing germplasm individuals
  VILA_coordinates <- VILA_coordinates[which(VILA_coordinates$POPTYPE=='Wild'),]
  # Remove unnecessary columns (species, state, and poptype), reorder and rename columns
  VILA_coordinates <- VILA_coordinates[,c(1,4,3)]
  colnames(VILA_coordinates) <- c('sampleID','decimalLatitude','decimalLongitude')
  # Write resulting coordinates data.frame as a CSV to disk, for future runs
  write.csv(VILA_coordinates, file=paste0(VILA_filePath,'Geographic/VILA_coordinates.csv'), 
            row.names = FALSE)
}
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.shp')))
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = VILA_coordinates)
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage 
# (solely in the U.S.)
ecoregion_poly <-
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level4/us_eco_l4.shp')))
# Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
# exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
world_poly_clip_W <- wrap(world_poly_clip)
ecoregion_poly_W <- wrap(ecoregion_poly)

# ---- GENETIC MATRIX
# Read in the VCF file for Vitis labrusca using vcfR::read.vcfR
VILA_vcf <- read.vcfR(paste0(VILA_filePath,'Genetic/VILA.vcf'))
# Convert the vcf to a genlight object, then the genlight to a genind object (via dartR function)
# This extra conversion step is necessary because of the formatting of loci names in the VCF
VILA_genlight <- vcfR2genlight(VILA_vcf)
VILA_genind <- gl2gi(VILA_genlight)
# Individual names are duplicated in genind object (such that 104_5 is called 
# 104_5_104_5). To address this, first, use the sub command to replace names
indNames(VILA_genind) <- sub("(.*?_.*?)_.*", "\\1", indNames(VILA_genind))
# Match the individuals in the genind file to those in the coordinate dataframe.
# This essentially drops 25 germplasm individuals
VILA_genind <- VILA_genind[VILA_coordinates[,1], drop=TRUE]

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('VILA_coordinates','VILA_genind','num_reps','geo_buffSize', 'eco_buffSize',
                              'world_poly_clip_W', 'ecoregion_poly_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres', 
                              'eco.intersectBuff','eco.compareBuff','gen.getAlleleCategories', 
                              'gen.buildDistMat', 'gen.calcGenDistCov', 'eco.totalEcoregionCount',
                              'calculateCoverage','exSituResample.Par', 'geo.gen.Resample.Par'))
# Specify file path, for saving resampling array
arrayDir <- paste0(VILA_filePath, 'resamplingData/VILA_SMBO2_5r_resampArr.Rdata')

# Run resampling (in parallel)
VILA_demoArray_Par <- 
  geo.gen.Resample.Par(genObj=VILA_genind, geoFlag=TRUE, genDistFlag=FALSE, coordPts=VILA_coordinates, 
                       geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag=TRUE, ecoBuff=eco_buffSize, 
                       ecoRegions=ecoregion_poly_W, ecoLayer='US', reps=num_reps, arrayFilepath=arrayDir, cluster=cl)
# Close cores
stopCluster(cl)

# # Run resampling not in parallel (for function testing purposes)
# VILA_coordinates_test <- VILA_coordinates[1:5,]
# VILA_genind_test <- VILA_genind[VILA_coordinates_test[,1], drop=TRUE]
# VILA_demoArray_IND <-
#   geo.gen.Resample(genObj = VILA_genind_test, geoFlag = TRUE, coordPts = VILA_coordinates_test,
#                    geoBuff = geo_buffSize, boundary = world_poly_clip, ecoFlag = TRUE,
#                    ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly, ecoLayer = "US", reps = 1)

# %%% ANALYZE DATA %%% ----
# Specify filepath for VILA geographic and genetic data, including resampling array
VILA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/VILA/')
arrayDir <- paste0(VILA_filePath, 'resamplingData/VILA_1km_GE_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
VILA_demoArray_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values
VILA_DF <- resample.array2dataframe(VILA_demoArray_Par)
# Calculate normalized root mean square value
VILA_nrmse_geo <- nrmse_func(obs=VILA_DF$Geo, pred=VILA_DF$Total) ; VILA_nrmse_geo
VILA_nrmse_eco <- nrmse_func(obs=VILA_DF$Eco, pred=VILA_DF$Total) ; VILA_nrmse_eco

# ---- PLOTTING ----
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(VILA_demoArray_Par, allValues = TRUE)
# Subset matrix of all average values to just Total allelic, geographic, and ecological coverage
averageValueMat_TEG <- averageValueMat[,c(1,6,7)]
# Calculate the absolute difference between genetic and geographic/ecological, and add to data.frame
averageValueMat_TEG <- cbind(averageValueMat_TEG, abs(averageValueMat_TEG$Total-averageValueMat_TEG$Geo))
averageValueMat_TEG <- cbind(averageValueMat_TEG, abs(averageValueMat_TEG$Total-averageValueMat_TEG$Eco))
names(averageValueMat_TEG) <- c(names(averageValueMat_TEG)[1:3], 'Geo_Difference', 'Eco_Difference')

# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','darkblue', 'purple')
plotColors_Fade <- alpha(plotColors, 0.65)
plotColors_Sub <- plotColors_Fade[-(2:4)]
# Two plots in a single window
par(mfrow=c(2,1))
# ---- CORRELATION PLOTS
plot(averageValueMat_TEG$Geo, averageValueMat_TEG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     main='V. labrusca: Geographic by genetic coverage',xlab='', ylab='', col=plotColors_Fade[[5]])
mtext(text='220 Individuals; 1 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Geographic/Ecological coverage (%)', side=1, line=3, cex=1.6)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add points for ecological coverage
points(x=averageValueMat$Eco, y=averageValueMat$Total, pch=20, col=plotColors_Fade[[6]])
# Add NRMSE values for each comparison
text(x = 76, y = 35, labels = paste0('NRMSE: ', VILA_nrmse_geo), col='darkblue', cex=0.9)
text(x = 76, y = 20, labels = paste0('NRMSE: ', VILA_nrmse_eco), col='purple', cex=0.9)
# Add legend
legend(x=58, y=162, inset = 0.05, xpd=TRUE,
       legend = c('Geographic', 'Ecological'),
       col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 0.08)
# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TEG[,1:3], ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='V. labrusca: Geo-Eco-Gen Coverage', line=1.5)
mtext(text='220 Individuals; 1 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Add legend
legend(x=85, y=180, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (1 km buffer)',
                  'Ecological coverage (1 km buffer, EPA Level IV)'),
       col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n',
       y.intersp = 0.08)

# %%%% SMBO: MULTIPLE BUFFER SIZES ----
# Specify filepath for VILA geographic and genetic data, including resampling array
VILA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/VILA/')
arrayDir <- paste0(VILA_filePath, 'resamplingData/VILA_SMBO2_5r_resampArr.Rdata')
# Read in array and build a data.frame of values
VILA_SMBO2_array <- readRDS(arrayDir)
# Specify geographic buffer size in meters (used above)
geo_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500))

# ---- CALCULATIONS ----
# Build a data.frame from array values
VILA_SMBO2_DF <- resample.array2dataframe(VILA_SMBO2_array)
# Build a matrix to capture NRMSE values
VILA_NRMSE_Mat <- matrix(NA, nrow=length(geo_buffSize), ncol=2)
# The names of this matrix match the different parts of the dataframe names
colnames(VILA_NRMSE_Mat) <- c('Geo_Buff','Eco_Buff')
rownames(VILA_NRMSE_Mat) <- paste0(geo_buffSize/1000, 'km')
# Loop through the dataframe columns. The first two columns are skipped, as they're sampleNumber and the
# predictve variable (genetic coverages)
for(i in 3:ncol(VILA_SMBO2_DF)){
  # Calculate NRMSE for the current column in the dataframe
  VILA_NRMSEvalue <- nrmse.func(VILA_SMBO2_DF[,i], pred = VILA_SMBO2_DF$Total)
  # Get the name of the current dataframe column
  dataName <- unlist(strsplit(names(VILA_SMBO2_DF)[[i]],'_'))
  # Match the data name to the relevant rows/columns of the receiving matrix
  matRow <- which(rownames(VILA_NRMSE_Mat) == dataName[[3]])
  matCol <- which(colnames(VILA_NRMSE_Mat) == paste0(dataName[[1]],'_',dataName[[2]]))
  # Locate the NRMSE value accordingly
  VILA_NRMSE_Mat[matRow,matCol] <- VILA_NRMSEvalue
}
print(VILA_NRMSE_Mat)
# Store the matrix as a CSV to disk
write.table(VILA_NRMSE_Mat,
            file=paste0(VILA_filePath, 'resamplingData/VILA_SMBO2_NRMSE.csv'), sep=',')

# %%%% SMBO3 ----
# Specify filepath for VILA geographic and genetic data, including resampling array
VILA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/VILA/')
arrayDir <- paste0(VILA_filePath, 'resamplingData/VILA_SMBO3_5r_resampArr.Rdata')
# Read in array
VILA_SMBO3_array <- readRDS(arrayDir)

# ---- CALCULATIONS ----
# Build a data.frame from array values
VILA_SMBO3_DF <- resample.array2dataframe(VILA_SMBO3_array)
# Build tables of NRSMSE values, calculated based on data.frame
VILA_NRMSE_Mat_CV <- buildNRMSEmatrix(resampDF=VILA_SMBO3_DF, genCovType='CV', sdmFlag=FALSE)
VILA_NRMSE_Mat_GD <- buildNRMSEmatrix(resampDF=VILA_SMBO3_DF, genCovType='GD', sdmFlag=FALSE)
# Combine the results of the NRMSE values calculated using allelic coverages and using
# genetic distances, and then rename the columns accordingly
VILA_NRMSE_Mat <- cbind(VILA_NRMSE_Mat_CV, VILA_NRMSE_Mat_GD)
# Store the matrix as a CSV to disk
write.table(VILA_NRMSE_Mat,
            file=paste0(VILA_filePath, 'resamplingData/VILA_SMBO3_NRMSE.csv'), sep=',')

# SMBO2: OPTIMAL BUFFER SIZES ----
# Read in VILA SMBO2 resampling array amd convert to data.frame
VILA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/VILA/')
VILA_arrayDir <- paste0(VILA_filePath, 'resamplingData/VILA_SMBO2_5r_resampArr.Rdata')
# From VILA resampling array, return a matrix of average coverage values for optimal buffer sizes
VILA_optCovMat <- extractOptCovs(VILA_arrayDir)
# Calculate MSSEs: minimum number of samples for 95% of each coverage type
VILA_Gen_MSSE <- min(which(VILA_optCovMat[,1] > 95)) ; VILA_Gen_MSSE
VILA_GeoBuff_MSSE <- min(which(VILA_optCovMat[,2] > 95)) ; VILA_GeoBuff_MSSE
VILA_Eco_MSSE <- min(which(VILA_optCovMat[,3] > 95)) ; VILA_Eco_MSSE

# PLOTTING
# Specify plot colors
plotColors <- c('red', 'darkblue','purple')
plotColors_Fade <- alpha(plotColors, c(0.45, rep(0.85, length(plotColors)-1)))
# Set plotting window to stack 3 graphs vertically
par(mfcol=c(2,1), oma=rep(0.2,4), mar=c(2,4,3,1)+0.1)
# Geo Buff
matplot(VILA_optCovMat[,c(1,2)], ylim=c(0,110), col=plotColors_Fade[c(1, 2)], pch=16, ylab='')
abline(h=95, col="black", lty=3)
abline(v=VILA_Gen_MSSE, col="red") ; abline(v=VILA_GeoBuff_MSSE, col="darkblue")
mtext(text=paste0('MSSE: ', VILA_Gen_MSSE), side=1, line=-1, at=VILA_Gen_MSSE+8, cex=0.8, col='red')
mtext(text=paste0(' MSSE: ', VILA_GeoBuff_MSSE), line=-1.5, side=1, at=VILA_GeoBuff_MSSE-8, cex=0.8, col='darkblue')
title('Vitis labrusca: Coverages at Optimal Buffer Sizes', cex.sub=1.2, line = 2)
mtext(text='Geographic (Total Buffer): 250 km', side=3, at=15, cex=0.9)
# Legend
legend(x=160, y=95, xpd=TRUE, cex=1, pch=rep(19,3),
       col=c('red','darkblue','purple'),
       legend=c('Genetic', 'Geographic (Total Buffer)', 'Ecological'),
       y.intersp = 0.3, bty='n')
# Eco Buff
par(mar=c(3,4,2,1)+0.1)
matplot(VILA_optCovMat[,c(1,3)], ylim=c(0,110), col=plotColors_Fade[c(1, 3)], pch=16, ylab='')
abline(h=95, col="black", lty=3)
abline(v=VILA_Gen_MSSE, col="red") ; abline(v=VILA_Eco_MSSE, col="purple")
mtext(text=paste0('MSSE: ', VILA_Gen_MSSE), side=1, line=-1, at=VILA_Gen_MSSE+8, cex=0.8, col='red')
mtext(text=paste0(' MSSE: ', VILA_Eco_MSSE), line=-1.5, side=1, at=VILA_Eco_MSSE-8, cex=0.8, col='purple')
mtext(text='Ecological: 190 km', side=3, at=15, cex=0.9)
