# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: CONRADINA GLABRA USING HAPLOTYPES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage
# for Conradina glabra, a dataset originally provided by Lauren Eserman. This script explores
# the impact of using haplotypes of different lengths on genetic resampling curves. To do this,
# resampling is performed multiple times using different input files, with the results of each
# haplotype length stored to a different array object.

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
geo_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500))
# Specify ecological buffer size in meters 
eco_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500))
# ---- SHAPEFILES
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.shp')))
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage 
# (solely in the U.S.)
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
# Specify filepath for COGL geographic and genetic data
COGL_filePath <- paste0(GeoGenCorr_wd, 'Datasets/COGL/')

# ---- GENETIC MATRICES
# For each haplotype length (ranging from a single SNP to 10 concatenated SNPs), the genetic
# data is stored as a STRUCTURE input file. Each STRUCTURE file is read in and converted into a
# genind object. 
COGL_haplos_filePath <- paste0(COGL_filePath, 'Genetic/haplotypicData')
COGL_haplos_input <- list.files(COGL_haplos_filePath, full.names = TRUE)
# Declare a list of integers, storing the number of loci for each STRUCTURE input file
COGL_lociLevels <- c(16291,7045,4001,2509,1627,1083,727,473,330,212)
# Predeclare an empty list to store genind objects
COGL_genList <- list()
# For each STRUCTURE input file, 
for(i in 1:length(COGL_haplos_input)){
  # Read the file in and convert it into a genind object
  COGL_genList[[i]] <- read.structure(file=COGL_haplos_input[[i]], n.ind=564, 
                                    n.loc=COGL_lociLevels[[i]],
                                    onerowperind=TRUE, col.lab=3, col.pop=1,
                                    col.others=2, row.marknames=0, NA.char=-9,
                                    sep=" ", ask=TRUE, quiet=TRUE)
  # Rename the individuals in the file. This is for sample names to match the coordinate file
  indNames(COGL_genList[[i]]) <- gsub("_sorted", "", indNames(COGL_genList[[i]]))
  indNames(COGL_genList[[i]]) <- paste0("Congla", indNames(COGL_genList[[i]]))
  # Remove two duplicate samples (with "_rep" in sample name), leaving 562 samples
  COGL_genList[[i]] <- COGL_genList[[i]][-grep('_rep', indNames(COGL_genList[[i]])), drop=TRUE]
  # Specify the population of every sample to be 'wild'
  pop(COGL_genList[[i]]) <- rep('wild', nInd(COGL_genList[[i]]))
}

# ---- COORDINATE POINTS
# CSV was provided by Lauren Eserman. Because individuals were sampled very close 
# together, coordinates were determined by obtaining coordinates for the center of a
# "clump" of individuals, and then individual coordinates were calculated based on the
# separation and direction of each individual from the center. Some samples have identical
# coordinate values.
# NOTE: these coordinates cannot be shared externally, due to the rare status of this species!
COGL_coordinates <- 
  read.csv(file=paste0(COGL_filePath, 'Geographic/COGL_coordinates.csv'), header = TRUE)
# Rename the columns of the geographic coordinates dataframe (because geo.compareBuff function expects certain strings)
colnames(COGL_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')
# Remove two duplicate samples (with "_rep" in sample name), leaving 562 samples
COGL_coordinates<- COGL_coordinates[-grep('_rep', COGL_coordinates$Sample.Name),]

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('COGL_coordinates','COGL_genList','num_reps','geo_buffSize', 'eco_buffSize',
                              'world_poly_clip_W', 'ecoregion_poly_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres', 
                              'eco.intersectBuff','eco.compareBuff','gen.getAlleleCategories',
                              'eco.totalEcoregionCount','calculateCoverage','exSituResample.Par',
                              'geo.gen.Resample.Par'))
# Specify file paths, for saving resampling arrays (if/else statement is for leading zeros)
COGL_haplos_output <- dirname(gsub("Genetic", "resamplingData", COGL_haplos_input))
for(i in 1:length(COGL_haplos_output)){
  if(i>10){
    COGL_haplos_output[[i]] <- paste0(COGL_haplos_output[[i]],'/SMBO2/COGL_SMBO2_GE_5r_Hap-0', i, '_resampArr.Rdata')
  } else{
    COGL_haplos_output[[i]] <- paste0(COGL_haplos_output[[i]],'/SMBO2/COGL_SMBO2_GE_5r_Hap-', i, '_resampArr.Rdata')
  }
}

# Run resampling (in parallel). Use a loop to iterate through the different haplotype lengths
for(i in 1:length(COGL_haplos_output)){
  print(paste0('%%% HAPLOTYPE LENGTH: ', i))
  # Run resampling
  COGL_demoArray_Par <- 
    geo.gen.Resample.Par(genObj = COGL_genList[[i]], geoFlag = TRUE, coordPts = COGL_coordinates, 
                         geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag = TRUE, 
                         ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer = 'US', 
                         reps = num_reps, arrayFilepath = COGL_haplos_output[[i]], cluster = cl)
}
# Close cores
stopCluster(cl)

# Run resampling not in parallel (for function testing purposes)
# COGL_demoArray_IND <-
#   geo.gen.Resample(gen_obj = unlist(COGL_genList[[9]]), geoFlag = TRUE, coordPts = COGL_coordinates,
#                    geoBuff = geo_buffSize, boundary = world_poly_clip, ecoFlag = FALSE, reps = 1)

# %%% ANALYZE DATA %%% ----
# Specify filepath for COGL geographic and genetic data, including resampling data
COGL_filePath <- paste0(GeoGenCorr_wd, 'Datasets/COGL/')
arrayDir <- paste0(COGL_filePath, 'resamplingData/haplotypicData/')
# Read in the resampling array .Rdata objects, saved to disk. Then, calculate the average
# value matrix (across resampling replicates) for each resampling array, and add that
# matrix as an item to a list.
COGL_haplos_output <- list.files(arrayDir, full.names = TRUE); COGL_haplos_averageValueMats <- list()
for(i in 1:length(COGL_haplos_output)){
  COGL_haplos_averageValueMats[[i]] <- meanArrayValues(readRDS(COGL_haplos_output[[i]]))
}
# Build a summary matrix that consists of the Total allelic representation values for each
# haplotype dataset (10 columns) and an average of all Geographic resampling values (1 column).
COGL_haplos_summaryMat <- array(data=NA, dim=c(nrow(COGL_haplos_averageValueMats[[1]]),11))
colnames(COGL_haplos_summaryMat) <- c(paste0(rep('Total_Hap', 10),seq(1:10)),'Geo')
# Total allelic representation values: pull the "Total" columns from each average value matrix,
# and cbind these together
COGL_haplos_summaryMat[,1:10] <- do.call(cbind,lapply(COGL_haplos_averageValueMats, function(x) x$Total))
# For the geographic coverages: build a matrix of the geographic coverage values
# from each average value matrix (using sapply). Then, calculate the means across the
# rows of this matrix (using rowMeans), and pass this to the last column of the summary matrix
COGL_haplos_summaryMat[,11] <- rowMeans(sapply(COGL_haplos_averageValueMats, function(df) df[, 2]))
# Convert the summary data into a data.frame
COGL_haplos_summaryMat <- as.data.frame(COGL_haplos_summaryMat)
# %%% NORMALIZED ROOT MEAN SQUARE ERROR: calculate the NRMSE for each haplotype length
COGL_NRMSE_Values <- vector(length=length(COGL_haplos_output))
for(i in 1:length(COGL_NRMSE_Values)){
  COGL_NRMSE_Values[i] <- nrmse_func(COGL_haplos_summaryMat$Geo, pred=COGL_haplos_summaryMat[,i])
}

# ---- PLOTTING ----
hapColors <- c('wheat2','tan','gold2','orange','salmon',
               'darkorange2','tomato3','red','red4','salmon4','darkblue')
hapColors_Fade <- alpha(hapColors, 0.5)
legText <- c(paste0(rep('Hap. length: ', 10), seq(1:10)), 'Geographic coverage (1 km buffer)')
# ---- COVERAGE PLOTS
matplot(COGL_haplos_summaryMat, ylim=c(0,110), col=hapColors_Fade, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='C. glabra: Haplotypic Coverages', line=1.5)
mtext(text='562 Individuals; 1 km Geographic buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.2)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.2, srt=90)
# Add legend
legend(x=400, y=87, inset = 0.05, legend = legText, col=hapColors, pch = c(19,19),
       cex=1.2, pt.cex = 2, bty='n', y.intersp = 0.5)

# %%%% SMBO2: MULTIPLE BUFFER SIZES ----
# Specify filepath for SMBO2 resampling arrays
COGL_filePath <- paste0(GeoGenCorr_wd, 'Datasets/COGL/resamplingData/haplotypicData/SMBO2/')
# Specify geographic buffer size in meters (used above)
geo_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500))
# Declare a list of the resampling array objects stored in the data directory
COGL_SMBO2_arrays <- grep('.Rdata', list.files(COGL_filePath, full.names = TRUE), value = TRUE, )
# ---- CALCULATIONS ----
# Loop through the list of arrays
for(i in 1:length(COGL_SMBO2_arrays)){
  # Read in the array
  COGL_SMBO2_array <- readRDS(COGL_SMBO2_arrays[[i]])
  # Build a dataframe from array values
  COGL_SMBO2_DF <- resample.array2dataframe(COGL_SMBO2_array)
  # Build a matrix to capture NRMSE values
  COGL_NRMSE_Mat <- matrix(NA, nrow=length(geo_buffSize), ncol=2)
  # The names of this matrix match the different parts of the dataframe names
  colnames(COGL_NRMSE_Mat) <- c('Geo_Buff','Eco_Buff')
  rownames(COGL_NRMSE_Mat) <- paste0(geo_buffSize/1000, 'km')
  # Loop through the dataframe columns. The first two columns are skipped, as they're sampleNumber and the
  # predictve variable (genetic coverages)
  for(j in 3:ncol(COGL_SMBO2_DF)){
    # Calculate NRMSE for the current column in the dataframe
    COGL_NRMSEvalue <- nrmse.func(COGL_SMBO2_DF[,j], pred = COGL_SMBO2_DF$Total)
    # Get the name of the current dataframe column
    dataName <- unlist(strsplit(names(COGL_SMBO2_DF)[[j]],'_'))
    # Match the data name to the relevant rows/columns of the receiving matrix
    matRow <- which(rownames(COGL_NRMSE_Mat) == dataName[[3]])
    matCol <- which(colnames(COGL_NRMSE_Mat) == paste0(dataName[[1]],'_',dataName[[2]]))
    # Locate the NRMSE value accordingly
    COGL_NRMSE_Mat[matRow,matCol] <- COGL_NRMSEvalue
  }
  cat(paste0('\n','%%% HAPLOTYPE LENGTH: ',i,' %%%','\n'))
  print(COGL_NRMSE_Mat)
  # Store the matrix as a CSV to disk
  write.table(COGL_NRMSE_Mat,
              file=paste0(COGL_filePath, 'COGL_SMBO2_NRMSE_Hap-', i, '.csv'), sep=',')
}

# %%%% SMBO2: HAPLOTYPES: MULTIPLE BUFFER SIZES ----
# Specify filepath for SMBO2 resampling arrays
COGL_filePath <- paste0(GeoGenCorr_wd, 'Datasets/COGL/resamplingData/haplotypicData/SMBO2/')
# Specify geographic buffer size in meters (used above)
geo_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500))
# Declare a list of the resampling array objects stored in the data directory
COGL_SMBO2_arrays <- grep('.Rdata', list.files(COGL_filePath, full.names = TRUE), value = TRUE, )
# # Decalre a list of NRMSE matrices, to store the results at each haplotype length
COGL_NRMSE_MatList <- list(Hap1=NA, Hap2=NA, Hap3=NA, Hap4=NA, Hap5=NA, 
                           Hap6=NA, Hap7=NA, Hap8=NA, Hap9=NA, Hap10=NA)

# ---- CALCULATIONS ----
# Loop through the list of arrays
for(i in 1:length(COGL_SMBO2_arrays)){
  # Read in the array
  COGL_SMBO2_array <- readRDS(COGL_SMBO2_arrays[[i]])
  # Build a dataframe from array values
  COGL_SMBO2_DF <- resample.array2dataframe(COGL_SMBO2_array)
  # Build a matrix to capture NRMSE values
  COGL_NRMSE_Mat <- matrix(NA, nrow=length(geo_buffSize), ncol=2)
  # The names of this matrix match the different parts of the dataframe names
  colnames(COGL_NRMSE_Mat) <- c('Geo_Buff','Eco_Buff')
  rownames(COGL_NRMSE_Mat) <- paste0(geo_buffSize/1000, 'km')
  # Loop through the dataframe columns. The first two columns are skipped, as they're sampleNumber and the
  # predictve variable (genetic coverages)
  for(j in 3:ncol(COGL_SMBO2_DF)){
    # Calculate NRMSE for the current column in the dataframe
    COGL_NRMSEvalue <- nrmse.func(COGL_SMBO2_DF[,j], pred = COGL_SMBO2_DF$Total)
    # Get the name of the current dataframe column
    dataName <- unlist(strsplit(names(COGL_SMBO2_DF)[[j]],'_'))
    # Match the data name to the relevant rows/columns of the receiving matrix
    matRow <- which(rownames(COGL_NRMSE_Mat) == dataName[[3]])
    matCol <- which(colnames(COGL_NRMSE_Mat) == paste0(dataName[[1]],'_',dataName[[2]]))
    # Locate the NRMSE value accordingly
    COGL_NRMSE_Mat[matRow,matCol] <- COGL_NRMSEvalue
  }
  cat(paste0('\n','%%% HAPLOTYPE LENGTH: ',i,' %%%','\n'))
  print(COGL_NRMSE_Mat)
  # Pass the current NRMSE matrix as an item in a list
  COGL_NRMSE_MatList[[i]] <- COGL_NRMSE_Mat
  # Store the matrix as a CSV to disk, if it isn't already written
  if(!file.exists(paste0(COGL_filePath, 'COGL_SMBO2_NRMSE_Hap-', i, '.csv'))){
    write.table(COGL_NRMSE_Mat,
                file=paste0(COGL_filePath, 'COGL_SMBO2_NRMSE_Hap-', i, '.csv'),
                sep=',')
  }
}

# ---- PLOTTING ----
# %%% PLOT HAPLOTYPE-WISE ----
# Specify plot colors
plotColors <- colorRampPalette(c("darkred","azure4","lightgray"))(129)
# Loop through the list of arrays
for(i in 1:length(COGL_SMBO2_arrays)){
  # Build a matrix of mean values (based on array; custom function returns a data.frame)
  COGL_SMBO2_meanValues <- as.matrix(meanArrayValues(readRDS(COGL_SMBO2_arrays[[i]])))
  # Subset the matrix of mean values according to each coverage type
  COGL_SMBO2_GeoBuffMeans <- COGL_SMBO2_meanValues[,grep('Geo_Buff',colnames(COGL_SMBO2_meanValues))]
  COGL_SMBO2_EcoBuffMeans <- COGL_SMBO2_meanValues[,grep('Eco_Buff',colnames(COGL_SMBO2_meanValues))]
  # Create colors based on the NRMSE values in matrix. Make all points transparent (alpha)
  GeoBuffCols <- 
    alpha(plotColors[as.numeric(cut(COGL_NRMSE_MatList[[i]][,1], breaks = length(plotColors)))], 0.15)
  EcoBuffCols <- 
    alpha(plotColors[as.numeric(cut(COGL_NRMSE_MatList[[i]][,3], breaks = length(plotColors)))], 0.15)
  # Determine each coverage type's optimal buffer size for current haplotype length, for plotting
  optBuff_GeoBuff <- names(which.min(COGL_NRMSE_MatList[[i]][,1]))
  optBuff_EcoBuff <- names(which.min(COGL_NRMSE_MatList[[i]][,3]))
  # For each color vector, decrease transparency of points corresponding to the lowest NRMSE
  GeoBuffCols[[which.min(COGL_NRMSE_MatList[[i]][,1])]] <- 
    alpha(GeoBuffCols[[which.min(COGL_NRMSE_MatList[[i]][,1])]], 0.65)
  EcoBuffCols[[which.min(COGL_NRMSE_MatList[[i]][,3])]] <- 
    alpha(EcoBuffCols[[which.min(COGL_NRMSE_MatList[[i]][,3])]], 0.65)
  
  # Use matplot to plot values for different coverages. Stack 3 plots vertically
  par(mfrow=c(3,1), mar = c(3,3.5,2,1)+0.1, cex.main=1.9)
  # GeoBuff
  matplot(COGL_SMBO2_GeoBuffMeans, ylim=c(0,100), col=GeoBuffCols, pch=16, 
          ylab='', xlab='')
  title(main=paste0('M. guttatus: Haplotype length: ',i))
  mtext(text='Number of individuals', side=1, line=2.5, at=140, cex=0.9)
  # Add points for genetic values, subtitle, optimal buffer size, and legend
  points(COGL_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
  mtext(text=paste0('*Optimal buffer size: ',optBuff_GeoBuff), side=1, line=-1.8, at=85, cex=1.1)
  legend(x=170, y=63, inset = 0.05, xpd=TRUE, cex=1.1, fill=c('darkred','darkgray','cyan4'), 
         legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Genetic values'),
         y.intersp = 0.75, border='white')
  # GeoSDM
  matplot(COGL_SMBO2_GeoSDMMeans, ylim=c(0,100), col=GeoBuffCols, pch=16, 
          ylab='', xlab='')
  mtext(text='Number of individuals', side=1, line=2.5, at=140, cex=0.9)
  # Add points for genetic values, subtitle, optimal buffer size, and legend
  points(COGL_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
  mtext(text=paste0('*Optimal buffer size: ',optBuff_GeoSDM), side=1, line=-1.8, at=85, cex=1.1)
  legend(x=170, y=63, inset = 0.05, xpd=TRUE, cex=1.1, fill=c('darkred','darkgray','cyan4'), 
         legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Genetic values'),
         y.intersp = 0.75, border='white')
  mtext(text='Coverage (%)', side=2, line=2.3, cex=1.1, srt=90)
  # EcoBuff
  matplot(COGL_SMBO2_EcoBuffMeans, ylim=c(0,100), col=EcoBuffCols, pch=16, 
          ylab='', xlab='')
  mtext(text='Number of individuals', side=1, line=2, at=140, cex=0.9)
  # Add points for genetic values, subtitle, optimal buffer size, and legend
  points(COGL_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
  mtext(text=paste0('*Optimal buffer size: ',optBuff_EcoBuff), side=1, line=-1.8, at=85, cex=1.1)
  legend(x=170, y=63, inset = 0.05, xpd=TRUE, cex=1.1, fill=c('darkred','darkgray','cyan4'), 
         legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Genetic values'),
         y.intersp = 0.75, border='white')
}

# %%% ALL PLOTS, 1 WINDOW ----
# Specify a graph layout, with coverage types as rows and haplotype lengths as columns
graphMat <- matrix(1:15,nrow=3,ncol=5)
layout(graphMat) ; par(mar = c(2,4,1,2)+0.1)
# Loop through the list of arrays
for(i in 1:length(COGL_SMBO2_arrays)){
  # Build a matrix of mean values (based on array; custom function returns a data.frame)
  COGL_SMBO2_meanValues <- as.matrix(meanArrayValues(readRDS(COGL_SMBO2_arrays[[i]])))
  # Subset the matrix of mean values according to each coverage type
  COGL_SMBO2_GeoBuffMeans <- COGL_SMBO2_meanValues[,grep('Geo_Buff',colnames(COGL_SMBO2_meanValues))]
  COGL_SMBO2_EcoBuffMeans <- COGL_SMBO2_meanValues[,grep('Eco_Buff',colnames(COGL_SMBO2_meanValues))]
  # Create colors based on the NRMSE values in matrix. Make all points transparent (alpha)
  GeoBuffCols <- 
    alpha(plotColors[as.numeric(cut(COGL_NRMSE_MatList[[i]][,1], breaks = length(plotColors)))], 0.15)
  EcoBuffCols <- 
    alpha(plotColors[as.numeric(cut(COGL_NRMSE_MatList[[i]][,3], breaks = length(plotColors)))], 0.15)
  # Determine each coverage type's optimal buffer size for current haplotype length, for plotting
  optBuff_GeoBuff <- names(which.min(COGL_NRMSE_MatList[[i]][,1]))
  optBuff_EcoBuff <- names(which.min(COGL_NRMSE_MatList[[i]][,3]))
  # For each color vector, decrease transparency of points corresponding to the lowest NRMSE
  GeoBuffCols[[which.min(COGL_NRMSE_MatList[[i]][,1])]] <- 
    alpha(GeoBuffCols[[which.min(COGL_NRMSE_MatList[[i]][,1])]], 0.65)
  EcoBuffCols[[which.min(COGL_NRMSE_MatList[[i]][,3])]] <- 
    alpha(EcoBuffCols[[which.min(COGL_NRMSE_MatList[[i]][,3])]], 0.65)
  
  # Adjust plotting, after first column of graphs
  if(i==2){
    par(mar = c(1,1,1,2)+0.1)
  }
  # GeoBuff
  matplot(COGL_SMBO2_GeoBuffMeans, ylim=c(0,100), col=GeoBuffCols, pch=16, 
          ylab='', xlab='')
  title(main=paste0('COGL: Hap ',i))
  if(i==1){
    mtext(text='Geo Buff', side=2, line=2.3, cex=1.1, srt=90)
  }
  # Add points for genetic values, subtitle, optimal buffer size, and legend
  points(COGL_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
  # Add points for genetic values, subtitle, optimal buffer size, and legend
  points(COGL_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
  # EcoBuff
  matplot(COGL_SMBO2_EcoBuffMeans, ylim=c(0,100), col=EcoBuffCols, pch=16, 
          ylab='', xlab='')
  # Add points for genetic values, subtitle, optimal buffer size, and legend
  points(COGL_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
  if(i==1){
    mtext(text='Eco Buff', side=2, line=2.3, cex=1.1, srt=90)
  }
}
