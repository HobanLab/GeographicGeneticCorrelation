# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: MIMULUS GUTTATUS USING HAPLOTYPES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage
# for Mimulus guttatus. This script explores the impact of using haplotypes of different 
# lengths on genetic resampling curves. To do this, resampling is performed multiple times 
# using different input files, with the results of each haplotype length stored to a 
# different array object.

pacman::p_load(adegenet, terra, parallel, RColorBrewer, viridis, scales, vcfR, usedist)

# Read in relevant functions
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/Code/'
setwd(GeoGenCorr_wd)
source('Scripts/functions_GeoGenCoverage.R')

# ---- VARIABLES ----
# Specify number of resampling replicates
num_reps <- 5
# ---- BUFFER SIZES
# Specify geographic buffer size in meters
geo_buffSize <- 1000*(c(2,3,4,5,seq(10,100,5),seq(110,250,10),500))
# Specify ecological buffer size in meters 
eco_buffSize <- 1000*(c(2,3,4,5,seq(10,100,5),seq(110,250,10),500))

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
  read.csv(file=paste0(MIGU_filePath, 'Geographic/MIGU_coordinates.csv'), header = TRUE)
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

# ---- GENETIC MATRICES
# For each haplotype length (ranging from a single SNP to 5 concatenated SNPs), the genetic
# data is stored as a STRUCTURE input file. Each STRUCTURE file is read in and converted into a
# genind object.
MIGU_haplos_filePath <- paste0(MIGU_filePath, 'Genetic/haplotypicData')
MIGU_haplos_input <- list.files(MIGU_haplos_filePath, full.names = TRUE)
# Declare a list of integers, storing the number of loci for each STRUCTURE input file
MIGU_lociLevels <- c(1498,608,326,190,115)

# Predeclare an empty list to store genind objects
MIGU_genList <- list()
# For each STRUCTURE input file, 
for(i in 1:length(MIGU_haplos_input)){
  # Read the file in and convert it into a genind object
  MIGU_genind_global <- read.structure(file=MIGU_haplos_input[[i]], n.ind=474, 
                                       n.loc=MIGU_lociLevels[[i]],onerowperind=TRUE, 
                                       col.lab=3, col.pop=1, col.others=2, row.marknames=0, 
                                       NA.char=-9, sep=" ", ask=TRUE, quiet=TRUE)
  
  # REMOVE INTRODUCED POPULATIONS: Subset global genind object to only contain individuals from native range. 
  # The 'drop' argument removes alleles no longer present in the dataset.
  MIGU_genList[[i]] <- MIGU_genind_global[MIGU_coordinates[,1], drop=TRUE]
  # Specify the population of every sample to be 'wild'
  pop(MIGU_genList[[i]]) <- rep('wild', nInd(MIGU_genList[[i]]))
}

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('MIGU_coordinates','MIGU_genList','num_reps','geo_buffSize', 'eco_buffSize',
                              'world_poly_clip_W', 'ecoregion_poly_W', 'MIGU_sdm_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres', 
                              'eco.intersectBuff','eco.compareBuff','gen.getAlleleCategories',
                              'eco.totalEcoregionCount','calculateCoverage','exSituResample.Par',
                              'geo.gen.Resample.Par'))
# Specify file paths, for saving resampling arrays 
MIGU_haplos_output <- dirname(gsub("Genetic", "resamplingData", MIGU_haplos_input))
for(i in 1:length(MIGU_haplos_output)){
    MIGU_haplos_output[[i]] <- paste0(MIGU_haplos_output[[i]],'/SMBO2_G2E/MIGU_SMBO2_G2E_5r_2-500_Hap-', i, '_resampArr.Rdata')
}

# Run resampling (in parallel). Use a loop to iterate through the different haplotype lengths
for(i in 1:length(MIGU_haplos_output)){
  print(paste0('%%% HAPLOTYPE LENGTH: ', i))
  # Run resampling
  MIGU_demoArray_Par <- 
    geo.gen.Resample.Par(genObj = MIGU_genList[[i]], geoFlag = TRUE, coordPts = MIGU_coordinates, 
                         geoBuff = geo_buffSize, SDMrast = MIGU_sdm_W, boundary=world_poly_clip_W, 
                         ecoFlag = TRUE, ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, 
                         ecoLayer = 'NA', reps = num_reps, arrayFilepath = MIGU_haplos_output[[i]], 
                         cluster = cl)
}
# Close cores
stopCluster(cl)

# %%% ANALYZE DATA %%% ----
# Specify filepath for MIGU geographic and genetic data, including resampling data
MIGU_filePath <- paste0(GeoGenCorr_wd, 'Datasets/MIGU/')
arrayDir <- paste0(MIGU_filePath, 'resamplingData/haplotypicData/SMBO2_G2E')
# Read in the resampling array .Rdata objects, saved to disk. Then, calculate the average
# value matrix (across resampling replicates) for each resampling array, and add that
# matrix as an item to a list.
MIGU_haplos_output <- list.files(arrayDir, full.names = TRUE); MIGU_haplos_averageValueMats <- list()
for(i in 1:length(MIGU_haplos_output)){
  MIGU_haplos_averageValueMats[[i]] <- meanArrayValues(readRDS(MIGU_haplos_output[[i]]))
}
# Build a summary matrix that consists of the Total allelic representation values for each
# haplotype dataset (5 columns), average Geographic representation values (1 column), and average
# Ecological representation values (1 column)
MIGU_haplos_summaryMat <- array(data=NA, dim=c(nrow(MIGU_haplos_averageValueMats[[1]]),7))
colnames(MIGU_haplos_summaryMat) <- c(paste0(rep('Total_Hap', 5),seq(1:5)),'Geo', 'Eco')
# Total allelic representation values: pull the "Total" columns from each average value matrix,
# and cbind these together
MIGU_haplos_summaryMat[,1:5] <- do.call(cbind,lapply(MIGU_haplos_averageValueMats, function(x) x$Total))
# For the geographic/ecological coverages: build a matrix of the coverage values
# from each average value matrix (using sapply). Then, calculate the means across the
# rows of this matrix (using rowMeans), and pass this to the last column of the summary matrix
MIGU_haplos_summaryMat[,6] <- rowMeans(sapply(MIGU_haplos_averageValueMats, function(df) df[, 2]))
MIGU_haplos_summaryMat[,7] <- rowMeans(sapply(MIGU_haplos_averageValueMats, function(df) df[, 3]))
# Convert the summary data into a data.frame
MIGU_haplos_summaryMat <- as.data.frame(MIGU_haplos_summaryMat)
# %%% NORMALIZED ROOT MEAN SQUARE ERROR: calculate the NRMSE for each haplotype length, for both
# geographic and ecological coverage
MIGU_NRMSE_Values <- matrix(nrow=length(MIGU_haplos_output), ncol=2)

for(i in 1:nrow(MIGU_NRMSE_Values)){
  MIGU_NRMSE_Values[i,1] <- nrmse_func(MIGU_haplos_summaryMat$Geo, pred=MIGU_haplos_summaryMat[,i])
  MIGU_NRMSE_Values[i,2] <- nrmse_func(MIGU_haplos_summaryMat$Eco, pred=MIGU_haplos_summaryMat[,i])
}

# ---- PLOTTING ----
hapColors <- c('gold2','orange','salmon','darkorange2','tomato3','darkblue','purple')
hapColors_Fade <- alpha(hapColors, 0.5)
legText <- c(paste0(rep('Hap. length: ', 5), seq(1:5)), 'Geographic coverage (50 km buffer)',
             'Ecological coverage (50 km buffer, EPA Level III)')
# ---- COVERAGE PLOTS
matplot(MIGU_haplos_summaryMat, ylim=c(0,110), col=hapColors_Fade, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='M. guttatus: Haplotypic Coverages', line=1.5)
mtext(text='255 Individuals; 50 km Geographic buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.2)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.2, srt=90)
# Add legend
legend(x=120, y=60, inset = 0.05, legend = legText, col=hapColors, pch = c(19,19),
       cex=1.2, pt.cex = 2, bty='n', y.intersp = 0.5)

# %%%% SMBO2: MULTIPLE BUFFER SIZES ----
# Specify filepath for SMBO2 resampling arrays
MIGU_filePath <- paste0(GeoGenCorr_wd, 'Datasets/MIGU/resamplingData/haplotypicData/SMBO2_G2E/')
# Specify geographic buffer size in meters (used above)
geo_buffSize <- 1000*(c(2,3,4,5,seq(10,100,5),seq(110,250,10),500))
# Declare a list of the resampling array objects stored in the data directory
MIGU_SMBO2_arrays <- grep('.Rdata', list.files(MIGU_filePath, full.names = TRUE), value = TRUE, )
# # Decalre a list of NRMSE matrices, to store the results at each haplotype length
MIGU_NRMSE_MatList <- list(Hap1=NA, Hap2=NA, Hap3=NA, Hap4=NA, Hap5=NA)

# ---- CALCULATIONS ----
# Loop through the list of arrays
for(i in 1:length(MIGU_SMBO2_arrays)){
  # Read in the array
  MIGU_SMBO2_array <- readRDS(MIGU_SMBO2_arrays[[i]])
  # Build a dataframe from array values
  MIGU_SMBO2_DF <- resample.array2dataframe(MIGU_SMBO2_array)
  # Build a matrix to capture NRMSE values
  MIGU_NRMSE_Mat <- matrix(NA, nrow=length(geo_buffSize), ncol=3)
  # The names of this matrix match the different parts of the dataframe names
  colnames(MIGU_NRMSE_Mat) <- c('Geo_Buff','Geo_SDM','Eco_Buff')
  rownames(MIGU_NRMSE_Mat) <- paste0(geo_buffSize/1000, 'km')
  # Loop through the dataframe columns. The first two columns are skipped, as they're sampleNumber and the
  # predictve variable (genetic coverages)
  for(j in 3:ncol(MIGU_SMBO2_DF)){
    # Calculate NRMSE for the current column in the dataframe
    MIGU_NRMSEvalue <- nrmse.func(MIGU_SMBO2_DF[,j], pred = MIGU_SMBO2_DF$Total)
    # Get the name of the current dataframe column
    dataName <- unlist(strsplit(names(MIGU_SMBO2_DF)[[j]],'_'))
    # Match the data name to the relevant rows/columns of the receiving matrix
    matRow <- which(rownames(MIGU_NRMSE_Mat) == dataName[[3]])
    matCol <- which(colnames(MIGU_NRMSE_Mat) == paste0(dataName[[1]],'_',dataName[[2]]))
    # Locate the NRMSE value accordingly
    MIGU_NRMSE_Mat[matRow,matCol] <- MIGU_NRMSEvalue
  }
  cat(paste0('\n','%%% HAPLOTYPE LENGTH: ',i,' %%%','\n'))
  print(MIGU_NRMSE_Mat)
  # Pass the current NRMSE matrix as an item in a list
  MIGU_NRMSE_MatList[[i]] <- MIGU_NRMSE_Mat
  # Store the matrix as a CSV to disk, if it isn't already written
  if(!file.exists(paste0(MIGU_filePath, 'MIGU_SMBO2_NRMSE_Hap-', i, '.csv'))){
    write.table(MIGU_NRMSE_Mat,
                file=paste0(MIGU_filePath, 'MIGU_SMBO2_NRMSE_Hap-', i, '.csv'))
  }
}

# ---- PLOTTING ----
# %%% PLOT HAPLOTYPE-WISE ----
# Specify plot colors
plotColors <- colorRampPalette(c("darkred","azure4","lightgray"))(129)
# Loop through the list of arrays
for(i in 1:length(MIGU_SMBO2_arrays)){
  # Build a matrix of mean values (based on array; custom function returns a data.frame)
  MIGU_SMBO2_meanValues <- as.matrix(meanArrayValues(readRDS(MIGU_SMBO2_arrays[[i]])))
  # Subset the matrix of mean values according to each coverage type
  MIGU_SMBO2_GeoBuffMeans <- MIGU_SMBO2_meanValues[,grep('Geo_Buff',colnames(MIGU_SMBO2_meanValues))]
  MIGU_SMBO2_GeoSDMMeans <- MIGU_SMBO2_meanValues[,grep('Geo_SDM',colnames(MIGU_SMBO2_meanValues))]
  MIGU_SMBO2_EcoBuffMeans <- MIGU_SMBO2_meanValues[,grep('Eco_Buff',colnames(MIGU_SMBO2_meanValues))]
  # Create colors based on the NRMSE values in matrix. Make all points transparent (alpha)
  GeoBuffCols <- 
    alpha(plotColors[as.numeric(cut(MIGU_NRMSE_MatList[[i]][,1], breaks = length(plotColors)))], 0.15)
  GeoSDMCols <- 
    alpha(plotColors[as.numeric(cut(MIGU_NRMSE_MatList[[i]][,2], breaks = length(plotColors)))], 0.15)
  EcoBuffCols <- 
    alpha(plotColors[as.numeric(cut(MIGU_NRMSE_MatList[[i]][,3], breaks = length(plotColors)))], 0.15)
  # Determine each coverage type's optimal buffer size for current haplotype length, for plotting
  optBuff_GeoBuff <- names(which.min(MIGU_NRMSE_MatList[[i]][,1]))
  optBuff_GeoSDM <- names(which.min(MIGU_NRMSE_MatList[[i]][,2]))
  optBuff_EcoBuff <- names(which.min(MIGU_NRMSE_MatList[[i]][,3]))
  # For each color vector, decrease transparency of points corresponding to the lowest NRMSE
  GeoBuffCols[[which.min(MIGU_NRMSE_MatList[[i]][,1])]] <- 
    alpha(GeoBuffCols[[which.min(MIGU_NRMSE_MatList[[i]][,1])]], 0.65)
  GeoSDMCols[[which.min(MIGU_NRMSE_MatList[[i]][,2])]] <- 
    alpha(GeoSDMCols[[which.min(MIGU_NRMSE_MatList[[i]][,2])]], 0.65)
  EcoBuffCols[[which.min(MIGU_NRMSE_MatList[[i]][,3])]] <- 
    alpha(EcoBuffCols[[which.min(MIGU_NRMSE_MatList[[i]][,3])]], 0.65)
  
  # Use matplot to plot values for different coverages. Stack 3 plots vertically
  par(mfrow=c(3,1), mar = c(3,3.5,2,1)+0.1, cex.main=1.9)
  # GeoBuff
  matplot(MIGU_SMBO2_GeoBuffMeans, ylim=c(0,100), col=GeoBuffCols, pch=16, 
          ylab='', xlab='')
  title(main=paste0('M. guttatus: Haplotype length: ',i))
  mtext(text='Number of individuals', side=1, line=2.5, at=140, cex=0.9)
  # Add points for genetic values, subtitle, optimal buffer size, and legend
  points(MIGU_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
  mtext(text=paste0('*Optimal buffer size: ',optBuff_GeoBuff), side=1, line=-1.8, at=85, cex=1.1)
  legend(x=170, y=63, inset = 0.05, xpd=TRUE, cex=1.1, fill=c('darkred','darkgray','cyan4'), 
         legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Genetic values'),
         y.intersp = 0.75, border='white')
  # GeoSDM
  matplot(MIGU_SMBO2_GeoSDMMeans, ylim=c(0,100), col=GeoBuffCols, pch=16, 
          ylab='', xlab='')
  mtext(text='Number of individuals', side=1, line=2.5, at=140, cex=0.9)
  # Add points for genetic values, subtitle, optimal buffer size, and legend
  points(MIGU_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
  mtext(text=paste0('*Optimal buffer size: ',optBuff_GeoSDM), side=1, line=-1.8, at=85, cex=1.1)
  legend(x=170, y=63, inset = 0.05, xpd=TRUE, cex=1.1, fill=c('darkred','darkgray','cyan4'), 
         legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Genetic values'),
         y.intersp = 0.75, border='white')
  mtext(text='Coverage (%)', side=2, line=2.3, cex=1.1, srt=90)
  # EcoBuff
  matplot(MIGU_SMBO2_EcoBuffMeans, ylim=c(0,100), col=EcoBuffCols, pch=16, 
          ylab='', xlab='')
  mtext(text='Number of individuals', side=1, line=2, at=140, cex=0.9)
  # Add points for genetic values, subtitle, optimal buffer size, and legend
  points(MIGU_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
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
for(i in 1:length(MIGU_SMBO2_arrays)){
  # Build a matrix of mean values (based on array; custom function returns a data.frame)
  MIGU_SMBO2_meanValues <- as.matrix(meanArrayValues(readRDS(MIGU_SMBO2_arrays[[i]])))
  # Subset the matrix of mean values according to each coverage type
  MIGU_SMBO2_GeoBuffMeans <- MIGU_SMBO2_meanValues[,grep('Geo_Buff',colnames(MIGU_SMBO2_meanValues))]
  MIGU_SMBO2_GeoSDMMeans <- MIGU_SMBO2_meanValues[,grep('Geo_SDM',colnames(MIGU_SMBO2_meanValues))]
  MIGU_SMBO2_EcoBuffMeans <- MIGU_SMBO2_meanValues[,grep('Eco_Buff',colnames(MIGU_SMBO2_meanValues))]
  # Create colors based on the NRMSE values in matrix. Make all points transparent (alpha)
  GeoBuffCols <- 
    alpha(plotColors[as.numeric(cut(MIGU_NRMSE_MatList[[i]][,1], breaks = length(plotColors)))], 0.15)
  GeoSDMCols <- 
    alpha(plotColors[as.numeric(cut(MIGU_NRMSE_MatList[[i]][,2], breaks = length(plotColors)))], 0.15)
  EcoBuffCols <- 
    alpha(plotColors[as.numeric(cut(MIGU_NRMSE_MatList[[i]][,3], breaks = length(plotColors)))], 0.15)
  # Determine each coverage type's optimal buffer size for current haplotype length, for plotting
  optBuff_GeoBuff <- names(which.min(MIGU_NRMSE_MatList[[i]][,1]))
  optBuff_GeoSDM <- names(which.min(MIGU_NRMSE_MatList[[i]][,2]))
  optBuff_EcoBuff <- names(which.min(MIGU_NRMSE_MatList[[i]][,3]))
  # For each color vector, decrease transparency of points corresponding to the lowest NRMSE
  GeoBuffCols[[which.min(MIGU_NRMSE_MatList[[i]][,1])]] <- 
    alpha(GeoBuffCols[[which.min(MIGU_NRMSE_MatList[[i]][,1])]], 0.65)
  GeoSDMCols[[which.min(MIGU_NRMSE_MatList[[i]][,2])]] <- 
    alpha(GeoSDMCols[[which.min(MIGU_NRMSE_MatList[[i]][,2])]], 0.65)
  EcoBuffCols[[which.min(MIGU_NRMSE_MatList[[i]][,3])]] <- 
    alpha(EcoBuffCols[[which.min(MIGU_NRMSE_MatList[[i]][,3])]], 0.65)
  
  # Adjust plotting, after first column of graphs
  if(i==2){
    par(mar = c(1,1,1,2)+0.1)
  }
  # GeoBuff
  matplot(MIGU_SMBO2_GeoBuffMeans, ylim=c(0,100), col=GeoBuffCols, pch=16, 
          ylab='', xlab='')
  title(main=paste0('MIGU: Hap ',i))
  if(i==1){
    mtext(text='Geo Buff', side=2, line=2.3, cex=1.1, srt=90)
  }
  # Add points for genetic values, subtitle, optimal buffer size, and legend
  points(MIGU_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
  # GeoSDM
  matplot(MIGU_SMBO2_GeoSDMMeans, ylim=c(0,100), col=GeoBuffCols, pch=16, 
          ylab='', xlab='')
  # Add points for genetic values, subtitle, optimal buffer size, and legend
  points(MIGU_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
  if(i==1){
    mtext(text='Geo SDM', side=2, line=2.3, cex=1.1, srt=90)
  }
  # EcoBuff
  matplot(MIGU_SMBO2_EcoBuffMeans, ylim=c(0,100), col=EcoBuffCols, pch=16, 
          ylab='', xlab='')
  # Add points for genetic values, subtitle, optimal buffer size, and legend
  points(MIGU_SMBO2_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
  if(i==1){
    mtext(text='Eco Buff', side=2, line=2.3, cex=1.1, srt=90)
  }
}
