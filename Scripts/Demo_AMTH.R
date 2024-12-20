# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: AMSONIA THARPII %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage
# for Amsonia tharpii. Uses data files provided by Dylan Cohen (at Chicago Botanic Garden):
# a VCF, for genetic data, and a spreadsheet that has coordinates for each population. The
# coordinate values need to be processed in order to build a data.frame with coordinates
# for each individual; these steps are taken below.

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
# Make sure libraries (adegenet, terra, etc.) are on cluster (but avoid printing output)
invisible(clusterEvalQ(cl, library('adegenet')))
invisible(clusterEvalQ(cl, library('terra')))
invisible(clusterEvalQ(cl, library('parallel')))
invisible(clusterEvalQ(cl, library('usedist')))
invisible(clusterEvalQ(cl, library('ape')))

# %%% CONDUCT RESAMPLING %%% ----
# ---- READ IN DATA ----
# Specify filepath for AMTH geographic and genetic data
AMTH_filePath <- paste0(GeoGenCorr_wd, 'Datasets/AMTH/')

# ---- GENETIC MATRIX
# Read in the VCF file for Amsonia tharpii using vcfR::read.vcfR
AMTH_vcf <- read.vcfR(paste0(AMTH_filePath,'Genetic/Tharpii_only.vcf'))
# Convert the vcf to a genind; the return.alleles TRUE value is suggested in the function's help file
AMTH_genind <- vcfR2genind(AMTH_vcf, return.alleles = TRUE)

# ---- GEOGRAPHIC/ECOLOGICAL DATA FILES
# The original coordinates file for this Amsonia dataset needs to be processed such that coordinates
# are listed for individuals, rather than populations.
# Check if the processed file (called AMTH_coordinates.csv) already exists; if not, then 
# run the necessary processing steps.
# NOTE: these coordinates cannot be shared externally, due to the rare status of this species!
if(file.exists(paste0(AMTH_filePath, 'Geographic/AMTH_coordinates.csv'))){
  # Read in the CSV of processed coordinates. The first column contains row numbers
  AMTH_coordinates <- read.csv(
    paste0(AMTH_filePath, 'Geographic/AMTH_coordinates.csv'), header=TRUE)[2:4]
} else {
  # The CSV of geographic coordinates contains the coordinates of each population for each 
  # Amsonia species analyzed (of which there are 3). The relevant populations (for Amsonia tharpii)
  # are at the top, so only the first 5 rows are extracted. This file has coordinates 
  # for each population; the population of each sample is represented by the letter at 
  # the front of each individual's name (see indNames(AMTH_genind))
  AMTH_popCoordinates <- read.csv2(paste0(AMTH_filePath, 'Geographic/AMTH_coordinates_Original.csv'),
                                   header = TRUE,sep = ',',nrows = 5)
  # Build a data.frame for individual coordinates, based on the sample names, and specify column names
  AMTH_coordinates <- data.frame(indNames(AMTH_genind))
  AMTH_coordinates <- cbind(AMTH_coordinates,NA,NA)
  colnames(AMTH_coordinates) <- c('sampleNames', 'decimalLatitude', 'decimalLongitude')
  # The coordinates file needs to list the coordinates for each individual (rather than each population)
  # To achieve this, detect the first letter of each sample name, and based on that letter, list
  # the coordinates for the corresponding population
  for(i in 1:length(indNames(AMTH_genind))){
    if(substring(AMTH_coordinates[i,1],1,1) == 'P'){
      AMTH_coordinates[i,2:3] <- AMTH_popCoordinates[1,4:5]
    }
    if(substring(AMTH_coordinates[i,1],1,1) == 'R'){
      AMTH_coordinates[i,2:3] <- AMTH_popCoordinates[2,4:5]
    }
    if(substring(AMTH_coordinates[i,1],1,1) == 'B'){
      AMTH_coordinates[i,2:3] <- AMTH_popCoordinates[3,4:5]
    }
    if(substring(AMTH_coordinates[i,1],1,1) == 'C'){
      AMTH_coordinates[i,2:3] <- AMTH_popCoordinates[4,4:5]
    }
    if(substring(AMTH_coordinates[i,1],1,1) == 'T'){
      AMTH_coordinates[i,2:3] <- AMTH_popCoordinates[5,4:5]
    }
  }
  # Ensure coordinate values are numeric, not strings
  AMTH_coordinates[,2] <- as.numeric(AMTH_coordinates[,2])
  AMTH_coordinates[,3] <- as.numeric(AMTH_coordinates[,3])
  # Write resulting coordinates data.frame as CSV to disk, for future runs
  write.csv(AMTH_coordinates, file=paste0(AMTH_filePath,'Geographic/AMTH_coordinates.csv'), 
            row.names = FALSE)
}
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/world_countries_10m/world_countries_10m.shp')))
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = AMTH_coordinates)
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage 
# (solely in the U.S.)
ecoregion_poly <-
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level4/us_eco_l4.shp')))
# Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
# exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
world_poly_clip_W <- wrap(world_poly_clip)
ecoregion_poly_W <- wrap(ecoregion_poly)

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('AMTH_coordinates','AMTH_genind','num_reps','geo_buffSize', 'eco_buffSize',
                              'world_poly_clip_W', 'ecoregion_poly_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres', 
                              'eco.intersectBuff','eco.compareBuff','gen.getAlleleCategories', 
                              'gen.buildDistMat', 'gen.calcGenDistCov', 'eco.totalEcoregionCount',
                              'calculateCoverage','exSituResample.Par', 'geo.gen.Resample.Par'))
# Specify file path, for saving resampling array
arrayDir <- paste0(AMTH_filePath, 'resamplingData/AMTH_SMBO3_G2GE_5r_resampArr.Rdata')

# Run resampling (in parallel)
AMTH_demoArray_Par <- 
  geo.gen.Resample.Par(genObj=AMTH_genind, genDistFlag=TRUE, geoFlag=TRUE, coordPts=AMTH_coordinates, 
                       geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag=TRUE, ecoBuff=eco_buffSize, 
                       ecoRegions=ecoregion_poly_W, ecoLayer='US', reps=num_reps, arrayFilepath=arrayDir, cluster=cl)
# Close cores
stopCluster(cl)

# Run resampling not in parallel (for function testing purposes)
AMTH_demoArray_IND <-
  geo.gen.Resample(genObj = AMTH_genind, genDistFlag=TRUE, geoFlag = TRUE, coordPts = AMTH_coordinates, 
                   geoBuff = geo_buffSize, boundary = world_poly_clip, ecoFlag = FALSE, 
                   ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly, ecoLayer = "US", reps = 1)

# # %%% ANALYZE DATA %%% ----
# # Specify filepath for AMTH geographic and genetic data, including resampling array
# AMTH_filePath <- paste0(GeoGenCorr_wd, 'Datasets/AMTH/')
# arrayDir <- paste0(AMTH_filePath, 'resamplingData/AMTH_1km_GE_5r_resampArr.Rdata')
# # Read in the resampling array .Rdata object, saved to disk
# AMTH_demoArray_Par <- readRDS(arrayDir)
# 
# # ---- CORRELATION ----
# # Build a data.frame from array values
# AMTH_DF <- resample.array2dataframe(AMTH_demoArray_Par)
# # Calculate normalized root mean square value
# AMTH_nrmse_geo <- nrmse_func(obs=AMTH_DF$Geo, pred=AMTH_DF$Total) ; AMTH_nrmse_geo
# AMTH_nrmse_eco <- nrmse_func(obs=AMTH_DF$Eco, pred=AMTH_DF$Total) ; AMTH_nrmse_eco
# 
# # ---- PLOTTING ----
# # Generate the average values (across replicates) for all proportions
# # This function has default arguments for returning just Total allelic geographic proportions
# averageValueMat <- meanArrayValues(AMTH_demoArray_Par, allValues = TRUE)
# # Subset matrix of all average values to just Total allelic, geographic, and ecological coverage
# averageValueMat_TEG <- averageValueMat[,c(1,6,7)]
# # Calculate the absolute difference between genetic and geographic/ecological, and add to data.frame
# averageValueMat_TEG <- cbind(averageValueMat_TEG, abs(averageValueMat_TEG$Total-averageValueMat_TEG$Geo))
# averageValueMat_TEG <- cbind(averageValueMat_TEG, abs(averageValueMat_TEG$Total-averageValueMat_TEG$Eco))
# names(averageValueMat_TEG) <- c(names(averageValueMat_TEG)[1:3], 'Geo_Difference', 'Eco_Difference')
# 
# # Specify plot colors
# plotColors <- c('red','red4','darkorange3','coral','darkblue', 'purple')
# plotColors_Fade <- alpha(plotColors, 0.65)
# plotColors_Sub <- plotColors_Fade[-(2:4)]
# # Two plots in a single window
# par(mfrow=c(2,1))
# # ---- CORRELATION PLOTS
# plot(averageValueMat_TEG$Geo, averageValueMat_TEG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
#      main='A. tharpii: Geographic by genetic coverage',xlab='', ylab='', col=plotColors_Fade[[5]])
# mtext(text='140 Individuals; 1 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
# mtext(text='Geographic/Ecological coverage (%)', side=1, line=3, cex=1.6)
# mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# # Add points for ecological coverage
# points(x=averageValueMat$Eco, y=averageValueMat$Total, pch=20, col=plotColors_Fade[[6]])
# # Add NRMSE values for each comparison
# text(x = 76, y = 35, labels = paste0('NRMSE: ', AMTH_nrmse_geo), col='darkblue', cex=0.9)
# text(x = 76, y = 20, labels = paste0('NRMSE: ', AMTH_nrmse_eco), col='purple', cex=0.9)
# # Add legend
# legend(x=58, y=162, inset = 0.05, xpd=TRUE,
#        legend = c('Geographic', 'Ecological'),
#        col=c('darkblue', 'purple'), pch = c(20,20), cex=0.9, pt.cex = 2, bty='n',
#        y.intersp = 0.08)
# # ---- COVERAGE PLOTS
# # Use the matplot function to plot the matrix of average values, with specified settings
# matplot(averageValueMat_TEG[,1:3], ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# # Add title and x-axis labels to the graph
# title(main='A. tharpii: Geo-Eco-Gen Coverage', line=1.5)
# mtext(text='140 Individuals; 1 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
# mtext(text='Number of individuals', side=1, line=2.4, cex=1.6)
# mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# # Add legend
# legend(x=85, y=180, inset = 0.05,
#        legend = c('Genetic coverage', 'Geographic coverage (1 km buffer)',
#                   'Ecological coverage (1 km buffer, EPA Level IV)'),
#        col=c('red', 'darkblue', 'purple'), pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n',
#        y.intersp = 0.08)
# 
# # %%%% SMBO: MULTIPLE BUFFER SIZES ----
# # Specify filepath for AMTH geographic and genetic data, including resampling array
# AMTH_filePath <- paste0(GeoGenCorr_wd, 'Datasets/AMTH/')
# arrayDir <- paste0(AMTH_filePath, 'resamplingData/AMTH_SMBO2_GE_5r_resampArr.Rdata')
# # Read in array and build a data.frame of values
# AMTH_SMBO2_array <- readRDS(arrayDir)
# # Specify geographic buffer size in meters (used above)
# geo_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500))
# 
# # ---- CALCULATIONS ----
# # Build a data.frame from array values
# AMTH_SMBO2_DF <- resample.array2dataframe(AMTH_SMBO2_array)
# # Build a matrix to capture NRMSE values
# AMTH_NRMSE_Mat <- matrix(NA, nrow=length(geo_buffSize), ncol=2)
# # The names of this matrix match the different parts of the dataframe names
# colnames(AMTH_NRMSE_Mat) <- c('Geo_Buff','Eco_Buff')
# rownames(AMTH_NRMSE_Mat) <- paste0(geo_buffSize/1000, 'km')
# # Loop through the dataframe columns. The first two columns are skipped, as they're sampleNumber and the
# # predictve variable (genetic coverages)
# for(i in 3:ncol(AMTH_SMBO2_DF)){
#   # Calculate NRMSE for the current column in the dataframe
#   AMTH_NRMSEvalue <- nrmse.func(AMTH_SMBO2_DF[,i], pred = AMTH_SMBO2_DF$Total)
#   # Get the name of the current dataframe column
#   dataName <- unlist(strsplit(names(AMTH_SMBO2_DF)[[i]],'_'))
#   # Match the data name to the relevant rows/columns of the receiving matrix
#   matRow <- which(rownames(AMTH_NRMSE_Mat) == dataName[[3]])
#   matCol <- which(colnames(AMTH_NRMSE_Mat) == paste0(dataName[[1]],'_',dataName[[2]]))
#   # Locate the NRMSE value accordingly
#   AMTH_NRMSE_Mat[matRow,matCol] <- AMTH_NRMSEvalue
# }
# print(AMTH_NRMSE_Mat)
# # Store the matrix as a CSV to disk
# write.table(AMTH_NRMSE_Mat,
#             file=paste0(AMTH_filePath, 'resamplingData/AMTH_SMBO2_NRMSE.csv'), sep=',')
# 
# # SMBO2: OPTIMAL BUFFER SIZES ----
# # Read in AMTH SMBO2 resampling array amd convert to data.frame
# AMTH_filePath <- paste0(GeoGenCorr_wd, 'Datasets/AMTH/')
# AMTH_arrayDir <- paste0(AMTH_filePath, 'resamplingData/AMTH_SMBO2_GE_5r_resampArr.Rdata')
# # From AMTH resampling array, return a matrix of average coverage values for optimal buffer sizes
# AMTH_optCovMat <- extractOptCovs(AMTH_arrayDir)
# # Calculate MSSEs: minimum number of samples for 95% of each coverage type
# AMTH_Gen_MSSE <- min(which(AMTH_optCovMat[,1] > 95)) ; AMTH_Gen_MSSE
# AMTH_GeoBuff_MSSE <- min(which(AMTH_optCovMat[,2] > 95)) ; AMTH_GeoBuff_MSSE
# AMTH_Eco_MSSE <- min(which(AMTH_optCovMat[,3] > 95)) ; AMTH_Eco_MSSE
