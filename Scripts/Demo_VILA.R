# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: VITIS LABRUSCA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage
# for Vitis labrusca. Uses data files provided by Dr. Zoe Migicovsky

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

# %%% CONDUCT RESAMPLING %%% ----
# ---- READ IN DATA ----
# Specify filepath for VILA geographic and genetic data
VILA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/VILA/')

# ---- GEOGRAPHIC/ECOLOGICAL DATA FILES
# A text file was provided which specifies the sample name, latitude, longitude, 
# source, and species for each individual. Read this in (using read.delim for TSV)
VILA_coordinates <- 
  read.delim(paste0(VILA_filePath, 'Geographic/VILA_coordinates.tsv'), header = TRUE)
# Subset coordinate dataframe to strictly wild samples by removing germplasm individuals
VILA_coordinates <- VILA_coordinates[which(VILA_coordinates$POPTYPE=='Wild'),]
# Remove unnecessary columns (species, state, and poptype), and rename columns
VILA_coordinates <- VILA_coordinates[,c(2,4,5)]
colnames(VILA_coordinates) <- c('sampleID','decimalLatitude','decimalLongitude')
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
# Read in the VCF file for Amsonia tharpii using vcfR::read.vcfR
VILA_vcf <- read.vcfR(paste0(VILA_filePath,'Genetic/VILA.vcf'))
# Convert the vcf to a genind; return.alleles TRUE value suggested in help file
VILA_genind <- vcfR2genind(VILA_vcf, return.alleles = TRUE)
# Individual names are duplicated in genind object (such that 104_5 is called 
# 104_5_104_5). To address this, first, use the sub command to replace names
indNames(VILA_genind) <- sub("(.*?_.*?)_.*", "\\1", indNames(VILA_genind))
# Match the individuals in the genind file to those in the coordinate dataframe.
# This essentially drops germplasm individuals
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
arrayDir <- paste0(VILA_filePath, 'resamplingData/VILA_SMBO3_5r_resampArr.Rdata')

# Run resampling (in parallel)
VILA_demoArray_Par <- 
  geo.gen.Resample.Par(genObj=VILA_genind, geoFlag=TRUE, coordPts=VILA_coordinates, 
                       geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag=TRUE, ecoBuff=eco_buffSize, 
                       ecoRegions=ecoregion_poly_W, ecoLayer='US', reps=num_reps, arrayFilepath=arrayDir, cluster=cl)
# Close cores
stopCluster(cl)

# # Run resampling not in parallel (for function testing purposes)
# VILA_demoArray_IND <-
#   geo.gen.Resample(genObj = VILA_genind, geoFlag = TRUE, coordPts = VILA_coordinates, 
#                    geoBuff = geo_buffSize, boundary = world_poly_clip, ecoFlag = FALSE, 
#                    ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly, ecoLayer = "US", reps = 1)

# # %%% ANALYZE DATA %%% ----
# # Specify filepath for VILA geographic and genetic data, including resampling array
# VILA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/VILA/')
# arrayDir <- paste0(VILA_filePath, 'resamplingData/VILA_1km_GE_5r_resampArr.Rdata')
# # Read in the resampling array .Rdata object, saved to disk
# VILA_demoArray_Par <- readRDS(arrayDir)
# 
# # ---- CORRELATION ----
# # Build a data.frame from array values
# VILA_DF <- resample.array2dataframe(VILA_demoArray_Par)
# # Calculate normalized root mean square value
# VILA_nrmse_geo <- nrmse_func(obs=VILA_DF$Geo, pred=VILA_DF$Total) ; VILA_nrmse_geo
# VILA_nrmse_eco <- nrmse_func(obs=VILA_DF$Eco, pred=VILA_DF$Total) ; VILA_nrmse_eco
# 
# # ---- PLOTTING ----
# # Generate the average values (across replicates) for all proportions
# # This function has default arguments for returning just Total allelic geographic proportions
# averageValueMat <- meanArrayValues(VILA_demoArray_Par, allValues = TRUE)
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
# text(x = 76, y = 35, labels = paste0('NRMSE: ', VILA_nrmse_geo), col='darkblue', cex=0.9)
# text(x = 76, y = 20, labels = paste0('NRMSE: ', VILA_nrmse_eco), col='purple', cex=0.9)
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
# # Specify filepath for VILA geographic and genetic data, including resampling array
# VILA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/VILA/')
# arrayDir <- paste0(VILA_filePath, 'resamplingData/VILA_SMBO2_GE_5r_resampArr.Rdata')
# # Read in array and build a data.frame of values
# VILA_SMBO2_array <- readRDS(arrayDir)
# # Specify geographic buffer size in meters (used above)
# geo_buffSize <- 1000*(c(0.5,1,2,3,4,5,seq(10,100,5),seq(110,250,10),500))
# 
# # ---- CALCULATIONS ----
# # Build a data.frame from array values
# VILA_SMBO2_DF <- resample.array2dataframe(VILA_SMBO2_array)
# # Build a matrix to capture NRMSE values
# VILA_NRMSE_Mat <- matrix(NA, nrow=length(geo_buffSize), ncol=2)
# # The names of this matrix match the different parts of the dataframe names
# colnames(VILA_NRMSE_Mat) <- c('Geo_Buff','Eco_Buff')
# rownames(VILA_NRMSE_Mat) <- paste0(geo_buffSize/1000, 'km')
# # Loop through the dataframe columns. The first two columns are skipped, as they're sampleNumber and the
# # predictve variable (genetic coverages)
# for(i in 3:ncol(VILA_SMBO2_DF)){
#   # Calculate NRMSE for the current column in the dataframe
#   VILA_NRMSEvalue <- nrmse.func(VILA_SMBO2_DF[,i], pred = VILA_SMBO2_DF$Total)
#   # Get the name of the current dataframe column
#   dataName <- unlist(strsplit(names(VILA_SMBO2_DF)[[i]],'_'))
#   # Match the data name to the relevant rows/columns of the receiving matrix
#   matRow <- which(rownames(VILA_NRMSE_Mat) == dataName[[3]])
#   matCol <- which(colnames(VILA_NRMSE_Mat) == paste0(dataName[[1]],'_',dataName[[2]]))
#   # Locate the NRMSE value accordingly
#   VILA_NRMSE_Mat[matRow,matCol] <- VILA_NRMSEvalue
# }
# print(VILA_NRMSE_Mat)
# # Store the matrix as a CSV to disk
# write.table(VILA_NRMSE_Mat,
#             file=paste0(VILA_filePath, 'resamplingData/VILA_SMBO2_NRMSE.csv'), sep=',')
