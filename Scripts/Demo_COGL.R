# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: CONRADINA GLABRA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script calculating the correlation between genetic, geographic, and ecological coverage
# for Conradina glabra. Uses data files from Lauren Eserman--a VCF containing filtered
# SNPs, for genetic data, and a CSV of coordinates. 

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
# Specify filepath for COGL geographic and genetic data
COGL_filePath <- paste0(GeoGenCorr_wd, 'Datasets/COGL/')

# ---- GENETIC MATRIX
# Read in the VCF of filtered SNPs provided by Lauren Eserman. These reads were trimmed, 
# and then a "de novo reference" was created by collecting 10 million reads in total
# evenly from all samples. Reads with counts fewer than 6 and more than 700 were excluded.
# Alleles with a population frequency of less than 3% were removed.
COGL_vcf <- 
  read.vcfR(file=paste0(COGL_filePath, 'Genetic/conradina.vcf'))
# Convert the vcf to a genind; the return.alleles TRUE value is suggested in the function's help file
# This genind file is made up of 564 individuals and 16,291 loci
COGL_genind <- vcfR2genind(COGL_vcf, sep = "/", return.alleles = TRUE)
# Rename individuals to match the sample names in coordinate CSV (below)
indNames(COGL_genind) <- gsub("_sorted", "", indNames(COGL_genind))
indNames(COGL_genind) <- paste0("Congla", indNames(COGL_genind))

# ---- COORDINATE POINTS
# CSV was provided by Lauren Eserman. Because individuals were sampled very close 
# together, coordinates were determined by obtaining coordinates for the center of a
# "clump" of individuals, and then individual coordinates were calculated based on the
# separation and direction of each individual from the center. Some samples have identical
# coordinate values.
# NOTE: these coordinates cannot be shared externally, due to the rare status of this species!
COGL_coordinates <- 
  read.csv(file=paste0(COGL_filePath, 'Geographic/COGL_coordinates.csv'), header = TRUE)
# Rename the columns of the geographic coordinates data.frame (because geo.compareBuff function expects certain strings)
colnames(COGL_coordinates)[2:3] <- c('decimalLatitude', 'decimalLongitude')
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- grabWorldAdmin(GeoGenCorr_wd = GeoGenCorr_wd, fileExtentsion = ".gpkg", overwrite = FALSE)
# Perform geographic filter on the admin layer. 
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = COGL_coordinates)
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage 
# (solely in the U.S.)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level4/us_eco_l4.shp')))
# Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
# exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
world_poly_clip_W <- wrap(world_poly_clip)
ecoregion_poly_W <- wrap(ecoregion_poly)

# ---- REMOVE DUPLICATE INDIVIDUALS ----
# Remove two duplicate samples (with "_rep" in sample name), leaving 562 samples
COGL_genind <- COGL_genind[-grep('_rep', indNames(COGL_genind)), drop=TRUE]
COGL_coordinates<- COGL_coordinates[-grep('_rep', COGL_coordinates$Sample.Name),]

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('COGL_coordinates','COGL_genind','num_reps','geo_buffSize','eco_buffSize',
                              'world_poly_clip_W','ecoregion_poly_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres', 
                              'eco.intersectBuff','eco.compareBuff','eco.totalEcoregionCount',
                              'gen.getAlleleCategories','gen.buildDistMat', 'gen.calcGenDistCov', 
                              'calculateCoverage', 'exSituResample.Par', 'geo.gen.Resample.Par'))
# Specify file path, for saving resampling array
arrayDir <- paste0(COGL_filePath, 'resamplingData/COGL_SMBO3_G2GE_5r_resampArr.Rdata')

# Run resampling (in parallel)
COGL_demoArray_Par <- 
  geo.gen.Resample.Par(genObj=COGL_genind, genDistFlag=TRUE, geoFlag=TRUE, coordPts=COGL_coordinates, 
                       geoBuff=geo_buffSize, boundary=world_poly_clip_W, ecoFlag=TRUE, 
                       ecoBuff=eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer='NA', 
                       reps=num_reps, arrayFilepath=arrayDir, cluster=cl)

# Close cores
stopCluster(cl)

# # Run resampling not in parallel (for function testing purposes)
# COGL_demoArray_IND <-
#   geo.gen.Resample(gen_obj = COGL_genind, geoFlag = TRUE, coordPts = COGL_coordinates,
#                    geoBuff = geo_buffSize, boundary = world_poly_clip, ecoFlag = FALSE, reps = 1)

# %%% ANALYZE DATA %%% ----
# Specify filepath for COGL geographic and genetic data, including resampling data
COGL_filePath <- paste0(GeoGenCorr_wd, 'Datasets/COGL/')
arrayDir <- paste0(COGL_filePath, 'resamplingData/COGL_50km_GE_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
COGL_array <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values
COGL_DF <- resample.array2dataframe(COGL_array)
# # Calculate Spearman's r for geographic coverage
# COGL_spearR_geo <- round(cor(COGL_DF$Geo, COGL_DF$Total, method = 'spearman'),3) ; COGL_spearR_geo
# Calculate normalized root mean square value
COGL_nrmse_geo <- nrmse_func(obs=COGL_DF$Geo, pred=COGL_DF$Total) ; COGL_nrmse_geo

# ---- PLOTTING ----
# ---- CALCULATE 95% MSSE AND AVERAGE VALUES
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(COGL_array, allValues = TRUE)
# Subset matrix of all average values to just Total allelic and geographic coverage
averageValueMat_TG <- averageValueMat[,c(1,6)]
# Calculate the absolute difference between genetic and geographic, and add this as a column to the data.frame
averageValueMat_TG <- cbind(averageValueMat_TG, abs(averageValueMat_TG$Total-averageValueMat_TG$Geo))
names(averageValueMat_TG) <- c(names(averageValueMat_TG)[1:2], "Difference")

# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue', 'purple')
plotColors <- alpha(plotColors, 0.45)
plotColors_Sub <- plotColors[-(2:5)]

# ---- CORRELATION PLOTS
par(mfrow=c(2,1), mar=c(4,4,3,2)+0.1)
# ---- GEOGRAPHIC-GENETIC
plot(averageValueMat_TG$Geo, averageValueMat_TG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     xlab='', ylab='', col=plotColors[[6]])
title(main='C. glabra: Geographic by genetic coverage', line=1.5)
mtext(text='562 Individuals; 1 km buffer; 5 replicates', side=3, line=0.1, cex=1.3)
mtext(text='Geographic coverage (%)', side=1, line=3, cex=1.2)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.2, srt=90)
# Add NRMSE values for each comparison
text(x = 76, y = 43, labels = paste0('NRMSE: ', COGL_nrmse_geo), col='darkblue', cex=1.2)
# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TG[,1:2], ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='C. glabra: Coverage values', line=1.5)
mtext(text='562 Individuals; 1 km buffer; 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.2)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.2, srt=90)
# Add legend
legend(x=350, y=80, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (1 km buffer)'),
       col=c('red', 'darkblue'), pch = c(20,20), cex=1.2, pt.cex = 2, bty='n',
       y.intersp = 0.8)

# %%%% SMBO3 ----
# Specify filepath for COGL geographic and genetic data, including resampling array
COGL_filePath <- paste0(GeoGenCorr_wd, 'Datasets/COGL/')
arrayDir <- paste0(COGL_filePath, 'resamplingData/COGL_SMBO3_G2GE_5r_resampArr.Rdata')
# Read in array
COGL_SMBO3_array <- readRDS(arrayDir)

# ---- CALCULATIONS ----
# Build a data.frame from array values
COGL_SMBO3_DF <- resample.array2dataframe(COGL_SMBO3_array)
# Build tables of NRSMSE values, calculated based on data.frame
COGL_NRMSE_Mat_CV <- buildNRMSEmatrix(resampDF=COGL_SMBO3_DF, genCovType='CV', sdmFlag=FALSE)
COGL_NRMSE_Mat_GD <- buildNRMSEmatrix(resampDF=COGL_SMBO3_DF, genCovType='GD', sdmFlag=FALSE)
# Combine the results of the NRMSE values calculated using allelic coverages and using
# genetic distances, and then rename the columns accordingly
COGL_NRMSE_Mat <- cbind(COGL_NRMSE_Mat_CV, COGL_NRMSE_Mat_GD)
# Store the matrix as a CSV to disk
write.table(COGL_NRMSE_Mat,
            file=paste0(COGL_filePath, 'resamplingData/COGL_SMBO3_NRMSE.csv'), sep=',')

# SMBO2: OPTIMAL BUFFER SIZES ----
# Read in COGL SMBO2 resampling array amd convert to data.frame
COGL_filePath <- paste0(GeoGenCorr_wd, 'Datasets/COGL/')
COGL_arrayDir <- paste0(COGL_filePath, 'resamplingData/COGL_SMBO2_GE_5r_resampArr.Rdata')
# From COGL resampling array, return a matrix of average coverage values for optimal buffer sizes
COGL_optCovMat <- extractOptCovs(COGL_arrayDir)
# Calculate MSSEs: minimum number of samples for 95% of each coverage type
COGL_Gen_MSSE <- min(which(COGL_optCovMat[,1] > 95)) ; COGL_Gen_MSSE
COGL_GeoBuff_MSSE <- min(which(COGL_optCovMat[,2] > 95)) ; COGL_GeoBuff_MSSE
COGL_Eco_MSSE <- min(which(COGL_optCovMat[,3] > 95)) ; COGL_Eco_MSSE

# PLOTTING
# Specify plot colors
plotColors <- c('red', 'darkblue','purple')
plotColors_Fade <- alpha(plotColors, c(0.45, rep(0.85, length(plotColors)-1)))
# Set plotting window to stack 3 graphs vertically
par(mfcol=c(2,1), oma=rep(0.2,4), mar=c(2,4,3,1)+0.1)
# Geo Buff
matplot(COGL_optCovMat[,c(1,2)], ylim=c(0,110), col=plotColors_Fade[c(1, 2)], pch=16, ylab='')
abline(h=95, col="black", lty=3)
abline(v=COGL_Gen_MSSE, col="red") ; abline(v=COGL_GeoBuff_MSSE, col="darkblue")
mtext(text=paste0('MSSE: ', COGL_Gen_MSSE), side=1, line=-1, at=COGL_Gen_MSSE+18, cex=0.8, col='red')
mtext(text=paste0(' MSSE: ', COGL_GeoBuff_MSSE), line=-1.5, side=1, at=COGL_GeoBuff_MSSE-20, cex=0.8, col='darkblue')
title('Conradina glabra: Coverages at Optimal Buffer Sizes', cex.sub=1.2, line = 2)
mtext(text='Geographic (Total Buffer): 15 km', side=3, at=80, cex=0.9)
# Legend
legend(x=240, y=95, xpd=TRUE, cex=1, pch=rep(19,3),
       col=c('red','darkblue','purple'),
       legend=c('Genetic', 'Geographic (Total Buffer)', 'Ecological'),
       y.intersp = 0.3, bty='n')
# Eco Buff
par(mar=c(3,4,2,1)+0.1)
matplot(COGL_optCovMat[,c(1,3)], ylim=c(0,110), col=plotColors_Fade[c(1, 3)], pch=16, ylab='')
abline(h=95, col="black", lty=3)
abline(v=COGL_Gen_MSSE, col="red") ; abline(v=COGL_Eco_MSSE, col="purple")
mtext(text=paste0('MSSE: ', COGL_Gen_MSSE), side=1, line=-1, at=COGL_Gen_MSSE+18, cex=0.8, col='red')
mtext(text=paste0(' MSSE: ', COGL_Eco_MSSE), line=-1.5, side=1, at=COGL_Eco_MSSE-14, cex=0.8, col='purple')
mtext(text='Ecological: 240 km', side=3, at=80, cex=0.9)
