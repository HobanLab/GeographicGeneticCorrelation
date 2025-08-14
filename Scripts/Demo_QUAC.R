# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: QUERCUS ACERIFOLIA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script demonstrating approach for calculating the correlation between genetic, geographic, 
# and ecological coverage. Uses a genind file from Quercus acerifolia (SNP loci, Complete dataset), 
# as well as a CSV containing sample names, latitudes and longitudes, to iteratively resample 
# wild points and measure coverage. Based on the parFlag value, the code is branched to 
# allow for procesing in parallel or not in parallel.

pacman::p_load(adegenet, terra, parallel, RColorBrewer, viridis, scales, vcfR, usedist, rnaturalearth)

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

# %%%% INDIVIDUAL-LEVEL GEOGRAPHIC COORDINATES %%%% ----
# In this analysis, we utilize a CSV file of lat/longs that specify the location of each individual

# ---- READ IN DATA ----
# Specify filepath for QUAC geographic and genetic data
QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')

# ---- GENETIC MATRIX
# Read in genind file: Reference aligment; R80, min-maf=0, first SNP/locus, only wild individuals,
# no Kessler individuals.
# Wild sample names/order must match those in the sample name column of the CSV (above)
QUAC_genind <- read.genepop(paste0(QUAC_filePath,'Genetic/QUAC_REF_Wild_R80_NOMAF_1SNP_NoK.gen'))
# Set all popnames for this genind object to 'wild'
pop(QUAC_genind) <- rep('wild', nInd(QUAC_genind))

# ---- GEOGRAPHIC/ECOLOGICAL DATA FILES
# Read in wild occurrence points. This CSV has 3 columns: sample name, latitude, and longitude. 
# The sample names (and order) have to match the sample names/order of the genind object 
# (rownams of the genetic matrix) read in below.
QUAC_coordinates <- read.csv(paste0(QUAC_filePath, 'Geographic/QUAC_coordinates.csv'), header=TRUE)
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- grabWorldAdmin(GeoGenCorr_wd = GeoGenCorr_wd, fileExtentsion = ".gpkg", overwrite = TRUE)
# Perform geographic filter on the admin layer
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip, wildPoints = QUAC_coordinates) 
# Read in raster data, for SDM
QUAC_sdm <- terra::rast(paste0(GeoGenCorr_wd,'/Datasets/QUAC/Geographic/QUAC_91inds_rast.tif'))
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage (solely in the U.S.)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level4/us_eco_l4.shp')))

# ---- PARALLELIZATION
# Flag for running resampling steps in parallel
parFlag <- FALSE

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
  invisible(clusterEvalQ(cl, library('ape')))
  # Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
  # exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
  world_poly_clip_W <- wrap(world_poly_clip)
  ecoregion_poly_W <- wrap(ecoregion_poly)
  QUAC_sdm_W <- wrap(QUAC_sdm)
}

# ---- RESAMPLING ----
if(parFlag==TRUE){
  # Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
  clusterExport(cl, varlist = c('QUAC_coordinates','QUAC_genind','num_reps','geo_buffSize','eco_buffSize',
                                'world_poly_clip_W','ecoregion_poly_W','QUAC_sdm_W'))
  # Export necessary functions (for calculating geographic and ecological coverage) to the cluster
  clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres', 
                                'eco.intersectBuff','eco.compareBuff','gen.getAlleleCategories', 
                                'gen.buildDistMat', 'gen.calcGenDistCov', 'eco.totalEcoregionCount',
                                'calculateCoverage','exSituResample.Par', 'geo.gen.Resample.Par'))
  # Specify file path, for saving resampling array
  # arrayDir <- paste0(QUAC_filePath, 'resamplingData/QUAC_SMBO3_G2G2E_5r_resampArr.Rdata')
  arrayDir <- paste0(QUAC_filePath, 'resamplingData/QUAC_TEST_5r_resampArr.Rdata')
  # Run resampling in parallel
  QUAC_demoArray_TEST <- 
    geo.gen.Resample.Par(genObj=QUAC_genind, genDistFlag=FALSE, geoFlag=TRUE, coordPts=QUAC_coordinates, 
                         SDMrast=QUAC_sdm_W, geoBuff=geo_buffSize, boundary=world_poly_clip_W, 
                         ecoFlag=TRUE, ecoBuff=eco_buffSize, ecoRegions=ecoregion_poly_W, ecoLayer='US',
                         reps=num_reps, arrayFilepath=arrayDir, cluster=cl)
  # Close cores
  stopCluster(cl)
} else {
  # Run resampling not in parallel (for function testing purposes)
  QUAC_demoArray_IND <-
    geo.gen.Resample(genObj=QUAC_genind,  genDistFlag=FALSE, geoFlag=TRUE, 
                     coordPts=QUAC_coordinates, SDMrast=NA, geoBuff=geo_buffSize, 
                     boundary=world_poly_clip, ecoFlag=FALSE, ecoBuff=eco_buffSize, 
                     ecoRegions=ecoregion_poly, ecoLayer='US', reps=1)
  arrayDir <- paste0(QUAC_filePath, 'resamplingData/QUAC_G2G2_1r_resampArr.Rdata')
  saveRDS(QUAC_demoArray_IND, arrayDir)
}

# %%% ANALYZE DATA %%% ----
QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')
arrayDir <- paste0(QUAC_filePath, 'resamplingData/QUAC_1kmIND_G_5r_resampArr.Rdata')
# Read in the resampling array .Rdata object, saved to disk
QUAC_demoArray_IND_Par <- readRDS(arrayDir)

# ---- CORRELATION ----
# Build a data.frame from array values
QUAC_DF <- resample.array2dataframe(QUAC_demoArray_IND_Par)
# Calculate normalized root mean square value
QUAC_nrmse_geo <- nrmse_func(obs=QUAC_DF$Geo, pred=QUAC_DF$Total) ; QUAC_nrmse_geo

# ---- PLOTTING ----
# ---- CALCULATE 95% MSSE AND AVERAGE VALUES
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic geographic proportions
averageValueMat <- meanArrayValues(QUAC_demoArray_IND_Par, allValues = TRUE)
# Subset matrix of all average values to just Total allelic and geographic coverage
averageValueMat_TG <- averageValueMat[,c(1,6)]
# Calculate the absolute difference between genetic and geographic, and add this as a column to the data.frame
averageValueMat_TG <- cbind(averageValueMat_TG, abs(averageValueMat_TG$Total-averageValueMat_TG$Geo))
names(averageValueMat_TG) <- c(names(averageValueMat_TG)[1:2], "Difference")

# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue', 'purple')
plotColors_Fade <- alpha(plotColors, 0.45)
plotColors_Sub <- plotColors_Fade[-(2:5)]

# ---- CORRELATION PLOTS
par(mfrow=c(2,1), mar=c(4,4,3,2)+0.1)
# ---- GEOGRAPHIC-GENETIC
plot(averageValueMat_TG$Geo, averageValueMat_TG$Total, pch=20, xlim=c(0,100), ylim=c(0,110),
     xlab='', ylab='', col=plotColors_Fade[[6]])
title(main='Q. acerifolia: Geographic by genetic coverage', line=1.5)
mtext(text='91 Individuals; 1 km buffer; 3 replicates', side=3, line=0.1, cex=1.3)
mtext(text='Geographic coverage (%)', side=1, line=3, cex=1.2)
mtext(text='Genetic coverage (%)', side=2, line=2.3, cex=1.2, srt=90)
# Add NRMSE values for each comparison
text(x = 76, y = 30, labels = paste0('NRMSE: ', QUAC_nrmse_geo), col='darkblue', cex=1.2)
# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_TG[,1:2], ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='')
# Add title and x-axis labels to the graph
title(main='Q. acerifolia: Coverage values', line=1.5)
mtext(text='91 Individuals; 1 km buffer; 3 replicates', side=3, line=0.3, cex=1.3)
mtext(text='Number of individuals', side=1, line=2.4, cex=1.2)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.2, srt=90)
# Add legend
legend(x=60, y=75, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (1 km buffer)'),
       col=c('red', 'darkblue'), pch = c(20,20), cex=1.2, pt.cex = 2, bty='n',
       y.intersp = 0.8)

# %%%% SMBO3 ----
# Specify filepath for QUAC geographic and genetic data, including resampling array
QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')
arrayDir <- paste0(QUAC_filePath, 'resamplingData/QUAC_SMBO3_G2G2E_5r_resampArr.Rdata')
# Read in array and build a data.frame of values
QUAC_SMBO3_array <- readRDS(arrayDir)

# ---- CALCULATIONS ----
# Build a data.frame from array values
QUAC_SMBO3_DF <- resample.array2dataframe(QUAC_SMBO3_array)
# Build tables of NRSMSE values, calculated based on data.frame
QUAC_NRMSE_Mat_CV <- buildNRMSEmatrix(resampDF=QUAC_SMBO3_DF, genCovType='CV', sdmFlag=TRUE)
QUAC_NRMSE_Mat_GD <- buildNRMSEmatrix(resampDF=QUAC_SMBO3_DF, genCovType='GD', sdmFlag=TRUE)
# Combine the results of the NRMSE values calculated using allelic coverages and using
# genetic distances, and then rename the columns accordingly
QUAC_NRMSE_Mat <- cbind(QUAC_NRMSE_Mat_CV, QUAC_NRMSE_Mat_GD)
# Store the matrix as a CSV to disk
write.table(QUAC_NRMSE_Mat,
            file=paste0(QUAC_filePath, 'resamplingData/QUAC_SMBO3_NRMSE.csv'), sep=',')

# ---- PLOTTING ----
# ALLELIC AND GENETIC DISTANCE COVERAGE
# Build a matrix of mean values (based on array; custom function returns a data.frame)
QUAC_SMBO3_meanValues <- as.matrix(meanArrayValues(QUAC_SMBO3_array))
# Subset the total mean Values to just Total, GenDist, and GeoBuff_0.5km
QUAC_SMBO3_Gen2GeoValues <- QUAC_SMBO3_meanValues[,1:3]
# Plot allelic coverage, genetic distance coverage, and geographic coverage at optimal buffer size
matplot(QUAC_SMBO3_Gen2GeoValues, ylim=c(0,100), col=c('red','darkred','purple'),
        pch=16, ylab='Coverage (%)', xlab='Number of individuals',
        main='Q. acerifolia: Allelic and Genetic Distance Coverage')
mtext(text='91 Individuals; 41 buffer sizes (0.5km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
legend(x=45, y=25, inset = 0.05,
       legend = c('Allelic coverage', 'Genetic distance coverage',
                  'Geographic coverage ("Total buffer", 0.5 km)'),
       col=c('red', 'darkred', 'purple'), pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.8)

# ALL SMBO3 RESULTS
# Specify plot colors
plotColors <- colorRampPalette(c("darkred","azure4","lightgray"))(129)
# Build a matrix of mean values (based on array; custom function returns a data.frame)
QUAC_SMBO3_meanValues <- as.matrix(meanArrayValues(QUAC_SMBO3_array))
# Subset the matrix of mean values according to each coverage type
QUAC_SMBO3_GeoBuffMeans <- QUAC_SMBO3_meanValues[,grep('Geo_Buff',colnames(QUAC_SMBO3_meanValues))]
QUAC_SMBO3_GeoSDMMeans <- QUAC_SMBO3_meanValues[,grep('Geo_SDM',colnames(QUAC_SMBO3_meanValues))]
QUAC_SMBO3_EcoBuffMeans <- QUAC_SMBO3_meanValues[,grep('Eco_Buff',colnames(QUAC_SMBO3_meanValues))]
# Create colors based on the NRMSE values in matrix (allelic coverages). Make all points transparent (alpha)
GeoBuffCols <-
  alpha(plotColors[as.numeric(cut(QUAC_NRMSE_Mat[,1], breaks = length(plotColors)))], 0.15)
GeoSDMCols <-
  alpha(plotColors[as.numeric(cut(QUAC_NRMSE_Mat[,2], breaks = length(plotColors)))], 0.15)
EcoBuffCols <-
  alpha(plotColors[as.numeric(cut(QUAC_NRMSE_Mat[,3], breaks = length(plotColors)))], 0.15)
# For each color vector, decrease transparency of points corresponding to the lowest NRMSE
GeoBuffCols[[which.min(QUAC_NRMSE_Mat[,1])]] <-
  alpha(GeoBuffCols[[which.min(QUAC_NRMSE_Mat[,1])]], 0.65)
GeoSDMCols[[which.min(QUAC_NRMSE_Mat[,2])]] <-
  alpha(GeoSDMCols[[which.min(QUAC_NRMSE_Mat[,2])]], 0.65)
EcoBuffCols[[which.min(QUAC_NRMSE_Mat[,3])]] <-
  alpha(EcoBuffCols[[which.min(QUAC_NRMSE_Mat[,3])]], 0.65)

# Use matplot to plot values for different coverages
# GeoBuff
matplot(QUAC_SMBO3_GeoBuffMeans, ylim=c(0,100), col=GeoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Q. acerifolia: Geographic Coverages (Total buffer)')
# Add points for allelic coverage values, genetic distance values, subtitle, optimal buffer size, and legend
points(QUAC_SMBO3_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
points(QUAC_SMBO3_meanValues[,2], col=alpha('blue', 0.55), pch=20)
mtext(text='91 Individuals; 41 buffer sizes (0.5km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal geo. buffer: 0.5 km', side=1, line=-10.0, at=77, cex=1.0)
legend(x=59, y=35, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4','blue'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Allelic coverage', 'Genetic distance coverage'),
       y.intersp = 0.75)
# GeoSDM
matplot(QUAC_SMBO3_GeoSDMMeans, ylim=c(0,100), col=GeoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Q. acerifolia: Geographic Coverages (SDM)')
# Add points for genetic values, subtitle, optimal buffer size, and legend
points(QUAC_SMBO3_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
points(QUAC_SMBO3_meanValues[,2], col=alpha('blue', 0.55), pch=20)
mtext(text='91 Individuals; 41 buffer sizes (0.5km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal geo. buffer: 10 km', side=1, line=-10.0, at=77, cex=1.0)
legend(x=59, y=35, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4','blue'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Allelic coverage', 'Genetic distance coverage'),
       y.intersp = 0.75)
# EcoBuff
matplot(QUAC_SMBO3_EcoBuffMeans, ylim=c(0,100), col=EcoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Q. acerifolia: Ecological Coverages')
# Add points for genetic values, subtitle, optimal buffer size, and legend
points(QUAC_SMBO3_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
points(QUAC_SMBO3_meanValues[,2], col=alpha('blue', 0.55), pch=20)
mtext(text='91 Individuals; 41 buffer sizes (0.5km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal ecological buffer size: 35 km', side=1, line=-10.0, at=77, cex=1.0)
legend(x=59, y=35, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4','blue'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Allelic coverage', 'Genetic distance coverage'),
       y.intersp = 0.75)
# Add arrows
arrows(x0=75, y0=40, x1=25, y1=80)

# Update SDM plots: remove data for first few buffer sizes (<10km), to make plots clearer
QUAC_SMBO3_GeoSDMMeans <- QUAC_SMBO3_GeoSDMMeans[,-(1:6)]
# Create new color vector
GeoSDMCols <-
  alpha(plotColors[as.numeric(cut(QUAC_NRMSE_Mat[,2], breaks = length(plotColors)))], 0.15)
GeoSDMCols[[which.min(QUAC_NRMSE_Mat[,2])]] <-
  alpha(GeoSDMCols[[which.min(QUAC_NRMSE_Mat[,2])]], 0.65)

# Call new plot
matplot(QUAC_SMBO3_GeoSDMMeans, ylim=c(0,100), col=GeoBuffCols, pch=16,
        ylab='Coverage (%)', xlab='Number of individuals',
        main='Q. acerifolia: Geographic Coverages (SDM)')
# Add points for genetic values, subtitle, optimal buffer size, and legend
points(QUAC_SMBO3_meanValues[,1], col=alpha('cyan4', 0.55), pch=20)
mtext(text='91 Individuals; 35 buffer sizes (10km -- 500km); 5 replicates', side=3, line=0.3, cex=1.3)
mtext(text='*Optimal geographic buffer size: 35 km', side=1, line=-1.7, at=200, cex=1.1)
legend(x=300, y=55, inset = 0.05, xpd=TRUE, cex=0.9, fill=c('darkred','darkgray','cyan4'),
       legend=c('Low NRMSE (better match)', 'High NRMSE (worse match)','Genetic values'),
       y.intersp = 0.75)

# SMBO2: OPTIMAL BUFFER SIZES ----
# Read in QUAC SMBO2 resampling array amd convert to data.frame
QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')
QUAC_arrayDir <- paste0(QUAC_filePath, 'resamplingData/QUAC_SMBO2_G2E_5r_resampArr.Rdata')
# From QUAC resampling array, return a matrix of average coverage values for optimal buffer sizes
QUAC_optCovMat <- extractOptCovs(QUAC_arrayDir)
# Calculate MSSEs: minimum number of samples for 95% of each coverage type
QUAC_Gen_MSSE <- min(which(QUAC_optCovMat[,1] > 95)) ; QUAC_Gen_MSSE
QUAC_GeoBuff_MSSE <- min(which(QUAC_optCovMat[,2] > 95)) ; QUAC_GeoBuff_MSSE
QUAC_GeoSDM_MSSE <- min(which(QUAC_optCovMat[,3] > 95)) ; QUAC_GeoSDM_MSSE # This will return Inf, as max coverage is 83%
QUAC_Eco_MSSE <- min(which(QUAC_optCovMat[,4] > 95)) ; QUAC_Eco_MSSE

# PLOTTING
# Specify plot colors
plotColors <- c('red', 'darkblue','darkorange3', 'purple')
plotColors_Fade <- alpha(plotColors, c(0.45, rep(0.85, length(plotColors)-1)))
# Set plotting window to stack 3 graphs vertically
par(mfcol=c(3,1), oma=rep(0.2,4), mar=c(2,4,3,1)+0.1)
# Geo Buff
matplot(QUAC_optCovMat[,c(1,2)], ylim=c(0,110), col=plotColors[c(1, 2)], pch=16, ylab='')
abline(h=95, col="black", lty=3)
abline(v=QUAC_Gen_MSSE, col="red") ; abline(v=QUAC_GeoBuff_MSSE, col="darkblue")
mtext(text=paste0('MSSE: ', QUAC_Gen_MSSE), side=1, line=-1, at=QUAC_Gen_MSSE+3, cex=0.8, col='red')
mtext(text=paste0(' MSSE: ', QUAC_GeoBuff_MSSE), line=-1.5, side=1, at=QUAC_GeoBuff_MSSE-3, cex=0.8, col='darkblue')
title('Quercus acerifolia: Coverages at Optimal Buffer Sizes', cex.sub=1.2, line = 2)
mtext(text='Geographic (Total Buffer): 0.5 km', side=3, at=15, cex=0.9)
# Geo SDM
par(mar=c(2,4,2,1)+0.1)
matplot(QUAC_optCovMat[,c(1,3)], ylim=c(0,110), col=plotColors[c(1, 3)], pch=16, ylab='')
abline(h=95, col="black", lty=3)
abline(v=QUAC_Gen_MSSE, col="red")
mtext(text=paste0('MSSE: ', QUAC_Gen_MSSE), side=1, line=-1, at=QUAC_Gen_MSSE+3, cex=0.8, col='red')
mtext(text='Geographic (SDM): 15 km', side=3, at=15, cex=0.9)
mtext(text="Coverage (%)", side=2, line=2.6, cex=1.2, srt=90)
# Legend
legend(x=45, y=95, xpd=TRUE, cex=1.4, pch=rep(19,4),
       col=c('red','darkblue','darkorange3', 'purple'),
       legend=c('Genetic', 'Geographic (Total Buffer)','Geographic (SDM)', 'Ecological'),
       y.intersp = 0.3, bty='n')
# Eco Buff
par(mar=c(3,4,2,1)+0.1)
matplot(QUAC_optCovMat[,c(1,4)], ylim=c(0,110), col=plotColors[c(1, 4)], pch=16, ylab='')
abline(h=95, col="black", lty=3)
abline(v=QUAC_Gen_MSSE, col="red") ; abline(v=QUAC_Eco_MSSE, col="purple")
mtext(text=paste0('MSSE: ', QUAC_Gen_MSSE), side=1, line=-1, at=QUAC_Gen_MSSE+3, cex=0.8, col='red')
mtext(text=paste0(' MSSE: ', QUAC_Eco_MSSE), line=-1.5, side=1, at=QUAC_Eco_MSSE-3, cex=0.8, col='purple')
mtext(text='Ecological: 120 km', side=3, at=15, cex=0.9)
