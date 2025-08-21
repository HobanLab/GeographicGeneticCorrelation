# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION: POINT SUMMARY METRICS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script calculates geographic metrics based on the collection of coordinate points from 
# each dataset used in the geographic-genetic correlation analyses. The goal of this analysis 
# is to examine any possible trends between these point-based metrics and the optimal geographic
# buffer sizes.

# The key geo.calc.pointSummaries function used in this script is declared in the 
# functions_GeoGenCoverage.R script. 

# Load packages 
pacman::p_load(redlistr, spatialEco, sf, terra, tmap, dplyr, sfdep, sfheaders,
               deldir, readr, Hmisc, corrplot, parallel)
tmap_mode("view")

# Read in relevant functions
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/Code/'
setwd(GeoGenCorr_wd)
source('Scripts/functions_GeoGenCoverage.R')

# CALCULATING POINT SUMMARY VALUES ----
# Create a list of the relevant CSVs in the Datasets repository 
pointsData <- list.files(paste0(getwd(), "/Datasets"), pattern = ".csv", recursive = TRUE, full.names = TRUE)
pointsData <- pointsData[grep('Geographic', pointsData)]
# Select one coordinate file per species
amth <- pointsData[grepl(pattern ="AMTH", pointsData)][2] # Necessary due to 2nd (original) coordinate file
arth <- pointsData[grepl(pattern ="ARTH", pointsData)][2] # Necessary due to 2nd (original) coordinate file
cogl <- pointsData[grepl(pattern ="COGL", pointsData)]
hiwa <- pointsData[grepl(pattern ="HIWA", pointsData)][2] # Necessary due to 2nd (original) coordinate file
migu <- pointsData[grepl(pattern ="MIGU", pointsData)][2] # Necessary because of a 2nd file of global coordinates
pico <- pointsData[grepl(pattern ="PICO", pointsData)]
quac <- pointsData[grepl(pattern ="QUAC", pointsData)]
qulo <- pointsData[grepl(pattern ="QULO", pointsData)]
vila <- pointsData[grepl(pattern ="VILA", pointsData)]
yubr <- pointsData[grepl(pattern ="YUBR", pointsData)][2] # Necessary due to 2nd (original) coordinate file

# Build a list of the relevant files, and fields within those files, for each dataset.
# Note that the column names between different coordinate files differ
pointsDataList <- list(
  amth = read_csv(amth) |>
    dplyr::mutate(taxon = "AMTH")|>
    dplyr::select(taxon, lat = decimalLatitude, lon = decimalLongitude),
  arth = read_csv(arth) |>
      plyr::mutate(taxon = "ARTH")|>
      dplyr::select(taxon, lat = decimalLatitude, lon = decimalLongitude),
  cogl = read_csv(cogl) |> 
    dplyr::mutate(taxon = "COGL")|>
    dplyr::select(taxon, lat = Latitude , lon = Longitude),
  hiwa = read_csv(hiwa)|> 
    dplyr::mutate(taxon = "HIWA")|>
    dplyr::select(taxon, lat = decimalLatitude  , lon = decimalLongitude),
  migu = read_csv(migu) |> 
    dplyr::mutate(taxon = "MIGU")|>
    dplyr::select(taxon, lat = Latitude , lon = Longitude), 
  pico = read_csv(pico) |> 
    dplyr::mutate(taxon = "PICO")|>
    dplyr::select(taxon, lat = Latitude , lon = Longitude),
  quac = read_csv(quac)|> 
    dplyr::mutate(taxon = "QUAC")|>
    dplyr::select(taxon, lat = decimalLatitude  , lon = decimalLongitude), 
  qulo = read_csv(qulo)|> 
    dplyr::mutate(taxon = "QULO")|>
    dplyr::select(taxon, lat = decimalLatitude  , lon = decimalLongitude),
  vila = read_csv(vila)|> 
    dplyr::mutate(taxon = "VILA")|>
    dplyr::select(taxon, lat = decimalLatitude  , lon = decimalLongitude),
  yubr = read_csv(yubr)|>
    dplyr::mutate(taxon = "YUBR")|>
    dplyr::select(taxon, lat = decimalLatitude  , lon = decimalLongitude)
)
# Apply function which calculates multiple point summary metrics to list of datasets
# pointSummaries <- lapply(pointsDataList, geo.calc.pointSummaries)
pointSummaries <- readRDS(paste0(GeoGenCorr_wd,'Datasets/Outputs/pointSummariesList.Rdata'))
# Transform the point summary values into a list, where columns are the species and rows are the metrics
pointSummariesMat <- matrix(unlist(pointSummaries), ncol = length(pointSummaries), byrow = FALSE)
colnames(pointSummariesMat) <- toupper(names(pointSummaries))
rownames(pointSummariesMat) <- c("EOO", "AOO", "ANN", "VOR", "STD", "ELA", "ELP")
# Write the matrix of point values to disk
# write_csv(x = as.data.frame(pointSummariesMat), file = paste0(GeoGenCorr_wd,"Datasets/pointSummaryMeasures.csv"))

# EXTRACTING OPTIMAL BUFFER SIZES ----
# To measure any possible correlation between optimal buffer sizes and the point summary statistics,
# the optimal buffer size for each dataset (and for each relevant coverage type) needs to be appended
# to the matrix of point summary values.

# Specify the filepaths to the resampling array for each dataset
resampArrList <- list(
  AMTH=paste0(GeoGenCorr_wd, 'Datasets/AMTH/resamplingData/AMTH_SMBO2_GE_5r_resampArr.Rdata'),
  ARTH=paste0(GeoGenCorr_wd, 'Datasets/ARTH/resamplingData/ARTH_SMBO2_GE_5r_resampArr.Rdata'),
  COGL=paste0(GeoGenCorr_wd, 'Datasets/COGL/resamplingData/COGL_SMBO2_GE_5r_resampArr.Rdata'),
  HIWA=paste0(GeoGenCorr_wd, 'Datasets/HIWA/resamplingData/HIWA_SMBO2_GE_5r_resampArr.Rdata'),
  MIGU=paste0(GeoGenCorr_wd, 'Datasets/MIGU/resamplingData/SMBO2_G2E/MIGU_SMBO2_G2E_5r_resampArr.Rdata'),
  PICO=paste0(GeoGenCorr_wd, 'Datasets/PICO/resamplingData/SMBO2_G2E/PICO_SMBO2_G2E_5r_resampArr.Rdata'),
  QUAC=paste0(GeoGenCorr_wd, 'Datasets/QUAC/resamplingData/QUAC_SMBO2_G2E_5r_resampArr.Rdata'),
  QULO=paste0(GeoGenCorr_wd, 'Datasets/QULO/resamplingData/SMBO2/QULO_SMBO2_G2E_5r_resampArr.Rdata'),
  VILA=paste0(GeoGenCorr_wd, 'Datasets/VILA/resamplingData/VILA_SMBO2_5r_resampArr.Rdata'),
  YUBR=paste0(GeoGenCorr_wd, 'Datasets/YUBR/resamplingData/YUBR_SMBO2_G2E_resampArr.Rdata')
)
# Based on data in resampling arrays, extract the optimal buffer sizes for each species
optBuffs <- lapply(resampArrList, extractOptBuffs)
# For datasets without SDM values, add a column (in order to match dimensions with other datasets)
optBuffs$AMTH <- c(optBuffs$AMTH[[1]],NA,optBuffs$AMTH[[2]])
optBuffs$ARTH <- c(optBuffs$ARTH[[1]],NA,optBuffs$ARTH[[2]])
optBuffs$COGL <- c(optBuffs$COGL[[1]],NA,optBuffs$COGL[[2]])
optBuffs$HIWA <- c(optBuffs$HIWA[[1]],NA,optBuffs$HIWA[[2]])
optBuffs$VILA <- c(optBuffs$VILA[[1]],NA,optBuffs$VILA[[2]])
names(optBuffs$AMTH) <- names(optBuffs$ARTH) <- names(optBuffs$COGL)<- names(optBuffs$HIWA) <- 
  names(optBuffs$VILA) <- names(optBuffs$QULO)
# Convert the list of optimal buffer size values to a matrix
optBuffsMat <- matrix(unlist(optBuffs), ncol = length(optBuffs), byrow = FALSE)
colnames(optBuffsMat) <- names(optBuffs)
rownames(optBuffsMat) <- c('Opt. Geo. Buff', 'Opt. Geo. SDM', 'Opt. Eco.')
# Combine the point summary matrix to the optimal buffer size matrix. Transpose such that 
# rows are datasets and columns are summary metrics, and order by optimal GeoBuff size
SMBO_Mat <- t(rbind(pointSummariesMat, optBuffsMat))
SMBO_Mat <- SMBO_Mat[order(SMBO_Mat[,'Opt. Geo. Buff']),]

# BUILDING AND PLOTTING CORRELATION MATRICES ----
# Setting the correlation type to Spearman, since we don't know whether the relationship
# between the optimal buffer sizes and each point statistic is linear (and data likely
# isn't Normally distributed)
corType <- 'spearman'
# GEO/ECO COVERAGES
# Build a correlation matrix based off of values
corMat_SMBO <- rcorr(SMBO_Mat, type=corType)
# Replace NAs in diagonal of p-value matrix with 0s, to match dimensions
corMat_SMBO$P[which(is.na(corMat_SMBO$P))] <- 0
# CORRECTION FOR MULTIPLE TESTING
# Get the unique (non-NA, off-diagonal) p-values
pValues <- corMat_SMBO$P[upper.tri(corMat_SMBO$P)]
# Adjust p-values (using BH, Benjamini–Hochberg: tests are not independent, we want to
# reduce false positives but retain power...)
adj_pValues <- p.adjust(pValues, method = "BH")
# Put adjusted p-values back into a symmetric matrix
adj_pMat <- matrix(NA, nrow = ncol(corMat_SMBO$P), ncol = ncol(corMat_SMBO$P))
rownames(adj_pMat) <- colnames(adj_pMat) <- colnames(corMat_SMBO$P)
adj_pMat[upper.tri(adj_pMat)] <- adj_pValues
adj_pMat[lower.tri(adj_pMat)] <- t(adj_pMat)[lower.tri(adj_pMat)]
diag(adj_pMat) <- 0  # Set diagonal to 0, for plotting purposes
corMat_SMBO$P <- adj_pMat # Reassign corrected p values in correlation matrix
# Plot correlation matrix using corrplot. Label significant correlations using
# asterisks
corrplot(corMat_SMBO$r, type="upper", order="original", p.mat = corMat_SMBO$P, 
         sig.level = 0.01, insig = "label_sig", diag = FALSE)

# Write image to disc
imageOutDir <- 
  '/home/akoontz/Documents/GeoGenCorr/Documentation/Images/20250813_MANUSCRIPT_DRAFT4/corMat_adjPvalues.png'
png(filename=imageOutDir, width=900, height=760)
par(oma=c(0,0,3,0), mar=c(5,4,7,2)+0.1)
corrplot(corMat_SMBO$r, type="upper", order="original", p.mat = corMat_SMBO$P, 
         sig.level = 0.01, insig = "label_sig", diag = FALSE, cl.cex = 1.4, tl.cex=1.2)
title('Correlations: Spatial statistics', line = 5.6, cex.main=1.3)
dev.off() # Turn off plotting device

# # SDM COVERAGES
# # Build a correlation matrix based off of values
# corMat_SMBO_SDM <- rcorr(SMBO_SDM_Mat, type=corType)
# # Replace NAs in diagonal of p-value matrix with 0s, to match dimensions
# corMat_SMBO_SDM$P[which(is.na(corMat_SMBO_SDM$P))] <- 0
# # Plot correlation matrix using corrplot. Label significant correlations using
# # asterisks
# corrplot(corMat_SMBO_SDM$r, type="upper", order="original", p.mat = corMat_SMBO_SDM$P, 
#          sig.level = 0.01, insig = "label_sig", diag = FALSE)
# mtext('Spearman correlations: Points-based statistics and Geo/SDM/Eco coverages', side=3, line=1.2, adj=0.6, cex=1.2)

# PLOTTING OPTIMAL BUFFER SIZES VERSUS POINTS BASED METRICS
plot(SMBO_Mat[,'ANN'], SMBO_Mat[,'Opt_Geo-Buff'], pch=16, col='black',
     ylab='Optimal Geographic Buffer Sizes', xlab='Average Nearest Neighbor Values',
     main='SMBO2: Buffer sizes across ANN metrics', ylim=c(-10,600))

plot(SMBO_Mat[,'EOO'], SMBO_Mat[,'Opt_Geo-Buff'], pch=16, col='black',
     ylab='Optimal Geographic Buffer Sizes', xlab='Extent of Occurrence Values',
     main='SMBO2: Buffer sizes across EOO metrics')

plot(SMBO_Mat[,'StDevEllP'], SMBO_Mat[,'Opt_Geo-Buff'], pch=16, col='black',
     ylab='Optimal Geographic Buffer Sizes', xlab='Standard Deviation Ellipses Perimeter Values',
     main='SMBO2: Buffer sizes across St. Dev. Ellipse Perimeter metrics')
