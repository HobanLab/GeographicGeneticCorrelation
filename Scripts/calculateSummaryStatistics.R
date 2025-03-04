# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION: CALCULATE SUMMARY STATISTICS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script analyzes the SMBO2 results (resampling arrays) to build matrices of correlation metrics, measuring
# the relationship between the genetic coverage values and the geographic and ecological coverage values. It does
# this by calculating NRMSE values (rounded to 2 different values of decimal points), Spearman's rho, and Pearson's r
# (the last of these we don't really utilize). 

# Some custom functions, used for calculating MSSE values across datasets and for average Spearman rho values, are
# declared up front. The end of the script contains commands for plotting a matrix of Spearman rho correlation values.
pacman::p_load(adegenet, terra, parallel, RColorBrewer, scales, vcfR, usedist, DescTools, ggplot2, reshape2)

# Read in relevant functions
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/Code/'
setwd(GeoGenCorr_wd)
source('Scripts/functions_GeoGenCoverage.R')
# Specify the directory to save outputs to
outputDir <- paste0(GeoGenCorr_wd,'Datasets/Outputs/')

# %%%% FUNCTIONS %%% ----
# Declare a function that, given a resampling array filepath, will return a vector of MSSE
# values for a specified threshold size.
calcMSSEs <- function(arrayDir='~/resamp.Rdata', thresh=95){
  # Read in the resampling array
  resampArr <- readRDS(arrayDir)
  # Build a matrix of mean values from the resampling array
  meanValuesMat <- meanArrayValues(resampArr, allValues=FALSE)
  # Use an apply statement to calculate the 95% minimum sample size. The lambda function 
  # calculates the row of the array which is greater than 95; the +1 adds one to the value,
  # since the first row in the mean values array corresponds to 2 samples
  MSSEs <- apply(meanValuesMat, 2, function(x) which(x >= thresh)[1]+1)
  # Return MSSE values
  return(MSSEs)
}

# Function which, given a vector of correlation values, will calculate the average correlation value. Because correlation
# coefficients are not additive, this requires converting the coefficient values using Fisher's z transform, then averaging,
# and then performing the inverse Fisher z transform.
averageCorrValue <- function(corrValues){
  # Perform Fisher's z transform on correlation coefficients. Then, average and perform the inverse transform
  # Remove NA values before calculating the mean
  averages <- FisherZInv(mean(FisherZ(corrValues), na.rm=TRUE))
  return(averages)
}

# %%% READ IN DATASETS %%% ----
# Specify relevant resampling arrays, for each dataset
QUAC_arrayDir <- 
  paste0(GeoGenCorr_wd, 'Datasets/QUAC/resamplingData/QUAC_SMBO2_G2E_5r_resampArr.Rdata')
YUBR_arrayDir <- 
  paste0(GeoGenCorr_wd, 'Datasets/YUBR/resamplingData/YUBR_SMBO2_G2E_resampArr.Rdata')
COGL_arrayDir <-
  paste0(GeoGenCorr_wd, 'Datasets/COGL/resamplingData/COGL_SMBO2_GE_5r_resampArr.Rdata')
QULO_arrayDir <- 
  paste0(GeoGenCorr_wd, 'Datasets/QULO/resamplingData/SMBO2/QULO_SMBO2_G2E_5r_resampArr.Rdata')
PICO_arrayDir <- 
  paste0(GeoGenCorr_wd, 'Datasets/PICO/resamplingData/SMBO2_G2E/PICO_SMBO2_G2E_5r_resampArr.Rdata')
MIGU_arrayDir <- 
  paste0(GeoGenCorr_wd, 'Datasets/MIGU/resamplingData/SMBO2_G2E/MIGU_SMBO2_G2E_5r_resampArr.Rdata')
AMTH_arrayDir <-
  paste0(GeoGenCorr_wd, 'Datasets/AMTH/resamplingData/AMTH_SMBO2_GE_5r_resampArr.Rdata')
HIWA_arrayDir <-
  paste0(GeoGenCorr_wd, 'Datasets/HIWA/resamplingData/HIWA_SMBO2_GE_5r_resampArr.Rdata')
ARTH_arrayDir <- 
  paste0(GeoGenCorr_wd, 'Datasets/ARTH/resamplingData/ARTH_SMBO2_GE_5r_resampArr.Rdata')
VILA_arrayDir <- 
  paste0(GeoGenCorr_wd, 'Datasets/VILA/resamplingData/VILA_SMBO2_5r_resampArr.Rdata')

# Build a single vector of relevant array directories
SMBO2_values <- c(QUAC_arrayDir, YUBR_arrayDir, COGL_arrayDir, AMTH_arrayDir, HIWA_arrayDir,
                  QULO_arrayDir, PICO_arrayDir, MIGU_arrayDir, ARTH_arrayDir, VILA_arrayDir)

# %%% CALCULATE SUMMARY STATISTICS %%% ----
# Apply function calculating correlation values (NRMSE unrounded/rounded, Spearman, and Pearson) to list of SMBO2 datasets
NRMSEs <- lapply(SMBO2_values, buildCorrelationMat, corMetric = 'NRMSE', NRMSEdigits=5)
NRMSEs_Round <- lapply(SMBO2_values, buildCorrelationMat, corMetric = 'NRMSE')
corSps <- lapply(SMBO2_values, buildCorrelationMat, corMetric = 'corSp')
corPes <- lapply(SMBO2_values, buildCorrelationMat, corMetric = 'corPe')
names(NRMSEs) <- names(NRMSEs_Round) <- names(corSps) <- names(corPes) <- 
  c('QUAC','YUBR','COGL','AMTH','HIWA','QULO','PICO','MIGU','ARTH','VILA')

# Export resutling lists to CSV in Datasets directory
lapply(NRMSEs, function(x) write.table(data.frame(x), paste0(outputDir,'SMBO2_NRMSEs.csv'), append= T, sep=',' ))
lapply(NRMSEs_Round, function(x) write.table(data.frame(x), paste0(outputDir,'SMBO2_NRMSEs_Round.csv'), append= T, sep=',' ))
lapply(corSps, function(x) write.table(data.frame(x), paste0(outputDir,'SMBO2_corSps.csv'), append= T, sep=',' ))
lapply(corPes, function(x) write.table(data.frame(x), paste0(outputDir,'SMBO2_corPes.csv'), append= T, sep=',' ))

# Apply function calculating MSSEs from resampling arrays to list of SMBO2 datasets
MSSEs <- lapply(SMBO2_values, calcMSSEs)
names(MSSEs) <- c('QUAC','YUBR','COGL','AMTH','HIWA','QULO','PICO','MIGU','ARTH','VILA')
lapply(MSSEs, function(x) write.table(data.frame(x), paste0(outputDir,'SMBO2_MSSEs.csv'), append= T, sep=',' ))

# %%% PLOT SPEARMAN CORRELATION TABLE %%% ----
# Declare a matrix for capturing average correlation coefficient values over all buffer sizes for each species
corMat_Sps <- matrix(NA, ncol = 3, nrow = length(corSps))  
rownames(corMat_Sps) <- names(corSps) ; colnames(corMat_Sps) <- c('Geographic (Total buffer)', 'Geographic (SDM)', 'Ecological')
# Loop through the data.frames of Spearman rho values, calculating average correlation coefficients across buffer sizes for each species
for(i in 1:length(corSps)){
  # Need to conditionalize based on whether or not SDM geographic coverages are present
  if(ncol(corSps[[i]]) == 3){
    corMat_Sps[i,] <- apply(corSps[[i]], 2, averageCorrValue)
  } else{
    # If SDM values are not present, then include an NA value in the middle
    aveargeValues <- apply(corSps[[i]], 2, averageCorrValue)
    corMat_Sps[i,] <- c(aveargeValues[[1]], NA, aveargeValues[[2]])
  }
}
# Create a column which is an average across datasets, for each coverage type
corMat_Sps <- rbind(corMat_Sps, apply(corMat_Sps, 2, averageCorrValue))
rownames(corMat_Sps)[[11]] <- 'Coverage (Mean)'
# Create a row which is an average for each dataset
corMat_Sps <- cbind(corMat_Sps, apply(corMat_Sps, 1, averageCorrValue))
colnames(corMat_Sps)[[4]] <- 'Datasets (Mean)'

# PLOTTING IN GGPLOT
# Transpose the matrix, and convert it into long format
cor_matrix <- t(corMat_Sps)
cor_melted <- melt(cor_matrix)
# Preserve the original row order
cor_melted$Var1 <- factor(cor_melted$Var1, levels = rev(rownames(cor_matrix))) 
# Plot heatmap with text labels
ggplot(cor_melted, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white", width = 0.8, height = 0.8) +  # Reduce tile size
  geom_text(aes(label = round(value, 2)), color = "black", size = 5) +  # Increase text size
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Correlation") +
  scale_x_discrete(position = "top") +  # Move column names to the top
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 14),  # Larger x-axis text
        axis.text.y = element_text(size = 14),  # Larger y-axis text
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        aspect.ratio = nrow(cor_matrix) / ncol(cor_matrix),
        plot.title = element_text(vjust = -10, hjust=0.5)) +  # Control aspect ratio for smaller cells
  labs(title = "Mean Spearman Rho Values")
