# Script looking through SMBO2 results to find buffer sizes which lead to the closest
# match between the 95% MSSEs

pacman::p_load(adegenet, terra, parallel, RColorBrewer, scales, vcfR, usedist, DescTools)

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

# Function for building a matrix of correlation values (NRMSE, Pearson correlation, or Spearman correlation)
# based on an array path (resampling array to read in) and the genetic coverage type used for the predictor variable.
# Given these arguments, a matrix is generated which has a NRMSE value for each set of geographic/ecological
# coverage values compared to the genetic coverage metric. This command expects the resampling data.frame
# to have certain row names and column names, and makes those checks
buildCorrelationMat <- function(arrayDir='~/resamp.Rdata', genCovType=c('CV', 'GD'), 
                                corMetric=c('NRMSE', 'corSp', 'corPe'), sdmFlag=FALSE){
  # Match arguments for type of correlation metric to calculate and response variable (genetic coverage approach) to use
  corMetric <- match.arg(corMetric) ; genCovType <- match.arg(genCovType)
  # Read in the specified resampling array
  resampArr <- readRDS(arrayDir)
  # Convert the resampling array to a data.frame
  resampDF <- resample.array2dataframe(resampArr)
  # If CV is chosen, set the predictive variable argument to the 'Total' column in the data.frame
  if(genCovType=='CV'){
    predVar <- resampDF$Total
    # Otherwise, set the predictve variable to the 'GenDist' column
  } else {
    predVar <- resampDF$GenDist
  }
  # Determine whether or not geographic coverages using an SDM are present, and set a flag accordingly
  if(length(grep('SDM', colnames(resampDF))) > 1){
    sdmFlag <- TRUE
  }
  # Extract buffer sizes through column names, for the geographic (total buffer) approach. This is 
  # a complicated command, but essentially relies on the column names of the resampling data.frame
  # to refer to the geographic buffer size used for that row of data.
  buffSizes <- 
    1000*(as.numeric(sub('km','', sapply(strsplit(grep('Geo_Buff_', colnames(resampDF), value=TRUE), '_'),'[',3))))
  # Check that the resampling data.frame has the same number of columns for Geo and Eco coverage
  # (We can't implement this check for SDM, because some datasets have fewer SDM buffer sizes)
  if(length(grep('Geo_Buff_', colnames(resampDF)))!=length(grep('Eco_Buff_', colnames(resampDF)))){
    stop('Different number of Geo_Buff and Eco_Buff columns in data.frame!')
  }
  # Build matrix of NRMSE values. Depending on the sdmFlag, this matrix will have 3 columns or just 2. 
  if(sdmFlag==TRUE){
    # Name matrix columns accordingly
    corMat <- matrix(NA, nrow=length(buffSizes), ncol=3)
    colnames(corMat) <- c('Geo_Buff','Geo_SDM','Eco_Buff')
  } else {
    # Name matrix columns accordingly
    corMat <- matrix(NA, nrow=length(buffSizes), ncol=2)
    colnames(corMat) <- c('Geo_Buff','Eco_Buff')
  }
  # Name matrix rows according to buffer sizes
  rownames(corMat) <- paste0(buffSizes/1000, 'km')
  # Not all resampling data.frames have the geographic coverage results starting in the same column.
  # To address this, make the loop start at an index according to the column names
  startCol <- min(grep('Geo_Buff_', colnames(resampDF)))
  # Loop through the dataframe columns. The first three columns are skipped, as they're sample number and the
  # predictive variables (allelic coverage and genetic distance proportions)
  for(i in startCol:ncol(resampDF)){
    # Syntax below allows the appropriate correlation metric to be calculated, according to function argument
    if (corMetric=='NRMSE') {
      # Calculate NRMSE value
      corValue <- nrmse.func(resampDF[,i], pred = predVar)
    } else if (corMetric=='corSp') {
      # Calculate Spearman correlation value
      corValue <- cor.test(resampDF[,i], predVar, method='spearman')$estimate
    } else {
      # Calculate Pearson correlation value
      corValue <- cor.test(resampDF[,i], predVar, method='pearson')$estimate
    }
    # Get the name of the current dataframe column
    dataName <- unlist(strsplit(names(resampDF)[[i]],'_'))
    # Match the data name to the relevant rows/columns of the receiving matrix
    matRow <- which(rownames(corMat) == dataName[[3]])
    matCol <- which(colnames(corMat) == paste0(dataName[[1]],'_',dataName[[2]]))
    # Locate the NRMSE value accordingly
    corMat[matRow,matCol] <- corValue
  }
  # Rename the columns based on the genetic coverage type used
  if(genCovType=='CV'){
    colnames(corMat) <- paste0(colnames(corMat),'_CV')
  } else {
    colnames(corMat) <- paste0(colnames(corMat),'_GD')
  }
  # Return the resulting matrix
  return(corMat)
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
# Apply function calculating correlation values (NRMSE, Spearman, and Pearson) from resampling arrays to list of SMBO2 datasets
NRMSEs <- lapply(SMBO2_values, buildCorrelationMat, corMetric = 'NRMSE')
corSps <- lapply(SMBO2_values, buildCorrelationMat, corMetric = 'corSp')
corPes <- lapply(SMBO2_values, buildCorrelationMat, corMetric = 'corPe')
names(NRMSEs) <- names(corSps) <- names(corPes) <- c('QUAC','YUBR','COGL','AMTH','HIWA','QULO','PICO','MIGU','ARTH','VILA')

# Export resutling lists to CSV in Datasets directory
lapply(NRMSEs, function(x) write.table(data.frame(x), paste0(outputDir,'SMBO2_NRMSEs.csv'), append= T, sep=',' ))
lapply(corSps, function(x) write.table(data.frame(x), paste0(outputDir,'SMBO2_corSps.csv'), append= T, sep=',' ))
lapply(corPes, function(x) write.table(data.frame(x), paste0(outputDir,'SMBO2_corPes.csv'), append= T, sep=',' ))

# Apply function calculating MSSEs from resampling arrays to list of SMBO2 datasets
MSSEs <- lapply(SMBO2_values, calcMSSEs)
names(MSSEs) <- c('QUAC','YUBR','COGL','AMTH','HIWA','QULO','PICO','MIGU','ARTH','VILA')
lapply(MSSEs, function(x) write.table(data.frame(x), paste0(outputDir,'SMBO2_MSSEs.csv'), append= T, sep=',' ))

# %%% PLOT SPEARMAN CORRELATION TABLE %%% ----
# Declare a matrix for capturing average correlation coefficient values over all buffer sizes for each species
corMat_Sps <- matrix(NA, ncol = 3, nrow = length(corSps))  
rownames(corMat_Sps) <- names(corSps) ; colnames(corMat_Sps) <- colnames(corSps$QUAC)
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

# NEED TO FIND A WAY TO PLOT THIS DATA, EITHER USING CORRPLOT OR SOMETHING ELSE...
