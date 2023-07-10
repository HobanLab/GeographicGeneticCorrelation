# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FUNCTIONS FOR GEOGRAPHIC-GENETIC COVERAGE CALCULATIONS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script declares the functions used for the analysis of correlation between geographic
# and genetic coverage. The createBuffers and compareBuffArea functions were sourced from 
# the gap analysis code developed by Emily Beckman Bruns.

library(adegenet)
library(terra)
library(parallel)

# Create buffers around points, using specified projection
createBuffers <- function(df, radius=1000, pt_proj="+proj=longlat +datum=WGS84", 
                          buff_proj="+proj=eqearth +datum=WGS84", boundary){
  # Turn occurrence point data into a SpatVector
  spat_pts <- vect(df, geom=c("decimalLongitude", "decimalLatitude"), crs=pt_proj)
  # Reproject to specified projection
  proj_df <- project(spat_pts, buff_proj)
  # Place buffer around each point, then dissolve into one polygon
  buffers <- buffer(proj_df,width=radius)
  buffers <- aggregate(buffers,dissolve = TRUE)
  # Clip by boundary so they don't extend into the water
  boundary <- project(boundary,buff_proj)
  buffers_clip <- crop(buffers,boundary)
  # Return buffer polygons
  return(buffers_clip)
}

# Given coordinate points and vector of sample names, calculate geographic coverage
compareBuffArea <- function(insitu, samp_vect, radius, pt_proj, buff_proj, boundary, parFlag=FALSE){
  # If running in parallel: world polygon shapefile needs to be "unwrapped", after being exported to cluster
  if(parFlag==TRUE){
    boundary <- unwrap(boundary)
  }
  # Build ex situ points by subseting complete in situ points data.frame, according to samp_vect
  exsitu <- insitu[sort(match(samp_vect, insitu[,1])),]
  # Create buffers
  buffer_insitu <- createBuffers(insitu,radius,pt_proj,buff_proj,boundary)
  buffer_exsitu <- createBuffers(exsitu,radius,pt_proj,buff_proj,boundary)
  # Calculate buffer area
  area_exsitu <- expanse(buffer_exsitu)/1000000
  area_insitu <- expanse(buffer_insitu)/1000000
  # Calculate difference between in situ and ex situ buffer areas (% coverage)
  area_diff_percent <- (area_exsitu/area_insitu)*100
  return(area_diff_percent)
}

# Function for reporting representation rates, using a vector of allele frequencies and a sample matrix.
# Assumes that freqVector represents the absolute allele frequencies for the population of interest 
# (typically, entire wild population). Allele names between frequency vector and sample matrix must match! 
# 1. The length of matches between garden and wild alleles is calculated (numerator). 
# 2. The complete number of wild alleles of that category (denominator) is calculated. 
# 3. From these 2 values, a percentage is calculated. 
# This function returns the numerators, denominators, and the proportion (representation rates) in a matrix.
getAlleleCategories <- function(freqVector, sampleMat){
  # Determine how many Total alleles in the sample matrix are found in the frequency vector 
  exSitu_totalAlleles <- length(which(names(freqVector) %in% colnames(sampleMat)))
  inSitu_totalAlleles <- length(freqVector)
  totalPercentage <- (exSitu_totalAlleles/inSitu_totalAlleles)*100
  # Very common alleles (greater than 10%)
  exSitu_vComAlleles <- length(which(names(which(freqVector > 10)) %in% colnames(sampleMat)))
  inSitu_vComAlleles <- length(which(freqVector > 10))
  vComPercentage <- (exSitu_vComAlleles/inSitu_vComAlleles)*100
  # Common alleles (greater than 5%)
  exSitu_comAlleles <- length(which(names(which(freqVector > 5)) %in% colnames(sampleMat)))
  inSitu_comAlleles <- length(which(freqVector > 5))
  comPercentage <- (exSitu_comAlleles/inSitu_comAlleles)*100
  # Low frequency alleles (between 1% and 10%)
  exSitu_lowFrAlleles <- length(which(names(which(freqVector < 10 & freqVector > 1)) %in% colnames(sampleMat)))
  inSitu_lowFrAlleles <- length(which(freqVector < 10 & freqVector > 1))
  lowFrPercentage <- (exSitu_lowFrAlleles/inSitu_lowFrAlleles)*100
  # Rare alleles (less than 1%)
  exSitu_rareAlleles <- length(which(names(which(freqVector < 1)) %in% colnames(sampleMat)))
  inSitu_rareAlleles <- length(which(freqVector < 1))
  rarePercentage <- (exSitu_rareAlleles/inSitu_rareAlleles)*100
  # Concatenate values to vectors
  exSituAlleles <- c(exSitu_totalAlleles, exSitu_vComAlleles, exSitu_comAlleles, exSitu_lowFrAlleles, exSitu_rareAlleles)
  inSituAlleles <- c(inSitu_totalAlleles, inSitu_vComAlleles, inSitu_comAlleles, inSitu_lowFrAlleles, inSitu_rareAlleles)
  repRates <- c(totalPercentage,vComPercentage,comPercentage,lowFrPercentage,rarePercentage) 
  # Bind vectors to a matrix, name dimensions, and return
  exSituValues <- cbind(exSituAlleles, inSituAlleles, repRates)
  rownames(exSituValues) <- c("Total","V. common","Common", "Low freq.","Rare")
  colnames(exSituValues) <- c("Ex situ", "In situ", "Rate (%)")
  return(exSituValues)
}

# Wrapper of getAlleleCategories and compareBuffArea. Given a genetic matrix (rows are samples, columns are alleles)
# and a data.frame of coordinates (3 columns: sample names, latitudes, and longitudes), it calculates the genetic
# and geographic coverage from a random draw of some amount of samples (numSamples). 
calculateCoverage <- function(gen_mat, geo_coordPts, geo_buff, geo_ptProj, geo_buffProj, 
                              geo_boundary, parFlag=FALSE, numSamples){
  # Check that sample names in genetic matrix match the column of sample names in the coordinate data.frame
  if(!identical(rownames(gen_mat), geo_coordPts[,1])){
    stop("Error: Sample names between the genetic matrix and the first 
         column of the coordinate point data.frame do not match.")
  }
  # GENETIC PROCESSING
  # Calculate a vector of allele frequencies, based on the total sample matrix
  freqVector <- colSums(gen_mat, na.rm = TRUE)/(nrow(gen_mat)*2)*100
  # Remove any missing alleles (those with frequencies of 0) from the frequency vector
  freqVector <- freqVector[which(freqVector != 0)]
  # From a matrix of individuals, select a set of random individuals (rows)
  samp <- gen_mat[sample(nrow(gen_mat), size=numSamples, replace = FALSE),]
  # Remove any missing alleles (those with colSums of 0) from the sample matrix
  samp <- samp[,which(colSums(samp, na.rm = TRUE) != 0)]
  # Genetic coverage: calculate sample's allelic representation
  gen_rates <- getAlleleCategories(freqVector, samp)
  # Subset matrix returned by getAlleleCategories to just 3rd column (representation rates), and return
  gen_rates <- gen_rates[,3]
  # GEOGRAPHIC PROCESSING
  # Geographic coverage: calculate sample's geographic representation, by passing sample names to compareBuffArea
  geo_rate <- compareBuffArea(geo_coordPts, rownames(samp), geo_buff, geo_ptProj, geo_buffProj, geo_boundary, parFlag)
  names(geo_rate) <- "Geo"
  # Combine genetic and geographic coverage rates into a vector, and return
  cov_rates <- c(gen_rates, geo_rate)
  return(cov_rates)
}

# Wrapper of exSitu_Sample: iterates that function over the entire sample matrix
exSituResample <- function(gen_obj, geo_coordPts, geo_buff=1000, geo_ptProj="+proj=longlat +datum=WGS84", 
                           geo_buffProj="+proj=eqearth +datum=WGS84", geo_boundary, parFlag){
  # Check populations of samples: if NULL, provide all samples with this popname (i.e. assume all are wild)
  if(is.null(pop(gen_obj))){
    pop(gen_obj) <- rep("wild", nInd(gen_obj))
  }
  # Create a matrix of wild individuals (those with population "wild") from genind object
  gen_mat <- gen_obj@tab[which(pop(gen_obj) == "wild"),]
  # Apply the exSitu_Sample function to all rows of the wild matrix
  # (except row 1, because we need at least 2 individuals to sample)
  # The resulting matrix needs to be transposed, in order to keep columns as different allele categories
  cov_matrix <- t(sapply(2:nrow(gen_mat), 
                         function(x) calculateCoverage(gen_mat=gen_mat, geo_coordPts=geo_coordPts, geo_buff=geo_buff,
                                                       geo_ptProj=geo_ptProj, geo_buffProj=geo_buffProj, 
                                                       geo_boundary=geo_boundary, parFlag=parFlag, numSamples=x)))
  # Return the matrix of representation values
  return(cov_matrix)
}

# Wrapper of geo.gen.Resample: runs resampling in parallel over a specified cluster
# Results are saved to a specified file path.
geo.gen.Resample.Parallel <- function(gen_obj, geo_coordPts, geo_buff=1000, geo_ptProj="+proj=longlat +datum=WGS84",
                                      geo_buffProj="+proj=eqearth +datum=WGS84", geo_boundary, reps=5,
                                      arrayFilepath="~/resamplingArray.Rdata", cluster){
  # Run resampling in parallel, capturing results to an array
  resamplingArray <- parSapply(cluster, 1:reps, 
                               function(x) {exSituResample(gen_obj=gen_obj, geo_coordPts=geo_coordPts, geo_buff=geo_buff,
                                                          geo_ptProj=geo_ptProj, geo_buffProj=geo_buffProj, 
                                                          geo_boundary=geo_boundary, parFlag=TRUE)}, 
                               simplify = "array")
  # Save the resampling array object to disk, for later usage
  saveRDS(resamplingArray, file = arrayFilepath)
  # Return resampling array to global environment
  return(resamplingArray)
}

# Wrapper for exSituResample, which will generate an array of values from a single genind object
geo.gen.Resample <- function(gen_obj, geo_coordPts, geo_buff=1000, 
                             geo_ptProj="+proj=longlat +datum=WGS84", 
                             geo_buffProj="+proj=eqearth +datum=WGS84", geo_boundary, reps=5){
  # Run resampling for all replicates, using sapply and lambda function
  resamplingArray <- sapply(1:reps, 
                            function(x) exSituResample(gen_obj=gen_obj, geo_coordPts=geo_coordPts, geo_buff=geo_buff,
                                                       geo_ptProj=geo_ptProj, geo_buffProj=geo_buffProj, 
                                                       geo_boundary=geo_boundary, parFlag=FALSE), 
           simplify = "array")
  # Return array
  return(resamplingArray)
}

# From resampling array, calculate the mean minimum sample size to represent 95% of the Total wild diversity
gen_min95Mean <- function(resamplingArray){
  # resampling array[,1,]: returns the Total column values for each replicate (3rd array dimension)
  # apply(resamplingArray[,1,],1,mean): calculates the average across replicates for each row
  # which(apply(resamplingArray[,1,],1,mean) > 95): returns the rows with averages greater than 95
  # min(which(apply(resamplingArray[,1,],1,mean) > 95)): the lowest row with an average greater than 95
  meanValue <- min(which(apply(resamplingArray[,1,],1, mean, na.rm=TRUE) > 95))
  return(meanValue)
}

# From resampling array, calculate the standard deviation, at the mean 95% value
gen_min95SD <- function(resamplingArray){
  # Determine the mean value for representing 95% of allelic diversity
  meanValue <- gen_min95Mean(resamplingArray)
  # Calculate the standard deviation, at that mean value, and return
  sdValue <- apply(resamplingArray[,1,],1,sd)[meanValue]
  return(sdValue)
}

# From resampling array, calculate the mean minimum sample size to represent 95% of the Total wild diversity
geo_min95Mean <- function(resamplingArray){
  meanValue <- min(which(apply(resamplingArray[,6,],1, mean, na.rm=TRUE) > 95))
  return(meanValue)
}

# From resampling array, calculate the standard deviation, at the mean 95% value
geo_min95SD <- function(resamplingArray){
  # Determine the mean value for representing 95% of allelic diversity
  meanValue <- geo_min95Mean(resamplingArray)
  # Calculate the standard deviation, at that mean value, and return
  sdValue <- apply(resamplingArray[,1,],1,sd)[meanValue]
  return(sdValue)
}

# From resampling array, calculate the mean values (across replicates) for each allele frequency category
meanArrayValues <- function(resamplingArray, allValues=FALSE){
  # Declare a matrix to receive average values
  meanValue_mat <- matrix(nrow=nrow(resamplingArray), ncol=ncol(resamplingArray))
  # For each column in the array, average results across replicates (3rd array dimension)
  for(i in 1:ncol(resamplingArray)){
    meanValue_mat[,i] <- apply(resamplingArray[,i,], 1, mean, na.rm=TRUE)
  }
  # Give names to all meanValue_mat columns
  colnames(meanValue_mat) <- c("Total","V. common","Common","Low freq.","Rare", "Geo")
  if(allValues==FALSE){
    # Unless returning all matrix values is specified, return just the first ("Total") and last ("Geo") columns
    meanValue_mat <- meanValue_mat[,-(2:5)]
  } 
  # Reformat the matrix as a data.frame, and return
  meanValue_mat <- as.data.frame(meanValue_mat)
  return(meanValue_mat)
}

# DATA EXPLORATION FUNCTIONS ----
# Function for generating a vector of wild allele frequencies from a genind object
getWildFreqs <- function(gen.obj){
  # Build a vector of rows corresponding to wild individuals (those that do not have a population of "garden")
  wildRows <- which(pop(gen.obj)!="garden")
  # Build the wild allele frequency vector: colSums of alleles (removing NAs), divided by number of haplotypes (Ne*2)
  wildFreqs <- colSums(gen.obj@tab[wildRows,], na.rm = TRUE)/(length(wildRows)*2)*100
  return(wildFreqs)
}

# Function for generating a vector of total allele frequencies from a genind object
getTotalFreqs <- function(gen.obj){
  # Build allele frequency vector: colSums of alleles (removing NAs), divided by number of haplotypes (Ne*2)
  totalFreqs <- colSums(gen.obj@tab, na.rm = TRUE)/(nInd(gen.obj)*2)*100
  return(totalFreqs)
}

# Exploratory function for reporting the proprtion of alleles of each category, from a (wild) frequency vector
getWildAlleleFreqProportions <- function(gen.obj){
  # Build the wild allele frequency vector, using the getWildFreqs function
  wildFreqs <- getWildFreqs(gen.obj)
  # Very common
  veryCommonAlleles <- wildFreqs[which(wildFreqs > 10)]
  veryCommon_prop <- (length(veryCommonAlleles)/length(wildFreqs))*100
  # Low frequency
  lowFrequencyAlleles <- wildFreqs[which(wildFreqs < 10 & wildFreqs > 1)]
  lowFrequency_prop <- (length(lowFrequencyAlleles)/length(wildFreqs))*100
  # Rare
  rareAlleles <- wildFreqs[which(wildFreqs < 1)]
  rare_prop <- (length(rareAlleles)/length(wildFreqs))*100
  # Build list of proportions, and return
  freqProportions <- c(veryCommon_prop, lowFrequency_prop, rare_prop)
  names(freqProportions) <- c("Very common (>10%)","Low frequency (1% -- 10%)","Rare (<1%)")
  return(freqProportions)
}

# Exploratory function for reporting the proprtion of alleles of each category, 
# from a frequency vector (of ALL alleles--garden AND wild)
getTotalAlleleFreqProportions <- function(gen.obj){
  # Build the wild allele frequency vector, using the getWildFreqs function
  totalFreqs <- getTotalFreqs(gen.obj)
  # Very common
  veryCommonAlleles <- totalFreqs[which(totalFreqs > 10)]
  veryCommon_prop <- (length(veryCommonAlleles)/length(totalFreqs))*100
  # Low frequency
  lowFrequencyAlleles <- totalFreqs[which(totalFreqs < 10 & totalFreqs > 1)]
  lowFrequency_prop <- (length(lowFrequencyAlleles)/length(totalFreqs))*100
  # Rare
  rareAlleles <- totalFreqs[which(totalFreqs < 1)]
  rare_prop <- (length(rareAlleles)/length(totalFreqs))*100
  # Build list of proportions, and return
  freqProportions <- c(veryCommon_prop, lowFrequency_prop, rare_prop)
  names(freqProportions) <- c("Very common (>10%)","Low frequency (1% -- 10%)","Rare (<1%)")
  return(freqProportions)
}
