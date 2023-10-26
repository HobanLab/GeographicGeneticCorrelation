# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FUNCTIONS FOR GEOGRAPHIC-GENETIC COVERAGE CALCULATIONS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script declares the functions used for the analysis of correlation between geographic
# and genetic coverage. The createBuffers and compareBuffArea functions were sourced from 
# the gap analysis code developed by Emily Beckman Bruns.

library(adegenet)
library(terra)
library(parallel)

# ---- BUILDING THE RESAMPLING ARRAY ----
# Create buffers around points, using specified projection
createBuffers <- function(df, radius=1000, ptProj='+proj=longlat +datum=WGS84', 
                          buffProj='+proj=eqearth +datum=WGS84', boundary){
  # Turn occurrence point data into a SpatVector
  spat_pts <- vect(df, geom=c('decimalLongitude', 'decimalLatitude'), crs=ptProj)
  # Reproject to specified projection
  proj_df <- project(spat_pts, buffProj)
  # Place buffer around each point, then dissolve into one polygon
  buffers <- buffer(proj_df,width=radius)
  buffers <- aggregate(buffers,dissolve = TRUE)
  # Clip by boundary so they don't extend into the water
  boundary <- project(boundary,buffProj)
  buffers_clip <- crop(buffers,boundary)
  # Return buffer polygons
  return(buffers_clip)
}

# WORKER FUNCTION: Given coordinate points and vector of sample names, calculate geographic coverage
geo.compareBuff <- function(inSitu, sampVect, radius, ptProj, buffProj, boundary, parFlag=FALSE){
  # If running in parallel: world polygon shapefile needs to be 'unwrapped', 
  # after being exported to cluster
  if(parFlag==TRUE){
    boundary <- unwrap(boundary)
  }
  # Build ex situ points by subseting complete in situ points data.frame, according to sampVect
  exSitu <- inSitu[sort(match(sampVect, inSitu[,1])),]
  # Create buffers
  geo_exSitu <- createBuffers(exSitu, radius, ptProj, buffProj, boundary)
  geo_inSitu <- createBuffers(inSitu, radius, ptProj, buffProj, boundary)
  # Calculate buffer area. 1,000,000 value is to convert area to km²
  geo_exSituArea <- expanse(geo_exSitu)/1000000
  geo_inSituArea <- expanse(geo_inSitu)/1000000
  # Calculate difference between in situ and ex situ buffer areas (% coverage)
  geo_Coverage <- (geo_exSituArea/geo_inSituArea)*100
  return(geo_Coverage)
}

# Create a data.frame with ecoregion data extracted for area covered by buffers
eco.intersectBuff <- function(df, radius, ptProj, buffProj, ecoRegion, boundary, parFlag=FALSE){
  # If running in parallel: world polygon and ecoregions shapefiles need to be 'unwrapped', 
  # after being exported to cluster
  if(parFlag==TRUE){
    ecoRegion <- unwrap(ecoRegion)
    boundary <- unwrap(boundary)
  }
  # Create buffers
  buffers <- createBuffers(df, radius, ptProj, buffProj, boundary)
  # Make sure ecoregions are in same projection as buffers
  ecoProj <- project(ecoRegion, buffProj)
  # Intersect buffers with ecoregions, and return
  ecoBuffJoin <- intersect(buffers, ecoProj)
  return(ecoBuffJoin)
}

# WORKER FUNCTION: Create a data.frame with ecoregion data extracted for area covered by buffers
# surrounding all points and sample points. Then, compare the ecoregions count of the
# sample to the total. The layerType argument allows for 3 possible values: US (EPA Level 4),
# NA (EPA Level 3), and GL (TNC Global Terrestrial) ecoregions. These should correspond with 
# the ecoRegion argument (which specifies the ecoregion shapefile).
eco.compareBuff <- function(inSitu, sampVect, radius, ptProj, buffProj, 
                            ecoRegion, layerType=c('US','NA','GL'), boundary, parFlag=FALSE){
  # Match layerType argument, which specifies which ecoregion data type to extract (below)
  layerType <- match.arg(layerType)
  # If running in parallel: world polygon and ecoregions shapefiles need to be 'unwrapped', 
  # after being exported to cluster
  if(parFlag==TRUE){
    ecoRegion <- unwrap(ecoRegion)
    boundary <- unwrap(boundary)
  }
  # Build sample ex situ points by subseting complete in situ points data.frame, according to sampVect
  exSitu <- inSitu[sort(match(sampVect, inSitu[,1])),]
  # Create data.frame of ecoregion-buffer intersection
  eco_exSitu <- eco.intersectBuff(exSitu, radius, ptProj, buffProj, ecoRegion, boundary)
  eco_inSitu <- eco.intersectBuff(inSitu, radius, ptProj, buffProj, ecoRegion, boundary)
  # Based on the ecoRegion shapefile and the specified layer type, count the number of ecoregions in
  # the random sample (exSitu) and all of the data points (inSitu)
  if(layerType=='US'){
    # Extract the number of EPA Level IV ('U.S. Only') ecoregions
    eco_exSituCount <- length(unique(eco_exSitu$US_L4CODE))
    eco_inSituCount <- length(unique(eco_inSitu$US_L4CODE))
  } else {
    # Extract the number of EPA Level III ('North America') ecoregions 
    if(layerType=='NA'){
      eco_exSituCount <- length(unique(eco_exSitu$NA_L3CODE))
      eco_inSituCount <- length(unique(eco_inSitu$NA_L3CODE))
    } else {
      # Extract the number of Nature Conservancy ('Global Terrestrial') ecoregions
      eco_exSituCount <- length(unique(eco_exSitu$ECO_ID_U))
      eco_inSituCount <- length(unique(eco_inSitu$ECO_ID_U))
    }
  }
  # Calculate difference in number of ecoregions between the sample and all data points, and return
  eco_Coverage <- (eco_exSituCount/eco_inSituCount)*100
  return(eco_Coverage)
}

# WORKER FUNCTION: Function for reporting representation rates, using a vector of allele frequencies 
# and a sample matrix. Assumes that freqVector represents the absolute allele frequencies 
# for the population of interest (typically, entire wild population). Allele names between 
# frequency vector and sample matrix must match! 
# 1. The length of matches between garden and wild alleles is calculated (numerator). 
# 2. The complete number of wild alleles of that category (denominator) is calculated. 
# 3. From these 2 values, a percentage is calculated. 
# This function returns the numerators, denominators, and the proportion (representation rates) in a matrix.
gen.getAlleleCategories <- function(freqVector, sampleMat){
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
  rownames(exSituValues) <- c('Total','V. common','Common', 'Low freq.','Rare')
  colnames(exSituValues) <- c('Ex situ', 'In situ', 'Rate (%)')
  return(exSituValues)
}

# Wrapper of gen.getAlleleCategories, geo.compareBuff, and eco.compareBuff worker functions. 
# Given a genetic matrix (rows are samples, columns are alleles) and a data.frame of coordinates 
# (3 columns: sample names, latitudes, and longitudes), it calculates the genetic,
# geographic (if flagged), and ecologcial (if flagged) coverage from a random draw of some amount of 
# samples (numSamples). The sample names between the genind object and the coordinate data.frame need
# to match (in order to properly subset across genetic, geographic, and ecological datasets). This is
# the core function of the gen-geo-eco resampling workflow.
calculateCoverage <- function(gen_mat, geoFlag=TRUE, coordPts, geoBuff=50000, 
                              ptProj='+proj=longlat +datum=WGS84',
                              buffProj='+proj=eqearth +datum=WGS84', boundary,
                              ecoFlag=TRUE, ecoBuff=50000, ecoRegions, ecoLayer=c('US','NA','GL'),
                              parFlag=FALSE, numSamples){
  
  # Check that sample names in genetic matrix match the column of sample names in the coordinate data.frame
  if(!identical(rownames(gen_mat), coordPts[,1])){
    stop('Error: Sample names between the genetic matrix and the first 
         column of the coordinate point data.frame do not match.')
  }
  
  # GENETIC PROCESSING
  # Calculate a vector of allele frequencies, based on the total sample matrix
  freqVector <- colSums(gen_mat, na.rm = TRUE)/(nrow(gen_mat)*2)*100
  # Remove any missing alleles (those with frequencies of 0) from the frequency vector
  freqVector <- freqVector[which(freqVector != 0)]
  # From a matrix of individuals, select a set of random individuals (rows). This is the set of individuals that will
  # be used for all downstream coverage calculations within this function (genetic, geographic, and ecological)
  samp <- gen_mat[sample(nrow(gen_mat), size=numSamples, replace = FALSE),]
  # Remove any missing alleles (those with colSums of 0) from the sample matrix
  samp <- samp[,which(colSums(samp, na.rm = TRUE) != 0)]
  # Genetic coverage: calculate sample's allelic representation
  genRates <- gen.getAlleleCategories(freqVector, samp)
  # Subset matrix returned by getAlleleCategories to just 3rd column (representation rates), and return
  genRates <- genRates[,3]
  
  # GEOGRAPHIC PROCESSING
  if(geoFlag==TRUE){
    # Check for the required arguments (ptProj and buffProj will use defaults, if not specified)
    if(missing(coordPts)) stop('For geographic coverage, a data.frame of wild coordinates (coordPts) is required')
    if(missing(geoBuff)) stop('For geographic coverage, an integer specifying the geographic buffer size (geoBuff) is required')
    if(missing(boundary)) stop('For geographic coverage, a SpatVector object of country boundaries (boundary) is required')
    # Check that the names of the latitude and longitude columns are properly written (this is unfortunately hard-coded)
    if(!identical(colnames(coordPts)[2:3], c('decimalLatitude', 'decimalLongitude'))){
      stop('The column names of the geographic coordinates data.frame (coordPts) need to be 
           decimalLatitude and decimalLongitude. Please rename your data.frame of geographic coordinates!')
    }
    # Geographic coverage: calculate sample's geographic representation, by passing all points (coordPts) and 
    # the random subset of points (rownames(samp)) to the geo.compareBuff worker function, which will calculate
    # the proportion of area covered in the random sample
    geoRate <- geo.compareBuff(inSitu=coordPts, sampVect=rownames(samp), radius=geoBuff, 
                               ptProj=ptProj, buffProj=buffProj, boundary=boundary, parFlag=parFlag)
  } else {
    geoRate <- NA
  }
  names(geoRate) <- 'Geo'
  
  # ECOLOGICAL PROCESSING
  if(ecoFlag==TRUE){
    # Match ecoLayer argument, to ensure it is 1 of 3 possible values ('US', 'NA', 'GL)
    ecoLayer <- match.arg(ecoLayer)
    # Check for the required arguments (ptProj, buffProj, and ecoLayer will use defaults, if not specified)
    if(missing(coordPts)) stop('For ecological coverage, a data.frame of wild coordinates (coordPts) is required')
    if(missing(ecoBuff)) stop('For ecological coverage, an integer specifying the ecological buffer size (ecoBuff) is required')
    if(missing(ecoRegions)) stop('For ecological coverage, a SpatVector object of ecoregions (ecoregions) is required')
    if(missing(boundary)) stop('For ecological coverage, a SpatVector object of country boundaries (boundary) is required')
    # Ecological coverage: calculate sample's ecological representation, by passing all points (coordPts) and 
    # the random subset of points (rownames(samp)) to the eco.compareBuff worker function, which will calculate
    # the proportion of ecoregions covered in the random sample
    ecoRate <- eco.compareBuff(inSitu=coordPts, sampVect=rownames(samp), radius=ecoBuff, 
                               ptProj=ptProj, buffProj=buffProj, ecoRegion=ecoRegions, 
                               layerType=ecoLayer, boundary=boundary, parFlag=parFlag)
  } else{
    ecoRate <- NA
  }
  names(ecoRate) <- 'Eco'
  
  # Combine genetic, geographic, and ecological coverage rates into a vector, and return
  covRates <- c(genRates, geoRate, ecoRate)
  return(covRates)
}

# !!! WAYS FOR CALCULATECOVERAGE TO IMPROVE !!!
# 1. Move calculation of ecoregions for all sample points out of the innermost function
# 2. Move calculation of all allele frequencies out of the innermost function

# Wrapper of exSitu_Sample: iterates calculateCoverage over the entire matrix of samples
exSituResample <- function(gen_obj, geoFlag=TRUE, coordPts, geoBuff=1000, ptProj='+proj=longlat +datum=WGS84', 
                           buffProj='+proj=eqearth +datum=WGS84', boundary, ecoFlag=TRUE, ecoBuff=1000,
                           ecoRegions, ecoLayer='US', parFlag){
  # Check populations of samples: if NULL, provide all samples with the popname 'wild' 
  # This means that if no populations are specified in the genind object, all samples will be used!
  if(is.null(pop(gen_obj))){
    pop(gen_obj) <- rep('wild', nInd(gen_obj))
  }
  # Create a matrix of wild individuals (those with population 'wild') from genind object
  gen_mat <- gen_obj@tab[which(pop(gen_obj) == 'wild'),]
  # Apply the calculateCoverage function to all rows of the wild matrix
  # (except row 1, because we need at least 2 individuals to sample)
  # The resulting matrix needs to be transposed, in order to keep columns as different coverage categories
  cov_matrix <- t(sapply(2:nrow(gen_mat), 
                         function(x) calculateCoverage(gen_mat=gen_mat, geoFlag=geoFlag, coordPts=coordPts, 
                                                       geoBuff=geoBuff, ptProj=ptProj, buffProj=buffProj, 
                                                       boundary=boundary, ecoFlag=ecoFlag, ecoBuff=ecoBuff,
                                                       ecoRegions=ecoRegions, ecoLayer=ecoLayer,
                                                       parFlag=parFlag, numSamples=x)))
  # Return the matrix of coverage values
  return(cov_matrix)
}

# Wrapper of exSituResample: runs resampling in parallel over a specified cluster
# Results are saved to a specified file path.
geo.gen.Resample.Parallel <- function(gen_obj, geoFlag=TRUE, coordPts, geoBuff=1000, ptProj='+proj=longlat +datum=WGS84',
                                      buffProj='+proj=eqearth +datum=WGS84', boundary, ecoFlag=TRUE, ecoBuff=1000, 
                                      ecoRegions, ecoLayer=c('US','NA','GL'), reps=5,
                                      arrayFilepath='~/resamplingArray.Rdata', cluster){
  # Print starting time
  startTime <- Sys.time() 
  print(paste0('%%% RUNTIME START: ', startTime))
  # Run resampling in parallel, capturing results to an array
  resamplingArray <- parSapply(cluster, 1:reps, 
                               function(x) {exSituResample(gen_obj=gen_obj, geoFlag=geoFlag, coordPts=coordPts, 
                                                           geoBuff=geoBuff, ptProj= ptProj, buffProj=buffProj,
                                                           boundary=boundary, ecoFlag=ecoFlag, ecoBuff=ecoBuff,
                                                           ecoRegions=ecoRegions, ecoLayer=ecoLayer, parFlag=TRUE)}, 
                               simplify = 'array')
  # Print ending time and total runtime
  endTime <- Sys.time() 
  print(paste0('%%% RUNTIME END: ', endTime))
  cat(paste0('\n', '%%% TOTAL RUNTIME: ', endTime-startTime))
  # Save the resampling array object to disk, for later usage
  saveRDS(resamplingArray, file = arrayFilepath)
  # Return resampling array to global environment
  return(resamplingArray)
}

# Wrapper for exSituResample, which will generate an array of values from a single genind object
# This function doesn't run in parallel, so it's primarily used for testing/demonstration purposes
geo.gen.Resample <- function(gen_obj, geoFlag=TRUE, coordPts, geoBuff=1000, ptProj='+proj=longlat +datum=WGS84', 
                             buffProj='+proj=eqearth +datum=WGS84', boundary, ecoFlag=TRUE, ecoBuff=1000, 
                             ecoRegions, ecoLayer=c('US','NA','GL'), reps=5){
  # Run resampling for all replicates, using sapply and lambda function
  resamplingArray <- sapply(1:reps, 
                            function(x) exSituResample(gen_obj=gen_obj, geoFlag=geoFlag, coordPts=coordPts, geoBuff=geoBuff,
                                                       ptProj=ptProj, buffProj=buffProj, boundary=boundary, ecoFlag=ecoFlag,
                                                       ecoBuff=ecoBuff, ecoRegions=ecoRegions, ecoLayer=ecoLayer, parFlag=FALSE), 
                            simplify = 'array')
  # Return array
  return(resamplingArray)
}

# ---- PROCESSING THE RESAMPLING ARRAY ----
# From resampling array, calculate the mean minimum sample size to represent 95% of the Total wild diversity
gen.min95Mean <- function(resamplingArray){
  # resampling array[,1,]: returns the Total column values for each replicate (3rd array dimension)
  # apply(resamplingArray[,1,],1,mean): calculates the average across replicates for each row
  # which(apply(resamplingArray[,1,],1,mean) > 95): returns the rows with averages greater than 95
  # min(which(apply(resamplingArray[,1,],1,mean) > 95)): the lowest row with an average greater than 95
  meanValue <- min(which(apply(resamplingArray[,1,],1, mean, na.rm=TRUE) > 95))
  return(meanValue)
}

# From resampling array, calculate the standard deviation, at the mean 95% value
gen.min95SD <- function(resamplingArray){
  # Determine the mean value for representing 95% of allelic diversity
  meanValue <- gen_min95Mean(resamplingArray)
  # Calculate the standard deviation, at that mean value, and return
  sdValue <- apply(resamplingArray[,1,],1,sd)[meanValue]
  return(sdValue)
}

# From resampling array, calculate the mean minimum sample size to represent 95% of the Total wild diversity
geo.min95Mean <- function(resamplingArray){
  meanValue <- min(which(apply(resamplingArray[,6,],1, mean, na.rm=TRUE) > 95))
  return(meanValue)
}

# From resampling array, calculate the standard deviation, at the mean 95% value
geo.min95SD <- function(resamplingArray){
  # Determine the mean value for representing 95% of allelic diversity
  meanValue <- geo_min95Mean(resamplingArray)
  # Calculate the standard deviation, at that mean value, and return
  sdValue <- apply(resamplingArray[,6,],1,sd)[meanValue]
  return(sdValue)
}

# From resampling array, calculate the mean minimum sample size to represent 95% of the Total wild diversity
eco.min95Mean <- function(resamplingArray){
  meanValue <- min(which(apply(resamplingArray[,7,],1, mean, na.rm=TRUE) > 95))
  return(meanValue)
}

# From resampling array, calculate the standard deviation, at the mean 95% value
eco.min95SD <- function(resamplingArray){
  # Determine the mean value for representing 95% of allelic diversity
  meanValue <- eco_min95Mean(resamplingArray)
  # Calculate the standard deviation, at that mean value, and return
  sdValue <- apply(resamplingArray[,7,],1,sd)[meanValue]
  return(sdValue)
}

# From resampling array, calculate the mean values (across replicates) for each allele frequency category
# The allValues flag indicates whether or not to return the coverage metrics for alleles of different 
# categories (by default, the function will only return Total allelic coverage)
meanArrayValues <- function(resamplingArray, allValues=FALSE){
  # Declare a matrix to receive average values
  meanValues_mat <- matrix(nrow=nrow(resamplingArray), ncol=ncol(resamplingArray))
  # Name columns in the mean value matrix according to columns from input array
  colnames(meanValues_mat) <- colnames(resamplingArray)
  # For each column in the array, average results across replicates (3rd array dimension)
  for(i in 1:ncol(resamplingArray)){
    meanValues_mat[,i] <- apply(resamplingArray[,i,], 1, mean, na.rm=TRUE)
  }
  # Unless returning all matrix values is specified, return just the first ('Total') and last ('Geo') columns
  if(allValues==FALSE){
    meanValues_mat <- meanValues_mat[,-(2:5)]
  }
  # Reformat the matrix as a data.frame, and return
  meanValues <- as.data.frame(meanValues_mat)
  return(meanValues)
}

# From resampling array, generate a data.frame by collapsing values across replicates into vectors
# allValues flag indicates whether or not to include categories of alleles other that 'Total'
resample.array2dataframe <- function(resamplingArray, allValues=FALSE){
  # Create a vector of sample numbers. The values in this vector range from 2:total number
  # of samples (at least 2 samples are required in order for sample function to work; see above).
  # These values are repeated for the number of replicates in the resampling array (3rd dimension)
  sampleNumbers <- rep(2:(nrow(resamplingArray)+1), dim(resamplingArray)[[3]])
  # Pass sample number vector to data.frame, which will be the final output of the function
  resamp_DF <- data.frame(sampleNumbers=sampleNumbers)
  # Loop through the array by colunms (variables)
  for(i in 1:ncol(resamplingArray)){
    # For each, collapse the column into one long vector, and add that vector to the data.frame
    resamp_DF <- cbind(resamp_DF, c(resamplingArray[,i,]))
  }
  # Rename the data.frame values according to the column names of the array
  names(resamp_DF) <- c('sampleNumbers', colnames(resamplingArray))
  # If allValues flag is FALSE, remove the allele categories other than 'Total'
  if(allValues==FALSE){
    resamp_DF <- resamp_DF[,-(3:6)]
  }
  return(resamp_DF)
}

# ---- DATA EXPLORATION FUNCTIONS ----
# Function for generating a vector of wild allele frequencies from a genind object
getWildFreqs <- function(gen.obj){
  # Build a vector of rows corresponding to wild individuals (those that do not have a population of 'garden')
  wildRows <- which(pop(gen.obj)!='garden')
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
  names(freqProportions) <- c('Very common (>10%)','Low frequency (1% -- 10%)','Rare (<1%)')
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
  names(freqProportions) <- c('Very common (>10%)','Low frequency (1% -- 10%)','Rare (<1%)')
  return(freqProportions)
}
