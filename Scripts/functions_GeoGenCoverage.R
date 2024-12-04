# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FUNCTIONS FOR GENETIC-GEOGRAPHIC-ECOLOGICAL COVERAGE CALCULATIONS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script declares the functions used for correlation analyses between genetic, geographic, and
# ecological coverage metrics. The functions here are used to generate resampling arrays, in which the
# rows are numbers, the columns are the various coverage values (genetic, geographic, ecological), and the
# slices are different resampling replicates. In general, these resampling replicates are built using a 
# genetic matrix (part of a genind file) and a data.frame of decimal latitude and longitude values.

# Many of the functions here are wrappers, which sapply other lower-level functions. The advantage of this 
# nested function approach is that it allows for a single function to be called at the "upper-most" level of the code
# (i.e. the level at which data is read in). For instance, the geo.gen.Resample function (and its parallelized version,
# geo.gen.Resample.Par) is the only resampling function called in scripts analyzing species. These functions are 
# wrappers of exSituResample (or exSituResample.Par) sapplied over resampling replicates; these functions, in turn,
# are wrappers of calculateCoverage (sapplied over different sample sizes).

# The functions in this script are divided into sections based on their role in the workflow: the majority
# of the most relevant functions are within the 'BUILDING THE RESAMPLING ARRAY" section. 

library(adegenet)
library(terra)
library(parallel)

# ---- BUILDING THE RESAMPLING ARRAY ----
# WORKER FUNCTION: Create buffers around points, using specified projection. This function is used for 
# calculations of both geographic and ecological buffers; it does not include any area calculations. 
# The default point and buffer projections are Web Mercator 84 (WGS84), also used in many gap analysis workflows. 
createBuffers <- function(df, radius=1000, ptProj='+proj=longlat +datum=WGS84', 
                          buffProj='+proj=eqearth +datum=WGS84', boundary){
  # Turn occurrence point data into a SpatVector
  spat_pts <- vect(df, geom=c('decimalLongitude', 'decimalLatitude'), crs=ptProj)
  # Reproject spatial vector to the specified projection
  proj_df <- terra::project(spat_pts, buffProj)
  # Place buffer around each point, then dissolve into one polygon
  buffers <- terra::buffer(proj_df, width=radius)
  buffers <- terra::aggregate(buffers, dissolve = TRUE)
  # Clip by boundary, so buffers don't extend into the water
  boundary <- terra::project(boundary, buffProj)
  buffers_clip <- terra::crop(buffers, boundary)
  # Return buffer polygons
  return(buffers_clip)
}

# WORKER FUNCTION: Given a data.frame of geographic coordinates, a vector of sample names, a specified
# buffer size, and projections for the coordinate points and the buffers, calculate geographic coverage.
# The sampVect argument represents a vector of sample names, which is used to subset the totalWildPoints
# dataframe to create a separate "ex situ" spatial object. Then, the createBuffers function is used to 
# place buffers around all wild points and the sample, and then the proportion of the total area covered 
# is calculated
geo.compareBuff <- function(totalWildPoints, sampVect, buffSize, ptProj, buffProj, boundary, parFlag=FALSE){
  # If running in parallel: world polygon shapefile needs to be 'unwrapped', 
  # after being exported to cluster
  if(parFlag==TRUE){
    boundary <- unwrap(boundary)
  }
  # Select "ex situ" coordinates by subsetting totalWildPoints data.frame, according to sampVect
  exSitu <- totalWildPoints[sort(match(sampVect, totalWildPoints[,1])),]
  # Create buffers around selected (exSitu) wild points and around all (total) occurrences 
  geo_exSitu <- createBuffers(exSitu, buffSize, ptProj, buffProj, boundary)
  geo_total <- createBuffers(totalWildPoints, buffSize, ptProj, buffProj, boundary)
  # Calculate the area under both buffers. The 1,000,000 value converts values to km²
  geo_exSituArea <- expanse(geo_exSitu)/1000000
  geo_totalArea <- expanse(geo_total)/1000000
  # Calculate the proportion of the ex situ buffer areas to the total buffer area (percent geographic coverage)
  geo_Coverage <- (geo_exSituArea/geo_totalArea)*100
  return(geo_Coverage)
}

# WORKER FUNCTION: Analogous to geo.compareBuff, but instead of comparing buffered area around a random
# sample of points to the total buffered area, compares it instead to the total area under SDM, which
# is passed to this function as a raster argument (model).The sampVect argument represents vector 
# of sample names, which is used to subset the totalWildPoints data.frame to create a separate 
# "ex situ" spatial object. Then, the createBuffers function is used to place buffers around sampled
# points, and the proportion of the total area covered is calculated. Authored by Dan Carver
geo.compareBuffSDM <- function(totalWildPoints, sampVect, buffSize, model, ptProj, buffProj, 
                               boundary, parFlag=FALSE){
  # If running in parallel, unwrap spatial features 
  if(parFlag==TRUE){
    boundary <- unwrap(boundary)
    model <- unwrap(model)
  }
  # Generate a mask of the model layer by converting all 0 values to NA  
  m <- c(0, 0, NA)
  mask <- model |>
    classify(m)
  # Select "ex situ" coordinates by subsetting totalWildPoints data.frame, according to sampVect
  exSitu <- totalWildPoints[sort(match(sampVect, totalWildPoints[,1])),]
  # Create buffers around selected (exSitu) wild points and around all (total) occurrences 
  geo_exSitu <- createBuffers(exSitu, buffSize, ptProj, buffProj, boundary) |>
    # Reproject to match crs of smd object 
      terra::project(mask)
  # Rasterize 
  buffRast <- geo_exSitu |>
    terra::rasterize(mask)
  # Apply mask
  buffMask <- buffRast * mask
  # Calculate the area under both rasters in km²
  geo_exSituArea <- terra::cellSize(buffMask, mask= TRUE, unit = "km") |>
    terra::values()|>
    sum(na.rm = TRUE)
  geo_totalArea <- terra::cellSize(mask, mask= TRUE, unit = "km") |>
    terra::values()|>
    sum(na.rm = TRUE)
  # Calculate the proportion of the ex situ buffer areas to the total buffer area (percent geographic coverage)
  geo_Coverage <- (geo_exSituArea/geo_totalArea)*100
  return(geo_Coverage)
}

# WORKER FUNCTION: This function compares the resolution of a raster argument to the geographic buffer size 
# being used. If the buffer size is smaller than the raster resolution, a warning is given, and the SDM 
# is resampled (scaled) to a lower resolution. (Without this scaling, geographic coverage commands will 
# throw an error). Authored by Dan Carver.
geo.checkSDMres <- function(buffSize, raster, parFlag=FALSE){
  # If running in parallel, unwrap the raster argument
  if(parFlag==TRUE){
    raster <- unwrap(raster)
  }
  # Convert buffer from meters to degrees (assuming 1 m = 0.000012726903908907691 degrees)
  buffSizeDegree <- buffSize * 0.000012726903908907691
  # Extract the raster resolution
  sdmDegree <- terra::res(raster)[1]
  # If geographic buffer less than SDM resolution, give a warning and resample SDM to smaller resolution
  if(buffSizeDegree < sdmDegree){
    warning(paste0('SDM provided has a resolution larger than geographic buffer size (',
                buffSize, ' km). SDM will be resampled to a smaller resolution.'))
    # Resample the raster to a smaller cell size, based on scaling factor
    scaleFactor <- ceiling(sdmDegree / buffSizeDegree)
    raster <- terra::disagg(raster, scaleFactor) 
  }
  # If running in parallel, rewrap the raster argument
  if(parFlag==TRUE){
    raster <- wrap(raster)
  }
  # Return the raster
  return(raster)
}

# WORKER FUNCTION: Create a data.frame with ecoregion data extracted for area covered by buffers
eco.intersectBuff <- function(df, buffSize, ptProj, buffProj, ecoRegion, boundary, parFlag=FALSE){
  # If running in parallel: world polygon and ecoregions shapefiles need to be unwrapped, 
  # after being exported to cluster
  if(parFlag==TRUE){
    ecoRegion <- unwrap(ecoRegion) ; boundary <- unwrap(boundary)
  }
  # Create buffers
  buffers <- createBuffers(df, buffSize, ptProj, buffProj, boundary)
  # Make sure ecoregions are in same projection as buffers
  ecoProj <- terra::project(ecoRegion, buffProj)
  # Intersect buffers with ecoregions, and return
  ecoBuffJoin <- intersect(buffers, ecoProj)
  return(ecoBuffJoin)
}

# WORKER FUNCTION: Create a dataframe with ecoregion data extracted for area covered by buffers
# surrounding the sample ("exSitu") points. Then, compare the ecoregions count of the
# sample to the total, which is passed as an argument (ecoTotalCount) to this function. The layerType 
# argument allows for 3 possible values: US (EPA Level 4), NA (EPA Level 3), and GL 
# (TNC Global Terrestrial) ecoregions. These should correspond with the ecoRegion argument 
# (which specifies the ecoregion shapefile).
eco.compareBuff <- function(totalWildPoints, sampVect, buffSize, ecoTotalCount, ptProj, buffProj, 
                            ecoRegion, layerType=c('US','NA','GL'), boundary, parFlag=FALSE){
  # Match layerType argument, which specifies which ecoregion data type to extract (below)
  layerType <- match.arg(layerType)
  # If running in parallel: world polygon and ecoregions shapefiles need to be 'unwrapped', 
  # after being exported to cluster
  if(parFlag==TRUE){
    ecoRegion <- unwrap(ecoRegion) ; boundary <- unwrap(boundary)
  }
  # Build sample ex situ points by subseting totalWildPoints dataframe, according to sampVect
  exSitu <- totalWildPoints[sort(match(sampVect, totalWildPoints[,1])),]
  # Create dataframe of ecoregion-buffer intersections for ex situ points
  eco_exSitu <- eco.intersectBuff(exSitu, buffSize, ptProj, buffProj, ecoRegion, boundary)
  # Based on the ecoRegion shapefile and the specified layer type, count the number of ecoregions in
  # the random sample (exSitu) 
  if(layerType=='US'){
    # Extract the number of EPA Level IV ('U.S. Only') ecoregions
    eco_exSituCount <- length(unique(eco_exSitu$US_L4CODE))
  } else {
    # Extract the number of EPA Level III ('North America') ecoregions 
    if(layerType=='NA'){
      eco_exSituCount <- length(unique(eco_exSitu$NA_L3CODE))
    } else {
      # Extract the number of Nature Conservancy ('Global Terrestrial') ecoregions
      eco_exSituCount <- length(unique(eco_exSitu$ECO_ID_U))
    }
  }
  # Calculate difference in number of ecoregions between the sample and total (all data points), and return
  eco_Coverage <- (eco_exSituCount/ecoTotalCount)*100
  return(eco_Coverage)
}

# WORKER FUNCTION: Create a dataframe with ecoregion data extracted for area covered by buffers
# surrounding all points, and then count the number of ecoregions. This function is run once (per buffer size), 
# and it's final result (eco_totalCount) is passed down to a lower level function (eco.compareBuff) 
# to use as the denominator for ecological coverage calculations. The layerType argument allows 
# for 3 possible values: US (EPA Level 4), NA (EPA Level 3), and GL (TNC Global Terrestrial) 
# ecoregions. These should correspond with the ecoRegion argument (which specifies the ecoregion shapefile).
eco.totalEcoregionCount <- function(totalWildPoints, buffSize, ptProj, buffProj, 
                                    ecoRegion, layerType=c('US','NA','GL'), boundary, parFlag=FALSE){
  # Match layerType argument, which specifies which ecoregion data type to extract (below)
  layerType <- match.arg(layerType)
  # If running in parallel: world polygon and ecoregions shapefiles need to be 'unwrapped', 
  # after being exported to cluster
  if(parFlag==TRUE){
    ecoRegion <- unwrap(ecoRegion) ; boundary <- unwrap(boundary)
  }
  # Create dataframe of ecoregion-buffer intersections for all points
  eco_total <- eco.intersectBuff(totalWildPoints, buffSize, ptProj, buffProj, ecoRegion, boundary)
  # Based on the ecoRegion shapefile and the specified layer type, count the number of ecoregions in
  # all of the data points
  if(layerType=='US'){
    # Extract the number of EPA Level IV ('U.S. Only') ecoregions
    eco_totalCount <- length(unique(eco_total$US_L4CODE))
  } else {
    # Extract the number of EPA Level III ('North America') ecoregions 
    if(layerType=='NA'){
      eco_totalCount <- length(unique(eco_total$NA_L3CODE))
    } else {
      # Extract the number of Nature Conservancy ('Global Terrestrial') ecoregions
      eco_totalCount <- length(unique(eco_total$ECO_ID_U))
    }
  }
  # Return the total number of ecoregions
  return(eco_totalCount)
}

# WORKER FUNCTION: Function for reporting representation rates, using a vector of allele frequencies 
# and a sample matrix. Assumes that freqVector represents the absolute allele frequencies 
# for the population of interest (the entire wild population). 
# 1. The length of matches between garden and wild alleles is calculated (numerator). 
# 2. The complete number of wild alleles of that category (denominator) is calculated. 
# 3. From these 2 values, a percentage is calculated. 
# This function returns the numerators, denominators, and the proportion (representation rates) in a matrix.
gen.getAlleleCategories <- function(freqVector, sampleMat){
  # Determine how many Total alleles in the sample matrix are found in the frequency vector 
  exSitu_allAlleles <- length(which(names(freqVector) %in% colnames(sampleMat)))
  total_allAlleles <- length(freqVector)
  allPercentage <- (exSitu_allAlleles/total_allAlleles)*100
  # Very common alleles (greater than 10%)
  exSitu_vComAlleles <- length(which(names(which(freqVector > 10)) %in% colnames(sampleMat)))
  total_vComAlleles <- length(which(freqVector > 10))
  vComPercentage <- (exSitu_vComAlleles/total_vComAlleles)*100
  # Common alleles (greater than 5%)
  exSitu_comAlleles <- length(which(names(which(freqVector > 5)) %in% colnames(sampleMat)))
  total_comAlleles <- length(which(freqVector > 5))
  comPercentage <- (exSitu_comAlleles/total_comAlleles)*100
  # Low frequency alleles (between 1% and 10%)
  exSitu_lowFrAlleles <- length(which(names(which(freqVector < 10 & freqVector > 1)) %in% colnames(sampleMat)))
  total_lowFrAlleles <- length(which(freqVector < 10 & freqVector > 1))
  lowFrPercentage <- (exSitu_lowFrAlleles/total_lowFrAlleles)*100
  # Rare alleles (less than 1%)
  exSitu_rareAlleles <- length(which(names(which(freqVector < 1)) %in% colnames(sampleMat)))
  total_rareAlleles <- length(which(freqVector < 1))
  rarePercentage <- (exSitu_rareAlleles/total_rareAlleles)*100
  # Concatenate values to vectors
  exSituAlleles <- c(exSitu_allAlleles, exSitu_vComAlleles, exSitu_comAlleles, exSitu_lowFrAlleles, exSitu_rareAlleles)
  totalWildAlleles <- c(total_allAlleles, total_vComAlleles, total_comAlleles, total_lowFrAlleles, total_rareAlleles)
  repRates <- c(allPercentage,vComPercentage,comPercentage,lowFrPercentage,rarePercentage) 
  # Bind vectors to a matrix, name dimensions, and return
  exSituValues <- cbind(exSituAlleles, totalWildAlleles, repRates)
  rownames(exSituValues) <- c('Total','V. common','Common', 'Low freq.','Rare')
  colnames(exSituValues) <- c('Ex situ', 'Total', 'Rate (%)')
  return(exSituValues)
}

# Declare a function which, given a genind object, builds a matrix of Euclidean distances between every individual
gen.buildDistMat <- function(genObj){
  browser()
  # Convert genind object to data.frame, ignoring any population designations
  df <- genind2df(genObj, usepop=FALSE)
  # Convert data.frame to genlight object
  genLight <- as.genlight(df)
  # Build a genetic distance matrix from the genlight object. The default method (Euclidean) is used.
  distMat <- dist(genLight)
  # Check that total genetic distance is a number; if it isn't, genetic distance calculations will fail
  if(is.na(sum(c(distMat)))){
    stop('Total genetic distance is NA, so genetic distance coverage cannot be calculated!')}
  return(distMat)
}

# Given a genetic distance matrix and a vector of sample names, calculate the sum total of the genetic 
# distances between all pairs of individuals (denominator) and the sum of the genetic distances 
# strictly between sampled individuals (numerator), and returns a coverage metric
gen.calcGenDistCov <- function(distMat, sampVect){
  # Calculate the total of all the genetic distances (denominator)
  distTotal <- sum(c(distMat))
  # Subset the genetic distance matrix to strictly samples included within the vector of sample names
  distMatSamp <- dist_subset(distMat, sampVect)
  # Calculate the total genetic distance in the subset matrix (numerator)
  distSample <- sum(c(distMatSamp))
  # Calculate the percent coverage of genetic distance, and return
  distCov <- (distSample/distTotal)*100
  return(distCov)
}

# CORE FUNCTION: Wrapper of gen.getAlleleCategories, geo.compareBuff, and eco.compareBuff worker functions. 
# Given a genetic matrix (rows are samples, columns are alleles) and a dataframe of coordinates 
# (3 columns: sample names, latitudes, and longitudes), it calculates the genetic,
# geographic (if flagged), and ecologcial (if flagged) coverage from a random draw of some amount of 
# samples (numSamples). The sample names between the genind object and the coordinate dataframe need
# to match (in order to properly subset across genetic, geographic, and ecological datasets). 

# The SDMrast argument allows users to calculate geographic coverage using a 
# rasterized SDM provided to the function (in addition to the standard approach, 
# which involves using buffered areas). Both geoBuff and ecoBuff can be single values or
# a vector of values, in which case multiple geographic/ecological coverages are calculated
# according to each buffer value.
calculateCoverage <- function(genMat, genDistMat=NA, geoFlag=TRUE, coordPts, geoBuff, 
                              SDMrast=NA, ptProj='+proj=longlat +datum=WGS84', 
                              buffProj='+proj=eqearth +datum=WGS84', boundary, 
                              ecoFlag=FALSE, ecoBuff, ecoTotalCount, ecoRegions, 
                              ecoLayer=c('US','NA','GL'), parFlag=FALSE, numSamples){
  # GENETIC PROCESSING
  # Calculate a vector of allele frequencies, based on the total sample matrix
  freqVector <- colSums(genMat, na.rm = TRUE)/(nrow(genMat)*2)*100
  # Remove any missing alleles (those with frequencies of 0) from the frequency vector
  freqVector <- freqVector[which(freqVector != 0)]
  # From matrix of individuals, select a random set (rows). This is the set of individuals that will
  # be used for all downstream coverage calculations within this function (genetic, geographic, and ecological)
  samp <- genMat[sample(nrow(genMat), size=numSamples, replace = FALSE),]
  # Remove any missing alleles (those with colSums of 0) from the sample matrix
  samp <- samp[,which(colSums(samp, na.rm = TRUE) != 0)]
  # Genetic coverage: calculate sample's allelic representation
  genRates <- gen.getAlleleCategories(freqVector, samp)
  # Check if genetic distance matrix was passed down by upper level functions;
  # if so, calculate coverages using a distance metric (in addition to allelic coverage)
  if(class(genDistMat)=='logical'){
    # Subset matrix returned by getAlleleCategories to just 3rd column (representation rates), and return
    genRates <- genRates[,3]
  } else {
    # Pass distance matrix and sample name vector to function calculating  
    # proportion of total pairwise genetic distances represented in sample
    genDistCov <- gen.calcGenDistCov(distMat=genDistMat, sampVect=rownames(samp))
    # Append the resulting coverage value to the allelic coverages
    genRates <- c(genRates[,3], genDistCov)
    names(genRates)[6] <- 'GenDist'
  }
  
  # GEOGRAPHIC PROCESSING
  if(geoFlag==TRUE){
    # Check that sample names in genetic matrix match the column of sample names in the coordinate dataframe
    if(!identical(rownames(genMat), coordPts[,1])){
      stop('Error: Sample names between the genetic matrix and the first 
         column of the coordinate point data.frame do not match.')
    }
    # Geographic coverage: for each buffer size, calculate sample's geograhpic representation, by 
    # passing all points (coordPts) and the random subset of points (rownames(samp)) to the 
    # geo.compareBuff worker function, which will calculate the proportion of area covered in the 
    # random sample. 
    geoRates <- 
      lapply(geoBuff, function(x) geo.compareBuff(totalWildPoints=coordPts, sampVect=rownames(samp),
                                                  buffSize=x, ptProj=ptProj, buffProj=buffProj, 
                                                  boundary=boundary, parFlag=parFlag))
    # If no rasterized SDM is provided, calculate geographic coverage using just the buffer approach (default)
    if(class(SDMrast)=='logical'){
      names(geoRates) <- paste0(rep('Geo_Buff_',), geoBuff/1000, 'km')
    } else {
      # If rasterized SDM provided, calculate geo. coverage using SDM approach, and append coverages. 
      # Using mapply to iterate over multiple lists (buffer sizes and SDM rasters)
      geoRates_SDM <- 
        mapply(function(b,r) geo.compareBuffSDM(totalWildPoints=coordPts, sampVect=rownames(samp),
                                                buffSize=b, model=r, ptProj=ptProj, 
                                                buffProj=buffProj, boundary=boundary, 
                                                parFlag=parFlag), b=geoBuff, r=SDMrast)
      # Combine the geographic coverage values (total buffered area approach and SDM approach)
      geoRates <- c(geoRates, geoRates_SDM)
      # Name geographic coverage values, according to buffer size
      names(geoRates) <- c(paste0(rep('Geo_Buff_',), geoBuff/1000, 'km'),
                           paste0(rep('Geo_SDM_',), geoBuff/1000, 'km'))
    }
  } else {
    # If geographic processing is not occurring, make coverage values NA
    geoRates <- NA ; names(geoRates) <- 'Geo'
  }
  
  # ECOLOGICAL PROCESSING
  if(ecoFlag==TRUE){
    # Ecological coverage: for each buffer size, calculate sample's ecological representation, by passing all 
    # points (coordPts) and the random subset of points (rownames(samp)) to the eco.compareBuff worker function, 
    # which will calculate the proportion of ecoregions covered in the random sample. The ecoTotalCount argument
    # provides the denominator (ecoregions covered by all points) used for coverage calculations.
    ecoRates <- 
      mapply(function(b,t) eco.compareBuff(totalWildPoints=coordPts, sampVect=rownames(samp),
                                           buffSize=b, ecoTotalCount=t, ptProj=ptProj, buffProj=buffProj, 
                                           ecoRegion=ecoRegions, layerType=ecoLayer, boundary=boundary, 
                                           parFlag=parFlag), b=ecoBuff, t=ecoTotalCount)
    # Name ecological coverage values, according to buffer size
    names(ecoRates) <- paste0(rep('Eco_Buff_',), ecoBuff/1000, 'km')
  } else {
    # If ecological processing is not occurring, make coverage values NA
    ecoRates <- NA ; names(ecoRates) <- 'Eco'
  }
  
  # Combine genetic, geographic, and ecological coverage rates into a vector (not a list), and return
  covRates <- unlist(c(genRates, geoRates, ecoRates))
  return(covRates)
}

# WRAPPER FUNCTION: iterates calculateCoverage over the entire matrix of samples
exSituResample <- function(genMat, genDistMat=NA, geoFlag=TRUE, coordPts, geoBuff=50000, SDMrast=NA, 
                           ptProj='+proj=longlat +datum=WGS84', buffProj='+proj=eqearth +datum=WGS84', 
                           boundary, ecoFlag=FALSE, ecoBuff=50000, ecoTotalCount, ecoRegions, ecoLayer='US', parFlag){
  # Apply the calculateCoverage function to all rows of the wild matrix
  # (except row 1, because we need at least 2 individuals to sample)
  # The resulting matrix needs to be transposed, in order to keep columns as different coverage categories
  cov_matrix <- 
    t(sapply(2:nrow(genMat), 
             function(x) calculateCoverage(genMat=genMat, genDistMat=genDistMat, geoFlag=geoFlag, coordPts=coordPts, 
                                           geoBuff=geoBuff, SDMrast=SDMrast, ptProj=ptProj, 
                                           buffProj=buffProj, boundary=boundary, ecoFlag=ecoFlag, 
                                           ecoBuff=ecoBuff, ecoTotalCount=ecoTotalCount, ecoRegions=ecoRegions, 
                                           ecoLayer=ecoLayer, parFlag=FALSE, numSamples=x), simplify = 'array'))
  # Return the matrix of coverage values
  return(cov_matrix)
}

# WRAPPER FUNCTION: iterates calculateCoverage over the entire matrix of samples, in parallel
exSituResample.Par <- function(genMat, genDistMat=NA, geoFlag=TRUE, coordPts, geoBuff=50000, SDMrast=NA, 
                               ptProj='+proj=longlat +datum=WGS84', buffProj='+proj=eqearth +datum=WGS84', 
                               boundary, ecoFlag=FALSE, ecoBuff=50000, ecoTotalCount, ecoRegions, ecoLayer='US', 
                               parFlag=TRUE, cluster){
  # Apply the calculateCoverage function to all rows of the wild matrix using parSapply
  # (except row 1, because we need at least 2 individuals to sample)
  # The resulting matrix needs to be transposed, in order to keep columns as different coverage categories
  cov_matrix <-
    t(parSapply(cluster, 2:nrow(genMat), 
                function(x) calculateCoverage(genMat=genMat, genDistMat=genDistMat, geoFlag=geoFlag, 
                                              coordPts=coordPts, geoBuff=geoBuff, SDMrast=SDMrast, ptProj=ptProj, 
                                              buffProj=buffProj, boundary=boundary, ecoFlag=ecoFlag, 
                                              ecoBuff=ecoBuff, ecoTotalCount=ecoTotalCount, ecoRegions=ecoRegions, 
                                              ecoLayer=ecoLayer, parFlag=parFlag, numSamples=x), simplify = 'array'))
  # Return the matrix of coverage values
  return(cov_matrix)
}

# WRAPPER FUNCTION: iterates exSituResample, which will generate an array of values from a single genind object.
# Checks for arguments are also performed; if ecological coverage is being calculated, then the total number
# of ecoregions will . This function doesn't run in parallel, so it's  primarily used
# for testing/demonstration purposes
geo.gen.Resample <- function(genObj, genDistFlag=FALSE, geoFlag=TRUE, coordPts, geoBuff=50000, SDMrast=NA,
                             ptProj='+proj=longlat +datum=WGS84', buffProj='+proj=eqearth +datum=WGS84',
                             boundary, ecoFlag=FALSE, ecoBuff=50000, ecoRegions,
                             ecoLayer=c('US', 'NA', 'GL'), reps=5){
  # Extract the genetic matrix from the genind object
  genMat <- genObj@tab
  # If genetic distance flag is set to TRUE, build the matrix of genetic distances to pass down to lower functions
  if(genDistFlag==TRUE){
    cat(paste0('\n', '- genDistFlag ON: will calculate genetic distance coverage -'))
    genDistMat <- gen.buildDistMat(genObj=genObj)
    # Otherwise, set genDistMat as NA
  } else {
    genDistMat <- NA
  } 
  # Otherwise, set genDistMat as NA
  # If calculating geographic coverage, check for arguments
  if(geoFlag==TRUE){
    # Check for the required arguments (ptProj and buffProj will use defaults, if not specified)
    if(missing(coordPts)) stop('For geographic coverage, a data.frame of wild coordinates (coordPts) is required')
    if(missing(geoBuff)) stop('For geographic coverage, an integer (or vector of integers) specifying the geographic 
                              buffer size(s) is required (geoBuff argument)')
    if(missing(boundary)) stop('For geographic coverage, a SpatVector object of country boundaries (boundary) is required')
    # Check that the names of the latitude and longitude columns are properly written (this is unfortunately hard-coded)
    if(!identical(colnames(coordPts)[2:3], c('decimalLatitude', 'decimalLongitude'))){
    stop('The column names of the geographic coordinates dataframe (coordPts) need to be
         decimalLatitude and decimalLongitude. Please rename your dataframe of geographic coordinates!')
      }
    # Print out message stating what coverages are being calculated, and how many buffer sizes
    cat('\n', '- geoFlag ON: will calculate geographic coverage (total buffer) -')
    cat(paste0('\n', '--- Number of buffer sizes (Geo, Total buffer): ', length(geoBuff), ' ---'))
    # If SDM is provided (meaning it's not NA, or class logical): 
    if(!class(SDMrast)=='logical'){
      # Print out message stating what coverages are being calculated, and how many buffer sizes
      cat(paste0('\n', '- SDM provided: will calculate geographic coverage (SDM) -'))
      cat(paste0('\n', '--- Number of buffer sizes (Geo, SDM): ', length(geoBuff), ' ---'))
      # Check that geographic buffer size is greater than SDM raster resolution,
      # and fix if not. If multiple buffer sizes are used, resample the resolution of the SDM according
      # multiple times, and return a list
      SDMrast <- lapply(geoBuff, function(x) geo.checkSDMres(buffSize=x, raster=SDMrast, parFlag=FALSE))
    }
  }
  # If calculating ecological coverage, check for arguments
  if(ecoFlag==TRUE){
    # Match ecoLayer argument, to ensure it is 1 of 3 possible values ('US', 'NA', 'GL)
    ecoLayer <- match.arg(ecoLayer)
    # Check for the required arguments (ptProj, buffProj, and ecoLayer will use defaults, if not specified)
    if(missing(coordPts)) stop('For ecological coverage, a data.frame of wild coordinates (coordPts) is required')
    if(missing(ecoBuff)) stop('For ecological coverage, an integer (or vector of integers) specifying
                            the ecological buffer size(s)  is required (ecoBuff argument)')
    if(missing(ecoRegions)) stop('For ecological coverage, a SpatVector object of ecoregions (ecoregions) is required')
    if(missing(boundary)) stop('For ecological coverage, a SpatVector object of country boundaries (boundary) is required')
    # Print out message stating what coverages are being calculated, and how many buffer sizes
    cat(paste0('\n', '- ecoFlag ON: will calculate ecological coverage -'))
    cat(paste0('\n', '--- Number of buffer sizes (Eco): ', length(ecoBuff), ' ---'))
    # CALCULATE TOTAL ECOLOGICAL COVERAGE: calculate the number of ecoregions found under all samples for all
    # buffer sizes. The resulting variable is passed down to lower level functions to optimize processing
    cat(paste0('\n', '--- CALCULATING TOTAL ECOREGION COVERAGE... ---'))
    ecoTotalCount <- lapply(ecoBuff,
                            function(x) eco.totalEcoregionCount(totalWildPoints=coordPts, buffSize=x,
                                                                ptProj=ptProj, buffProj=buffProj,
                                                                ecoRegion=ecoRegions, layerType=ecoLayer,
                                                                boundary=boundary, parFlag=FALSE))
  }
  # Print starting time
  startTime <- Sys.time()
  cat(paste0('\n', '%%% RESAMPLING START: ', startTime, '\n'))
  # Run resampling for all replicates, using sapply and lambda function
  resamplingArray <-
    sapply(1:reps, function(x) exSituResample(genMat=genMat, genDistMat=genDistMat, geoFlag=geoFlag,
                                              coordPts=coordPts, geoBuff=geoBuff,
                                              SDMrast=SDMrast, ptProj=ptProj,
                                              buffProj=buffProj,boundary=boundary,
                                              ecoFlag=ecoFlag, ecoBuff=ecoBuff,
                                              ecoTotalCount=ecoTotalCount, ecoRegions=ecoRegions,
                                              ecoLayer=ecoLayer, parFlag=FALSE), simplify = 'array')
    # Print ending time and total runtime
    endTime <- Sys.time()
    cat(paste0('\n', '%%% RESAMPLING END: ', endTime))
    cat(paste0('\n', '%%% TOTAL RUNTIME: ', endTime-startTime))
    # Return array
    return(resamplingArray)
}

# WRAPPER FUNCTION: iterates exSituResample.Par, which will generate an array of values from a single genind object
# This function iterates the parallelized version of exSituResample, such that each different sample size for a 
# single resampling replicate is processed on a single core. Results (resampling array) are saved to a specified file path.
geo.gen.Resample.Par <- function(genObj, genDistFlag=FALSE, geoFlag=TRUE, coordPts, geoBuff=50000, 
                                 SDMrast=NA, ptProj='+proj=longlat +datum=WGS84', 
                                 buffProj='+proj=eqearth +datum=WGS84', boundary, ecoFlag=FALSE, 
                                 ecoBuff=50000, ecoRegions, ecoLayer=c('US','NA','GL'), reps=5,
                                 arrayFilepath='~/resamplingArray.Rdata', cluster){
  # Extract the genetic matrix from the genind object
  genMat <- genObj@tab
  # If genetic distance flag is set to TRUE, build the matrix of genetic distances to pass down to lower functions
  if(genDistFlag==TRUE){
    cat(paste0('\n', '- genDistFlag ON: will calculate genetic distance coverage -'))
    genDistMat <- gen.buildDistMat(genObj=genObj)
    # Otherwise, set genDistMat as NA
  } else {
    genDistMat <- NA
  }
  # If calculating geographic coverage, check for arguments
  if(geoFlag==TRUE){
    # Check for the required arguments (ptProj and buffProj will use defaults, if not specified)
    if(missing(coordPts)) stop('For geographic coverage, a data.frame of wild coordinates (coordPts) is required')
    if(missing(geoBuff)) stop('For geographic coverage, an integer (or vector of integers) specifying the geographic 
                              buffer size(s) is required (geoBuff argument)')
    if(missing(boundary)) stop('For geographic coverage, a SpatVector object of country boundaries (boundary) is required')
    # Check that the names of the latitude and longitude columns are properly written (this is unfortunately hard-coded)
    if(!identical(colnames(coordPts)[2:3], c('decimalLatitude', 'decimalLongitude'))){
      stop('The column names of the geographic coordinates dataframe (coordPts) need to be 
           decimalLatitude and decimalLongitude. Please rename your dataframe of geographic coordinates!')
    }
    # Print out message stating what coverages are being calculated, and how many buffer sizes
    cat('\n', '- geoFlag ON: will calculate geographic coverage (total buffer) -')
    cat(paste0('\n', '--- Number of buffer sizes (Geo, Total buffer): ', length(geoBuff), ' ---'))
    # If SDM is provided (meaning it's not NA, or class logical):
    if(!class(SDMrast)=='logical'){
      # Print out message stating what coverages are being calculated, and how many buffer sizes
      cat(paste0('\n', '- SDM provided: will calculate geographic coverage (SDM) -'))
      cat(paste0('\n', '--- Number of buffer sizes (Geo, SDM): ', length(geoBuff), ' ---'))
      # Check that geographic buffer size is greater than SDM raster resolution, 
      # and fix if not. If multiple buffer sizes are used, resample the resolution of the SDM according
      # multiple times, and return a list
      SDMrast <- lapply(geoBuff, function(x) geo.checkSDMres(buffSize=x, raster=SDMrast, parFlag=TRUE))
      # Export list of raster objects to the cluster
      clusterExport(cl=cluster, varlist='SDMrast', envir=environment())
    }
  }
  # If calculating ecological coverage, check for arguments
  if(ecoFlag==TRUE){
    # Match ecoLayer argument, to ensure it is 1 of 3 possible values ('US', 'NA', 'GL)
    ecoLayer <- match.arg(ecoLayer)
    # Check for the required arguments (ptProj, buffProj, and ecoLayer will use defaults, if not specified)
    if(missing(coordPts)) stop('For ecological coverage, a data.frame of wild coordinates (coordPts) is required')
    if(missing(ecoBuff)) stop('For ecological coverage, an integer (or vector of integers) specifying 
                              the ecological buffer size(s)  is required (ecoBuff argument)')
    if(missing(ecoRegions)) stop('For ecological coverage, a SpatVector object of ecoregions (ecoregions) is required')
    if(missing(boundary)) stop('For ecological coverage, a SpatVector object of country boundaries (boundary) is required')
    # Print out message stating what coverages are being calculated, and how many buffer sizes
    cat(paste0('\n', '- ecoFlag ON: will calculate ecological coverage -'))
    cat(paste0('\n', '--- Number of buffer sizes (Eco): ', length(ecoBuff), ' ---'))
    # CALCULATE TOTAL ECOLOGICAL COVERAGE: calculate the number of ecoregions found under all samples for all 
    # buffer sizes, and pass this down to lower level functions, to optimize processing
    cat(paste0('\n', '--- CALCULATING TOTAL ECOREGION COVERAGE... ---'))
    ecoTotalCount <- lapply(ecoBuff, 
                            function(x) eco.totalEcoregionCount(totalWildPoints=coordPts, buffSize=x,
                                                                ptProj=ptProj, buffProj=buffProj, 
                                                                ecoRegion=ecoRegions, layerType=ecoLayer,
                                                                boundary=boundary, parFlag=TRUE))
  }
  # Print starting time
  startTime <- Sys.time() 
  cat(paste0('\n', '%%% RESAMPLING START: ', startTime, '\n'))
  # Run resampling for all replicates, using sapply and lambda function
  resamplingArray <- 
    sapply(1:reps, function(x) exSituResample.Par(genMat=genMat, genDistMat=genDistMat, geoFlag=geoFlag, 
                                                  coordPts=coordPts, geoBuff=geoBuff, 
                                                  SDMrast=SDMrast, ptProj=ptProj, 
                                                  buffProj=buffProj, boundary=boundary, 
                                                  ecoFlag=ecoFlag, ecoBuff=ecoBuff,
                                                  ecoTotalCount=ecoTotalCount, ecoRegions=ecoRegions, 
                                                  ecoLayer=ecoLayer, parFlag=TRUE, cluster), simplify = 'array')
  # Print ending time and total runtime
  endTime <- Sys.time() 
  cat(paste0('\n', '%%% RESAMPLING END: ', endTime))
  cat(paste0('\n', '%%% TOTAL RUNTIME: ', endTime-startTime))
  # Save the resampling array object to disk, for later usage
  saveRDS(resamplingArray, file = arrayFilepath)
  cat(paste0('\n', '%%% Resampling array object saved to: ', arrayFilepath, '\n'))
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
  meanValue <- min(which(apply(resamplingArray[,2,],1, mean, na.rm=TRUE) > 95))
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

# Function for calculating normalized root mean square error. Takes two vectors of equal length
# (for this project, typically geographic or ecological coverages and allelic representation values)
# and calculates the root mean square error between them, and then normalizes that value based on the 
# 'type' argument. A lower value indicates similarity between values. In the context of this project, 
# 'obs_var' is the explanatory variable (what we're using as a 'proxy' for genetic coverage), 
# and 'pred_var' is what we're trying to predict (genetic coverage).
nrmse.func <-  function(obs_var, pred_var, norm_type='mean'){
  # Check that lengths of observed and predicted variables match, and error if not
  if(length(obs_var) != length(pred_var)) stop('Lengths of observed and predicted variables do not match!')
  # Calculate root mean square error 
  squared_sums <- sum((obs_var - pred_var)^2)
  mse <- squared_sums/length(obs_var)
  rmse <- sqrt(mse)
  # Check that type argument matches preset options; if not, return the non-normalized RMSE value
  if (!norm_type %in% c('mean', 'sd', 'maxmin', 'iq')){
    message('Wrong type argument for how to normalize! Non-normalized root mean square error value returned.')
    rmse <- round(rmse, 3)
    return(rmse)
  } else {
    # Normalize RMSE, based on type argument
    if (norm_type == 'sd') nrmse <- rmse/sd(obs_var)
    if (norm_type == 'mean') nrmse <- rmse/mean(obs_var)
    if (norm_type == 'maxmin') nrmse <- rmse/(max(obs_var) - min(obs_var))
    if (norm_type == 'iq') nrmse <- rmse/(quantile(obs_var, 0.75) - quantile(obs_var, 0.25))
    # Round result and return
    nrmse <- round(nrmse, 5)
    return(nrmse)
  }
}

# Function for building a matrix of normalized root mean square error (NRMSE) values based on 
# a data.frame of resampling values and the genetic coverage type used for the predictor variable.
# Given these arguments, a matrix is generated which has a NRMSE value for each set of geographic/ecological
# coverage values compared to the genetic coverage metric. This command expects the resampling data.frame
# to have certain row names and column names, and makes those checks
buildNRMSEmatrix <- function(resampDF, genCovType=c('CV', 'GD'), sdmFlag=TRUE){
  # Match argument for type of response variable (genetic coverage approach) to use
  genCovType <- match.arg(genCovType)
  # If CV is chosen, set the predictive variable argument to the 'Total' column in the data.frame
  if(genCovType=='CV'){
    predVar <- resampDF$Total
    # Otherwise, set the predictve variable to the 'GenDist' column
  } else {
    predVar <- resampDF$GenDist
  }
  # Extract buffer sizes through column names, for the geographic (total buffer) approach. This is 
  # a complicated command, but essentially relies on the column names of the resampling data.frame
  # to refer to the geographic buffer size used for that row of data.
  buffSizes <- 
    1000*(as.numeric(sub('km','', sapply(strsplit(grep('Geo_Buff_', colnames(resampDF), value=TRUE), '_'),'[',3))))
  # Build matrix of NRMSE values. Depending on the sdmFlag, this matrix will either have
  # 3 columns or just 2. 
  if(sdmFlag==TRUE){
    # Check that the resampling data.frame has the same number of columns for each coverage type
    if(length(grep('Geo_Buff_', colnames(resampDF)))!=length(grep('Geo_SDM_', colnames(resampDF)))){
      stop('Different number of Geo_Buff and Geo_SDM columns in data.frame!')
    }
    if(length(grep('Geo_Buff_', colnames(resampDF)))!=length(grep('Eco_Buff_', colnames(resampDF)))){
      stop('Different number of Geo_Buff and Eco_Buff columns in data.frame!')
    }
    # Name matrix columns accordingly
    NRMSEmat <- matrix(NA, nrow=length(buffSizes), ncol=3)
    colnames(NRMSEmat) <- c('Geo_Buff','Geo_SDM','Eco_Buff')
  } else {
    # Check that the resampling data.frame has the same number of columns for each coverage type
    if(length(grep('Geo_Buff_', colnames(resampDF)))!=length(grep('Eco_Buff_', colnames(resampDF)))){
      stop('Different number of Geo_Buff and Eco_Buff columns in data.frame!')
    }
    # Name matrix columns accordingly
    NRMSEmat <- matrix(NA, nrow=length(buffSizes), ncol=2)
    colnames(NRMSEmat) <- c('Geo_Buff','Eco_Buff')
  }
  # Name matrix rows according to buffer sizes
  rownames(NRMSEmat) <- paste0(buffSizes/1000, 'km')
  # Not all resampling data.frames have the geographic coverage results starting in the same column.
  # To address this, make the loop start at an index according to the column names
  startCol <- min(grep('Geo_Buff_', colnames(resampDF)))
  # Loop through the dataframe columns. The first three columns are skipped, as they're sample number and the
  # predictive variables (allelic coverage and genetic distance proportions)
  for(i in startCol:ncol(resampDF)){
    # Calculate NRMSE for the current column in the dataframe
    NRMSEvalue <- nrmse.func(resampDF[,i], pred = predVar)
    # Get the name of the current dataframe column
    dataName <- unlist(strsplit(names(resampDF)[[i]],'_'))
    # Match the data name to the relevant rows/columns of the receiving matrix
    matRow <- which(rownames(NRMSEmat) == dataName[[3]])
    matCol <- which(colnames(NRMSEmat) == paste0(dataName[[1]],'_',dataName[[2]]))
    # Locate the NRMSE value accordingly
    NRMSEmat[matRow,matCol] <- NRMSEvalue
  }
  # Rename the columns based on the genetic coverage type used
  if(genCovType=='CV'){
    colnames(NRMSEmat) <- paste0(colnames(NRMSEmat),'_CV')
  } else {
    colnames(NRMSEmat) <- paste0(colnames(NRMSEmat),'_GD')
  }
  # Return the resulting matrix
  return(NRMSEmat)
}

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

# ---- SDM GEOGRAPHIC COVERAGE FUNCTIONS ----
# The functions in this section are associated with the SDM approach to calculated geographic coverage,
# but are not used in the geographic genetic resampling process. Instead, they are used to process the
# data layers (e.g. country borders) prior to resampling, as well as to explore the results of geographic 
# resampling analyses (using both the total buffer and the SDM approach). These functions were written by
# Dan Carver.

#' grabWorldAdmin -- Dan Carver
#'
#' @param GeoGenCorr_wd : current working directory define by either variable or getwd()
#' @param fileExtentsion : choice between shp and gpkg.Prefer gpkg for single file. 
#' @param overwrite : defaults to false; True will force a redownload of the file.
#'
#' @description
#' Checks to see if the file exists. If true it is loaded. If false it is downloaded.
#' Allow uses to specify the file extention. .shp or .gpkg are options
#'
#' @return terra vect object of the world admin layer
#' 
grabWorldAdmin <- function(GeoGenCorr_wd, fileExtentsion, overwrite=FALSE){
  path <- file.path(paste0(GeoGenCorr_wd,
                           'GIS_shpFiles/world_countries_10m/world_countries_10m',
                           fileExtentsion))
  # test for presence of file. 
  if(file.exists(path) | overwrite == TRUE){
    # read in for return 
    world_poly_clip <- terra::vect(path)
  }else{
    # download data from natural earth 
    download <- rnaturalearth::ne_countries(scale = 10,
                                            type = "countries",
                                            returnclass = "sf")
    # convert to terra vect
    world_poly_clip <- download |>
      terra::vect()
    
    # export the file. 
    ## create folder structure 
    dir <- file.path(paste0(GeoGenCorr_wd,'GIS_shpFiles/world_countries_10m'))
    if(!dir.exists(dir)){
      dir.create(dir, recursive = TRUE)
    }
    # export the file
    terra::writeVector(x = world_poly_clip, filename = path)
  }
  return(world_poly_clip)
}

#' prepWorldAdmin -- Dan Carver
#'
#' @param worldAdmin : full world admin feature
#' @param occurranceData : tabular records with lat long of know species occurrances 
#' 
#' @description
#' We buffer a convex hull of the occurrance data for the species to 50km. This is used to 
#' Select the countries that overlap this area. These countries are filterer out of the larger
#' admin file and dissolve so there are no internal political boundaries.  
#' 
#' @return A geographic limits and disvolved verison of the world admin layer that only includes elements that intersect the occurrences
#' 
prepWorldAdmin <- function(world_poly_clip, wildPoints){
  
  # buffer the wild points to max expected buffer distance 
  bbox <- vect(wildPoints, 
               geom = c("decimalLongitude", "decimalLatitude"),
               crs = crs(world_poly_clip))|> # setting wgs1984
    terra::convHull() |>
    terra::buffer(width = 50000)
  crs(bbox) <- crs(world_poly_clip)
  
  # Unit is meter if x has a longitude/latitude CRS
  # intersect with the admin feature   
  uniqueLocs <- terra::intersect(world_poly_clip, bbox) 
  # pull the iso3 for a filter later on 
  counties <- unique(uniqueLocs$adm0_a3)
  # filter world admin to countries on overlap
  admin <- world_poly_clip |>
    sf::st_as_sf()|>
    dplyr::filter(adm0_a3 %in% counties)|>
    terra::vect()|> 
    terra::aggregate()
  
  # return the simplified object 
  return(admin)
}

#' makeAMap -- Dan Carver
#'
#' @param points : sf point object  
#' @param raster : terra raster object... converted to raster within the function. 
#' @param buffer : optional sf point buffered feature. 
#'
#' @return one of two maps depending on if a buffer input object was defined or not.
makeAMap <- function(points,raster,buffer=NA){
  # Create the centroid
  centroid <- points |>
    dplyr::mutate(group = 1)|>
    group_by(group) |>
    summarize(geometry = st_union(geometry)) |>
    st_centroid()|>
    st_coordinates()
  
  # Define base map 
  map1 <- leaflet(options = leafletOptions(minZoom = 4)) |>
    # Set zoom levels
    setView(lng = centroid[1]
            , lat = centroid[2]
            , zoom = 6) |>
    # Tile providers 
    addProviderTiles("OpenStreetMap", group = "OpenStreetMap") |>
    # Add point features 
    addCircleMarkers(
      data = points,
      color = "#4287f5",
      radius = 0.2,
      group = "Points",
      # Add highlight options to make labels a bit more intuitive 
    ) |> 
    addRasterImage(
      x = raster::raster(raster),
      colors = c("#ffffff95", "#8bed80"),
      group = "Raster"
    )|>
    addLegend(
      position = "topright",
      colors = c("#ffffff95", "#8bed80"),
      labels = c("potential area", "predicted area"),
      group = "Raster"
    )|>
    addLayersControl(
      overlayGroups = c(
        "Points",
        "Raster"
      ),
      position = "bottomleft",
      options = layersControlOptions(collapsed = FALSE),
    ) 
  
  # Add buffers if those are specified
  if(!is.na(buffer)){
    buffs <- st_as_sf(buffer)
    # Build a new map 
    map2 <- map1 |>
      addPolygons(
        data = buffs,
        weight = 1,
        fill = FALSE,
        opacity = 0.6,
        color = "#181b54",
        dashArray = "3",
        group = "Buffer"
      )|>
      addLayersControl(
        overlayGroups = c(
          "Buffer",
          "Points",
          "Raster"
        ),
        position = "bottomleft",
        options = layersControlOptions(collapsed = FALSE),
      ) 
    map <- map2
  }else{
    map <- map1 
  }
  return(map)
}

# ---- POINT SUMMARY FUNCTIONS ----
# The functions in this section are associated with the spatial metrics for buffer optimization (SMBO)
# analyses, which seek to determine whether an optimal buffer size for matching genetic and geographic
# coverage can be estimated using summaries of the geographic coordinates of each dataset. The functions
# were authored by Dan Carver.

#' geo.generateSpatialObject -- Dan Carver
#' Build a spatial object from coordinate data
geo.generateSpatialObject <- function(data){
  # Clean 
  data1 <- data |>
    dplyr::filter(!is.na(lon))|>
    dplyr::filter(!is.na(lat))
  # Lat long data
  sp1 <- sf::st_as_sf(x = data1,
                      coords = c("lon","lat"),
                      crs = CRS("+proj=longlat +datum=WGS84 +no_defs +type=crs"),
                      remove = FALSE)
  # Projected data 
  sp1_proj <- sf::st_transform(x = sp1, crs ="+proj=moll")
  return(sp1_proj)
}

#' geo.calc.EOO
#' Calculate the extent of occurence (EOO) value. Units: square meters
geo.calc.EOO <- function(data){
  # Convert to a multipoint object
  bb <- st_convex_hull(st_union(data))
  # Export value 
  eooArea <- st_area(bb,)
  return(eooArea)
}

#' geo.calc.AOO
#' Calculate the area of occurence (EOO) value. Units: square meters
geo.calc.AOO <- function(data){
  bb <- st_convex_hull(st_union(data))
  # Create the gridded feature 
  allGrids <- st_make_grid(bb,
                           square = T, 
                           cellsize = 2000)
  # Test for intersetion with the conver hull 
  intersectionCells <- sf::st_intersects(x = allGrids,y = bb,
                                         sparse= FALSE) # the grid, covering bounding box
  # Filter the full polygon feature
  selectAreas <- allGrids[intersectionCells]
  # Determine the number of areas that have an observation 
  intersectionPoints <- sf::st_intersects(x = selectAreas, y = data, sparse= TRUE)
  # Count of all features with at least one observation 
  testIntersect <- function(area){
    if(length(area)>0){
      return(1)
    }else{
      0
    }
  }
  areasWithPoints <- purrr::map(.x = intersectionPoints, .f = testIntersect )|> 
    unlist()|> 
    sum()
  
  # Calculate A00 and return
  aoo <- areasWithPoints / nrow(intersectionPoints) *100
  return(aoo)
}

#' geo.calc.averageNearestNeighbor 
#' Calculate the average nearest neighbor metric. Units: meters
geo.calc.averageNearestNeighbor <- function(data){
  # Convert data back to lat lon 
  refSP <- data |>
    sf::st_transform(crs = 4326)
  # Remove duplicated locations 
  refSP_noDup <- refSP[!duplicated(refSP),]
  # Calculate average nearest neighbor
  aveNearNeighbor <- refSP_noDup |> 
    mutate(
      nb = st_dist_band(geometry),
      dists = st_nb_dists(geometry, nb,longlat = TRUE),
      avg_dist = purrr::map_dbl(dists, mean),
      group_average = mean(avg_dist, na.rm = TRUE)
    ) 
  return(aveNearNeighbor$group_average[1])
}

#' geo.calc.voronoiAreas
#' Calculate the area of Voronoi cells. Unitless
geo.calc.voronoiAreas <- function(data){
  # Produce the tesselations 
  tesselation <- deldir(data$lon, data$lat)
  # Creates the spatial representation on the areas 
  tiles <- tile.list(tesselation)
  # Declare a helper function for indexing out of the tiles object
  selectArea <- function(tile){tile[6]}
  # Go through each tile and determine the overall average area
  aveVoronoiArea <- purrr::map(.x = tiles, .f = selectArea ) |> unlist() |> mean(na.rm=TRUE)
  return(aveVoronoiArea)
}

#' geo.calc.stdDistance
#' Calculate the standard distance metric. Units: meters
geo.calc.stdDistance <- function(data){
  stdDist <- std_distance(geometry = data)
  return(stdDist)
}

#' geo.calc.stdDistanceEllipseArea
#' Calculate the standard distance of the ellipse area. Units: meters
geo.calc.stdDistanceEllipseArea <- function(data){
  # Determine the standard distance 
  stdDist <- std_distance(geometry = data)
  # Determine the mean center 
  meanCenter <- sfdep::center_mean(geometry = data)
  # Grab coords from the mean center
  meanCoords <- as.data.frame(sf::st_coordinates(meanCenter))
  # Produces a list of points within the ellipse 
  stdDistEllipse <- ellipse(x = meanCoords$X,y = meanCoords$Y , sx = stdDist, sy = stdDist)
  # Convert to a polygon and project 
  poly <- stdDistEllipse |>
    as.data.frame()|>
    sfheaders::sf_polygon(x = "x", y = "y")
  sf::st_crs(poly) <- crs(data)
  # Make valid and reproject 
  validPoly <- st_make_valid(poly)
  # Determine area of the polygon, and return
  stdEllipseArea <- st_area(poly)
  return(stdEllipseArea)
}

#' geo.calc.stdDeviationEllipseArea
#' Calculate the standard deviation of the ellipse area. Units: meters
geo.calc.stdDevationEllipseArea <- function(data){
  stdDevElli <- std_dev_ellipse(geometry = data)
  # Points for the ellipse 
  stdDevEllipse <- st_ellipse(geometry =stdDevElli,
                              sx = stdDevElli$sx,
                              sy = stdDevElli$sy,
                              rotation = -stdDevElli$theta)
  # Determine length of the perimeter  
  stdDevationEllipseArea <- st_length(stdDevEllipse)
  return(stdDevationEllipseArea)
}

# Wrapper function of the above point summary functions, which will calculate
# each point summary statistic for a given dataset, and return 
geo.calc.pointSummaries <- function(geoData){
  # Begin by converting data into spatial object
  geoSpat <- geo.generateSpatialObject(geoData)
  # Run through the different spatial metric calculations, and store results to a data.frame
  EOO <- geo.calc.EOO(geoSpat)
  AOO <- geo.calc.AOO(geoSpat)
  ANN <- geo.calc.averageNearestNeighbor(geoSpat)
  VOR <- geo.calc.voronoiAreas(geoSpat)
  StdDist <- geo.calc.stdDistance(geoSpat)
  StdDistEA <- geo.calc.stdDistanceEllipseArea(geoSpat)
  StdDevEA <- geo.calc.stdDevationEllipseArea(geoSpat)
  # Combine metrics into a data.frame, and return
  ptSummaryDF <- data.frame(EOO, AOO, ANN, VOR, StdDist, StdDistEA, StdDevEA)
  return(ptSummaryDF)
}
