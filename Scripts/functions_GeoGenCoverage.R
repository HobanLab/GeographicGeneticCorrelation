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
library(ape)
library(terra)
library(parallel)

# ---- WORKER FUNCTIONS USED TO CALCULATE COVERAGES ----
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
  # Clip by boundary, so buffers don't extend into the water. Skip reprojection if boundary is already 
  # in the target projection -- avoids redundant reprojection when the caller (e.g. calculateCoverage) 
  # has already projected boundary once for the whole call, rather than on every buffer operation.
  if(!terra::same.crs(boundary, buffProj)){
    boundary <- terra::project(boundary, buffProj)
  }
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
geo.compareBuff <- function(totalWildPoints, sampVect, buffSize, ptProj, buffProj, boundary, 
                            geoTotalArea=NULL, parFlag=FALSE){
  # If running in parallel: world polygon shapefile needs to be 'unwrapped', 
  # after being exported to cluster
  if(parFlag==TRUE){
    boundary <- unwrap(boundary)
  }
  # Select "ex situ" coordinates by subsetting totalWildPoints data.frame, according to sampVect
  exSitu <- totalWildPoints[sort(match(sampVect, totalWildPoints[,1])),]
  # Create buffers around selected (exSitu) wild points 
  geo_exSitu <- createBuffers(exSitu, buffSize, ptProj, buffProj, boundary)
  # Calculate the area under the exSitu buffer. The 1,000,000 value converts values to km²
  geo_exSituArea <- expanse(geo_exSitu)/1000000
  # If a precomputed total buffer area wasn't provided, compute it here (original behavior). The total 
  # area (all wild points) is constant for a given buffSize -- see geo.totalBuffArea(), which precomputes 
  # this once per buffer size, rather than it being redundantly recomputed on every call here.
  if(is.null(geoTotalArea)){
    geo_total <- createBuffers(totalWildPoints, buffSize, ptProj, buffProj, boundary)
    geoTotalArea <- expanse(geo_total)/1000000
  }
  # Calculate the proportion of the ex situ buffer areas to the total buffer area (percent geographic coverage)
  geo_Coverage <- (geo_exSituArea/geoTotalArea)*100
  return(geo_Coverage)
}

# WORKER FUNCTION: Analogous to geo.compareBuff, but instead of comparing buffered area around a random
# sample of points to the total buffered area, compares it instead to the total area under SDM, which
# is passed to this function as a raster argument (model).The sampVect argument represents vector 
# of sample names, which is used to subset the totalWildPoints data.frame to create a separate 
# "ex situ" spatial object. Then, the createBuffers function is used to place buffers around sampled
# points, and the proportion of the total area covered is calculated. Authored by Dan Carver
geo.compareBuffSDM <- function(totalWildPoints, sampVect, buffSize, model, ptProj, buffProj, 
                               boundary, geoTotalArea_SDM=NULL, parFlag=FALSE){
  # If running in parallel, unwrap spatial features 
  if(parFlag==TRUE){
    boundary <- unwrap(boundary)
    model <- unwrap(model)
  }
  # Select "ex situ" coordinates by subsetting totalWildPoints data.frame, according to sampVect
  exSitu <- totalWildPoints[sort(match(sampVect, totalWildPoints[,1])),]
  # Create buffers around selected (exSitu) wild points 
  geo_exSitu <- createBuffers(exSitu, buffSize, ptProj, buffProj, boundary) |>
    # Reproject to match crs of smd object 
    terra::project(model)
  # Crop the (potentially very large, e.g. disaggregated global-extent) raster down to just the buffer's 
  # local extent BEFORE doing any pixel-level operations. A single buffer covers a tiny fraction of a 
  # global raster; classify/rasterize/multiply/cellSize on the full raster for every sample was the actual 
  # dominant per-call cost -- not redundant repetition (already fixed), but legitimate work sized to the 
  # whole raster instead of the small area actually needed. geoTotalArea_SDM (the denominator) is computed 
  # separately, from the uncropped raster, so this does not change what's being measured -- only how much 
  # irrelevant area gets processed to compute the ex situ buffer's area.
  modelCrop <- terra::crop(model, geo_exSitu)
  # Generate a mask of the (cropped) model layer by converting all 0 values to NA  
  m <- c(0, 0, NA)
  mask <- modelCrop |>
    classify(m)
  # Rasterize 
  buffRast <- geo_exSitu |>
    terra::rasterize(mask)
  # Apply mask
  buffMask <- buffRast * mask
  # Calculate the area under the exSitu raster in km²
  geo_exSituArea <- terra::cellSize(buffMask, mask= TRUE, unit = "km") |>
    terra::values()|>
    sum(na.rm = TRUE)
  # If a precomputed total SDM area wasn't provided, compute it here (original behavior) -- from the FULL, 
  # uncropped raster, since the total area denominator must reflect the entire SDM extent. This total area 
  # depends only on the model raster -- not on buffSize or sample -- see geo.totalSDMArea(), which 
  # precomputes this once for the whole run, rather than it being redundantly recomputed on every call here
  # (this was especially costly given the size of the SDM raster).
  if(is.null(geoTotalArea_SDM)){
    geoTotalArea_SDM <- terra::cellSize(model |> classify(m), mask= TRUE, unit = "km") |>
      terra::values()|>
      sum(na.rm = TRUE)
  }
  # Calculate the proportion of the ex situ buffer areas to the total buffer area (percent geographic coverage)
  geo_Coverage <- (geo_exSituArea/geoTotalArea_SDM)*100
  return(geo_Coverage)
}

# DIAGNOSTIC HELPER: prints a timestamped system memory snapshot (via `free -h`) to stdout. Used to bracket 
# major steps in geo.gen.Resample.Par, to identify exactly which step memory usage climbs during. Read-only; 
# does not affect any computation.
memCheck <- function(label){
  cat(paste0('\n@@@ MEM CHECK [', label, '] ', Sys.time(), ' @@@\n'))
  cat(paste(system('free -h', intern=TRUE), collapse='\n'))
  cat('\n')
}

# WORKER FUNCTION: writes a timestamped progress+memory checkpoint to a per-worker log file, every 
# `every` samples. Worker stdout is NOT visible in the master's nohup.out during parSapply execution (PSOCK 
# workers' output is not forwarded by default), so this is the only way to see what a worker was doing --  
# and what memory looked like -- at the moment it stops responding. logDir must be a location visible to 
# all workers (e.g. shared with the master, same as SDMrast_files); defaults to tempdir().
workerProgressLog <- function(sampleIndex, every=1, logDir=tempdir()){
  if(sampleIndex %% every != 0) return(invisible(NULL))
  logFile <- file.path(logDir, paste0('worker_', Sys.getpid(), '_progress.log'))
  memInfo <- paste(system('free -h', intern=TRUE), collapse=' | ')
  cat(paste0(Sys.time(), ' | sample ', sampleIndex, ' | ', memInfo, '\n'), file=logFile, append=TRUE)
}

# WORKER FUNCTION: This function compares the resolution of a raster argument to the geographic buffer size 
# being used, and returns the disaggregation scale factor needed (1 = no resampling needed). Does NOT return
# a raster copy itself -- see geo.buildSDMscaleFactors()/geo.getSDMrast(), which use this to avoid building 
# or broadcasting a resampled raster per buffer size (most buffer sizes typically need no resampling, or 
# share the same scale factor).
# Authored by Dan Carver; modified to return only a scale factor.
geo.checkSDMres <- function(buffSize, raster, parFlag=FALSE){
  # If running in parallel, unwrap the raster argument
  if(parFlag==TRUE){
    raster <- unwrap(raster)
  }
  # Convert buffer from meters to degrees (assuming 1 m = 0.000012726903908907691 degrees)
  buffSizeDegree <- buffSize * 0.000012726903908907691
  # Extract the raster resolution
  sdmDegree <- terra::res(raster)[1]
  # If geographic buffer less than SDM resolution, give a warning and flag that resampling is needed
  if(buffSizeDegree < sdmDegree){
    scaleFactor <- ceiling(sdmDegree / buffSizeDegree)
    warning(paste0('SDM provided has a resolution larger than geographic buffer size (',
                   buffSize, ' m). SDM will be resampled to a smaller resolution (scale factor: ', 
                   scaleFactor, ').'))
  } else {
    scaleFactor <- 1
  }
  # Return the scale factor (not a raster)
  return(scaleFactor)
}

# WORKER FUNCTION: Given a vector of geographic buffer sizes and a single SDM raster, determine the 
# disaggregation scale factor needed for EACH buffer size. Unlike the earlier geo.buildSDMrastList, this 
# does NOT build or return any resampled raster -- it only returns the small numeric scale-factor vector. 
# Building the (potentially much larger, disaggregated) raster is deferred to each worker individually, 
# via geo.getSDMrast() below -- this avoids the master building a large disaggregated raster and 
# broadcasting a full copy of it to every cluster worker via clusterExport, which was very costly for a 
# global-extent SDM raster.
geo.buildSDMscaleFactors <- function(geoBuff, raster, parFlag=FALSE){
  sapply(geoBuff, function(x) geo.checkSDMres(buffSize=x, raster=raster, parFlag=parFlag))
}

# WORKER FUNCTION: Returns a (locally cached, per-worker) version of the SDM raster, resampled to the given 
# scale factor. If a pre-written file path is available for this scale factor (SDMrast_files -- written 
# ONCE by the master via geo.buildSDMrastFiles, for scale factors > 1), the raster is read from disk, 
# which is disk-backed/lightweight rather than a full in-memory disaggregation. This avoids every worker 
# independently disaggregating the same raster in memory at once (a "thundering herd" at the start of 
# resampling, when all workers begin their first task simultaneously). If no file is available (scale 
# factor 1, or SDMrast_files not provided), falls back to unwrapping/disaggregating SDMrast_orig directly. 
# Either way, the result is cached in a persistent cache environment so all subsequent calls -- across 
# every sample and every replicate processed by that worker -- reuse the cached version. The cache 
# environment (.sdmRastCache) is created lazily, in whichever session's global environment this function 
# is running in, the first time it's called -- this makes the function self-contained: it does NOT rely on 
# the cache being separately shipped to or initialized on each worker (a function's free variables, when 
# the function is defined at top level, resolve in whichever R session's global environment it is actually 
# running in -- so a cache created only on the master would not be visible to workers; each worker needs, 
# and now gets, its own).
geo.getSDMrast <- function(scaleFactor, SDMrast_orig, SDMrast_files=NULL){
  if(!exists('.sdmRastCache', envir=.GlobalEnv)){
    assign('.sdmRastCache', new.env(), envir=.GlobalEnv)
  }
  cache <- get('.sdmRastCache', envir=.GlobalEnv)
  key <- as.character(scaleFactor)
  if(!exists(key, envir=cache)){
    if(scaleFactor > 1 && !is.null(SDMrast_files) && key %in% names(SDMrast_files)){
      r <- terra::rast(SDMrast_files[[key]])
    } else {
      r <- unwrap(SDMrast_orig)
      if(scaleFactor > 1){
        r <- terra::disagg(r, scaleFactor)
      }
    }
    assign(key, r, envir=cache)
  }
  return(get(key, envir=cache))
}

# MASTER-SIDE FUNCTION: For each DISTINCT scale factor > 1 required across a set of buffer sizes, 
# disaggregate the SDM raster ONCE and write the result to a GeoTIFF file (in a location accessible to 
# both the master and all cluster workers, since a local PSOCK cluster's workers share the master's 
# filesystem). Returns a named vector of file paths (keyed by scale factor, as a string), which is small 
# and cheap to clusterExport -- unlike broadcasting the disaggregated raster object(s) themselves. Workers 
# then read these files via geo.getSDMrast(), rather than each independently disaggregating in memory.
geo.buildSDMrastFiles <- function(scaleFactors, raster, parFlag=FALSE, outDir=tempdir()){
  uniqueFactors <- unique(scaleFactors)
  uniqueFactors <- uniqueFactors[uniqueFactors > 1]
  if(length(uniqueFactors) == 0){
    return(NULL)
  }
  files <- sapply(uniqueFactors, function(sf){
    r <- if(parFlag==TRUE) unwrap(raster) else raster
    r <- terra::disagg(r, sf)
    fp <- file.path(outDir, paste0('SDMrast_scaleFactor', sf, '.tif'))
    terra::writeRaster(r, fp, overwrite=TRUE)
    fp
  })
  names(files) <- as.character(uniqueFactors)
  return(files)
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
  # Make sure ecoregions are in same projection as buffers. Skip reprojection if already in the target 
  # projection -- avoids redundantly reprojecting the (large, global) ecoregions layer on every call.
  if(!terra::same.crs(ecoRegion, buffProj)){
    ecoProj <- terra::project(ecoRegion, buffProj)
  } else {
    ecoProj <- ecoRegion
  }
  # Crop the (large, global) ecoregion polygon down to just the buffers' extent BEFORE intersecting. 
  # A buffer around a handful of points covers a small fraction of a global ecoregion layer's extent; 
  # intersecting against the full layer for every sample is the same class of oversized-per-call-work 
  # problem already fixed for the SDM raster (geo.compareBuffSDM crop) -- this had not yet been applied 
  # here. Cropping first can only remove features that couldn't intersect anyway, so this does not 
  # change the result.
  ecoProjCrop <- terra::crop(ecoProj, buffers)
  # Intersect buffers with (cropped) ecoregions, and return
  ecoBuffJoin <- terra::intersect(buffers, ecoProjCrop)
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

# WORKER FUNCTION: Analogous to eco.totalEcoregionCount, but for the geographic buffer approach: computes 
# the TOTAL buffered area (all wild points) for a given buffer size, ONCE. This value is constant across 
# all samples/reps for a given buffer size, so precomputing it (rather than recomputing it inside every 
# geo.compareBuff call) avoids substantial redundant buffer/aggregate/crop computation.
geo.totalBuffArea <- function(totalWildPoints, buffSize, ptProj, buffProj, boundary, parFlag=FALSE){
  # If running in parallel: world polygon shapefile needs to be 'unwrapped', after being exported to cluster
  if(parFlag==TRUE){
    boundary <- unwrap(boundary)
  }
  # Create buffer around all (total) occurrences, and calculate area in km²
  geo_total <- createBuffers(totalWildPoints, buffSize, ptProj, buffProj, boundary)
  geo_totalArea <- expanse(geo_total)/1000000
  return(geo_totalArea)
}

# WORKER FUNCTION: Analogous to geo.totalBuffArea, but for the SDM-based approach: computes the total 
# masked SDM area, ONCE. This value does not depend on buffSize or sample -- only on the SDM raster 
# (model) -- so it only needs to be computed a single time per run, rather than recomputed inside every 
# geo.compareBuffSDM call (which is especially costly for a large raster).
geo.totalSDMArea <- function(model, parFlag=FALSE){
  # If running in parallel, unwrap the raster argument
  if(parFlag==TRUE){
    model <- unwrap(model)
  }
  # Generate a mask of the model layer by converting all 0 values to NA
  m <- c(0, 0, NA)
  mask <- model |> classify(m)
  # Calculate the total area under the mask, in km²
  geo_totalArea <- terra::cellSize(mask, mask=TRUE, unit="km") |>
    terra::values() |>
    sum(na.rm=TRUE)
  return(geo_totalArea)
}

# WORKER FUNCTION: Function for reporting representation rates, given (1) a genetic matrix,
# which includes ALL samples, and (2) a vector of sample names, which represents the samples
# of interest to calculate genetic coverage for. The function first builds a vector of allele
# frequencies from the entire genetic matrix, before performing the following steps:
# 1. The length of matches between garden and wild alleles is calculated (numerator). 
# 2. The complete number of wild alleles of that category (denominator) is calculated. 
# 3. From these 2 values, a percentage is calculated. 
# This function returns the numerators, denominators, and the proportion (representation rates) 
# in a matrix. Individual samples are allowed.
gen.getAlleleCategories <- function(genMat, sampNames){
  # Calculate allele frequency vector from complete genetic matrix, removing missing alleles
  freqVector <- colSums(genMat, na.rm = TRUE)/(nrow(genMat)*2)*100
  freqVector <- freqVector[which(freqVector != 0)]
  # Based on sample names vector, extract the relevant rows from complete genetic matrix
  sampMat <- genMat[which(rownames(genMat) %in% sampNames),]
  # Remove absent alleles. Conditional is to accommodate sample sizes of 1
  if(length(sampNames) == 1){
    # Remove any missing alleles from the sample "matrix" (vector)
    sampMat <- sampMat[!is.na(sampMat) & sampMat != 0]
    # sampMat <- sampMat[which(sampMat != 0)]
    # Capture names of present alleles in sample "matrix" (vector)
    sampAlleleNames <- names(sampMat)
  } else {
    # Remove any missing alleles (those with colSums of 0) from the sample matrix
    sampMat <- sampMat[, which(colSums(sampMat, na.rm = TRUE) != 0)]
    # Capture names of present alleles in sample matrix
    sampAlleleNames <- colnames(sampMat)
  }
  
  # CALCULATE GENETIC COVERAGES
  # Determine how many Total alleles in the sample matrix are found in the frequency vector 
  exSitu_allAlleles <- length(which(names(freqVector) %in% sampAlleleNames))
  total_allAlleles <- length(freqVector)
  allPercentage <- (exSitu_allAlleles/total_allAlleles)*100
  # Very common alleles (greater than 10%)
  exSitu_vComAlleles <- length(which(names(which(freqVector > 10)) %in% sampAlleleNames))
  total_vComAlleles <- length(which(freqVector > 10))
  vComPercentage <- (exSitu_vComAlleles/total_vComAlleles)*100
  # Common alleles (greater than 5%)
  exSitu_comAlleles <- length(which(names(which(freqVector > 5)) %in% sampAlleleNames))
  total_comAlleles <- length(which(freqVector > 5))
  comPercentage <- (exSitu_comAlleles/total_comAlleles)*100
  # Low frequency alleles (between 1% and 10%)
  exSitu_lowFrAlleles <- length(which(names(which(freqVector < 10 & freqVector > 1)) %in% sampAlleleNames))
  total_lowFrAlleles <- length(which(freqVector < 10 & freqVector > 1))
  lowFrPercentage <- (exSitu_lowFrAlleles/total_lowFrAlleles)*100
  # Rare alleles (less than 1%)
  exSitu_rareAlleles <- length(which(names(which(freqVector < 1)) %in% sampAlleleNames))
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

# WORKER FUNCTION: given a genind object, builds a matrix of Euclidean distances 
# between every individual. Exploratory function testing alternative approaches 
# for measuring genetic coverage
gen.buildDistMat <- function(genObj){
  # Convert genind object to data.frame, ignoring any population designations
  df <- genind2df(genObj, usepop=FALSE)
  # Calculate distance matrix, removing loci with at least one missing value, and return
  distMat <- dist.gene(df, method='pairwise', pairwise.deletion=TRUE)
  return(distMat)
}

# WORKER FUNCTION: given a genetic distance matrix and a vector of sample names, 
# calculate the sum total of the genetic distances between all pairs of individuals 
# (denominator) and the sum of the genetic distances strictly between sampled 
# individuals (numerator), and returns a coverage metric
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

# ---- BUILDING THE RESAMPLING ARRAY ----
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
calculateCoverage <- function(genMat, genType=c('CV','DI','EN','AN','EE','HE'), 
                              genDistMat=NA, genCHGeno=NA, geoFlag=TRUE, coordPts, 
                              geoBuff, SDMrast=NA, SDMrast_scaleFactor=NULL, SDMrast_files=NULL, geoTotalArea=NULL, 
                              geoTotalArea_SDM=NULL, ptProj='+proj=longlat +datum=WGS84', 
                              buffProj='+proj=eqearth +datum=WGS84', boundary, 
                              ecoFlag=FALSE, ecoBuff, ecoTotalCount, ecoRegions, 
                              ecoLayer=c('US','NA','GL'), parFlag=FALSE, numSamples){
  # DRAW RANDOM SAMPLES
  # From matrix of individuals (geMat), select a random set (rows). This is the set of individuals that will
  # be used for all downstream coverage calculations within this function. The sampNames object
  # is simply the vector of sample names.
  sampNames <- sample(rownames(genMat), numSamples)
  # Log progress (every 20 samples) to a per-worker file -- worker stdout isn't visible in nohup.out during 
  # parSapply execution, so this is the only way to see how far a worker got, and what memory looked like, 
  # if it stops responding.
  if(parFlag==TRUE){
    workerProgressLog(sampleIndex=numSamples)
  }
  
  # GENETIC PROCESSING
  # Conditionals below capture different gen. metrics; first is "standard" (default) allelic coverage
  if(genType=='CV'){
    # Genetic coverage: calculate sample's allelic representation
    genRates <- gen.getAlleleCategories(genMat, sampNames)
    # Subset matrix returned by getAlleleCategories to just 3rd column (representation rates), and return
    genRates <- genRates[,3]
  }
  # Calculate total genetic distance
  if(genType=='DI'){
    # Ensure the genetic distance matrix has been provided
    if(class(genDistMat)=='logical'){
      stop('Genetic distance matrix must be provided when specifying DI genType!')
    }
    # Pass distance matrix and sample name vector to function calculating  
    # proportion of total pairwise genetic distances represented in sample
    genDistCov <- gen.calcGenDistCov(distMat=genDistMat, sampVect=sampNames)
    # Append the resulting coverage value to the allelic coverages
    genRates <- c(genRates[,3], genDistCov)
    names(genRates)[6] <- 'GenDist'
  }
  # Calculate CoreHunter metrics
  if (genType %in% c('EN', 'AN', 'EE', 'HE')){
    # Ensure the CoreHunter genotype object has been provided
    if(length(class(genCHGeno))==1){
      stop('CoreHunter genotype object must be provided when specifying EN, AN, EE, or HE genType!')
    }
    # Calculate genRate value using corehunter::evaluateCore function
    genRates <- evaluateCore(sampNames, genCHGeno, objective = objective(type = genType))
    names(genRates) <- paste0('CH_',genType)
  }
  
  # If running in parallel, unwrap spatial objects ONCE per calculateCoverage call (i.e. once per sample),
  # rather than once per buffer size inside each geo.compareBuff/geo.compareBuffSDM/eco.compareBuff call below.
  # Previously, a single calculateCoverage call could re-unwrap boundary/ecoRegions/SDM rasters from their 
  # packed form dozens of times (once per buffer size, per geo/eco function) -- across many samples, buffer 
  # sizes, and reps, this adds up to hundreds of thousands of redundant unwrap() calls cluster-wide, which was 
  # causing memory to climb quickly. Downstream calls receive parFlag=FALSE, since the objects passed to them
  # below are now already live (unwrapped) -- this matches how the non-parallel path already calls them.
  if(parFlag==TRUE){
    if(geoFlag==TRUE || ecoFlag==TRUE){
      boundary <- unwrap(boundary)
    }
    if(ecoFlag==TRUE){
      ecoRegions <- unwrap(ecoRegions)
    }
    parFlag <- FALSE
  }
  # Note: SDMrast is intentionally NOT unwrapped here. It stays as the single, small, wrapped original 
  # raster; geo.getSDMrast() (called below, in the SDM branch) unwraps and -- if needed -- disaggregates 
  # it, caching the result in a persistent, worker-local cache so the (potentially expensive) 
  # disaggregation only happens once per worker for the life of the run, rather than once per call, and 
  # so the disaggregated raster is never itself broadcast over the network.
  # Reproject boundary/ecoRegions ONCE per calculateCoverage call, rather than inside every downstream 
  # buffer/intersect operation (createBuffers/eco.intersectBuff now skip their own internal reprojection 
  # when the object they receive is already in the target projection -- see the terra::same.crs() checks 
  # there). For the GLOBAL ecoregions layer in particular, reprojection is expensive; this was previously 
  # being repeated on every one of the ~41 buffer-size iterations, for every sample, for every rep. This 
  # step runs regardless of parFlag, so the non-parallel path benefits too.
  if(geoFlag==TRUE || ecoFlag==TRUE){
    boundary <- terra::project(boundary, buffProj)
  }
  if(ecoFlag==TRUE){
    ecoRegions <- terra::project(ecoRegions, buffProj)
  }
  
  # GEOGRAPHIC PROCESSING
  if(geoFlag==TRUE){
    # Check that sample names in genetic matrix match the column of sample names in the coordinate dataframe
    if(!identical(rownames(genMat), coordPts[,1])){
      stop('Error: Sample names between the genetic matrix and the first 
         column of the coordinate point data.frame do not match.')
    }
    # Geographic coverage: for each buffer size, calculate sample's geograhpic representation, by 
    # passing all points (coordPts) and the random subset of points (sampNames) to the 
    # geo.compareBuff worker function, which will calculate the proportion of area covered in the 
    # random sample. 
    geoRates <- 
      lapply(seq_along(geoBuff), function(i) geo.compareBuff(totalWildPoints=coordPts, sampVect=sampNames,
                                                             buffSize=geoBuff[i], ptProj=ptProj, buffProj=buffProj, 
                                                             boundary=boundary, 
                                                             geoTotalArea=if(is.null(geoTotalArea)) NULL else geoTotalArea[i],
                                                             parFlag=parFlag))
    # If no rasterized SDM is provided, calculate geographic coverage using just the buffer approach (default)
    if(class(SDMrast)=='logical'){
      names(geoRates) <- paste0(rep('Geo_Buff_',), geoBuff/1000, 'km')
    } else {
      # If rasterized SDM provided, calculate geo. coverage using SDM approach, and append coverages. 
      # SDMrast is the single, original (undisaggregated) raster; SDMrast_scaleFactor (same length/order 
      # as geoBuff) gives the disaggregation scale factor needed for each buffer size. geo.getSDMrast() 
      # builds (and worker-locally caches) the resampled raster version for each scale factor as needed, 
      # rather than indexing into a pre-built list that was broadcast to every worker up front.
      geoRates_SDM <- 
        mapply(function(b,sf) geo.compareBuffSDM(totalWildPoints=coordPts, sampVect=sampNames,
                                                 buffSize=b, model=geo.getSDMrast(scaleFactor=sf, SDMrast_orig=SDMrast, SDMrast_files=SDMrast_files), 
                                                 ptProj=ptProj, buffProj=buffProj, boundary=boundary, 
                                                 geoTotalArea_SDM=geoTotalArea_SDM,
                                                 parFlag=FALSE), b=geoBuff, sf=SDMrast_scaleFactor)
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
    # points (coordPts) and the random subset of points (sampNames) to the eco.compareBuff worker function, 
    # which will calculate the proportion of ecoregions covered in the random sample. The ecoTotalCount argument
    # provides the denominator (ecoregions covered by all points) used for coverage calculations.
    ecoRates <- 
      mapply(function(b,t) eco.compareBuff(totalWildPoints=coordPts, sampVect=sampNames,
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
exSituResample <- function(genMat, genType=c('CV','DI','EN','AN','EE','HE'),
                           genDistMat=NA, genCHGeno=NA, geoFlag=TRUE, coordPts, geoBuff=50000, SDMrast=NA, 
                           SDMrast_scaleFactor=NULL, SDMrast_files=NULL, ptProj='+proj=longlat +datum=WGS84', buffProj='+proj=eqearth +datum=WGS84', 
                           boundary, ecoFlag=FALSE, ecoBuff=50000, ecoTotalCount, ecoRegions, ecoLayer='US', parFlag){
  # Apply the calculateCoverage function to all rows of the wild matrix.
  # The resulting matrix needs to be transposed, in order to keep columns as different coverage categories
  cov_matrix <- 
    t(sapply(1:nrow(genMat), 
             function(x) calculateCoverage(genMat=genMat, genType=genType, genDistMat=genDistMat, 
                                           genCHGeno=genCHGeno, geoFlag=geoFlag, coordPts=coordPts, 
                                           geoBuff=geoBuff, SDMrast=SDMrast, SDMrast_scaleFactor=SDMrast_scaleFactor, SDMrast_files=SDMrast_files, ptProj=ptProj, 
                                           buffProj=buffProj, boundary=boundary, ecoFlag=ecoFlag, 
                                           ecoBuff=ecoBuff, ecoTotalCount=ecoTotalCount, ecoRegions=ecoRegions, 
                                           ecoLayer=ecoLayer, parFlag=FALSE, numSamples=x), simplify = 'array'))
  # Return the matrix of coverage values
  return(cov_matrix)
}

# WRAPPER FUNCTION: iterates calculateCoverage over the entire matrix of samples, in parallel
exSituResample.Par <- function(genMat, genType=c('CV','DI','EN','AN','EE','HE'),
                               genDistMat=NA, genCHGeno=NA, geoFlag=TRUE, coordPts, geoBuff=50000, SDMrast=NA, 
                               SDMrast_scaleFactor=NULL, SDMrast_files=NULL, geoTotalArea=NULL, geoTotalArea_SDM=NULL, 
                               ptProj='+proj=longlat +datum=WGS84', buffProj='+proj=eqearth +datum=WGS84', 
                               boundary, ecoFlag=FALSE, ecoBuff=50000, ecoTotalCount, ecoRegions, ecoLayer='US', 
                               parFlag=TRUE, cluster){
  # Use parSapply (load balanced) to iterate calculate coverage for each sample size. Transpose the resulting matrix.
  cov_matrix <-
    t(parSapply(cluster, 1:nrow(genMat),
                function(x) calculateCoverage(genMat=genMat, genType=genType, genDistMat=genDistMat, genCHGeno=genCHGeno,
                                              geoFlag=geoFlag, coordPts=coordPts, geoBuff=geoBuff, SDMrast=SDMrast,
                                              SDMrast_scaleFactor=SDMrast_scaleFactor, SDMrast_files=SDMrast_files, geoTotalArea=geoTotalArea, 
                                              geoTotalArea_SDM=geoTotalArea_SDM, ptProj=ptProj, buffProj=buffProj, boundary=boundary, ecoFlag=ecoFlag,
                                              ecoBuff=ecoBuff, ecoTotalCount=ecoTotalCount, ecoRegions=ecoRegions,
                                              ecoLayer=ecoLayer, parFlag=parFlag, numSamples=x), simplify = 'array'))
  # Return the matrix of coverage values
  return(cov_matrix)
}

# WRAPPER FUNCTION: iterates exSituResample, which will generate an array of values from a single genind object.
# Checks for arguments are also performed; if ecological coverage is being calculated, then the total number
# of ecoregions will . This function doesn't run in parallel, so it's  primarily used
# for testing/demonstration purposes
geo.gen.Resample <- function(genObj, genType='CV', geoFlag=TRUE, coordPts, geoBuff=50000, SDMrast=NA,
                             ptProj='+proj=longlat +datum=WGS84', buffProj='+proj=eqearth +datum=WGS84',
                             boundary, ecoFlag=FALSE, ecoBuff=50000, ecoRegions,
                             ecoLayer=c('US', 'NA', 'GL'), reps=5){
  # Check that genType argument is an allowed value; if not, notify user
  if(!(genType %in% c('CV','DI','EN', 'AN', 'EE', 'HE'))){
    stop('The genType argument can only be equal to CV, DI, EN, AN, EE, or HE!')
  }
  # Extract the genetic matrix from the genind object
  genMat <- genObj@tab
  # If genetic distances are being calculated, build matrix of genetic distances to pass down to lower functions
  if(genType=='DI'){
    cat(paste0('\n', '-- genType set to DI: will calculate genetic distance coverage --'))
    genDistMat <- gen.buildDistMat(genObj=genObj)
  } else {
    # Otherwise, set genDistMat as NA
    genDistMat <- NA
  } 
  # If CoreHunter metrics are being calculated, build CoreHunter genotype object to pass down to lower functions
  if(genType %in% c('EN', 'AN', 'EE', 'HE')){
    cat(paste0('\n', '-- genType set to ', genType, ': will use CoreHunter metrics --'))
    genCHGeno <- genotypes(genMat, format="biparental")
  } else {
    # Otherwise, set genCHGeno as NA
    genCHGeno <- NA
  } 
  # Initialize SDMrast_scaleFactor as NULL; only set if geoFlag=TRUE and an SDM raster is provided (below)
  SDMrast_scaleFactor <- NULL
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
      # Check that geographic buffer size is greater than SDM raster resolution, and fix if not. Rather 
      # than building one (potentially large) resampled raster per buffer size up front, only the 
      # per-buffer-size scale factor is computed here; geo.getSDMrast() builds (and locally caches) the 
      # actual resampled raster for a given scale factor lazily, the first time it's needed.
      SDMrast_scaleFactor <- geo.buildSDMscaleFactors(geoBuff=geoBuff, raster=SDMrast, parFlag=FALSE)
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
    sapply(1:reps, function(x) exSituResample(genMat=genMat, genType=genType, genDistMat=genDistMat, 
                                              genCHGeno=genCHGeno, geoFlag=geoFlag, coordPts=coordPts, 
                                              geoBuff=geoBuff, SDMrast=SDMrast, SDMrast_scaleFactor=SDMrast_scaleFactor, ptProj=ptProj,
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
geo.gen.Resample.Par <- function(genObj, genType='CV', genCHGeno=NA, geoFlag=TRUE, coordPts, 
                                 geoBuff=50000, SDMrast=NA, ptProj='+proj=longlat +datum=WGS84', 
                                 buffProj='+proj=eqearth +datum=WGS84', boundary, ecoFlag=FALSE, 
                                 ecoBuff=50000, ecoRegions, ecoLayer=c('US','NA','GL'), reps=5,
                                 arrayFilepath='~/resamplingArray.Rdata', cluster){
  # Check that genType argument is an allowed value; if not, notify user
  if(!(genType %in% c('CV','DI','EN', 'AN', 'EE', 'HE'))){
    stop('The genType argument can only be equal to CV, DI, EN, AN, EE, or HE!')
  }
  # Extract the genetic matrix from the genind object, and export it to the cluster
  genMat <- genObj@tab
  clusterExport(cl = cluster, varlist = 'genMat', envir = environment())
  # If genetic distance flag is set to TRUE, build the matrix of genetic distances to pass down to lower functions
  if(genType=='DI'){
    cat(paste0('\n', '-- genType set to DI: will calculate genetic distance coverage --'))
    genDistMat <- gen.buildDistMat(genObj=genObj)
  } else {
    # Otherwise, set genDistMat as NA
    genDistMat <- NA
  }
  # If CoreHunter metrics are being calculated, make sure CoreHunter genotype object has been specified
  if(genType %in% c('EN', 'AN', 'EE', 'HE')){
    # Ensure the genetic distance matrix has been provided
    if(length(class(genCHGeno))==1){
      stop('CoreHunter Genotype object must be provided when specifying genType EE, AN, EE, or HE.')
    } else{
      cat(paste0('\n', '-- genType set to ', genType, ': expecting CoreHunter genotype to be exported to cluster! --'))
    }
  } 
  # Initialize SDMrast_scaleFactor/SDMrast_files as NULL; only set if geoFlag=TRUE and an SDM raster is 
  # provided (below)
  SDMrast_scaleFactor <- NULL
  SDMrast_files <- NULL
  # Initialize precomputed geographic total-area denominators as NULL; only set if geoFlag=TRUE (below).
  # These mirror ecoTotalCount (below): the "total" area used as the coverage denominator does not depend 
  # on the sample or replicate -- only on buffer size (or, for the SDM approach, not even on buffer size) 
  # -- so it's computed ONCE here, rather than being redundantly recomputed inside every geo.compareBuff / 
  # geo.compareBuffSDM call (previously happening on the order of hundreds of thousands of times).
  geoTotalArea <- NULL
  geoTotalArea_SDM <- NULL
  memCheck('function start')
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
    # Unwrap and project boundary ONCE here (locally), for use in the total-area precomputation below. This 
    # local, live copy (boundary_live) is only used here on the master for precomputation -- it is never 
    # exported to the cluster (only the resulting numeric areas are).
    boundary_live <- unwrap(boundary)
    boundary_live <- terra::project(boundary_live, buffProj)
    memCheck('boundary unwrapped+projected')
    # CALCULATE TOTAL GEOGRAPHIC (BUFFER) COVERAGE: the total buffered area across all wild points, for 
    # each buffer size. This is constant across samples/reps for a given buffer size.
    cat(paste0('\n', '--- CALCULATING TOTAL GEOGRAPHIC BUFFER COVERAGE... ---'))
    geoTotalArea <- sapply(geoBuff, function(x) geo.totalBuffArea(totalWildPoints=coordPts, buffSize=x,
                                                                  ptProj=ptProj, buffProj=buffProj,
                                                                  boundary=boundary_live, parFlag=FALSE))
    memCheck('after geoTotalArea precompute')
    # If SDM is provided (meaning it's not NA, or class logical):
    if(!class(SDMrast)=='logical'){
      # Print out message stating what coverages are being calculated, and how many buffer sizes
      cat(paste0('\n', '- SDM provided: will calculate geographic coverage (SDM) -'))
      cat(paste0('\n', '--- Number of buffer sizes (Geo, SDM): ', length(geoBuff), ' ---'))
      # Determine the disaggregation scale factor needed for each buffer size.
      SDMrast_scaleFactor <- geo.buildSDMscaleFactors(geoBuff=geoBuff, raster=SDMrast, parFlag=TRUE)
      cat(paste0('\n', '--- Distinct SDM scale factors needed: ', length(unique(SDMrast_scaleFactor)), 
                 ' (of ', length(geoBuff), ' buffer sizes) ---'))
      memCheck('after SDM scale factor check')
      # For any scale factor > 1, disaggregate ONCE here (on the master) and write the result to disk, 
      # rather than broadcasting the disaggregated raster object to every worker (very costly -- ~64GB 
      # spike in earlier testing) OR letting every worker independently disaggregate in memory the moment 
      # resampling starts (a "thundering herd" of simultaneous heavy operations -- the cause of the most 
      # recent failure). Workers read the small file path and load the raster lazily/disk-backed via 
      # geo.getSDMrast(), instead of holding a full in-memory copy built by each worker separately.
      SDMrast_files <- geo.buildSDMrastFiles(scaleFactors=SDMrast_scaleFactor, raster=SDMrast, parFlag=TRUE)
      memCheck('after SDM raster file(s) written')
      # CALCULATE TOTAL SDM COVERAGE: the total masked SDM area does not depend on buffer size or sample, 
      # and does not depend on the disaggregation scale factor either (disaggregating doesn't change total 
      # area) -- so it's computed ONCE here, directly from the original (undisaggregated) raster.
      geoTotalArea_SDM <- geo.totalSDMArea(model=SDMrast, parFlag=TRUE)
      memCheck('after geoTotalArea_SDM precompute')
      # Export the small original raster, the scale-factor vector, the (small) file path lookup, and the 
      # precomputed total areas to the cluster. SDMrast here is still just the single, original, small 
      # (undisaggregated) raster -- the disaggregated raster itself is never sent over the network.
      clusterExport(cl=cluster, varlist=c('SDMrast','SDMrast_scaleFactor','SDMrast_files','geoTotalArea','geoTotalArea_SDM'), envir=environment())
      memCheck('after clusterExport (SDMrast/geoTotalArea)')
    } else {
      # Export the precomputed total buffer areas to the cluster
      clusterExport(cl=cluster, varlist='geoTotalArea', envir=environment())
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
    # Unwrap and project boundary/ecoRegions ONCE here (locally; reusing boundary_live if it was already 
    # built above), for use in the total-ecoregion-count precomputation below. This avoids the unwrap + 
    # reproject that used to happen inside every single one of the 41 eco.totalEcoregionCount calls below 
    # -- which, for the (large, global) ecoregions layer, was very costly, and was a major contributor to 
    # memory climbing steeply even before the parallel resampling loop started.
    if(!exists('boundary_live')){
      boundary_live <- unwrap(boundary)
      boundary_live <- terra::project(boundary_live, buffProj)
    }
    ecoRegions_live <- unwrap(ecoRegions)
    ecoRegions_live <- terra::project(ecoRegions_live, buffProj)
    memCheck('ecoRegions unwrapped+projected')
    # CALCULATE TOTAL ECOLOGICAL COVERAGE: calculate the number of ecoregions found under all samples for all 
    # buffer sizes, and pass this down to lower level functions, to optimize processing
    cat(paste0('\n', '--- CALCULATING TOTAL ECOREGION COVERAGE... ---'))
    ecoTotalCount <- lapply(ecoBuff, 
                            function(x) eco.totalEcoregionCount(totalWildPoints=coordPts, buffSize=x,
                                                                ptProj=ptProj, buffProj=buffProj, 
                                                                ecoRegion=ecoRegions_live, layerType=ecoLayer,
                                                                boundary=boundary_live, parFlag=FALSE))
    memCheck('after ecoTotalCount precompute')
  }
  # Print starting time
  startTime <- Sys.time() 
  memCheck('immediately before resampling loop')
  cat(paste0('\n', '%%% RESAMPLING START: ', startTime, '\n'))
  # Run resampling for all replicates, using sapply and lambda function
  resamplingArray <- 
    sapply(1:reps, function(x) exSituResample.Par(genMat=genMat, genType=genType, genDistMat=genDistMat, 
                                                  genCHGeno=genCHGeno, geoFlag=geoFlag, 
                                                  coordPts=coordPts, geoBuff=geoBuff, 
                                                  SDMrast=SDMrast, SDMrast_scaleFactor=SDMrast_scaleFactor, SDMrast_files=SDMrast_files, 
                                                  geoTotalArea=geoTotalArea, geoTotalArea_SDM=geoTotalArea_SDM, ptProj=ptProj, 
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
nrmse.func <-  function(obs_var, pred_var, norm_type='mean', digits=NA){
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
    # If digits argument is specified, round result to specified number of digits, and return
    if( is.numeric(digits) == TRUE){
      rmse <- round(rmse, digits)
    }
    return(rmse)
  } else {
    # Normalize RMSE, based on type argument
    if (norm_type == 'sd') nrmse <- rmse/sd(obs_var)
    if (norm_type == 'mean') nrmse <- rmse/mean(obs_var)
    if (norm_type == 'maxmin') nrmse <- rmse/(max(obs_var) - min(obs_var))
    if (norm_type == 'iq') nrmse <- rmse/(quantile(obs_var, 0.75) - quantile(obs_var, 0.25))
    # If digits argument is specified, round result to specified number of digits, and return
    if( is.numeric(digits) == TRUE){
      nrmse <- round(nrmse, digits)
    }
    return(nrmse)
  }
}

# Function for building a matrix of normalized root mean square error (NRMSE) values based on 
# a data.frame of resampling values and the genetic coverage type used for the predictor variable.
# Given these arguments, a matrix is generated which has a NRMSE value for each set of geographic/ecological
# coverage values compared to the genetic coverage metric. This command expects the resampling data.frame
# to have certain row names and column names, and makes those checks
buildNRMSEmatrix <- function(resampDF, genCovType=c('CV', 'GD'), sdmFlag=FALSE, NRMSEdigits=2){
  # Match argument for type of response variable (genetic coverage approach) to use
  genCovType <- match.arg(genCovType)
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
    NRMSEmat <- matrix(NA, nrow=length(buffSizes), ncol=3)
    colnames(NRMSEmat) <- c('Geo_Buff','Geo_SDM','Eco_Buff')
  } else {
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
    NRMSEvalue <- nrmse.func(resampDF[,i], pred_var = predVar, digits=NRMSEdigits)
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

# Function which, given a matrix of NRMSE values for each given coverage type, will 
# return the optimal buffer sizes (those with the minimum NRMSE values) as a vector
getOptBuffers <- function(nrmseMat){
  # Create vector capturing optimal buffer sizes; name it according to the NRMSE matrix column names
  optBuffs <- vector(length=ncol(nrmseMat)) ; names(optBuffs) <- colnames(nrmseMat)
  # Loop through the columns of the NRMSE matrix. For each column, find the minimum value, and return
  # that columns name (which is the buffer size value) as a numeric
  for(i in 1:length(optBuffs)){
    optBuffs[i] <- as.numeric(unlist(strsplit(names(which.min(nrmseMat[,i])),'km')))
  }
  # Return the list of optimal buffer sizes
  return(optBuffs)
}

# Wrapper function which, given a filepath to a resampling array, will return the optimal 
# buffer sizes according to the data in that array
extractOptBuffs <- function(arrayDir='~/resamp.Rdata', genCovType=c('CV', 'GD'), sdmFlag=FALSE, NRMSEdigits=2){
  # Match genetic coverage type argument
  genCovType <- match.arg(genCovType)
  # Read in the array, then convert it into a data.frame
  array <- readRDS(arrayDir)
  df <- resample.array2dataframe(array)
  # From the data.frame, calculate a matrix of NRMSE values for each coverage type. Pass arguments
  # along for the genetic coverage type (allelic/haplotypic coverage vs. genetic distance; whether
  # or not SDM is included)
  nrmseMat <- buildNRMSEmatrix(resampDF=df, genCovType=genCovType, sdmFlag=sdmFlag, NRMSEdigits=NRMSEdigits)
  # From the NRMSE matrix, extract the optimal buffer sizes, and return
  optBuffs <- getOptBuffers(nrmseMat)
  return(optBuffs)
}

# Function for building a matrix of average coverage values for optimal buffer sizes
extractOptCovs <- function(arrayDir='~/resamp.Rdata', genCovType=c('CV', 'GD'), sdmFlag=FALSE, NRMSEdigits=2){
  # Match genetic coverage type argument
  genCovType <- match.arg(genCovType) 
  # Read in the array, then convert it into a data.frame
  array <- readRDS(arrayDir)
  df <- resample.array2dataframe(array)
  # From the data.frame, calculate a matrix of NRMSE values for each coverage type. Pass arguments
  # along for the genetic coverage type (allelic/haplotypic coverage vs. genetic distance; whether
  # or not SDM is included)
  nrmseMat <- buildNRMSEmatrix(resampDF=df, genCovType=genCovType, sdmFlag=sdmFlag, NRMSEdigits=NRMSEdigits)
  # From the NRMSE matrix, extract the optimal buffer sizes
  optBuffs <- getOptBuffers(nrmseMat)
  # From array, generate average value matrix (across resampling replicates)
  averageValMat <- meanArrayValues(array)
  # Subset the average value matrix to only optimal buffer size values
  if(genCovType=='CV'){
    optCovMat <- averageValMat[,paste0(strsplit(names(optBuffs),split = '_CV'),'_', optBuffs,'km'),]
  } else {
    optCovMat <- averageValMat[,paste0(strsplit(names(optBuffs),split = '_GD'),'_', optBuffs,'km'),]
  }
  # Add, as the first column, the average genetic coverage value (Total), and return
  optCovMat <- cbind(averageValMat[,'Total'], optCovMat); colnames(optCovMat)[[1]] <- 'Gen_Total'
  return(optCovMat)
}

# Function for building a matrix of correlation values (NRMSE, Pearson correlation, or Spearman correlation)
# based on an array path (resampling array to read in) and the genetic coverage type used for the predictor variable.
# Given these arguments, a matrix is generated which has a NRMSE value for each set of geographic/ecological
# coverage values compared to the genetic coverage metric. This command expects the resampling data.frame
# to have certain row names and column names, and makes those checks
buildCorrelationMat <- function(arrayDir='~/resamp.Rdata', genCovType=c('CV', 'GD'), 
                                corMetric=c('NRMSE', 'corSp', 'corPe'), sdmFlag=FALSE, NRMSEdigits=2){
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
      corValue <- nrmse.func(resampDF[,i], pred = predVar, digits=NRMSEdigits)
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
#' Build a spatial object from coordinate data (columns lat, lon in decimal degrees, WGS84).
#' Points are projected to Mollweide, an equal-area projection, so that the area based metrics
#' (EOO, AOO, ELA) are measured in an equal-area projection as recommended by the IUCN Red List
#' guidelines. Distances (ANN, STD, ELP) are handled separately within each function.
geo.generateSpatialObject <- function(data){
  # Clean 
  data1 <- data |>
    dplyr::filter(!is.na(lon))|>
    dplyr::filter(!is.na(lat))
  # Lat long data
  sp1 <- sf::st_as_sf(x = data1,
                      coords = c("lon","lat"),
                      crs = 4326,
                      remove = FALSE)
  # Projected data 
  sp1_proj <- sf::st_transform(x = sp1, crs = "+proj=moll")
  return(sp1_proj)
}

#' geo.calc.EOO
#' Calculate the extent of occurrence (EOO): the area of the minimum convex polygon (convex hull)
#' around all locations, following the IUCN Red List guidelines. The hull is drawn and measured
#' in the equal-area projection of the input (see geo.generateSpatialObject). Units: km2
geo.calc.EOO <- function(data){
  # Minimum convex polygon around all locations
  mcp <- st_convex_hull(st_union(data))
  # Area in km2
  eooArea <- as.numeric(st_area(mcp)) / 1e6
  return(eooArea)
}

#' geo.calc.AOO
#' Calculate the area of occupancy (AOO). A grid of 2 km x 2 km cells (the IUCN Red List reference
#' scale) is placed over the minimum convex polygon, in the equal-area projection of the input, and
#' the cells containing at least one location are counted. Two values are returned:
#'   AOO_km2: number of occupied cells x 4 km2, the IUCN Red List AOO
#'   AOO_pct: occupied cells as a percentage of all cells intersecting the minimum convex polygon
#'            (the value reported in earlier drafts; a measure of occupancy relative to EOO)
#' Note that the cell counts depend on where the grid is placed; the grid is anchored at the lower
#' left corner of the bounding box of the locations.
geo.calc.AOO <- function(data, cellsize = 2000){
  # Minimum convex polygon around all locations
  mcp <- st_convex_hull(st_union(data))
  # Grid of 2 km cells covering the bounding box of the minimum convex polygon
  allGrids <- st_make_grid(mcp, square = TRUE, cellsize = cellsize)
  # Keep the cells intersecting the minimum convex polygon
  mcpCells <- allGrids[sf::st_intersects(x = allGrids, y = mcp, sparse = FALSE)]
  # Count the cells with at least one location
  occupied <- sf::st_intersects(x = mcpCells, y = data, sparse = TRUE) |>
    lengths()
  occupiedCells <- sum(occupied > 0)
  # Area of occupancy (km2) and percent of cells occupied
  cellArea <- (cellsize / 1000)^2
  output <- data.frame(AOO_km2 = occupiedCells * cellArea,
                       AOO_pct = occupiedCells / length(mcpCells) * 100)
  # Specifying only the return of the km^2 metric, to keep downstream analyses clean
  return(output)
}

#' geo.calc.averageNearestNeighbor 
#' Calculate the average nearest neighbor distance (Clark and Evans 1954): for each location, the
#' distance to the closest other location, averaged over all locations. Distances are geodesic
#' (great circle distances on a sphere, via sf::st_distance with s2) so they do not depend on the projection. Individuals 
#' sharing the same coordinates are collapsed to a single location first; otherwise site level 
#' datasets (e.g. AMTH, PICO, VILA) would have nearest neighbor distances of zero. Units: km
geo.calc.averageNearestNeighbor <- function(data){
  # Unique locations, in lat lon for geodesic distances
  uniqueLoc <- data[!duplicated(sf::st_coordinates(data)), ] |>
    sf::st_transform(crs = 4326)
  # Pairwise distance matrix (m); a location is not its own neighbor
  distMat <- sf::st_distance(uniqueLoc) |>
    units::drop_units()
  diag(distMat) <- NA
  # Distance from each location to its nearest neighbor
  nearestDist <- apply(distMat, 1, min, na.rm = TRUE)
  # Average, in km
  ANN <- mean(nearestDist) / 1000
  return(ANN)
}

#' geo.calc.voronoiAreas
#' Calculate the evenness of the spatial sampling from a Voronoi tessellation. Each unique location
#' is assigned the area closer to it than to any other location; the tessellation is built in the
#' equal-area projection of the input and clipped to the minimum convex polygon (EOO). The metric
#' is the coefficient of variation (sd / mean) of the cell areas: evenly spaced locations produce
#' cells of similar size (values near 0), clustered sampling produces a few large and many small
#' cells (larger values). Note that the mean cell area itself is always EOO / number of locations,
#' regardless of the arrangement of the points, so it is not reported. Unitless
geo.calc.voronoiAreas <- function(data){
  # Unique locations (duplicated coordinates produce empty cells)
  uniqueLoc <- data[!duplicated(sf::st_coordinates(data)), ]
  # Minimum convex polygon, used to bound the outer cells
  mcp <- sf::st_convex_hull(sf::st_union(uniqueLoc))
  # Voronoi tessellation, clipped to the minimum convex polygon
  cells <- sf::st_voronoi(sf::st_union(uniqueLoc)) |>
    sf::st_collection_extract("POLYGON") |>
    sf::st_intersection(mcp)
  # Cell areas (km2) and their coefficient of variation
  cellAreas <- as.numeric(sf::st_area(cells)) / 1e6
  vorCV <- sd(cellAreas) / mean(cellAreas)
  return(vorCV)
}

#' geo.calc.stdDistance
#' Calculate the standard distance: the root mean square distance of the unique locations from
#' their mean center, a measure of the overall dispersion of the sampling (Bachi 1963). Calculated
#' in the equal-area projection of the input. Units: km
geo.calc.stdDistance <- function(data){
  # Unique locations
  uniqueLoc <- data[!duplicated(sf::st_coordinates(data)), ]
  # Standard distance (m) to km
  stdDist <- sfdep::std_distance(geometry = uniqueLoc) / 1000
  return(stdDist)
}

#' geo.calc.ellipseElongation
#' Calculate the elongation of the standard deviational ellipse (Yuill 1971) of the unique
#' locations: the ratio of the major to the minor axis (one standard deviation along each axis).
#' A value of 1 indicates an isotropic (circular) arrangement of the locations, larger values an
#' increasingly linear arrangement. Note that the area and perimeter of the ellipse, and the area
#' of the standard distance circle (pi * STD^2), are all size measures which rank the datasets in 
#' the same order as STD; the elongation is the shape information the ellipse adds. Unitless
geo.calc.ellipseElongation <- function(data){
  # Unique locations
  uniqueLoc <- data[!duplicated(sf::st_coordinates(data)), ]
  # Standard deviational ellipse (sx, sy: axis lengths; theta: rotation)
  sde <- sfdep::std_dev_ellipse(geometry = uniqueLoc)
  # Ratio of the major to the minor axis
  elongation <- max(sde$sx, sde$sy) / min(sde$sx, sde$sy)
  return(elongation)
}

#' #' geo.calc.stdDistanceEllipseArea <-- NO LONGER USED
#' #' Calculate the standard distance of the ellipse area. Units: meters
#' geo.calc.stdDistanceEllipseArea <- function(data){
#'   # Determine the standard distance 
#'   stdDist <- std_distance(geometry = data)
#'   # Determine the mean center 
#'   meanCenter <- sfdep::center_mean(geometry = data)
#'   # Grab coords from the mean center
#'   meanCoords <- as.data.frame(sf::st_coordinates(meanCenter))
#'   # Produces a list of points within the ellipse 
#'   stdDistEllipse <- ellipse(x = meanCoords$X,y = meanCoords$Y , sx = stdDist, sy = stdDist)
#'   # Convert to a polygon and project 
#'   poly <- stdDistEllipse |>
#'     as.data.frame()|>
#'     sfheaders::sf_polygon(x = "x", y = "y")
#'   sf::st_crs(poly) <- crs(data)
#'   # Make valid and reproject 
#'   validPoly <- st_make_valid(poly)
#'   # Determine area of the polygon and return
#'   stdEllipseArea <- st_area(poly)
#'   return(stdEllipseArea)
#' }
#' 
#' #' geo.calc.stdDeviationEllipseArea <-- NO LONGER USED
#' #' Calculate the standard deviation of the ellipse perimeter. Units: meters
#' geo.calc.stdDevationEllipseArea <- function(data){
#'   stdDevElli <- std_dev_ellipse(geometry = data)
#'   # Points for the ellipse 
#'   stdDevEllipse <- st_ellipse(geometry =stdDevElli,
#'                               sx = stdDevElli$sx,
#'                               sy = stdDevElli$sy,
#'                               rotation = -stdDevElli$theta)
#'   # Determine length of the perimeter and return
#'   stdDevationEllipseArea <- st_length(stdDevEllipse)
#'   return(stdDevationEllipseArea)
#' }

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
  STD <- geo.calc.stdDistance(geoSpat)
  ELG <- geo.calc.ellipseElongation(geoSpat)
  # Combine metrics into a data.frame, and return
  ptSummaryDF <- data.frame(EOO = EOO, AOO = AOO$AOO_km2, AOO_pct = AOO$AOO_pct, 
                            ANN = ANN, VOR = VOR, STD = STD, ELG = ELG)
  return(ptSummaryDF)
}
