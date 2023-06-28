# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FUNCTIONS FOR GEOGRAPHIC-GENETIC COVERAGE CALCULATIONS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script declares the functions used for the analysis of correlation between geographic
# and genetic coverage. The createBuffers and compareBuffArea functions were sourced from 
# the gap analysis code developed by Emily Beckman Bruns.

# Load adegenet library, since some functions use the pop() and nInd() accessors
library(adegenet)
# Load parallel library, for parallelized resampling functions
library(parallel)

# Create buffers around points, using specified projection
createBuffers <- function(df, radius, pt_proj, buff_proj, boundary){
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
compareBuffArea <- function(insitu, samp_vect, radius, pt_proj, buff_proj, boundary){
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
exSitu_Sample <- function(gen_mat, geo_coordPts, geo_buff, geo_ptProj, geo_buffProj, geo_boundary, numSamples){
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
  geo_rate <- compareBuffArea(geo_coordPts, rownames(samp), geo_buff, geo_ptProj, geo_buffProj, geo_boundary)
  names(geo_rate) <- "Geo. coverage"
  # Combine genetic and geographic coverage rates into a vector, and return
  cov_rates <- c(gen_rates, geo_rate)
  return(cov_rates)
}

# Wrapper of exSitu_Sample: iterates that function over the entire sample matrix
exSitu_Resample <- function(gen_obj, geo_coordPts, geo_buff, geo_ptProj, geo_buffProj, geo_boundary){
  # Create a matrix of wild individuals (those with population "wild") from genind object
  gen_mat <- gen_obj@tab[which(pop(gen_obj) == "wild"),]
  # Apply the exSitu_Sample function to all rows of the wild matrix
  # (except row 1, because we need at least 2 individuals to sample)
  # The resulting matrix needs to be transposed, in order to keep columns as different allele categories
  cov_matrix <- t(sapply(2:nrow(gen_mat), 
                        function(x) exSitu_Sample(gen_mat, geo_coordPts, geo_buff, 
                                                  geo_ptProj, geo_buffProj, geo_boundary, x)))
  # Return the matrix of representation values
  return(cov_matrix)
}

# Wrapper of exSitu_Resample: runs resampling in parallel over a specified cluster. Results
# are saved to a specified file path.
exSitu_Resample_Parallel <- function(gen.obj, geo_coordPts, geo_buff, geo_ptProj, geo_buffProj, geo_boundary, 
                                     cluster, reps, arrayFilepath="~/resamplingArray.Rdata"){
  # Run resampling in parallel, capturing results to an array
  resamplingArray <- 
    parSapply(cluster, 1:reps, function(a) exSitu_Resample(gen.obj = gen.obj), simplify = "array")
  # Save the resampling array object to disk, for later usage
  saveRDS(resamplingArray, file = arrayFilepath)
  # Return resampling array to global environment
  return(resamplingArray)
}

# From resampling array, calculate the mean minimum sample size to represent 95% of the Total wild diversity
resample_min95_mean <- function(resamplingArray){
  # resampling array[,1,]: returns the Total column values for each replicate (3rd array dimension)
  # apply(resamplingArray[,1,],1,mean): calculates the average across replicates for each row
  # which(apply(resamplingArray[,1,],1,mean) > 95): returns the rows with averages greater than 95
  # min(which(apply(resamplingArray[,1,],1,mean) > 95)): the lowest row with an average greater than 95
  meanValue <- min(which(apply(resamplingArray[,1,],1, mean, na.rm=TRUE) > 95))
  return(meanValue)
}

# From resampling array, calculate the standard deviation, at the mean 95% value
resample_min95_sd <- function(resamplingArray){
  # Determine the mean value for representing 95% of allelic diversity
  meanValue <- resample_min95_mean(resamplingArray)
  # Calculate the standard deviation, at that mean value, and return
  sdValue <- apply(resamplingArray[,1,],1,sd)[meanValue]
  return(sdValue)
}

# From resampling array, calculate the mean values (across replicates) for each allele frequency category
resample_meanValues <- function(resamplingArray){
  # Declare a matrix to receive average values
  meanValue_mat <- matrix(nrow=nrow(resamplingArray), ncol=ncol(resamplingArray))
  # For each column in the array, average results across replicates (3rd array dimension)
  meanValue_mat[,1] <- apply(resamplingArray[,1,], 1, mean, na.rm=TRUE)
  meanValue_mat[,2] <- apply(resamplingArray[,2,], 1, mean, na.rm=TRUE)
  meanValue_mat[,3] <- apply(resamplingArray[,3,], 1, mean, na.rm=TRUE)
  meanValue_mat[,4] <- apply(resamplingArray[,4,], 1, mean, na.rm=TRUE)
  meanValue_mat[,5] <- apply(resamplingArray[,5,], 1, mean, na.rm=TRUE)
  # Give names to meanValue_mat columns, and return
  colnames(meanValue_mat) <- c("Total","Very common","Common","Low frequency","Rare")
  return(meanValue_mat)
}

# GRAPHING FUNCTIONS ----

# From resampling array, plot the results of the resampling analysis and save to a PDF file
resample_singlePlot_PDF <- function(arrayPath, imagePath="~/", colors, xLeg, yLeg, minSampleLineDist, mainText){
  # Create two vectors for colors. This is to show points on the graph and in the legend clearly
  fullColors <- colors
  fadedColors <- c(colors[1], alpha(colors[2:5], 0.4))
  # Read in the resampling array, based on the array path argument
  resamplingArray <- readRDS(file=arrayPath)
  # Generate the average values (across replicates) for each allele frequency category 
  averageValueMat <- resample_meanValues(resamplingArray)
  # Generate the minimum sample size to represent 95% of allelic diversity (across replicates)
  min95_Value <- resample_min95_mean(resamplingArray)
  # Call pdf command, to save resampling plot to disk. 
  pdf(file = imagePath, width = 9, height = 7.5)
  # Use the matplot function to plot the matrix of average values, with specified settings
  matplot(averageValueMat, ylim=c(0,110), col=fadedColors, pch=16, ylab="Allelic Representation (%)")
  # Add title and x-axis labels to the graph
  title(main=mainText, line=0.5)
  mtext(text="Number of individuals", side=1, line=2.4)
  # Mark the 95% threshold line, as well as the 95% minimum sampling size
  abline(h=95, col="black", lty=3); abline(v=min95_Value, col="black")
  # Add text for the minimum sampling size line. Location based on min 95 value and function argument
  mtext(text=paste0("Minimum sampling size (95%) = ", min95_Value),
        side=1, line=-1.5, at=min95_Value-minSampleLineDist, cex=1)
  # Add legend
  legend(x=xLeg, y=yLeg, inset = 0.05,
         legend = c("Total","Very common","Common","Low frequency", "Rare"),
         col=fullColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty="n", y.intersp = 1)
  # Turn off plotting device
  dev.off()
}

# Plots the results of two different resampling analyses (usually R0 and R80), and saves to a PDF file
resample_doublePlot_PDF <- function(arrayPath1, arrayPath2, imagePath="~/",
                                    colors, xLeg, yLeg, minSampleLineDist, mainText1, mainText2){
  # Create two vectors for colors. This is to show points on the graph and in the legend clearly
  fullColors <- colors
  fadedColors <- c(colors[1], alpha(colors[2:5], 0.4))
  # Read in and process resampling arrays
  # %%% FIRST ARRAY
  resamplingArray1 <- readRDS(file=arrayPath1)
  averageValueMat1 <- resample_meanValues(resamplingArray1)
  min95_Value1 <- resample_min95_mean(resamplingArray1)
  # %%% SECOND ARRAY
  resamplingArray2 <- readRDS(file=arrayPath2)
  averageValueMat2 <- resample_meanValues(resamplingArray2)
  min95_Value2 <- resample_min95_mean(resamplingArray2)
  # Call pdf command, to save resampling plot to disk. 
  pdf(file = imagePath, width = 9, height = 7.5)
  # Set plotting window to stack 2 graphs vertically
  par(mfcol=c(2,1), oma=rep(0.1,4), mar=c(3,4,2,1))
  
  # %%% FIRST ARRAY
  # Use the matplot function to plot the matrix of average values, with specified settings
  matplot(averageValueMat1, ylim=c(0,110), col=fadedColors, pch=16, ylab="Allelic Representation (%)")
  # Add title and x-axis labels to the graph
  title(main=mainText1, line=0.5)
  mtext(text="Number of individuals", side=1, line=1.8)
  # Mark the 95% threshold line, as well as the 95% minimum sampling size
  abline(h=95, col="black", lty=3); abline(v=min95_Value1, col="black")
  # Add text for the minimum sampling size line. Location based on min 95 value and function argument
  mtext(text=paste0("Minimum sampling size (95%) = ", min95_Value1),
        side=1, line=-1.5, at=min95_Value1-minSampleLineDist, cex=1)
  # Add legend
  legend(x=xLeg, y=yLeg, inset = 0.05,
         legend = c("Total","Very common","Common","Low frequency", "Rare"),
         col=fullColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty="n", y.intersp = 1)
  
  # %%% SECOND ARRAY
  # Use the matplot function to plot the matrix of average values, with specified settings
  matplot(averageValueMat2, ylim=c(0,110), col=fadedColors, pch=16, ylab="Allelic Representation (%)")
  # Add title and x-axis labels to the graph
  title(main=mainText2, line=0.5)
  mtext(text="Number of individuals", side=1, line=1.8)
  # Mark the 95% threshold line, as well as the 95% minimum sampling size
  abline(h=95, col="black", lty=3); abline(v=min95_Value2, col="black")
  # Add text for the minimum sampling size line. Location based on min 95 value and function argument
  mtext(text=paste0("Minimum sampling size (95%) = ", min95_Value2),
        side=1, line=-1.5, at=min95_Value2-minSampleLineDist, cex=1)
  # Add legend
  legend(x=xLeg, y=yLeg, inset = 0.05,
         legend = c("Total","Very common","Common","Low frequency", "Rare"),
         col=fullColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty="n", y.intersp = 1)
  
  # Turn off plotting device
  dev.off()
}
