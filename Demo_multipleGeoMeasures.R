# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEOGRAPHIC-GENETIC CORRELATION DEMO: QUERCUS ACERIFOLIA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script demonstrating first draft approach for calculating the correlation 
# between genetic and geographic coverage. Uses a genind file from Quercus
# acerifolia (SNP loci, complete dataset), as well as a CSV containing sample 
# names, latitudes and longitudes, to iteratively resample wild points and measure
# coverage.

library(adegenet)
library(terra)
library(parallel)

# Read in relevant functions
GeoGenCorr.wd <- "/home/akoontz/Documents/GeoGenCorr/Code/"
setwd(GeoGenCorr.wd)
source("functions_GeoGenCoverage.R")

# ---- VARIABLES ----
# Declare a directory within which to store .Rdata objects of resampling arrays
resamplingDataDir <- paste0(GeoGenCorr.wd, "resamplingData/")
# Specify number of resampling replicates. 
num_reps <- 5
# ---- GEOGRAPHIC VARIABLES
# Specify geographic buffer size in meters (i.e. 1 km)
buffSize <- 1000    
# Define projections for points (WGS84) and for calculations (which use units in meters and equal area)
ptProj <- "+proj=longlat +datum=WGS84"
calcProj <- "+proj=eqearth +datum=WGS84"
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- 
  vect(file.path("~/Documents/GeoGenCorr/QUAC_demo/gis_layers/world_countries_10m/world_countries_10m.shp"))
# ---- PARALLELIZATION
# Set up relevant cores 
num_cores <- detectCores() - 8 
cl <- makeCluster(num_cores)
# Make sure libraries (adegenet + terra) are on cluster
clusterEvalQ(cl, library("adegenet"))
clusterEvalQ(cl, library("terra"))
clusterEvalQ(cl, library("parallel"))

# ---- READ IN GEOGRAPHIC AND GENETIC DATA ----
# ---- GEOGRAPHIC
# Read in wild occurrence points. This CSV has 3 columns: sample name,latitude, and longitude. 
# The sample names (and order) have to match 
# the sample names/order of the genind object (rownams of the genetic matrix) read in below.
wildPoints <- read.csv("~/Documents/GeoGenCorr/QUAC_demo/Quercus_acerifolia.csv", header = TRUE)

# ---- GENETIC
# Read in genind file: Optimized de novo assembly; R80, min-maf=0, 
# first SNP/locus, 2 populations (garden and wild), no Kessler individuals.
# Wild sample names/order must match those in the sample name column of the CSV (above)
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops_NoK/"
QUAC.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames of genind. For this analysis, we'll only utilize wild samples (i.e. those in pop "wild")
pop(QUAC.genind) <- 
  factor(read.table(paste0(genpop.filePath, "QUAC_popmap_GardenWild_NoK"), header=FALSE)[,2])

# ---- RESAMPLING ----
# Export the coordinate points data.frame and genind object to the cluster
# clusterExport(cl, varlist = c("wildPoints", "QUAC.genind"))
clusterExport(cl, varlist = c("wildPoints","QUAC.genind","num_reps","buffSize","ptProj","calcProj","world_poly_clip"))
clusterExport(cl, varlist = c("createBuffers", "compareBuffArea", "getAlleleCategories","calculateCoverage",
                              "exSituResample", "geo.gen.Resample", "geo.gen.Resample.Parallel"))
# Specify file path, for saving resampling array
arrayDir <- paste0(resamplingDataDir, "QUACdemoResampArr2.Rdata")
# Run resampling
print("%%% BEGIN RESAMPLING %%%")

# Below line generates error:
# "Error in evaluating the argument 'x' in selecting a method for function 't': NULL value passed as symbol address"
QUAC_demoArray_Par <- geo.gen.Resample.Parallel(gen_obj=QUAC.genind, geo_coordPts=wildPoints,
                                                geo_ptProj=ptProj, geo_buffProj=calcProj,
                                                geo_boundary=world_poly_clip, reps=5,
                                                arrayFilepath=arrayDir, cluster=cl)

# Working: Run resampling (single replicate)
exSituResample(QUAC.genind, wildPoints, geo_boundary=world_poly_clip)
# Working: Run resampling (5 replicates, not in parallel), and save array
QUAC_demoArray <- geo.gen.Resample(gen_obj = QUAC.genind, geo_coordPts = wildPoints,
                                   geo_buff = buffSize, geo_ptProj = ptProj, geo_buffProj = calcProj,
                                   geo_boundary = world_poly_clip, reps = 5)
saveRDS(QUAC_demoArray, file=paste0(resamplingDataDir, "QUACdemoResampArr.Rdata"))

# ---- PLOTTING ----
# Specify plot colors
plotColors <- c("red","red4","darkorange3","coral","purple", "darkblue")
plotColors[2:5] <- alpha(plotColors[2:5], 0.35)
plotColors_Sub <- plotColors[-(2:5)]

# Calculate minimum 95% sample size for genetic and geographic values
gen_min95Value <- gen_min95Mean(QUAC_demoArray) ; gen_min95Value
gen_min95SD(QUAC_demoArray)
geo_min95Value <- geo_min95Mean(QUAC_demoArray) ; geo_min95Value
geo_min95SD(QUAC_demoArray)
# Generate the average values (across replicates) for all proportions 
averageValueMat <- meanArrayValues(QUAC_demoArray)
averageValueMat_Sub <- averageValueMat[,-(2:5)]

# %%%% TOTAL ALLELIC AND GEOGRAPHIC COLORS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat_Sub, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab="Coverage (%)")
# Add title and x-axis labels to the graph
title(main="Quercus acerifolia: Geo-Gen Coverage (5 resampling replicates)", line=0.5)
mtext(text="Number of individuals (91 Total)", side=1, line=2.4)
# Mark the 95% threshold line, and the genetic/geographic points
abline(h=95, col="black", lty=3) 
abline(v=gen_min95Value, col="red")
abline(v=geo_min95Value, col="darkblue")
# Add text for the minimum sampling size lines
mtext(text=paste0("Gen 95% MSSE = ", gen_min95Value),
      side=1, line=-1.5, at=70, cex=1)
mtext(text=paste0("Geo 95% MSSE = ", geo_min95Value),
      side=1, line=-1.5, at=15, cex=1)
# Add legend
legend(x=65, y=80, inset = 0.05,
       legend = c("Genetic coverage (Total)", "Geographic coverage (1km buffer)"),
       col=plotColors_Sub, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty="n", y.intersp = 0.3)

# %%%% ALL ALLELE CATEGORIES
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat, ylim=c(0,100), col=plotColors, pch=16, ylab="Coverage (%)")
# Add title and x-axis labels to the graph
title(main="Quercus acerifolia: Geographic and Genetic Coverage", line=0.5)
mtext(text="Number of individuals", side=1, line=2.4)
# Mark the 95% threshold line, as well as the 95% minimum sampling size
abline(h=95, col="black", lty=3); abline(v=min95_Value, col="black")
# Add text for the minimum sampling size line. Location based on min 95 value and function argument
mtext(text=paste0("Minimum sampling size (95%) = ", min95_Value),
      side=1, line=-1.5, at=50, cex=1)
# Add legend
legend(x=70, y=80, inset = 0.05,
       legend = c("Total","Very common","Common","Low frequency", "Rare", "Geo. coverage"),
       col=plotColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty="n", y.intersp = 0.3)

# Close cores
stopCluster(cl)
