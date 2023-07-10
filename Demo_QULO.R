# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEOGRAPHIC-GENETIC CORRELATION DEMO: QUERCUS LOBATA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script demonstrating first draft approach for calculating the correlation 
# between genetic and geographic coverage. Uses data files from Gugger et al. 2021
# to pull in genetic data (as a table) and geographic coordinates (in a CSV)
# to conduct geographic-genetic correlation analyses.

library(adegenet)
library(terra)
library(parallel)
library(RColorBrewer)
library(scales)

# Read in relevant functions
GeoGenCorr.wd <- "/home/akoontz/Documents/GeoGenCorr/"
setwd(GeoGenCorr.wd)
source("Code/functions_GeoGenCoverage.R")

# ---- VARIABLES ----
# Declare a directory within which to store .Rdata objects of resampling arrays
resamplingDataDir <- paste0(GeoGenCorr.wd, "Code/resamplingData/")
# Specify number of resampling replicates. 
num_reps <- 5
# ---- GEOGRAPHIC VARIABLES
# Specify geographic buffer size in meters (i.e. 1 km, or 50 km)
# buffSize <- 1000
buffSize <- 50000
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- 
  vect(file.path(paste0(GeoGenCorr.wd, "GIS_shpFiles/world_countries_10m/world_countries_10m.shp")))
# This shapefile is by default a "non-exportable" object, which means it must be processed before it can be
# exported to the cluster (for geographic calculations). The terra::wrap function is used to do this.
world_poly_clip_W <- wrap(world_poly_clip)

# ---- PARALLELIZATION
# Set up relevant cores 
num_cores <- detectCores() - 8 
cl <- makeCluster(num_cores)
# Make sure libraries (adegenet + terra) are on cluster
clusterEvalQ(cl, library("adegenet"))
clusterEvalQ(cl, library("terra"))
clusterEvalQ(cl, library("parallel"))

# ---- READ IN GEOGRAPHIC AND GENETIC DATA ----
# Specify filepath for QUAC geographic and genetic data
QULO.filePath <- paste0(GeoGenCorr.wd, "Datasets/QULO/")

# ---- GEOGRAPHIC
# Read in wild occurrence points. This CSV has 3 columns: sample name, latitude, and longitude. 
# The sample names (and order) have to match the sample names/order of the genind object 
# (rownams of the genetic matrix) read in below.
wildPoints <- read.csv(paste0(QULO.filePath, "Quercus_lobata.csv"), header=TRUE)

# ---- GENETIC
# Processing input file
QULO.tab <- read.table(paste0(QULO.filePath, "SNPs80.forR"), header = TRUE)
# Make sample names the row names
rownames(QULO.tab) <- QULO.tab[,1] ; QULO.tab <- QULO.tab[,-1]
# Convert to genind. ncode based off similar practice for other geninds (values of 1, 2, and 3 generate identical results)
QULO.genind <- df2genind(QULO.tab, ncode = 3, ploidy = 2)

# ---- RESAMPLING ----
# Export the coordinate points data.frame and genind object to the cluster
clusterExport(cl, varlist = c("wildPoints","QULO.genind","num_reps","buffSize","world_poly_clip_W"))
clusterExport(cl, varlist = c("createBuffers", "compareBuffArea", "getAlleleCategories","calculateCoverage",
                              "exSituResample", "geo.gen.Resample.Parallel"))
# Specify file path, for saving resampling array
arrayDir <- paste0(QULO.filePath, "resamplingData/QULO_50km_5r_resampArr.Rdata")
# Run resampling (in parallel)
QULO_demoArray_Par <- geo.gen.Resample.Parallel(gen_obj=QULO.genind, geo_coordPts=wildPoints,
                                                geo_buff=buffSize,
                                                geo_boundary=world_poly_clip_W, reps=5,
                                                arrayFilepath=arrayDir, cluster=cl)
# Close cores
stopCluster(cl)

# ---- CORRELATION AND PLOTTTING ----
# Specify plot colors
plotColors <- c("red","red4","darkorange3","coral","purple", "darkblue")
plotColors[2:5] <- alpha(plotColors[2:5], 0.35)
plotColors_Sub <- plotColors[-(2:5)]

# Calculate minimum 95% sample size for genetic and geographic values
gen_min95Value <- gen_min95Mean(QULO_demoArray_Par) ; gen_min95Value
gen_min95SD(QULO_demoArray_Par)
geo_min95Value <- geo_min95Mean(QULO_demoArray_Par) ; geo_min95Value
geo_min95SD(QULO_demoArray_Par)
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic and geographic proportions
averageValueMat <- meanArrayValues(QULO_demoArray_Par)
# Generate linear model between both coverage values
QULO_model <- lm (Total ~ Geo, data=averageValueMat)
QULO_model_summary <- summary(QULO_model) ; QULO_model_summary
# Pull R-squared and p-value estimates from model
QULO_model_rSquared <- QULO_model_summary$adj.r.squared
QULO_model_pValue <- QULO_model_summary$coefficients[2, 4]

# %%%% GEOGRAPHIC-GENETIC CORRELATION
plot(averageValueMat$Geo, averageValueMat$Total, pch=20, main="Q. lobata: Geographic by genetic coverage",
     xlab="Geographic coverage (%)", ylab="Genetic coverage (%)")
mtext(text="436 Individuals; 50 km buffer; 5 replicates", side=3, line=0.3)
mylabel = bquote(italic(R)^2 == .(format(QULO_model_rSquared, digits = 3)))
text(x = 15, y = 85, labels = mylabel)

# %%%% TOTAL ALLELIC AND GEOGRAPHIC COVERAGE
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab="Coverage (%)")
# Add title and x-axis labels to the graph
title(main="Quercus lobata: Geo-Gen Coverage", line=1.5)
mtext(text="436 Individuals; 50 km buffer; 5 replicates", side=3, line=0.3)
mtext(text="Number of individuals", side=1, line=2.4)
# Mark the 95% threshold line, and the genetic/geographic points
abline(h=95, col="black", lty=3) 
abline(v=gen_min95Value, col="red")
abline(v=geo_min95Value, col="darkblue")
# Add text for the minimum sampling size lines
mtext(text=paste0("Gen 95% MSSE = ", gen_min95Value),
      side=1, line=-1.5, at=95, cex=1)
mtext(text=paste0("Geo 95% MSSE = ", geo_min95Value),
      side=1, line=-1.5, at=200, cex=1)
# Add legend
legend(x=205, y=60, inset = 0.05,
       legend = c("Genetic coverage (Total)", "Geographic coverage (50 km buffer)"),
       col=plotColors_Sub, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty="n", y.intersp = 0.5)
