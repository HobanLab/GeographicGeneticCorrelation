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
# num_reps <- 5
num_reps <- 3
# ---- GEOGRAPHIC VARIABLES
# Specify geographic buffer size in meters 
buffSize <- 1000
# buffSize <- 50000
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
# Make sure libraries (adegenet, terra, and parallel) are on cluster
clusterEvalQ(cl, library("adegenet"))
clusterEvalQ(cl, library("terra"))
clusterEvalQ(cl, library("parallel"))

# %%%% INDIVIDUAL-LEVEL GEOGRAPHIC COORDINATES %%%% ----
# In this analysis, we utilize a CSV file of lat/longs that specify the location of each individual
# ---- READ IN GEOGRAPHIC AND GENETIC DATA ----
# Specify filepath for QUAC geographic and genetic data
QUAC.filePath <- paste0(GeoGenCorr.wd, "Datasets/QUAC/")

# ---- GEOGRAPHIC
# Read in wild occurrence points. This CSV has 3 columns: sample name, latitude, and longitude. 
# The sample names (and order) have to match the sample names/order of the genind object 
# (rownams of the genetic matrix) read in below.
wildPoints <- read.csv(paste0(QUAC.filePath, "QUAC_coord_ind.csv"), header=TRUE)

# ---- GENETIC
# Read in genind file: Optimized de novo assembly; R80, min-maf=0, 
# first SNP/locus, 2 populations (garden and wild), no Kessler individuals.
# Wild sample names/order must match those in the sample name column of the CSV (above)
QUAC.genind <- read.genepop(paste0(QUAC.filePath,"QUAC_DNFA_populations_R80_NOMAF_1SNP_2Pops_NoK.gen"))
# Correct popNames of genind. For this analysis, we'll only utilize wild samples (i.e. those in pop "wild")
pop(QUAC.genind) <- 
  factor(read.table(paste0(QUAC.filePath, "QUAC_popmap_GardenWild_NoK"), header=FALSE)[,2])

# ---- RESAMPLING ----
# Export the coordinate points data.frame and genind object to the cluster
clusterExport(cl, varlist = c("wildPoints","QUAC.genind","num_reps","buffSize","world_poly_clip_W"))
clusterExport(cl, varlist = c("createBuffers", "compareBuffArea", "getAlleleCategories","calculateCoverage",
                              "exSituResample", "geo.gen.Resample.Parallel"))
# Specify file path, for saving resampling array
# arrayDir <- paste0(QUAC.filePath, "resamplingData/QUAC_1kmIND_5r_resampArr.Rdata")
arrayDir <- paste0(QUAC.filePath, "resamplingData/QUAC_1kmIND_3r_resampArr.Rdata")
# Run resampling
QUAC_demoArray_Par <- geo.gen.Resample.Parallel(gen_obj=QUAC.genind, geo_coordPts=wildPoints,
                                                geo_buff=buffSize,
                                                geo_boundary=world_poly_clip_W, reps=num_reps,
                                                arrayFilepath=arrayDir, cluster=cl)
# Close cores
stopCluster(cl)

# ---- CORRELATION AND PLOTTTING ----
# Specify plot colors
plotColors <- c("red","red4","darkorange3","coral","purple", "darkblue")
plotColors[2:5] <- alpha(plotColors[2:5], 0.35)
plotColors_Sub <- plotColors[-(2:5)]

# Calculate minimum 95% sample size for genetic and geographic values
gen_min95Value <- gen_min95Mean(QUAC_demoArray_Par) ; gen_min95Value
gen_min95SD(QUAC_demoArray_Par)
geo_min95Value <- geo_min95Mean(QUAC_demoArray_Par) ; geo_min95Value
geo_min95SD(QUAC_demoArray_Par)
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic and geographic proportions
averageValueMat <- meanArrayValues(QUAC_demoArray_Par)
# Generate linear model between both coverage values
QUAC_model <- lm (Total ~ Geo, data=averageValueMat)
QUAC_model_summary <- summary(QUAC_model) ; QUAC_model_summary
# Pull R-squared and p-value estimates from model
QUAC_model_rSquared <- QUAC_model_summary$adj.r.squared; QUAC_model_rSquared
QUAC_model_pValue <- QUAC_model_summary$coefficients[2, 4]

# %%%% GEOGRAPHIC-GENETIC CORRELATION
plot(averageValueMat$Geo, averageValueMat$Total, pch=20, main="Q. acerifolia: Geographic by genetic coverage",
     xlab="Geographic coverage (%)", ylab="Genetic coverage (%)")
mtext(text="91 Individuals; 1 km buffer (individuals); 5 replicates", side=3, line=0.3)
mylabel = bquote(italic(R)^2 == .(format(QUAC_model_rSquared, digits = 3)))
text(x = 45, y = 95, labels = mylabel)

# %%%% TOTAL ALLELIC AND GEOGRAPHIC COVERAGE
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab="Coverage (%)")
# Add title and x-axis labels to the graph
title(main="Quercus acerifolia: Geo-Gen Coverage", line=1.5)
mtext(text="91 Individuals; 1 km buffer (individuals); 5 replicates", side=3, line=0.3)
mtext(text="Number of individuals", side=1, line=2.4)
# Mark the 95% threshold line, and the genetic/geographic points
abline(h=95, col="black", lty=3) 
abline(v=gen_min95Value, col="red")
abline(v=geo_min95Value, col="darkblue")
# Add text for the minimum sampling size lines
mtext(text=paste0("Gen 95% MSSE = ", gen_min95Value),
      side=1, line=-1.5, at=76, cex=1)
mtext(text=paste0("Geo 95% MSSE = ", geo_min95Value),
      side=1, line=-1.5, at=10, cex=1)
# Add legend
legend(x=65, y=80, inset = 0.05,
       legend = c("Genetic coverage (Total)", "Geographic coverage (1 km buffer IND)"),
       col=plotColors_Sub, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty="n", y.intersp = 0.75)

# %%%% TOTAL ALLELIC AND GEOGRAPHIC COVERAGE: 3 SAMPLE EMPHASIS
# (For IMLS NLG subgroup presentation, 2023-08-17)
# Use the matplot function to 3 average values, with specified settings
matplot(averageValueMat[1:3,], col=plotColors_Sub, pch=16, xlim=c(0,100), ylim=c(0,100), ylab = "")
# Add title and x-axis labels to the graph
title(main="Quercus acerifolia: Geo-Gen Coverage", line=1.5)
mtext(text="91 Individuals; 1 km buffer (individuals); 5 replicates", side=3, line=0.3, cex=1.2)
mtext(text="Number of individuals", side=1, line=3, cex=1.6)
mtext(text="Coverage (%)", side=2, line=2.3, cex=1.6, srt=90)
# Mark the 95% threshold line, and the genetic/geographic points
abline(h=95, col="black", lty=3) 
abline(v=3, col="black") 
# Add text for 3 sample example
mtext(text="COVERAGE VALUES AT 3 (RANDOM) SAMPLES", side=1, line=-4.5, at=24.8, cex=1.2)
mtext(text="Genetic coverage: 65.68%", side=1, line=-2.5, at=15.7, cex=1.2)
mtext(text="Geographic coverage: 51.05%", side=1, line=-1.5, at=17.5, cex=1.2)
# Add legend
legend(x=58, y=70, inset = 0.05,
       legend = c("Genetic coverage", "Geographic coverage (1 km buffer)"),
       col=plotColors_Sub, pch = c(20,20,20), cex=1.1, pt.cex = 2, bty="n", y.intersp = 0.75)

# Use the matplot function to plot the entire matrix of average values, with specified settings
matplot(averageValueMat, col=plotColors_Sub, pch=16, xlim=c(0,100), ylim=c(0,100), ylab = "")
# Add title and x-axis labels to the graph
title(main="Quercus acerifolia: Geo-Gen Coverage", line=1.5)
mtext(text="91 Individuals; 1 km buffer (individuals); 5 replicates", side=3, line=0.3, cex=1.2)
mtext(text="Number of individuals", side=1, line=3, cex=1.6)
mtext(text="Coverage (%)", side=2, line=2.3, cex=1.6, srt=90)
# Mark the 95% threshold line, and the genetic/geographic points
abline(h=95, col="black", lty=3) 
# Add legend
legend(x=58, y=70, inset = 0.05,
       legend = c("Genetic coverage", "Geographic coverage (1 km buffer)"),
       col=plotColors_Sub, pch = c(20,20,20), cex=1.1, pt.cex = 2, bty="n", y.intersp = 0.75)

# %%%% POPULATION-LEVEL GEOGRAPHIC COORDINATES %%%% ----
# In this analysis, we utilize a CSV file of lat/longs that uses the same value for each individual 
# in a given population. Essentially, there are only 4 unique combinations of latitude and longitude
# ---- READ IN GEOGRAPHIC AND GENETIC DATA ----
# Specify filepath for QUAC geographic and genetic data
QUAC.filePath <- paste0(GeoGenCorr.wd, "Datasets/QUAC/")

# ---- GEOGRAPHIC
# Read in wild occurrence points. This CSV has 3 columns: sample name, latitude, and longitude. 
# The sample names (and order) have to match the sample names/order of the genind object 
# (rownams of the genetic matrix) read in below.
wildPoints <- read.csv(paste0(QUAC.filePath, "QUAC_coord_pop.csv"), header=TRUE)

# ---- GENETIC
# Read in genind file: Optimized de novo assembly; R80, min-maf=0, 
# first SNP/locus, 2 populations (garden and wild), no Kessler individuals.
# Wild sample names/order must match those in the sample name column of the CSV (above)
QUAC.genind <- read.genepop(paste0(QUAC.filePath,"QUAC_DNFA_populations_R80_NOMAF_1SNP_2Pops_NoK.gen"))
# Correct popNames of genind. For this analysis, we'll only utilize wild samples (i.e. those in pop "wild")
pop(QUAC.genind) <- 
  factor(read.table(paste0(QUAC.filePath, "QUAC_popmap_GardenWild_NoK"), header=FALSE)[,2])

# ---- RESAMPLING ----
# Export the coordinate points data.frame and genind object to the cluster
clusterExport(cl, varlist = c("wildPoints","QUAC.genind","num_reps","buffSize","world_poly_clip_W"))
clusterExport(cl, varlist = c("createBuffers", "compareBuffArea", "getAlleleCategories","calculateCoverage",
                              "exSituResample", "geo.gen.Resample.Parallel"))
# Specify file path, for saving resampling array
arrayDir <- paste0(QUAC.filePath, "resamplingData/QUAC_50kmPOP_5r_resampArr.Rdata")
# Run resampling
QUAC_demoArray_Par <- geo.gen.Resample.Parallel(gen_obj=QUAC.genind, geo_coordPts=wildPoints,
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
gen_min95Value <- gen_min95Mean(QUAC_demoArray_Par) ; gen_min95Value
gen_min95SD(QUAC_demoArray_Par)
geo_min95Value <- geo_min95Mean(QUAC_demoArray_Par) ; geo_min95Value
geo_min95SD(QUAC_demoArray_Par)
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic and geographic proportions
averageValueMat <- meanArrayValues(QUAC_demoArray_Par)
# Generate linear model between both coverage values
QUAC_model <- lm (Total ~ Geo, data=averageValueMat)
QUAC_model_summary <- summary(QUAC_model) ; QUAC_model_summary
# Pull R-squared and p-value estimates from model
QUAC_model_rSquared <- QUAC_model_summary$adj.r.squared ; QUAC_model_rSquared
QUAC_model_pValue <- QUAC_model_summary$coefficients[2, 4]

# %%%% GEOGRAPHIC-GENETIC CORRELATION
plot(averageValueMat$Geo, averageValueMat$Total, pch=20, main="Q. acerifolia: Geographic by genetic coverage",
     xlab="Geographic coverage (%)", ylab="Genetic coverage (%)")
mtext(text="91 Individuals; 50 km buffer (populations); 5 replicates", side=3, line=0.3)
mylabel = bquote(italic(R)^2 == .(format(QUAC_model_rSquared, digits = 3)))
text(x = 45, y = 95, labels = mylabel)

# %%%% TOTAL ALLELIC AND GEOGRAPHIC COVERAGE
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab="Coverage (%)")
# Add title and x-axis labels to the graph
title(main="Quercus acerifolia: Geo-Gen Coverage", line=1.5)
mtext(text="91 Individuals; 50 km buffer (populations); 5 replicates", side=3, line=0.3)
mtext(text="Number of individuals", side=1, line=2.4)
# Mark the 95% threshold line, and the genetic/geographic points
abline(h=95, col="black", lty=3) 
abline(v=gen_min95Value, col="red")
abline(v=geo_min95Value, col="darkblue")
# Add text for the minimum sampling size lines
mtext(text=paste0("Gen 95% MSSE = ", gen_min95Value),
      side=1, line=-1.5, at=76, cex=1)
mtext(text=paste0("Geo 95% MSSE = ", geo_min95Value),
      side=1, line=-1.5, at=10, cex=1)
# Add legend
legend(x=65, y=80, inset = 0.05,
       legend = c("Genetic coverage (Total)", "Geographic coverage (50 km buffer POP)"),
       col=plotColors_Sub, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty="n", y.intersp = 0.75)
