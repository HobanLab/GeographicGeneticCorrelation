# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEO-ECO-GEN CORRELATION DEMO: QUERCUS LOBATA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script demonstrating first draft approach for calculating the correlation 
# between genetic, geographic, and ecological coverage. 
# Uses data files from Gugger et al. 2021 to pull in genetic data 
# (as a table) and geographic coordinates (in a CSV) to conduct correlation analyses.

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
# ---- BUFFER SIZES
# Specify geographic buffer size in meters 
# geo_buffSize <- 1000
geo_buffSize <- 50000
# Specify ecological buffer size in meters 
# eco_buffSize <- 1000
eco_buffSize <- 50000
# ---- SHAPEFILES
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- 
  vect(file.path(paste0(GeoGenCorr.wd, "GIS_shpFiles/world_countries_10m/world_countries_10m.shp")))
# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage (solely in the U.S.)
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr.wd, "GIS_shpFiles/ecoregions_EPA_level4/us_eco_l4.shp")))
# Shapefiles are by default a "non-exportable" object, which means the must be processed before being
# exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
world_poly_clip_W <- wrap(world_poly_clip)
ecoregion_poly_W <- wrap(ecoregion_poly)

# ---- PARALLELIZATION
# Set up relevant cores 
num_cores <- detectCores() - 4 
cl <- makeCluster(num_cores)
# Make sure libraries (adegenet + terra) are on cluster
clusterEvalQ(cl, library("adegenet"))
clusterEvalQ(cl, library("terra"))
clusterEvalQ(cl, library("parallel"))

# ---- READ IN DATA ----
# Specify filepath for QUAC geographic and genetic data
QULO.filePath <- paste0(GeoGenCorr.wd, "Datasets/QULO/")

# ---- COORDINATE POINTS
# Read in wild occurrence points. This CSV has 3 columns: sample name, latitude, and longitude. 
# The sample names (and order) have to match the sample names/order of the genind object 
# (rownams of the genetic matrix) read in below.
wildPoints <- read.csv(paste0(QULO.filePath, "Quercus_lobata.csv"), header=TRUE)

# ---- GENETIC MATRIX
# Processing input file
QULO.tab <- read.table(paste0(QULO.filePath, "SNPs80.forR"), header = TRUE)
# Make sample names the row names
rownames(QULO.tab) <- QULO.tab[,1] ; QULO.tab <- QULO.tab[,-1]
# Convert to genind. ncode based off similar practice for other geninds (values of 1, 2, and 3 generate identical results)
QULO.genind <- df2genind(QULO.tab, ncode = 3, ploidy = 2)

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c("wildPoints","QULO.genind","num_reps","geo_buffSize", "eco_buffSize",
                              "world_poly_clip_W", "ecoregion_poly_W"))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c("createBuffers", "geo_compareBuff", "eco_intersectBuff", "eco_compareBuff",
                              "getAlleleCategories","calculateCoverage", "exSituResample", 
                              "geo.gen.Resample.Parallel"))
# Specify file path, for saving resampling array
# arrayDir <- paste0(QUAC.filePath, "resamplingData/QULO_1km_GE_5r_resampArr.Rdata")
arrayDir <- paste0(QULO.filePath, "resamplingData/QULO_50km_GE_5r_resampArr.Rdata")
# Run resampling (in parallel)
QULO_demoArray_Par <- 
  geo.gen.Resample.Parallel(gen_obj = QULO.genind, geoFlag = TRUE, coordPts = wildPoints, 
                            geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag = TRUE, 
                            ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer = "US", 
                            reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
# Close cores
stopCluster(cl)

# ---- CORRELATION AND PLOTTTING ----
# Specify plot colors
# plotColors <- c("red","red4","darkorange3","coral","purple", "darkblue")
# plotColors[2:5] <- alpha(plotColors[2:5], 0.35)
# plotColors_Sub <- plotColors[-(2:5)]
# 
# # Calculate minimum 95% sample size for genetic and geographic values
# gen_min95Value <- gen_min95Mean(QULO_demoArray_Par) ; gen_min95Value
# gen_min95SD(QULO_demoArray_Par)
# geo_min95Value <- geo_min95Mean(QULO_demoArray_Par) ; geo_min95Value
# geo_min95SD(QULO_demoArray_Par)
# # Generate the average values (across replicates) for all proportions
# # This function has default arguments for returning just Total allelic and geographic proportions
# averageValueMat <- meanArrayValues(QULO_demoArray_Par)
# # Generate linear model between both coverage values
# QULO_model <- lm (Total ~ Geo, data=averageValueMat)
# QULO_model_summary <- summary(QULO_model) ; QULO_model_summary
# # Pull R-squared and p-value estimates from model
# QULO_model_rSquared <- QULO_model_summary$adj.r.squared
# QULO_model_pValue <- QULO_model_summary$coefficients[2, 4]
# 
# # %%%% GEOGRAPHIC-GENETIC CORRELATION
# plot(averageValueMat$Geo, averageValueMat$Total, pch=20, 
#      main="Q. lobata: Geographic by genetic coverage",xlab="", ylab="")
# mtext(text="436 Individuals; 50 km buffer; 5 replicates", side=3, line=0.3, cex=1.3)
# mtext(text="Geographic coverage (%)", side=1, line=3, cex=1.6)
# mtext(text="Genetic coverage (%)", side=2, line=2.3, cex=1.6, srt=90)
# mylabel = bquote(italic(R)^2 == .(format(QULO_model_rSquared, digits = 3)))
# text(x = 20, y = 85, labels = mylabel, cex=1.4)
# 
# # %%%% TOTAL ALLELIC AND GEOGRAPHIC COVERAGE
# # Use the matplot function to plot the matrix of average values, with specified settings
# matplot(averageValueMat, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab="")
# # Add title and x-axis labels to the graph
# title(main="Quercus lobata: Geo-Gen Coverage", line=1.5)
# mtext(text="436 Individuals; 50 km buffer; 5 replicates", side=3, line=0.3, cex=1.3)
# mtext(text="Number of individuals", side=1, line=2.4, cex=1.6)
# mtext(text="Coverage (%)", side=2, line=2.3, cex=1.6, srt=90)
# # Mark the 95% threshold line, and the genetic/geographic points
# abline(h=95, col="black", lty=3) 
# abline(v=gen_min95Value, col="red")
# abline(v=geo_min95Value, col="darkblue")
# # Add text for the minimum sampling size lines
# mtext(text=paste0("Gen 95% MSSE = ", gen_min95Value),
#       side=1, line=-1.5, at=95, cex=1.3)
# mtext(text=paste0("Geo 95% MSSE = ", geo_min95Value),
#       side=1, line=-1.5, at=200, cex=1.3)
# # Add legend
# legend(x=205, y=60, inset = 0.05,
#        legend = c("Genetic coverage (Total)", "Geographic coverage (50 km buffer)"),
#        col=plotColors_Sub, pch = c(20,20,20), cex=1.2, pt.cex = 2, bty="n", y.intersp = 0.8)
# 
# # %%%% BOTH PLOTS (For IMLS NLG subgroup presentation, 2023-08-17) ----
# par(mfrow=c(2,1))
# 
# # %%%% TOTAL ALLELIC AND GEOGRAPHIC COVERAGE
# # Use the matplot function to plot the matrix of average values, with specified settings
# matplot(averageValueMat, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab="")
# # Add title and x-axis labels to the graph
# title(main="Quercus lobata: Geo-Gen Coverage", line=1.5)
# mtext(text="436 Individuals; 50 km buffer; 5 replicates", side=3, line=0.3, cex=1.3)
# mtext(text="Number of individuals", side=1, line=2.4, cex=1.6)
# mtext(text="Coverage (%)", side=2, line=2.3, cex=1.6, srt=90)
# # Mark the 95% threshold line, and the genetic/geographic points
# abline(h=95, col="black", lty=3) 
# # Add legend
# legend(x=215, y=65, inset = 0.05,
#        legend = c("Genetic coverage (Total)", "Geographic coverage (50 km buffer)"),
#        col=plotColors_Sub, pch = c(20,20,20), cex=1.2, pt.cex = 2, bty="n", y.intersp = 0.6)
# 
# # %%%% GEOGRAPHIC-GENETIC CORRELATION
# plot(averageValueMat$Geo, averageValueMat$Total, pch=20, main="",xlab="", ylab="")
# mtext(text="Geographic coverage (%)", side=1, line=3, cex=1.6)
# mtext(text="Genetic coverage (%)", side=2, line=2.3, cex=1.6, srt=90)
# mylabel = bquote(italic(R)^2 == .(format(QULO_model_rSquared, digits = 3)))
# text(x = 20, y = 85, labels = mylabel, cex=1.4)
