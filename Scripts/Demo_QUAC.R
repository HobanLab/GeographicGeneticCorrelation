# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION DEMO: QUERCUS ACERIFOLIA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script demonstrating approach for calculating the correlation between genetic, geographic, 
# and ecological coverage. Uses a genind file from Quercus acerifolia (SNP loci, Complete dataset), 
# as well as a CSV containing sample names, latitudes and longitudes, to iteratively resample 
# wild points and measure coverage. Based on the parFlag value, the code is branched to 
# allow for procesing in parallel or not in parallel.

library(adegenet)
library(terra)
library(parallel)
library(RColorBrewer)
library(scales)
library(rnaturalearth)

# Read in relevant functions
GeoGenCorr_wd <- '~/Documents/GeographicGeneticCorrelation/'
setwd(GeoGenCorr_wd)
source('Scripts/functions_GeoGenCoverage.R')
source('Scripts/worldAdmin.R')

# ---- VARIABLES ----
# Specify number of resampling replicates
num_reps <- 3
# ---- BUFFER SIZES
# Specify geographic buffer size in meters 
geo_buffSize <- 1000
# Specify ecological buffer size in meters 
eco_buffSize <- 1000
# ---- SHAPEFILES
# Read in world countries layer (created as part of the gap analysis workflow)
# This layer is used to clip buffers, to make sure they're not in the water
world_poly_clip <- grabWorldAdmin(GeoGenCorr_wd = GeoGenCorr_wd,
                                  fileExtentsion = ".gpkg",
                                  overwrite = FALSE)


wildPoints <- read.csv(paste0(GeoGenCorr_wd,
                              '/Datasets/QUAC/Geographic/QUAC_coord_ind.csv'),
                       header=TRUE)

# perfor geographic filter on the admin layer
world_poly_clip <- prepWorldAdmin(world_poly_clip = world_poly_clip,
                        wildPoints = wildPoints) 

# ---- rasters 
sdm <- terra::rast(paste0(GeoGenCorr_wd,
                          '/Datasets/QUAC/Geographic/prj_threshold.tif'))

### test resolution of buffer against raster features 
## assuming that 1 meter == 0.000012726903908907691 degrees 
geoBuffDegree <- geo_buffSize * 0.000012726903908907691
sdmDegree <- terra::res(sdm)[1]
if(geoBuffDegree < sdmDegree){
  warning("Your input SDM has a resolution that is larger then your suggested buffer.
          The SDM will be resample to a smaller resolution.")
  # resample the raster to a smaller cell size 
  scaleFactor <- ceiling(sdmDegree /geoBuffDegree)
  sdm <- terra::disagg(sdm, scaleFactor) 
}


# defined again in the primary workflow. 
rm(wildPoints)




# Read in the EPA Level IV ecoregion shapefile, which is used for calculating ecological coverage (solely in the U.S.)
## no direct way to download so these need to be grad from epa website https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states
## seems like the data with states was used originally, I don't that that reference is utialize so it might be 
## better to use the smaller file. 
ecoregion_poly <- 
  vect(file.path(paste0(GeoGenCorr_wd, 'GIS_shpFiles/ecoregions_EPA_level4/us_eco_l4.shp')))


# ---- PARALLELIZATION
# Flag for running resampling steps in parallel
parFlag <- FALSE

# If running in parallel, set up cores and export required libraries
if(parFlag==TRUE){
  # Set up relevant cores 
  num_cores <- detectCores() - 4 
  cl <- makeCluster(num_cores)
  # Make sure libraries (adegenet, terra, and parallel) are on cluster
  clusterEvalQ(cl, library('adegenet'))
  clusterEvalQ(cl, library('terra'))
  clusterEvalQ(cl, library('parallel'))
  
  # Shapefiles are by default a 'non-exportable' object, which means the must be processed before being
  # exported to the cluster (for parallelized calculations). The terra::wrap function is used to do this.
  world_poly_clip_W <- wrap(world_poly_clip)
  ecoregion_poly_W <- wrap(ecoregion_poly)
  sdm_W <- wrap(sdm)
}

###
# I wonder if these parameters above are going to be consistent across species. Maybe best to move it into a 
# config files that source across the inputs. Also might not be such a big deal to duplicate between species. 
#
###



# %%%% INDIVIDUAL-LEVEL GEOGRAPHIC COORDINATES %%%% ----
# In this analysis, we utilize a CSV file of lat/longs that specify the location of each individual

# ---- READ IN DATA ----
# Specify filepath for QUAC geographic and genetic data
QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')

# ---- GENETIC MATRIX
# Read in genind file: Optimized de novo assembly; R80, min-maf=0, 
# first SNP/locus, 2 populations (garden and wild), no Kessler individuals.
# Wild sample names/order must match those in the sample name column of the CSV (above)
QUAC_genind <- read.genepop(paste0(QUAC_filePath,'Genetic/QUAC_DNFA_populations_R80_NOMAF_1SNP_2Pops_NoK.gen'))
# Correct popNames of genind. For this analysis, we'll only utilize wild samples (i.e. those in pop 'wild')
pop(QUAC_genind) <- 
  factor(read.table(paste0(QUAC_filePath, 'Genetic/QUAC_popmap_GardenWild_NoK'), header=FALSE)[,2])

# ---- GEOGRAPHIC COORDINATES
# Read in wild occurrence points. This CSV has 3 columns: sample name, latitude, and longitude. 
# The sample names (and order) have to match the sample names/order of the genind object 
# (rownams of the genetic matrix) read in below.
wildPoints <- read.csv(paste0(QUAC_filePath, 'Geographic/QUAC_coord_ind.csv'), header=TRUE)

# ---- RESAMPLING ----
if(parFlag==TRUE){
  # Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
  clusterExport(cl, varlist = c('wildPoints','QUAC_genind','num_reps','geo_buffSize', 'eco_buffSize',
                                'world_poly_clip_W', 'ecoregion_poly_W', 'sdm_W'))
  # Export necessary functions (for calculating geographic and ecological coverage) to the cluster
  clusterExport(cl, varlist = c('createBuffers', 'geo.compareBuff', "geo.compareBuffSDM",
                                'eco.intersectBuff', 'eco.compareBuff',
                                'gen.getAlleleCategories','calculateCoverage', 'exSituResample.Par', 
                                'geo.gen.Resample.Par'))
  # Specify file path, for saving resampling array
  arrayDir <- paste0(QUAC_filePath, 'resamplingData/QUAC_1kmIND_GE_3r_resampArr.Rdata')
  # Run resampling in parallel
  QUAC_demoArray_IND_Par <- 
    geo.gen.Resample.Par(
      gen_obj = QUAC_genind,
      geoFlag = TRUE,
      coordPts = wildPoints,
      SDMrast = sdm_W,
      geoBuff = geo_buffSize,
      boundary = world_poly_clip_W,
      ecoFlag = TRUE,
      ecoBuff = eco_buffSize,
      ecoRegions = ecoregion_poly_W,
      ecoLayer = 'US',
      reps = num_reps,
      arrayFilepath = arrayDir,
      cluster = cl
    )
  # Close cores
  stopCluster(cl)
} else{
  # Run resampling not in parallel (for function testing purposes)
  QUAC_demoArray_IND <-
    geo.gen.Resample(
      gen_obj = QUAC_genind,
      SDMrast = sdm,
      geoFlag = TRUE,
      coordPts = wildPoints,
      geoBuff = geo_buffSize,
      boundary = world_poly_clip,
      ecoFlag = FALSE,
      reps = num_reps
    )
}

# ---- CORRELATION ----
# Build a data.frame from array values, to pass to linear models
QUAC_DF <- resample.array2dataframe(QUAC_demoArray_IND_Par)

# ---- LINEAR MODELS
# Generate linear models, using Total allelic coverage as the response variable
# GEOGRAPHIC COVERAGE AS PREDICTOR VARIABLE
QUAC_geoModel <- lm (Total ~ Geo, data=QUAC_DF)
QUAC_geoModel_summary <- summary(QUAC_geoModel) ; QUAC_geoModel_summary
# Pull R-squared and p-value estimates from model
QUAC_geoModel_rSquared <- round(QUAC_geoModel_summary$adj.r.squared,2)
QUAC_geoModel_pValue <- QUAC_geoModel_summary$coefficients[2, 4]

# ---- PLOTTING ----
# ---- CALCULATE 95% MSSE AND AVERAGE VALUES
# Calculate minimum 95% sample size for genetic and geographic values
gen_min95Value <- gen.min95Mean(QUAC_demoArray_Par) ; gen_min95Value
gen_min95SD(QUAC_demoArray_Par)
geo_min95Value <- geo.min95Mean(QUAC_demoArray_Par) ; geo_min95Value
geo_min95SD(QUAC_demoArray_Par)
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic and geographic proportions
averageValueMat <- meanArrayValues(QUAC_demoArray_Par)

# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue')
plotColors[2:5] <- alpha(plotColors[2:5], 0.35)
plotColors_Sub <- plotColors[-(2:5)]

# ---- CORRELATION PLOTS
plot(averageValueMat$Geo, averageValueMat$Total, pch=20, main='Q. acerifolia: Geographic by genetic coverage',
     xlab='Geographic coverage (%)', ylab='Genetic coverage (%)')
mtext(text='91 Individuals; 1 km buffer (individuals); 5 replicates', side=3, line=0.3)
mylabel = bquote(italic(R)^2 == .(format(QUAC_model_rSquared, digits = 3)))
text(x = 45, y = 95, labels = mylabel)

# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='Coverage (%)')
# Add title and x-axis labels to the graph
title(main='Quercus acerifolia: Geo-Gen Coverage', line=1.5)
mtext(text='91 Individuals; 1 km buffer (individuals); 5 replicates', side=3, line=0.3)
mtext(text='Number of individuals', side=1, line=2.4)
# Mark the 95% threshold line, and the genetic/geographic points
abline(h=95, col='black', lty=3) 
abline(v=gen_min95Value, col='red')
abline(v=geo_min95Value, col='darkblue')
# Add text for the minimum sampling size lines
mtext(text=paste0('Gen 95% MSSE = ', gen_min95Value),
      side=1, line=-1.5, at=76, cex=1)
mtext(text=paste0('Geo 95% MSSE = ', geo_min95Value),
      side=1, line=-1.5, at=10, cex=1)
# Add legend
legend(x=65, y=80, inset = 0.05,
       legend = c('Genetic coverage (Total)', 'Geographic coverage (1 km buffer IND)'),
       col=plotColors_Sub, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.75)

# %%%% 2023-08-17 TOTAL ALLELIC AND GEOGRAPHIC COVERAGE: 3 SAMPLE EMPHASIS
# (For IMLS NLG subgroup presentation, 2023-08-17)
# Use the matplot function to 3 average values, with specified settings
matplot(averageValueMat[1:3,], col=plotColors_Sub, pch=16, xlim=c(0,100), ylim=c(0,100), ylab = '')
# Add title and x-axis labels to the graph
title(main='Quercus acerifolia: Geo-Gen Coverage', line=1.5)
mtext(text='91 Individuals; 1 km buffer (individuals); 5 replicates', side=3, line=0.3, cex=1.2)
mtext(text='Number of individuals', side=1, line=3, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Mark the 95% threshold line, and the genetic/geographic points
abline(h=95, col='black', lty=3) 
abline(v=3, col='black') 
# Add text for 3 sample example
mtext(text='COVERAGE VALUES AT 3 (RANDOM) SAMPLES', side=1, line=-4.5, at=24.8, cex=1.2)
mtext(text='Genetic coverage: 65.68%', side=1, line=-2.5, at=15.7, cex=1.2)
mtext(text='Geographic coverage: 51.05%', side=1, line=-1.5, at=17.5, cex=1.2)
# Add legend
legend(x=58, y=70, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (1 km buffer)'),
       col=plotColors_Sub, pch = c(20,20,20), cex=1.1, pt.cex = 2, bty='n', y.intersp = 0.75)

# Use the matplot function to plot the entire matrix of average values, with specified settings
matplot(averageValueMat, col=plotColors_Sub, pch=16, xlim=c(0,100), ylim=c(0,100), ylab = '')
# Add title and x-axis labels to the graph
title(main='Quercus acerifolia: Geo-Gen Coverage', line=1.5)
mtext(text='91 Individuals; 1 km buffer (individuals); 5 replicates', side=3, line=0.3, cex=1.2)
mtext(text='Number of individuals', side=1, line=3, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Mark the 95% threshold line, and the genetic/geographic points
abline(h=95, col='black', lty=3) 
# Add legend
legend(x=58, y=70, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (1 km buffer)'),
       col=plotColors_Sub, pch = c(20,20,20), cex=1.1, pt.cex = 2, bty='n', y.intersp = 0.75)

# %%%% 2023-09-27 TOTAL ALLELIC AND GEOGRAPHIC COVERAGE: 3 SAMPLE EMPHASIS
# (For IMLS NLG subgroup presentation, 2023-09-27)
# Alter the values in the averageValueMat, to correspond with the presentation
averageValueMat[1:3,1] <- c(70.2, 75.3, 80.9)
# Use the matplot function to 3 average values, with specified settings
matplot(averageValueMat[1:3,], col=plotColors_Sub, pch=16, xlim=c(0,100), ylim=c(0,100), ylab = '')
# Add title and x-axis labels to the graph
title(main='Quercus acerifolia: Geo-Gen Coverage', line=1.5)
mtext(text='91 Individuals; 1 km buffer (individuals); 5 replicates', side=3, line=0.3, cex=1.2)
mtext(text='Number of individuals', side=1, line=3, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Mark the 95% threshold line, and the genetic/geographic points
abline(h=95, col='black', lty=3) 
abline(v=3, col='black') 
# Add text for 3 sample example
mtext(text='COVERAGE VALUES AT 3 (RANDOM) SAMPLES', side=1, line=-4.5, at=24.8, cex=1.2)
mtext(text='Genetic coverage: 80.9%', side=1, line=-2.5, at=15.7, cex=1.2)
mtext(text='Geographic coverage: 51.05%', side=1, line=-1.5, at=17.5, cex=1.2)
# Add legend
legend(x=58, y=70, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (1 km buffer)'),
       col=plotColors_Sub, pch = c(20,20,20), cex=1.1, pt.cex = 2, bty='n', y.intersp = 0.75)

# Use the matplot function to plot the entire matrix of average values, with specified settings
matplot(averageValueMat, col=plotColors_Sub, pch=16, xlim=c(0,100), ylim=c(0,100), ylab = '')
# Add title and x-axis labels to the graph
title(main='Quercus acerifolia: Geo-Gen Coverage', line=1.5)
mtext(text='91 Individuals; 1 km buffer (individuals); 5 replicates', side=3, line=0.3, cex=1.2)
mtext(text='Number of individuals', side=1, line=3, cex=1.6)
mtext(text='Coverage (%)', side=2, line=2.3, cex=1.6, srt=90)
# Mark the 95% threshold line, and the genetic/geographic points
abline(h=95, col='black', lty=3) 
# Add legend
legend(x=58, y=70, inset = 0.05,
       legend = c('Genetic coverage', 'Geographic coverage (1 km buffer)'),
       col=plotColors_Sub, pch = c(20,20,20), cex=1.1, pt.cex = 2, bty='n', y.intersp = 0.75)

# %%%% POPULATION-LEVEL GEOGRAPHIC COORDINATES %%%% ----
# In this analysis, we utilize a CSV file of lat/longs that uses the same value for each individual 
# in a given population. Essentially, there are only 4 unique combinations of latitude and longitude

# ---- READ IN DATA ----
# Specify filepath for QUAC geographic and genetic data
QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')

# ---- GENETIC MATRIX
# Read in genind file: Optimized de novo assembly; R80, min-maf=0, 
# first SNP/locus, 2 populations (garden and wild), no Kessler individuals.
# Wild sample names/order must match those in the sample name column of the CSV (above)
QUAC.genind <- read.genepop(paste0(QUAC.filePath,'Genetic/QUAC_DNFA_populations_R80_NOMAF_1SNP_2Pops_NoK.gen'))
# Correct popNames of genind. For this analysis, we'll only utilize wild samples (i.e. those in pop 'wild')
pop(QUAC.genind) <- 
  factor(read.table(paste0(QUAC.filePath, 'Genetic/QUAC_popmap_GardenWild_NoK'), header=FALSE)[,2])

# ---- GEOGRAPHIC COORDINATES
# Read in wild occurrence points. This CSV has 3 columns: sample name, latitude, and longitude. 
# The sample names (and order) have to match the sample names/order of the genind object 
# (rownams of the genetic matrix) read in below.
wildPoints <- read.csv(paste0(QUAC.filePath, 'Geographic/QUAC_coord_pop.csv'), header=TRUE)

# ---- RESAMPLING ----
# Export necessary objects (genind, coordinate points, buffer size variables, polygons) to the cluster
clusterExport(cl, varlist = c('wildPoints','QUAC.genind','num_reps','geo_buffSize', 'eco_buffSize',
                              'world_poly_clip_W', 'ecoregion_poly_W'))
# Export necessary functions (for calculating geographic and ecological coverage) to the cluster
clusterExport(cl, varlist = c('createBuffers', 'geo_compareBuff', 'eco_intersectBuff', 'eco_compareBuff',
                              'gen.getAlleleCategories','calculateCoverage', 'exSituResample', 
                              'geo.gen.Resample.Parallel'))
# Specify file path, for saving resampling array
arrayDir <- paste0(QUAC.filePath, 'resamplingData/QUAC_1kmPOP_GE_5r_resampArr.Rdata')
# Run resampling
QUAC_demoArray_POP_Par <- 
  geo.gen.Resample.Parallel(gen_obj = QUAC.genind, geoFlag = TRUE, coordPts = wildPoints, 
                            geoBuff = geo_buffSize, boundary=world_poly_clip_W, ecoFlag = TRUE, 
                            ecoBuff = eco_buffSize, ecoRegions = ecoregion_poly_W, ecoLayer = 'US', 
                            reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
# Close cores
stopCluster(cl)

# ---- CORRELATION ----
# Build a data.frame from array values, to pass to linear models
QUAC_DF <- resample.array2dataframe(QUAC_demoArray_POP_Par)

# ---- LINEAR MODELS
# Generate linear models, using Total allelic coverage as the response variable
# GEOGRAPHIC COVERAGE AS PREDICTOR VARIABLE
QUAC_geoModel <- lm (Total ~ Geo, data=QUAC_DF)
QUAC_geoModel_summary <- summary(QUAC_geoModel) ; QUAC_geoModel_summary
# Pull R-squared and p-value estimates from model
QUAC_geoModel_rSquared <- round(QUAC_geoModel_summary$adj.r.squared,2)
QUAC_geoModel_pValue <- QUAC_geoModel_summary$coefficients[2, 4]

# ---- PLOTTING ----
# ---- CALCULATE 95% MSSE AND AVERAGE VALUES
# Calculate minimum 95% sample size for genetic and geographic values
gen_min95Value <- gen.min95Mean(QUAC_demoArray_Par) ; gen_min95Value
gen_min95SD(QUAC_demoArray_Par)
geo_min95Value <- geo.min95Mean(QUAC_demoArray_Par) ; geo_min95Value
geo_min95SD(QUAC_demoArray_Par)
# Generate the average values (across replicates) for all proportions
# This function has default arguments for returning just Total allelic and geographic proportions
averageValueMat <- meanArrayValues(QUAC_demoArray_Par)

# Specify plot colors
plotColors <- c('red','red4','darkorange3','coral','purple', 'darkblue')
plotColors[2:5] <- alpha(plotColors[2:5], 0.35)
plotColors_Sub <- plotColors[-(2:5)]

# ---- CORRELATION PLOTS
plot(averageValueMat$Geo, averageValueMat$Total, pch=20, main='Q. acerifolia: Geographic by genetic coverage',
     xlab='Geographic coverage (%)', ylab='Genetic coverage (%)')
mtext(text='91 Individuals; 50 km buffer (populations); 5 replicates', side=3, line=0.3)
mylabel = bquote(italic(R)^2 == .(format(QUAC_model_rSquared, digits = 3)))
text(x = 45, y = 95, labels = mylabel)

# ---- COVERAGE PLOTS
# Use the matplot function to plot the matrix of average values, with specified settings
matplot(averageValueMat, ylim=c(0,100), col=plotColors_Sub, pch=16, ylab='Coverage (%)')
# Add title and x-axis labels to the graph
title(main='Quercus acerifolia: Geo-Gen Coverage', line=1.5)
mtext(text='91 Individuals; 50 km buffer (populations); 5 replicates', side=3, line=0.3)
mtext(text='Number of individuals', side=1, line=2.4)
# Mark the 95% threshold line, and the genetic/geographic points
abline(h=95, col='black', lty=3) 
abline(v=gen_min95Value, col='red')
abline(v=geo_min95Value, col='darkblue')
# Add text for the minimum sampling size lines
mtext(text=paste0('Gen 95% MSSE = ', gen_min95Value),
      side=1, line=-1.5, at=76, cex=1)
mtext(text=paste0('Geo 95% MSSE = ', geo_min95Value),
      side=1, line=-1.5, at=10, cex=1)
# Add legend
legend(x=65, y=80, inset = 0.05,
       legend = c('Genetic coverage (Total)', 'Geographic coverage (50 km buffer POP)'),
       col=plotColors_Sub, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty='n', y.intersp = 0.75)
