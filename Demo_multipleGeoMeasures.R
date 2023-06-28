# Demo geographic resampling script

# load packages
my.packages <- c('tidyverse','textclean','terra','leaflet','rnaturalearth','Polychrome', 'adegenet')
lapply(my.packages, require, character.only=TRUE); rm(my.packages)

# Read in relevant functions
GeoGenCorr.wd <- "/home/akoontz/Documents/GeoGenCorr/Code/"
setwd(GeoGenCorr.wd)
source("functions_GeoGenCoverage.R")
# Declare a directory within which to store .Rdata objects of resampling arrays, and plots
resamplingDataDir <- paste0(GeoGenCorr.wd, "resamplingData/")
# Dimensions for resmpling plots
plotWidth <- 1262 ; plotHeight <- 734
# Plotting colors (for all plots!), with transparency for values other than Total (alpha = 1.0 is transparent)
plotColors <- c("red","red4","darkorange3","coral","purple")
plotColors[2:5] <- alpha(plotColors[2:5], 0.6)
# Set up relevant cores, and make sure adegenet library is present on cluster
num_cores <- detectCores() - 8 ; cl <- makeCluster(num_cores)
clusterEvalQ(cl, library("adegenet"))
# Specify number of replicates. This value is used uniformly, for ALL datasets (QUAC and QUBO)
num_reps <- 5000
# Export relevant functions and variables
clusterExport(cl, varlist = c("getAlleleCategories", "exSitu_Sample", "exSitu_Resample", "num_reps"))


################################################################################
# Set up standards to use throughout
################################################################################
# Specify buffer size in meters (i.e. 50 km)
med_buff <- 1000    

# define projections
  #	points will be WGS84
pt.proj <- "+proj=longlat +datum=WGS84"
  # for calculations, we need something with units in meters and equal area
calc.proj <- "+proj=eqearth +datum=WGS84"

################################################################################
# Read in and prep polygon data
################################################################################
# read in world countries layer created in 1-prep_gis_layers.R
# this will be used to clip buffers so they're not in the water
world_poly_clip <- vect(file.path("~/Documents/GeoGenCorr/QUAC_demo/gis_layers/world_countries_10m/world_countries_10m.shp"))

# Read in genind file: Optimized de novo assembly; R80, min-maf=0, first SNP/locus, 2 populations (garden and wild), no Kessler individuals
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops_NoK/"
QUAC.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R80.genind) <- 
  factor(read.table(paste0(genpop.filePath, "QUAC_popmap_GardenWild_NoK"), header=FALSE)[,2])

################################################################################
## Calculate & map geographic and ecological coverage of ex situ collections
################################################################################
## read in occurrence points (includes ex situ)
insitu_pt <- read.csv("~/Documents/GeoGenCorr/QUAC_demo/Quercus_acerifolia.csv", header = TRUE)

exSitu_Resample(QUAC.SNP.DN.R80.genind, insitu_pt, med_buff, pt.proj, calc.proj, world_poly_clip)

exS
