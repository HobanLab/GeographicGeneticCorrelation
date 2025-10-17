pacman::p_load(terra,dplyr,readr,furrr, tictoc, purrr)

# source and global vars  -------------------------------------------------
buffDists <- c(1000,5000,10000,25000,50000,100000,250000)
source("Scripts/geoSamplingFunctions.R")

# preprocessing data  -----------------------------------------------------
prepData()

# individual species processing  ------------------------------------------
future::plan(strategy = "multicore", workers = 8)
taxon <- c("MIGU", "PICO", "QUAC", "QULO", "YUBR")
# testing 
taxon <- c("PICO")

# Create a data frame of all parameter combinations
params <- tidyr::expand_grid(species = taxon, buffDist = buffDists)

# iterate over the two columns
furrr::future_map2(
  .x = params$buffDist,
  .y = params$species,
  .f = runGeoSelection
)

# troubleshooting -- more specific control
# for(i in 1:nrow(params)){
#   print(i)
#   runGeoSelection(buffDist = params$buffDist[i], species = params$species[i])
# }





