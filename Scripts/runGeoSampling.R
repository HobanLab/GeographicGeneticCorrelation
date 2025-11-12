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
taxon <- c("YUBR")

# Create a data frame of all parameter combinations
params <- tidyr::expand_grid(species = taxon, buffDist = buffDists)

# iterate over the two columns
furrr::future_map2(
  .x = params$buffDist,
  .y = params$species,
  .f = runGeoSelection,
  area = 0
)

# troubleshooting -- more specific control
p1 <- params
# p1 <- params[params$species == "QULO",]
for(i in 1:nrow(p1)){
  print(i)
  runGeoSelection(buffDist = p1$buffDist[i], species = p1$species[i],area = 0)
}

# DELETE ALL the THINGS or just a specific species 
## change doit to TRUE and set species 
## species == NA will include all species, species == "MIGU" with grep only files from the taxon 
## tried to make it safe as gone is forever in this case
# removeFiles(doit = TRUE,species = "YUBR")



