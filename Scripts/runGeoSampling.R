pacman::p_load(terra, dplyr, readr, furrr, tictoc, purrr)

# source and global vars  -------------------------------------------------
buffDists <- c(1000, 5000, 10000, 25000, 50000, 100000, 250000)
source("Scripts/geoSamplingFunctions.R")

# preprocessing data  -----------------------------------------------------
# allSpecies <- c("MIGU", "PICO", "QUAC", "QULO", "YUBR", "AMTH", "ARTH", "COGL", "HIWA", "VILA")
taxon <- c( "VILA")
# species with no coordinates at the moment  "AMTH", "COGL", "HIWA",
prepData(species = taxon)

# individual species processing  ------------------------------------------
# taxon <- c("ARTH")

# Create a data frame of all parameter combinations
params <- tidyr::expand_grid(species = taxon, buffDist = buffDists)

# trying just a for loop
testing = FALSE
if(testing){
  for(i in 1:nrow(params)){
    print(i)
    runGeoSelection(
      buffDist = params$buffDist[i],
      species = params$species[i],
      area = 0
    )
  }
}


  ## removing this startegry for now for troubleshooting 
future::plan(strategy = "multicore", workers = 8)
# iterate over the two columns
furrr::future_walk2(
  .x = params$buffDist,
  .y = params$species,
  .f = runGeoSelection,
  area = 0
)


# DELETE ALL the THINGS or just a specific species
## change doit to TRUE and set species
## species == NA will include all species, species == "MIGU" with grep only files from the taxon
## tried to make it safe as gone is forever in this case
# for (i in c("ARTH")) {
#   #  c("MIGU", "PICO", "QUAC", "YUBR")
#   removeFiles(doit = TRUE, species = i)
# }
