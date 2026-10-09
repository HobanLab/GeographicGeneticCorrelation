# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# %%% RUN GEOGRAPHIC MAXIMIZATION %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This script calls the functions declared in the geoSamplingFunctions.R script
# and iterates runGeoSelection over a vector of buffer sizes and a vector of 
# species acronyms. This results in the CSVs stored in the GeoCoreSets folder
# of each species.

pacman::p_load(terra, dplyr, readr, furrr, tictoc, purrr)

# source and global vars  -------------------------------------------------
# Specify buffer sizes (in meters)
buffDists <- c(1000, 5000, 10000, 25000, 50000, 100000, 250000)
# Read in functions for geographic maximization
source("Scripts/geoSamplingFunctions.R")

# preprocessing data  -----------------------------------------------------
taxon <- c("MIGU", "PICO", "QUAC", "QULO", "YUBR", "AMTH", "ARTH", "COGL", "HIWA", "VILA")
prepData(species = taxon)
# Create a data frame of all parameter combinations
params <- tidyr::expand_grid(species = taxon, buffDist = buffDists)

# Initialize parallelization
future::plan(strategy = "multicore", workers = 8)
# iterate over the two columns
furrr::future_walk2(
  .x = params$buffDist,
  .y = params$species,
  .f = runGeoSelection,
  area = 0
)
