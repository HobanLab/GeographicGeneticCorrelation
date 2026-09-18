# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION: UNIQUE COORDINATE SUMMARY %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script summarizes, for each dataset, how many individuals are used in the genetic
# analyses and how many unique geographic locations those individuals represent. Several
# datasets report coordinates at the site level (or share coordinates between individuals),
# which means the number of distinct buffers available to the geographic/ecological coverage
# calculations, and to the geographically maximized sampling, is lower than the number of individuals.

# The coordinate files are used as the record of individuals in the genetic data: the resampling
# functions throw an error unless the sample names/order match the genind object exactly.

pacman::p_load(dplyr, readr, purrr)

# datasetLookup provides the coordinate file and column names for each dataset
source("Scripts/geoSamplingFunctions.R")

# order datasets as in Table 1 of the manuscript (increasing EOO)
allSpecies <- c("HIWA", "COGL", "YUBR", "QUAC", "AMTH", "VILA", "QULO", "PICO", "MIGU", "ARTH")

# summarize a single dataset ----------------------------------------------
summarizeCoordinates <- function(species, lookup) {
  vals <- lookup[lookup$taxon == species, ]
  lonLat <- vals$latLonCol[[1]]
  # PICO ids are read as character to keep the specific formatting
  d1 <- read_csv(
    paste0("Datasets/", species, "/Geographic/", vals$pointPath),
    col_types = cols(.default = col_guess(), !!vals$idCol := col_character()),
    show_col_types = FALSE,
    name_repair = "unique_quiet" # AMTH file has an unnamed row number column
  ) |>
    dplyr::select(id = all_of(vals$idCol), lon = all_of(lonLat[1]), lat = all_of(lonLat[2]))
  # individuals with coordinates
  d2 <- d1 |>
    dplyr::filter(!is.na(lon), !is.na(lat))
  # number of individuals at each unique location
  perLocation <- d2 |>
    dplyr::count(lon, lat, name = "individuals")

  dplyr::tibble(
    dataset = species,
    individuals = nrow(d1),
    missingCoordinates = nrow(d1) - nrow(d2),
    uniqueCoordinates = nrow(perLocation),
    percentUnique = round(nrow(perLocation) / nrow(d2) * 100, 1),
    meanIndPerLocation = round(mean(perLocation$individuals), 1),
    maxIndPerLocation = max(perLocation$individuals)
  )
}

# apply to all datasets ---------------------------------------------------
lookup <- datasetLookup()
# datasets with protected coordinates (AMTH, COGL, HIWA) are only present locally
available <- allSpecies[file.exists(paste0("Datasets/", allSpecies, "/Geographic/", allSpecies, "_coordinates.csv"))]
coordinateSummary <- purrr::map(.x = available, .f = summarizeCoordinates, lookup = lookup) |>
  dplyr::bind_rows()

print(coordinateSummary)
write_csv(x = coordinateSummary, file = "Datasets/coordinateSummary.csv")
