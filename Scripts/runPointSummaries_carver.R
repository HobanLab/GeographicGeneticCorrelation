# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION: POINT SUMMARY METRICS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load packages 
pacman::p_load(sf, sp, terra, dplyr, purrr, readr, tidyr, tibble, sfdep, Hmisc, corrplot)

# Read in relevant functions
source('Scripts/functions_GeoGenCoverage.R')
source('Scripts/geoSamplingFunctions.R')
# Run from the root of the repository (the .Rproj location); paths below are built from here
GeoGenCorr_wd <- paste0(getwd(), '/')

# CALCULATING POINT SUMMARY VALUES ----
# Datasets are kept in alphabetical order, as the optimal buffer sizes (see below) are
# combined with the point summaries by column position
allSpecies <- c('AMTH','ARTH','COGL','HIWA','MIGU','PICO','QUAC','QULO','VILA','YUBR')
# Number of individuals used in the genetic analyses of each dataset (Table 1 of the manuscript). 
# This is a check that the coordinate file matches the genetic data, and not an earlier or
# unfiltered version of the dataset
expectedInds <- c(AMTH=140, ARTH=1010, COGL=562, HIWA=197, MIGU=255, 
                  PICO=929, QUAC=91, QULO=436, VILA=157, YUBR=319)
# One coordinate file per dataset (<species>_coordinates.csv). Files such as 
# ARTH_coordinates_Original.csv or MIGU_coordinates_Global.csv include individuals that are
# not part of the analyses (e.g. outside of the native range) and are not used
lookup <- datasetLookup()

# Read in the coordinates for a dataset, and standardize the column names (these differ between files)
readPointsData <- function(species, lookup, expectedInds){
  vals <- lookup[lookup$taxon == species, ]
  lonLat <- vals$latLonCol[[1]]
  coordFile <- paste0(GeoGenCorr_wd, 'Datasets/', species, '/Geographic/', vals$pointPath)
  # name_repair: the AMTH file has an unnamed row number column
  d1 <- read_csv(coordFile, show_col_types = FALSE, name_repair = 'unique_quiet') |>
    dplyr::mutate(taxon = species) |>
    dplyr::select(taxon, lat = all_of(lonLat[2]), lon = all_of(lonLat[1]))
  # Check the number of individuals against the genetic data
  if(nrow(d1) != expectedInds[[species]]){
    stop(paste0(species, ': ', nrow(d1), ' rows in ', vals$pointPath, ', but ', 
                expectedInds[[species]], ' individuals expected'))
  }
  # Individuals without coordinates can not contribute to the point summaries. These are dropped 
  # here (rather than within geo.generateSpatialObject) so the number removed is reported
  noCoords <- is.na(d1$lat) | is.na(d1$lon)
  if(any(noCoords)){
    message(paste0(species, ': ', sum(noCoords), ' individuals without coordinates removed'))
  }
  return(d1[!noCoords, ])
}

# Build a list of the coordinates for each dataset
pointsDataList <- purrr::map(.x = allSpecies, .f = readPointsData, 
                             lookup = lookup, expectedInds = expectedInds)
names(pointsDataList) <- tolower(allSpecies)

# Summary of the inputs: individuals with coordinates, and unique locations, per dataset
pointsDataSummary <- purrr::map_dfr(pointsDataList, function(d){
  data.frame(taxon = d$taxon[1], individuals = nrow(d), uniqueLocations = nrow(dplyr::distinct(d, lat, lon)))
})
print(pointsDataSummary)


# ---- POINT SUMMARY FUNCTIONS ----
# The functions in this section are associated with the spatial metrics for buffer optimization (SMBO)
# analyses, which seek to determine whether an optimal buffer size for matching genetic and geographic
# coverage can be estimated using summaries of the geographic coordinates of each dataset. The functions
# were authored by Dan Carver.

#' geo.generateSpatialObject -- Dan Carver
#' Build a spatial object from coordinate data (columns lat, lon in decimal degrees, WGS84).
#' Points are projected to Mollweide, an equal-area projection, so that the area based metrics
#' (EOO, AOO, ELA) are measured in an equal-area projection as recommended by the IUCN Red List
#' guidelines. Distances (ANN, STD, ELP) are handled separately within each function.
geo.generateSpatialObject <- function(data){
  # Clean 
  data1 <- data |>
    dplyr::filter(!is.na(lon))|>
    dplyr::filter(!is.na(lat))
  # Lat long data
  sp1 <- sf::st_as_sf(x = data1,
                      coords = c("lon","lat"),
                      crs = 4326,
                      remove = FALSE)
  # Projected data 
  sp1_proj <- sf::st_transform(x = sp1, crs = "+proj=moll")
  return(sp1_proj)
}

#' geo.calc.EOO
#' Calculate the extent of occurrence (EOO): the area of the minimum convex polygon (convex hull)
#' around all locations, following the IUCN Red List guidelines. The hull is drawn and measured
#' in the equal-area projection of the input (see geo.generateSpatialObject). Units: km2
geo.calc.EOO <- function(data){
  # Minimum convex polygon around all locations
  mcp <- st_convex_hull(st_union(data))
  # Area in km2
  eooArea <- as.numeric(st_area(mcp)) / 1e6
  return(eooArea)
}

#' geo.calc.AOO
#' Calculate the area of occupancy (AOO). A grid of 2 km x 2 km cells (the IUCN Red List reference
#' scale) is placed over the minimum convex polygon, in the equal-area projection of the input, and
#' the cells containing at least one location are counted. Two values are returned:
#'   AOO_km2: number of occupied cells x 4 km2, the IUCN Red List AOO
#'   AOO_pct: occupied cells as a percentage of all cells intersecting the minimum convex polygon
#'            (the value reported in earlier drafts; a measure of occupancy relative to EOO)
#' Note that the cell counts depend on where the grid is placed; the grid is anchored at the lower
#' left corner of the bounding box of the locations.
geo.calc.AOO <- function(data, cellsize = 2000){
  # Minimum convex polygon around all locations
  mcp <- st_convex_hull(st_union(data))
  # Grid of 2 km cells covering the bounding box of the minimum convex polygon
  allGrids <- st_make_grid(mcp, square = TRUE, cellsize = cellsize)
  # Keep the cells intersecting the minimum convex polygon
  mcpCells <- allGrids[sf::st_intersects(x = allGrids, y = mcp, sparse = FALSE)]
  # Count the cells with at least one location
  occupied <- sf::st_intersects(x = mcpCells, y = data, sparse = TRUE) |>
    lengths()
  occupiedCells <- sum(occupied > 0)
  # Area of occupancy (km2) and percent of cells occupied
  cellArea <- (cellsize / 1000)^2
  output <- data.frame(AOO_km2 = occupiedCells * cellArea,
             AOO_pct = occupiedCells / length(mcpCells) * 100)
  return(output)
}

#' geo.calc.averageNearestNeighbor 
#' Calculate the average nearest neighbor distance (Clark and Evans 1954): for each location, the
#' distance to the closest other location, averaged over all locations. Distances are geodesic
#' (great circle distances on a sphere, via sf::st_distance with s2) so they do not depend on the projection. Individuals 
#' sharing the same coordinates are collapsed to a single location first; otherwise site level 
#' datasets (e.g. AMTH, PICO, VILA) would have nearest neighbor distances of zero. Units: km
geo.calc.averageNearestNeighbor <- function(data){
  # Unique locations, in lat lon for geodesic distances
  uniqueLoc <- data[!duplicated(sf::st_coordinates(data)), ] |>
    sf::st_transform(crs = 4326)
  # Pairwise distance matrix (m); a location is not its own neighbor
  distMat <- sf::st_distance(uniqueLoc) |>
    units::drop_units()
  diag(distMat) <- NA
  # Distance from each location to its nearest neighbor
  nearestDist <- apply(distMat, 1, min, na.rm = TRUE)
  # Average, in km
  ANN <- mean(nearestDist) / 1000
  return(ANN)
}

#' geo.calc.voronoiAreas
#' Calculate the evenness of the spatial sampling from a Voronoi tessellation. Each unique location
#' is assigned the area closer to it than to any other location; the tessellation is built in the
#' equal-area projection of the input and clipped to the minimum convex polygon (EOO). The metric
#' is the coefficient of variation (sd / mean) of the cell areas: evenly spaced locations produce
#' cells of similar size (values near 0), clustered sampling produces a few large and many small
#' cells (larger values). Note that the mean cell area itself is always EOO / number of locations,
#' regardless of the arrangement of the points, so it is not reported. Unitless
geo.calc.voronoiAreas <- function(data){
  # Unique locations (duplicated coordinates produce empty cells)
  uniqueLoc <- data[!duplicated(sf::st_coordinates(data)), ]
  # Minimum convex polygon, used to bound the outer cells
  mcp <- sf::st_convex_hull(sf::st_union(uniqueLoc))
  # Voronoi tessellation, clipped to the minimum convex polygon
  cells <- sf::st_voronoi(sf::st_union(uniqueLoc)) |>
    sf::st_collection_extract("POLYGON") |>
    sf::st_intersection(mcp)
  # Cell areas (km2) and their coefficient of variation
  cellAreas <- as.numeric(sf::st_area(cells)) / 1e6
  vorCV <- sd(cellAreas) / mean(cellAreas)
  return(vorCV)
}

#' geo.calc.stdDistance
#' Calculate the standard distance: the root mean square distance of the unique locations from
#' their mean center, a measure of the overall dispersion of the sampling (Bachi 1963). Calculated
#' in the equal-area projection of the input. Units: km
geo.calc.stdDistance <- function(data){
  # Unique locations
  uniqueLoc <- data[!duplicated(sf::st_coordinates(data)), ]
  # Standard distance (m) to km
  stdDist <- sfdep::std_distance(geometry = uniqueLoc) / 1000
  return(stdDist)
}

#' geo.calc.ellipseElongation
#' Calculate the elongation of the standard deviational ellipse (Yuill 1971) of the unique
#' locations: the ratio of the major to the minor axis (one standard deviation along each axis).
#' A value of 1 indicates an isotropic (circular) arrangement of the locations, larger values an
#' increasingly linear arrangement. Note that the area and perimeter of the ellipse, and the area
#' of the standard distance circle (pi * STD^2), are all size measures which rank the datasets in 
#' the same order as STD; the elongation is the shape information the ellipse adds. Unitless
geo.calc.ellipseElongation <- function(data){
  # Unique locations
  uniqueLoc <- data[!duplicated(sf::st_coordinates(data)), ]
  # Standard deviational ellipse (sx, sy: axis lengths; theta: rotation)
  sde <- sfdep::std_dev_ellipse(geometry = uniqueLoc)
  # Ratio of the major to the minor axis
  elongation <- max(sde$sx, sde$sy) / min(sde$sx, sde$sy)
  return(elongation)
}

# Wrapper function of the above point summary functions, which will calculate
# each point summary statistic for a given dataset, and return a single row data.frame.
# This replaces the version in functions_GeoGenCoverage.R (different set of metrics)
geo.calc.pointSummaries <- function(geoData){
  # Begin by converting data into spatial object
  geoSpat <- geo.generateSpatialObject(geoData)
  # Run through the different spatial metric calculations, and store results to a data.frame
  EOO <- geo.calc.EOO(geoSpat)
  AOO <- geo.calc.AOO(geoSpat)
  ANN <- geo.calc.averageNearestNeighbor(geoSpat)
  VOR <- geo.calc.voronoiAreas(geoSpat)
  STD <- geo.calc.stdDistance(geoSpat)
  ELG <- geo.calc.ellipseElongation(geoSpat)
  # Combine metrics into a data.frame, and return
  ptSummaryDF <- data.frame(EOO = EOO, AOO = AOO$AOO_km2, AOO_pct = AOO$AOO_pct, 
                            ANN = ANN, VOR = VOR, STD = STD, ELG = ELG)
  return(ptSummaryDF)
}


# TESTING THE METRICS AGAINST THE MANUSCRIPT (DRAFT 5) ----
# Each block below runs one metric on all datasets and compares it to the values currently in
# the manuscript. Set runMetricTests to FALSE to skip these and go straight to the full calculation
runMetricTests <- TRUE
if(runMetricTests){

# starting with EOO let's calll the function and test all results 
## Values currently reported in the manuscript (Draft 5, Table 1 and Table S4; km2)
manuscriptEOO <- c(AMTH=60767.70, ARTH=73591778.47, COGL=10.05, HIWA=8.00, MIGU=6529925.52,
                   PICO=1308901.79, QUAC=3506.07, QULO=135430.40, VILA=88917.24, YUBR=186.90)
eooTest <- purrr::map_dfr(pointsDataList, function(d){
  eoo <- geo.calc.EOO(geo.generateSpatialObject(d))
  data.frame(taxon = d$taxon[1], individuals = nrow(d), EOO_km2 = round(eoo, 2), 
             manuscript_km2 = manuscriptEOO[[d$taxon[1]]])
}) |>
  dplyr::mutate(percentDifference = round((EOO_km2 - manuscript_km2) / manuscript_km2 * 100, 1),
                rank = rank(EOO_km2), manuscriptRank = rank(manuscript_km2))
print(eooTest)

# AOO: run on the five datasets with the smallest ranges (the 2 km grid over the convex hull of the
# larger datasets takes minutes; ARTH ~22 minutes)
## Values currently reported in the manuscript (Draft 5, Table S4; percent of cells occupied)
manuscriptAOO <- c(AMTH=0.09017, ARTH=0.004025, COGL=83.33, HIWA=85.57, MIGU=0.006113,
                   PICO=0.09157, QUAC=0.517, QULO=0.5327, VILA=0.1799, YUBR=38.81)
smallSpecies <- c('hiwa','cogl','yubr','quac','amth')
aooTest <- purrr::map_dfr(pointsDataList[smallSpecies], function(d){
  aoo <- geo.calc.AOO(geo.generateSpatialObject(d))
  data.frame(taxon = d$taxon[1], individuals = nrow(d), 
             uniqueLocations = nrow(dplyr::distinct(d, lat, lon)),
             AOO_km2 = aoo$AOO_km2, occupiedCells = aoo$AOO_km2 / 4,
             AOO_pct = round(aoo$AOO_pct, 3), manuscript_pct = manuscriptAOO[[d$taxon[1]]])
})
print(aooTest)

# ANN: all datasets (fast). Table S4 values were the mean distance, in decimal degrees, to all
# neighbors within a distance band (not the nearest neighbor), then divided by 1000
manuscriptANN <- c(AMTH=9.991E-4, ARTH=1.122E-2, COGL=3.198E-6, HIWA=9.932E-6, MIGU=2.411E-3,
                   PICO=2.085E-3, QUAC=5.178E-7, QULO=1.593E-4, VILA=1.028E-3, YUBR=5.157E-6)
annTest <- purrr::map_dfr(pointsDataList, function(d){
  data.frame(taxon = d$taxon[1], uniqueLocations = nrow(dplyr::distinct(d, lat, lon)),
             ANN_km = signif(geo.calc.averageNearestNeighbor(geo.generateSpatialObject(d)), 4),
             manuscript = manuscriptANN[[d$taxon[1]]])
}) |>
  dplyr::mutate(rank = rank(ANN_km), manuscriptRank = rank(manuscript))
print(annTest)

# VOR: all datasets. Table S4 values were the mean Voronoi cell area in square degrees, which
# equals the (expanded) bounding box area / number of locations and so carries no information
# on the arrangement of the points. The metric is now the CV of the cell areas (see the function)
manuscriptVOR <- c(AMTH=1.568, ARTH=24.47, COGL=4.260E-6, HIWA=3.302E-5, MIGU=26.80,
                   PICO=1.143, QUAC=9.254E-3, QULO=1.121E-1, VILA=2.538, YUBR=1.450E-4)
vorTest <- purrr::map_dfr(pointsDataList, function(d){
  data.frame(taxon = d$taxon[1], uniqueLocations = nrow(dplyr::distinct(d, lat, lon)),
             VOR_cv = round(geo.calc.voronoiAreas(geo.generateSpatialObject(d)), 3),
             manuscript = manuscriptVOR[[d$taxon[1]]])
}) |>
  dplyr::mutate(rank = rank(VOR_cv), manuscriptRank = rank(manuscript))
print(vorTest)

# STD and ELG: all datasets. Table S4 STD was calculated on all individuals (now unique locations).
# ELA (pi * STD^2) and ELP (ellipse perimeter) rank the datasets identically to STD and are
# replaced by the elongation of the standard deviational ellipse
manuscriptSTD <- c(AMTH=128.85, ARTH=3081.23, COGL=1.51, HIWA=1.93, MIGU=1257.47,
                   PICO=402.62, QUAC=54.09, QULO=207.83, VILA=643.05, YUBR=6.23)
ellipseTest <- purrr::map_dfr(pointsDataList, function(d){
  g <- geo.generateSpatialObject(d)
  data.frame(taxon = d$taxon[1], uniqueLocations = nrow(dplyr::distinct(d, lat, lon)),
             STD_km = round(geo.calc.stdDistance(g), 2), manuscriptSTD = manuscriptSTD[[d$taxon[1]]],
             ELG = round(geo.calc.ellipseElongation(g), 2))
}) |>
  dplyr::mutate(rankSTD = rank(STD_km), manuscriptRankSTD = rank(manuscriptSTD), rankELG = rank(ELG))
print(ellipseTest)




















} # end of runMetricTests

# CALCULATING ALL POINT SUMMARIES ----
# Apply the wrapper to every dataset. AOO on the 2 km grid is slow for the wide ranging datasets
# (MIGU ~2.5 minutes, ARTH ~22 minutes), so the results are written to disk and read back on
# later runs; delete the file (or set overwriteSummaries to TRUE) to recalculate
overwriteSummaries <- FALSE
summariesFile <- paste0(GeoGenCorr_wd, 'Datasets/pointSummaryMeasures_carver.csv')
summariesList <- paste0(GeoGenCorr_wd, 'Datasets/pointSummariesList_carver.Rdata')

if(!file.exists(summariesFile) | overwriteSummaries){
  # One row per dataset, with the time taken (seconds) for reference
  pointSummaries <- purrr::map(pointsDataList, function(d){
    t0 <- Sys.time()
    ps <- geo.calc.pointSummaries(d)
    message(paste0(d$taxon[1], ' done (', round(as.numeric(difftime(Sys.time(), t0, units = 'secs'))), ' s)'))
    ps
  })
  # List of one row data.frames (the structure used by runPointSummaries.R)
  saveRDS(pointSummaries, file = summariesList)
  # Long table with the dataset as a column, for the CSV
  pointSummariesDF <- purrr::map_dfr(pointSummaries, ~.x, .id = 'taxon') |>
    dplyr::mutate(taxon = toupper(taxon))
  write_csv(x = pointSummariesDF, file = summariesFile)
} else {
  pointSummariesDF <- read_csv(summariesFile, show_col_types = FALSE)
  pointSummaries <- split(pointSummariesDF[, -1], tolower(pointSummariesDF$taxon))
}
print(pointSummariesDF)

# Matrix used by the correlation section: rows are the metrics, columns are the datasets
# (the percent version of AOO is not carried forward; it is a ratio of AOO to EOO)
pointSummariesMat <- pointSummariesDF |>
  dplyr::select(-AOO_pct) |>
  tibble::column_to_rownames('taxon') |>
  as.matrix() |>
  t()

# Comparison with Table S4 of the manuscript (Draft 5), values and ranks, to review with coauthors
tableS4 <- data.frame(
  taxon = c('AMTH','ARTH','COGL','HIWA','MIGU','PICO','QUAC','QULO','VILA','YUBR'),
  EOO = c(60767.70, 73591778.47, 10.05, 8.00, 6529925.52, 1308901.79, 3506.07, 135430.40, 88917.24, 186.90),
  AOO_pct = c(0.09017, 0.004025, 83.33, 85.57, 0.006113, 0.09157, 0.517, 0.5327, 0.1799, 38.81),
  ANN = c(9.991E-4, 1.122E-2, 3.198E-6, 9.932E-6, 2.411E-3, 2.085E-3, 5.178E-7, 1.593E-4, 1.028E-3, 5.157E-6),
  VOR = c(1.568, 24.47, 4.260E-6, 3.302E-5, 26.80, 1.143, 9.254E-3, 1.121E-1, 2.538, 1.450E-4),
  STD = c(128.85, 3081.23, 1.51, 1.93, 1257.47, 402.62, 54.09, 207.83, 643.05, 6.23),
  ELA = c(52120, 2.981E7, 7.3, 11.7, 4964000, 508900, 9186, 135600, 1298000, 121.7),
  ELP = c(861.5, 18117.73, 9.16, 11.32, 7916.43, 2432.1, 330.81, 1250.23, 3659.98, 38.65))
# Long format: one row per dataset and metric, with the new value, the manuscript value (where the
# metric existed), and the ranks under each
comparisonTable <- pointSummariesDF |>
  tidyr::pivot_longer(-taxon, names_to = 'metric', values_to = 'new') |>
  dplyr::left_join(tableS4 |> tidyr::pivot_longer(-taxon, names_to = 'metric', values_to = 'draft5'),
                   by = c('taxon', 'metric')) |>
  dplyr::group_by(metric) |>
  dplyr::mutate(newRank = rank(new), 
                draft5Rank = ifelse(is.na(draft5), NA, rank(draft5, na.last = 'keep'))) |>
  dplyr::ungroup() |>
  dplyr::arrange(metric, newRank)
write_csv(x = comparisonTable, file = paste0(GeoGenCorr_wd, 'Datasets/pointSummaryComparison_draft5.csv'))
print(comparisonTable, n = Inf)

# EXTRACTING OPTIMAL BUFFER SIZES ----
# To measure any possible correlation between optimal buffer sizes and the point summary statistics,
# the optimal buffer size for each dataset (and for each relevant coverage type) needs to be appended
# to the matrix of point summary values.

# Specify the filepaths to the resampling array for each dataset
resampArrList <- list(
  AMTH=paste0(GeoGenCorr_wd, 'Datasets/AMTH/resamplingData/AMTH_SMBO2_GE_5r_resampArr.Rdata'),
  ARTH=paste0(GeoGenCorr_wd, 'Datasets/ARTH/resamplingData/ARTH_SMBO2_GE_5r_resampArr.Rdata'),
  COGL=paste0(GeoGenCorr_wd, 'Datasets/COGL/resamplingData/COGL_SMBO2_GE_5r_resampArr.Rdata'),
  HIWA=paste0(GeoGenCorr_wd, 'Datasets/HIWA/resamplingData/HIWA_SMBO2_GE_5r_resampArr.Rdata'),
  MIGU=paste0(GeoGenCorr_wd, 'Datasets/MIGU/resamplingData/SMBO2_G2E/MIGU_SMBO2_G2E_5r_resampArr.Rdata'),
  PICO=paste0(GeoGenCorr_wd, 'Datasets/PICO/resamplingData/SMBO2_G2E/PICO_SMBO2_G2E_5r_resampArr.Rdata'),
  QUAC=paste0(GeoGenCorr_wd, 'Datasets/QUAC/resamplingData/QUAC_SMBO2_G2E_5r_resampArr.Rdata'),
  QULO=paste0(GeoGenCorr_wd, 'Datasets/QULO/resamplingData/SMBO2/QULO_SMBO2_G2E_5r_resampArr.Rdata'),
  VILA=paste0(GeoGenCorr_wd, 'Datasets/VILA/resamplingData/VILA_SMBO2_5r_resampArr.Rdata'),
  YUBR=paste0(GeoGenCorr_wd, 'Datasets/YUBR/resamplingData/YUBR_SMBO2_G2E_resampArr.Rdata')
)
# The resampling arrays are not tracked in the repository; stop here (the point summaries above
# are already written to disk) if they are not available on this machine
missingArrays <- unlist(resampArrList)[!file.exists(unlist(resampArrList))]
if(length(missingArrays) > 0){
  stop(paste0('Point summaries complete. Resampling arrays not found for the correlation section: ',
              paste(names(missingArrays), collapse = ', ')))
}
# Based on data in resampling arrays, extract the optimal buffer sizes for each species
optBuffs <- lapply(resampArrList, extractOptBuffs)
# For datasets without SDM values, add a column (in order to match dimensions with other datasets)
optBuffs$AMTH <- c(optBuffs$AMTH[[1]],NA,optBuffs$AMTH[[2]])
optBuffs$ARTH <- c(optBuffs$ARTH[[1]],NA,optBuffs$ARTH[[2]])
optBuffs$COGL <- c(optBuffs$COGL[[1]],NA,optBuffs$COGL[[2]])
optBuffs$HIWA <- c(optBuffs$HIWA[[1]],NA,optBuffs$HIWA[[2]])
optBuffs$VILA <- c(optBuffs$VILA[[1]],NA,optBuffs$VILA[[2]])
names(optBuffs$AMTH) <- names(optBuffs$ARTH) <- names(optBuffs$COGL)<- names(optBuffs$HIWA) <- 
  names(optBuffs$VILA) <- names(optBuffs$QULO)
# Convert the list of optimal buffer size values to a matrix
optBuffsMat <- matrix(unlist(optBuffs), ncol = length(optBuffs), byrow = FALSE)
colnames(optBuffsMat) <- names(optBuffs)
rownames(optBuffsMat) <- c('Opt_Geo-Buff', 'Opt_SDM-Buff', 'Opt_Eco-Buff')
# Combine the point summary matrix to the optimal buffer size matrix. Transpose such that 
# rows are datasets and columns are summary metrics, and order by optimal GeoBuff size
SMBO_Mat <- t(rbind(pointSummariesMat, optBuffsMat))
SMBO_Mat <- SMBO_Mat[order(SMBO_Mat[,'Opt_Geo-Buff']),]
# Create a separate matrix identical to the first, but only for species with SDMs
SMBO_SDM_Mat <- SMBO_Mat[-which(is.na(SMBO_Mat[,'Opt_SDM-Buff'])),]
SMBO_SDM_Mat <- SMBO_SDM_Mat[,-c(8,10)]
# Remove the row corresponding to SDM optimal buffer sizes from the original matrix
SMBO_Mat <- SMBO_Mat[,-9]

# BUILDING AND PLOTTING CORRELATION MATRICES ----
# Setting the correlation type to Spearman, since we don't know whether the relationship
# between the optimal buffer sizes and each point statistic is linear (and data likely
# isn't Normally distributed)
corType <- 'spearman'
# GEO/ECO COVERAGES
# Build a correlation matrix based off of values
corMat_SMBO <- rcorr(SMBO_Mat, type=corType)
# Replace NAs in diagonal of p-value matrix with 0s, to match dimensions
corMat_SMBO$P[which(is.na(corMat_SMBO$P))] <- 0
# Plot correlation matrix using corrplot. Label significant correlations using
# asterisks
corrplot(corMat_SMBO$r, type="upper", order="original", p.mat = corMat_SMBO$P, 
         sig.level = 0.01, insig = "label_sig", diag = FALSE)
mtext('Spearman correlation: Point stats and Geo/Eco coverages', side=3, line=1.2, adj=0.05, cex=1.2)

# SDM COVERAGES
# Build a correlation matrix based off of values
corMat_SMBO_SDM <- rcorr(SMBO_SDM_Mat, type=corType)
# Replace NAs in diagonal of p-value matrix with 0s, to match dimensions
corMat_SMBO_SDM$P[which(is.na(corMat_SMBO_SDM$P))] <- 0
# Plot correlation matrix using corrplot. Label significant correlations using
# asterisks
corrplot(corMat_SMBO_SDM$r, type="upper", order="original", p.mat = corMat_SMBO_SDM$P, 
         sig.level = 0.01, insig = "label_sig", diag = FALSE)
mtext('Spearman correlation: Point stats and Geo/SDM/Eco coverages', side=3, line=1.2, adj=0.05, cex=1.2)

# PLOTTING OPTIMAL BUFFER SIZES VERSUS POINTS BASED METRICS
plot(SMBO_Mat[,'ANN'], SMBO_Mat[,'Opt_Geo-Buff'], pch=16, col='black',
     ylab='Optimal Geographic Buffer Sizes', xlab='Average Nearest Neighbor Values',
     main='SMBO2: Buffer sizes across ANN metrics', ylim=c(-10,600))

plot(SMBO_Mat[,'EOO'], SMBO_Mat[,'Opt_Geo-Buff'], pch=16, col='black',
     ylab='Optimal Geographic Buffer Sizes', xlab='Extent of Occurrence Values',
     main='SMBO2: Buffer sizes across EOO metrics')

plot(SMBO_Mat[,'StDevEllP'], SMBO_Mat[,'Opt_Geo-Buff'], pch=16, col='black',
     ylab='Optimal Geographic Buffer Sizes', xlab='Standard Deviation Ellipses Perimeter Values',
     main='SMBO2: Buffer sizes across St. Dev. Ellipse Perimeter metrics')
