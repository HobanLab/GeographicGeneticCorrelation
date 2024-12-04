# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GEN-GEO-ECO CORRELATION: POINT SUMMARY METRICS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script calculates geographic metrics based on the collection of coordinate points from 
# each dataset used in the geographic-genetic correlation analyses. The goal of this analysis 
# is to examine any possible trends between these point-based metrics and the optimal geographic
# buffer sizes.

# The key geo.calc.pointSummaries function used in this script is declared in the functions_GeoGenCoverage.R script. 

# Load packages 
pacman::p_load(redlistr, spatialEco, sf, terra, tmap, dplyr, sfdep, sfheaders,
               deldir, readr)
tmap_mode("view")

# Read in relevant functions
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/Code/'
setwd(GeoGenCorr_wd)
source('Scripts/functions_GeoGenCoverage.R')

# Create a list of the relevant CSVs in the Datasets repository 
pointsData <- list.files(paste0(getwd(), "/Datasets"), pattern = ".csv", recursive = TRUE, full.names = TRUE)
pointsData <- pointsData[grep('Geographic', pointsData)]
# Select one coordinate file per species
amth <- pointsData[grepl(pattern ="AMTH", pointsData)]
arth <- pointsData[grepl(pattern ="ARTH", pointsData)]
cogl <- pointsData[grepl(pattern ="COGL", pointsData)]
hiwa <- pointsData[grepl(pattern ="HIWA", pointsData)]
migu <- pointsData[grepl(pattern ="MIGU", pointsData)][2] # Necessary because of a 2nd file of global coordinates
pico <- pointsData[grepl(pattern ="PICO", pointsData)]
quac <- pointsData[grepl(pattern ="QUAC", pointsData)]
qulo <- pointsData[grepl(pattern ="QULO", pointsData)]
yubr <- pointsData[grepl(pattern ="YUBR", pointsData)][2] # Necessary because of a 2nd (original) file of coordinates

# Build a list of the relevant files, and fields within those files, for each dataset.
# Note that the column names between different coordinate files differ
pointsDataList <- list(
  amth = read_csv(amth) |>
    dplyr::mutate(taxon = "AMTH")|>
    dplyr::select(taxon, lat = lat, lon = long),
  arth = read_csv(arth,col_names = FALSE) |>
      plyr::mutate(taxon = "ARTH")|>
      dplyr::select(taxon, lat = X6, lon = X7),
  cogl = read_csv(cogl) |> 
    dplyr::mutate(taxon = "COGL")|>
    dplyr::select(taxon, lat = Latitude , lon = Longitude),
  hiwa = read_csv(hiwa)|> 
    dplyr::mutate(taxon = "HIWA")|>
    dplyr::select(taxon, lat = decimalLatitude  , lon = decimalLongitude),
  migu = read_csv(migu) |> 
    dplyr::mutate(taxon = "MIGU")|>
    dplyr::select(taxon, lat = Latitude , lon = Longitude), 
  pico = read_csv(pico) |> 
    dplyr::mutate(taxon = "PICO")|>
    dplyr::select(taxon, lat = Latitude , lon = Longitude),
  quac = read_csv(quac)|> 
    dplyr::mutate(taxon = "QUAC")|>
    dplyr::select(taxon, lat = decimalLatitude  , lon = decimalLongitude), 
  qulo = read_csv(qulo)|> 
    dplyr::mutate(taxon = "QULO")|>
    dplyr::select(taxon, lat = decimalLatitude  , lon = decimalLongitude),
  yubr = read_csv(yubr)|>
    dplyr::mutate(taxon = "YUBR")|>
    dplyr::select(taxon, lat = decimalLatitude  , lon = decimalLongitude)
)

# Apply function which calculates multiple point summary metrics to list of datasets
pointSummaries <- lapply(pointsDataList, geo.calc.pointSummaries)
# Transform the point summary values into a list, where columns are the species and rows are the metrics
pointSummariesMat <- matrix(unlist(pointSummaries), ncol = length(pointSummaries), byrow = FALSE)
colnames(pointSummariesMat) <- toupper(names(pointSummaries))
rownames(pointSummariesMat) <- c("EOO", "AOO", "Average Nearest Neighbor", "Average Voroni Area",
                                 "Standard Distance", "Standard Distance Ellipse Area",
                                 "Standard Deviation Ellipse Perimeter")
# Write the matrix of point values to disk
write_csv(x = pointSummaries, file = paste0(GeoGenCorr_wd,"Datasets/pointSummaryMeasures.csv"))
