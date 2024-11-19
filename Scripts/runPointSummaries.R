###
# a script to render the point summaries of all available datasets
#
#
### 

# packages 
pacman::p_load(redlistr, spatialEco, sf, terra, tmap, dplyr, sfdep,sfheaders,
               deldir, readr)
tmap_mode("view")

# environment
setwd("~/Documents/GeographicGeneticCorrelation")
source(paste0(getwd(),"/Scripts/pointSummaryFunctions.R"))

# observation data 
data <- list.files(paste0(getwd(), "/Datasets"),
                   pattern = ".csv",
                   recursive = TRUE,
                   full.names = TRUE)
# manually selecting one file per species and standardizing 
amth <- data[grepl(pattern ="AMTH", data)]
arth <- data[grepl(pattern ="ARTH", data)]
migu <- data[grepl(pattern ="MIGU", data)][2]
pico <- data[grepl(pattern ="PICO", data)]
quac <- data[grepl(pattern ="QUAC", data)][2]
qulo <- data[grepl(pattern ="QULO", data)]
yubr <- data[grepl(pattern ="YUBR", data)][2]

dataList <- list(
  # not present in the data folder? 
  # amth = read_csv(amth) |> 
  #   dplyr::mutate(taxon = "AMTH")|>
  #   dplyr::select(taxon, lat = lat, lon = long), 
  arth = read_csv(arth,col_names = FALSE) |>
      plyr::mutate(taxon = "ARTH")|>
      dplyr::select(taxon, lat = X6, lon = X7),
  migu = read_csv(migu) |> 
    dplyr::mutate(taxon = "MIGU")|>
    dplyr::select(taxon, lat = Latitude , lon = Longitude), # maybe just North America? 
  pico = read_csv(pico) |> 
    dplyr::mutate(taxon = "PICO")|>
    dplyr::select(taxon, lat = Latitude , lon = Longitude),
  quac = read_csv(quac)|> 
    dplyr::mutate(taxon = "QUAC")|>
    dplyr::select(taxon, lat = decimalLatitude  , lon = decimalLongitude), # not sure what the difference is? 
  qulo = read_csv(qulo)|> 
    dplyr::mutate(taxon = "QULO")|>
    dplyr::select(taxon, lat = decimalLatitude  , lon = decimalLongitude)
  # some weird stuff going on with this dataset? 
  # yubr = read_csv(yubr)|> 
  #   dplyr::mutate(taxon = "YUBR")|>
  #   dplyr::select(taxon, lat = decimalLatitude  , lon = decimalLongitude)
)


# create spatial objects 
spList <- purrr::map(.x = dataList, .f = generateSpatailObject)

# eoo area -- square meters
eooList <- purrr::map(.x = spList, .f = boundingBoxAreas)|>
  as.data.frame()|>
  dplyr::mutate(across(everything(), as.numeric))

# aoo val -- takes a bit to run 
aooList<- purrr::map(.x = spList, .f = aooLikeMeasures)|>
  as.data.frame()

# average nearest neighbor - meters 
annList <- purrr::map(.x = spList, .f = averageNearestNeighbor)|>
  as.data.frame()

# average voroni areas -- unitless 
vorList <- purrr::map(.x = spList, .f = voriniAreas)|>
  as.data.frame()

# standard distance - meters 
standardDistanceList <- purrr::map(.x = spList, .f = standardDistance)|>
  as.data.frame()

# standard distance ellipse area in meters 
standardDistanceEllispeAreaList <- purrr::map(.x = spList, .f = standardDistanceEllipseArea)|>
  as.data.frame()|>
  dplyr::mutate(across(everything(), as.numeric))


# standard deviation ellipse length in meteres  
standardDeviationEllispeLengthList <- purrr::map(.x = spList, .f = stdDevationEllipseArea) |>
  as.data.frame()|>
  dplyr::mutate(across(everything(), as.numeric))
  

# store results into a dataframe 
df <- eooList |>
  bind_rows(aooList)|>
  bind_rows(annList)|>
  bind_rows(vorList)|>
  bind_rows(standardDistanceList)|>
  bind_rows(standardDistanceEllispeAreaList)|>
  bind_rows(standardDeviationEllispeLengthList)

df$measures <- c("EOO", "AOO", "Average Nearest Neighbor", "Average Voroni Area",
                 "Standard Distance", "Area of the Standard Distance Ellipse",
                 "Perimeter lenght of the Standard Deviation Ellipse")

df2 <- df |>
  dplyr::select(measures, everything())

# export
write_csv(x = df2, file = paste0(getwd(),"/outputs/pointSummaryMeasures.csv"))
