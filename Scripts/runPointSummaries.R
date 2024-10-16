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
## should replace with grelp indexing on the file path, as this hard code index will not work for long... 
dataList <- list(
  amth = read_csv(data[1]) |> 
    dplyr::mutate(taxon = "AMTH")|>
    dplyr::select(lat = lat, lon = long), 
  migu = read_csv(data[4]) |> 
    dplyr::mutate(taxon = "MIGU")|>
    dplyr::select(lat = Latitude , lon = Longitude), # maybe just North America? 
  pico = read_csv(data[5]) |> 
    dplyr::mutate(taxon = "PICO")|>
    dplyr::select(lat = Latitude , lon = Longitude),
  quac = read_csv(data[6])|> 
    dplyr::mutate(taxon = "QUAC")|>
    dplyr::select(lat = decimalLatitude  , lon = decimalLongitude), # not sure what the difference is? 
  qulo = read_csv(data[8])|> 
    dplyr::mutate(taxon = "QULO")|>
    dplyr::select( lat = decimalLatitude  , lon = decimalLongitude),
  yubr = read_csv(data[10])|> 
    dplyr::mutate(taxon = "YUBR")|>
    dplyr::select(lat = decimalLatitude  , lon = decimalLongitude)
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
