###
# quick map to show the spatial variability in the subset point locations 
#
#
###

pacman::p_load(leaflet, terra,googledrive,readr, RColorBrewer )
setwd("~/trueNAS/work/GeographicGeneticCorrelation")
# only need if your doing the google drive thing
# drive_auth() 
# read in MIGU data 
path <- paste0(getwd(), "/Datasets/MIGU")
# files 
files <- list.files(path = path, 
                    recursive = TRUE,
                    full.names = TRUE)
# grab the CSV for the coord and the distribution 
p1 <- read_csv(files[9])
r1 <- terra::rast(files[7])

# generate dataset from the subset elements 
file_id <- "https://docs.google.com/spreadsheets/d/1FWc2bEcOUPZ9s6bQ6NYCjuLgKZcBy6VhGpmLaZ6g9Jg/edit?gid=110657149#gid=110657149"
temp <- drive_download(as_id(file_id), type = "csv", overwrite = TRUE) 
d1 <- read.csv(temp$local_path)

# generate the spatial object 
points <- vect(p1, geom = c("Longitude","Latitude"), crs = "EPSG:4326") 


# produce subsets of the points 
s1 <- points[points$`Sample Name` %in% d1$Subset1, ]
s2 <- points[points$`Sample Name` %in% d1$Subset2, ]
s3 <- points[points$`Sample Name` %in% d1$Subset3, ]
s4 <- points[points$`Sample Name` %in% d1$Subset4, ]
s5 <- points[points$`Sample Name` %in% d1$Subset5, ]
s6 <- points[points$`Sample Name` %in% d1$Subset6, ]
s7 <- points[points$`Sample Name` %in% d1$Subset7, ]

# generate function to produce the map 
## buffer distance in meter 
bufferDistList <- c(5000,10000,50000)
# testing parameters 
points <- s1
raster <- r1
title <- "group1"
mapBuffers <- function(points, raster, title, bufferDistList){
  # Create buffers
  buffers <- lapply(bufferDistList, function(dist) {
    terra::buffer(points, width = dist)
  })
  
  # set palette 
  buffer_colors <- RColorBrewer::brewer.pal(n = length(bufferDistList), name = "Blues")
  
  # define control group names 
  controlGroups <- paste0("buffer dist: ", bufferDistList)
  
  # generate map with the point features 
  map <- leaflet() |>
    addTiles() |>  # Add a basemap
    addRasterImage(
      x = r1,
      colors = c("transparent", "#d8b365"),  # Set 0 to transparent, 1 to blue
      opacity = 0.8,
      group = "Distribution"
    )|>
    addControl(
      html = paste0("<h3>", title, "</h3>"), # Your title HTML
      position = "topleft"  # Position of the title
    ) 
    
    
  for(i in rev(seq_along(bufferDistList))){
    map <- map |>
      addPolygons(
        data = buffers[[i]],
        color = buffer_colors[i],
        fillOpacity = 0.4,
        group = controlGroups[i]
      )
  }
  # add control groups 
  map <- map |>
    # Add the points
    addCircleMarkers(
      data = points,
      radius = 5,
      color = "Black",
      stroke = FALSE,
      fillOpacity = ,
      group = "Points",
      popup = points$`Sample Name`)|>
    addLayersControl(
      overlayGroups = c("Points","Distribution", controlGroups),
      options = layersControlOptions(collapsed = FALSE)
    )|>
    addLegend(
      position = "bottomright",
      colors = buffer_colors,
      labels = controlGroups
    )
  
  return(map)
}

# single feature 
map <- mapBuffers(points = s1, title = "group1", raster = r1, bufferDistList = bufferDistList)

# name of the map 
titles <- paste0("Group ", 1:7)
# mutliple features 
maps <- purrr::map2(.x = c(s1,s2,s3,s4,s5,s6,s7), .y = titles, .f = mapBuffers,
            raster = r1, bufferDistList = bufferDistList )
