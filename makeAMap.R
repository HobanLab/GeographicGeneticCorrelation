
#### note some required libraries as these might not be loaded in the workflow. 
# pacman::p_load(leaflet,sf,terra,tmap,dplyr,raster)

# tmap_mode("view")
# 
# # testDataSets 
# raster <- terra::rast("~/Documents/GeographicGeneticCorrelation/Datasets/QUAC/Geographic/prj_threshold.tif")
# points <- read.csv("~/Documents/GeographicGeneticCorrelation/Datasets/QUAC/Geographic/QUAC_coord_ind.csv")|>
#   sf::st_as_sf(coords = c("decimalLongitude","decimalLatitude"),
#               crs= terra::crs(raster),
#               remove = FALSE) # Keeps the lat lon values in the dataframe. defaults to TRUE
# buffer <- terra::buffer(x = vect(points), width = 5000) # width in meters 
# 
# # quick look. 
# qtm(points)
# qtm(raster)
# 
# # test function 
# ## just points and rasters 
# m1 <- makeAMap(points = points, raster=raster)
# ## with buffer 
# m2 <- makeAMap(points = points, raster=raster, buffer = buffer)


#' make a map
#'
#' @param points : sf point object  
#' @param raster : terra raster object... converted to raster within the function. 
#' @param buffer : optional sf point buffered feature. 
#'
#' @return one of two maps depending on if a buffer input object was defined or not.
makeAMap <- function(points,raster,buffer=NA){
  
  #centroid
  ## obnoxious number of steps here but wil help with genealize this process... 
  centroid <- points |>
    dplyr::mutate(group = 1)|>
    group_by(group) |>
    summarize(geometry = st_union(geometry)) |>
    st_centroid()|>
    st_coordinates()
  
  # define base map 
  map1 <- leaflet(options = leafletOptions(minZoom = 4)) |>
    #set zoom levels
    setView(lng = centroid[1]
            , lat = centroid[2]
            , zoom = 6) |>
    ## tile providers ----------------------------------------------------------
    addProviderTiles("OpenStreetMap", group = "OpenStreetMap") |>
    ## add point features -----------------------------------------------------
    addCircleMarkers(
      data = points,
      color = "#4287f5",
      radius = 0.2,
      group = "Points",
      # add highlight options to make labels a bit more intuitive 
    ) |> 
    addRasterImage(
      x = raster::raster(raster),
      colors = c("#ffffff95", "#8bed80"),
      group = "Raster"
    )|>
    addLegend(
      position = "topright",
      colors = c("#ffffff95", "#8bed80"),
      labels = c("potential area", "predicted area"),
      group = "Raster"
    )|>
    addLayersControl(
      #baseGroups = c("OpenStreetMap", "Satellite"),
      overlayGroups = c(
        "Points",
        "Raster"
      ),
      position = "bottomleft",
      options = layersControlOptions(collapsed = FALSE),
    ) 
  
  # add buffers if those are wanted. 
  if(!is.na(buffer)){
    buffs <- st_as_sf(buffer)
    # new map 
    map2 <- map1 |>
      addPolygons(
        data = buffs,
        weight = 1,
        fill = FALSE,
        opacity = 0.6,
        color = "#181b54",
        dashArray = "3",
        group = "Buffer"
      )|>
  addLayersControl(
    #baseGroups = c("OpenStreetMap", "Satellite"),
    overlayGroups = c(
      "Buffer",
      "Points",
      "Raster"
    ),
    position = "bottomleft",
    options = layersControlOptions(collapsed = FALSE),
  ) 
  map <- map2
  
  }else{
    map <- map1 
  }
  
  return(map)
}
