##
# figure for publicaiton 

pacman::p_load(sf, dplyr, tmap,terra, leaflet, tigris)
tmap_mode("view")

# 1. Crop the buffers on the right side, such that buffers don't extend past the SDM area,
# 2. Include a representation of other buffer sizes. Maybe these could just be 
# wider circles of different colors
# 3. Values on the bottom, representing the actual geographic coverage values 
# for different approaches/buffer sizes. I think this point is less crucial, but could be nice.



# use QULO
r1 <- terra::rast("Datasets/QULO/Geographic/QULO_436inds_rast_Carver.tif")
# subset distribution 
r2 <- terra::subst(r1, 0, NA)
## querqus lobata 
p1 <- read.csv("Datasets/QULO/Geographic/QULO_coordinates.csv") |>
  dplyr::filter(!is.na(decimalLongitude))|>
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"))|>
  terra::vect()|>
  terra::set.crs(terra::crs(r1))
# select 5 random locations for selected points 
selected <- p1[sample(nrow(p1), 5), ]


# get country data to crop 
states <- tigris::nation() |>
  terra::vect() |> 
  terra::project(terra::crs(r1))


# produce buffer objects 
bufferCrop <- function(p1, dist, crop){
  terra::buffer(p1, dist)|>
    terra::aggregate()|> 
    terra::crop(crop)
}


km5 <- bufferCrop(p1 = p1, dist = 5000, crop = states)
km25 <-  bufferCrop(p1 = p1, dist = 25000, crop = states)
km50 <-  bufferCrop(p1 = p1, dist = 50000, crop = states)
km200 <-  bufferCrop(p1 = p1, dist = 200000, crop = states)
# 
s5 <- bufferCrop(p1 = selected, dist = 5000, crop = km5)
s25 <- bufferCrop(p1 = selected, dist = 25000, crop = km25)
s50 <- bufferCrop(p1 = selected, dist = 50000, crop = km50)
s200 <- bufferCrop(p1 = selected, dist = 200000, crop = km200)


# initialize map area map 
m <- leaflet() |>
  addTiles() 

# try buffer options first 
addBuffer <- function(map, object, color){
  map2 <- map |>
    addPolygons(data = object,
                color = color,
                )
}
addPoint <- function(map, object, color){
  
}



mb <- m |> 
  # addBuffer(object = km200, color = "Green")|>
  addBuffer(object = km50, color = "Green")
  # addBuffer(object = km25, color = "Green")|>
  # addBuffer(object = km5, color = "Green") |>


# add the selected areas buffers 
mba <- mb |>
  addBuffer(object = s50, color = "purple")|>
  addBuffer(object = s25, color = "purple")|>
  addBuffer(object = s5, color = "purple")

mba |>   addCircles(data = p1)
  