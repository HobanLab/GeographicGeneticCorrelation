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

# crop all points buffer to state 
km5 <- bufferCrop(p1 = p1, dist = 5000, crop = states)
km25 <-  bufferCrop(p1 = p1, dist = 25000, crop = states)
km50 <-  bufferCrop(p1 = p1, dist = 50000, crop = states)
km200 <-  bufferCrop(p1 = p1, dist = 200000, crop = states)

# crop to selected points to state cropped buffer 
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
                fillOpacity = 0.8
                # stroke = NA
                )
}
addPoint <- function(map, object, color){
  
}



mb <- m |> 
  # addBuffer(object = km200, color = "Green")|>
  addBuffer(object = km50, color = "#31a354")
  # addBuffer(object = km25, color = "Green")|>
  # addBuffer(object = km5, color = "Green") |>


# add the selected areas buffers 
mba <- mb |>
  addBuffer(object = s50, color =  "#54278f")|>
  addBuffer(object = s25, color = "#bcbddc")|>
  addBuffer(object = s5, color = "#dadaeb")

mba |>   
  addCircles(data = p1, 
             color =  "#c7e9b4",
             fillColor = "black",
             opacity = 0.5,
             stroke = 0.5,
             weight = 3)|>
  addCircles(data = selected, color =  "#fe9929") 
  


# SDM method  -------------------------------------------------------------
## crop buffers to the sdm 
poly1 <- as.polygons(r1, aggregate=TRUE, values=TRUE, na.rm=TRUE)
poly2 <- poly1[poly1$Threshold == 1, ]

# crop all points buffer to sdm  
r5 <- bufferCrop(p1 = p1, dist = 5000, crop = poly2)
r25 <-  bufferCrop(p1 = p1, dist = 25000, crop = poly2)
r50 <-  bufferCrop(p1 = p1, dist = 50000, crop = poly2)
r200 <-  bufferCrop(p1 = p1, dist = 200000, crop = poly2)

# crop to selected points to state cropped buffer 
sr5 <- bufferCrop(p1 = selected, dist = 5000, crop = r5)
sr25 <- bufferCrop(p1 = selected, dist = 25000, crop = r25)
sr50 <- bufferCrop(p1 = selected, dist = 50000, crop = r50)
sr200 <- bufferCrop(p1 = selected, dist = 200000, crop = r200)


mb <- m |> 
  addPolygons(data = r50,
              color = "#31a354",
              stroke = NA,
              opacity = 1,
              fillOpacity = 1
  )
  # addBuffer(object = km200, color = "Green")|>
  # addBuffer(object = r50, color = "Green")
# addBuffer(object = km25, color = "Green")|>
# addBuffer(object = km5, color = "Green") |>


# add the selected areas buffers 
mba<- mb |>
  addBuffer(object = sr50, color = "#54278f")|>
  addBuffer(object = sr25, color = "#bcbddc")|>
  addBuffer(object = sr5, color = "#dadaeb")

mba |>   
  addCircles(data = p1, 
             # color = "#d0d1e6",
             stroke = NA,
             radius = 1,
             fillColor = "#d0d1e6") |>
  addCircles(data = selected, 
             color = "#045a8d",
             radius = 1)

### todo 
# both buffer and the sdm map 
# all points 
## sample are bold 
## others are light 
## match color scheme 








