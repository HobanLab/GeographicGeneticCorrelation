###
# generalized method for determining the least number of locations require to get maximum coverage of a sdm 
###
# environment 
pacman::p_load(terra, dplyr,readr)

# read in datasets 
r1 <- terra::rast("Datasets/MIGU/Geographic/MIGU_255inds_rast_Carver.tif")
# point data 
p1 <- read_csv("Datasets/MIGU/Geographic/MIGU_coordinates.csv")

## transform the datasets 
# convert raster to vector 
r2 <- terra::as.polygons(r1)
# select the predicted presence area
r2 <- r2[r2$Threshold == 1, ]

# produce a point object
s1 <- terra::vect(p1, geom = c("Longitude", "Latitude"), crs = r1)

# buffer all locations 
# Buffer the points by 50km
bufferRadius <- 50000  # in meters
buffers <- buffer(s1, width = bufferRadius)

# get the original area 
buffers$fullArea <-expanse(buffers, unit = "km")

# crop buffers to model area
croppedBuffers <- crop(buffers, r2)
# new area 
croppedBuffers$cropArea <- expanse(croppedBuffers, unit = "km") 

# generate a dataframe for holding change in area results 
croppedBuffers$equalArea <- croppedBuffers$fullArea - croppedBuffers$cropArea
