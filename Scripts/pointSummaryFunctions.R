###
# a grouping of all point summary calculations 
# 
###

degreeToMeterWithLat <- function(degrees, latitude){
  40075000*cos(latitude)/360
}

generateSpatailObject <- function(data){
  ## clean 
  data1 <- data |>
    dplyr::filter(!is.na(lon))|>
    dplyr::filter(!is.na(lat))
  
  ## lat long data
  sp1 <- sf::st_as_sf(x = data1,
                      coords = c("lon","lat"),
                      crs = CRS("+proj=longlat +datum=WGS84 +no_defs +type=crs"),
                      remove = FALSE)
  ## projected data #alberts equal area 8857
  sp1_proj <- sf::st_transform(x = sp1, crs ="+proj=moll")
  return(sp1_proj)
}

# redlistR  ---------------------------------------------------------------
## EOO ----
boundingBoxAreas <- function(data){
  ## convert to a multipoint object
  bb <- st_convex_hull(st_union(data))
  ## export value 
  eooArea <- st_area(bb,)
  return(eooArea)
}


## AOO ---------------------------------------------------------------------
aooLikeMeasures <- function(data){
  bb <- st_convex_hull(st_union(data))
  # created the gridded feature 
  allGrids <- st_make_grid(bb,
                           square = T, 
                           cellsize = 2000)
  # test for intersetion with the conver hull 
  intersectionCells <- sf::st_intersects(x = allGrids,
                                         y = bb,
                                         sparse= FALSE) # the grid, covering bounding box
  # filter the full polygon feature
  selectAreas <- allGrids[intersectionCells]
  # Determine the number of areas that have an observation 
  intersectionPoints <- sf::st_intersects(x = selectAreas,
                                          y = data,
                                          sparse= TRUE)
  # count of all features with at least one observation 
  testIntersect <- function(area){
    if(length(area)>0){
      return(1)
    }else{
      0
    }
  }
  areasWithPoints <- purrr::map(.x = intersectionPoints, .f = testIntersect )|> 
    unlist()|> 
    sum()
  
  # calculate A00 --- 
  aoo <- areasWithPoints / nrow(intersectionPoints) *100
  return(aoo)
}


# average nearest neighbor ------------------------------------------------

averageNearestNeighbor <- function(data){
  # convert back to lat lon 
  refSP <- data |>
    sf::st_transform(crs = 4326)
  # remove duplicated locations 
  refSP_noDup <- refSP[!duplicated(refSP),]
  # average nearest neighbor
  aveNearNeighbor <- refSP_noDup |> 
    mutate(
      nb = st_dist_band(geometry),
      dists = st_nb_dists(geometry, nb,longlat = TRUE),
      avg_dist = purrr::map_dbl(dists, mean),
      group_average = mean(avg_dist, na.rm = TRUE)
    ) 
  return(aveNearNeighbor$group_average[1])
}


# vironi diagrams  --------------------------------------------------------

voriniAreas <- function(data){
  # produce the tesselations 
  tesselation <- deldir(data$lon, data$lat)
  # creates the spatial representation on the areas 
  tiles <- tile.list(tesselation)
  ## helper for indexing out of the tiles object
  selectArea <- function(tile){
    tile[6]
  }
  # go through each tile and determine the overall average area
  aveVoroniArea <- purrr::map(.x = tiles,.f = selectArea ) |> unlist() |> mean(na.rm=TRUE)
  return(aveVoroniArea)
}


# Standard Distance  ------------------------------------------------------
standardDistance <- function(data){
  stdDist <- std_distance(geometry = data)
  return(stdDist)
}

# standard deviation ellipse area 

standardDistanceEllipseArea <- function(data){
  
  # determine the standard distance 
  stdDist <- std_distance(geometry = data)
  
  # determine the mean center 
  meanCenter <- sfdep::center_mean(geometry = data)
  
  # grap coords from the mean center
  meanCoords <- as.data.frame(sf::st_coordinates(meanCenter))
  
  ## produces a list of points within the ellipse 
  stdDistEllipse <- ellipse(x = meanCoords$X,y = meanCoords$Y , sx = stdDist, sy = stdDist)
  # convert to a ploygon and project 
  poly <- stdDistEllipse |>
    as.data.frame()|>
    sfheaders::sf_polygon(x = "x", y = "y")
  sf::st_crs(poly) <- crs(data)
  # make valid and reproject 
  validPoly <- st_make_valid(poly)
    
  # determine area of the polygon 
  stdEllipseArea <- st_area(poly)
  return(stdEllipseArea)
  
}

stdDevationEllipseArea <- function(data){
  stdDevElli <- std_dev_ellipse(geometry = data)
  ## points for the ellipse 
  stdDevEllipse <- st_ellipse(geometry =stdDevElli,
                              sx = stdDevElli$sx,
                              sy = stdDevElli$sy,
                              rotation = -stdDevElli$theta)
  # determine length of the perimater  
  stdDevationEllipseArea <- st_length(stdDevEllipse)
  return(stdDevationEllipseArea)
  }

