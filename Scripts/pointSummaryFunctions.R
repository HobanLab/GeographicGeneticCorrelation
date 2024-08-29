###
# a grouping of all point summary calculations 
# 
###

pacman::p_load(redlistr, spatialEco, sf, terra, tmap, dplyr, sfdep,sfheaders,
               deldir)
tmap_mode("view")

# reference 
degreeToMeter <- 111139
degreeToMeterWithLat <- function(degrees, latitude){
  40075000*cos(latitude)/360
}




wd <- "~/Documents/GeographicGeneticCorrelation"
tempData <- read.csv(paste0(wd, "/Datasets/QULO/Geographic/Quercus_lobata.csv"))
# genereate spatial Data 
## lat long data
sp1 <- sf::st_as_sf(x = tempData,
                    coords = c("decimalLongitude","decimalLatitude"),
                    crs = CRS("+proj=longlat +datum=WGS84 +no_defs +type=crs"))
## projected data #alberts equal area 8857
sp1_proj <- sf::st_transform(x = sp1, crs ="+proj=moll")

vect1 <- terra::vect(x = tempData, 
                     geom=c("decimalLongitude","decimalLatitude"), 
                     crs= crs("+proj=longlat +datum=WGS84 +no_defs +type=crs"))


# redlistR  ---------------------------------------------------------------
## redlistr
### these functions seem to require a raster object, which is weird so I'm rewriting them 
# EOO 
## need to convert to a multipoint object
bb <- st_convex_hull(st_union(sp1_proj))
## export value 
eooArea <- st_area(bb)


# AOO ---------------------------------------------------------------------
## should really be using a distribution with this so we might want to call it something difference
# created the gridded feature 
allGrids <- st_make_grid(bb,
                         square = T, 
                         cellsize = 2000)
# test for intersetion with the conver hull 
intersectionCells <- sf::st_intersects(x = allGrids, y = bb, sparse= FALSE) # the grid, covering bounding box
# filter the full polygon feature
selectAreas <- allGrids[intersectionCells]
# Determine the number of areas that have an observation 
intersectionPoints <- sf::st_intersects(x = selectAreas, y = sp1_proj, sparse= TRUE)
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


# reference spatial feature -----------------------------------------------
refSP <- tempData |>
  sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
               crs = CRS("epsg:4326"),
               remove = FALSE)
## break though with sfdep package 


# meanCenter --------------------------------------------------------------
meanCenter <- sfdep::center_mean(geometry = refSP)
meanCoords <- as.data.frame(sf::st_coordinates(meanCenter))
medianCenter <- sfdep::center_median(geometry = refSP)
euclideanMedian <- sfdep::euclidean_median(geometry = refSP)


# Standard Distance  ------------------------------------------------------
stdDist <- std_distance(geometry = refSP)
## produces a list of points 
stdDistEllipse <- ellipse(x = meanCoords$X,y = meanCoords$Y , sx = stdDist, sy = stdDist)
poly <- stdDistEllipse |>
  as.data.frame()|>
  sfheaders::sf_polygon(x = "x", y = "y")
sf::st_crs(poly) <- 4326


# standard deviation ellipse  ---------------------------------------------
stdDevElli <- std_dev_ellipse(geometry = refSP)
## points for the ellipse 
stdDevEllipse <- st_ellipse(geometry =stdDevElli,
                            sx = stdDevElli$sx,
                            sy = stdDevElli$sy,
                            rotation = -stdDevElli$theta
                            )



# global density
# Points/Area

# kernal density 
# https://cran.r-project.org/web/packages/eks/vignettes/tidysf_kde.html
## good options for producing this dataset, just not 100% sure how to translate 
## the produce to a comparable number 

# average nearest neighbor
aveNearNeighbor <- refSP |> 
  mutate(
    nb = st_dist_band(geometry),
    dists = st_nb_dists(geometry, nb,longlat = TRUE),
    avg_dist = purrr::map_dbl(dists, mean),
    group_average = mean(avg_dist)
  ) 


## library unclear 
# voronoi diagram 

tesselation <- deldir(refSP$decimalLongitude, refSP$decimalLatitude)
tiles <- tile.list(tesselation)
# get the average area of all the tiles 
## this is basically a unitless measure that is valuable for comparision only 
selectArea <- function(tile){
  tile[6]
  }
aveVoroniArea <- purrr::map(.x = tiles,.f = selectArea ) |> unlist() |> mean(na.rm=TRUE)

# comparison against random 

