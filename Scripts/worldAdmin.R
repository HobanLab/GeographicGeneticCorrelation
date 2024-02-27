###
# helper function to pull input world data if the content does not exist
# carverd@colostate.edu
# 2023-12-22
### 


#' grabWorldAdmin
#'
#' @param GeoGenCorr_wd : current working directory define by either variable or getwd()
#' @param fileExtentsion : choice between shp and gpkg.Prefer gpkg for single file. 
#' @param overwrite : defaults to false; True will force a redownload of the file.
#'
#' @description
#' Checks to see if the file exists. If true it is loaded. If false it is downloaded.
#' Allow uses to specify the file extention. .shp or .gpkg are options
#' 
#'
#' @return terra vect object of the world admin layer
#' 
grabWorldAdmin <- function(GeoGenCorr_wd, fileExtentsion, overwrite=FALSE){
  path <- file.path(paste0(GeoGenCorr_wd,
                           'GIS_shpFiles/world_countries_10m/world_countries_10m',
                           fileExtentsion))
  # test for presence of file. 
  if(file.exists(path) | overwrite == TRUE){
    # read in for return 
    world_poly_clip <- terra::vect(path)
  }else{
    # download data from natural earth 
    download <- rnaturalearth::ne_countries(scale = 10,
                                            type = "countries",
                                            returnclass = "sf")
    # convert to terra vect and simply 
    # world_poly_clip <- download |>
    #   terra::vect() |>
    #   #this select is option but I don't think we need much of any attribute data on this file. 
    #   dplyr::select("admin", "adm0_a3", "name","name_long" )
    
    # ACK: the dplyr::select command errors for me 
    # (no applicable method for 'select' applied to an object of class "SpatVector")
    # Command seems optional anyway, so commented out for now
    world_poly_clip <- download |>
      terra::vect()
    
    # export the file. 
    ## create folder structure 
    dir <- file.path(paste0(GeoGenCorr_wd,'GIS_shpFiles/world_countries_10m'))
    if(!dir.exists(dir)){
      dir.create(dir, recursive = TRUE)
    }
    # export the file
    terra::writeVector(x = world_poly_clip, filename = path)
  }
  return(world_poly_clip)
}




#' prepWorldAdmin
#'
#' @param worldAdmin : full world admin feature
#' @param occurranceData : tabular records with lat long of know species occurrances 
#' 
#' @description
#' We buffer a convex hull of the occurrance data for the species to 50km. This is used to 
#' Select the countries that overlap this area. These countries are filterer out of the larger
#' admin file and dissolve so there are no internal political boundaries.  
#' 
#' @return A geographic limits and disvolved verison of the world admin layer that only includes elements that intersect the occurrences
#' 
prepWorldAdmin <- function(world_poly_clip, wildPoints){

  # buffer the wild points to max expected buffer distance 
  bbox <- vect(wildPoints, 
               geom = c("decimalLongitude", "decimalLatitude"),
               crs = "epsg:4326")|> # setting wgs1984
    terra::convHull()|>
    terra::buffer(width = 50000)# Unit is meter if x has a longitude/latitude CRS
  # intersect with the admin feature   
  uniqueLocs <- terra::extract(world_poly_clip, bbox) |> 
    dplyr::select(adm0_a3)|>
    dplyr::pull()
  # filter world admin to countries on overlap
  ### this is silly... but the dplyr filter is not working on the vect object... 
  admin <- world_poly_clip |>
    sf::st_as_sf()|>
    dplyr::filter(adm0_a3 %in% uniqueLocs)|>
    terra::vect()|> 
    terra::aggregate()
                      
  # return the simplified object 
  return(admin)
}
