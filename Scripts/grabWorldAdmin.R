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
    world_poly_clip <- download |>
      terra::vect() |>
      #this select is option but I don't think we need much of any attribute data on this file. 
      dplyr::select("admin", "adm0_a3", "name","name_long" )
    
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
