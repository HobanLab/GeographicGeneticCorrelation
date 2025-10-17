# environment 
pacman::p_load(terra, dplyr,readr, tmap, furrr, purrr, tictoc)
tmap_mode("view")

source("Scripts/optimalSpatialSampling.R")



# functions  --------------------------------------------------------------
runSampling <- function(rastPath, pointPath, latLonCol, exportPath,idColumn,trueRandom = FALSE, furrr = TRUE){
  # read in datasets for MIGU 
  r1 <- terra::rast(rastPath)
  # # point data 
  p1 <- read_csv(pointPath) |>
    dplyr::select(all_of(c(idColumn, latLonCol)))
  s1 <- terra::vect(p1, geom = latLonCol, crs = r1)
  # standardize names of the point obect
  ## should probably do this inside the function
  names(s1) <- 'id'
  # # call process
  buffDist <- c(5000,10000,25000,50000,100000,250000,500000)
  if(trueRandom == TRUE){
    furrr::future_map(.x = buffDist,.f = produceGeoOptimizationTrueRandom,
                      raster = r1,
                      points = s1,
                      exportPath = exportPath)
  }else{
    if(furrr == TRUE){
      furrr::future_map(.x = buffDist,.f = produceGeoOptimization,
                        raster = r1,
                        points = s1,
                        exportPath = exportPath,
                        random = TRUE)
    }else{
      purrr::map(.x = buffDist,.f = produceGeoOptimization,
                 raster = r1,
                 points = s1,
                 exportPath = exportPath,
                 random = TRUE)
    }
  }
  

  
}

constructPaths <- function(name){
  core <- paste0("Datasets/", name)
  geo <- paste0(core, "/Geographic")
  opt <- paste0(core, "/geographicSampingOptimization")
  if(!dir.exists(opt)){
    dir.create(opt)
  }
  return(
    paths = list(
     core = core,
     geo = geo,
     opt = opt
    )
  )
}

# run results for the species of interest  --------------------------------
# set the furrr environmen
future::plan("multicore", workers = 8)
fullList <- c("MIGU", "PICO", "QUAC", "QULO","YUBR")
runList <- fullList[5:5]


## MIGU --------------------------------------------------------------------
paths <- constructPaths("MIGU")
core <- paths$core
rastPath <- paste0(paths$geo, "/MIGU_255inds_rast_Carver.tif")
pointPath <- paste0(paths$geo, "/MIGU_coordinates.csv")
# look at names to assign lat lon column 
names(read_csv(pointPath))
latLonCol <- c("Longitude", "Latitude")
idColumn <- "Sample Name"
exportPath <- paths$opt
# run method 
# run time : 299.081 sec elapsed
# tic()
# runSampling(rastPath, pointPath, latLonCol, exportPath, furrr = FALSE)
# toc()
# run time : 67.444 sec elapsed
if("MIGU" %in% runList){
  # optimized 
  # runSampling(rastPath, pointPath, latLonCol, exportPath, idColumn, furrr = TRUE)
  # true random
  runSampling(rastPath, pointPath, latLonCol, exportPath, idColumn, trueRandom = TRUE, furrr = TRUE)
}


## QUAC --------------------------------------------------------------------
paths <- constructPaths("QUAC")
core <- paths$core
rastPath <- paste0(paths$geo, "/QUAC_91inds_rast.tif")
pointPath <- paste0(paths$geo, "/QUAC_coordinates.csv")
# look at names to assign lat lon column 
names(read_csv(pointPath))
latLonCol <- c("decimalLongitude", "decimalLatitude")
idColumn <- "sampleID"
exportPath <- paths$opt
if("QUAC" %in% runList){
  # runSampling(rastPath, pointPath, latLonCol, exportPath,idColumn, furrr = TRUE)
  runSampling(rastPath, pointPath, latLonCol, exportPath, idColumn, trueRandom = TRUE, furrr = TRUE)
}



runSampling
## QULO --------------------------------------------------------------------
paths <- constructPaths("QULO")
core <- paths$core
rastPath <- paste0(paths$geo, "/QULO_436inds_rast_Carver.tif")
pointPath <- paste0(paths$geo, "/QULO_coordinates.csv")
# look at names to assign lat lon column 
names(read_csv(pointPath))
latLonCol <- c("decimalLongitude", "decimalLatitude")
idColumn <- "sampleID"
exportPath <- paths$opt
if("QULO" %in% runList){
  runSampling(rastPath, pointPath, latLonCol, exportPath,idColumn, furrr = TRUE)
}
## YUBR  --------------------------------------------------------------------
paths <- constructPaths("YUBR")
core <- paths$core
rastPath <- paste0(paths$geo, "/YUBR_319inds_rast_Carver.tif")
pointPath <- paste0(paths$geo, "/YUBR_coordinates.csv")
# look at names to assign lat lon column 
names(read_csv(pointPath))
latLonCol <- c("decimalLongitude", "decimalLatitude")
idColumn <- "Sample"
exportPath <- paths$opt
if("YUBR" %in% runList){
  # runSampling(rastPath, pointPath, latLonCol, exportPath,idColumn, furrr = TRUE)
  runSampling(rastPath, pointPath, latLonCol, exportPath, idColumn, trueRandom = TRUE, furrr = TRUE)
  
}
## PICO --------------------------------------------------------------------
paths <- constructPaths("PICO")
core <- paths$core
rastPath <- paste0(paths$geo, "/PICO_929inds_rast_Carver.tif")
pointPath <- paste0(paths$geo, "/PICO_coordinates.csv")
# look at names to assign lat lon column 
names(read_csv(pointPath))
latLonCol <- c("Longitude", "Latitude")
idColumn <- "Internal_ID"
exportPath <- paths$opt
if("PICO" %in% runList){
  # runSampling(rastPath, pointPath, latLonCol, idColumn, exportPath, furrr = TRUE)
  runSampling(rastPath, pointPath, latLonCol, exportPath, idColumn, trueRandom = TRUE, furrr = TRUE)
  
}


