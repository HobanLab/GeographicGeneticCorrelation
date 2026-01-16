# standardize all CSV and rast objects  -----------------------------------
standardizeDataSets <- function(
  species,
  rastName,
  pointName,
  idCol,
  latLonCol
) {
  path <- paste0("Datasets/", species, "/Geographic")
  files <- list.files(path = path, full.names = TRUE)
  # MIGU export paths
  export1 <- paste0(path, "/sdmGO.gpkg")
  export2 <- paste0(path, "/pointGO.gpkg")
  # raster
  if (!file.exists(export1)) {
    print("exporting sdm")
    # read in
    r1 <- terra::rast(paste0(path, "/", rastName))
    # convert raster to vector
    r2 <- terra::as.polygons(r1)
    # select the predicted presence area
    r2 <- r2[r2$Threshold == 1, ]
    # export
    terra::writeVector(x = r2, filename = export1)
  }
  # points
  if (!file.exists(export2)) {
    print("exporting points")
    # read in rast for CRS
    r1 <- terra::rast(paste0(path, "/", rastName))
    # condition for the id column of pico with specific formatting
    if (species == "PICO") {
      d1 <- read_csv(
        paste0(path, "/", pointName),
        col_types = "cdd"
      ) |>
        dplyr::select(
          id = idCol,
          all_of(latLonCol)
        ) |>
        terra::vect(geom = latLonCol, crs = crs(r1))
      # export
      terra::writeVector(x = d1, filename = export2)
    } else {
      # read in csv and convert to vect
      d1 <- read_csv(
        paste0(path, "/", pointName)
      ) |>
        dplyr::select(
          id = idCol,
          all_of(latLonCol)
        ) |>
        terra::vect(geom = latLonCol, crs = crs(r1))
      # export
      terra::writeVector(x = d1, filename = export2)
    }
  }
  print(paste0(species, " data prepped"))
}

# call the standardization function  --------------------------------------
# Updated prepData function, for most recent 5 datasets
prepData <- function() {
  ## read in raster and point data and format
  ### assigning unique values as a vect to pass to the processing function
  species <- c("AMTH", "ARTH", "COGL", "HIWA", "VILA")
  rastName <- c(
    "AMTH_thresh.tif",
    "ARTH_thresh.tif",
    "COGL_thresh.tif",
    "HIWA_thresh.tif",
    "VILA_thresh.tif"
  )
  pointName <- c(
    "AMTH_coordinates.csv",
    "ARTH_coordinates.csv",
    "COGL_coordinates.csv",
    "HIWA_coordinates.csv",
    "VILA_coordinates.csv"
  )
  idCol <- c("sampleNames", "Acc_ID", "Sample Name", "SampleName", "sampleID")
  latLonCol <- list(
    c("decimalLongitude", "decimalLatitude"),
    c("decimalLongitude", "decimalLatitude"),
    c("Longitude", "Latitude"),
    c("decimalLongitude", "decimalLatitude"),
    c("decimalLongitude", "decimalLatitude")
  )
  for (i in seq_along(species)) {
    standardizeDataSets(
      species = species[i],
      rastName = rastName[i],
      pointName = pointName[i],
      idCol = idCol[i],
      latLonCol = latLonCol[[i]]
    )
  }
}

# Old prepData function, based on original 5 datasets
# prepData <- function() {
#   ## read in raster and point data and format
#   ### assigning unique values as a vect to pass to the processing function
#   species <- c("MIGU", "PICO", "QUAC", "QULO", "YUBR")
#   rastName <- c(
#     "MIGU_255inds_rast_Carver.tif",
#     "PICO_929inds_rast_Carver.tif",
#     "QUAC_91inds_rast.tif",
#     "QULO_436inds_rast_Carver.tif",
#     "YUBR_319inds_rast_Carver.tif"
#   )
#   pointName <- c(
#     "MIGU_coordinates.csv",
#     "PICO_coordinates.csv",
#     "QUAC_coordinates.csv",
#     "QULO_coordinates.csv",
#     "YUBR_coordinates.csv"
#   )
#   idCol <- c("Sample Name", "Internal_ID", "sampleID", "sampleID", "Sample")
#   latLonCol <- list(
#     c("Longitude", "Latitude"),
#     c("decimalLongitude", "decimalLatitude"),
#     c("decimalLongitude", "decimalLatitude"),
#     c("decimalLongitude", "decimalLatitude"),
#     c("decimalLongitude", "decimalLatitude")
#   )
#   for (i in seq_along(species)) {
#     standardizeDataSets(
#       species = species[i],
#       rastName = rastName[i],
#       pointName = pointName[i],
#       idCol = idCol[i],
#       latLonCol = latLonCol[[i]]
#     )
#   }
# }

# construct path for species ----------------------------------------------
constructPaths <- function(name) {
  core <- paste0("Datasets/", name)
  geo <- paste0(core, "/Geographic")
  opt <- paste0(core, "/geographicSampingOptimization")
  if (!dir.exists(opt)) {
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


# inital buffer  ----------------------------------------------------------
initBuffer <- function(buffDist, points, sdm) {
  # buffer features
  buffers <- terra::buffer(points, width = buffDist)
  # get the initial area of each buffer element
  buffers$fullArea <- round(expanse(buffers, unit = "m"), digits = 6)
  # crop buffers to model area
  croppedBuffers <- crop(buffers, sdm)
  return(croppedBuffers)
}


# Random Sample  ----------------------------------------------------------

randomSelect <- function(buffDist, buffer, paths, area) {
  print("stating random area selection ")
  # get original area to test percent change
  originalArea <- aggregate(buffer) |>
    terra::expanse(unit = "m") |>
    round(, digits = 6)
  # empty vectors for stroing information
  coverageOrder <- c()
  areaCoverage <- c()
  # export path
  exportPathPartial <- paste0(paths$opt, "/randomSelection_", buffDist)
  # set specific exports
  pdfPath <- paste0(exportPathPartial, ".pdf")
  csvPath <- paste0(exportPathPartial, ".csv")

  if (FALSE %in% file.exists(c(pdfPath, csvPath))) {
    # pdf of results
    pdf(paste0(exportPathPartial, ".pdf"))

    buffer <- buffer
    for (i in 1:nrow(buffer)) {
      # randomly select a buffer element
      draw <- sample(x = buffer, size = 1)
      # assign selection
      coverageOrder <- c(coverageOrder, draw$id)
      # erase the area from buffer objects
      buffer <- terra::erase(x = buffer, y = draw)
      # condition for capturing all areas
      if (nrow(buffer) == 0) {
        areaCoverage <- c(areaCoverage, 0)
        # generate the plot
        print(paste0("Sampled ", i, " with ", 0, " area left"))
        break
      } else {
        # Get the new area of coverage for each buffer
        buffer$area <- round(terra::expanse(buffer, unit = "km"), digits = 6)
        # drop all features with 0 area
        buffer <- buffer[buffer$area > 0, ]
        # calculate the total area change
        newArea <- aggregate(buffer) |>
          terra::expanse(unit = "km") |>
          round(, digits = 2)
        # calc the percent area left
        areaLeft <- (newArea / originalArea) * 100
        # assign area coverage
        areaCoverage <- c(areaCoverage, areaLeft)
        # generate the plot
        if (i %% 10 == 0) {
          print(paste0("Sampled ", i, " with ", areaLeft, " area left"))
          title <- paste0("Removed:", i)
          try(terra::plot(buffer, main = title))
        }
        if (areaLeft <= area) {
          break
          print(i)
        }
      }
    }
    # construct a dataframe of the results
    df <- data.frame(
      id = coverageOrder,
      area = areaCoverage
    )
    # export results
    write_csv(x = df, file = paste0(exportPathPartial, ".csv"))
    # Close the PDF device
    dev.off()
    # print
    print(paste0(
      "Random selection at ",
      buffDist,
      " required ",
      nrow(df),
      " selections"
    ))
  } else {
    print(paste0("Random selection at ", buffDist, " previously generated"))
  }
}


# max area sample ---------------------------------------------------------
maxSelect <- function(buffDist, buffer, points, paths, area) {
  print("stating Max area selection ")
  # get original area to test percent change
  originalArea <- aggregate(buffer) |>
    terra::expanse(unit = "km") |>
    round(, digits = 2)
  # empty vectors for stroing information
  coverageOrder <- c()
  areaCoverage <- c()
  # export path
  exportPathPartial <- paste0(paths$opt, "/maxSelection_", buffDist)
  # set specific exports
  pdfPath <- paste0(exportPathPartial, ".pdf")
  csvPath <- paste0(exportPathPartial, ".csv")
  # pull all potential ids
  ids <- buffer$id

  if (FALSE %in% file.exists(c(pdfPath, csvPath))) {
    # pdf of results
    pdf(pdfPath)
    # loop over rows of
    buffer2 <- buffer
    for (i in 1:nrow(buffer2)) {
      #nrow(buffer2)
      print(nrow(buffer2))
      # fix geomentry
      buffer2 <- terra::makeValid(buffer2)
      # calculate the current area of each buffer
      val1 <- round(terra::expanse(buffer2, unit = "km"), digits = 6)
      print(length(val1))
      buffer2$area <- val1
      # select the max area
      max <- buffer2[buffer2$area == max(buffer2$area), ] |>
        terra::makeValid()
      if (nrow(max) > 1) {
        max <- sample(x = max, 1)
      }
      # assign selection
      coverageOrder <- c(coverageOrder, max$id)
      # erase the area from buffer objects
      erasedBuff <- terra::erase(x = buffer2, y = max)

      # condition for capturing all areas
      if (nrow(erasedBuff) == 0) {
        areaCoverage <- c(areaCoverage, 0)
        # generate the plot
        print(paste0("Sampled ", i, " with ", 0, " area left"))
        break
      } else {
        # Get the new area of coverage for each buffer
        val2 <- round(terra::expanse(erasedBuff, unit = "km"), digits = 3)
        # if(length(buffer2) == length(val2)){
        erasedBuff$area <- val2
        # }
        # drop all features with 0 area
        buffer2 <- erasedBuff[erasedBuff$area > 0, ]

        # calculate the total area change
        newArea <- aggregate(buffer2) |>
          terra::expanse(unit = "km") |>
          round(, digits = 6)
        # calc the percent area left
        areaLeft <- (newArea / originalArea) * 100
        # assign area coverage
        areaCoverage <- c(areaCoverage, areaLeft)
        # generate the plot
        if (i %% 10 == 0) {
          print(paste0("Sampled ", i, " with ", areaLeft, " captured"))
          title <- paste0("Removed:", i)
          try(terra::plot(buffer2, main = title))
        }
        if (areaLeft <= area) {
          break
          print(i)
        }
      }
    }
    # construct a dataframe of the results
    df <- data.frame(
      id = coverageOrder,
      area = areaCoverage
    )
    # Add all the unsampled ids to the datafrom
    ## grabbing these ids from the original point object to include buffered features
    ## that were outside of the SDM boundaries
    pIds <- points$id
    id2 <- pIds[!pIds %in% df$id]
    # randomly distribute the data
    randomIDs <- sample(id2)
    # construct data frame
    df2 <- data.frame(
      id = randomIDs,
      area = NA
    )
    # combine data
    df <- dplyr::bind_rows(df, df2)

    # export results
    write_csv(x = df, file = paste0(exportPathPartial, ".csv"))
    dev.off()
    # print
    print(paste0(
      "Max selection at ",
      buffDist,
      " required ",
      nrow(df),
      " selections"
    ))
  } else {
    print(paste0("Max selection at ", buffDist, " previously generated"))
  }
}


# adding a while loop version  --------------------------------------------
maxSelectWhile <- function(buffDist, buffer, points, paths, area) {
  print("stating Max area selection ")
  # get original area to test percent change
  originalArea <- aggregate(buffer) |>
    terra::expanse(unit = "m") |>
    round(, digits = 6)
  # empty vectors for storing information
  coverageOrder <- c()
  areaCoverage <- c()
  areaRemoved <- c()
  # export path
  exportPathPartial <- paste0(paths$opt, "/maxSelection_", buffDist)
  # set specific exports
  csvPath <- paste0(exportPathPartial, ".csv")

  # skip is data exists
  if (!file.exists(csvPath)) {
    buffer2 <- buffer
    iteration <- 1 # counter for while loop
    #standard loop
    while (nrow(buffer2) > 0) {
      #testing loop
      # while (iteration <= 170) {
      # Try to repair any broken geometries
      buffer2 <- terra::makeValid(buffer2)

      # Still having issues with invalid geoms so adding this measure. Could effect the outcome
      # but it at least gives us an option
      valid_geoms <- terra:::is.valid(buffer2)
      if (any(!valid_geoms)) {
        print(paste("removing", sum(!valid_geoms), "that are not valid"))
        buffer2 <- buffer2[valid_geoms, ]
      }

      # exclude no area buffers
      buffer2 <- buffer2[!terra::is.empty(buffer2), ]

      # End loop if nothing left
      if (nrow(buffer2) == 0) {
        print("No valid geometries left to process.")
        break
      }
      print(paste("Iteration", iteration))

      # calculate the current area of each buffer
      val1 <- round(terra::expanse(buffer2, unit = "m"), digits = 6)

      # Safety check
      if (length(val1) != nrow(buffer2)) {
        stop(paste(
          "Mismatch found: buffer has",
          nrow(buffer2),
          "rows, but expanse returned",
          length(val1),
          "values."
        ))
      }
      buffer2$area <- val1

      # select the max area
      max_features <- buffer2[buffer2$area == max(buffer2$area), ] |>
        terra::makeValid() # Good to validate the selection

      if (nrow(max_features) > 1) {
        max_feature <- sample(x = max_features, 1)
      } else {
        max_feature <- max_features
      }
      # add area removed to the vector
      areaRemoved <- c(areaRemoved, max_feature$area)

      # assign selection
      coverageOrder <- c(coverageOrder, max_feature$id)

      # erase the area from buffer objects
      erasedBuff <- terra::erase(x = buffer2, y = max_feature)

      # aggressive cleaning steps to attempt to get rid on non valid geoms
      erasedBuff <- terra::makeValid(erasedBuff)
      valid_geoms_erase <- terra::is.valid(erasedBuff)
      if (any(!valid_geoms_erase)) {
        print(paste("Removing", sum(!valid_geoms_erase), "invalid geoms."))
        erasedBuff <- erasedBuff[valid_geoms_erase, ]
      }
      erasedBuff <- erasedBuff[!terra::is.empty(erasedBuff), ]

      # condition for capturing all areas
      if (nrow(erasedBuff) == 0) {
        areaCoverage <- c(areaCoverage, 0)
        print(paste0("Sampled ", iteration, " with ", 0, " area left"))
        break
      } else {
        # Filter out any nonpolygon features
        ## weird nastyness to work this one out... basically the crop was generating points and line
        ## features that were resulting in some errorsd
        erasedBuff <- erasedBuff[
          terra::geomtype(erasedBuff) %in% c("polygons", "multipolygon"),
        ]

        # If the filter removed everything, end while loop
        if (nrow(erasedBuff) == 0) {
          areaCoverage <- c(areaCoverage, 0)
          print(paste0(
            "Sampled ",
            iteration,
            " with ",
            0,
            " area left (post-geomtype filter)."
          ))
          break
        }

        # Get the new area of coverage for each buffer
        val2 <- round(terra::expanse(erasedBuff, unit = "m"), digits = 3)

        if (length(val2) != nrow(erasedBuff)) {
          stop(paste(
            "Mismatch post-erase: buffer has",
            nrow(erasedBuff),
            "rows, but expanse returned",
            length(val2),
            "values."
          ))
        }

        erasedBuff$area <- val2

        # drop all features with 0 area
        buffer2 <- erasedBuff[erasedBuff$area > 0, ]

        # Check if filtering removed everything
        if (nrow(buffer2) == 0) {
          areaCoverage <- c(areaCoverage, 0)
          print(paste0(
            "Sampled ",
            iteration,
            " with ",
            0,
            " area left (post-area filter)."
          ))
          break
        }

        # calculate the total area change
        newArea <- aggregate(buffer2) |>
          terra::expanse(unit = "m") |>
          round(, digits = 6)

        # calc the percent area left
        areaLeft <- (newArea / originalArea) * 100

        # assign area coverage
        areaCoverage <- c(areaCoverage, areaLeft)

        # Print progress
        if (iteration %% 10 == 0) {
          print(paste0("Sampled ", iteration, " with ", areaLeft, " captured"))
        }

        if (areaLeft <= area) {
          print(paste("Target area reached at iteration", iteration))
          break
        }
      }
      # add to counter
      iteration <- iteration + 1
    } # End of while loop

    # results df
    df <- data.frame(
      id = coverageOrder,
      area = areaCoverage,
      area_removed = areaRemoved
    )

    # Add all the unsampled ids to the dataframe using the point object
    pIds <- points$id
    id2 <- pIds[!pIds %in% df$id]

    # Check if there are any unsampled IDs before sampling
    if (length(id2) > 0) {
      # randomly distribute the data
      randomIDs <- sample(id2)
      # construct data frame
      df2 <- data.frame(
        id = randomIDs,
        area = NA,
        area_removed = NA
      )
      # combine data
      df <- dplyr::bind_rows(df, df2)
    } else {
      print("All point IDs were sampled.")
    }

    # export results
    write_csv(x = df, file = csvPath)

    # print
    print(paste0(
      "Max selection at ",
      buffDist,
      " required ",
      length(coverageOrder),
      " selections"
    ))
  } else {
    print(paste0("Max selection at ", buffDist, " previously generated"))
  }
}


# run geo optimized sampling method  --------------------------------------
runGeoSelection <- function(buffDist, species, area) {
  paths <- constructPaths(species)
  points <- terra::vect(paste0(paths$geo, "/pointGO.gpkg"))
  sdm <- terra::vect(paste0(paths$geo, "/sdmGO.gpkg"))
  buffer <- initBuffer(points = points, sdm = sdm, buffDist = buffDist)
  # ran1 <- randomSelect(buffer = buffer, buffDist = buffDist, paths = paths,area = area)
  max1 <- maxSelectWhile(
    buffer = buffer,
    buffDist = buffDist,
    points = points,
    paths = paths,
    area = area
  )
}


# remove files for clean reruns
removeFiles <- function(doit = FALSE, species = NA) {
  # max files
  max <- list.files(
    recursive = TRUE,
    pattern = "maxSelection",
    full.names = TRUE
  )
  # random files
  random <- list.files(
    recursive = TRUE,
    pattern = "randomSelection",
    full.names = TRUE
  )

  files <- c(max, random)
  if (!is.na(species)) {
    files <- files[grepl(pattern = species, x = files)]
  }
  # test condition
  if (doit == TRUE) {
    print("Deleting all results")
    result <- file.remove(files)
  } else {
    print("set doit to true to remove files")
  }
}
