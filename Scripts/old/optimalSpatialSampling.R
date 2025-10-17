###
# generalized method for determining the least number of locations require to get maximum coverage of a sdm 
###

## why is is the true random sample not rendering to 5% coverage? 
## QULO: missing the true random results 

# environment 
pacman::p_load(terra, dplyr,readr, tmap)
tmap_mode("view")

# function for testing erase method 
defineTopArea <-function(c1, coverageOrder, areaCoverage, originalArea,
                         randomSelection, random = FALSE, init = TRUE){
  # 022025 change 
  ## recalcualting the crop area after this fuction had removed data 
  c1$cropArea <- round(terra::expanse(c1, unit = "km"), digits = 2)
  
  # drop all features with cropArea == 0 
  c1 <- c1[c1$cropArea >0, ]
  
  
  # true random 
  
  # filter to the max area of coverage 
  ## selecting the buffer object that covers the greatest area of the distribution 
  if(random == TRUE){
    # select one of the top 10 feautres
    s1 <- c1[order(c1$cropArea, decreasing = TRUE),]
    # draw sample 
    if(init == TRUE){
      b1 <- s1[randomSelection,] 
    }else{
      if(nrow(s1) >= 10){
        b1 <- s1[sample(1:10, 1),] 
      }else{
        b1 <- s1[sample(1:nrow(s1), 1),1]
      }
    }
  }else{
    b1 <- c1[c1$cropArea == max(c1$cropArea), ]
  }

  
  
  # select the feature of interest 
  ## if two or more features share the same area of coverage, randomly select one 
  if(nrow(b1) > 1){
    # random seleciton 
    selection <- b1[sample(nrow(b1))[1],]
  }else{
    selection <- b1
  }
  
  
  
  # assign name to the feautre 
  ## this is how we track what is removed 
  coverageOrder <- c(coverageOrder, selection$id)
  # Erase the area
  ## remove all the area covered by the selected buffer from the full buffer object 
  c2 <- terra::erase(x = c1, y = selection)
  # remove all zeros agian... 
  c2 <- c2[c2$cropArea >0, ]
  # Get the new area of coverage for each buffer 
  ca <- round(terra::expanse(c2, unit = "km"), digits = 2)
  length(ca)
  c2$cropArea <- ca
  # calculate new total Area 
  newArea <- round(terra::expanse(x = terra::aggregate(c2), unit="km"), digits = 2)
  # get the percent of original area present
  ## some cases this is returning two values which is causing an error in the constrution of the dataframe object
  areaCoverage <- c(areaCoverage, (newArea/originalArea)*100)
  # returns the erase spatial object, the name of the feature used and the areacover percentage
  return(list(spatialObject = terra::makeValid(c2),
              orderList = coverageOrder,
              areaList = areaCoverage))
}


# function to render the results  -----------------------------------------
produceGeoOptimization <- function(buffDist, raster, points, exportPath,
                                   random = FALSE){
  # create path if it doesn't exist 
  if(!dir.exists(exportPath)){
    dir.create(exportPath)
  }
  
  # convert raster to vector 
  r2 <- terra::as.polygons(raster)
  # select the predicted presence area
  r2 <- r2[r2$Threshold == 1, ]
  
  # testing different buffer areas 
    bufferRadius <- buffDist  # in meters
    buffers <- buffer(points, width = bufferRadius)
    # get the original area of the buffered objects 
    ## adding a rounding step
    buffers$fullArea <- round(expanse(buffers, unit = "km"), digits = 2)
    # crop buffers to model area
    croppedBuffers <- crop(buffers, r2) 
    # new area 
    # croppedBuffers$cropArea <- round(expanse(croppedBuffers, unit = "km"), digits = 2) 
    
    # assign some variable names for the function
    c1 <- croppedBuffers
    # dissolve the vect feature to get the total area coverd by all buffers 
    originalArea <- aggregate(c1) |> terra::expanse(unit = "km") |>
      round(,digits = 2)
    # empty vectors for stroing information 
    coverageOrder <- c()
    areaCoverage <- c()
    
    # repeating this with different random seeds to 
    # export plots
    # Create a PDF file
    if(random==TRUE){
      pdf(paste0(exportPath, "/bufferEval_",buffDist,"_random.pdf")) 
    }else{
      pdf(paste0(exportPath, "/bufferEval_",buffDist,".pdf")) 
    }

    
    outputDF <- list()
    # test 1 max area draw and 5 random selection
    for(seed in 1:6){
      if(seed == 1){
        terra::plot(c1, main =paste0("seed:", seed, " removed: 0"))
        for(element in 1:nrow(c1)){
          print(element)
          if(element ==1){
            out1 <- defineTopArea(c1 = c1,
                                  coverageOrder = coverageOrder,
                                  areaCoverage = areaCoverage,
                                  originalArea = originalArea,
                                  random = FALSE,
                                  randomSelection = NA,
                                  init = FALSE)
          }else{
            out1 <- defineTopArea(c1 = out1$spatialObject, 
                                  coverageOrder = out1$orderList,
                                  areaCoverage = out1$areaList,
                                  originalArea = originalArea,
                                  random = FALSE,
                                  randomSelection = NA,
                                  init=FALSE)
          }
          # generate the plot 
          if (element %% 10 == 0) {
            title <- paste0("Seed: ",seed, " removed:", element)
            try(terra::plot(out1$spatialObject, main = title))
          }
          # this stops the workflow when the remaining area of buffered objects is less then 5% of when it started 
          ## implemented to address some un resolved errors which were breaking the workflow 
          if(min(out1$areaList) < 5){
            break  
          }
        }
        # combine the end results into a dataframe 
        print(paste0(element, " bind df"))
        df <- data.frame(
          siteID = out1$orderList,
          areaCoverage = round(out1$areaList, digits = 4)
        )
        names(df) <- c(paste0("siteID_MAX",seed), paste0("areaCoverage_MAX",seed))
        
        outputDF[[seed]] <- df
      }else{
        # random element 
        set.seed(seed)
        terra::plot(c1, main =paste0("seed:", seed, " removed: 0"))
        randomSelection <- sample(x = 1:10, 5)
        # loop thought the features -- this works.... 
        for(element in 1:nrow(c1)){
          print(element)
          if(element ==1){
            out1 <- defineTopArea(c1 = c1,
                                  coverageOrder = coverageOrder,
                                  areaCoverage = areaCoverage,
                                  originalArea = originalArea,
                                  random = random,
                                  randomSelection = randomSelection[seed],
                                  init = TRUE)
          }else{
            out1 <- defineTopArea(c1 = out1$spatialObject, 
                                  coverageOrder = out1$orderList,
                                  areaCoverage = out1$areaList,
                                  originalArea = originalArea,
                                  random = TRUE,
                                  randomSelection = randomSelection[seed],
                                  init=FALSE)
          }
          # generate the plot 
          if (element %% 10 == 0) {
            title <- paste0("Seed: ",seed, " removed:", element)
            try(terra::plot(out1$spatialObject, main = title))
          }
          # this stops the workflow when the remaining area of buffered objects is less then 5% of when it started 
          ## implemented to address some un resolved errors which were breaking the workflow 
          if(min(out1$areaList) < 5){
            break  
          }
        }
        # combine the end results into a dataframe 
        print(paste0(element, " bind df"))
        df <- data.frame(
          siteID = out1$orderList,
          areaCoverage = round(out1$areaList, digits = 4)[1:length(out1$orderList)]
        )
        names(df) <- c(paste0("siteID_Random",seed), paste0("areaCoverage_Random",seed))
        
        outputDF[[seed]] <- df
      }
    }
    # Close the PDF device
    dev.off() 
    totalRows <- max(lapply(X = outputDF, FUN = nrow) |> unlist())
    # add additional rows if needed 
    for(val in 1:length(outputDF)){
      d1 <- outputDF[[val]]
      dif <- totalRows - nrow(d1)
      print(dif)
      if(dif != 0){
        print(val)
        # Create a data frame of NAs with the correct number of rows and columns
        padding_df <- data.frame(matrix(NA, nrow = dif, ncol = ncol(d1)))
        colnames(padding_df) <- colnames(df) # Important: Set column names!
        
        # reassign values
        outputDF[[val]] <- bind_rows(d1, padding_df)
      }
    }
    ### wierd column names 
    allData <- bind_cols(outputDF) |> 
      select(where(~!is.na(.[1])))
    names(allData) <- c(
      "siteID_MAX1"              ,"areaCoverage_MAX1",       
      "siteID_Random2"           ,"areaCoverage_Radom2",     
      "siteID_Random3"           ,"areaCoverage_Radom3",     
      "siteID_Random4"           ,"areaCoverage_Radom4",     
      "siteID_Random5"           ,"areaCoverage_Radom5",     
       "siteID_Random6"          ,"areaCoverage_Radom6"
    )

      
    
    
    
    #Export 
    if(random == TRUE){
      write_csv(allData, 
                file = paste0(exportPath, "/bufferEval_",as.character(buffDist),"_random.csv"))
    }else{
      write_csv(allData, 
                file = paste0(exportPath, "/bufferEval_",as.character(buffDist),".csv"))
    }
  # }
}


# true random 
defineTopAreaTrueRandom <-function(c1, coverageOrder, areaCoverage, originalArea){
  # 022025 change 
  ## recalcualting the crop area after this fuction had removed data 
  c1$cropArea <- round(terra::expanse(c1, unit = "km"), digits = 2)
  
  # drop all features with cropArea == 0 
  c1 <- c1[c1$cropArea >0, ]
  
  
  # true random 
  b1 <- c1[sample(1:nrow(c1), 1), ]
  
  
  # select the feature of interest 
  ## if two or more features share the same area of coverage, randomly select one 
  if(nrow(b1) > 1){
    # random seleciton 
    selection <- b1[sample(nrow(b1))[1],]
  }else{
    selection <- b1
  }
  
  # assign name to the feautre 
  ## this is how we track what is removed 
  coverageOrder <- c(coverageOrder, selection$id)
  # Erase the area
  ## remove all the area covered by the selected buffer from the full buffer object 
  c2 <- terra::erase(x = c1, y = selection)
  # remove all zeros agian... 
  c2 <- c2[c2$cropArea >0, ]
  # Get the new area of coverage for each buffer 
  ca <- round(terra::expanse(c2, unit = "km"), digits = 2)
  length(ca)
  c2$cropArea <- ca
  # calculate new total Area 
  newArea <- round(terra::expanse(x = terra::aggregate(c2), unit="km"), digits = 2)
  # get the percent of original area present
  ## some cases this is returning two values which is causing an error in the constrution of the dataframe object
  areaCoverage <- c(areaCoverage, (newArea/originalArea)*100)
  # returns the erase spatial object, the name of the feature used and the areacover percentage
  return(list(spatialObject = terra::makeValid(c2),
              orderList = coverageOrder,
              areaList = areaCoverage))
}

produceGeoOptimizationTrueRandom <- function(buffDist,
                                             raster,
                                         points,
                                         exportPath
                                         ){
  # create path if it doesn't exist 
  if(!dir.exists(exportPath)){
    dir.create(exportPath)
  }
  
  # convert raster to vector 
  r2 <- terra::as.polygons(raster)
  # select the predicted presence area
  r2 <- r2[r2$Threshold == 1, ]
  
  # testing different buffer areas 
  # for(dist in c(5000, 10000, 50000,250000,500000)){
  # Buffer the points by 50km
  bufferRadius <- buffDist  # in meters
  buffers <- buffer(points, width = bufferRadius)
  # get the original area of the buffered objects 
  ## adding a rounding step
  buffers$fullArea <- round(expanse(buffers, unit = "km"), digits = 2)
  # crop buffers to model area
  croppedBuffers <- crop(buffers, r2) 
  # new area 
  # croppedBuffers$cropArea <- round(expanse(croppedBuffers, unit = "km"), digits = 2) 
  
  # assign some variable names for the function
  c1 <- croppedBuffers
  # dissolve the vect feature to get the total area coverd by all buffers 
  originalArea <- aggregate(c1) |> terra::expanse(unit = "km") |>
    round(,digits = 2)
  # empty vectors for stroing information 
  coverageOrder <- c()
  areaCoverage <- c()
  
  # repeating this with different random seeds to 
  # export plots
  # Create a PDF file
  pdf(paste0(exportPath, "/bufferEval_",buffDist,"_trueRandom.pdf")) 
  
  outputDF <- list()
  for(seed in 1:5){
    set.seed(seed)
    terra::plot(c1, main =paste0("seed:", seed, " removed: 0"))
    # loop thought the features -- this works.... 
    for(i in 1:nrow(c1)){
      print(i)
      if(i ==1){
        out1 <- defineTopAreaTrueRandom(c1 = c1,
                              coverageOrder = coverageOrder,
                              areaCoverage = areaCoverage,
                              originalArea = originalArea)
      }else{
        out1 <- defineTopArea(c1 = out1$spatialObject, 
                              coverageOrder = out1$orderList,
                              areaCoverage = out1$areaList,
                              originalArea = originalArea)
      }
      # generate the plot 
      if (i %% 10 == 0) {
        title <- paste0("Seed: single rep - removed:", i)
        try(terra::plot(out1$spatialObject, main = title))
      }
      # this stops the workflow when the remaining area of buffered objects is less then 5% of when it started 
      ## implemented to address some un resolved errors which were breaking the workflow 
      if(min(out1$areaList) < 10){
        break  
      }
    }
    
    # combine the end results into a dataframe 
    print(paste0(i, " bind df"))
    df <- data.frame(
      siteID = out1$orderList,
      areaCoverage = round(out1$areaList, digits = 4)
    )
    names(df) <- c(paste0("siteID_"), paste0("areaCoverage_"))
    
    outputDF[[seed]] <- df
  }
  # Close the PDF device
  dev.off()
  totalRows <- max(lapply(X = outputDF, FUN = nrow) |> unlist())
  # add additional rows if needed 
  for(val in 1:length(outputDF)){
    d1 <- outputDF[[val]]
    dif <- totalRows - nrow(d1)
    print(dif)
    if(dif != 0){
      print(val)
      # Create a data frame of NAs with the correct number of rows and columns
      padding_df <- data.frame(matrix(NA, nrow = dif, ncol = ncol(d1)))
      colnames(padding_df) <- colnames(df) # Important: Set column names!
      
      # reassign values
      outputDF[[val]] <- bind_rows(d1, padding_df)
    }
  }
  ### wierd column names
  allData <- bind_cols(outputDF) |> 
    select(where(~!is.na(.[1])))
  names(allData) <- c(
    "siteID_trueRandom1"              ,"areaCoverage_trueRandom1",       
    "siteID_trueRandom2"           ,"areaCoverage_trueRandom2",     
    "siteID_trueRandom3"           ,"areaCoverage_trueRandom3",     
    "siteID_trueRandom4"           ,"areaCoverage_trueRandom4",     
    "siteID_trueRandom5"           ,"areaCoverage_trueRandom5"
  )
  
  
  
  #Export 
  write_csv(allData, 
            file = paste0(exportPath, "/bufferEval_",as.character(buffDist),"_trueRandom.csv"))
  return(allData)
}



