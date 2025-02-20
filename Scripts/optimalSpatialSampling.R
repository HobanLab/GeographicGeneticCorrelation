###
# generalized method for determining the least number of locations require to get maximum coverage of a sdm 
###
# environment 
pacman::p_load(terra, dplyr,readr)


# Second attempt  ---------------------------------------------------------
## here I just want to brute for the process 
## start with the highest coverage points 
## if multiple objects are persent at the same rank, select randomly 
## erase that area from the full distribution 
## recalculate the area and rank then select the next highest coverage features 



# function for testing erase method 
defineTopArea <-function(c1, coverageOrder, areaCoverage, originalArea, random = FALSE){
  # 022025 change 
  ## recalcualting the crop area after this fuction had removed data 
  c1$cropArea <- round(terra::expanse(c1, unit = "km"), digits = 2)
  
  # drop all features with cropArea == 0 
  c1 <- c1[c1$cropArea >0, ]
  
  
  # filter to the max area of coverage 
  ## selecting the buffer object that covers the greatest area of the distribution 
  if(random == TRUE){
    # select one of the top 10 feautres
    s1 <- c1[order(c1$cropArea, decreasing = TRUE),]
    b1 <- s1[sample(1:10, 1), ]
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





produceGeoOptimization <- function(raster, points, buffDist, exportPath, random = FALSE){
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
    buffers <- buffer(s1, width = bufferRadius)
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
    for(seed in 1:10){
      set.seed(seed)
      terra::plot(c1, main =paste0("seed:", seed, " removed: 0"))
      # loop thought the features -- this works.... 
      for(i in 1:nrow(c1)){
        print(i)
        if(i ==1){
          out1 <- defineTopArea(c1 = c1,
                                coverageOrder = coverageOrder,
                                areaCoverage = areaCoverage,
                                originalArea = originalArea,
                                random = random)
        }else{
          out1 <- defineTopArea(c1 = out1$spatialObject, 
                                coverageOrder = out1$orderList,
                                areaCoverage = out1$areaList,
                                originalArea = originalArea,
                                random = random)
          print(paste0(i, "_out1"))
        }
        # generate the plot 
        if (i %% 10 == 0) {
          title <- paste0("Seed: ",seed, " removed:", i)
          try(terra::plot(out1$spatialObject, main = title))
        }
        # this stops the workflow when the remaining area of buffered objects is less then 5% of when it started 
        ## implemented to address some un resolved errors which were breaking the workflow 
        if(min(out1$areaList) < 5){
          break  
        }
      }
      
      # combine the end results into a dataframe 
      print(paste0(i, " bind df"))
      df <- data.frame(
        siteID = out1$orderList,
        areaCoverage = round(out1$areaList, digits = 4)
      )
      names(df) <- c(paste0("siteID_",seed), paste0("areaCoverage_",seed))
      
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
    
    allData <- bind_cols(outputDF)
    
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




# read in datasets for MIGU 
r1 <- terra::rast("Datasets/MIGU/Geographic/MIGU_255inds_rast_Carver.tif")
# # point data 
p1 <- read_csv("Datasets/MIGU/Geographic/MIGU_coordinates.csv")
s1 <- terra::vect(p1, geom = c("Longitude", "Latitude"), crs = r1)
# standardize names of the point obect
## should probably do this inside the function
names(s1) <- 'id'
# # call process
for(i in c(5000,10000,25000,50000,100000,250000,500000)){
  produceGeoOptimization(raster = r1,
                         points = s1,
                         buffDist = i,
                         exportPath = "Datasets/MIGU/geographicSampingOptimization",
                         random = TRUE)
}

# 

# run records for QULO ----------------------------------------------------
r1 <- terra::rast("Datasets/QULO/Geographic/QULO_436inds_rast_Carver.tif")
p1 <- read_csv("Datasets/QULO/Geographic/QULO_coordinates.csv")
s1 <- terra::vect(p1, geom = c("decimalLongitude", "decimalLatitude"), crs = r1)
# standardize names of the point obect 
names(s1) <- 'id'
# call process
produceGeoOptimization(raster = r1,
                       points = s1, 
                       buffDist = 50000, 
                       exportPath = "Datasets/QULO/geographicSampingOptimization",
                       random = TRUE)
 


# 
# 
# # View(p1)
# ## transform the datasets 
# # convert raster to vector 
# r2 <- terra::as.polygons(r1)
# # select the predicted presence area
# r2 <- r2[r2$Threshold == 1, ]
# 
# # testing different buffer areas 
# for(dist in c(5000, 10000, 50000,250000,500000)){
#   # Buffer the points by 50km
#   bufferRadius <- dist  # in meters
#   buffers <- buffer(s1, width = bufferRadius)
#   # get the original area of the buffered objects 
#   ## adding a rounding step
#   buffers$fullArea <- round(expanse(buffers, unit = "km"), digits = 2)
#   # crop buffers to model area
#   croppedBuffers <- crop(buffers, r2) 
#   # new area 
#   # croppedBuffers$cropArea <- round(expanse(croppedBuffers, unit = "km"), digits = 2) 
#   
#   # assign some variable names for the function
#   c1 <- croppedBuffers
#   # dissolve the vect feature to get the total area coverd by all buffers 
#   originalArea <- aggregate(c1) |> terra::expanse(unit = "km") |>
#     round(,digits = 2)
#   # empty vectors for stroing information 
#   coverageOrder <- c()
#   areaCoverage <- c()
#   
#   # repeating this with different random seeds to 
#   # export plots
#   # Create a PDF file
#   pdf(paste0("Datasets/MIGU/geographicSampingOptimization/bufferEval_",dist,".pdf")) 
# 
#   outputDF <- list()
#   for(seed in 1:10){
#     set.seed(seed)
#     terra::plot(c1, main =paste0("seed:", seed, " removed: 0"))
#     # loop thought the features -- this works.... 
#     for(i in 1:nrow(c1)){
#       print(i)
#       if(i ==1){
#         out1 <- defineTopArea(c1 = c1,
#                               coverageOrder = coverageOrder,
#                               areaCoverage = areaCoverage,
#                               originalArea = originalArea,
#                               )
#       }else{
#         out1 <- defineTopArea(c1 = out1$spatialObject, 
#                               coverageOrder = out1$orderList,
#                               areaCoverage = out1$areaList,
#                               originalArea = originalArea)
#       }
#       # generate the plot 
#       if (i %% 10 == 0) {
#         title <- paste0("Seed: ",seed, " removed:", i)
#         terra::plot(out1$spatialObject, main = title)
#       }
#       # this stops the workflow when the remaining area of buffered objects is less then 5% of when it started 
#       ## implemented to address some un resolved errors which were breaking the workflow 
#       if(min(out1$areaList) < 5){
#         break  
#       }
#     }
#     
#     # combine the end results into a dataframe 
#     df <- data.frame(
#       siteID = out1$orderList,
#       areaCoverage = round(out1$areaList, digits = 4)
#     )
#     names(df) <- c(paste0("siteID_",seed), paste0("areaCoverage_",seed))
#     
#     outputDF[[seed]] <- df
#   }
#   # Close the PDF device
#   dev.off() 
#   totalRows <- max(lapply(X = outputDF, FUN = nrow) |> unlist())
#   # add additional rows if needed 
#   for(val in 1:length(outputDF)){
#     d1 <- outputDF[[val]]
#     dif <- totalRows - nrow(d1)
#     print(dif)
#     if(dif != 0){
#       print(val)
#       # Create a data frame of NAs with the correct number of rows and columns
#       padding_df <- data.frame(matrix(NA, nrow = dif, ncol = ncol(d1)))
#       colnames(padding_df) <- colnames(df) # Important: Set column names!
# 
#       # reassign values
#       outputDF[[val]] <- bind_rows(d1, padding_df)
#     }
#   }
#   ### wierd column names 
#   
#   allData <- bind_cols(outputDF)
#   
#   #Export 
#   write_csv(allData, 
#             file = paste0("Datasets/MIGU/geographicSampingOptimization/bufferEval_",as.character(dist),".csv"))
#   
# }
# 
# 
# r1 <- read_csv( paste0("Datasets/MIGU/geographicSampingOptimization/bufferEval_",as.character(5000),".csv"))
# 
# 
# 
# 
# 
# 
# # first attempt -----------------------------------------------------------
# 
# 
# # # read in datasets 
# # r1 <- terra::rast("Datasets/MIGU/Geographic/MIGU_255inds_rast_Carver.tif")
# # # point data 
# # p1 <- read_csv("Datasets/MIGU/Geographic/MIGU_coordinates.csv")
# # 
# # ## transform the datasets 
# # # convert raster to vector 
# # r2 <- terra::as.polygons(r1)
# # # select the predicted presence area
# # r2 <- r2[r2$Threshold == 1, ]
# # 
# # # produce a point object
# # s1 <- terra::vect(p1, geom = c("Longitude", "Latitude"), crs = r1)
# # terra::writeVector(x = s1, filename = "Datasets/MIGU/Geographic/MIGU.gpkg")
# # # buffer all locations 
# # # Buffer the points by 50km
# # bufferRadius <- 50000  # in meters
# # buffers <- buffer(s1, width = bufferRadius)
# # terra::writeVector(x = buffers, filename = "Datasets/MIGU/Geographic/MIGU_50kmBuffer.gpkg")
# # 
# # # get the original area 
# # buffers$fullArea <-expanse(buffers, unit = "km")
# # 
# # # crop buffers to model area
# # croppedBuffers <- crop(buffers, r2) 
# # # new area 
# # croppedBuffers$cropArea <- expanse(croppedBuffers, unit = "km") 
# # 
# # # generate a dataframe for holding change in area results 
# # croppedBuffers$equalArea <- croppedBuffers$fullArea - croppedBuffers$cropArea
# # 
# # # Assign a rank order 
# # croppedBuffers$rankAreaCovered <- rank(croppedBuffers$equalArea)
# # # Sort the data frame by rank
# # croppedBuffers <- croppedBuffers[order(croppedBuffers$rankAreaCovered), ]
# # 
# # 
# # # construct a dataframe to store data
# # df <- croppedBuffers |>
# #   as.data.frame() |>
# #   dplyr::mutate(erasedAreaChange = NA)
# # 
# # for(i in seq_along(croppedBuffers$`Sample Name`)){
# #   # select feature of interest 
# #   feat <- croppedBuffers[i,]
# #   name <- feat$`Sample Name`
# #   # mask the croppedBuffer features based on the 
# #   cb1 <- terra::erase(croppedBuffers, feat)
# #   cb1$area2 <- expanse(cb1, unit = "km") 
# #   # convert to df then compare areas 
# #   cb2 <- cb1 |> 
# #     as.data.frame()|>
# #     dplyr::mutate(
# #       areaDiff2 = cropArea - area2
# #     )
# #   # total area changed 
# #   totalArea <- sum(cb2$areaDiff2)
# #   #assign the values 
# #   df[i, "erasedAreaChange"] <- totalArea
# # }
# # 
# # # assign rank based on total area changes 
# # df2 <- df |> 
# #   dplyr::mutate(erasedAreaRank = rank(erasedAreaChange))
# # df_sorted <- df2[order(df2$rankAreaCovered, df2$erasedAreaRank), ]  
# # 
# 
# 
# 
