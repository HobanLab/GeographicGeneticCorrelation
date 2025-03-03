# Script for plotting SMBO datasets, using ggplot
pacman::p_load(adegenet, terra, parallel, RColorBrewer, viridis, scales, vcfR, usedist, ggplot2, reshape2)

# Read in relevant functions
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/Code/'
setwd(GeoGenCorr_wd)
source('Scripts/functions_GeoGenCoverage.R')

# %%%% FUNCTIONS %%%% ----
# Function to normalize and adjust saturation. 
adjust_saturation <- function(saturation_vector) {
  scaled <- rescale(abs(saturation_vector), to = c(0.5, 0.1))  
  scaled[which.min(saturation_vector)] <- 1  # Full saturation for lowest value column
  return(scaled)
}

# Function to process data frames with separate colors and saturations
process_df <- function(df, panel_name, color_palette, saturation_scaled, genSat=0.25) {
  df_long <- melt(df, variable.name = "Series", value.name = "Value")
  df_long$Index <- rep(1:nrow(df), times = ncol(df))
  df_long$Series <- as.character(df_long$Series)
  
  gradient_colors <- colorRampPalette(color_palette)(100)
  color_mapping <- gradient_colors[as.numeric(cut(saturation_scaled, breaks = 100))]
  
  color_values <- c("red", color_mapping)
  df_long$Color <- color_values[match(df_long$Series, colnames(df))]
  
  alpha_values <- c(genSat, saturation_scaled)
  df_long$Alpha <- alpha_values[match(df_long$Series, colnames(df))]
  
  df_long$PlotOrder <- ifelse(df_long$Series == colnames(df)[1], 2, 1)
  df_long <- df_long[order(df_long$PlotOrder), ]
  
  df_long$Panel <- panel_name
  
  return(df_long)
}

# Define legend colors/labels. This is relict code, but is required for ggplot to get colors right...
legend_colors <- c("red", "darkblue", "brown", "purple", "gray")
legend_labels <- c("genetic", "Geographic (Total buffer)", "Geographic (SDM)", "Ecological", "Low match")

# %%%% QULO %%%% ----
# Read in QULO SMBO2 resampling array amd convert to data.frame
QULO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QULO/')
QULO_arrayDir <- paste0(QULO_filePath, 'resamplingData/SMBO2/QULO_SMBO2_G2E_5r_resampArr.Rdata')
QULO_array <- readRDS(QULO_arrayDir)
# Create average value matrix
QULO_avMat <- meanArrayValues(QULO_array)
# Divide into several average value matrices, depending on coverage type
QULO_avMat_GeoBuff <- QULO_avMat[,c(1,grep('Geo_Buff', colnames(QULO_avMat)))]
QULO_avMat_GeoSDM <- QULO_avMat[,c(1,grep('Geo_SDM', colnames(QULO_avMat)))]
# Remove first 7 SDM columns--outliers making graph illegible
QULO_avMat_GeoSDM <- QULO_avMat_GeoSDM[,-(2:8)] 
QULO_avMat_Eco <- QULO_avMat[,c(1,grep('Eco_Buff', colnames(QULO_avMat)))]
# Calculate NRMSE values
QULO_NRMSEs <-  buildCorrelationMat(QULO_arrayDir, corMetric = 'NRMSE', sdmFlag = TRUE)
# Create separate saturation vectors, based on each set of NRMSEs; filter based on columns
# present in average value matrix
satVect_GeoBuff <- QULO_NRMSEs[,1]
# Remove first 7 SDM columns--outliers making graph illegible
satVect_GeoSDM <- QULO_NRMSEs[,2][-(1:7)] 
satVect_Eco <- QULO_NRMSEs[,3]

# Specify saturations for each coverage type
satScaled_GeoBuff <- adjust_saturation(satVect_GeoBuff)
satScaled_GeoSDM <- adjust_saturation(satVect_GeoSDM)
satScaled_Eco <- adjust_saturation(satVect_Eco)

# Process individual dataframes to specified color saturation
QULO_avMat_GeoBuff_long <- process_df(QULO_avMat_GeoBuff, "Geo. (Total buffer)", c("gray", "darkblue"), satScaled_GeoBuff)
QULO_avMat_GeoSDM_long <- process_df(QULO_avMat_GeoSDM, "Geo. (SDM)", c("gray", "brown"), satScaled_GeoSDM)
QULO_avMat_Eco_long <- process_df(QULO_avMat_Eco, "Eco.", c("gray", "purple"), satScaled_Eco)
# Combine coverages into single, large dataframe
QULO_combined <- rbind(QULO_avMat_GeoBuff_long, QULO_avMat_GeoSDM_long, QULO_avMat_Eco_long)
# Specify order of panels
QULO_combined$Panel <- 
  factor(QULO_combined$Panel, levels = c("Geo. (Total buffer)", "Geo. (SDM)", "Eco."))

# Plot with facets
ggplot(QULO_combined, aes(x = Index, y = Value, group = Series, color = Color, alpha = Alpha)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  scale_size_identity() +  # Ensures manual size mapping
  scale_color_identity(guide = "legend", breaks = legend_colors, labels = legend_labels) +
  scale_alpha_identity(guide = "none") +
  expand_limits(y = c(0, 105)) +
  guides(color = guide_legend(title = "Legend", override.aes = list(size = 4))) +
  labs(
    title = "Q. lobata: Genetic, Geographic, and Ecological Coverages",
    subtitle = "436 individuals; 41 buffer sizes (0.5 km -- 500 km); 5 replicates",
    x = "Number of Samples",
    y = "Coverage (%)"
  ) +
  facet_grid(factor(Panel, levels = c("Geo. (Total buffer)", "Geo. (SDM)", "Eco.")) ~ ., scales = "free_y") +
  theme_minimal() +
  theme(
    legend.box.background = element_rect(color = "black", fill = "white"),
    legend.position = c(0.8, 0.1)  
  )

# %%%% PICO %%%% ----
# Read in PICO SMBO2 resampling array amd convert to data.frame
PICO_filePath <- paste0(GeoGenCorr_wd, 'Datasets/PICO/')
PICO_arrayDir <- paste0(PICO_filePath, 'resamplingData/SMBO2_G2E/PICO_SMBO2_G2E_5r_resampArr.Rdata')
PICO_array <- readRDS(PICO_arrayDir)
# Create average value matrix
PICO_avMat <- meanArrayValues(PICO_array)
# Divide into several average value matrices, depending on coverage type
PICO_avMat_GeoBuff <- PICO_avMat[,c(1,grep('Geo_Buff', colnames(PICO_avMat)))]
PICO_avMat_GeoSDM <- PICO_avMat[,c(1,grep('Geo_SDM', colnames(PICO_avMat)))]
# Remove first 10 SDM columns--outliers making graph illegible
PICO_avMat_GeoSDM <- PICO_avMat_GeoSDM[,-(2:11)] 
PICO_avMat_Eco <- PICO_avMat[,c(1,grep('Eco_Buff', colnames(PICO_avMat)))]
# Calculate NRMSE values
PICO_NRMSEs <-  buildCorrelationMat(PICO_arrayDir, corMetric = 'NRMSE', sdmFlag = TRUE)
# Create separate saturation vectors, based on each set of NRMSEs; filter based on columns
# present in average value matrix
satVect_GeoBuff <- PICO_NRMSEs[,1]
# Remove first 10 SDM columns--outliers making graph illegible
satVect_GeoSDM <- PICO_NRMSEs[,2][-(1:10)] 
satVect_Eco <- PICO_NRMSEs[,3]

# Specify saturations for each coverage type
satScaled_GeoBuff <- adjust_saturation(satVect_GeoBuff)
satScaled_GeoSDM <- adjust_saturation(satVect_GeoSDM)
satScaled_Eco <- adjust_saturation(satVect_Eco)

# Process individual dataframes to specified color saturation
PICO_avMat_GeoBuff_long <- process_df(PICO_avMat_GeoBuff, "Geo. (Total buffer)", c("gray", "darkblue"), satScaled_GeoBuff)
PICO_avMat_GeoSDM_long <- process_df(PICO_avMat_GeoSDM, "Geo. (SDM)", c("gray", "brown"), satScaled_GeoSDM)
PICO_avMat_Eco_long <- process_df(PICO_avMat_Eco, "Eco.", c("gray", "purple"), satScaled_Eco)
# Combine coverages into single, large dataframe
PICO_combined <- rbind(PICO_avMat_GeoBuff_long, PICO_avMat_GeoSDM_long, PICO_avMat_Eco_long)
# Specify order of panels
PICO_combined$Panel <- 
  factor(PICO_combined$Panel, levels = c("Geo. (Total buffer)", "Geo. (SDM)", "Eco."))

# Plot with facets
ggplot(PICO_combined, aes(x = Index, y = Value, group = Series, color = Color, alpha = Alpha)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  scale_size_identity() +  # Ensures manual size mapping
  scale_color_identity(guide = "legend", breaks = legend_colors, labels = legend_labels) +
  scale_alpha_identity(guide = "none") +
  expand_limits(y = c(0, 105)) +
  guides(color = guide_legend(title = "Legend", override.aes = list(size = 4))) +
  labs(
    title = "P. contorta: Genetic, Geographic, and Ecological Coverages",
    subtitle = "929 individuals; 41 buffer sizes (0.5 km -- 500 km); 5 replicates",
    x = "Number of Samples",
    y = "Coverage (%)"
  ) +
  facet_grid(factor(Panel, levels = c("Geo. (Total buffer)", "Geo. (SDM)", "Eco.")) ~ ., scales = "free_y") +
  theme_minimal() +
  theme(
    legend.box.background = element_rect(color = "black", fill = "white"),
    legend.position = c(0.8, 0.1)  
  )

# %%%% MIGU %%%% ----
# Read in MIGU SMBO2 resampling array amd convert to data.frame
MIGU_filePath <- paste0(GeoGenCorr_wd, 'Datasets/MIGU/')
MIGU_arrayDir <- paste0(MIGU_filePath, 'resamplingData/SMBO2_G2E/MIGU_SMBO2_G2E_5r_resampArr.Rdata')
MIGU_array <- readRDS(MIGU_arrayDir)
# Create average value matrix
MIGU_avMat <- meanArrayValues(MIGU_array)
# Divide into several average value matrices, depending on coverage type
MIGU_avMat_GeoBuff <- MIGU_avMat[,c(1,grep('Geo_Buff', colnames(MIGU_avMat)))]
MIGU_avMat_GeoSDM <- MIGU_avMat[,c(1,grep('Geo_SDM', colnames(MIGU_avMat)))]
# Remove first 13 SDM columns--outliers making graph illegible
MIGU_avMat_GeoSDM <- MIGU_avMat_GeoSDM[,-(2:14)] 
MIGU_avMat_Eco <- MIGU_avMat[,c(1,grep('Eco_Buff', colnames(MIGU_avMat)))]
# Calculate NRMSE values
MIGU_NRMSEs <-  buildCorrelationMat(MIGU_arrayDir, corMetric = 'NRMSE', sdmFlag = TRUE)
# Create separate saturation vectors, based on each set of NRMSEs; filter based on columns
# present in average value matrix
satVect_GeoBuff <- MIGU_NRMSEs[,1]
# Remove first 13 SDM columns--outliers making graph illegible
satVect_GeoSDM <- MIGU_NRMSEs[,2][-(1:13)] 
satVect_Eco <- MIGU_NRMSEs[,3]

# Specify saturations for each coverage type
satScaled_GeoBuff <- adjust_saturation(satVect_GeoBuff)
satScaled_GeoSDM <- adjust_saturation(satVect_GeoSDM)
satScaled_Eco <- adjust_saturation(satVect_Eco)

# Process individual dataframes to specified color saturation
MIGU_avMat_GeoBuff_long <- process_df(MIGU_avMat_GeoBuff, "Geo. (Total buffer)", c("gray", "darkblue"), satScaled_GeoBuff)
MIGU_avMat_GeoSDM_long <- process_df(MIGU_avMat_GeoSDM, "Geo. (SDM)", c("gray", "brown"), satScaled_GeoSDM)
MIGU_avMat_Eco_long <- process_df(MIGU_avMat_Eco, "Eco.", c("gray", "purple"), satScaled_Eco)
# Combine coverages into single, large dataframe
MIGU_combined <- rbind(MIGU_avMat_GeoBuff_long, MIGU_avMat_GeoSDM_long, MIGU_avMat_Eco_long)
# Specify order of panels
MIGU_combined$Panel <- 
  factor(MIGU_combined$Panel, levels = c("Geo. (Total buffer)", "Geo. (SDM)", "Eco."))

# Plot with facets
ggplot(MIGU_combined, aes(x = Index, y = Value, group = Series, color = Color, alpha = Alpha)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  scale_size_identity() +  # Ensures manual size mapping
  scale_color_identity(guide = "legend", breaks = legend_colors, labels = legend_labels) +
  scale_alpha_identity(guide = "none") +
  expand_limits(y = c(0, 105)) +
  guides(color = guide_legend(title = "Legend", override.aes = list(size = 4))) +
  labs(
    title = "M. guttatus: Genetic, Geographic, and Ecological Coverages",
    subtitle = "255 individuals; 41 buffer sizes (0.5 km -- 500 km); 5 replicates",
    x = "Number of Samples",
    y = "Coverage (%)"
  ) +
  facet_grid(factor(Panel, levels = c("Geo. (Total buffer)", "Geo. (SDM)", "Eco.")) ~ ., scales = "free_y") +
  theme_minimal() +
  theme(
    legend.box.background = element_rect(color = "black", fill = "white"),
    legend.position = c(0.8, 0.1)  
  )

# %%%% ARTH %%%% ----
# Read in ARTH SMBO2 resampling array amd convert to data.frame
ARTH_filePath <- paste0(GeoGenCorr_wd, 'Datasets/ARTH/')
ARTH_arrayDir <- paste0(ARTH_filePath, 'resamplingData/ARTH_SMBO2_GE_5r_resampArr.Rdata')
ARTH_array <- readRDS(ARTH_arrayDir)
# Create average value matrix
ARTH_avMat <- meanArrayValues(ARTH_array)
# Divide into several average value matrices, depending on coverage type
ARTH_avMat_GeoBuff <- ARTH_avMat[,c(1,grep('Geo_Buff', colnames(ARTH_avMat)))]
ARTH_avMat_Eco <- ARTH_avMat[,c(1,grep('Eco_Buff', colnames(ARTH_avMat)))]
# Calculate NRMSE values
ARTH_NRMSEs <-  buildCorrelationMat(ARTH_arrayDir, corMetric = 'NRMSE', sdmFlag = FALSE)
# Create separate saturation vectors, based on each set of NRMSEs; filter based on columns
# present in average value matrix
satVect_GeoBuff <- ARTH_NRMSEs[,1]
satVect_Eco <- ARTH_NRMSEs[,2]

# Specify saturations for each coverage type
satScaled_GeoBuff <- adjust_saturation(satVect_GeoBuff)
satScaled_Eco <- adjust_saturation(satVect_Eco)

# Process individual dataframes to specified color saturation
ARTH_avMat_GeoBuff_long <- process_df(ARTH_avMat_GeoBuff, "Geo. (Total buffer)", c("gray", "darkblue"), satScaled_GeoBuff)
ARTH_avMat_Eco_long <- process_df(ARTH_avMat_Eco, "Eco.", c("gray", "purple"), satScaled_Eco)
# Combine coverages into single, large dataframe
ARTH_combined <- rbind(ARTH_avMat_GeoBuff_long, ARTH_avMat_Eco_long)
# Specify order of panels
ARTH_combined$Panel <- 
  factor(ARTH_combined$Panel, levels = c("Geo. (Total buffer)", "Eco."))

# Plot with facets
ggplot(ARTH_combined, aes(x = Index, y = Value, group = Series, color = Color, alpha = Alpha)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  scale_size_identity() +  # Ensures manual size mapping
  scale_color_identity(guide = "legend", breaks = legend_colors, labels = legend_labels) +
  scale_alpha_identity(guide = "none") +
  expand_limits(y = c(0, 105)) +
  guides(color = guide_legend(title = "Legend", override.aes = list(size = 4))) +
  labs(
    title = "A. thaliana: Genetic, Geographic, and Ecological Coverages",
    subtitle = "1,010 individuals; 41 buffer sizes (0.5 km -- 500 km); 5 replicates",
    x = "Number of Samples",
    y = "Coverage (%)"
  ) +
  facet_grid(factor(Panel, levels = c("Geo. (Total buffer)", "Eco.")) ~ ., scales = "free_y") +
  theme_minimal() +
  theme(
    legend.box.background = element_rect(color = "black", fill = "white"),
    legend.position = c(0.8, 0.1)  
  )

# %%%% VILA %%%% ----
# Read in VILA SMBO2 resampling array amd convert to data.frame
VILA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/VILA/')
VILA_arrayDir <- paste0(VILA_filePath, 'resamplingData/VILA_SMBO2_5r_resampArr.Rdata')
VILA_array <- readRDS(VILA_arrayDir)
# Create average value matrix
VILA_avMat <- meanArrayValues(VILA_array)
# Divide into several average value matrices, depending on coverage type
VILA_avMat_GeoBuff <- VILA_avMat[,c(1,grep('Geo_Buff', colnames(VILA_avMat)))]
VILA_avMat_Eco <- VILA_avMat[,c(1,grep('Eco_Buff', colnames(VILA_avMat)))]
# Calculate NRMSE values
VILA_NRMSEs <-  buildCorrelationMat(VILA_arrayDir, corMetric = 'NRMSE', sdmFlag = FALSE)
# Create separate saturation vectors, based on each set of NRMSEs; filter based on columns
# present in average value matrix
satVect_GeoBuff <- VILA_NRMSEs[,1]
satVect_Eco <- VILA_NRMSEs[,2]

# Specify saturations for each coverage type
satScaled_GeoBuff <- adjust_saturation(satVect_GeoBuff)
satScaled_Eco <- adjust_saturation(satVect_Eco)

# Process individual dataframes to specified color saturation
VILA_avMat_GeoBuff_long <- process_df(VILA_avMat_GeoBuff, "Geo. (Total buffer)", c("gray", "darkblue"), satScaled_GeoBuff, genSat = 0.6)
VILA_avMat_Eco_long <- process_df(VILA_avMat_Eco, "Eco.", c("gray", "purple"), satScaled_Eco, genSat = 0.6)
# Combine coverages into single, large dataframe
VILA_combined <- rbind(VILA_avMat_GeoBuff_long, VILA_avMat_Eco_long)
# Specify order of panels
VILA_combined$Panel <- 
  factor(VILA_combined$Panel, levels = c("Geo. (Total buffer)", "Eco."))

# Plot with facets
ggplot(VILA_combined, aes(x = Index, y = Value, group = Series, color = Color, alpha = Alpha)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  scale_size_identity() +  # Ensures manual size mapping
  scale_color_identity(guide = "legend", breaks = legend_colors, labels = legend_labels) +
  scale_alpha_identity(guide = "none") +
  expand_limits(y = c(0, 105)) +
  guides(color = guide_legend(title = "Legend", override.aes = list(size = 4))) +
  labs(
    title = "V. labrusca: Genetic, Geographic, and Ecological Coverages",
    subtitle = "220 individuals; 41 buffer sizes (0.5 km -- 500 km); 5 replicates",
    x = "Number of Samples",
    y = "Coverage (%)"
  ) +
  facet_grid(factor(Panel, levels = c("Geo. (Total buffer)", "Eco.")) ~ ., scales = "free_y") +
  theme_minimal() +
  theme(
    legend.box.background = element_rect(color = "black", fill = "white"),
    legend.position = c(0.8, 0.1)  
  )

# %%%% QUAC %%%% ----
# Read in QUAC SMBO2 resampling array amd convert to data.frame
QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')
QUAC_arrayDir <- paste0(QUAC_filePath, 'resamplingData/QUAC_SMBO2_G2E_5r_resampArr.Rdata')
QUAC_array <- readRDS(QUAC_arrayDir)
# Create average value matrix
QUAC_avMat <- meanArrayValues(QUAC_array)
# Divide into several average value matrices, depending on coverage type
QUAC_avMat_GeoBuff <- QUAC_avMat[,c(1,grep('Geo_Buff', colnames(QUAC_avMat)))]
QUAC_avMat_GeoSDM <- QUAC_avMat[,c(1,grep('Geo_SDM', colnames(QUAC_avMat)))]
# Remove first 4 SDM columns--outliers making graph illegible
QUAC_avMat_GeoSDM <- QUAC_avMat_GeoSDM[,-(2:6)] 
QUAC_avMat_Eco <- QUAC_avMat[,c(1,grep('Eco_Buff', colnames(QUAC_avMat)))]
# Calculate NRMSE values
QUAC_NRMSEs <-  buildCorrelationMat(QUAC_arrayDir, corMetric = 'NRMSE', sdmFlag = TRUE)
# Create separate saturation vectors, based on each set of NRMSEs; filter based on columns
# present in average value matrix
satVect_GeoBuff <- QUAC_NRMSEs[,1]
# Remove first 4 SDM columns--outliers making graph illegible
satVect_GeoSDM <- QUAC_NRMSEs[,2][-(1:5)] 
satVect_Eco <- QUAC_NRMSEs[,3]

# Specify saturations for each coverage type
satScaled_GeoBuff <- adjust_saturation(satVect_GeoBuff)
satScaled_GeoSDM <- adjust_saturation(satVect_GeoSDM)
satScaled_Eco <- adjust_saturation(satVect_Eco)

# Process individual dataframes to specified color saturation
QUAC_avMat_GeoBuff_long <- process_df(QUAC_avMat_GeoBuff, "Geo. (Total buffer)", c("gray", "darkblue"), satScaled_GeoBuff, genSat = 0.7)
QUAC_avMat_GeoSDM_long <- process_df(QUAC_avMat_GeoSDM, "Geo. (SDM)", c("gray", "brown"), satScaled_GeoSDM, genSat = 0.7)
QUAC_avMat_Eco_long <- process_df(QUAC_avMat_Eco, "Eco.", c("gray", "purple"), satScaled_Eco, genSat = 0.5)
# Combine coverages into single, large dataframe
QUAC_combined <- rbind(QUAC_avMat_GeoBuff_long, QUAC_avMat_GeoSDM_long, QUAC_avMat_Eco_long)
# Specify order of panels
QUAC_combined$Panel <- 
  factor(QUAC_combined$Panel, levels = c("Geo. (Total buffer)", "Geo. (SDM)", "Eco."))

# Plot with facets
ggplot(QUAC_combined, aes(x = Index, y = Value, group = Series, color = Color, alpha = Alpha)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  scale_size_identity() +  # Ensures manual size mapping
  scale_color_identity(guide = "legend", breaks = legend_colors, labels = legend_labels) +
  scale_alpha_identity(guide = "none") +
  expand_limits(y = c(0, 105)) +
  guides(color = guide_legend(title = "Legend", override.aes = list(size = 4))) +
  labs(
    title = "Q. acerifolia: Genetic, Geographic, and Ecological Coverages",
    subtitle = "91 individuals; 41 buffer sizes (0.5 km -- 500 km); 5 replicates",
    x = "Number of Samples",
    y = "Coverage (%)"
  ) +
  facet_grid(factor(Panel, levels = c("Geo. (Total buffer)", "Geo. (SDM)", "Eco.")) ~ ., scales = "free_y") +
  theme_minimal() +
  theme(
    legend.box.background = element_rect(color = "black", fill = "white"),
    legend.position = c(0.8, 0.1)  
  )

# %%%% YUBR %%%% ----
# Read in YUBR SMBO2 resampling array amd convert to data.frame
YUBR_filePath <- paste0(GeoGenCorr_wd, 'Datasets/YUBR/')
YUBR_arrayDir <- paste0(YUBR_filePath, 'resamplingData/YUBR_SMBO2_G2E_resampArr.Rdata')
YUBR_array <- readRDS(YUBR_arrayDir)
# Create average value matrix
YUBR_avMat <- meanArrayValues(YUBR_array)
# Divide into several average value matrices, depending on coverage type
YUBR_avMat_GeoBuff <- YUBR_avMat[,c(1,grep('Geo_Buff', colnames(YUBR_avMat)))]
YUBR_avMat_GeoSDM <- YUBR_avMat[,c(1,grep('Geo_SDM', colnames(YUBR_avMat)))]
# Remove first 2 SDM columns--outliers making graph illegible
YUBR_avMat_GeoSDM <- YUBR_avMat_GeoSDM[,-(2:3)] 
YUBR_avMat_Eco <- YUBR_avMat[,c(1,grep('Eco_Buff', colnames(YUBR_avMat)))]
# Calculate NRMSE values
YUBR_NRMSEs <-  buildCorrelationMat(YUBR_arrayDir, corMetric = 'NRMSE', sdmFlag = TRUE)
# Create separate saturation vectors, based on each set of NRMSEs; filter based on columns
# present in average value matrix
satVect_GeoBuff <- YUBR_NRMSEs[,1]
# Remove first 2 SDM columns--outliers making graph illegible
satVect_GeoSDM <- YUBR_NRMSEs[,2][-(1:2)] 
satVect_Eco <- YUBR_NRMSEs[,3]

# Specify saturations for each coverage type
satScaled_GeoBuff <- adjust_saturation(satVect_GeoBuff)
satScaled_GeoSDM <- adjust_saturation(satVect_GeoSDM)
satScaled_Eco <- adjust_saturation(satVect_Eco)

# Process individual dataframes to specified color saturation
YUBR_avMat_GeoBuff_long <- process_df(YUBR_avMat_GeoBuff, "Geo. (Total buffer)", c("gray", "darkblue"), satScaled_GeoBuff, genSat = 0.5)
YUBR_avMat_GeoSDM_long <- process_df(YUBR_avMat_GeoSDM, "Geo. (SDM)", c("gray", "brown"), satScaled_GeoSDM, genSat = 0.5)
YUBR_avMat_Eco_long <- process_df(YUBR_avMat_Eco, "Eco.", c("gray", "purple"), satScaled_Eco, genSat = 0.5)
# Combine coverages into single, large dataframe
YUBR_combined <- rbind(YUBR_avMat_GeoBuff_long, YUBR_avMat_GeoSDM_long, YUBR_avMat_Eco_long)
# Specify order of panels
YUBR_combined$Panel <- 
  factor(YUBR_combined$Panel, levels = c("Geo. (Total buffer)", "Geo. (SDM)", "Eco."))

# Plot with facets
ggplot(YUBR_combined, aes(x = Index, y = Value, group = Series, color = Color, alpha = Alpha)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  scale_size_identity() +  # Ensures manual size mapping
  scale_color_identity(guide = "legend", breaks = legend_colors, labels = legend_labels) +
  scale_alpha_identity(guide = "none") +
  expand_limits(y = c(0, 105)) +
  guides(color = guide_legend(title = "Legend", override.aes = list(size = 4))) +
  labs(
    title = "Y. brevifolia: Genetic, Geographic, and Ecological Coverages",
    subtitle = "319 individuals; 41 buffer sizes (0.5 km -- 500 km); 5 replicates",
    x = "Number of Samples",
    y = "Coverage (%)"
  ) +
  facet_grid(factor(Panel, levels = c("Geo. (Total buffer)", "Geo. (SDM)", "Eco.")) ~ ., scales = "free_y") +
  theme_minimal() +
  theme(
    legend.box.background = element_rect(color = "black", fill = "white"),
    legend.position = c(0.8, 0.1)  
  )

# %%%% COGL %%%% ----
# Read in COGL SMBO2 resampling array amd convert to data.frame
COGL_filePath <- paste0(GeoGenCorr_wd, 'Datasets/COGL/')
COGL_arrayDir <- paste0(COGL_filePath, 'resamplingData/COGL_SMBO2_GE_5r_resampArr.Rdata')
COGL_array <- readRDS(COGL_arrayDir)
# Create average value matrix
COGL_avMat <- meanArrayValues(COGL_array)
# Divide into several average value matrices, depending on coverage type
COGL_avMat_GeoBuff <- COGL_avMat[,c(1,grep('Geo_Buff', colnames(COGL_avMat)))]
COGL_avMat_Eco <- COGL_avMat[,c(1,grep('Eco_Buff', colnames(COGL_avMat)))]
# Calculate NRMSE values
COGL_NRMSEs <-  buildCorrelationMat(COGL_arrayDir, corMetric = 'NRMSE', sdmFlag = FALSE)
# Create separate saturation vectors, based on each set of NRMSEs; filter based on columns
# present in average value matrix
satVect_GeoBuff <- COGL_NRMSEs[,1]
satVect_Eco <- COGL_NRMSEs[,2]

# Specify saturations for each coverage type
satScaled_GeoBuff <- adjust_saturation(satVect_GeoBuff)
satScaled_Eco <- adjust_saturation(satVect_Eco)

# Process individual dataframes to specified color saturation
COGL_avMat_GeoBuff_long <- process_df(COGL_avMat_GeoBuff, "Geo. (Total buffer)", c("gray", "darkblue"), satScaled_GeoBuff)
COGL_avMat_Eco_long <- process_df(COGL_avMat_Eco, "Eco.", c("gray", "purple"), satScaled_Eco)
# Combine coverages into single, large dataframe
COGL_combined <- rbind(COGL_avMat_GeoBuff_long, COGL_avMat_Eco_long)
# Specify order of panels
COGL_combined$Panel <- 
  factor(COGL_combined$Panel, levels = c("Geo. (Total buffer)", "Eco."))

# Plot with facets
ggplot(COGL_combined, aes(x = Index, y = Value, group = Series, color = Color, alpha = Alpha)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  scale_size_identity() +  # Ensures manual size mapping
  scale_color_identity(guide = "legend", breaks = legend_colors, labels = legend_labels) +
  scale_alpha_identity(guide = "none") +
  expand_limits(y = c(0, 105)) +
  guides(color = guide_legend(title = "Legend", override.aes = list(size = 4))) +
  labs(
    title = "C. glabra: Genetic, Geographic, and Ecological Coverages",
    subtitle = "562 individuals; 41 buffer sizes (0.5 km -- 500 km); 5 replicates",
    x = "Number of Samples",
    y = "Coverage (%)"
  ) +
  facet_grid(factor(Panel, levels = c("Geo. (Total buffer)", "Eco.")) ~ ., scales = "free_y") +
  theme_minimal() +
  theme(
    legend.box.background = element_rect(color = "black", fill = "white"),
    legend.position = c(0.8, 0.1)  
  )

# %%%% AMTH %%%% ----
# Read in AMTH SMBO2 resampling array amd convert to data.frame
AMTH_filePath <- paste0(GeoGenCorr_wd, 'Datasets/AMTH/')
AMTH_arrayDir <- paste0(AMTH_filePath, 'resamplingData/AMTH_SMBO2_GE_5r_resampArr.Rdata')
AMTH_array <- readRDS(AMTH_arrayDir)
# Create average value matrix
AMTH_avMat <- meanArrayValues(AMTH_array)
# Divide into several average value matrices, depending on coverage type
AMTH_avMat_GeoBuff <- AMTH_avMat[,c(1,grep('Geo_Buff', colnames(AMTH_avMat)))]
AMTH_avMat_Eco <- AMTH_avMat[,c(1,grep('Eco_Buff', colnames(AMTH_avMat)))]
# Calculate NRMSE values
AMTH_NRMSEs <-  buildCorrelationMat(AMTH_arrayDir, corMetric = 'NRMSE', sdmFlag = FALSE)
# Create separate saturation vectors, based on each set of NRMSEs; filter based on columns
# present in average value matrix
satVect_GeoBuff <- AMTH_NRMSEs[,1]
satVect_Eco <- AMTH_NRMSEs[,2]

# Specify saturations for each coverage type
satScaled_GeoBuff <- adjust_saturation(satVect_GeoBuff)
satScaled_Eco <- adjust_saturation(satVect_Eco)

# Process individual dataframes to specified color saturation
AMTH_avMat_GeoBuff_long <- process_df(AMTH_avMat_GeoBuff, "Geo. (Total buffer)", c("gray", "darkblue"), satScaled_GeoBuff, genSat = 0.5)
AMTH_avMat_Eco_long <- process_df(AMTH_avMat_Eco, "Eco.", c("gray", "purple"), satScaled_Eco, genSat = 0.5)
# Combine coverages into single, large dataframe
AMTH_combined <- rbind(AMTH_avMat_GeoBuff_long, AMTH_avMat_Eco_long)
# Specify order of panels
AMTH_combined$Panel <- 
  factor(AMTH_combined$Panel, levels = c("Geo. (Total buffer)", "Eco."))

# Plot with facets
ggplot(AMTH_combined, aes(x = Index, y = Value, group = Series, color = Color, alpha = Alpha)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  scale_size_identity() +  # Ensures manual size mapping
  scale_color_identity(guide = "legend", breaks = legend_colors, labels = legend_labels) +
  scale_alpha_identity(guide = "none") +
  expand_limits(y = c(0, 105)) +
  guides(color = guide_legend(title = "Legend", override.aes = list(size = 4))) +
  labs(
    title = "A. tharpii: Genetic, Geographic, and Ecological Coverages",
    subtitle = "140 individuals; 41 buffer sizes (0.5 km -- 500 km); 5 replicates",
    x = "Number of Samples",
    y = "Coverage (%)"
  ) +
  facet_grid(factor(Panel, levels = c("Geo. (Total buffer)", "Eco.")) ~ ., scales = "free_y") +
  theme_minimal() +
  theme(
    legend.box.background = element_rect(color = "black", fill = "white"),
    legend.position = c(0.8, 0.1)  
  )

# %%%% HIWA %%%% ----
# Read in HIWA SMBO2 resampling array amd convert to data.frame
HIWA_filePath <- paste0(GeoGenCorr_wd, 'Datasets/HIWA/')
HIWA_arrayDir <- paste0(HIWA_filePath, 'resamplingData/HIWA_SMBO2_GE_5r_resampArr.Rdata')
HIWA_array <- readRDS(HIWA_arrayDir)
# Create average value matrix
HIWA_avMat <- meanArrayValues(HIWA_array)
# Divide into several average value matrices, depending on coverage type
HIWA_avMat_GeoBuff <- HIWA_avMat[,c(1,grep('Geo_Buff', colnames(HIWA_avMat)))]
HIWA_avMat_Eco <- HIWA_avMat[,c(1,grep('Eco_Buff', colnames(HIWA_avMat)))]
# Calculate NRMSE values
HIWA_NRMSEs <-  buildCorrelationMat(HIWA_arrayDir, corMetric = 'NRMSE', sdmFlag = FALSE)
# Create separate saturation vectors, based on each set of NRMSEs; filter based on columns
# present in average value matrix
satVect_GeoBuff <- HIWA_NRMSEs[,1]
satVect_Eco <- HIWA_NRMSEs[,2]

# Specify saturations for each coverage type
satScaled_GeoBuff <- adjust_saturation(satVect_GeoBuff)
satScaled_Eco <- adjust_saturation(satVect_Eco)

# Process individual dataframes to specified color saturation
HIWA_avMat_GeoBuff_long <- process_df(HIWA_avMat_GeoBuff, "Geo. (Total buffer)", c("gray", "darkblue"), satScaled_GeoBuff, genSat = 0.5)
HIWA_avMat_Eco_long <- process_df(HIWA_avMat_Eco, "Eco.", c("gray", "purple"), satScaled_Eco, genSat = 0.5)
# Combine coverages into single, large dataframe
HIWA_combined <- rbind(HIWA_avMat_GeoBuff_long, HIWA_avMat_Eco_long)
# Specify order of panels
HIWA_combined$Panel <- 
  factor(HIWA_combined$Panel, levels = c("Geo. (Total buffer)", "Eco."))

# Plot with facets
ggplot(HIWA_combined, aes(x = Index, y = Value, group = Series, color = Color, alpha = Alpha)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  scale_size_identity() +  # Ensures manual size mapping
  scale_color_identity(guide = "legend", breaks = legend_colors, labels = legend_labels) +
  scale_alpha_identity(guide = "none") +
  expand_limits(y = c(0, 105)) +
  guides(color = guide_legend(title = "Legend", override.aes = list(size = 4))) +
  labs(
    title = "H. waimeae: Genetic, Geographic, and Ecological Coverages",
    subtitle = "197 individuals; 41 buffer sizes (0.5 km -- 500 km); 5 replicates",
    x = "Number of Samples",
    y = "Coverage (%)"
  ) +
  facet_grid(factor(Panel, levels = c("Geo. (Total buffer)", "Eco.")) ~ ., scales = "free_y") +
  theme_minimal() +
  theme(
    legend.box.background = element_rect(color = "black", fill = "white"),
    legend.position = c(0.8, 0.1)  
  )
