# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% AMSONIA THARPII : GEOGRAPHIC OPTIMIZATION ANALYSIS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat('---- GEOGRAPHIC OPTIMIZATION: ALT GENETIC METRICS RUN: AMTH ----')

# This script utilizes the Amsonia tharpii dataset from Cohen et al. 2024 and calculates
# genetic metrics in "optimized" subsets of individuals. It first calculates specified genetic metrics
# of geographically-optimized sample sets. Then, it calculates the same metrics 
# in genetically-optimized sample sets (determined by the CoreHunter package). Finally, it calculates
# the same metrics in random sets of individuals. The generated metrics are saved to a disk and plotted
pacman::p_load(adegenet, terra, parallel, RColorBrewer, viridis, 
               scales, vcfR, usedist, corehunter, dplyr, ggplot2, patchwork)
# Read in relevant functions
GeoGenCorr_wd <- '/home/akoontz/Documents/GeoGenCorr/Code/'
setwd(GeoGenCorr_wd)
source('Scripts/functions_GeoGenCoverage.R')
# Specify filepath for AMTH geographic and genetic data
AMTH_filePath <- paste0(GeoGenCorr_wd, 'Datasets/AMTH/')
# Specify filepath for genetic coverages of geographically optimized datasets
# If this .Rdata file exists, then the steps below to generate this file will
# be skipped, and the existing file will be read in.
AMTH_geoOptOutFile <- paste0(AMTH_filePath, 'resamplingData/AMTH_GeoOpt_AltGenMets.Rdata')

# ---- PARALLELIZATION
# Flag for running resampling steps in parallel
parFlag <- TRUE
# If running in parallel, set up cores and export required libraries
if(parFlag==TRUE){
  # Set up relevant cores
  num_cores <- detectCores() - 4
  cl <- makeCluster(num_cores)
  # Make sure libraries are on cluster (but avoid printing output)
  invisible(clusterEvalQ(cl, library('corehunter')))
}

# %%% GENETIC METRICS OF GEOGRAPHICALLY OPTIMIZED SUBSETS %%% ----
# %%% READ IN RELEVANT DATA ----
# ---- GENETIC MATRIX
# Read in the VCF file for Amsonia tharpii using vcfR::read.vcfR
AMTH_vcf <- read.vcfR(paste0(AMTH_filePath,'Genetic/Tharpii_only.vcf'))
# Convert the vcf to a genind; the return.alleles TRUE value is suggested in the function's help file
AMTH_genind <- vcfR2genind(AMTH_vcf, return.alleles = TRUE)

# ---- GEOGRAPHIC OPTIMIZATION FILES
cat(paste0('\n','--- OPTIMIZATION: Genetic coverages of geographically optimized datasets ---','\n'))
# Read in lists of AMTH geographic core sets. These are CSVs adapted from results provided by Dan Carver,
# with a unique group of core sets for each buffer size (1km, 5km, 10km, 25km, 50km, 100km, 250km)
AMTH_1km_geoCoreSets <-
  read.csv(file='/home/akoontz/Documents/GeoGenCorr/Datasets/Amsonia_tharpii/Geographic/GeoCoreSets/AMTH_geoCoreSets_1km.csv')
AMTH_5km_geoCoreSets <-
  read.csv(file='/home/akoontz/Documents/GeoGenCorr/Datasets/Amsonia_tharpii/Geographic/GeoCoreSets/AMTH_geoCoreSets_5km.csv')
AMTH_10km_geoCoreSets <-
  read.csv(file='/home/akoontz/Documents/GeoGenCorr/Datasets/Amsonia_tharpii/Geographic/GeoCoreSets/AMTH_geoCoreSets_10km.csv')
AMTH_25km_geoCoreSets <-
  read.csv(file='/home/akoontz/Documents/GeoGenCorr/Datasets/Amsonia_tharpii/Geographic/GeoCoreSets/AMTH_geoCoreSets_25km.csv')
AMTH_50km_geoCoreSets <-
  read.csv(file='/home/akoontz/Documents/GeoGenCorr/Datasets/Amsonia_tharpii/Geographic/GeoCoreSets/AMTH_geoCoreSets_50km.csv')
AMTH_100km_geoCoreSets <-
  read.csv(file='/home/akoontz/Documents/GeoGenCorr/Datasets/Amsonia_tharpii/Geographic/GeoCoreSets/AMTH_geoCoreSets_100km.csv')
AMTH_250km_geoCoreSets <-
  read.csv(file='/home/akoontz/Documents/GeoGenCorr/Datasets/Amsonia_tharpii/Geographic/GeoCoreSets/AMTH_geoCoreSets_250km.csv')
# Make a matrix geographic core sets (strictly the sample names), and name the columns according to buffer sizes
AMTH_geoCoreMat <- as.matrix(cbind(AMTH_250km_geoCoreSets[,1], AMTH_100km_geoCoreSets[,1], AMTH_50km_geoCoreSets[,1],
                                   AMTH_25km_geoCoreSets[,1], AMTH_10km_geoCoreSets[,1], AMTH_5km_geoCoreSets[,1],
                                   AMTH_1km_geoCoreSets[,1]))
colnames(AMTH_geoCoreMat) <- c('250km_1n','100km_4n', '50km_4n','25km_4n','10km_4n', '5km_4n', '1km_4n')

# %%% CALCULATE ALLELIC COVERAGES ----
cat(paste0('\n','--- GEO-OPT: Allelic coverage ','\n'))
# Declare list of different genetic metrics (Core Hunter "objectives"), besides allelic coverage
CH_objs <- list(
  objective(type="HE"),
  objective(type="EN")
)
# Build output array to store genetic values (row=sample number; column=core set; slice=genetic metric)
AMTH_geoOptGenCovs <- array(NA, dim=c(nrow(AMTH_geoCoreMat), ncol(AMTH_geoCoreMat), length(CH_objs)+1))
dimnames(AMTH_geoOptGenCovs)[[2]] <- colnames(AMTH_geoCoreMat)
dimnames(AMTH_geoOptGenCovs)[[3]] <- c('CV',unlist(lapply(CH_objs, function(x) x$type)))
# Loop through the geographic core sets, calculating the genetic metrics for each buffer size
for(i in 1:ncol(AMTH_geoCoreMat)){
  # Specify which cores set (buffer size) to analyze
  coreSet <- AMTH_geoCoreMat[,i]
  # Use sapply to iterate the gen.getAlleleCategores function for every sample size
  # This leads to entire columns of AMTH_geoOpt_genMat to be calculated at a single time
  AMTH_geoOptGenCovs[,i,1] <- sapply(seq_along(coreSet), function(i) {
    gen.getAlleleCategories(genMat=AMTH_genind@tab, samp=coreSet[1:i])[1,3]
  })
}

# %%% CALCULATE ALTERNATE METRICS ----
cat(paste0('\n','--- GEO-OPT: Alternate genetic metrics ','\n'))
# Flags for operating with and without parallelization
if (parFlag) {
  # Export relevant objects to the cluster
  clusterExport(cl, c("AMTH_geoCoreMat","AMTH_genind","evaluateCore"))
  # Build genotype object once per worker (necessary because of Java type object on cluster)
  clusterEvalQ(cl, {AMTH_chgeno <- genotypes(AMTH_genind@tab, format="biparental")})
  # Loop through CoreHunter objectives (genetic metrics)
  for (i in seq_along(CH_objs)) {
    CH_OB <- CH_objs[[i]]
    clusterExport(cl, "CH_OB")
    # Specify parallelized lapply, to iterate through columns (buffer sizes) of geo opt matrix
    res <- parLapplyLB(cl, seq_len(ncol(AMTH_geoCoreMat)), function(j){
      coreSet <- AMTH_geoCoreMat[, j]
      # Calculate genetic metric values for all sample sizes within that core set
      sapply(seq_along(coreSet), function(x){
        evaluateCore(coreSet[1:x], AMTH_chgeno, objective = CH_OB)
      })
    })
    # Populate the relevant slice of the array
    AMTH_geoOptGenCovs[,,i+1] <- do.call(cbind, res)
  }
  # End parallelization
  stopCluster(cl)
  
} else {
  # Create CoreHunter base object
  AMTH_chgeno <- genotypes(AMTH_genind@tab, format="biparental")
  # Loop through CoreHunter objectives (genetic metrics)
  for (i in seq_along(CH_objs)) {
    CH_OB <- CH_objs[[i]]
    # Iterate through columnds (buffer sizes) of geo opt matrix
    for (j in seq_len(ncol(AMTH_geoCoreMat))){
      coreSet <- AMTH_geoCoreMat[, j]
      # Calculate genetic metric values for all sample sizes within that core set
      AMTH_geoOptGenCovs[,j,i+1] <- sapply(seq_along(coreSet), function(x){
        evaluateCore(coreSet[1:x], AMTH_chgeno, objective = CH_OB)
      })
    }
  }
}

# %%% GENETIC METRICS OF RANDOMIZED SUBSETS %%% ----
cat(paste0('\n','--- OPTIMIZATION: Genetic coverages of randomized datasets ---','\n'))
# Specify resampling replicates used for randomized genetic coverages, and parallel flag
num_reps <- 5
parFlag <- FALSE
# Conditional steps based on parallel processing
if(parFlag==TRUE){
  # Set up relevant cores
  num_cores <- detectCores() - 4
  cl <- makeCluster(num_cores)
  # Make sure libraries (adegenet, terra, etc.) are on cluster (but avoid printing output)
  invisible(clusterEvalQ(cl, library('adegenet')))
  invisible(clusterEvalQ(cl, library('terra')))
  invisible(clusterEvalQ(cl, library('parallel')))
  invisible(clusterEvalQ(cl, library('usedist')))
  invisible(clusterEvalQ(cl, library('ape')))
  invisible(clusterEvalQ(cl, library('corehunter')))
  # Export necessary functions (for calculating geographic and ecological coverage) to the cluster
  clusterExport(cl, varlist = c('createBuffers','geo.compareBuff','geo.compareBuffSDM','geo.checkSDMres',
                                'eco.intersectBuff','eco.compareBuff','gen.getAlleleCategories',
                                'gen.buildDistMat', 'gen.calcGenDistCov', 'eco.totalEcoregionCount',
                                'calculateCoverage','exSituResample.Par', 'geo.gen.Resample.Par'))
  # Export relevant objects to the cluster
  clusterExport(cl, c("AMTH_genind"))
  # Build genotype object once per worker (necessary because of Java type object on cluster)
  clusterEvalQ(cl, {genCHGeno <- genotypes(AMTH_genind@tab, format="biparental")})
  # Specify array directory, and run resampling in parallel
  arrayDir <- paste0(AMTH_filePath, 'resamplingData/AMTH_GeoOpt_newGenMets.Rdata')
  # ERRORING: genCHGeno object not being found on the cluster, despite multiple methods taken.
  AMTH_randGenCovs <- geo.gen.Resample.Par(genObj = AMTH_genind, genType='EN', genCHGeno = genCHGeno,
                                           geoFlag = FALSE, ecoFlag = FALSE,
                                           reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
  # Close cores
  stopCluster(cl)
} else {
  # Run resampling, without geographic or ecological coverages. Subset to only necessary columns
  # Allelic coverage
  AMTH_randGenCovs_CV <- geo.gen.Resample(genObj=AMTH_genind, genType='CV', geoFlag=FALSE, ecoFlag=FALSE, reps=num_reps)
  AMTH_randGenCovs_CV <- AMTH_randGenCovs_CV[,1,]
  # Expected proportion of heterozygous loci
  AMTH_randGenCovs_HE <- geo.gen.Resample(genObj=AMTH_genind, genType='HE', geoFlag=FALSE, ecoFlag=FALSE, reps=num_reps)
  AMTH_randGenCovs_HE <- AMTH_randGenCovs_HE[,1,]
  # Average entry-to-nearest-entry distance
  AMTH_randGenCovs_EN <- geo.gen.Resample(genObj=AMTH_genind, genType='EN', geoFlag=FALSE, ecoFlag=FALSE, reps=num_reps)
  AMTH_randGenCovs_EN <- AMTH_randGenCovs_EN[,1,]
}

# %%% PROCESSING GENETIC COVERAGE MATRICES %%% ----
# Combine random sampling results for different genetic metrics
AMTH_randGenCovs <- list(AMTH_randGenCovs_CV, AMTH_randGenCovs_HE, AMTH_randGenCovs_EN)
# For each set of genetic metrics, calculate the average across replicates
AMTH_meanRandGenCovs <- lapply(AMTH_randGenCovs, function(x) apply(x, MARGIN=1, mean))
# Name objects according to Geo Opt array
names(AMTH_randGenCovs) <- names(AMTH_meanRandGenCovs) <- dimnames(AMTH_geoOptGenCovs)[[3]]
# Create matrices for each genetic metric
AMTH_GenCovs_CVmat <- cbind(AMTH_geoOptGenCovs[,,1], AMTH_meanRandGenCovs[[1]])
AMTH_GenCovs_HEmat <- cbind(AMTH_geoOptGenCovs[,,2], AMTH_meanRandGenCovs[[2]])
AMTH_GenCovs_ENmat <- cbind(AMTH_geoOptGenCovs[,,3], AMTH_meanRandGenCovs[[3]])
# Create an array to store the geo-optimized and randomized genetic coverages
AMTH_GenCovs <-
  array(c(AMTH_GenCovs_CVmat, AMTH_GenCovs_HEmat, AMTH_GenCovs_ENmat),
        dim=c(nrow(AMTH_geoOptGenCovs), ncol(AMTH_geoOptGenCovs)+1, length(AMTH_meanRandGenCovs)))
dimnames(AMTH_GenCovs)[[2]] <- c(colnames(AMTH_geoOptGenCovs), 'Random')
dimnames(AMTH_GenCovs)[[3]] <- names(AMTH_randGenCovs)
# Save this array to disk
saveRDS(AMTH_GenCovs, file = paste0(AMTH_filePath, 'resamplingData/AMTH_GeoOpt_AltGenMets.Rdata'))
cat('%%% GEOGRAPHIC OPTIMIZATION: ALT GENETIC METRICS RUN COMPLETE %%%')

# %%% PLOTTING: DELTA PLOTS %%% ----
# %%% READ DATA
AMTH_GenCovs <- 
  readRDS(file = paste0(AMTH_filePath, 'resamplingData/AMTH_GeoOpt_AltGenMets.Rdata'))
stopifnot(is.array(AMTH_GenCovs), length(dim(AMTH_GenCovs)) == 3)
# Mapping: metric abbreviation → full display name (with abbreviation)
# (Add or modify entries to match the metrics in your array)
metric_display <- c(
  CV = "Allelic coverage (CV)",
  HE = "Proportion of heterozygous loci (HE)",
  EN = "Entry-to-nearest-entry distance (EN)"
)
# Specify array dimensions
n_samples    <- dim(AMTH_GenCovs)[1]
strat_names  <- dimnames(AMTH_GenCovs)[[2]]
metric_names <- dimnames(AMTH_GenCovs)[[3]]
# Identify the "Random" column
rand_col <- which(strat_names == "Random")
if (length(rand_col) == 0) {
  rand_col <- length(strat_names)
  message("No column named 'Random'; using last column as baseline: ", strat_names[rand_col])
}
geo_cols <- setdiff(seq_along(strat_names), rand_col)

# %%% CLEAN STRATEGY LABELS
# Extract just the buffer distance (e.g. "250km_5n" → "250 km")
clean_strat <- function(x) {
  km_val <- gsub("^(\\d+)km.*", "\\1", x)
  paste(km_val, "km")
}
strat_labels <- clean_strat(strat_names[geo_cols])

# %%% COMPUTE DELTAS
delta_list <- list()
for (m in seq_along(metric_names)) {
  rand_vals <- AMTH_GenCovs[, rand_col, m]
  for (g in geo_cols) {
    delta_list[[length(delta_list) + 1]] <- data.frame(
      Samples  = seq_len(n_samples),
      Delta    = AMTH_GenCovs[, g, m] - rand_vals,
      Strategy = strat_labels[which(geo_cols == g)],
      Metric   = metric_names[m],
      stringsAsFactors = FALSE
    )
  }
}
df <- bind_rows(delta_list)
df$Strategy <- factor(df$Strategy, levels = strat_labels)

# Apply full display names as factor levels for facet labels
display_levels <- ifelse(metric_names %in% names(metric_display),
                         metric_display[metric_names], metric_names)
df$Metric <- factor(df$Metric, levels = metric_names, labels = display_levels)

# %%% AUTO X-AXIS LIMIT
metric_ranges <- sapply(metric_names, function(m) diff(range(AMTH_GenCovs[,,m], na.rm = TRUE)))
df <- df %>%
  left_join(
    data.frame(Metric = factor(metric_names, levels = metric_names, labels = display_levels),
               thresh = metric_ranges * 0.01),
    by = "Metric"
  )
convergence_n <- df %>%
  group_by(Strategy, Metric) %>%
  summarise(conv_n = max(which(abs(Delta) >= thresh), na.rm = TRUE), .groups = "drop") %>%
  pull(conv_n)
x_max <- min(n_samples, max(20, ceiling(max(convergence_n, na.rm = TRUE) * 1.3)))

# %%% COLORBLIND-FRIENDLY GRADIENT
# Strategies are ordered largest buffer → smallest buffer (fewest → most samples).
# viridis runs cool (purple) → warm (yellow), so large buffers / few samples
# appear cool, and small buffers / many samples appear warm.
n_strats <- length(geo_cols)
pal <- setNames(viridis::viridis(n_strats), levels(df$Strategy))

# %%% PLOT
x_breaks <- c(1, 2, 5, 10, 20, 50, 100, 200, 500)
x_breaks <- x_breaks[x_breaks <= x_max]

AMTH_plot <- ggplot(df, aes(x = Samples, y = Delta, color = Strategy)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_line(alpha = 0.75, linewidth = 0.8) +
  scale_color_manual(values = pal) +
  scale_x_continuous(
    trans = "log10",
    limits = c(1, x_max),
    breaks = x_breaks,
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  facet_wrap(~ Metric, ncol = 1, scales = "free_y") +
  labs(
    title = expression(italic("A. tharpii") * ": Genetic metric differences between geographically maximized and randomized strategies"),
    x     = "Sample size (log scale)",
    y     = "Difference from Randomized",
    color = "Geo. Max. Buffer Size"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title       = element_text(size = 14, face = "bold"),
    strip.text       = element_text(size = 12, face = "bold"),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )
# Display plot
print(AMTH_plot)
# Specify output file
outFile <- paste0(GeoGenCorr_wd, "../Documentation/Images/2026-08-14_GeoOptComps/AMTH_optGeoRand.png")
# Specify plot height, and save
plot_h <- 2.1 * length(metric_names)
ggsave(outFile, plot = AMTH_plot, width = 10, height = plot_h, dpi = 150)

# # %%% PLOTTING: COVERAGE PLOTS %%% ----
# # Read in array
# AMTH_GenCovs <- 
#   readRDS(file = paste0(AMTH_filePath, 'resamplingData/AMTH_GeoOpt_AltGenMets.Rdata'))
# # %%% REFORMAT ARRAY
# AMTH_GenCovs_df <- expand.grid(
#   Samples = 1:dim(AMTH_GenCovs)[1],
#   Dataset = 1:dim(AMTH_GenCovs)[2],
#   Metric  = 1:dim(AMTH_GenCovs)[3]
# )
# AMTH_GenCovs_df$Value <- as.vector(AMTH_GenCovs)
# # Factor labeling
# AMTH_GenCovs_df$Dataset <- factor(
#   AMTH_GenCovs_df$Dataset,
#   levels = 1:8,
#   labels = c(
#     "Geo Opt: 1 km",
#     "Geo Opt: 5 km",
#     "Geo Opt: 10 km",
#     "Geo Opt: 25 km",
#     "Geo Opt: 50 km",
#     "Geo Opt: 100 km",
#     "Geo Opt: 250 km",
#     "Randomized"
#   )
# )
# # Specify genetic metrics
# AMTH_GenCovs_df$Metric <- factor(
#   AMTH_GenCovs_df$Metric,
#   levels = 1:3,
#   labels = c(
#     "Allelic coverage (CV)",
#     "Proportion of heterozygous loci (HE)",
#     "Entry-to-nearest-entry distance (EN)"
#   )
# )
# 
# # Specify colors; non-hexadecimal for GenOpt and Randomized datasets
# cols <- c(
#   "Geo Opt: 1 km"   = "#0072B2",
#   "Geo Opt: 5 km"   = "#009E73",
#   "Geo Opt: 10 km"  = "#56B4E9",
#   "Geo Opt: 25 km"  = "#D55E00",
#   "Geo Opt: 50 km"  = "#CC79A7",
#   "Geo Opt: 100 km" = "#E69F8F",
#   "Geo Opt: 250 km" = "#661100",  
#   "Randomized"      = "black"
# )
# # 95% horizontal line (only for first panel)
# AMTH_GenCovs_df_95 <- AMTH_GenCovs_df %>% 
#   filter(Metric == "Allelic coverage (CV)") %>% 
#   group_by(Metric) %>% 
#   summarize(y95 = 95)
# 
# 
# # Common plotting function
# make_plot <- function(metric_name){
#   
#   ggplot(
#     AMTH_GenCovs_df %>% filter(Metric == metric_name),
#     aes(x = Samples, y = Value, color = Dataset, group = Dataset)
#   ) +
#     # GeoOpt datasets
#     geom_point(
#       data = \(x) subset(x, Dataset != "Randomized"),
#       size = 1.7,
#       alpha = 0.6
#     ) +
#     # Randomized dataset
#     geom_point(
#       data = \(x) subset(x, Dataset == "Randomized"),
#       size = 1.2,
#       shape = 4,
#       stroke = 1,
#       alpha = 0.6
#     ) +
#     scale_color_manual(values = cols) +
#     labs(
#       x = "Samples",
#       y = "Genetic Metric Values",
#       color = "Dataset"
#     ) +
#     theme_bw(base_size = 14) +
#     theme(
#       axis.title = element_text(size = 16),
#       axis.text = element_text(size = 12),
#       axis.title.y = element_text(size = 11),
#       strip.text = element_text(size = 14, face = "bold"),
#       legend.text = element_text(size = 13),
#       legend.title = element_text(size = 14),
#       legend.key.width = unit(2, "lines"),
#       legend.key.height = unit(1, "lines")
#     )
# }
# 
# # Build panels
# p1 <- make_plot("Allelic coverage (CV)") +
#   labs(y = "Allelic coverage (CV)") +
#   geom_hline(yintercept = 95, linetype = "dotted", color = "red")
# p2 <- make_plot("Proportion of heterozygous loci (HE)") +
#   labs(y = "Proportion of heterozygous loci (HE)")
# p3 <- make_plot("Entry-to-nearest-entry distance (EN)") +
#   labs(y = "Entry-to-nearest-entry distance (EN)")
# 
# # Layout
# AMTH_plot <- 
#   (p1 / (p2 | p3)) +
#   plot_layout(guides = "collect") &
#   theme(
#     legend.position = "right",
#     plot.title = element_text(size = 18, face = "bold")
#   )
# AMTH_plot <- 
#   AMTH_plot +
#   plot_annotation(
#     title = "AMTH: Genetic metrics across sampling strategies"
#   )
# 
# # Save image to PNG
# imageDir <- paste0(GeoGenCorr_wd, "../Documentation/Images/2026-06-18_GeoOptComps/")
# png(file = paste0(imageDir, 'AMTH_optGeoRand.png'), width = 1200, height = 795)
# print(AMTH_plot)
# dev.off()
