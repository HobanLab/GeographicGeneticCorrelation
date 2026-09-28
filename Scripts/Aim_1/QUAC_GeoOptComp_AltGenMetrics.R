# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUERCUS ACERIFOLIA : GEOGRAPHIC OPTIMIZATION ANALYSIS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat('---- GEOGRAPHIC OPTIMIZATION: ALT GENETIC METRICS RUN: QUAC ----')

# This script utilizes the Quercus acerifolia dataset from Koontz et al. 2023 and calculates
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
# Specify filepath for QUAC geographic and genetic data
QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')
# Specify filepath for genetic coverages of geographically optimized datasets
# If this .Rdata file exists, then the steps below to generate this file will
# be skipped, and the existing file will be read in.
QUAC_geoOptOutFile <- paste0(QUAC_filePath, 'resamplingData/QUAC_GeoOpt_AltGenMets.Rdata')

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

## %%% GENETIC METRICS OF GEOGRAPHICALLY OPTIMIZED SUBSETS %%% ----
# %%% READ IN RELEVANT DATA ----
# ---- GENETIC MATRIX
# Read in genind file: Reference aligment; R80, min-maf=0, first SNP/locus, only wild individuals,
# no Kessler individuals.
# Wild sample names/order must match those in the sample name column of the CSV (above)
QUAC_genind <- read.genepop(paste0(QUAC_filePath,'Genetic/QUAC_REF_Wild_R80_NOMAF_1SNP_NoK.gen'))
# Set all popnames for this genind object to 'wild'
pop(QUAC_genind) <- rep('wild', nInd(QUAC_genind))

# ---- GEOGRAPHIC OPTIMIZATION FILES
cat(paste0('\n','--- OPTIMIZATION: Genetic coverages of geographically optimized datasets ---','\n'))
# Read in lists of QUAC geographic core sets. These are CSVs adapted from results provided by Dan Carver,
# with a unique group of core sets for each buffer size (1km, 5km, 10km, 25km, 50km, 100km, 250km)
QUAC_1km_geoCoreSets <- read.csv(file=paste0(QUAC_filePath,'Geographic/GeoCoreSets/QUAC_geoCoreSets_1km.csv'))
QUAC_5km_geoCoreSets <- read.csv(file=paste0(QUAC_filePath,'Geographic/GeoCoreSets/QUAC_geoCoreSets_5km.csv'))
QUAC_10km_geoCoreSets <- read.csv(file=paste0(QUAC_filePath,'Geographic/GeoCoreSets/QUAC_geoCoreSets_10km.csv'))
QUAC_25km_geoCoreSets <- read.csv(file=paste0(QUAC_filePath,'Geographic/GeoCoreSets/QUAC_geoCoreSets_25km.csv'))
QUAC_50km_geoCoreSets <- read.csv(file=paste0(QUAC_filePath,'Geographic/GeoCoreSets/QUAC_geoCoreSets_50km.csv'))
QUAC_100km_geoCoreSets <- read.csv(file=paste0(QUAC_filePath,'Geographic/GeoCoreSets/QUAC_geoCoreSets_100km.csv'))
QUAC_250km_geoCoreSets <- read.csv(file=paste0(QUAC_filePath,'Geographic/GeoCoreSets/QUAC_geoCoreSets_250km.csv'))
# Make a matrix geographic core sets (strictly the sample names), and name the columns according to buffer sizes
QUAC_geoCoreMat <- as.matrix(cbind(QUAC_250km_geoCoreSets[,1], QUAC_100km_geoCoreSets[,1], QUAC_50km_geoCoreSets[,1],
                                   QUAC_25km_geoCoreSets[,1], QUAC_10km_geoCoreSets[,1], QUAC_5km_geoCoreSets[,1],
                                   QUAC_1km_geoCoreSets[,1]))
colnames(QUAC_geoCoreMat) <- c('250km_1n','100km_1n', '50km_4n','25km_4n','10km_10n', '5km_10n', '1km_31n')

# %%% CALCULATE ALLELIC COVERAGES ----
cat(paste0('\n','--- GEO-OPT: Allelic coverage ','\n'))
# Declare list of different genetic metrics (Core Hunter "objectives"), besides allelic coverage
CH_objs <- list(
  objective(type="HE"),
  objective(type="EN")
)
# Build output array to store genetic values (row=sample number; column=core set; slice=genetic metric)
QUAC_geoOptGenCovs <- array(NA, dim=c(nrow(QUAC_geoCoreMat), ncol(QUAC_geoCoreMat), length(CH_objs)+1))
dimnames(QUAC_geoOptGenCovs)[[2]] <- colnames(QUAC_geoCoreMat)
dimnames(QUAC_geoOptGenCovs)[[3]] <- c('CV',unlist(lapply(CH_objs, function(x) x$type)))
# Loop through the geographic core sets, calculating the genetic metrics for each buffer size
for(i in 1:ncol(QUAC_geoCoreMat)){
  # Specify which cores set (buffer size) to analyze
  coreSet <- QUAC_geoCoreMat[,i]
  # Use sapply to iterate the gen.getAlleleCategores function for every sample size
  # This leads to entire columns of QUAC_geoOpt_genMat to be calculated at a single time
  QUAC_geoOptGenCovs[,i,1] <- sapply(seq_along(coreSet), function(i) {
    gen.getAlleleCategories(genMat=QUAC_genind@tab, samp=coreSet[1:i])[1,3]
  })
}

# %%% CALCULATE ALTERNATE METRICS ----
cat(paste0('\n','--- GEO-OPT: Alternate genetic metrics ','\n'))
# Flags for operating with and without parallelization
if (parFlag) {
  # Export relevant objects to the cluster
  clusterExport(cl, c("QUAC_geoCoreMat","QUAC_genind","evaluateCore"))
  # Build genotype object once per worker (necessary because of Java type object on cluster)
  clusterEvalQ(cl, {QUAC_chgeno <- genotypes(QUAC_genind@tab, format="biparental")})
  # Loop through CoreHunter objectives (genetic metrics)
  for (i in seq_along(CH_objs)) {
    CH_OB <- CH_objs[[i]]
    clusterExport(cl, "CH_OB")
    # Specify parallelized lapply, to iterate through columns (buffer sizes) of geo opt matrix
    res <- parLapplyLB(cl, seq_len(ncol(QUAC_geoCoreMat)), function(j){
      coreSet <- QUAC_geoCoreMat[, j]
      # Calculate genetic metric values for all sample sizes within that core set
      sapply(seq_along(coreSet), function(x){
        evaluateCore(coreSet[1:x], QUAC_chgeno, objective = CH_OB)
      })
    })
    # Populate the relevant slice of the array
    QUAC_geoOptGenCovs[,,i+1] <- do.call(cbind, res)
  }
  # End parallelization
  stopCluster(cl)
  
} else {
  # Create CoreHunter base object
  QUAC_chgeno <- genotypes(QUAC_genind@tab, format="biparental")
  # Loop through CoreHunter objectives (genetic metrics)
  for (i in seq_along(CH_objs)) {
    CH_OB <- CH_objs[[i]]
    # Iterate through columnds (buffer sizes) of geo opt matrix
    for (j in seq_len(ncol(QUAC_geoCoreMat))){
      coreSet <- QUAC_geoCoreMat[, j]
      # Calculate genetic metric values for all sample sizes within that core set
      QUAC_geoOptGenCovs[,j,i+1] <- sapply(seq_along(coreSet), function(x){
        evaluateCore(coreSet[1:x], QUAC_chgeno, objective = CH_OB)
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
  clusterExport(cl, c("QUAC_genind"))
  # Build genotype object once per worker (necessary because of Java type object on cluster)
  clusterEvalQ(cl, {genCHGeno <- genotypes(QUAC_genind@tab, format="biparental")})
  # Specify array directory, and run resampling in parallel
  arrayDir <- paste0(QUAC_filePath, 'resamplingData/QUAC_GeoOpt_newGenMets.Rdata')
  # ERRORING: genCHGeno object not being found on the cluster, despite multiple methods taken.
  QUAC_randGenCovs <- geo.gen.Resample.Par(genObj = QUAC_genind, genType='EN', genCHGeno = genCHGeno,
                                           geoFlag = FALSE, ecoFlag = FALSE,
                                           reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
  # Close cores
  stopCluster(cl)
} else {
  # Run resampling, without geographic or ecological coverages. Subset to only necessary columns
  # Allelic coverage
  QUAC_randGenCovs_CV <- geo.gen.Resample(genObj=QUAC_genind, genType='CV', geoFlag=FALSE, ecoFlag=FALSE, reps=num_reps)
  QUAC_randGenCovs_CV <- QUAC_randGenCovs_CV[,1,]
  # Expected proportion of heterozygous loci
  QUAC_randGenCovs_HE <- geo.gen.Resample(genObj=QUAC_genind, genType='HE', geoFlag=FALSE, ecoFlag=FALSE, reps=num_reps)
  QUAC_randGenCovs_HE <- QUAC_randGenCovs_HE[,1,]
  # Average entry-to-nearest-entry distance
  QUAC_randGenCovs_EN <- geo.gen.Resample(genObj=QUAC_genind, genType='EN', geoFlag=FALSE, ecoFlag=FALSE, reps=num_reps)
  QUAC_randGenCovs_EN <- QUAC_randGenCovs_EN[,1,]
}

# %%% PROCESSING GENETIC COVERAGE MATRICES %%% ----
# Combine random sampling results for different genetic metrics
QUAC_randGenCovs <- list(QUAC_randGenCovs_CV, QUAC_randGenCovs_HE, QUAC_randGenCovs_EN)
# For each set of genetic metrics, calculate the average across replicates
QUAC_meanRandGenCovs <- lapply(QUAC_randGenCovs, function(x) apply(x, MARGIN=1, mean))
# Name objects according to objectives object (used for the geographically optimized sampling)
names(QUAC_randGenCovs) <- names(QUAC_meanRandGenCovs) <- dimnames(QUAC_geoOptGenCovs)[[3]]
# Create matrices for each genetic metric
QUAC_GenCovs_CVmat <- cbind(QUAC_geoOptGenCovs[,,1], QUAC_meanRandGenCovs[[1]])
QUAC_GenCovs_HEmat <- cbind(QUAC_geoOptGenCovs[,,2], QUAC_meanRandGenCovs[[2]])
QUAC_GenCovs_ENmat <- cbind(QUAC_geoOptGenCovs[,,3], QUAC_meanRandGenCovs[[3]])
# Create an array to store the geo-optimized and randomized genetic coverages
QUAC_GenCovs <-
  array(c(QUAC_GenCovs_CVmat, QUAC_GenCovs_HEmat, QUAC_GenCovs_ENmat),
        dim=c(nrow(QUAC_geoOptGenCovs), ncol(QUAC_geoOptGenCovs)+1, length(QUAC_meanRandGenCovs)))
dimnames(QUAC_GenCovs)[[2]] <- c(colnames(QUAC_geoOptGenCovs), 'Random')
dimnames(QUAC_GenCovs)[[3]] <- names(QUAC_randGenCovs)
# Save this array to disk
saveRDS(QUAC_GenCovs, file = paste0(QUAC_filePath, 'resamplingData/QUAC_GeoOpt_AltGenMets.Rdata'))
cat('%%% GEOGRAPHIC OPTIMIZATION: ALT GENETIC METRICS RUN COMPLETE %%%')

# %%% PLOTTING: DELTA PLOTS %%% ----
# %%% READ DATA
QUAC_GenCovs <- 
  readRDS(file = paste0(QUAC_filePath, 'resamplingData/QUAC_GeoOpt_AltGenMets.Rdata'))
stopifnot(is.array(QUAC_GenCovs), length(dim(QUAC_GenCovs)) == 3)
# Mapping: metric abbreviation → full display name (with abbreviation)
# (Add or modify entries to match the metrics in your array)
metric_display <- c(
  CV = "Allelic coverage (CV)",
  HE = "Proportion of heterozygous loci (HE)",
  EN = "Entry-to-nearest-entry distance (EN)"
)
# Specify array dimensions
n_samples    <- dim(QUAC_GenCovs)[1]
strat_names  <- dimnames(QUAC_GenCovs)[[2]]
metric_names <- dimnames(QUAC_GenCovs)[[3]]
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
  rand_vals <- QUAC_GenCovs[, rand_col, m]
  for (g in geo_cols) {
    delta_list[[length(delta_list) + 1]] <- data.frame(
      Samples  = seq_len(n_samples),
      Delta    = QUAC_GenCovs[, g, m] - rand_vals,
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
metric_ranges <- sapply(metric_names, function(m) diff(range(QUAC_GenCovs[,,m], na.rm = TRUE)))
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

QUAC_plot <- ggplot(df, aes(x = Samples, y = Delta, color = Strategy)) +
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
    title = expression(italic("Q. acerifolia") * ": Genetic metric differences between geographically maximized and randomized strategies"),
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
print(QUAC_plot)
# Specify output file
outFile <- paste0(GeoGenCorr_wd, "../Documentation/Images/2026-08-14_GeoOptComps/QUAC_optGeoRand.png")
# Specify plot height, and save
plot_h <- 2.1 * length(metric_names)
ggsave(outFile, plot = QUAC_plot, width = 10, height = plot_h, dpi = 150)

# # %%% PLOTTING: COVERAGE PLOTS %%% ----
# # Read in array
# QUAC_GenCovs <- 
#   readRDS(file = paste0(QUAC_filePath, 'resamplingData/QUAC_GeoOpt_AltGenMets.Rdata'))
# # %%% REFORMAT ARRAY
# QUAC_GenCovs_df <- expand.grid(
#   Samples = 1:dim(QUAC_GenCovs)[1],
#   Dataset = 1:dim(QUAC_GenCovs)[2],
#   Metric  = 1:dim(QUAC_GenCovs)[3]
# )
# QUAC_GenCovs_df$Value <- as.vector(QUAC_GenCovs)
# # Factor labeling
# QUAC_GenCovs_df$Dataset <- factor(
#   QUAC_GenCovs_df$Dataset,
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
# QUAC_GenCovs_df$Metric <- factor(
#   QUAC_GenCovs_df$Metric,
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
# QUAC_GenCovs_df_95 <- QUAC_GenCovs_df %>% 
#   filter(Metric == "Allelic coverage (CV)") %>% 
#   group_by(Metric) %>% 
#   summarize(y95 = 95)
# 
# 
# # Common plotting function
# make_plot <- function(metric_name){
#   
#   ggplot(
#     QUAC_GenCovs_df %>% filter(Metric == metric_name),
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
# QUAC_plot <- 
#   (p1 / (p2 | p3)) +
#   plot_layout(guides = "collect") &
#   theme(
#     legend.position = "right",
#     plot.title = element_text(size = 18, face = "bold")
#   )
# QUAC_plot <- 
#   QUAC_plot +
#   plot_annotation(
#     title = "QUAC: Genetic metrics across sampling strategies"
#   )
# 
# # Save image to PNG
# imageDir <- paste0(GeoGenCorr_wd, "../Documentation/Images/2026-06-18_GeoOptComps/")
# png(file = paste0(imageDir, 'QUAC_optGeoRand.png'), width = 1200, height = 795)
# print(QUAC_plot)
# dev.off()
