# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% CONRADINA GLABRA : GEOGRAPHIC OPTIMIZATION ANALYSIS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat('---- GEOGRAPHIC OPTIMIZATION: ALT GENETIC METRICS RUN: COGL ----')

# This script utilizes the Conradina glabra dataset from Eserman et al. 2024 and calculates
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
# Specify filepath for COGL geographic and genetic data
COGL_filePath <- paste0(GeoGenCorr_wd, 'Datasets/COGL/')
# Specify filepath for genetic coverages of geographically optimized datasets
# If this .Rdata file exists, then the steps below to generate this file will
# be skipped, and the existing file will be read in.
COGL_geoOptOutFile <- paste0(COGL_filePath, 'resamplingData/COGL_GeoOpt_AltGenMets.Rdata')

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
cat(paste0('\n','--- OPTIMIZATION: Genetic coverages of geographically optimized datasets ---','\n'))
# %%% READ IN RELEVANT DATA ----
# ---- GENETIC MATRIX
# Read in the VCF of filtered SNPs provided by Lauren Eserman. These reads were trimmed, 
# and then a "de novo reference" was created by collecting 10 million reads in total
# evenly from all samples. Reads with counts fewer than 6 and more than 700 were excluded.
# Alleles with a population frequency of less than 3% were removed.
COGL_vcf <- 
  read.vcfR(file=paste0(COGL_filePath, 'Genetic/conradina.vcf'))
# Convert the vcf to a genind; the return.alleles TRUE value is suggested in the function's help file
# This genind file is made up of 564 individuals and 16,291 loci
COGL_genind <- vcfR2genind(COGL_vcf, sep = "/", return.alleles = TRUE)
# Rename individuals to match the sample names in geographic optimization results (below)
indNames(COGL_genind) <- gsub("_sorted", "", indNames(COGL_genind))
indNames(COGL_genind) <- paste0("Congla", indNames(COGL_genind))
# Remove two duplicate samples (with "_rep" in sample name), leaving 562 samples
COGL_genind <- COGL_genind[-grep('_rep', indNames(COGL_genind)), drop=TRUE]

# ---- GEOGRAPHIC OPTIMIZATION FILES
cat(paste0('\n','--- OPTIMIZATION: Genetic coverages of geographically optimized datasets ---','\n'))
# Read in lists of COGL geographic core sets. These are CSVs adapted from results provided by Dan Carver,
# with a unique group of core sets for each buffer size (1km, 5km, 10km, 25km, 50km, 100km, 250km)
COGL_1km_geoCoreSets <- read.csv(file=paste0(COGL_filePath,'Geographic/GeoCoreSets/COGL_geoCoreSets_1km.csv'))
COGL_5km_geoCoreSets <- read.csv(file=paste0(COGL_filePath,'Geographic/GeoCoreSets/COGL_geoCoreSets_5km.csv'))
COGL_10km_geoCoreSets <- read.csv(file=paste0(COGL_filePath,'Geographic/GeoCoreSets/COGL_geoCoreSets_10km.csv'))
COGL_25km_geoCoreSets <- read.csv(file=paste0(COGL_filePath,'Geographic/GeoCoreSets/COGL_geoCoreSets_25km.csv'))
COGL_50km_geoCoreSets <- read.csv(file=paste0(COGL_filePath,'Geographic/GeoCoreSets/COGL_geoCoreSets_50km.csv'))
COGL_100km_geoCoreSets <- read.csv(file=paste0(COGL_filePath,'Geographic/GeoCoreSets/COGL_geoCoreSets_100km.csv'))
COGL_250km_geoCoreSets <- read.csv(file=paste0(COGL_filePath,'Geographic/GeoCoreSets/COGL_geoCoreSets_250km.csv'))
# Make a matrix geographic core sets (strictly the sample names), and name the columns according to buffer sizes
COGL_geoCoreMat <- as.matrix(cbind(COGL_250km_geoCoreSets[,1], COGL_100km_geoCoreSets[,1], COGL_50km_geoCoreSets[,1],
                                   COGL_25km_geoCoreSets[,1], COGL_10km_geoCoreSets[,1], COGL_5km_geoCoreSets[,1],
                                   COGL_1km_geoCoreSets[,1]))
colnames(COGL_geoCoreMat) <- c('250km_1n','100km_1n', '50km_1n','25km_1n','10km_1n', '5km_3n', '1km_26n')

# %%% CALCULATE ALLELIC COVERAGES ----
cat(paste0('\n','--- GEO-OPT: Allelic coverage ','\n'))
# Declare list of different genetic metrics (Core Hunter "objectives"), besides allelic coverage
CH_objs <- list(
  objective(type="HE"),
  objective(type="EN")
)
# Build output array to store genetic values (row=sample number; column=core set; slice=genetic metric)
COGL_geoOptGenCovs <- array(NA, dim=c(nrow(COGL_geoCoreMat), ncol(COGL_geoCoreMat), length(CH_objs)+1))
dimnames(COGL_geoOptGenCovs)[[2]] <- colnames(COGL_geoCoreMat)
dimnames(COGL_geoOptGenCovs)[[3]] <- c('CV',unlist(lapply(CH_objs, function(x) x$type)))
# Loop through the geographic core sets, calculating the genetic metrics for each buffer size
for(i in 1:ncol(COGL_geoCoreMat)){
  # Specify which cores set (buffer size) to analyze
  coreSet <- COGL_geoCoreMat[,i]
  # Use sapply to iterate the gen.getAlleleCategores function for every sample size
  # This leads to entire columns of COGL_geoOpt_genMat to be calculated at a single time
  COGL_geoOptGenCovs[,i,1] <- sapply(seq_along(coreSet), function(i) {
    gen.getAlleleCategories(genMat=COGL_genind@tab, samp=coreSet[1:i])[1,3]
  })
}

# %%% CALCULATE ALTERNATE METRICS ----
cat(paste0('\n','--- GEO-OPT: Alternate genetic metrics ','\n'))
# Flags for operating with and without parallelization
if (parFlag) {
  # Export relevant objects to the cluster
  clusterExport(cl, c("COGL_geoCoreMat","COGL_genind","evaluateCore"))
  # Build genotype object once per worker (necessary because of Java type object on cluster)
  clusterEvalQ(cl, {COGL_chgeno <- genotypes(COGL_genind@tab, format="biparental")})
  # Loop through CoreHunter objectives (genetic metrics)
  for (i in seq_along(CH_objs)) {
    CH_OB <- CH_objs[[i]]
    clusterExport(cl, "CH_OB")
    # Specify parallelized lapply, to iterate through columns (buffer sizes) of geo opt matrix
    res <- parLapplyLB(cl, seq_len(ncol(COGL_geoCoreMat)), function(j){
      coreSet <- COGL_geoCoreMat[, j]
      # Calculate genetic metric values for all sample sizes within that core set
      sapply(seq_along(coreSet), function(x){
        evaluateCore(coreSet[1:x], COGL_chgeno, objective = CH_OB)
      })
    })
    # Populate the relevant slice of the array
    COGL_geoOptGenCovs[,,i+1] <- do.call(cbind, res)
  }
  # End parallelization
  stopCluster(cl)
  
} else {
  # Create CoreHunter base object
  COGL_chgeno <- genotypes(COGL_genind@tab, format="biparental")
  # Loop through CoreHunter objectives (genetic metrics)
  for (i in seq_along(CH_objs)) {
    CH_OB <- CH_objs[[i]]
    # Iterate through columnds (buffer sizes) of geo opt matrix
    for (j in seq_len(ncol(COGL_geoCoreMat))){
      coreSet <- COGL_geoCoreMat[, j]
      # Calculate genetic metric values for all sample sizes within that core set
      COGL_geoOptGenCovs[,j,i+1] <- sapply(seq_along(coreSet), function(x){
        evaluateCore(coreSet[1:x], COGL_chgeno, objective = CH_OB)
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
  clusterExport(cl, c("COGL_genind"))
  # Build genotype object once per worker (necessary because of Java type object on cluster)
  clusterEvalQ(cl, {genCHGeno <- genotypes(COGL_genind@tab, format="biparental")})
  # Specify array directory, and run resampling in parallel
  arrayDir <- paste0(COGL_filePath, 'resamplingData/COGL_GeoOpt_newGenMets.Rdata')
  # ERRORING: genCHGeno object not being found on the cluster, despite multiple methods taken.
  COGL_randGenCovs <- geo.gen.Resample.Par(genObj = COGL_genind, genType='EN', genCHGeno = genCHGeno,
                                           geoFlag = FALSE, ecoFlag = FALSE,
                                           reps = num_reps, arrayFilepath = arrayDir, cluster = cl)
  # Close cores
  stopCluster(cl)
} else {
  # Run resampling, without geographic or ecological coverages. Subset to only necessary columns
  # Allelic coverage
  COGL_randGenCovs_CV <- geo.gen.Resample(genObj=COGL_genind, genType='CV', geoFlag=FALSE, ecoFlag=FALSE, reps=num_reps)
  COGL_randGenCovs_CV <- COGL_randGenCovs_CV[,1,]
  # Expected proportion of heterozygous loci
  COGL_randGenCovs_HE <- geo.gen.Resample(genObj=COGL_genind, genType='HE', geoFlag=FALSE, ecoFlag=FALSE, reps=num_reps)
  COGL_randGenCovs_HE <- COGL_randGenCovs_HE[,1,]
  # Average entry-to-nearest-entry distance
  COGL_randGenCovs_EN <- geo.gen.Resample(genObj=COGL_genind, genType='EN', geoFlag=FALSE, ecoFlag=FALSE, reps=num_reps)
  COGL_randGenCovs_EN <- COGL_randGenCovs_EN[,1,]
}

# %%% PROCESSING GENETIC COVERAGE MATRICES %%% ----
# Combine random sampling results for different genetic metrics
COGL_randGenCovs <- list(COGL_randGenCovs_CV, COGL_randGenCovs_HE, COGL_randGenCovs_EN)
# For each set of genetic metrics, calculate the average across replicates
COGL_meanRandGenCovs <- lapply(COGL_randGenCovs, function(x) apply(x, MARGIN=1, mean))
# Name objects according to Geo Opt array
names(COGL_randGenCovs) <- names(COGL_meanRandGenCovs) <- dimnames(COGL_geoOptGenCovs)[[3]]
# Create matrices for each genetic metric
COGL_GenCovs_CVmat <- cbind(COGL_geoOptGenCovs[,,1], COGL_meanRandGenCovs[[1]])
COGL_GenCovs_HEmat <- cbind(COGL_geoOptGenCovs[,,2], COGL_meanRandGenCovs[[2]])
COGL_GenCovs_ENmat <- cbind(COGL_geoOptGenCovs[,,3], COGL_meanRandGenCovs[[3]])
# Create an array to store the geo-optimized and randomized genetic coverages
COGL_GenCovs <-
  array(c(COGL_GenCovs_CVmat, COGL_GenCovs_HEmat, COGL_GenCovs_ENmat),
        dim=c(nrow(COGL_geoOptGenCovs), ncol(COGL_geoOptGenCovs)+1, length(COGL_meanRandGenCovs)))
dimnames(COGL_GenCovs)[[2]] <- c(colnames(COGL_geoOptGenCovs), 'Random')
dimnames(COGL_GenCovs)[[3]] <- names(COGL_randGenCovs)
# Save this array to disk
saveRDS(COGL_GenCovs, file = paste0(COGL_filePath, 'resamplingData/COGL_GeoOpt_AltGenMets.Rdata'))
cat('%%% GEOGRAPHIC OPTIMIZATION: ALT GENETIC METRICS RUN COMPLETE %%%')

# %%% PLOTTING: DELTA PLOTS %%% ----
# %%% READ DATA
COGL_GenCovs <- 
  readRDS(file = paste0(COGL_filePath, 'resamplingData/COGL_GeoOpt_AltGenMets.Rdata'))
stopifnot(is.array(COGL_GenCovs), length(dim(COGL_GenCovs)) == 3)
# Mapping: metric abbreviation → full display name (with abbreviation)
# (Add or modify entries to match the metrics in your array)
metric_display <- c(
  CV = "Allelic coverage (CV)",
  HE = "Proportion of heterozygous loci (HE)",
  EN = "Entry-to-nearest-entry distance (EN)"
)
# Specify array dimensions
n_samples    <- dim(COGL_GenCovs)[1]
strat_names  <- dimnames(COGL_GenCovs)[[2]]
metric_names <- dimnames(COGL_GenCovs)[[3]]
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
  rand_vals <- COGL_GenCovs[, rand_col, m]
  for (g in geo_cols) {
    delta_list[[length(delta_list) + 1]] <- data.frame(
      Samples  = seq_len(n_samples),
      Delta    = COGL_GenCovs[, g, m] - rand_vals,
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
metric_ranges <- sapply(metric_names, function(m) diff(range(COGL_GenCovs[,,m], na.rm = TRUE)))
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

COGL_plot <- ggplot(df, aes(x = Samples, y = Delta, color = Strategy)) +
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
    title = expression(italic("C. glabra") * ": Genetic metric differences between geographically maximized and randomized strategies"),
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
print(COGL_plot)
# Specify output file
outFile <- paste0(GeoGenCorr_wd, "../Documentation/Images/2026-08-14_GeoOptComps/COGL_optGeoRand.png")
# Specify plot height, and save
plot_h <- 2.1 * length(metric_names)
ggsave(outFile, plot = COGL_plot, width = 10, height = plot_h, dpi = 150)

# # %%% PLOTTING: COVERAGE PLOTS %%% ----
# # Read in array
# COGL_GenCovs <- 
#   readRDS(file = paste0(COGL_filePath, 'resamplingData/COGL_GeoOpt_AltGenMets.Rdata'))
# 
# # %%% REFORMAT ARRAY
# COGL_GenCovs_df <- expand.grid(
#   Samples = 1:dim(COGL_GenCovs)[1],
#   Dataset = 1:dim(COGL_GenCovs)[2],
#   Metric  = 1:dim(COGL_GenCovs)[3]
# )
# COGL_GenCovs_df$Value <- as.vector(COGL_GenCovs)
# # Factor labeling
# COGL_GenCovs_df$Dataset <- factor(
#   COGL_GenCovs_df$Dataset,
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
# COGL_GenCovs_df$Metric <- factor(
#   COGL_GenCovs_df$Metric,
#   levels = 1:3,
#   labels = c(
#     "Allelic coverage (CV)",
#     "Proportion of heterozygous loci (HE)",
#     "Entry-to-nearest-entry distance (EN)"
#   )
# )
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
# COGL_GenCovs_df_95 <- COGL_GenCovs_df %>% 
#   filter(Metric == "Allelic coverage (CV)") %>% 
#   group_by(Metric) %>% 
#   summarize(y95 = 95)
# 
# # Common plotting function
# make_plot <- function(metric_name){
#   
#   ggplot(
#     COGL_GenCovs_df %>% filter(Metric == metric_name),
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
# COGL_plot <- 
#   (p1 / (p2 | p3)) +
#   plot_layout(guides = "collect") &
#   theme(
#     legend.position = "right",
#     plot.title = element_text(size = 18, face = "bold")
#   )
# COGL_plot <- 
#   COGL_plot +
#   plot_annotation(
#     title = "COGL: Genetic metrics across sampling strategies"
#   )
# 
# # Save image to PNG
# imageDir <- paste0(GeoGenCorr_wd, "../Documentation/Images/2026-06-18_GeoOptComps/")
# png(file = paste0(imageDir, 'COGL_optGeoRand.png'), width = 1200, height = 795)
# print(COGL_plot)
# dev.off()
