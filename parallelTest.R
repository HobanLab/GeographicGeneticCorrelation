# # Load necessary packages
# pacman::p_load(
#   adegenet,
#   terra,
#   parallel,
#   RColorBrewer,
#   viridis,
#   scales,
#   vcfR,
#   usedist,
#   corehunter
# )
# # Specif file path to input files
# GeoGenCorr_wd <- '/home/dune/trueNAS/work/GeographicGeneticCorrelation/'
# setwd(GeoGenCorr_wd)
# QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')

# # ---- PARALLELIZATION
# # Set up relevant cores
# num_cores <- detectCores() - 20
# cl <- makeCluster(num_cores)
# # Make sure libraries are on cluster (but avoid printing output)
# invisible(clusterEvalQ(cl, library('corehunter')))

# # ---- GENETIC MATRIX
# # Read in genind
# QUAC_genind <- read.genepop(paste0(
#   QUAC_filePath,
#   'Genetic/QUAC_REF_Wild_R80_NOMAF_1SNP_NoK.gen'
# ))
# # Subset QUAC genind object (to accelerate testing)
# QUAC_genind <- QUAC_genind[1:5, , drop = TRUE]
# # Build the CoreHunter input file, specifying SNPs (format=biparental)
# QUAC_chgeno <- genotypes(QUAC_genind@tab, format = "biparental")
# # Specify allelic coverage CoreHunter objective, and export
# obj <- objective(type = "CV")
# ı
# clusterExport(cl = cl, varlist = c('obj', 'QUAC_chgeno'))
# # ---- CALL COMMAND
# # Using parSapplyLB (Generates error:
# # "java.lang.IllegalArgumentException: At least one type of data (genotypic, phenotypic, distances) should be defined." --- hitting this same error
# QUAC_genOptGenCovs_CV <- parSapplyLB(
#   cl = cl,
#   # Iterate over sample sizes from 2 to one less than the total number of samples
#   1:(nInd(QUAC_genind) - 2),
#   function(i) sampleCore(QUAC_chgeno, obj = obj, size = i + 1, steps = 10)$CV
# )
# # Close cluster
# stopCluster(cl)
# # Using standard sapply (works)
# # QUAC_genOptGenCovs_CV <- sapply(
# #   # Iterate over sample sizes from 2 to one less than the total number of samples
# #   1:(nInd(QUAC_genind) - 2),
# #   function(i) sampleCore(QUAC_chgeno, obj = obj, size = i + 1, steps = 10)$CV
# # )

# attempt 2
# Load necessary packages
pacman::p_load(
  adegenet,
  terra,
  parallel,
  RColorBrewer,
  viridis,
  scales,
  vcfR,
  usedist,
  corehunter
)
# Specif file path to input files
GeoGenCorr_wd <- '/home/dune/trueNAS/work/GeographicGeneticCorrelation/'
setwd(GeoGenCorr_wd)
QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')

# ---- PARALLELIZATION
# Set up relevant cores
num_cores <- detectCores() - 24
cl <- makeCluster(num_cores)
# Make sure libraries are on cluster (but avoid printing output)
invisible(clusterEvalQ(cl, library('corehunter')))

# ---- GENETIC MATRIX
# Read in genind
QUAC_genind <- read.genepop(paste0(
  QUAC_filePath,
  'Genetic/QUAC_REF_Wild_R80_NOMAF_1SNP_NoK.gen'
))
# Subset QUAC genind object (to accelerate testing)
QUAC_genind <- QUAC_genind[1:5, , drop = TRUE]
# Extract the base-R matrix (This serializes perfectly across clusters)
QUAC_raw_matrix <- QUAC_genind@tab
# Specify allelic coverage CoreHunter objective
obj <- objective(type = "CV")
# Export the raw matrix and the objective, NOT the Java-backed chgeno object
clusterExport(cl = cl, varlist = c('obj', 'QUAC_raw_matrix'))

# ---- CALL COMMAND
# Using parSapplyLB with local Java instantiation
QUAC_genOptGenCovs_CV <- parSapplyLB(
  cl = cl,
  # Iterate over sample sizes from 2 to one less than the total number of samples
  1:(nrow(QUAC_raw_matrix) - 2),
  function(i) {
    # 1. Instantiate the Java object LOCALLY on the worker node
    local_chgeno <- corehunter::genotypes(
      QUAC_raw_matrix,
      format = "biparental"
    )
    # 2. Execute sampleCore using the local object
    result <- corehunter::sampleCore(
      local_chgeno,
      obj = obj,
      size = i + 1,
      steps = 10
    )
    # 3. Return just the CV value
    return(result$CV)
  }
)

# Close cluster
stopCluster(cl)

# #furrr implimentation
# # Load necessary packages (adding future and furrr)
# pacman::p_load(
#   adegenet,
#   terra,
#   parallel,
#   RColorBrewer,
#   viridis,
#   scales,
#   vcfR,
#   usedist,
#   corehunter,
#   future,
#   furrr
# )

# # Specif file path to input files
# GeoGenCorr_wd <- '/home/dune/trueNAS/work/GeographicGeneticCorrelation/'
# setwd(GeoGenCorr_wd)
# QUAC_filePath <- paste0(GeoGenCorr_wd, 'Datasets/QUAC/')

# # ---- PARALLELIZATION SETUP ----
# # Set up relevant cores (using future's availableCores() is generally safer)
# num_cores <- availableCores() - 12
# # 20 cores was maxing out CPU, and about 110gb ram

# # Tell the future package to resolve asynchronously via background R sessions
# plan(multisession, workers = num_cores)

# # ---- GENETIC MATRIX ----

# # Read in genind
# QUAC_genind <- read.genepop(paste0(
#   QUAC_filePath,
#   'Genetic/QUAC_REF_Wild_R80_NOMAF_1SNP_NoK.gen'
# ))
# # Subset QUAC genind object (to accelerate testing)
# # QUAC_genind <- QUAC_genind[1:5, , drop = TRUE]

# # Extract the base-R matrix to avoid Java serialization issues
# QUAC_raw_matrix <- QUAC_genind@tab

# # Specify allelic coverage CoreHunter objective
# obj <- objective(type = "CV")

# # ---- CALL COMMAND ----
# # Using future_map_dbl to iterate and return a numeric vector
# # furrr automatically identifies and exports QUAC_raw_matrix and obj to the workers
# QUAC_genOptGenCovs_CV <- future_map_dbl(
#   .x = 1:(nrow(QUAC_raw_matrix) - 2),
#   .f = function(i) {
#     # 1. Instantiate the Java object LOCALLY on the worker node
#     local_chgeno <- corehunter::genotypes(
#       QUAC_raw_matrix,
#       format = "biparental"
#     )
#     # 2. Execute sampleCore using the local object
#     result <- corehunter::sampleCore(
#       local_chgeno,
#       obj = obj,
#       size = i + 1,
#       steps = 10
#     )
#     # 3. Return just the CV value (future_map_dbl expects a numeric output)
#     return(result$CV)
#   },
#   # Good practice: ensures safe random number generation across parallel sessions
#   .options = furrr_options(seed = TRUE)
# )

# # Optional but good practice: shut down the background workers when finished
# plan(sequential)
