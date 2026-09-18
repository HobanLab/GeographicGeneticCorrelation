# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FIGURE 3 FROM THE NEW POINT SUMMARIES AND TABLE 3 OPTIMA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Reproduces the "BUILDING AND PLOTTING CORRELATION MATRICES" section of runPointSummaries.R
# (main branch, 2025-08-13) using the recalculated spatial metrics (runPointSummaries_carver.R).
# The resampling arrays which that script reads the optimal buffer sizes from are not tracked in
# the repository, so the optimal buffer sizes are taken from Table 3 of the manuscript (Draft 5)
# instead. The correlation, multiple testing correction and plotting code below is otherwise
# unchanged from runPointSummaries.R.

pacman::p_load(dplyr, readr, tibble, Hmisc, corrplot)

# INPUTS ----
# Point summaries: rows are datasets, columns are the metrics (the percent version of AOO is not used)
pointSummariesMat <- read_csv('Datasets/pointSummaryMeasures_carver.csv', show_col_types = FALSE) |>
  dplyr::select(-AOO_pct) |>
  tibble::column_to_rownames('taxon') |>
  as.matrix()
# Optimal buffer sizes (km) for each coverage type, from Table 3 of the manuscript (Draft 5).
# Column names as in runPointSummaries.R; ARTH has no SDM value
optBuffsMat <- tibble::tribble(
  ~taxon, ~`Opt. Geo. Buff`, ~`Opt. Geo. SDM`, ~`Opt. Eco.`,
  'AMTH',   0.5,  65,  15,
  'ARTH', 240,    NA,  45,
  'COGL',  10,    25,   0.5,
  'HIWA',  15,    20,   0.5,
  'MIGU', 230,   240, 250,
  'PICO', 230,   210, 110,
  'QUAC',   0.5,  15,  20,
  'QULO', 240,   190, 200,
  'VILA', 210,    90, 150,
  'YUBR',   4,    15,  25) |>
  tibble::column_to_rownames('taxon') |>
  as.matrix()
# Combine the point summary matrix to the optimal buffer size matrix. Rows are datasets and
# columns are summary metrics; order by optimal GeoBuff size (as in runPointSummaries.R)
SMBO_Mat <- cbind(pointSummariesMat, optBuffsMat[rownames(pointSummariesMat), ])
SMBO_Mat <- SMBO_Mat[order(SMBO_Mat[,'Opt. Geo. Buff']),]

# BUILDING AND PLOTTING CORRELATION MATRICES ----
# Setting the correlation type to Spearman, since we don't know whether the relationship
# between the optimal buffer sizes and each point statistic is linear (and data likely
# isn't Normally distributed)
corType <- 'spearman'
# GEO/ECO COVERAGES
# Build a correlation matrix based off of values
corMat_SMBO <- rcorr(SMBO_Mat, type=corType)
# Replace NAs in diagonal of p-value matrix with 0s, to match dimensions
corMat_SMBO$P[which(is.na(corMat_SMBO$P))] <- 0
# CORRECTION FOR MULTIPLE TESTING
# Get the unique (non-NA, off-diagonal) p-values
pValues <- corMat_SMBO$P[upper.tri(corMat_SMBO$P)]
# Adjust p-values (using BH, Benjamini–Hochberg: tests are not independent, we want to
# reduce false positives but retain power...)
adj_pValues <- p.adjust(pValues, method = "BH")
# Put adjusted p-values back into a symmetric matrix
adj_pMat <- matrix(NA, nrow = ncol(corMat_SMBO$P), ncol = ncol(corMat_SMBO$P))
rownames(adj_pMat) <- colnames(adj_pMat) <- colnames(corMat_SMBO$P)
adj_pMat[upper.tri(adj_pMat)] <- adj_pValues
adj_pMat[lower.tri(adj_pMat)] <- t(adj_pMat)[lower.tri(adj_pMat)]
diag(adj_pMat) <- 0  # Set diagonal to 0, for plotting purposes
corMat_SMBO$P <- adj_pMat # Reassign corrected p values in correlation matrix
# Plot correlation matrix using corrplot. Label significant correlations using
# asterisks
corrplot(corMat_SMBO$r, type="upper", order="original", p.mat = corMat_SMBO$P, 
         sig.level = 0.01, insig = "label_sig", diag = FALSE)

# Write image to disc
imageOutDir <- 'Datasets/figure3_pointSummaries_carver.png'
png(filename=imageOutDir, width=900, height=760)
par(oma=c(0,0,3,0), mar=c(5,4,7,2)+0.1)
corrplot(corMat_SMBO$r, type="upper", order="original", p.mat = corMat_SMBO$P, 
         sig.level = 0.01, insig = "label_sig", diag = FALSE, cl.cex = 1.4, tl.cex=1.2)
title('Correlations: Spatial statistics', line = 5.6, cex.main=1.3)
dev.off() # Turn off plotting device

# TABLE OF VALUES (addition, not in runPointSummaries.R) ----
# The correlation coefficients, sample sizes and BH adjusted p values for each metric by
# optimal buffer size pair, for the text of the manuscript
metrics <- colnames(pointSummariesMat)
optima <- colnames(optBuffsMat)
corTable <- expand.grid(metric = metrics, optimum = optima, stringsAsFactors = FALSE) |>
  dplyr::mutate(rho = round(mapply(function(m, o) corMat_SMBO$r[m, o], metric, optimum), 3),
                n = mapply(function(m, o) corMat_SMBO$n[m, o], metric, optimum),
                p_BH = round(mapply(function(m, o) corMat_SMBO$P[m, o], metric, optimum), 3)) |>
  dplyr::arrange(optimum, desc(abs(rho)))
print(corTable, row.names = FALSE)
write_csv(corTable, 'Datasets/figure3_correlations_carver.csv')
