## ============================================================================
## Cluster Analysis of Weekly Entropy-Complexity Features (RSAM, WIZ, 2019)
## Steps: (1) correlation matrix of features, (2) PCA on the (scaled)
## features, (3) DBSCAN clustering of weeks in PCA space, (4) PCA biplot
## (weeks + feature loadings together, coloured by DBSCAN cluster).
## ============================================================================
##
## ASSUMPTION: DBSCAN clusters OBSERVATIONS (weeks), not features. The
## correlation matrix is a diagnostic on feature redundancy; PCA reduces
## the 13 correlated features to a handful of orthogonal components; DBSCAN
## then clusters weeks using those PCA scores.
##
## OTHER ASSUMPTIONS (edit CONFIG below if wrong):
##  1. Features standardised (mean 0, sd 1) before PCA and before DBSCAN --
##     the 13 features are on very different numeric scales.
##  2. Rows with any NA among the selected 13 features are dropped before
##     correlation/PCA/DBSCAN -- a warning is printed showing how many were
##     dropped.
##  3. dbscan_eps <- 1.8, read off the k-NN distance plot (k = 4) as the
##     point where the curve bends from a gradual rise to a steep climb
##     (around the 55th sorted point, y ~ 1.8-2.0). Since the bend is
##     fairly gradual rather than a sharp knee, it's worth re-running with
##     eps = 1.5 and eps = 2.2 as well and checking whether the number of
##     clusters/noise points stays stable across that range -- if it swings
##     a lot, that's worth flagging as a limitation of the clustering.
##  4. dbscan_minPts auto-set to 2 x n_PCs retained (Sander et al. 1998
##     heuristic: minPts >= dim + 1).
##  5. PCs retained for DBSCAN/biplot: enough to reach `pca_var_threshold`
##     (default 90%) cumulative variance, capped at a minimum of 2 PCs.
##  6. ALL plots use base_family/text family = "serif", and are print()'d
##     to the active graphics device before being saved with ggsave(), so
##     you see every figure interactively as the script runs, in addition
##     to the saved PDF.
## ============================================================================

library(readr)
library(tidyverse)
library(corrplot)
library(ggcorrplot)
library(dbscan)
library(factoextra)
library(ggrepel)

select <- dplyr::select   # guards against MASS::select() masking dplyr's,
# in case MASS was loaded earlier in this R
# session (e.g. from the LDA/MANOVA script) --
# this script doesn't load MASS itself, but base
# R doesn't "unload" packages between scripts
# within the same session


## ---- CONFIG ---------------------------------------------------------------
results_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"
out_dir     <- file.path(results_dir, "Cluster_Analysis_2019")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(results_dir, "weekly_entropy_complexity_rsam_avg_2019.csv")

# The 13 features you specified
feature_cols <- c(
  "H_Shannon", "C_Shannon", "Fisher_Info", "C_Fisher",
  "H_Renyi", "C_Renyi", "H_Tsallis", "C_Tsallis",
  "Var_H_Shannon", "Var_C_Shannon", "Var_H_Renyi", "Var_H_Tsallis",
  "Disequilibrium"
)

pca_var_threshold <- 0.90   # cumulative variance to retain for DBSCAN/biplot
min_pcs           <- 2      # always keep at least this many PCs

dbscan_minPts <- NA    # NA = auto-set to 2 x n_PCs_retained (see Step 3)
dbscan_eps    <- 1.8   # read off the k-NN plot elbow at k = 4 (see ASSUMPTION 3)

# Common serif theme addition -- appended to every ggplot object below so
# all figures (axis text, titles, legends) render in a consistent serif
# font for the paper.
serif_theme <- theme(
  text = element_text(family = "serif"),
  plot.title = element_text(family = "serif")
)

## ---- Load & select features ------------------------------------------------
df_raw <- read_csv(input_file, show_col_types = FALSE)

df <- df_raw %>%
  select(Week_Index, Week_Start, Week_End, all_of(feature_cols))

n_before <- nrow(df)
df <- df %>% drop_na(all_of(feature_cols))
n_after <- nrow(df)

if (n_after < n_before) {
  warning(sprintf(
    "%d of %d weeks dropped due to NA in one or more of the 13 features (likely edge-effect weeks). Remaining: %d weeks.",
    n_before - n_after, n_before, n_after
  ))
}

## ============================================================================
## STEP 1 -- Correlation matrix (diagnostic on feature redundancy)
## ============================================================================
feat_matrix <- df %>% select(all_of(feature_cols)) %>% as.matrix()

# Standardise: mean 0, sd 1 per feature (see ASSUMPTION 1)
feat_scaled <- scale(feat_matrix)

cor_mat <- cor(feat_scaled, use = "pairwise.complete.obs")

p_corr <- ggcorrplot(
  cor_mat,
  type = "lower",
  lab = TRUE,
  lab_size = 2.6,
  colors = c("#2c7bb6", "white", "#d7191c"),
  outline.color = "grey70",
  tl.cex = 8
) +
  labs(title = "Feature correlation matrix (RSAM, 2019)") +
  serif_theme

# The in-tile correlation coefficient labels (lab = TRUE) are a separate
# geom_text() layer added internally by ggcorrplot() with no `family`
# argument exposed -- theme(text = ...) above doesn't reach it, so it's
# patched directly on the layer.
text_layer <- which(sapply(p_corr$layers, function(l) inherits(l$geom, "GeomText")))
p_corr$layers[[text_layer]]$aes_params$family <- "serif"

print(p_corr)

ggsave(
  filename = file.path(out_dir, "correlation_matrix_rsam_2019.pdf"),
  plot = p_corr, width = 20, height = 18, units = "cm", device = "pdf"
)

# Print highly correlated pairs (|r| > 0.8) to the console
high_cor <- which(abs(cor_mat) > 0.8 & abs(cor_mat) < 1, arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  high_cor_pairs <- data.frame(
    Feature_1 = rownames(cor_mat)[high_cor[, 1]],
    Feature_2 = colnames(cor_mat)[high_cor[, 2]],
    r = cor_mat[high_cor]
  ) %>%
    filter(Feature_1 < Feature_2) %>%
    arrange(desc(abs(r)))
  message("Highly correlated feature pairs (|r| > 0.8):")
  print(high_cor_pairs)
}

## ============================================================================
## STEP 2 -- PCA on the scaled features
## ============================================================================
pca_fit <- prcomp(feat_matrix, center = TRUE, scale. = TRUE)

# ---- Scree plot -------------------------------------------------------
p_scree <- fviz_eig(
  pca_fit, addlabels = TRUE, barfill = "#2c7bb6", barcolor = "grey30"
) +
  labs(title = "PCA scree plot (RSAM, 2019)") +
  serif_theme

print(p_scree)

ggsave(
  filename = file.path(out_dir, "pca_scree_rsam_2019.pdf"),
  plot = p_scree, width = 16, height = 12, units = "cm", device = "pdf"
)

# ---- Variable loadings plot (PC1 vs PC2) -------------------------------
p_loadings <- fviz_pca_var(
  pca_fit, col.var = "contrib",
  gradient.cols = c("#2c7bb6", "#ffffbf", "#d7191c"),
  repel = TRUE
) +
  labs(title = "PCA variable loadings (PC1 vs PC2)") +
  serif_theme

print(p_loadings)

ggsave(
  filename = file.path(out_dir, "pca_loadings_rsam_2019.pdf"),
  plot = p_loadings, width = 16, height = 14, units = "cm", device = "pdf"
)

# ---- How many PCs to retain for DBSCAN/biplot (see ASSUMPTION 5) ------
var_explained     <- summary(pca_fit)$importance["Proportion of Variance", ]
cum_var_explained <- cumsum(var_explained)
n_pcs <- max(min_pcs, which(cum_var_explained >= pca_var_threshold)[1])

message(sprintf(
  "Retaining %d PCs for DBSCAN (%.1f%% cumulative variance explained).",
  n_pcs, cum_var_explained[n_pcs] * 100
))

pca_scores <- as.data.frame(pca_fit$x[, 1:n_pcs, drop = FALSE])
pca_scores$Week_Index <- df$Week_Index
pca_scores$Week_Start <- df$Week_Start

## ============================================================================
## STEP 3 -- DBSCAN on the PCA scores
## ============================================================================
dbscan_input <- as.matrix(pca_scores[, 1:n_pcs, drop = FALSE])

if (is.na(dbscan_minPts)) {
  dbscan_minPts <- 2 * n_pcs
  message("dbscan_minPts auto-set to ", dbscan_minPts, " (2 x n_PCs).")
}

# k-NN distance plot (base R graphics, not ggplot -- family set via par())
pdf(file.path(out_dir, "knn_distance_plot.pdf"), width = 7, height = 5, family = "serif")
par(family = "serif")
kNNdistplot(dbscan_input, k = dbscan_minPts)
abline(h = seq(0, 8, by = 0.5), col = "grey85", lty = 3)
abline(h = dbscan_eps, col = "red", lty = 2)  # marks the chosen eps for reference
title(main = sprintf("k-NN distance plot (k = %d) -- eps = %.2f", dbscan_minPts, dbscan_eps))
dev.off()

# Also draw it to the active/interactive device so you see it inline
par(family = "serif")
kNNdistplot(dbscan_input, k = dbscan_minPts)
abline(h = seq(0, 8, by = 0.5), col = "grey85", lty = 3)
abline(h = dbscan_eps, col = "red", lty = 2)
title(main = sprintf("k-NN distance plot (k = %d) -- eps = %.2f", dbscan_minPts, dbscan_eps))

if (is.na(dbscan_eps)) {
  stop(
    "dbscan_eps is not set. Inspect the k-NN distance plot just produced, ",
    "identify the 'elbow' (where distance sharply increases), set ",
    "dbscan_eps to that value in CONFIG, then re-run this script."
  )
}

db_fit <- dbscan(dbscan_input, eps = dbscan_eps, minPts = dbscan_minPts)

pca_scores$Cluster <- factor(db_fit$cluster)  # cluster 0 = DBSCAN "noise"

message(sprintf(
  "DBSCAN found %d cluster(s) plus %d noise point(s), out of %d weeks.",
  length(unique(db_fit$cluster[db_fit$cluster != 0])),
  sum(db_fit$cluster == 0),
  nrow(pca_scores)
))
## ---- Cluster plot in PC1-PC2 space, labelling ALL noise points -----------
## Instead of a hard-coded weeks_to_label vector, label every point DBSCAN
## assigned to cluster "0" (noise) -- this automatically includes whichever
## weeks are flagged as noise, whether or not you already knew about them
## (e.g. the extra unrecognised point you mentioned), rather than relying on
## a list you have to keep updating by hand.
label_df <- pca_scores %>% filter(Cluster == "0")

message(sprintf("Labelling %d noise point(s): weeks %s",
                nrow(label_df), paste(label_df$Week_Index, collapse = ", ")))

p_cluster <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_brewer(
    palette = "Set1",
    labels = function(x) ifelse(x == "0", "Noise", paste("Cluster", x))
  ) +
  geom_point(
    data = label_df, aes(x = PC1, y = PC2),
    inherit.aes = FALSE, shape = 21, size = 3.5, stroke = 1, color = "black", fill = NA
  ) +
  geom_text_repel(
    data = label_df, aes(x = PC1, y = PC2, label = paste0("week ", Week_Index)),
    inherit.aes = FALSE, color = "black", size = 3, fontface = "bold",
    min.segment.length = 0, seed = 42, family = "serif",
    max.overlaps = Inf   # ensures every noise label is drawn even if crowded --
    # without this, ggrepel can silently drop labels it
    # judges "too crowded to place," which is exactly how
    # you ended up with one unlabelled noise point before
  ) +
  labs(
    title = "DBSCAN clusters of weekly RSAM features (PCA space, 2019)",
    x = sprintf("PC1 (%.1f%% variance)", var_explained[1] * 100),
    y = sprintf("PC2 (%.1f%% variance)", var_explained[2] * 100)
  ) +
  theme_classic(base_size = 12, base_family = "serif") +
  theme(legend.position = "bottom") +
  serif_theme

print(p_cluster)
ggsave(
  filename = file.path(out_dir, "dbscan_clusters_rsam_2019.pdf"),
  plot = p_cluster, width = 18, height = 16, units = "cm", device = "pdf"
)

## ============================================================================
## STEP 4 -- PCA biplot (weeks + feature loadings together),
## coloured by DBSCAN cluster
## ============================================================================
## Combines the scree/loadings information with the DBSCAN result: weeks
## (points, coloured by cluster) AND the 13 feature loading vectors in the
## same PC1-PC2 space, so you can see which features pull which
## weeks/clusters in which direction.

p_biplot <- fviz_pca_biplot(
  pca_fit,
  axes = c(1, 2),
  geom.ind = "point",
  col.ind = pca_scores$Cluster,   # colour weeks by DBSCAN cluster -- relies
  # on pca_scores and pca_fit$x sharing row
  # order (true here since neither was
  # reordered/filtered after prcomp())
  palette = "Set1",
  addEllipses = FALSE,            # see note below before switching to TRUE
  col.var = "black",
  label = "var",                  # only label the 13 loading arrows, not
  # every individual week point
  repel = TRUE,
  pointsize = 2,
  legend.title = "Cluster"
) +
  labs(title = "PCA biplot: weeks (by DBSCAN cluster) and feature loadings") +
  serif_theme

print(p_biplot)

ggsave(
  filename = file.path(out_dir, "pca_biplot_clusters_rsam_2019.pdf"),
  plot = p_biplot, width = 18, height = 16, units = "cm", device = "pdf"
)
# Note: addEllipses = TRUE draws a confidence ellipse per cluster, but
# DBSCAN's "Cluster 0" (noise) usually isn't a real cluster -- an ellipse
# drawn around scattered noise points is misleading. Only switch this on
# if your clusters (excluding noise) each have a reasonable number of
# points; check with table(pca_scores$Cluster) first.

## ---- Save cluster assignments back to a CSV for further use --------------
write_csv(
  pca_scores %>% select(Week_Index, Week_Start, Cluster, everything()),
  file.path(out_dir, "dbscan_cluster_assignments_rsam_2019.csv")
)

message("Done. Outputs saved to: ", out_dir)