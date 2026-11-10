## ============================================================================
## Clustering of Weekly H-C (Entropy-Complexity) Features -- RSAM, 2007-2019
## Step 1: Correlation plot (all 13 features)
## Step 2: K-means clustering on reduced 7-feature set (+ elbow plot)
## Step 3: DBSCAN clustering on reduced 7-feature set (+ kNN distance plot)
## Step 4: Export cluster results to Excel
## Step 5: Diagnostic -- interactively identify specific points (easylabel)
## Step 6: Load VAL ground-truth data and assign True_Label per week
## Step 7: Build confusion matrix (DBSCAN cluster 0/1 vs VAL Unrest/Not unrest)
## Step 8: Isolate and visualize misclassified (FN/FP) weeks
## All axis/title text uses the serif font family throughout.
## All outputs are saved to OUT_DIR (set below).
## ============================================================================

library(readr)
library(readxl)
library(here)
library(dplyr)
library(ggplot2)
library(ggcorrplot)
library(factoextra)   # kmeans/dbscan cluster viz (PCA-based), elbow/silhouette
library(cluster)
library(dbscan)       # dbscan() + kNNdistplot()
library(writexl)      # write_xlsx() -- export cluster results to Excel
library(easylabel)    # interactive diagnostic labelling of scatter points

## ---- Output / input directories -------------------------------------------
#OUT_DIR <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results/Weekly_OP_Analysis_2007to2019/Cluster Analysis_2007to2019"
#DATA_DIR <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"
OUT_DIR <- here("results", "Weekly_OP_Analysis_2007to2019", "Cluster Analysis_2007to2019")
DATA_DIR <- here("data")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## Shared serif theme, reused on every ggplot
serif_theme <- theme(
  text        = element_text(family = "serif"),
  axis.text   = element_text(family = "serif"),
  axis.title  = element_text(family = "serif"),
  plot.title  = element_text(family = "serif", face = "bold", hjust = 0.5),
  legend.text = element_text(family = "serif"),
  legend.title = element_text(family = "serif")
)

## ---- Load data --------------------------------------------------------
weekly_HC_rsam_avg_2007to2019 <- read_csv(
  "results/Weekly_OP_Analysis_2007to2019/weekly_entropy_complexity_rsam_avg_2007to2019.csv"
)
#View(weekly_HC_rsam_avg_2007to2019)

df <- data.frame(weekly_HC_rsam_avg_2007to2019)
#print(names(df))

## Full feature set -- used ONLY for the correlation plot, to decide
## which features to keep for clustering.
features1 <- c("H_Shannon", "C_Shannon",
               "Fisher_Info", "C_Fisher",
               "H_Renyi", "C_Renyi",
               "H_Tsallis", "C_Tsallis",
               "Var_H_Shannon", "Var_C_Shannon",
               "Var_H_Tsallis", "Var_H_Renyi",
               "Disequilibrium"
)

## Reduced feature set -- chosen after inspecting the correlation plot
## below. Used for ALL clustering (k-means + DBSCAN).
features <- c(
  "Fisher_Info", "C_Fisher",
  "Var_H_Shannon", "Var_C_Shannon",
  "Var_H_Tsallis", "Var_H_Renyi",
  "Disequilibrium"
)

## Keep only complete cases on the FULL feature set (features1), so the
## same rows/weeks are used consistently through correlation, clustering,
## and the Excel export. `keep_idx` lets cluster labels be joined back
## onto the original weeks (Week_Index/Week_Start/Week_End) later.
keep_idx  <- complete.cases(df[, features1])
df_feat   <- df[keep_idx, features1]   # 13 features -- correlation plot only
df_scaled <- scale(df_feat)            # 13 features, standardized

## Reduced 7-feature subset used for clustering. Standardize separately
## since clustering can't handle NAs and features are on very different
## scales (entropy/complexity vs. variance terms).
df_feat_clust   <- df[keep_idx, features]   # 7 features
df_scaled_clust <- scale(df_feat_clust)     # 7 features, standardized

## ============================================================================
## 1. CORRELATION PLOT (all 13 features)
## ============================================================================

corr_mat <- cor(df_feat, use = "pairwise.complete.obs", method = "spearman")

p_corr <- ggcorrplot(
  corr_mat,
  method = "square",
  type = "lower",
  lab = TRUE,
  lab_size = 3,
  colors = c("#3B4CC0", "white", "#B40426"),
  outline.color = "grey90"
) +
  labs(title = "Correlation of Entropy Complexity Features (RSAM, 2007-2019)") +
  serif_theme

print(p_corr)
ggsave(file.path(OUT_DIR, "correlation_plot_rsam.pdf"), p_corr, width = 8, height = 7)

## ============================================================================
## 2. K-MEANS CLUSTERING (reduced 7-feature set)
## ============================================================================

## Elbow plot (within-cluster sum of squares vs k) -- inspect this, pick
## the k where the curve bends, then set `k` below.
p_elbow <- fviz_nbclust(df_scaled_clust, kmeans, method = "wss", k.max = 10) +
  labs(title = "Elbow Method for Choosing k") +
  serif_theme
print(p_elbow)
ggsave(file.path(OUT_DIR, "elbow_plot_rsam.pdf"), p_elbow, width = 8, height = 7)

set.seed(012345)
k <- 2   # <-- set this from the elbow plot above
km_res <- kmeans(df_scaled_clust, centers = k, nstart = 25)

## fviz_cluster projects onto the first two principal components when
## there are >2 features, so this works directly with all 7 features.
p_kmeans <- fviz_cluster(
  km_res, data = df_scaled_clust,
  geom = "point", ellipse.type = "convex",
  palette = "jco", ggtheme = theme_minimal()
) +
  labs(title = paste0("K-means Clustering (k = ", k, ")")) +
  serif_theme
print(p_kmeans)
ggsave(file.path(OUT_DIR, "kmeans_cluster_plot_rsam.pdf"), p_kmeans, width = 8, height = 7)

## ============================================================================
## 3. DBSCAN CLUSTERING (reduced 7-feature set)
## ============================================================================

## kNN distance plot to help pick eps: look for the "knee" (sharp bend)
## in the sorted k-distances. minPts is commonly set to ~2 * ncol(data).
minPts_val <- 2 * ncol(df_scaled_clust)

pdf(file.path(OUT_DIR, "knn_distance_plot_rsam.pdf"), width = 8, height = 7)
par(family = "serif")
kNNdistplot(df_scaled_clust, k = minPts_val)
title(main = "k-NN Distance Plot (for choosing eps)", family = "serif")
abline(h = 1.0, lty = 2, col = "red")  # move this line to match the knee you see
dev.off()

par(family = "serif")
kNNdistplot(df_scaled_clust, k = minPts_val)
title(main = "k-NN Distance Plot (for choosing eps)", family = "serif")
abline(h = 1.0, lty = 2, col = "red")

eps_val <- 1.1   # <-- set this by reading off the knee in the plot above

db_res <- dbscan(df_scaled_clust, eps = eps_val, minPts = minPts_val)

p_dbscan <- fviz_cluster(
  list(data = df_scaled_clust, cluster = db_res$cluster),
  geom = "point", ellipse = TRUE,
  palette = "jco", ggtheme = theme_minimal()
) +
  labs(title = paste0("DBSCAN Clustering (eps = ", eps_val,
                      ", minPts = ", minPts_val, ")")) +
  serif_theme
print(p_dbscan)
ggsave(file.path(OUT_DIR, "dbscan_cluster_plot_rsam.pdf"), p_dbscan, width = 8, height = 7)

## Note: DBSCAN labels noise points as cluster "0" -- check
## table(db_res$cluster) to see how many points were flagged as noise.

## ============================================================================
## 4. EXPORT CLUSTER RESULTS TO EXCEL
## ============================================================================

## Any identifying columns present in the original data (e.g. week labels
## or dates) are carried along so each row in the Excel file is traceable
## back to a specific week -- add/remove names here to match your CSV.
id_cols     <- intersect(c("Week_Index", "Week_Start", "Week_End"), names(df))
df_ids_kept <- df[keep_idx, id_cols, drop = FALSE]

results_df <- cbind(
  df_ids_kept,
  df_feat,
  KMeans_Cluster = km_res$cluster,
  DBSCAN_Cluster = db_res$cluster
)

## Quick per-cluster feature means, handy as a summary sheet alongside
## the full row-level results.
kmeans_summary <- aggregate(df_feat, by = list(KMeans_Cluster = km_res$cluster), FUN = mean)
dbscan_summary <- aggregate(df_feat, by = list(DBSCAN_Cluster = db_res$cluster), FUN = mean)

write_xlsx(
  list(
    Cluster_Results = results_df,
    KMeans_Summary  = kmeans_summary,
    DBSCAN_Summary  = dbscan_summary
  ),
  path = file.path(OUT_DIR, "cluster_results_rsam.xlsx")
)

## ============================================================================
## 5. DIAGNOSTIC: identify specific points via easylabel
## ============================================================================
## Interactive check -- e.g. a point visually inside the "wrong" cluster
## (such as a DBSCAN noise point sitting among a dense cluster). Hover to
## see its week, click to label it, then close the app when done.

## Explicit PCA on the same standardized 7-feature data DBSCAN/fviz_cluster used
pca_res <- prcomp(df_scaled_clust)

plot_df <- data.frame(
  Dim1 = pca_res$x[, 1],
  Dim2 = pca_res$x[, 2],
  cluster = factor(db_res$cluster),
  Week_Start = df_ids_kept$Week_Start   # swap for Week_Index if Week_Start is absent
)

easylabel(
  plot_df,
  x = 'Dim1',
  y = 'Dim2',
  col = 'cluster',
  labs = 'Week_Start',
  size = 8
)

## ============================================================================
## 6. LOAD VAL GROUND-TRUTH DATA AND ASSIGN TRUE LABEL PER WEEK
## ============================================================================
## Ground truth: VAL = 1 ("Signs of volcano unrest" / "Minor volcanic unrest")
## is treated as the long-running BASELINE state -> "Not unrest".
## VAL >= 2 (eruptive activity, moderate/heightened unrest, eruption) is
## treated as an elevated/documented event -> "Unrest".
## (This split was confirmed after checking that VAL records cover the
## entire 2007-2019 window with no gaps, so "no VAL record at all" cannot
## be used to define the Not-unrest class.)

val_df <- read_csv(file.path(DATA_DIR, "VAL_WhakaariWhiteIsland.csv"))

val_df <- val_df %>%
  mutate(
    PeriodStart_UTC = as.POSIXct(`PeriodStart UTC`, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    PeriodEnd_UTC   = as.POSIXct(`PeriodEnd UTC`,   format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )

## Cap the open-ended 9999 row at "now" so date arithmetic doesn't break
open_ended <- is.na(val_df$PeriodEnd_UTC) | format(val_df$`PeriodEnd UTC`, "%Y") == "9999"
val_df$PeriodEnd_UTC[open_ended] <- Sys.time()

val_df$true_label <- ifelse(val_df$VAL == 1, "Not unrest", "Unrest")

## For each week window, find which VAL period(s) it overlaps and take the
## one with the LONGEST overlap as that week's ground truth.
assign_true_label <- function(w_start, w_end, val_df) {
  ov_start <- pmax(val_df$PeriodStart_UTC, w_start)
  ov_end   <- pmin(val_df$PeriodEnd_UTC, w_end)
  ov_secs  <- as.numeric(ov_end - ov_start, units = "secs")
  ov_secs[ov_secs < 0] <- 0
  if (all(ov_secs == 0)) return(NA_character_)
  val_df$true_label[which.max(ov_secs)]
}

results_df$True_Label <- mapply(
  assign_true_label,
  results_df$Week_Start,
  results_df$Week_End,
  MoreArgs = list(val_df = val_df)
)

## Sanity check both classes are populated
table(results_df$True_Label)

## ============================================================================
## 7. CONFUSION MATRIX: DBSCAN CLUSTER vs VAL GROUND TRUTH
## ============================================================================
## Cluster 0 (DBSCAN noise) -> predicted "Unrest"
## Any non-zero cluster     -> predicted "Not unrest"

results_df$Predicted_Label <- ifelse(results_df$DBSCAN_Cluster == 0, "Unrest", "Not unrest")

conf_df <- results_df %>% filter(!is.na(True_Label))

conf_mat <- table(
  Observation = factor(conf_df$True_Label, levels = c("Unrest", "Not unrest")),
  Predicted   = factor(conf_df$Predicted_Label, levels = c("Unrest", "Not unrest"),
                       labels = c("cluster 0", "cluster 1"))
)
print(conf_mat)

TP <- conf_mat["Unrest", "cluster 0"]
FN <- conf_mat["Unrest", "cluster 1"]
FP <- conf_mat["Not unrest", "cluster 0"]
TN <- conf_mat["Not unrest", "cluster 1"]

accuracy    <- (TP + TN) / sum(conf_mat)
precision   <- TP / (TP + FP)
recall      <- TP / (TP + FN)     # sensitivity -- how many true unrest weeks were caught
specificity <- TN / (TN + FP)
f1          <- 2 * precision * recall / (precision + recall)

cat(sprintf(
  "Accuracy: %.3f\nPrecision: %.3f\nRecall (Sensitivity): %.3f\nSpecificity: %.3f\nF1: %.3f\n",
  accuracy, precision, recall, specificity, f1
))

metrics_df <- data.frame(
  Metric = c("Accuracy", "Precision", "Recall (Sensitivity)", "Specificity", "F1"),
  Value  = c(accuracy, precision, recall, specificity, f1)
)

write_xlsx(
  list(
    Weekly_Results   = conf_df,
    Confusion_Matrix = as.data.frame.matrix(conf_mat),
    Metrics          = metrics_df
  ),
  path = file.path(OUT_DIR, "dbscan_confusion_matrix_rsam.xlsx")
)

## ============================================================================
## 8. ISOLATE AND VISUALIZE MISCLASSIFIED (FN / FP) WEEKS
## ============================================================================

## False Negatives: TRUE unrest, but DBSCAN put them in the dense cluster (1)
##   -> DBSCAN missed these
false_negatives <- conf_df %>%
  filter(True_Label == "Unrest", Predicted_Label == "Not unrest") %>%
  arrange(Week_Start)

## False Positives: TRUE not-unrest, but DBSCAN flagged them as noise (0)
##   -> DBSCAN over-flagged these
false_positives <- conf_df %>%
  filter(True_Label == "Not unrest", Predicted_Label == "Unrest") %>%
  arrange(Week_Start)

cat("False Negatives (missed unrest):", nrow(false_negatives), "\n")
cat("False Positives (over-flagged):", nrow(false_positives), "\n")

## Pull in the VAL description for context on why each week was labelled
## "Unrest" or "Not unrest" in the first place
mismatches <- bind_rows(
  false_negatives %>% mutate(Error_Type = "False Negative"),
  false_positives %>% mutate(Error_Type = "False Positive")
) %>%
  rowwise() %>%
  mutate(
    VAL_Description = val_df$Description[
      which.max(pmin(val_df$PeriodEnd_UTC, Week_End) -
                  pmax(val_df$PeriodStart_UTC, Week_Start))
    ]
  ) %>%
  ungroup() %>%
  arrange(Week_Start)

print(mismatches[, c("Week_Start", "Week_End", "True_Label",
                     "Predicted_Label", "Error_Type", "VAL_Description")])

write_xlsx(
  list(Mismatches = mismatches),
  path = file.path(OUT_DIR, "dbscan_mismatches_rsam.xlsx")
)

## ---- Diagnostic plot: mismatches on the PCA-projected cluster space -------

plot_df2 <- data.frame(
  Dim1 = pca_res$x[, 1],
  Dim2 = pca_res$x[, 2],
  DBSCAN_Cluster = factor(db_res$cluster),
  Week_Start = df_ids_kept$Week_Start
) %>%
  left_join(
    conf_df %>% select(Week_Start, True_Label, Predicted_Label),
    by = "Week_Start"
  ) %>%
  mutate(
    Match_Status = case_when(
      is.na(True_Label)                                    ~ "No VAL data",
      True_Label == "Unrest" & Predicted_Label == "Unrest"  ~ "TP",
      True_Label == "Not unrest" & Predicted_Label == "Not unrest" ~ "TN",
      True_Label == "Unrest" & Predicted_Label == "Not unrest" ~ "FN",
      True_Label == "Not unrest" & Predicted_Label == "Unrest" ~ "FP"
    )
  )

## Define the color mapping once, reuse everywhere
match_colors <- c(
  TP = "#2ca02c", TN = "grey70",
  FN = "#d62728", FP = "#ff7f0e",
  `No VAL data` = "grey90"
)

p_mismatch <- ggplot(plot_df2, aes(Dim1, Dim2, color = Match_Status)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = match_colors) +
  labs(title = "DBSCAN Classification vs VAL Ground Truth",
       color = "Result") +
  serif_theme
print(p_mismatch)
ggsave(file.path(OUT_DIR, "dbscan_mismatch_diagnostic.pdf"), p_mismatch, width = 8, height = 7)

## Interactively click through the FN/FP points specifically, using the SAME colors
easylabel(
  plot_df2 %>% filter(Match_Status %in% c("FN", "FP")),
  x = 'Dim1', y = 'Dim2', col = 'Match_Status',
  labs = 'Week_Start', size = 8,
  colScheme = match_colors[c("FN", "FP")]
)

## Interactively click through the TP/TN points specifically, using the SAME colors
easylabel(
  plot_df2 %>% filter(Match_Status %in% c("TN", "TP")),
  x = 'Dim1', y = 'Dim2', col = 'Match_Status',
  labs = 'Week_Start', size = 8,
  colScheme = match_colors[c("TN", "TP")]
)




########################################################################################
############################################################################################
# Same analysis with 13 features 
## ============================================================================
## Clustering of Weekly H-C (Entropy-Complexity) Features -- RSAM, 2007-2019
## FULL 13-FEATURE VERSION
## Step 1: Correlation plot (all 13 features)
## Step 2: K-means clustering on FULL 13-feature set (+ elbow plot)
## Step 3: DBSCAN clustering on FULL 13-feature set (+ kNN distance plot)
## Step 4: Export cluster results to Excel
## Step 5: Diagnostic -- interactively identify specific points (easylabel)
## Step 6: Load VAL ground-truth data and assign True_Label per week
## Step 7: Build confusion matrix (DBSCAN cluster 0/1 vs VAL Unrest/Not unrest)
## Step 8: Isolate and visualize misclassified (FN/FP) weeks
## All axis/title text uses the serif font family throughout.
## All outputs are saved to OUT_DIR (set below), with a "_13features" suffix
## so they sit alongside the 7-feature run's outputs without overwriting them.
## ============================================================================

library(readr)
library(readxl)
library(here)
library(dplyr)
library(ggplot2)
library(ggcorrplot)
library(factoextra)   # kmeans/dbscan cluster viz (PCA-based), elbow/silhouette
library(cluster)
library(dbscan)       # dbscan() + kNNdistplot()
library(writexl)      # write_xlsx() -- export cluster results to Excel
library(easylabel)    # interactive diagnostic labelling of scatter points

## ---- Output / input directories -------------------------------------------
OUT_DIR <- here("results", "Weekly_OP_Analysis_2007to2019", "Cluster Analysis_2007to2019")
DATA_DIR <- here("data")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## Shared serif theme, reused on every ggplot
serif_theme <- theme(
  text        = element_text(family = "serif"),
  axis.text   = element_text(family = "serif"),
  axis.title  = element_text(family = "serif"),
  plot.title  = element_text(family = "serif", face = "bold", hjust = 0.5),
  legend.text = element_text(family = "serif"),
  legend.title = element_text(family = "serif")
)

## ---- Load data --------------------------------------------------------
weekly_HC_rsam_avg_2007to2019 <- read_csv(
  here("results", "Weekly_OP_Analysis_2007to2019", "weekly_entropy_complexity_rsam_avg_2007to2019.csv")
)

df <- data.frame(weekly_HC_rsam_avg_2007to2019)

## Full 13-feature set -- used for BOTH the correlation plot AND all
## clustering in this version of the script (no reduced subset here).
features1 <- c("H_Shannon", "C_Shannon",
               "Fisher_Info", "C_Fisher",
               "H_Renyi", "C_Renyi",
               "H_Tsallis", "C_Tsallis",
               "Var_H_Shannon", "Var_C_Shannon",
               "Var_H_Tsallis", "Var_H_Renyi",
               "Disequilibrium"
)

## Keep only complete cases on the full feature set, so the same rows/weeks
## are used consistently through correlation, clustering, and export.
keep_idx  <- complete.cases(df[, features1])
df_feat   <- df[keep_idx, features1]   # 13 features
df_scaled <- scale(df_feat)            # 13 features, standardized

## No reduced subset in this version -- clustering uses the full 13-feature
## standardized matrix directly.
df_feat_clust   <- df_feat
df_scaled_clust <- df_scaled

## ============================================================================
## 1. CORRELATION PLOT (all 13 features)
## ============================================================================

corr_mat <- cor(df_feat, use = "pairwise.complete.obs", method = "spearman")

p_corr <- ggcorrplot(
  corr_mat,
  method = "square",
  type = "lower",
  lab = TRUE,
  lab_size = 3,
  colors = c("#3B4CC0", "white", "#B40426"),
  outline.color = "grey90"
) +
  labs(title = "Correlation of Entropy Complexity Features (RSAM, 2007-2019)") +
  serif_theme

print(p_corr)
ggsave(file.path(OUT_DIR, "correlation_plot_rsam.pdf"), p_corr, width = 8, height = 7)

## ============================================================================
## 2. K-MEANS CLUSTERING (full 13-feature set)
## ============================================================================

p_elbow <- fviz_nbclust(df_scaled_clust, kmeans, method = "wss", k.max = 10) +
  labs(title = "Elbow Method for Choosing k (13 Features)") +
  serif_theme
print(p_elbow)
ggsave(file.path(OUT_DIR, "elbow_plot_rsam_13features.pdf"), p_elbow, width = 8, height = 7)

set.seed(012345)
k <- 2   # <-- set this from the elbow plot above
km_res <- kmeans(df_scaled_clust, centers = k, nstart = 25)

p_kmeans <- fviz_cluster(
  km_res, data = df_scaled_clust,
  geom = "point", ellipse.type = "convex",
  palette = "jco", ggtheme = theme_minimal()
) +
  labs(title = paste0("K-means Clustering, 13 Features (k = ", k, ")")) +
  serif_theme
print(p_kmeans)
ggsave(file.path(OUT_DIR, "kmeans_cluster_plot_rsam_13features.pdf"), p_kmeans, width = 8, height = 7)

## ============================================================================
## 3. DBSCAN CLUSTERING (full 13-feature set)
## ============================================================================

## minPts = 2 * ncol(data) -- now 2*13 = 26
minPts_val <- 2 * ncol(df_scaled_clust)

pdf(file.path(OUT_DIR, "knn_distance_plot_rsam_13features.pdf"), width = 8, height = 7)
par(family = "serif")
kNNdistplot(df_scaled_clust, k = minPts_val)
title(main = "k-NN Distance Plot, 13 Features (for choosing eps)", family = "serif")
abline(h = 1.0, lty = 2, col = "red")  # move this line to match the knee you see
dev.off()

par(family = "serif")
kNNdistplot(df_scaled_clust, k = minPts_val)
title(main = "k-NN Distance Plot, 13 Features (for choosing eps)", family = "serif")
abline(h = 1.0, lty = 2, col = "red")

eps_val <- 1   # <-- set this by reading off the knee in the plot above

db_res <- dbscan(df_scaled_clust, eps = eps_val, minPts = minPts_val)

p_dbscan <- fviz_cluster(
  list(data = df_scaled_clust, cluster = db_res$cluster),
  geom = "point", ellipse = TRUE,
  palette = "jco", ggtheme = theme_minimal()
) +
  labs(title = paste0("DBSCAN Clustering, 13 Features (eps = ", eps_val,
                      ", minPts = ", minPts_val, ")")) +
  serif_theme
print(p_dbscan)
ggsave(file.path(OUT_DIR, "dbscan_cluster_plot_rsam_13features.pdf"), p_dbscan, width = 8, height = 7)

## Note: DBSCAN labels noise points as cluster "0" -- check
## table(db_res$cluster) to see how many points were flagged as noise.

## ============================================================================
## 4. EXPORT CLUSTER RESULTS TO EXCEL
## ============================================================================

id_cols     <- intersect(c("Week_Index", "Week_Start", "Week_End"), names(df))
df_ids_kept <- df[keep_idx, id_cols, drop = FALSE]

results_df <- cbind(
  df_ids_kept,
  df_feat,
  KMeans_Cluster = km_res$cluster,
  DBSCAN_Cluster = db_res$cluster
)

kmeans_summary <- aggregate(df_feat, by = list(KMeans_Cluster = km_res$cluster), FUN = mean)
dbscan_summary <- aggregate(df_feat, by = list(DBSCAN_Cluster = db_res$cluster), FUN = mean)

write_xlsx(
  list(
    Cluster_Results = results_df,
    KMeans_Summary  = kmeans_summary,
    DBSCAN_Summary  = dbscan_summary
  ),
  path = file.path(OUT_DIR, "cluster_results_rsam_13features.xlsx")
)

## ============================================================================
## 5. DIAGNOSTIC: identify specific points via easylabel
## ============================================================================

## Explicit PCA on the same standardized 13-feature data DBSCAN/fviz_cluster used
pca_res <- prcomp(df_scaled_clust)

plot_df <- data.frame(
  Dim1 = pca_res$x[, 1],
  Dim2 = pca_res$x[, 2],
  cluster = factor(db_res$cluster),
  Week_Start = df_ids_kept$Week_Start
)

easylabel(
  plot_df,
  x = 'Dim1',
  y = 'Dim2',
  col = 'cluster',
  labs = 'Week_Start',
  size = 8
)

## ============================================================================
## 6. LOAD VAL GROUND-TRUTH DATA AND ASSIGN TRUE LABEL PER WEEK
## ============================================================================

val_df <- read_csv(file.path(DATA_DIR, "VAL_WhakaariWhiteIsland.csv"))

val_df <- val_df %>%
  mutate(
    PeriodStart_UTC = as.POSIXct(`PeriodStart UTC`, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    PeriodEnd_UTC   = as.POSIXct(`PeriodEnd UTC`,   format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )

open_ended <- is.na(val_df$PeriodEnd_UTC) | format(val_df$`PeriodEnd UTC`, "%Y") == "9999"
val_df$PeriodEnd_UTC[open_ended] <- Sys.time()

val_df$true_label <- ifelse(val_df$VAL == 1, "Not unrest", "Unrest")

assign_true_label <- function(w_start, w_end, val_df) {
  ov_start <- pmax(val_df$PeriodStart_UTC, w_start)
  ov_end   <- pmin(val_df$PeriodEnd_UTC, w_end)
  ov_secs  <- as.numeric(ov_end - ov_start, units = "secs")
  ov_secs[ov_secs < 0] <- 0
  if (all(ov_secs == 0)) return(NA_character_)
  val_df$true_label[which.max(ov_secs)]
}

results_df$True_Label <- mapply(
  assign_true_label,
  results_df$Week_Start,
  results_df$Week_End,
  MoreArgs = list(val_df = val_df)
)

table(results_df$True_Label)

## ============================================================================
## 7. CONFUSION MATRIX: DBSCAN CLUSTER vs VAL GROUND TRUTH (13 features)
## ============================================================================

results_df$Predicted_Label <- ifelse(results_df$DBSCAN_Cluster == 0, "Unrest", "Not unrest")

conf_df <- results_df %>% filter(!is.na(True_Label))

conf_mat <- table(
  Observation = factor(conf_df$True_Label, levels = c("Unrest", "Not unrest")),
  Predicted   = factor(conf_df$Predicted_Label, levels = c("Unrest", "Not unrest"),
                       labels = c("cluster 0", "cluster 1"))
)
print(conf_mat)

TP <- conf_mat["Unrest", "cluster 0"]
FN <- conf_mat["Unrest", "cluster 1"]
FP <- conf_mat["Not unrest", "cluster 0"]
TN <- conf_mat["Not unrest", "cluster 1"]

accuracy    <- (TP + TN) / sum(conf_mat)
precision   <- TP / (TP + FP)
recall      <- TP / (TP + FN)
specificity <- TN / (TN + FP)
f1          <- 2 * precision * recall / (precision + recall)

cat(sprintf(
  "13-feature results (eps=%.2f, minPts=%d)\nAccuracy: %.3f\nPrecision: %.3f\nRecall (Sensitivity): %.3f\nSpecificity: %.3f\nF1: %.3f\n",
  eps_val, minPts_val, accuracy, precision, recall, specificity, f1
))

metrics_df <- data.frame(
  Metric = c("Accuracy", "Precision", "Recall (Sensitivity)", "Specificity", "F1"),
  Value  = c(accuracy, precision, recall, specificity, f1)
)

write_xlsx(
  list(
    Weekly_Results   = conf_df,
    Confusion_Matrix = as.data.frame.matrix(conf_mat),
    Metrics          = metrics_df
  ),
  path = file.path(OUT_DIR, "dbscan_confusion_matrix_rsam_13features.xlsx")
)

## ============================================================================
## 8. ISOLATE AND VISUALIZE MISCLASSIFIED (FN / FP) WEEKS
## ============================================================================

false_negatives <- conf_df %>%
  filter(True_Label == "Unrest", Predicted_Label == "Not unrest") %>%
  arrange(Week_Start)

false_positives <- conf_df %>%
  filter(True_Label == "Not unrest", Predicted_Label == "Unrest") %>%
  arrange(Week_Start)

cat("False Negatives (missed unrest):", nrow(false_negatives), "\n")
cat("False Positives (over-flagged):", nrow(false_positives), "\n")

mismatches <- bind_rows(
  false_negatives %>% mutate(Error_Type = "False Negative"),
  false_positives %>% mutate(Error_Type = "False Positive")
) %>%
  rowwise() %>%
  mutate(
    VAL_Description = val_df$Description[
      which.max(pmin(val_df$PeriodEnd_UTC, Week_End) -
                  pmax(val_df$PeriodStart_UTC, Week_Start))
    ]
  ) %>%
  ungroup() %>%
  arrange(Week_Start)

print(mismatches[, c("Week_Start", "Week_End", "True_Label",
                     "Predicted_Label", "Error_Type", "VAL_Description")])

write_xlsx(
  list(Mismatches = mismatches),
  path = file.path(OUT_DIR, "dbscan_mismatches_rsam_13features.xlsx")
)

## ---- True Positives: weeks correctly flagged as Unrest --------------------

true_positives <- conf_df %>%
  filter(True_Label == "Unrest", Predicted_Label == "Unrest") %>%
  arrange(Week_Start) %>%
  rowwise() %>%
  mutate(
    VAL_Description = val_df$Description[
      which.max(pmin(val_df$PeriodEnd_UTC, Week_End) -
                  pmax(val_df$PeriodStart_UTC, Week_Start))
    ]
  ) %>%
  ungroup()

cat("True Positives (correctly flagged unrest weeks):", nrow(true_positives), "\n")
print(true_positives[, c("Week_Start", "Week_End", "VAL_Description")])

write_xlsx(
  list(True_Positives = true_positives),
  path = file.path(OUT_DIR, "dbscan_true_positives_rsam_13features.xlsx")
)

## ---- Diagnostic plot: mismatches on the PCA-projected cluster space -------

plot_df2 <- data.frame(
  Dim1 = pca_res$x[, 1],
  Dim2 = pca_res$x[, 2],
  DBSCAN_Cluster = factor(db_res$cluster),
  Week_Start = df_ids_kept$Week_Start
) %>%
  left_join(
    conf_df %>% select(Week_Start, True_Label, Predicted_Label),
    by = "Week_Start"
  ) %>%
  mutate(
    Match_Status = case_when(
      is.na(True_Label)                                    ~ "No VAL data",
      True_Label == "Unrest" & Predicted_Label == "Unrest"  ~ "TP",
      True_Label == "Not unrest" & Predicted_Label == "Not unrest" ~ "TN",
      True_Label == "Unrest" & Predicted_Label == "Not unrest" ~ "FN",
      True_Label == "Not unrest" & Predicted_Label == "Unrest" ~ "FP"
    )
  )

match_colors <- c(
  TP = "#2ca02c", TN = "grey70",
  FN = "#d62728", FP = "#ff7f0e",
  `No VAL data` = "grey90"
)

p_mismatch <- ggplot(plot_df2, aes(Dim1, Dim2, color = Match_Status)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = match_colors) +
  labs(title = "DBSCAN Classification vs VAL Ground Truth (13 Features)",
       color = "Result") +
  serif_theme
print(p_mismatch)
ggsave(file.path(OUT_DIR, "dbscan_mismatch_diagnostic_13features.pdf"), p_mismatch, width = 8, height = 7)

easylabel(
  plot_df2 %>% filter(Match_Status %in% c("FN", "FP")),
  x = 'Dim1', y = 'Dim2', col = 'Match_Status',
  labs = 'Week_Start', size = 8,
  colScheme = match_colors[c("FN", "FP")]
)

easylabel(
  plot_df2 %>% filter(Match_Status %in% c("TN", "TP")),
  x = 'Dim1', y = 'Dim2', col = 'Match_Status',
  labs = 'Week_Start', size = 8,
  colScheme = match_colors[c("TN", "TP")]
)


#Final version of cluster analysis based on the experience from different 
# epsilon values from DBSCAN clustering. In this code block I analyzed DBSCAN clustering 
#with epsilon 1.1 and SVM with different kernel functions and cost parameters are 
#considered.
###################################################################
library(readr)
library(dplyr)
library(ggplot2)
library(ggcorrplot)
library(dbscan)
library(factoextra)
library(easylabel)
library(caret)
library(e1071)
library(writexl)

serif_theme <- theme(
  text = element_text(family = "serif"),
  axis.text = element_text(family = "serif"),
  axis.title = element_text(family = "serif"),
  plot.title = element_text(
    family = "serif",
    face = "bold",
    hjust = 0.5
  ),
  legend.text = element_text(family = "serif"),
  legend.title = element_text(family = "serif")
)

OUT_DIR <- "results/Weekly_OP_Analysis_2007to2019/Cluster Analysis_2007to2019"
DATA_FILE <- "results/Weekly_OP_Analysis_2007to2019/weekly_entropy_complexity_rsam_avg_2007to2019.csv"
VAL_FILE <- "data/VAL_WhakaariWhiteIsland.csv"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(DATA_FILE, show_col_types = FALSE)

features1 <- c(
  "H_Shannon", "C_Shannon", "Fisher_Info", "C_Fisher",
  "H_Renyi", "C_Renyi", "H_Tsallis", "C_Tsallis",
  "Var_H_Shannon", "Var_C_Shannon", "Var_H_Tsallis",
  "Var_H_Renyi", "Disequilibrium"
)

# Features retained from the previous analysis
features <- c(
  "Fisher_Info", "C_Fisher", "Var_H_Shannon",
  "Var_C_Shannon", "Var_H_Tsallis", "Var_H_Renyi",
  "Disequilibrium"
)

keep <- complete.cases(df[, features1])
analysis_df <- df[keep, ]
X <- analysis_df[, features]
X_scaled <- scale(X)

# Correlation analysis using all features
corr_mat <- cor(analysis_df[, features1], method = "spearman")
print(round(corr_mat, 2))

p_corr <- ggcorrplot(
  corr_mat,
  type = "lower",
  lab = TRUE,
  lab_size = 2,
  colors = c("#3B4CC0", "white", "#B40426")
) +
  labs(title = "Correlation of Entropy-Complexity Features") +
  serif_theme

print(p_corr)
ggsave(file.path(OUT_DIR, "correlation_plot.pdf"), p_corr, width = 9, height = 8)

# Display pairs with an absolute correlation greater than 0.80
high_cor <- which(abs(corr_mat) > 0.80 & upper.tri(corr_mat), arr.ind = TRUE)
high_cor_pairs <- data.frame(
  Feature_1 = rownames(corr_mat)[high_cor[, 1]],
  Feature_2 = colnames(corr_mat)[high_cor[, 2]],
  Correlation = corr_mat[high_cor]
)
print(high_cor_pairs)

# DBSCAN
minPts_val <- 2 * ncol(X_scaled)
kNNdistplot(X_scaled, k = minPts_val)
title(main = "k-NN Distance Plot", family = "serif")
abline(h = 1.1, col = "red", lty = 2)

eps_val <- 1.1  # Change this after examining the kNN-distance plot
db_res <- dbscan(X_scaled, eps = eps_val, minPts = minPts_val)
print(table(db_res$cluster))

p_dbscan <- fviz_cluster(
  list(data = X_scaled, cluster = db_res$cluster),
  geom = "point",
  ellipse = TRUE,
  ggtheme = theme_minimal(base_family = "serif")
) +
  labs(title = "DBSCAN Clustering") +
  serif_theme

print(p_dbscan)
ggsave(file.path(OUT_DIR, "dbscan_plot.pdf"), p_dbscan, width = 8, height = 7)

# PCA coordinates are used only to display and label the DBSCAN points
pca_res <- prcomp(X_scaled)
plot_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Cluster = factor(db_res$cluster),
  Week_Start = analysis_df$Week_Start
)

easylabel(
  plot_df,
  x = "PC1", y = "PC2", col = "Cluster",
  labs = "Week_Start", size = 8
)

# Assign the VAL true label to each week
analysis_df <- analysis_df %>%
  mutate(
    Week_Start = as.POSIXct(Week_Start, tz = "UTC"),
    Week_End = as.POSIXct(Week_End, tz = "UTC")
  )

val_df <- read_csv(VAL_FILE, show_col_types = FALSE) %>%
  mutate(
    PeriodStart_UTC = as.POSIXct(
      `PeriodStart UTC`, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"
    ),
    PeriodEnd_UTC = as.POSIXct(
      `PeriodEnd UTC`, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"
    ),
    PeriodEnd_UTC = if_else(
      is.na(PeriodEnd_UTC) | grepl("^9999", `PeriodEnd UTC`),
      Sys.time(), PeriodEnd_UTC
    ),
    True_Label = if_else(VAL == 1, "Not unrest", "Unrest")
  )

assign_true_label <- function(week_start, week_end) {
  overlap <- as.numeric(
    pmin(val_df$PeriodEnd_UTC, week_end) -
      pmax(val_df$PeriodStart_UTC, week_start),
    units = "secs"
  )
  overlap[is.na(overlap) | overlap < 0] <- 0
  if (all(overlap == 0)) return(NA_character_)
  val_df$True_Label[which.max(overlap)]
}

analysis_df$True_Label <- mapply(
  assign_true_label,
  analysis_df$Week_Start,
  analysis_df$Week_End
)

results_df <- analysis_df %>%
  mutate(
    DBSCAN_Cluster = db_res$cluster,
    DBSCAN_Label = if_else(DBSCAN_Cluster == 0, "Unrest", "Not unrest")
  ) %>%
  filter(!is.na(True_Label))

# DBSCAN confusion matrix: cluster 0 is treated as an anomaly/unrest prediction
dbscan_cm <- confusionMatrix(
  factor(results_df$DBSCAN_Label, levels = c("Unrest", "Not unrest")),
  factor(results_df$True_Label, levels = c("Unrest", "Not unrest")),
  positive = "Unrest",
  mode = "everything"
)
print(dbscan_cm)

# ---------------------------------------------------------------------------
# SVM kernel/cost/gamma selection under severe class imbalance
#
# With ~20 positives out of ~792 rows, a single 70/30 (or 80/20) holdout
# split is unreliable: the test set may end up with very few anomalies,
# making sensitivity/F1 estimates noisy and sometimes undefined. Instead we
# use repeated stratified k-fold CV, which keeps the class ratio consistent
# in every fold and averages performance over many resamples so every
# anomaly contributes to both training and evaluation across repeats.
#
# Model selection uses Balanced Accuracy and F1 rather than raw accuracy,
# since with this imbalance a trivial "always majority class" model scores
# ~97% accuracy while missing every anomaly.
# ---------------------------------------------------------------------------

model_df <- results_df %>%
  select(all_of(features), True_Label) %>%
  mutate(True_Label = factor(True_Label, levels = c("Not unrest", "Unrest")))

cat("\nClass balance in modeling data:\n")
print(table(model_df$True_Label))

# ---- Hyperparameter grid across kernels ----
costs   <- c(0.1, 1, 10, 100)
gammas  <- c(0.01, 0.05, 0.1, 0.5, 1)
degrees <- c(2, 3)

grid_linear <- expand.grid(
  kernel = "linear", cost = costs, gamma = NA_real_, degree = NA_real_
)
grid_radial <- expand.grid(
  kernel = "radial", cost = costs, gamma = gammas, degree = NA_real_
)
grid_sigmoid <- expand.grid(
  kernel = "sigmoid", cost = costs, gamma = gammas, degree = NA_real_
)
grid_poly <- expand.grid(
  kernel = "polynomial", cost = costs, gamma = gammas, degree = degrees
)

param_grid <- bind_rows(grid_linear, grid_radial, grid_sigmoid, grid_poly)
cat(sprintf("\nTotal hyperparameter combinations to evaluate: %d\n", nrow(param_grid)))

# ---- Repeated stratified k-fold CV settings ----
n_folds   <- 5   # small enough that each fold still contains several anomalies
n_repeats <- 20  # averages out noise from the tiny minority class

evaluate_params <- function(kernel, cost, gamma, degree, data,
                            n_folds = 5, n_repeats = 20) {
  bal_acc_rep <- numeric(n_repeats)
  f1_rep      <- numeric(n_repeats)
  sens_rep    <- numeric(n_repeats)
  spec_rep    <- numeric(n_repeats)
  
  for (r in seq_len(n_repeats)) {
    set.seed(1000 + r)
    folds <- createFolds(data$True_Label, k = n_folds, list = TRUE)
    
    fold_bal_acc <- rep(NA_real_, n_folds)
    fold_f1      <- rep(NA_real_, n_folds)
    fold_sens    <- rep(NA_real_, n_folds)
    fold_spec    <- rep(NA_real_, n_folds)
    
    for (f in seq_len(n_folds)) {
      test_idx   <- folds[[f]]
      train_fold <- data[-test_idx, ]
      test_fold  <- data[test_idx, ]
      
      # Skip degenerate folds where a class is entirely absent from training
      if (length(unique(train_fold$True_Label)) < 2) next
      if (length(unique(test_fold$True_Label)) < 2) {
        # still usable for prediction, confusionMatrix handles single-class
        # truth as long as both levels are declared
      }
      
      class_counts_fold <- table(train_fold$True_Label)
      cw <- sum(class_counts_fold) / (length(class_counts_fold) * class_counts_fold)
      
      svm_args <- list(
        formula = True_Label ~ .,
        data = train_fold,
        kernel = kernel,
        cost = cost,
        scale = TRUE,
        class.weights = cw
      )
      if (!is.na(gamma))  svm_args$gamma  <- gamma
      if (!is.na(degree)) svm_args$degree <- degree
      
      fit <- tryCatch(do.call(svm, svm_args), error = function(e) NULL)
      if (is.null(fit)) next
      
      pred <- tryCatch(predict(fit, test_fold), error = function(e) NULL)
      if (is.null(pred)) next
      
      cm <- tryCatch(
        confusionMatrix(
          factor(pred, levels = c("Not unrest", "Unrest")),
          factor(test_fold$True_Label, levels = c("Not unrest", "Unrest")),
          positive = "Unrest",
          mode = "everything"
        ),
        error = function(e) NULL
      )
      if (is.null(cm)) next
      
      fold_bal_acc[f] <- cm$byClass["Balanced Accuracy"]
      fold_f1[f]      <- ifelse(is.na(cm$byClass["F1"]), 0, cm$byClass["F1"])
      fold_sens[f]    <- cm$byClass["Sensitivity"]
      fold_spec[f]    <- cm$byClass["Specificity"]
    }
    
    bal_acc_rep[r] <- mean(fold_bal_acc, na.rm = TRUE)
    f1_rep[r]      <- mean(fold_f1, na.rm = TRUE)
    sens_rep[r]    <- mean(fold_sens, na.rm = TRUE)
    spec_rep[r]    <- mean(fold_spec, na.rm = TRUE)
  }
  
  data.frame(
    kernel = kernel, cost = cost, gamma = gamma, degree = degree,
    Mean_Balanced_Accuracy = mean(bal_acc_rep, na.rm = TRUE),
    SD_Balanced_Accuracy   = sd(bal_acc_rep, na.rm = TRUE),
    Mean_F1                = mean(f1_rep, na.rm = TRUE),
    SD_F1                  = sd(f1_rep, na.rm = TRUE),
    Mean_Sensitivity       = mean(sens_rep, na.rm = TRUE),
    Mean_Specificity       = mean(spec_rep, na.rm = TRUE)
  )
}

cat("\nRunning repeated stratified CV grid search (this may take a while)...\n")

results_grid <- do.call(rbind, lapply(seq_len(nrow(param_grid)), function(i) {
  row <- param_grid[i, ]
  out <- evaluate_params(
    kernel = row$kernel, cost = row$cost, gamma = row$gamma, degree = row$degree,
    data = model_df, n_folds = n_folds, n_repeats = n_repeats
  )
  cat(sprintf(
    "[%d/%d] kernel=%-10s cost=%-6s gamma=%-6s degree=%-3s -> BalAcc=%.3f F1=%.3f\n",
    i, nrow(param_grid), row$kernel, row$cost,
    ifelse(is.na(row$gamma), "-", row$gamma),
    ifelse(is.na(row$degree), "-", row$degree),
    out$Mean_Balanced_Accuracy, out$Mean_F1
  ))
  out
}))

# Rank primarily by Balanced Accuracy, tie-break with F1
results_grid <- results_grid %>%
  arrange(desc(Mean_Balanced_Accuracy), desc(Mean_F1))

cat("\nTop 10 hyperparameter combinations:\n")
print(head(results_grid, 10))

best_params <- results_grid[1, ]
cat("\nSelected best configuration:\n")
print(best_params)

# ---- Confusion matrix for the selected best model ----
# Built from ONE representative stratified 5-fold CV pass using the winning
# hyperparameters: every row gets exactly one out-of-fold prediction (made
# by a model that never saw that row during training), so this is a fair
# confusion matrix rather than the model grading its own training data.
# The Mean_* columns in results_grid (averaged over 20 repeats) remain the
# more robust performance estimate; this confusion matrix is for inspection
# and for direct comparison against dbscan_cm's table below.
set.seed(1001)
cv_folds_best <- createFolds(model_df$True_Label, k = n_folds, list = TRUE)
oof_pred <- factor(rep(NA_character_, nrow(model_df)), levels = c("Not unrest", "Unrest"))

for (f in seq_len(n_folds)) {
  test_idx   <- cv_folds_best[[f]]
  train_fold <- model_df[-test_idx, ]
  test_fold  <- model_df[test_idx, ]
  
  class_counts_fold <- table(train_fold$True_Label)
  cw <- sum(class_counts_fold) / (length(class_counts_fold) * class_counts_fold)
  
  svm_args <- list(
    formula = True_Label ~ ., data = train_fold,
    kernel = best_params$kernel, cost = best_params$cost,
    scale = TRUE, class.weights = cw
  )
  if (!is.na(best_params$gamma))  svm_args$gamma  <- best_params$gamma
  if (!is.na(best_params$degree)) svm_args$degree <- best_params$degree
  
  fit <- do.call(svm, svm_args)
  oof_pred[test_idx] <- predict(fit, test_fold)
}

svm_cm <- confusionMatrix(
  oof_pred,
  model_df$True_Label,
  positive = "Unrest",
  mode = "everything"
)

cat("\nConfusion matrix for best SVM model (out-of-fold, single 5-fold CV pass):\n")
print(svm_cm)

# ---- Refit the best model on the full dataset for deployment ----
full_class_counts <- table(model_df$True_Label)
full_cw <- sum(full_class_counts) / (length(full_class_counts) * full_class_counts)

final_svm_args <- list(
  formula = True_Label ~ .,
  data = model_df,
  kernel = best_params$kernel,
  cost = best_params$cost,
  scale = TRUE,
  class.weights = full_cw
)
if (!is.na(best_params$gamma))  final_svm_args$gamma  <- best_params$gamma
if (!is.na(best_params$degree)) final_svm_args$degree <- best_params$degree

final_svm_model <- do.call(svm, final_svm_args)
print(final_svm_model)

saveRDS(final_svm_model, file.path(OUT_DIR, "final_svm_model.rds"))
saveRDS(results_grid, file.path(OUT_DIR, "svm_grid_search_results.rds"))

# ---- Compare against DBSCAN using the same style of summary table ----
# Two SVM rows are reported: the 20-repeat average (most robust estimate)
# and the single out-of-fold pass shown in svm_cm above (directly comparable
# in form to dbscan_cm, since both are single confusion-matrix snapshots).
performance <- data.frame(
  Model = c("DBSCAN", "Best SVM (20-repeat CV average)", "Best SVM (single OOF pass)"),
  Balanced_Accuracy = c(
    dbscan_cm$byClass["Balanced Accuracy"],
    best_params$Mean_Balanced_Accuracy,
    svm_cm$byClass["Balanced Accuracy"]
  ),
  Sensitivity = c(
    dbscan_cm$byClass["Sensitivity"],
    best_params$Mean_Sensitivity,
    svm_cm$byClass["Sensitivity"]
  ),
  Specificity = c(
    dbscan_cm$byClass["Specificity"],
    best_params$Mean_Specificity,
    svm_cm$byClass["Specificity"]
  ),
  F1 = c(
    dbscan_cm$byClass["F1"],
    best_params$Mean_F1,
    svm_cm$byClass["F1"]
  )
)
print(performance)

write_xlsx(
  list(
    Grid_Search_Results = results_grid,
    Best_Parameters     = as.data.frame(best_params),
    Model_Comparison    = performance,
    DBSCAN_Confusion    = as.data.frame.matrix(dbscan_cm$table),
    SVM_Confusion        = as.data.frame.matrix(svm_cm$table)
  ),
  file.path(OUT_DIR, "SVM_kernel_tuning_and_results.xlsx")
)

cat(sprintf(
  "\nDone. Best model: kernel=%s, cost=%s, gamma=%s, degree=%s\n",
  best_params$kernel, best_params$cost,
  ifelse(is.na(best_params$gamma), "NA", best_params$gamma),
  ifelse(is.na(best_params$degree), "NA", best_params$degree)
))
cat("Saved: final_svm_model.rds, svm_grid_search_results.rds, SVM_kernel_tuning_and_results.xlsx\n")

