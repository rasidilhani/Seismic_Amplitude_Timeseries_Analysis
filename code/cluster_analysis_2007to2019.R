## ============================================================================
## Clustering of Weekly H-C (Entropy-Complexity) Features -- RSAM, 2007-2019
## Step 1: Correlation plot
## Step 2: K-means clustering (+ elbow plot to help pick k)
## Step 3: DBSCAN clustering (+ kNN distance plot to help pick eps)
## All axis/title text uses the serif font family throughout.
## ============================================================================

library(readr)
library(dplyr)
library(ggplot2)
library(ggcorrplot)
library(factoextra)   # kmeans/dbscan cluster viz (PCA-based), elbow/silhouette
library(cluster)
library(dbscan)       # dbscan() + kNNdistplot()
library(writexl)      # write_xlsx() -- export cluster results to Excel

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

#features <- c("H_Shannon", "C_Shannon",
#              "Fisher_Info", "C_Fisher",
#              "H_Renyi", "C_Renyi",
#              "H_Tsallis", "C_Tsallis",
#              "Var_H_Shannon", "Var_C_Shannon",
#              "Var_H_Tsallis", "Var_H_Renyi", 
#                "Disequilibrium"
#)

features <- c(
              "Fisher_Info", "C_Fisher",
              "Var_H_Shannon", "Var_C_Shannon",
              "Var_H_Tsallis", "Var_H_Renyi", 
              "Disequilibrium"
)

## Keep only complete cases on the chosen features (clustering can't
## handle NAs), and standardize -- required since these features are on
## very different scales (entropy/complexity vs. variance terms).
## `keep_idx` is kept so the cluster labels can be joined back onto the
## original weeks (Week_Index/Week_Start/Week_End) later for the Excel export.
keep_idx  <- complete.cases(df[, features])
df_feat   <- df[keep_idx, features]
df_scaled <- scale(df_feat)

## ============================================================================
## 1. CORRELATION PLOT
## ============================================================================

corr_mat <- cor(df_feat, use = "pairwise.complete.obs")

p_corr <- ggcorrplot(
  corr_mat,
  method = "square",
  type = "lower",
  lab = TRUE,
  lab_size = 3,
  colors = c("#3B4CC0", "white", "#B40426"),
  outline.color = "grey90"
) +
  labs(title = "Correlation of H-C Features (RSAM, 2007-2019)") +
  serif_theme

print(p_corr)
ggsave("correlation_plot_rsam.pdf", p_corr, width = 8, height = 7)

## ============================================================================
## 2. K-MEANS CLUSTERING
## ============================================================================

## Elbow plot (within-cluster sum of squares vs k) -- inspect this, pick
## the k where the curve bends, then set `k` below.
p_elbow <- fviz_nbclust(df_scaled, kmeans, method = "wss", k.max = 10) +
  labs(title = "Elbow Method for Choosing k") +
  serif_theme
print(p_elbow)

set.seed(0123456)
k <- 2   # <-- set this from the elbow plot above
km_res <- kmeans(df_scaled, centers = k, nstart = 25)

## fviz_cluster projects onto the first two principal components when
## there are >2 features, so this works directly with all 13 features.
p_kmeans <- fviz_cluster(
  km_res, data = df_scaled,
  geom = "point", ellipse.type = "convex",
  palette = "jco", ggtheme = theme_minimal()
) +
  labs(title = paste0("K-means Clustering (k = ", k, ")")) +
  serif_theme
print(p_kmeans)
ggsave("kmeans_cluster_plot_rsam.pdf", p_kmeans, width = 8, height = 7)

## ============================================================================
## 3. DBSCAN CLUSTERING
## ============================================================================

## kNN distance plot to help pick eps: look for the "knee" (sharp bend)
## in the sorted k-distances. minPts is commonly set to ~2 * ncol(data).
minPts_val <- 2 * ncol(df_scaled)

par(family = "serif")
kNNdistplot(df_scaled, k = minPts_val)
title(main = "k-NN Distance Plot (for choosing eps)", family = "serif")
abline(h = 1.0, lty = 2, col = "red")  # move this line to match the knee you see

eps_val <- 1   # <-- set this by reading off the knee in the plot above

db_res <- dbscan(df_scaled, eps = eps_val, minPts = minPts_val)

p_dbscan <- fviz_cluster(
  list(data = df_scaled, cluster = db_res$cluster),
  geom = "point", ellipse = FALSE,
  palette = "jco", ggtheme = theme_minimal()
) +
  labs(title = paste0("DBSCAN Clustering (eps = ", eps_val,
                      ", minPts = ", minPts_val, ")")) +
  serif_theme
print(p_dbscan)
ggsave("dbscan_cluster_plot_rsam.pdf", p_dbscan, width = 8, height = 7)

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
  path = "cluster_results_rsam.xlsx"
)

###########################################################################
