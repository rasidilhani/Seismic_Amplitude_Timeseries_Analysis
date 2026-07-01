# ==============================================================================
# WIZ Ordinal Pattern Analysis: PCA + Clustering + Binary Classification
# Purpose: identify which months show eruption-driven amplitude peaks, using
# pooled (all-years) features rather than per-year clustering.
# ==============================================================================

library(here)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(factoextra)   # PCA / cluster visualisation
library(cluster)      # silhouette
library(randomForest) # classification

# ------------------------------------------------------------------------
# 1) USER SWITCHES -- change these two lines and re-run the whole script
# ------------------------------------------------------------------------
VARIABLE    <- "RSAM"        # options: "Displacement" or "RSAM"
PERIOD_TYPE <- "Monthly"     # options: "Monthly" or "Quarterly"

# Monthly results live in sheet 1, Quarterly in sheet 2 -- picked automatically
sheet_to_read <- if (PERIOD_TYPE == "Monthly") 1 else 2

# The 13 ordinal-pattern features we analyse (same names in both sheets)
features <- c("H_Shannon","H_Tsallis","H_Renyi","C_Shannon","H_Fisher","C_Fisher",
              "C_Tsallis","C_Renyi","Disequilibrium",
              "Var_H_Shannon","Var_H_Tsallis","Var_H_Renyi","Var_C_Shannon")

# ------------------------------------------------------------------------
# 2) LOAD + FILTER DATA
# ------------------------------------------------------------------------
input_path <- here("results", "WIZ_OrdinalPatterns_D5_Monthly_Quarterly.xlsx")
out_dir    <- here("results", "WIZ_Clustering_Classification", VARIABLE, PERIOD_TYPE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_xlsx(input_path, sheet = sheet_to_read) %>%
  filter(Variable == VARIABLE, Period_Type == PERIOD_TYPE) %>%
  mutate(
    Period = factor(Period, levels = if (PERIOD_TYPE == "Monthly")
      month.abb else c("Q1","Q2","Q3","Q4")),
    Year = as.numeric(Year)
  ) %>%
  arrange(Year, Period)

# Pool ALL years together (this is the key fix vs. per-year clustering --
# 12 points/year is too few; pooling gives ~130+ rows so clusters reflect
# a real long-term baseline, and eruption months stand out against it)
cat("Rows pooled across all years:", nrow(df), "\n")

# ------------------------------------------------------------------------
# 3) SCALE FEATURES (PCA/clustering need comparable units)
# ------------------------------------------------------------------------
X <- df %>% select(all_of(features)) %>% drop_na()
keep_rows <- complete.cases(df %>% select(all_of(features)))  # track which rows survived NA removal
df_clean  <- df[keep_rows, ]
X_scaled  <- scale(X)

# ------------------------------------------------------------------------
# 4) PCA -- reduces 13 correlated features down to a few components, and
#    tells us (via loadings) which original features matter most
# ------------------------------------------------------------------------
pca <- prcomp(X_scaled, center = FALSE, scale. = FALSE)  # already scaled above

print(summary(pca))  # check % variance explained by PC1, PC2, ...

# Save scree plot (how many components are worth keeping)
png(file.path(out_dir, "pca_scree.pdf"), width = 800, height = 600)
print(fviz_eig(pca, addlabels = TRUE))
dev.off()

# Save loadings (which features drive PC1/PC2 -- use this to justify a
# reduced feature set, e.g. keep only features with high |loading|)
loadings <- as.data.frame(pca$rotation[, 1:2])
loadings$Feature <- rownames(loadings)
write.csv(loadings, file.path(out_dir, "pca_loadings.csv"), row.names = FALSE)
print(loadings[order(-abs(loadings$PC1)), ])  # features ranked by PC1 importance

# Attach PC scores back to the data for plotting / clustering
df_clean$PC1 <- pca$x[, 1]
df_clean$PC2 <- pca$x[, 2]

# ------------------------------------------------------------------------
# 5) CLUSTERING on the pooled PCA scores
#    Goal: see if months naturally split into seismic "regimes"
#    (e.g. quiet / elevated / eruptive), not just hunt for outliers.
# ------------------------------------------------------------------------
# Choose k via silhouette width (tries k = 2..5, picks best average score)
sil_scores <- sapply(2:5, function(k) {
  km <- kmeans(df_clean[, c("PC1","PC2")], centers = k, nstart = 25)
  mean(silhouette(km$cluster, dist(df_clean[, c("PC1","PC2")]))[, 3])
})
best_k <- which.max(sil_scores) + 1
cat("Best k by silhouette score:", best_k, "\n")

km <- kmeans(df_clean[, c("PC1","PC2")], centers = best_k, nstart = 25)
df_clean$Cluster <- factor(km$cluster)

# Plot: points coloured by cluster, labelled with Year-Period so you can
# visually confirm whether known eruption months (e.g. 2019) fall into
# their own cluster
p <- ggplot(df_clean, aes(PC1, PC2, color = Cluster,
                          label = paste0(Year, "-", Period))) +
  geom_point(size = 3) +
  geom_text(vjust = -0.6, size = 2.5, check_overlap = TRUE) +
  theme_minimal() +
  labs(title = paste(VARIABLE, PERIOD_TYPE, "- Pooled PCA + k-means clusters"))
ggsave(file.path(out_dir, "cluster_plot.pdf"), p, width = 9, height = 6)

write.csv(df_clean, file.path(out_dir, "clustered_results.csv"), row.names = FALSE)

# ------------------------------------------------------------------------
# 5b) HC-PLANE PLOTS (Shannon and Fisher) -- the classic entropy vs.
#     complexity view used in ordinal-pattern literature. Points sitting
#     away from the main "cloud" (esp. high complexity / low-to-mid
#     entropy) are the ones worth checking against eruption timing.
#     Coloured by cluster and labelled by Year-Period for easy spotting.
# ------------------------------------------------------------------------
plot_hc_plane <- function(data, h_col, c_col, title, fname) {
  p <- ggplot(data, aes(x = .data[[h_col]], y = .data[[c_col]],
                        color = Cluster, label = paste0(Year, "-", Period))) +
    geom_point(size = 3) +
    geom_text(vjust = -0.6, size = 2.5, check_overlap = TRUE) +
    theme_minimal() +
    labs(title = title, x = h_col, y = c_col)
  ggsave(file.path(out_dir, fname), p, width = 9, height = 6)
  p
}

# Shannon HC plane
plot_hc_plane(df_clean, "H_Shannon", "C_Shannon",
              paste(VARIABLE, PERIOD_TYPE, "- Shannon HC plane"),
              "hc_plane_shannon.pdf")

# Fisher HC plane (note: scale/shape differs from Shannon -- Fisher
# information is more sensitive to abrupt local changes, often useful
# for picking up sharp eruption-related transitions)
plot_hc_plane(df_clean, "H_Fisher", "C_Fisher",
              paste(VARIABLE, PERIOD_TYPE, "- Fisher HC plane"),
              "hc_plane_fisher.pdf")

# ------------------------------------------------------------------------
# 6) BINARY CLASSIFICATION (optional -- needs your visually-identified
#    eruption months first)
# ------------------------------------------------------------------------
# Fill in the eruption months you identified by eye, e.g.:
# eruption_months <- tibble(Year = c(2019, 2019), Period = c("Apr","May"))
# Leave empty and skip this section until you have your list.

eruption_months <- tibble(
  Year   = c(2019),     # <- add more years here
  Period = c("Apr")     # <- add matching months/quarters here
)

if (nrow(eruption_months) > 0) {
  
  df_clean <- df_clean %>%
    mutate(Eruption = if_else(
      paste(Year, Period) %in% paste(eruption_months$Year, eruption_months$Period),
      "Yes", "No"
    ) %>% factor())
  
  # Random forest: simple, handles correlated features, gives importance
  # ranking for free (answers "which features matter most")
  rf <- randomForest(
    x = df_clean[, features],
    y = df_clean$Eruption,
    importance = TRUE,
    ntree = 500
  )
  
  print(rf)                          # confusion matrix / OOB error
  print(importance(rf))              # feature importance scores
  
  pdf(file.path(out_dir, "rf_importance.pdf"), width = 800, height = 600)
  varImpPlot(rf, main = paste(VARIABLE, PERIOD_TYPE, "- Feature importance"))
  dev.off()
  
  saveRDS(rf, file.path(out_dir, "rf_model.rds"))
  
} else {
  cat("No eruption months supplied yet -- skipping classification step.\n",
      "Fill in 'eruption_months' above once you've identified them visually.\n")
}

cat("\nDone. Outputs saved to:", out_dir, "\n")