# ============================================================
# RSAM yearly monthly feature analysis
# Each year contains 12 months
# Outputs saved separately
# ============================================================

library(readxl)
library(tidyverse)
library(here)
library(cluster)
library(factoextra)
library(caret)
library(randomForest)
library(corrplot)
library(patchwork)

# ============================================================
# Chunk 1: Load data
# ============================================================

D <- 5

input_path <- here("results", "WIZ_OrdinalPatterns_D5_Monthly_Quarterly.xlsx")

out_dir <- here("results", "Stepwise_RSAM_Yearly_Monthly_Clustering")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df_all <- read_xlsx(input_path, sheet = 1)

df_all <- df_all %>%
  filter(
    Period_Type == "Monthly",
    Year >= 2011,
    Year <= 2021
  ) %>%
  mutate(
    Period = factor(Period, levels = month.abb),
    Year = factor(Year)
  )

unique(df_all$Variable)

# ============================================================
# Chunk 2: Select RSAM and reduced features
# ============================================================

selected_variable <- "RSAM"

feature_cols <- c(
  "H_Shannon",
  "C_Shannon",
  "H_Fisher",
  "C_Fisher",
  "Disequilibrium",
  "Var_H_Shannon",
  "Var_C_Shannon"
)

missing_features <- setdiff(feature_cols, names(df_all))

if (length(missing_features) > 0) {
  stop("Missing columns: ", paste(missing_features, collapse = ", "))
}

df_rsam <- df_all %>%
  filter(Variable == selected_variable) %>%
  select(
    Year,
    Period,
    Variable,
    all_of(feature_cols)
  ) %>%
  drop_na(all_of(feature_cols))

summary_stats <- df_rsam %>%
  summarise(across(
    all_of(feature_cols),
    list(
      Min = ~ min(.x, na.rm = TRUE),
      Max = ~ max(.x, na.rm = TRUE),
      Mean = ~ mean(.x, na.rm = TRUE),
      SD = ~ sd(.x, na.rm = TRUE)
    )
  ))

write_csv(
  summary_stats,
  file.path(out_dir, "01_RSAM_Feature_Summary_Statistics.csv")
)

print(summary_stats)

# ============================================================
# Chunk 3: Correlation heatmaps by year
# Output 1: one faceted-style saved set, one plot per year
# ============================================================

cor_dir <- file.path(out_dir, "02_Correlation_By_Year")
dir.create(cor_dir, recursive = TRUE, showWarnings = FALSE)

for (yr in levels(df_rsam$Year)) {
  
  year_df <- df_rsam %>%
    filter(Year == yr)
  
  cor_mat <- year_df %>%
    select(all_of(feature_cols)) %>%
    cor(use = "pairwise.complete.obs", method = "spearman")
  
  write_csv(
    as.data.frame(cor_mat) %>% rownames_to_column("Feature"),
    file.path(cor_dir, paste0("RSAM_Correlation_Matrix_", yr, ".csv"))
  )
  
  corrplot(
    cor_mat,
    method = "color",
    type = "upper",
    order = "hclust",
    tl.col = "black",
    tl.cex = 0.8,
    addCoef.col = "black",
    number.cex = 0.6,
    title = paste("RSAM feature correlation:", yr),
    mar = c(0, 0, 2, 0)
  )
  
  pdf(file.path(cor_dir, paste0("RSAM_Correlation_Heatmap_", yr, ".pdf")), width = 8, height = 7)
  corrplot(
    cor_mat,
    method = "color",
    type = "upper",
    order = "hclust",
    tl.col = "black",
    tl.cex = 0.8,
    addCoef.col = "black",
    number.cex = 0.6,
    title = paste("RSAM feature correlation:", yr),
    mar = c(0, 0, 2, 0)
  )
  dev.off()
}

# ============================================================
# Chunk 4: Log-scale feature boxplots
# Output 2: faceted by feature, months on x-axis, all years together
# ============================================================

eps <- 1e-10

box_df <- df_rsam %>%
  pivot_longer(
    cols = all_of(feature_cols),
    names_to = "Feature",
    values_to = "Value"
  ) %>%
  mutate(
    Log_Value = log10(Value + eps)
  )

p_box <- ggplot(box_df, aes(x = Period, y = Log_Value, fill = Period)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_grid(Feature ~ Year, scales = "free_y") +
  labs(
    title = "RSAM ordinal-pattern features by month within each year",
    x = "Month",
    y = expression(log[10](value + epsilon))
  ) +
  theme_bw(base_size = 9, base_family = "serif") +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 7),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_box)

ggsave(
  file.path(out_dir, "03_RSAM_Log_Boxplots_Months_By_Year.pdf"),
  p_box,
  width = 18,
  height = 12,
  units = "in"
)

# ============================================================
# Chunk 5: Helper function
# ============================================================

scale_features <- function(df) {
  df %>%
    select(all_of(feature_cols)) %>%
    scale() %>%
    as.data.frame()
}

# ============================================================
# Chunk 6: K-means silhouette evaluation by year
# Output 3: one faceted silhouette plot
# ============================================================

silhouette_results <- list()

for (yr in levels(df_rsam$Year)) {
  
  year_df <- df_rsam %>%
    filter(Year == yr)
  
  x_scaled <- scale_features(year_df)
  
  max_k <- min(5, nrow(year_df) - 1)
  
  if (nrow(year_df) < 4 || max_k < 2) next
  
  sil_df <- map_dfr(2:max_k, function(k) {
    
    set.seed(123)
    
    km <- kmeans(
      x_scaled,
      centers = k,
      nstart = 50
    )
    
    sil <- silhouette(km$cluster, dist(x_scaled))
    
    tibble(
      Year = yr,
      k = k,
      Average_Silhouette = mean(sil[, 3])
    )
  })
  
  silhouette_results[[yr]] <- sil_df
}

silhouette_all <- bind_rows(silhouette_results)

write_csv(
  silhouette_all,
  file.path(out_dir, "04_RSAM_KMeans_Silhouette_By_Year.csv")
)

p_sil <- ggplot(silhouette_all, aes(x = k, y = Average_Silhouette)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ Year, ncol = 4) +
  labs(
    title = "K-means silhouette evaluation by year",
    x = "Number of clusters, k",
    y = "Average silhouette width"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black")
  )

print(p_sil)

ggsave(
  file.path(out_dir, "04_RSAM_KMeans_Silhouette_By_Year.pdf"),
  p_sil,
  width = 10,
  height = 7,
  units = "in"
)

# ============================================================
# Chunk 7: K-means clustering and fviz plots by year
# Output 4: fviz cluster plots
# Output 5: cluster assignment table
# ============================================================

fviz_dir <- file.path(out_dir, "05_fviz_KMeans_By_Year")
dir.create(fviz_dir, recursive = TRUE, showWarnings = FALSE)

cluster_assignments <- list()
fviz_plots <- list()

for (yr in levels(df_rsam$Year)) {
  
  year_df <- df_rsam %>%
    filter(Year == yr)
  
  sil_y <- silhouette_all %>%
    filter(Year == yr)
  
  if (nrow(year_df) < 4 || nrow(sil_y) == 0) next
  
  best_k <- sil_y %>%
    arrange(desc(Average_Silhouette)) %>%
    slice(1) %>%
    pull(k)
  
  x_scaled <- scale_features(year_df)
  
  set.seed(123)
  
  km <- kmeans(
    x_scaled,
    centers = best_k,
    nstart = 50
  )
  
  year_cluster_df <- year_df %>%
    mutate(
      Best_k = best_k,
      Cluster = factor(km$cluster)
    )
  
  cluster_assignments[[yr]] <- year_cluster_df
  
  p_fviz <- fviz_cluster(
    list(
      data = x_scaled,
      cluster = km$cluster
    ),
    geom = "point",
    ellipse.type = "convex",
    palette = "jco",
    ggtheme = theme_bw(base_size = 11, base_family = "serif"),
    main = paste("RSAM K-means clusters:", yr, ", k =", best_k)
  ) +
    labs(color = "Cluster") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black")
    )
  
  print(p_fviz)
  
  ggsave(
    file.path(fviz_dir, paste0("RSAM_fviz_KMeans_", yr, ".pdf")),
    p_fviz,
    width = 7,
    height = 5,
    units = "in"
  )
  
  fviz_plots[[yr]] <- p_fviz
}

cluster_all <- bind_rows(cluster_assignments)

write_csv(
  cluster_all,
  file.path(out_dir, "05_RSAM_KMeans_Cluster_Assignments_By_Year.csv")
)

cluster_all %>%
  select(Year, Period, Best_k, Cluster)

# ============================================================
# Chunk 8: Cluster membership summary
# Output 6: month-to-cluster table
# ============================================================

cluster_summary <- cluster_all %>%
  count(Year, Cluster, Period) %>%
  arrange(Year, Cluster, Period)

write_csv(
  cluster_summary,
  file.path(out_dir, "06_RSAM_Cluster_Membership_Summary_By_Year.csv")
)

print(cluster_summary)

p_cluster_months <- cluster_all %>%
  ggplot(aes(x = Period, y = Year, fill = Cluster)) +
  geom_tile(color = "white") +
  labs(
    title = "RSAM monthly cluster membership by year",
    x = "Month",
    y = "Year",
    fill = "Cluster"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

print(p_cluster_months)

ggsave(
  file.path(out_dir, "06_RSAM_Cluster_Membership_Heatmap_By_Year.pdf"),
  p_cluster_months,
  width = 9,
  height = 5,
  units = "in"
)

# ============================================================
# Chunk 9: PCA plots by year
# Output 7: PCA plots coloured by cluster, labelled by month
# ============================================================

pca_dir <- file.path(out_dir, "07_PCA_KMeans_By_Year")
dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)

pca_all_scores <- list()

for (yr in levels(cluster_all$Year)) {
  
  year_df <- cluster_all %>%
    filter(Year == yr) %>%
    mutate(
      Cluster = factor(Cluster),
      Period = factor(Period, levels = month.abb)
    )
  
  if (nrow(year_df) < 4 || n_distinct(year_df$Cluster) < 2) next
  
  x_scaled <- year_df %>%
    select(all_of(feature_cols)) %>%
    scale()
  
  pca_fit <- prcomp(x_scaled, center = FALSE, scale. = FALSE)
  
  variance_explained <- summary(pca_fit)$importance[2, 1:2] * 100
  
  pca_scores <- as.data.frame(pca_fit$x[, 1:2]) %>%
    mutate(
      Year = year_df$Year,
      Period = year_df$Period,
      Cluster = year_df$Cluster,
      PC1_Percent = variance_explained[1],
      PC2_Percent = variance_explained[2]
    )
  
  pca_all_scores[[yr]] <- pca_scores
  
  p_pca <- ggplot(
    pca_scores,
    aes(x = PC1, y = PC2, color = Cluster, label = Period)
  ) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
    stat_ellipse(
      aes(group = Cluster),
      linewidth = 0.5,
      linetype = "dashed",
      show.legend = FALSE
    ) +
    labs(
      title = paste("PCA of RSAM ordinal-pattern features:", yr),
      subtitle = "Points labelled by month and coloured by K-means cluster",
      x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
      color = "Cluster"
    ) +
    theme_bw(base_size = 11, base_family = "serif") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text = element_text(color = "black"),
      legend.position = "bottom"
    )
  
  print(p_pca)
  
  ggsave(
    file.path(pca_dir, paste0("RSAM_PCA_KMeans_Clusters_", yr, ".pdf")),
    p_pca,
    width = 7,
    height = 5,
    units = "in"
  )
}

pca_all_scores <- bind_rows(pca_all_scores)

write_csv(
  pca_all_scores,
  file.path(out_dir, "07_RSAM_PCA_Scores_By_Year.csv")
)

# ============================================================
# Chunk 10: Random Forest feature importance by year
# Output 8: RF accuracy table
# Output 9: RF feature importance plots
# ============================================================

rf_dir <- file.path(out_dir, "08_Random_Forest_By_Year")
dir.create(rf_dir, recursive = TRUE, showWarnings = FALSE)

rf_results <- list()
rf_importance_all <- list()

for (yr in levels(cluster_all$Year)) {
  
  year_df <- cluster_all %>%
    filter(Year == yr) %>%
    mutate(Cluster = factor(Cluster))
  
  if (nrow(year_df) < 8 || n_distinct(year_df$Cluster) < 2) next
  
  x <- year_df %>%
    select(all_of(feature_cols))
  
  y <- year_df$Cluster
  
  ctrl <- trainControl(
    method = "LOOCV",
    classProbs = FALSE,
    savePredictions = "final"
  )
  
  set.seed(123)
  
  rf_fit <- train(
    x = x,
    y = y,
    method = "rf",
    trControl = ctrl,
    ntree = 500,
    importance = TRUE
  )
  
  acc_df <- tibble(
    Year = yr,
    Accuracy = max(rf_fit$results$Accuracy, na.rm = TRUE),
    Kappa = max(rf_fit$results$Kappa, na.rm = TRUE),
    Best_mtry = rf_fit$bestTune$mtry
  )
  
  rf_results[[yr]] <- acc_df
  
  imp_df <- varImp(rf_fit)$importance %>%
    as.data.frame() %>%
    rownames_to_column("Feature")
  
  if (!"Overall" %in% names(imp_df)) {
    imp_df <- imp_df %>%
      mutate(
        Overall = rowMeans(
          select(., where(is.numeric)),
          na.rm = TRUE
        )
      )
  }
  
  imp_df <- imp_df %>%
    arrange(desc(Overall)) %>%
    mutate(Year = yr)
  
  rf_importance_all[[yr]] <- imp_df
  
  p_imp <- ggplot(imp_df, aes(x = reorder(Feature, Overall), y = Overall)) +
    geom_col(fill = "grey35") +
    coord_flip() +
    labs(
      title = paste("Random Forest feature importance:", yr),
      x = NULL,
      y = "Importance"
    ) +
    theme_bw(base_size = 11, base_family = "serif") +
    theme(
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  print(p_imp)
  
  ggsave(
    file.path(rf_dir, paste0("RSAM_RF_Feature_Importance_", yr, ".pdf")),
    p_imp,
    width = 7,
    height = 5,
    units = "in"
  )
}

rf_accuracy_all <- bind_rows(rf_results)
rf_importance_all <- bind_rows(rf_importance_all)

write_csv(
  rf_accuracy_all,
  file.path(out_dir, "08_RSAM_RF_Accuracy_By_Year.csv")
)

write_csv(
  rf_importance_all,
  file.path(out_dir, "09_RSAM_RF_Feature_Importance_By_Year.csv")
)

print(rf_accuracy_all)