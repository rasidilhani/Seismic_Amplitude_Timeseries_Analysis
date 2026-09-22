# =============================================================================
# WIZ and WSRZ seismic activity classification - final, reduced method set
#
# Based on the full comparison done earlier, this keeps only the methods
# that are actually worth reporting:
#
#   - K-means and DBSCAN: unsupervised clustering.
#   - SVM: the only supervised classifier that caught any real eruption
#     weeks.
#   - Feature sets compared: 13 original features, 7 selected
#     features, and a 4-feature variance-only set.
#   - Class-imbalance methods compared: no balancing, oversampling,
#     undersampling, SMOTE.
#
# =============================================================================

library(tidyverse)
library(corrplot)
library(dbscan)
library(e1071)         # SVM
library(caret)
library(cluster)       # for silhouette

if (!requireNamespace("mclust", quietly = TRUE)) {
  stop("Package 'mclust' is needed for the Adjusted Rand Index. Install it with: install.packages('mclust')")
}
# mclust is used only via mclust::adjustedRandIndex() below, not library()'d in full,
# because it overwrites some tidyverse/purrr function names (e.g. map)

set.seed(123456789, kind = "Mersenne-Twister")

# --- plot style: serif font, article-style look, applied to every plot -----

theme_article <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = "serif") +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      axis.title   = element_text(face = "bold", size = base_size),
      axis.text    = element_text(colour = "black", size = base_size - 1),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text  = element_text(size = base_size - 1),
      strip.text   = element_text(face = "bold", size = base_size - 1),
      panel.grid   = element_blank()
    )
}

theme_set(theme_article())

# --- paths ---------------------------------------------------------------

DATA_DIR    <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"
RESULTS_DIR <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"
CLUSTER_DIR <- file.path(RESULTS_DIR, "cluster_analysis_R3")
dir.create(CLUSTER_DIR, recursive = TRUE, showWarnings = FALSE)

VAL_FILE <- file.path(DATA_DIR, "VAL_WhakaariWhiteIsland.csv")

station_files <- c(
  WIZ  = file.path(DATA_DIR, "HC_RSAM_WIZ_20140701_to_20191231.csv"),
  WSRZ = file.path(DATA_DIR, "HC_RSAM_WSRZ_20140701_to_20191231.csv")
)

ANALYSIS_START <- as.POSIXct("2014-07-01 00:00:00", tz = "UTC")
ANALYSIS_END   <- as.POSIXct("2019-12-31 23:59:59", tz = "UTC")
VAL_SYSTEM_VERSION <- 3

# --- feature sets ----------------------------------------------------------

all_features <- c(
  "H_Shannon", "C_Shannon", "Fisher_Info", "C_Fisher",
  "H_Renyi", "C_Renyi", "H_Tsallis", "C_Tsallis",
  "Var_H_Shannon", "Var_C_Shannon", "Var_H_Renyi", "Var_H_Tsallis",
  "Disequilibrium"
)

feature_cols <- c(
  "Fisher_Info", "C_Fisher",
  "Var_H_Shannon", "Var_C_Shannon",
  "Var_H_Tsallis", "Var_H_Renyi",
  "Disequilibrium"
)

variance_features <- c("Var_H_Shannon", "Var_C_Shannon", "Var_H_Renyi", "Var_H_Tsallis")

# feature sets compared for SVM classification accuracy
feature_sets <- list(
  "13_Features"       = all_features,
  "7_Features"        = feature_cols,
  "4_Features" = variance_features
)

# feature sets compared for clustering quality (kept to two, since DBSCAN's
# grid search is slow to re-run)
clustering_feature_sets <- list(
  "7_Features"        = feature_cols,
  "4_Features" = variance_features
)

balance_methods <- c("No_Balancing", "oversample", "undersample", "smote")

model_label_levels <- c(
  "1_Background_unrest",
  "2_Moderate_to_heightened_unrest",
  "3_Minor_to_moderate_eruption"
)

# =============================================================================
# VAL activity data: read it, keep version 3, collapse to 3 classes
#
# Uses the UTC period columns (not NZ local-time), to match the OP feature
# data, which was computed in UTC.
#
# read_csv() already detects "PeriodStart UTC" and "PeriodEnd UTC" as proper
# datetime (POSIXct) columns, so we just rename them - no need to strip a
# trailing "UTC" and re-parse as text. Re-parsing was the source of the
# midnight-timestamp bug described at the top of this script.
# =============================================================================

val_df <- read_csv(VAL_FILE, show_col_types = FALSE) %>%
  rename(
    PeriodStart_UTC = `PeriodStart UTC`,
    PeriodEnd_UTC   = `PeriodEnd UTC`
  )

open_ended <- is.na(val_df$PeriodEnd_UTC) | val_df$PeriodEnd_UTC >= as.POSIXct("9999-01-01", tz = "UTC")
val_df$PeriodEnd_UTC[open_ended] <- ANALYSIS_END

val_df <- val_df %>%
  filter(
    `VAL System Version` == VAL_SYSTEM_VERSION,
    PeriodEnd_UTC > ANALYSIS_START,
    PeriodStart_UTC < ANALYSIS_END
  ) %>%
  mutate(
    Model_Label_Num = case_when(
      Description == "Minor volcanic unrest" ~ 1L,
      Description == "Moderate to heightened volcanic unrest" ~ 2L,
      Description %in% c("Minor volcanic eruption", "Moderate volcanic eruption") ~ 3L,
      TRUE ~ NA_integer_
    ),
    Model_Label = model_label_levels[Model_Label_Num]
  ) %>%
  filter(!is.na(Model_Label_Num))

assign_val_class <- function(week_start, week_end) {
  overlap <- val_df %>% filter(PeriodStart_UTC < week_end, PeriodEnd_UTC > week_start)
  if (nrow(overlap) == 0) return(NA_integer_)
  max(overlap$Model_Label_Num, na.rm = TRUE)
}

# =============================================================================
# Helper: 13-feature correlogram (displayed and saved as its own PDF)
# =============================================================================

make_correlogram <- function(data, station_name) {
  
  missing_features <- setdiff(all_features, names(data))
  if (length(missing_features) > 0) {
    stop(station_name, " is missing: ", paste(missing_features, collapse = ", "))
  }
  
  cor_matrix <- cor(data[, all_features], use = "pairwise.complete.obs", method = "spearman")
  
  draw_plot <- function() {
    par(family = "serif")
    corrplot(cor_matrix, method = "color", type = "upper", order = "hclust",
             tl.col = "black", tl.srt = 45, tl.cex = 0.75,
             addCoef.col = "black", number.cex = 0.45, diag = FALSE,
             title = paste(station_name, "- Correlation of 13 Original Features"),
             mar = c(0, 0, 2, 0))
  }
  
  draw_plot()   # on-screen
  
  pdf(file.path(CLUSTER_DIR, paste0(station_name, "_13_Feature_Correlogram.pdf")), width = 11, height = 9)
  draw_plot()
  dev.off()
  
  cor_matrix
}

# =============================================================================
# Helper: build one PCA-based cluster plot, for K-means or DBSCAN
# =============================================================================
#
# K-means and DBSCAN run on the same scaled feature matrix X, so they share
# one PCA projection. Axis labels show the variance each PC captures.
# =============================================================================

make_cluster_plot <- function(pca, cluster_labels, station_name, method_name) {
  
  var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)
  pc_scores <- as_tibble(pca$x[, 1:2]) %>% mutate(Cluster = factor(cluster_labels))
  
  ggplot(pc_scores, aes(PC1, PC2, colour = Cluster)) +
    geom_point(alpha = 0.7, size = 2) +
    labs(
      title = paste0(station_name, ": ", method_name, " Clustering"),
      x = paste0("PC1 (", var_explained[1], "% of variance)"),
      y = paste0("PC2 (", var_explained[2], "% of variance)")
    )
}

# =============================================================================
# Helper: DBSCAN, searched to find exactly 3 clusters with least noise
# =============================================================================

find_three_cluster_dbscan <- function(X) {
  
  grid <- expand_grid(eps = seq(0.05, 5, by = 0.05), minPts = 3:20) %>%
    mutate(n_clusters = NA_integer_, noise = NA_integer_)
  
  for (i in seq_len(nrow(grid))) {
    fit <- dbscan(X, eps = grid$eps[i], minPts = grid$minPts[i])
    grid$n_clusters[i] <- length(setdiff(unique(fit$cluster), 0))
    grid$noise[i]      <- sum(fit$cluster == 0)
  }
  
  candidates <- grid %>% filter(n_clusters == 3)
  if (nrow(candidates) == 0) stop("No DBSCAN solution gives exactly 3 clusters.")
  
  best <- candidates %>% arrange(noise, abs(eps - 1.1), minPts) %>% slice(1)
  final_fit <- dbscan(X, eps = best$eps, minPts = best$minPts)
  
  list(model = final_fit, eps = best$eps, minPts = best$minPts)
}

# =============================================================================
# Class-imbalance methods
# =============================================================================
#
# All three work on one training fold at a time and return a more balanced
# version of it. The test fold is never touched by these functions.
# =============================================================================

pick_one <- function(x) {
  if (length(x) == 1) return(x)
  sample(x, 1)
}

oversample_data <- function(X, y) {
  target_n <- max(table(y))
  idx <- unlist(lapply(levels(y), function(cls) {
    cls_idx <- which(y == cls)
    sample(cls_idx, target_n, replace = TRUE)
  }))
  list(X = X[idx, , drop = FALSE], y = y[idx])
}

undersample_data <- function(X, y) {
  target_n <- min(table(y))
  idx <- unlist(lapply(levels(y), function(cls) {
    cls_idx <- which(y == cls)
    sample(cls_idx, target_n, replace = FALSE)
  }))
  list(X = X[idx, , drop = FALSE], y = y[idx])
}

smote_data <- function(X, y, k = 5) {
  
  target_n <- max(table(y))
  X_new <- X
  y_new <- as.character(y)
  
  for (cls in levels(y)) {
    
    cls_idx <- which(y == cls)
    n_cls <- length(cls_idx)
    n_needed <- target_n - n_cls
    if (n_needed <= 0) next
    
    cls_X <- X[cls_idx, , drop = FALSE]
    dist_matrix <- as.matrix(dist(cls_X))
    
    synthetic <- matrix(NA_real_, nrow = n_needed, ncol = ncol(X))
    colnames(synthetic) <- colnames(X)
    
    for (i in seq_len(n_needed)) {
      
      point_idx <- pick_one(seq_len(n_cls))
      point <- cls_X[point_idx, ]
      
      neighbour_order <- order(dist_matrix[point_idx, ])
      neighbour_order <- neighbour_order[neighbour_order != point_idx]
      neighbours <- head(neighbour_order, min(k, length(neighbour_order)))
      
      neighbour_idx <- pick_one(neighbours)
      neighbour <- cls_X[neighbour_idx, ]
      
      gap <- runif(1)
      synthetic[i, ] <- point + gap * (neighbour - point)
    }
    
    X_new <- rbind(X_new, synthetic)
    y_new <- c(y_new, rep(cls, n_needed))
  }
  
  list(X = X_new, y = factor(y_new, levels = levels(y)))
}

balance_training_data <- function(X, y, method) {
  switch(method,
         No_Balancing = list(X = X, y = y),
         oversample   = oversample_data(X, y),
         undersample  = undersample_data(X, y),
         smote        = smote_data(X, y),
         stop("Unknown balancing method: ", method)
  )
}

# =============================================================================
# Helper: per-class metrics, with raw counts (not just percentages) - and
# balanced accuracy
# =============================================================================
#
# Plain accuracy is misleading when one class has far more rows than the
# others. Balanced accuracy averages the RECALL of each class, so all three
# count equally. Correct/Total are the raw counts behind each recall value -
# with only 3-4 true eruption weeks, "caught 2 of 4" is more honest than "50%".
# =============================================================================

class_metrics <- function(true, predicted) {
  
  classes <- levels(true)
  per_class <- tibble(Class = classes, Correct = NA_integer_, Total = NA_integer_,
                      Recall = NA_real_, Precision = NA_real_, F1 = NA_real_)
  
  for (i in seq_along(classes)) {
    cls <- classes[i]
    
    true_positive   <- sum(true == cls & predicted == cls)
    actual_count    <- sum(true == cls)
    predicted_count <- sum(predicted == cls)
    
    recall    <- if (actual_count > 0) true_positive / actual_count else NA_real_
    precision <- if (predicted_count > 0) true_positive / predicted_count else NA_real_
    f1 <- if (!is.na(recall) && !is.na(precision) && (recall + precision) > 0) {
      2 * recall * precision / (recall + precision)
    } else {
      NA_real_
    }
    
    per_class$Correct[i]   <- true_positive
    per_class$Total[i]     <- actual_count
    per_class$Recall[i]    <- recall
    per_class$Precision[i] <- precision
    per_class$F1[i]        <- f1
  }
  
  list(per_class = per_class, balanced_accuracy = mean(per_class$Recall, na.rm = TRUE))
}

# turns the 3-row per-class table into one row of columns, e.g.
# Eruption_Correct, Eruption_Total, Eruption_Recall, Eruption_Precision ...
# so it can sit alongside Accuracy/Balanced_Accuracy in one wide table
flatten_per_class <- function(per_class) {
  
  short_names <- c(
    "1_Background_unrest"             = "Background",
    "2_Moderate_to_heightened_unrest" = "Moderate",
    "3_Minor_to_moderate_eruption"    = "Eruption"
  )
  per_class$Class <- short_names[per_class$Class]
  
  per_class %>%
    pivot_wider(
      names_from = Class,
      values_from = c(Correct, Total, Recall, Precision, F1),
      names_glue = "{Class}_{.value}"
    )
}

warn_if_class_too_small <- function(y, min_safe_size = 10) {
  counts <- table(y)
  small_classes <- counts[counts < min_safe_size]
  if (length(small_classes) > 0) {
    cat("\n  NOTE: these classes have very few observations:\n")
    print(small_classes)
    cat("  Cross-validation folds and SMOTE neighbours will be limited for them -\n")
    cat("  check their individual confusion-matrix rows, not just overall accuracy.\n")
  }
}

# =============================================================================
# Helper: 5-fold (or fewer) SVM cross-validation, with optional balancing
# =============================================================================

run_svm_cv <- function(X, y, balance = "No_Balancing") {
  
  y <- droplevels(y)
  k <- 4   # fixed at 4 folds; with 4 eruption weeks, each fold gets exactly one
  if (min(table(y)) < k) {
    stop("Smallest class has fewer than ", k, " observations - cannot do ", k,
         "-fold cross-validation without an empty fold for that class.")
  }
  cat("  [run_svm_cv] using k =", k, "folds; smallest class size =", min(table(y)), "\n")
  
  folds <- createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  prediction <- factor(rep(NA_character_, length(y)), levels = levels(y))
  
  for (fold in folds) {
    train_idx <- setdiff(seq_along(y), fold)
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- droplevels(y[train_idx])
    
    balanced <- balance_training_data(X_train, y_train, balance)
    X_train <- balanced$X
    y_train <- balanced$y
    
    counts <- table(y_train)
    weights <- setNames(sum(counts) / (length(counts) * counts), names(counts))
    
    model <- svm(x = X_train, y = y_train,
                 kernel = "radial", cost = 1, gamma = 1 / ncol(X_train),
                 class.weights = weights, scale = FALSE)
    
    prediction[fold] <- as.character(predict(model, X[fold, , drop = FALSE]))
  }
  
  metrics <- class_metrics(y, prediction)
  
  list(prediction = prediction, accuracy = mean(prediction == y),
       balanced_accuracy = metrics$balanced_accuracy, per_class = metrics$per_class, folds = k)
}

# =============================================================================
# Helper: SVM across feature sets and balancing methods
# =============================================================================

compare_svm_across_settings <- function(data, station_name) {
  
  cat("\n---- SVM across feature sets and balancing methods ----\n")
  
  results <- tibble()
  
  for (set_name in names(feature_sets)) {
    
    cols <- feature_sets[[set_name]]
    X_set <- scale(data %>% select(all_of(cols)))
    
    for (method in balance_methods) {
      
      cat("  Feature set:", set_name, " | Balance method:", method, "\n")
      
      svm_res <- run_svm_cv(X_set, data$Model_Label, balance = method)
      
      row <- tibble(Station = station_name, Feature_Set = set_name, Balance_Method = method,
                    Accuracy = svm_res$accuracy, Balanced_Accuracy = svm_res$balanced_accuracy)
      
      results <- bind_rows(results, bind_cols(row, flatten_per_class(svm_res$per_class)))
    }
  }
  
  results <- results %>%
    mutate(Feature_Set = factor(Feature_Set, levels = names(feature_sets)),
           Balance_Method = factor(Balance_Method, levels = balance_methods))
  
  # the metric that matters most here: does the model actually catch real
  # eruption weeks, not just score well on average
  plot_data <- results %>%
    select(Feature_Set, Balance_Method, Accuracy, Balanced_Accuracy, Eruption_Recall) %>%
    pivot_longer(cols = c(Accuracy, Balanced_Accuracy, Eruption_Recall),
                 names_to = "Metric", values_to = "Value")
  
  svm_plot <- ggplot(plot_data, aes(x = Balance_Method, y = Value, fill = Feature_Set)) +
    geom_col(position = "dodge") +
    facet_wrap(~ Metric, ncol = 1) +
    labs(title = paste0(station_name, ": SVM Performance by Feature Set and Balancing Method"),
         x = "Balancing Method", y = NULL) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  print(svm_plot)
  plot_path <- file.path(CLUSTER_DIR, paste0(station_name, "_SVM_FeatureSet_Balance_Plot.pdf"))
  ggsave(plot_path, svm_plot, width = 9, height = 9)
  
  list(results = results, plot_path = plot_path)
}

# =============================================================================
# Helper: K-means and DBSCAN quality across feature sets
# =============================================================================
#
# Scores each result two ways: average silhouette width (an INTERNAL check -
# are the clusters tight and well separated?) and Adjusted Rand Index (an
# EXTERNAL check against the true VAL classes, correcting for chance
# agreement: 1 = perfect match, 0 = no better than random, negative = worse).
# =============================================================================

run_clustering_feature_set_comparison <- function(data, station_name) {
  
  cat("\n---- Clustering quality across feature sets ----\n")
  
  metrics <- tibble()
  
  for (set_name in names(clustering_feature_sets)) {
    
    cat("  Feature set:", set_name, "\n")
    
    cols <- clustering_feature_sets[[set_name]]
    X_set <- scale(data %>% select(all_of(cols)))
    
    set.seed(123456789, kind = "Mersenne-Twister")
    km <- kmeans(X_set, centers = 3, nstart = 50)
    km_sil <- mean(silhouette(km$cluster, dist(X_set))[, "sil_width"])
    km_ari <- mclust::adjustedRandIndex(km$cluster, data$Model_Label)
    
    db <- find_three_cluster_dbscan(X_set)
    db_sil <- mean(silhouette(db$model$cluster, dist(X_set))[, "sil_width"])
    db_ari <- mclust::adjustedRandIndex(db$model$cluster, data$Model_Label)
    
    metrics <- bind_rows(
      metrics,
      tibble(Station = station_name, Feature_Set = set_name, Method = "Kmeans",
             Avg_Silhouette = km_sil, Adjusted_Rand_Index = km_ari),
      tibble(Station = station_name, Feature_Set = set_name, Method = "DBSCAN",
             Avg_Silhouette = db_sil, Adjusted_Rand_Index = db_ari)
    )
  }
  
  metrics <- metrics %>% mutate(Feature_Set = factor(Feature_Set, levels = names(clustering_feature_sets)))
  
  metrics_long <- metrics %>%
    pivot_longer(cols = c(Avg_Silhouette, Adjusted_Rand_Index), names_to = "Metric", values_to = "Value")
  
  clustering_plot <- ggplot(metrics_long, aes(x = Method, y = Value, fill = Feature_Set)) +
    geom_col(position = "dodge") +
    facet_wrap(~ Metric, scales = "free_y") +
    labs(title = paste0(station_name, ": Clustering Quality by Feature Set"),
         x = "Clustering Method", y = NULL)
  
  print(clustering_plot)
  plot_path <- file.path(CLUSTER_DIR, paste0(station_name, "_Clustering_FeatureSet_Comparison.pdf"))
  ggsave(plot_path, clustering_plot, width = 9, height = 6)
  
  list(metrics = metrics, plot_path = plot_path)
}

# =============================================================================
# Helper: stack the three confusion tables (Kmeans, DBSCAN, SVM baseline)
# into one long table - combined across stations at the bottom of the script
# =============================================================================

combine_confusion_tables <- function(station_name, kmeans_table, dbscan_table, svm_table) {
  
  to_long <- function(tbl, method_name) {
    df <- as.data.frame(tbl)
    names(df)[2] <- "Predicted_Group"
    df$Method <- method_name
    df$Station <- station_name
    df
  }
  
  bind_rows(
    to_long(kmeans_table, "Kmeans"),
    to_long(dbscan_table, "DBSCAN"),
    to_long(svm_table,    "SVM")
  )
}

# =============================================================================
# Main per-station analysis
# =============================================================================

analyse_station <- function(station_name, file_name) {
  
  cat("\n\n==================== STATION:", station_name, "====================\n")
  
  data <- read_csv(file_name, show_col_types = FALSE) %>%
    mutate(Week_Start = as.POSIXct(Week_Start, tz = "UTC"),
           Week_End   = as.POSIXct(Week_End, tz = "UTC")) %>%
    filter(Week_Start >= ANALYSIS_START, Week_Start <= ANALYSIS_END)
  
  cat("Weekly observations:", nrow(data), "\n")
  
  cat("\nCorrelation of 13 original features...\n")
  make_correlogram(data, station_name)
  
  data$Model_Label_Num <- mapply(assign_val_class, data$Week_Start, data$Week_End)
  data <- data %>%
    mutate(Model_Label = factor(model_label_levels[Model_Label_Num], levels = model_label_levels)) %>%
    filter(!is.na(Model_Label), if_all(all_of(all_features), is.finite))
  
  cat("\nTrue VAL classes:\n"); print(table(data$Model_Label))
  warn_if_class_too_small(data$Model_Label)
  
  X <- scale(data %>% select(all_of(feature_cols)))
  pca <- prcomp(X)   # shared by both cluster plots below
  
  # ---- K-means -----------------------------------------------------------
  cat("\n---- K-means ----\n")
  
  set.seed(123456789, kind = "Mersenne-Twister")
  km_model <- kmeans(X, centers = 3, nstart = 50)
  data$Kmeans_Cluster <- factor(km_model$cluster)
  
  kmeans_table <- table(True_Class = data$Model_Label, Kmeans_Cluster = data$Kmeans_Cluster)
  cat("\nK-means vs true VAL class:\n"); print(kmeans_table)
  
  betweenss_ratio <- km_model$betweenss / km_model$totss
  km_sil <- mean(silhouette(km_model$cluster, dist(X))[, "sil_width"])
  cat("Between-SS / total-SS:", round(betweenss_ratio, 3), " Avg silhouette:", round(km_sil, 3), "\n")
  
  kmeans_plot <- make_cluster_plot(pca, data$Kmeans_Cluster, station_name, "K-means")
  print(kmeans_plot)
  ggsave(file.path(CLUSTER_DIR, paste0(station_name, "_Kmeans_Cluster_Plot.pdf")), kmeans_plot, width = 8, height = 6)
  
  # ---- DBSCAN --------------------------------------------------------------
  cat("\n---- DBSCAN ----\n")
  
  db_result <- find_three_cluster_dbscan(X)
  db_model <- db_result$model
  cat("eps =", db_result$eps, " minPts =", db_result$minPts, "\n")
  
  data$DBSCAN_Cluster <- factor(db_model$cluster, levels = 0:3,
                                labels = c("Noise", "Cluster 1", "Cluster 2", "Cluster 3"))
  
  dbscan_table <- table(True_Class = data$Model_Label, DBSCAN_Cluster = data$DBSCAN_Cluster)
  cat("\nDBSCAN vs true VAL class:\n"); print(dbscan_table)
  
  db_sil <- mean(silhouette(db_model$cluster, dist(X))[, "sil_width"])
  cat("Avg silhouette:", round(db_sil, 3), " Noise:", sum(db_model$cluster == 0), "of", nrow(data), "\n")
  
  dbscan_plot <- make_cluster_plot(pca, data$DBSCAN_Cluster, station_name, "DBSCAN")
  print(dbscan_plot)
  ggsave(file.path(CLUSTER_DIR, paste0(station_name, "_DBSCAN_Cluster_Plot.pdf")), dbscan_plot, width = 8, height = 6)
  
  # ---- clustering quality across feature sets ------------------------------
  
  clustering_comparison <- run_clustering_feature_set_comparison(data, station_name)
  
  # ---- SVM baseline (7 features, no balancing) -----------------------------
  cat("\n---- SVM baseline (7 features, no balancing) ----\n")
  
  svm_result <- run_svm_cv(X, data$Model_Label)
  data$SVM_Prediction <- svm_result$prediction
  
  cat("Accuracy:", round(svm_result$accuracy, 4),
      " Balanced accuracy:", round(svm_result$balanced_accuracy, 4), "\n")
  cat("Per-class (Correct / Total):\n"); print(svm_result$per_class %>% select(Class, Correct, Total, Recall, Precision))
  
  svm_table <- table(True_Class = data$Model_Label, SVM_Prediction = data$SVM_Prediction)
  cat("\nSVM confusion matrix:\n"); print(svm_table)
  
  # ---- SVM across feature sets and balancing methods -----------------------
  
  svm_comparison <- compare_svm_across_settings(data, station_name)
  
  # ---- results returned for combining across stations at the bottom -------
  
  list(
    data = data,
    kmeans_table = kmeans_table, dbscan_table = dbscan_table, svm_table = svm_table,
    clustering_comparison = clustering_comparison, svm_comparison = svm_comparison,
    summary = tibble(
      Station = station_name,
      Number_of_Weeks = nrow(data),
      Kmeans_BetweenSS_Ratio = betweenss_ratio,
      Kmeans_Avg_Silhouette = km_sil,
      DBSCAN_eps = db_result$eps,
      DBSCAN_minPts = db_result$minPts,
      DBSCAN_Avg_Silhouette = db_sil,
      DBSCAN_Noise_Percent = 100 * mean(db_model$cluster == 0),
      SVM_Baseline_Accuracy = svm_result$accuracy,
      SVM_Baseline_Balanced_Accuracy = svm_result$balanced_accuracy
    )
  )
}

# =============================================================================
# Run both stations
# =============================================================================

WIZ_result  <- analyse_station("WIZ",  station_files["WIZ"])
WSRZ_result <- analyse_station("WSRZ", station_files["WSRZ"])

# =============================================================================
# Combine everything into a small set of result files
# =============================================================================

cat("\n\n==================== Combining results ====================\n")

# 1. true class counts, both stations
true_class_distribution <- bind_rows(
  WIZ_result$data  %>% count(Model_Label, name = "Number_of_Weeks") %>% mutate(Station = "WIZ"),
  WSRZ_result$data %>% count(Model_Label, name = "Number_of_Weeks") %>% mutate(Station = "WSRZ")
)
write_csv(true_class_distribution, file.path(CLUSTER_DIR, "WIZ_WSRZ_True_Class_Distribution.csv"))

# 2. one summary row per station (K-means, DBSCAN, SVM baseline)
model_summary <- bind_rows(WIZ_result$summary, WSRZ_result$summary)
write_csv(model_summary, file.path(CLUSTER_DIR, "WIZ_WSRZ_Model_Summary.csv"))

# 3. every confusion table (Kmeans, DBSCAN, SVM baseline), both stations
combined_confusion <- bind_rows(
  combine_confusion_tables("WIZ",  WIZ_result$kmeans_table,  WIZ_result$dbscan_table,  WIZ_result$svm_table),
  combine_confusion_tables("WSRZ", WSRZ_result$kmeans_table, WSRZ_result$dbscan_table, WSRZ_result$svm_table)
)
write_csv(combined_confusion, file.path(CLUSTER_DIR, "WIZ_WSRZ_Confusion_Tables_Combined.csv"))

# 4. SVM across feature sets and balancing methods, WITH per-class recall/precision -
#    this is the main table to sort by Eruption_Correct / Eruption_Recall
svm_comparison_combined <- bind_rows(WIZ_result$svm_comparison$results, WSRZ_result$svm_comparison$results)
write_csv(svm_comparison_combined, file.path(CLUSTER_DIR, "WIZ_WSRZ_SVM_FeatureSet_Balance_Comparison.csv"))

# 5. clustering quality (Kmeans, DBSCAN) across feature sets, both stations
clustering_comparison_combined <- bind_rows(WIZ_result$clustering_comparison$metrics, WSRZ_result$clustering_comparison$metrics)
write_csv(clustering_comparison_combined, file.path(CLUSTER_DIR, "WIZ_WSRZ_Clustering_FeatureSet_Comparison.csv"))

# 6. full weekly results (true class, cluster labels, SVM prediction), both stations
weekly_results <- bind_rows(
  WIZ_result$data  %>% mutate(Station = "WIZ"),
  WSRZ_result$data %>% mutate(Station = "WSRZ")
) %>%
  select(Station, Week_Start, Week_End, Model_Label, Kmeans_Cluster, DBSCAN_Cluster, SVM_Prediction, all_of(all_features))
write_csv(weekly_results, file.path(CLUSTER_DIR, "WIZ_WSRZ_Weekly_Results.csv"))

cat("\nResult files written to:", CLUSTER_DIR, "\n")
cat(" - WIZ_WSRZ_True_Class_Distribution.csv\n")
cat(" - WIZ_WSRZ_Model_Summary.csv\n")
cat(" - WIZ_WSRZ_Confusion_Tables_Combined.csv\n")
cat(" - WIZ_WSRZ_SVM_FeatureSet_Balance_Comparison.csv   <- sort by Eruption_Correct / Eruption_Recall\n")
cat(" - WIZ_WSRZ_Clustering_FeatureSet_Comparison.csv\n")
cat(" - WIZ_WSRZ_Weekly_Results.csv\n")
cat("\nPlots (each its own PDF, per station):\n")
cat(" - <station>_13_Feature_Correlogram.pdf\n")
cat(" - <station>_Kmeans_Cluster_Plot.pdf\n")
cat(" - <station>_DBSCAN_Cluster_Plot.pdf\n")
cat(" - <station>_Clustering_FeatureSet_Comparison.pdf\n")
cat(" - <station>_SVM_FeatureSet_Balance_Plot.pdf\n")

cat("\nAnalysis complete.\n")
