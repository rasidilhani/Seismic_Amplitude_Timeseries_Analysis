# =============================================================================
# WIZ and WSRZ seismic activity classification
#
# Steps:
#   - correlation check across all 13 original features
#   - K-means, DBSCAN, and ensemble (consensus) clustering
#   - SVM, Random Forest, and Naive Bayes classification, with cross-validation
#   - class-imbalance handling: oversampling, undersampling, SMOTE
#   - comparison across three feature sets: 13 features, 7 features,
#     and a variance-only 4-feature set
#   - all result tables and all plots combined into single files
# =============================================================================

library(tidyverse)
library(corrplot)
library(dbscan)
library(e1071)      # SVM and Naive Bayes
library(caret)
library(cluster)      # for silhouette
library(randomForest)
library(qpdf)         # for combining PDF plots into one file

if (!requireNamespace("mclust", quietly = TRUE)) {
  stop("Package 'mclust' is needed for the Adjusted Rand Index. Install it with: install.packages('mclust')")
}
# mclust is used only via mclust::adjustedRandIndex() below, not library()'d in full,
# because it overwrites some tidyverse/purrr function names (e.g. map)

set.seed(123456789, kind = "Mersenne-Twister")

# --- plot style: serif font, article-style look, applied to every plot -----
#
# theme_set() makes this the default for every ggplot() call from here on,
# so individual plots don't need to repeat it.
# =============================================================================

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
CLUSTER_DIR <- file.path(RESULTS_DIR, "cluster_analysis_R2")
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

# all 13 original features
all_features <- c(
  "H_Shannon", "C_Shannon", "Fisher_Info", "C_Fisher",
  "H_Renyi", "C_Renyi", "H_Tsallis", "C_Tsallis",
  "Var_H_Shannon", "Var_C_Shannon", "Var_H_Renyi", "Var_H_Tsallis",
  "Disequilibrium"
)

# the 7 features chosen after the correlation check (used for clustering)
feature_cols <- c(
  "Fisher_Info", "C_Fisher",
  "Var_H_Shannon", "Var_C_Shannon",
  "Var_H_Tsallis", "Var_H_Renyi",
  "Disequilibrium"
)

# just the variance-based features, to see if they carry most of the signal
variance_features <- c("Var_H_Shannon", "Var_C_Shannon", "Var_H_Renyi", "Var_H_Tsallis")

# the three feature sets we compare for classification accuracy
feature_sets <- list(
  "13_Features"       = all_features,
  "7_Features"        = feature_cols,
  "Variance_Features" = variance_features
)

# the feature sets we compare for CLUSTERING quality (kept smaller than the
# classification comparison above, since clustering is much slower to run -
# DBSCAN in particular re-searches its whole eps/minPts grid every time)
clustering_feature_sets <- list(
  "7_Features"        = feature_cols,
  "Variance_Features" = variance_features
)

# the class-imbalance methods we compare
balance_methods <- c("none", "oversample", "undersample", "smote")

model_label_levels <- c(
  "1_Background_unrest",
  "2_Moderate_to_heightened_unrest",
  "3_Minor_to_moderate_eruption"
)

# =============================================================================
# VAL activity data: read it, keep version 3, collapse to 3 classes
#
# Using the UTC period columns here (not the NZ local-time ones) so this
# lines up with the OP feature data, which was computed in UTC.
# =============================================================================

val_df <- read_csv(VAL_FILE, show_col_types = FALSE) %>%
  mutate(
    PeriodStart_UTC = as.POSIXct(str_remove(`PeriodStart UTC`, "\\s*UTC$"),
                                 format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    PeriodEnd_UTC   = as.POSIXct(str_remove(`PeriodEnd UTC`, "\\s*UTC$"),
                                 format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  )

# open-ended VAL periods just run to the end of the analysis window
open_ended <- is.na(val_df$PeriodEnd_UTC) | str_starts(val_df$`PeriodEnd UTC`, "9999")
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

# if a week overlaps more than one VAL period, use the more severe one
assign_val_class <- function(week_start, week_end) {
  overlap <- val_df %>% filter(PeriodStart_UTC < week_end, PeriodEnd_UTC > week_start)
  if (nrow(overlap) == 0) return(NA_integer_)
  max(overlap$Model_Label_Num, na.rm = TRUE)
}

# =============================================================================
# Helper: 13-feature correlogram (plot + saved PDF + CSV)
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
  
  pdf_path <- file.path(CLUSTER_DIR, paste0(station_name, "_13_Feature_Correlogram.pdf"))
  pdf(pdf_path, width = 11, height = 9)
  draw_plot()
  dev.off()
  
  #write.csv(cor_matrix, file.path(CLUSTER_DIR, paste0(station_name, "_13_Feature_Correlation_Matrix.csv")))
  
  list(matrix = cor_matrix, pdf_path = pdf_path)
}

# =============================================================================
# Helper: build one PCA-based cluster plot, for any clustering method
# =============================================================================
#
# All three clustering methods (K-means, DBSCAN, ensemble) run on the same
# scaled feature matrix X, so they can all share one PCA projection. The
# axis labels show how much of the total variance each PC captures, so the
# plot is honest about how much of the picture a 2D view actually shows.
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
  
  list(model = final_fit, eps = best$eps, minPts = best$minPts, search = grid)
}

# =============================================================================
# Helper: ensemble ("consensus") clustering
# =============================================================================
#
# A single K-means run can land slightly differently depending on its random
# starting points. Ensemble clustering makes the result more stable by:
#
#   1. Running K-means many times, each with a different random seed.
#   2. For every pair of weeks, counting how often the two runs put them in
#      the same cluster -> a co-occurrence score between 0 and 1 per pair.
#   3. Treating (1 - co-occurrence) as a distance and running one final
#      hierarchical clustering pass on it, cut into 3 groups.
#
# The result is the consensus: the grouping most of the individual K-means
# runs agreed on, rather than the outcome of just one run.
# =============================================================================

run_ensemble_clustering <- function(X, k = 3, n_runs = 50) {
  
  n <- nrow(X)
  co_occurrence <- matrix(0, n, n)
  
  for (run in seq_len(n_runs)) {
    set.seed(run)
    fit <- kmeans(X, centers = k, nstart = 1)
    same_cluster <- outer(fit$cluster, fit$cluster, `==`)
    co_occurrence <- co_occurrence + same_cluster
  }
  
  co_occurrence <- co_occurrence / n_runs
  
  consensus_dist <- as.dist(1 - co_occurrence)
  consensus_tree <- hclust(consensus_dist, method = "average")
  consensus_cluster <- cutree(consensus_tree, k = k)
  
  list(cluster = consensus_cluster, co_occurrence = co_occurrence, n_runs = n_runs)
}

# =============================================================================
# Class-imbalance methods
# =============================================================================
#
# All three work on one training fold at a time: X (a feature matrix) and
# y (the true class for each row). They return a new, more balanced
# X and y to train on. The test fold is never touched by these functions.
# =============================================================================

# pick one value safely - sample() misbehaves if you give it a vector of length 1
pick_one <- function(x) {
  if (length(x) == 1) return(x)
  sample(x, 1)
}

# duplicate rows from the smaller classes until every class matches the
# largest class
oversample_data <- function(X, y) {
  target_n <- max(table(y))
  idx <- unlist(lapply(levels(y), function(cls) {
    cls_idx <- which(y == cls)
    sample(cls_idx, target_n, replace = TRUE)
  }))
  list(X = X[idx, , drop = FALSE], y = y[idx])
}

# drop rows from the bigger classes until every class matches the smallest
undersample_data <- function(X, y) {
  target_n <- min(table(y))
  idx <- unlist(lapply(levels(y), function(cls) {
    cls_idx <- which(y == cls)
    sample(cls_idx, target_n, replace = FALSE)
  }))
  list(X = X[idx, , drop = FALSE], y = y[idx])
}

# SMOTE: instead of copying minority rows, invent new ones by interpolating
# between a minority point and one of its nearest neighbours in the same class
smote_data <- function(X, y, k = 5) {
  
  target_n <- max(table(y))
  X_new <- X
  y_new <- as.character(y)
  
  for (cls in levels(y)) {
    
    cls_idx <- which(y == cls)
    n_cls <- length(cls_idx)
    n_needed <- target_n - n_cls
    if (n_needed <= 0) next   # this class is already the biggest one
    
    cls_X <- X[cls_idx, , drop = FALSE]
    dist_matrix <- as.matrix(dist(cls_X))
    
    synthetic <- matrix(NA_real_, nrow = n_needed, ncol = ncol(X))
    colnames(synthetic) <- colnames(X)
    
    for (i in seq_len(n_needed)) {
      
      point_idx <- pick_one(seq_len(n_cls))
      point <- cls_X[point_idx, ]
      
      # nearest neighbours of this point, within the same class, excluding itself
      neighbour_order <- order(dist_matrix[point_idx, ])
      neighbour_order <- neighbour_order[neighbour_order != point_idx]
      neighbours <- head(neighbour_order, min(k, length(neighbour_order)))
      
      neighbour_idx <- pick_one(neighbours)
      neighbour <- cls_X[neighbour_idx, ]
      
      # a new point somewhere between "point" and its neighbour
      gap <- runif(1)
      synthetic[i, ] <- point + gap * (neighbour - point)
    }
    
    X_new <- rbind(X_new, synthetic)
    y_new <- c(y_new, rep(cls, n_needed))
  }
  
  list(X = X_new, y = factor(y_new, levels = levels(y)))
}

# picks the right balancing function by name
balance_training_data <- function(X, y, method) {
  switch(method,
         none        = list(X = X, y = y),
         oversample  = oversample_data(X, y),
         undersample = undersample_data(X, y),
         smote       = smote_data(X, y),
         stop("Unknown balancing method: ", method)
  )
}

# =============================================================================
# Helper: per-class metrics and balanced accuracy
# =============================================================================
#
# Plain accuracy is misleading when one class has far more rows than the
# others (here: 320 / 12 / 3). A model that just guesses the biggest class
# every time would score ~95% accuracy while being useless for the rare
# classes. Balanced accuracy instead averages the RECALL of each class, so
# all three classes count equally regardless of how many rows they have.
# =============================================================================

class_metrics <- function(true, predicted) {
  
  classes <- levels(true)
  per_class <- tibble(Class = classes, Recall = NA_real_, Precision = NA_real_, F1 = NA_real_)
  
  for (i in seq_along(classes)) {
    cls <- classes[i]
    
    true_positive  <- sum(true == cls & predicted == cls)
    actual_count   <- sum(true == cls)
    predicted_count <- sum(predicted == cls)
    
    recall    <- if (actual_count > 0) true_positive / actual_count else NA_real_
    precision <- if (predicted_count > 0) true_positive / predicted_count else NA_real_
    f1 <- if (!is.na(recall) && !is.na(precision) && (recall + precision) > 0) {
      2 * recall * precision / (recall + precision)
    } else {
      NA_real_
    }
    
    per_class$Recall[i]    <- recall
    per_class$Precision[i] <- precision
    per_class$F1[i]        <- f1
  }
  
  list(per_class = per_class, balanced_accuracy = mean(per_class$Recall, na.rm = TRUE))
}

# warns when a class is so small that cross-validation and SMOTE only have
# a handful of examples to work with
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

run_svm_cv <- function(X, y, balance = "none") {
  
  y <- droplevels(y)
  k <- min(5, min(table(y)))
  if (k < 2) stop("Not enough observations per class for cross-validation.")
  
  folds <- createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  prediction <- factor(rep(NA_character_, length(y)), levels = levels(y))
  
  for (fold in folds) {
    train_idx <- setdiff(seq_along(y), fold)
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- droplevels(y[train_idx])
    
    # balance the TRAINING fold only - the test fold stays untouched
    balanced <- balance_training_data(X_train, y_train, balance)
    X_train <- balanced$X
    y_train <- balanced$y
    
    # inverse-frequency class weights, on top of whatever balancing was done
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
# Helper: 5-fold (or fewer) Random Forest cross-validation, with balancing
# =============================================================================

run_rf_cv <- function(X, y, balance = "none", ntree = 500) {
  
  y <- droplevels(y)
  k <- min(5, min(table(y)))
  if (k < 2) stop("Not enough observations per class for cross-validation.")
  
  folds <- createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  prediction <- factor(rep(NA_character_, length(y)), levels = levels(y))
  importance_sum <- NULL
  
  for (fold in folds) {
    train_idx <- setdiff(seq_along(y), fold)
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- droplevels(y[train_idx])
    
    balanced <- balance_training_data(X_train, y_train, balance)
    X_train <- balanced$X
    y_train <- balanced$y
    
    if (balance == "none") {
      # data is still imbalanced here, so make each tree sample every
      # class evenly instead of letting the majority class dominate
      min_class_size <- min(table(y_train))
      model <- randomForest(x = X_train, y = y_train, ntree = ntree,
                            sampsize = rep(min_class_size, length(levels(y_train))),
                            strata = y_train, importance = TRUE)
    } else {
      # already balanced above, so just fit normally
      model <- randomForest(x = X_train, y = y_train, ntree = ntree, importance = TRUE)
    }
    
    prediction[fold] <- as.character(predict(model, X[fold, , drop = FALSE]))
    
    imp <- importance(model)[, "MeanDecreaseGini"]
    importance_sum <- if (is.null(importance_sum)) imp else importance_sum + imp
  }
  
  metrics <- class_metrics(y, prediction)
  
  list(
    prediction = prediction,
    accuracy = mean(prediction == y),
    balanced_accuracy = metrics$balanced_accuracy,
    per_class = metrics$per_class,
    folds = k,
    importance = sort(importance_sum / length(folds), decreasing = TRUE)
  )
}

# =============================================================================
# Helper: 5-fold (or fewer) Naive Bayes cross-validation, with balancing
# =============================================================================
#
# Naive Bayes is a GENERATIVE model: it learns, for each class, what a
# "typical" feature vector looks like (assuming the features are roughly
# Gaussian and independent of each other), then classifies a new week by
# asking "which class's typical pattern does this look most like?" via
# Bayes' rule. This is the opposite approach to SVM and Random Forest, which
# are DISCRIMINATIVE - they only learn the boundary between classes, without
# modelling what each class's data actually looks like.
# =============================================================================

run_nb_cv <- function(X, y, balance = "none") {
  
  y <- droplevels(y)
  k <- min(5, min(table(y)))
  if (k < 2) stop("Not enough observations per class for cross-validation.")
  
  folds <- createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  prediction <- factor(rep(NA_character_, length(y)), levels = levels(y))
  
  for (fold in folds) {
    train_idx <- setdiff(seq_along(y), fold)
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- droplevels(y[train_idx])
    
    balanced <- balance_training_data(X_train, y_train, balance)
    X_train <- balanced$X
    y_train <- balanced$y
    
    model <- naiveBayes(x = X_train, y = y_train)
    
    prediction[fold] <- as.character(predict(model, X[fold, , drop = FALSE]))
  }
  
  metrics <- class_metrics(y, prediction)
  
  list(prediction = prediction, accuracy = mean(prediction == y),
       balanced_accuracy = metrics$balanced_accuracy, per_class = metrics$per_class, folds = k)
}

# =============================================================================
# Helper: compare feature sets x balancing methods, for all three classifiers
# =============================================================================

compare_feature_sets_and_balancing <- function(data, station_name) {
  
  cat("\n---- Comparing feature sets and class-imbalance methods ----\n")
  
  results <- tibble()
  
  for (set_name in names(feature_sets)) {
    
    cols <- feature_sets[[set_name]]
    X_set <- scale(data %>% select(all_of(cols)))
    
    for (method in balance_methods) {
      
      cat("  Feature set:", set_name, " | Balance method:", method, "\n")
      
      svm_res <- run_svm_cv(X_set, data$Model_Label, balance = method)
      rf_res  <- run_rf_cv(X_set, data$Model_Label, balance = method)
      nb_res  <- run_nb_cv(X_set, data$Model_Label, balance = method)
      
      results <- bind_rows(
        results,
        tibble(Station = station_name, Feature_Set = set_name, Balance_Method = method,
               Classifier = "SVM", Accuracy = svm_res$accuracy,
               Balanced_Accuracy = svm_res$balanced_accuracy),
        tibble(Station = station_name, Feature_Set = set_name, Balance_Method = method,
               Classifier = "Random_Forest", Accuracy = rf_res$accuracy,
               Balanced_Accuracy = rf_res$balanced_accuracy),
        tibble(Station = station_name, Feature_Set = set_name, Balance_Method = method,
               Classifier = "Naive_Bayes", Accuracy = nb_res$accuracy,
               Balanced_Accuracy = nb_res$balanced_accuracy)
      )
    }
  }
  
  results <- results %>%
    mutate(Feature_Set = factor(Feature_Set, levels = names(feature_sets)),
           Balance_Method = factor(Balance_Method, levels = balance_methods))
  
  # look at Balanced_Accuracy, not Accuracy, to judge which method actually
  # helps the rare classes - plain accuracy barely moves because the biggest
  # class (320 of 335 rows) dominates it either way
  results_long <- results %>%
    pivot_longer(cols = c(Accuracy, Balanced_Accuracy), names_to = "Metric", values_to = "Value")
  
  accuracy_plot <- ggplot(results_long, aes(x = Balance_Method, y = Value, fill = Classifier)) +
    geom_col(position = "dodge") +
    facet_grid(Metric ~ Feature_Set) +
    labs(title = paste0(station_name, ": Accuracy vs Balanced Accuracy by Feature Set and Balancing Method"),
         x = "Balancing Method", y = NULL) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  print(accuracy_plot)
  
  plot_path <- file.path(CLUSTER_DIR, paste0(station_name, "_FeatureSet_Balance_Accuracy.pdf"))
  ggsave(plot_path, accuracy_plot, width = 12, height = 7)
  
  write_csv(results, file.path(CLUSTER_DIR, paste0(station_name, "_FeatureSet_Balance_Accuracy.csv")))
  
  list(results = results, plot_path = plot_path)
}

# =============================================================================
# Helper: compare clustering quality across feature sets
# =============================================================================
#
# Runs K-means, DBSCAN, and Ensemble clustering separately for each feature
# set, and scores each result two ways:
#
#   - Average silhouette width: an INTERNAL check - are the clusters
#     themselves tight and well separated? Doesn't use the true VAL classes.
#   - Adjusted Rand Index (ARI): an EXTERNAL check against the true VAL
#     classes, correcting for chance agreement. 1 = perfect match to the
#     true classes, 0 = no better than a random grouping, negative = worse
#     than random.
#
# ARI is the more direct answer to "does this feature set improve
# clustering", since it's the one that actually checks against ground truth.
# =============================================================================

run_clustering_feature_set_comparison <- function(data, station_name) {
  
  cat("\n---- Comparing clustering quality across feature sets ----\n")
  
  metrics <- tibble()
  
  for (set_name in names(clustering_feature_sets)) {
    
    cat("  Feature set:", set_name, "\n")
    
    cols <- clustering_feature_sets[[set_name]]
    X_set <- scale(data %>% select(all_of(cols)))
    
    # K-means
    set.seed(123456789, kind = "Mersenne-Twister")
    km <- kmeans(X_set, centers = 3, nstart = 50)
    km_sil <- mean(silhouette(km$cluster, dist(X_set))[, "sil_width"])
    km_ari <- mclust::adjustedRandIndex(km$cluster, data$Model_Label)
    
    # DBSCAN
    db <- find_three_cluster_dbscan(X_set)
    db_sil <- mean(silhouette(db$model$cluster, dist(X_set))[, "sil_width"])
    db_ari <- mclust::adjustedRandIndex(db$model$cluster, data$Model_Label)
    
    # Ensemble
    ens <- run_ensemble_clustering(X_set, k = 3, n_runs = 50)
    ens_sil <- mean(silhouette(ens$cluster, dist(X_set))[, "sil_width"])
    ens_ari <- mclust::adjustedRandIndex(ens$cluster, data$Model_Label)
    
    metrics <- bind_rows(
      metrics,
      tibble(Station = station_name, Feature_Set = set_name, Method = "Kmeans",
             Avg_Silhouette = km_sil, Adjusted_Rand_Index = km_ari),
      tibble(Station = station_name, Feature_Set = set_name, Method = "DBSCAN",
             Avg_Silhouette = db_sil, Adjusted_Rand_Index = db_ari),
      tibble(Station = station_name, Feature_Set = set_name, Method = "Ensemble",
             Avg_Silhouette = ens_sil, Adjusted_Rand_Index = ens_ari)
    )
  }
  
  metrics <- metrics %>%
    mutate(Feature_Set = factor(Feature_Set, levels = names(clustering_feature_sets)))
  
  metrics_long <- metrics %>%
    pivot_longer(cols = c(Avg_Silhouette, Adjusted_Rand_Index), names_to = "Metric", values_to = "Value")
  
  clustering_plot <- ggplot(metrics_long, aes(x = Method, y = Value, fill = Feature_Set)) +
    geom_col(position = "dodge") +
    facet_wrap(~ Metric, scales = "free_y") +
    labs(title = paste0(station_name, ": Clustering Quality by Feature Set"),
         x = "Clustering Method", y = NULL)
  
  print(clustering_plot)
  
  plot_path <- file.path(CLUSTER_DIR, paste0(station_name, "_Clustering_FeatureSet_Comparison.pdf"))
  ggsave(plot_path, clustering_plot, width = 10, height = 6)
  
  write_csv(metrics, file.path(CLUSTER_DIR, paste0(station_name, "_Clustering_FeatureSet_Comparison.csv")))
  
  list(metrics = metrics, plot_path = plot_path)
}

# =============================================================================
# Helper: stack every confusion table into one long, combined CSV
# =============================================================================

combine_confusion_tables <- function(station_name, kmeans_table, dbscan_table,
                                     ensemble_table, svm_table, rf_table, nb_table) {
  
  to_long <- function(tbl, method_name) {
    df <- as.data.frame(tbl)
    names(df)[2] <- "Predicted_Group"   # 2nd column name differs per method
    df$Method <- method_name
    df
  }
  
  combined <- bind_rows(
    to_long(kmeans_table,   "Kmeans"),
    to_long(dbscan_table,   "DBSCAN"),
    to_long(ensemble_table, "Ensemble"),
    to_long(svm_table,      "SVM"),
    to_long(rf_table,       "Random_Forest"),
    to_long(nb_table,       "Naive_Bayes")
  )
  
  out_path <- file.path(CLUSTER_DIR, paste0(station_name, "_All_Confusion_Tables_Combined.csv"))
  write_csv(combined, out_path)
  
  list(combined = combined, csv_path = out_path)
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
  corr_result <- make_correlogram(data, station_name)
  
  # true VAL class per week
  data$Model_Label_Num <- mapply(assign_val_class, data$Week_Start, data$Week_End)
  data <- data %>%
    mutate(Model_Label = factor(model_label_levels[Model_Label_Num], levels = model_label_levels)) %>%
    filter(!is.na(Model_Label), if_all(all_of(all_features), is.finite))
  
  cat("\nTrue VAL classes:\n"); print(table(data$Model_Label))
  warn_if_class_too_small(data$Model_Label)
  
  cluster_data <- data %>% select(all_of(feature_cols))
  X <- scale(cluster_data)
  
  # one PCA, shared by every cluster plot below (K-means, DBSCAN, Ensemble)
  pca <- prcomp(X)
  
  # ---- K-means -----------------------------------------------------------
  cat("\n---- K-means ----\n")
  
  set.seed(123456789, kind = "Mersenne-Twister")
  km_model <- kmeans(X, centers = 3, nstart = 50)
  data$Kmeans_Cluster <- factor(km_model$cluster)
  
  cat("Cluster sizes:\n"); print(table(data$Kmeans_Cluster))
  
  kmeans_table <- table(True_Class = data$Model_Label, Kmeans_Cluster = data$Kmeans_Cluster)
  cat("\nK-means vs true VAL class:\n"); print(kmeans_table)
  
  cat("\nCluster centres (scaled units):\n"); print(km_model$centers)
  
  centres_original <- sweep(km_model$centers, 2, attr(X, "scaled:scale"), `*`)
  centres_original <- sweep(centres_original, 2, attr(X, "scaled:center"), `+`)
  cat("\nCluster centres (original units):\n"); print(round(centres_original, 4))
  
  betweenss_ratio <- km_model$betweenss / km_model$totss
  sil <- silhouette(km_model$cluster, dist(X))
  avg_sil <- mean(sil[, "sil_width"])
  
  cat("\nBetween-SS / total-SS:", round(betweenss_ratio, 3), "\n")
  cat("Average silhouette width:", round(avg_sil, 3), "\n")
  
  kmeans_plot <- make_cluster_plot(pca, data$Kmeans_Cluster, station_name, "K-means")
  print(kmeans_plot)
  kmeans_plot_path <- file.path(CLUSTER_DIR, paste0(station_name, "_Kmeans_Cluster_Plot.pdf"))
  ggsave(kmeans_plot_path, kmeans_plot, width = 8, height = 6)
  
  write.csv(round(centres_original, 4),
            file.path(CLUSTER_DIR, paste0(station_name, "_Kmeans_Cluster_Centres.csv")))
  
  # ---- DBSCAN --------------------------------------------------------------
  cat("\n---- DBSCAN ----\n")
  
  db_result <- find_three_cluster_dbscan(X)
  db_model <- db_result$model
  cat("eps =", db_result$eps, " minPts =", db_result$minPts, "\n")
  
  data$DBSCAN_Cluster <- factor(db_model$cluster, levels = 0:3,
                                labels = c("Noise", "Cluster 1", "Cluster 2", "Cluster 3"))
  
  cat("\nCluster sizes:\n"); print(table(data$DBSCAN_Cluster))
  
  dbscan_table <- table(True_Class = data$Model_Label, DBSCAN_Cluster = data$DBSCAN_Cluster)
  cat("\nDBSCAN vs true VAL class:\n"); print(dbscan_table)
  
  dbscan_plot <- make_cluster_plot(pca, data$DBSCAN_Cluster, station_name, "DBSCAN")
  print(dbscan_plot)
  dbscan_plot_path <- file.path(CLUSTER_DIR, paste0(station_name, "_DBSCAN_Cluster_Plot.pdf"))
  ggsave(dbscan_plot_path, dbscan_plot, width = 8, height = 6)
  
  # ---- Ensemble (consensus) clustering ------------------------------------
  cat("\n---- Ensemble clustering (50 K-means runs, consensus) ----\n")
  
  ensemble_result <- run_ensemble_clustering(X, k = 3, n_runs = 50)
  data$Ensemble_Cluster <- factor(ensemble_result$cluster)
  
  cat("Cluster sizes:\n"); print(table(data$Ensemble_Cluster))
  
  ensemble_table <- table(True_Class = data$Model_Label, Ensemble_Cluster = data$Ensemble_Cluster)
  cat("\nEnsemble clustering vs true VAL class:\n"); print(ensemble_table)
  
  ensemble_sil <- silhouette(ensemble_result$cluster, dist(X))
  ensemble_avg_sil <- mean(ensemble_sil[, "sil_width"])
  cat("Average silhouette width:", round(ensemble_avg_sil, 3), "\n")
  
  ensemble_plot <- make_cluster_plot(pca, data$Ensemble_Cluster, station_name, "Ensemble")
  print(ensemble_plot)
  ensemble_plot_path <- file.path(CLUSTER_DIR, paste0(station_name, "_Ensemble_Cluster_Plot.pdf"))
  ggsave(ensemble_plot_path, ensemble_plot, width = 8, height = 6)
  
  # ---- clustering quality across feature sets (7 features vs variance-only) --
  
  clustering_comparison <- run_clustering_feature_set_comparison(data, station_name)
  
  # ---- SVM (baseline, 7 features, no balancing) ---------------------------
  cat("\n---- SVM (3-class, 5-fold CV) ----\n")
  
  svm_result <- run_svm_cv(X, data$Model_Label)
  data$SVM_Prediction <- svm_result$prediction
  
  cat("Folds:", svm_result$folds,
      " Accuracy:", round(svm_result$accuracy, 4),
      " Balanced accuracy:", round(svm_result$balanced_accuracy, 4), "\n")
  cat("\nPer-class recall / precision / F1:\n"); print(svm_result$per_class)
  
  svm_table <- table(True_Class = data$Model_Label, SVM_Prediction = data$SVM_Prediction)
  cat("\nSVM confusion matrix:\n"); print(svm_table)
  
  # ---- Random Forest (baseline, 7 features, no balancing) -----------------
  cat("\n---- Random Forest (3-class, 5-fold CV) ----\n")
  
  rf_result <- run_rf_cv(X, data$Model_Label)
  data$RF_Prediction <- rf_result$prediction
  
  cat("Folds:", rf_result$folds,
      " Accuracy:", round(rf_result$accuracy, 4),
      " Balanced accuracy:", round(rf_result$balanced_accuracy, 4), "\n")
  cat("\nPer-class recall / precision / F1:\n"); print(rf_result$per_class)
  cat("\nFeature importance (Mean Decrease Gini, averaged over folds):\n")
  print(round(rf_result$importance, 2))
  
  rf_table <- table(True_Class = data$Model_Label, RF_Prediction = data$RF_Prediction)
  cat("\nRandom Forest confusion matrix:\n"); print(rf_table)
  
  # ---- Naive Bayes (baseline, 7 features, no balancing) --------------------
  cat("\n---- Naive Bayes (3-class, 5-fold CV) ----\n")
  
  nb_result <- run_nb_cv(X, data$Model_Label)
  data$NB_Prediction <- nb_result$prediction
  
  cat("Folds:", nb_result$folds,
      " Accuracy:", round(nb_result$accuracy, 4),
      " Balanced accuracy:", round(nb_result$balanced_accuracy, 4), "\n")
  cat("\nPer-class recall / precision / F1:\n"); print(nb_result$per_class)
  
  nb_table <- table(True_Class = data$Model_Label, NB_Prediction = data$NB_Prediction)
  cat("\nNaive Bayes confusion matrix:\n"); print(nb_table)
  
  # ---- feature set x balancing method comparison --------------------------
  
  feature_comparison <- compare_feature_sets_and_balancing(data, station_name)
  
  # ---- combine all the confusion tables into one CSV ----------------------
  
  combined_tables <- combine_confusion_tables(
    station_name, kmeans_table, dbscan_table, ensemble_table, svm_table, rf_table, nb_table
  )
  
  # ---- save the individual result tables (kept for detailed checking) -----
  
  write.csv(as.data.frame.matrix(kmeans_table),   file.path(CLUSTER_DIR, paste0(station_name, "_Kmeans_vs_VAL.csv")))
  write.csv(as.data.frame.matrix(dbscan_table),   file.path(CLUSTER_DIR, paste0(station_name, "_DBSCAN_vs_VAL.csv")))
  write.csv(as.data.frame.matrix(ensemble_table), file.path(CLUSTER_DIR, paste0(station_name, "_Ensemble_vs_VAL.csv")))
  write.csv(as.data.frame.matrix(svm_table),      file.path(CLUSTER_DIR, paste0(station_name, "_SVM_vs_VAL.csv")))
  write.csv(as.data.frame.matrix(rf_table),       file.path(CLUSTER_DIR, paste0(station_name, "_RF_vs_VAL.csv")))
  write.csv(as.data.frame.matrix(nb_table),       file.path(CLUSTER_DIR, paste0(station_name, "_NB_vs_VAL.csv")))
  write.csv(rf_result$importance,                 file.path(CLUSTER_DIR, paste0(station_name, "_RF_Feature_Importance.csv")))
  write_csv(svm_result$per_class %>% mutate(Classifier = "SVM") %>%
              bind_rows(rf_result$per_class %>% mutate(Classifier = "Random_Forest")) %>%
              bind_rows(nb_result$per_class %>% mutate(Classifier = "Naive_Bayes")),
            file.path(CLUSTER_DIR, paste0(station_name, "_PerClass_Metrics.csv")))
  write_csv(data, file.path(CLUSTER_DIR, paste0(station_name, "_Weekly_Cluster_Results.csv")))
  
  # ---- combine every plot for this station into one PDF --------------------
  
  station_plot_paths <- c(corr_result$pdf_path, kmeans_plot_path, dbscan_plot_path,
                          ensemble_plot_path, clustering_comparison$plot_path,
                          feature_comparison$plot_path)
  combined_pdf_path <- file.path(CLUSTER_DIR, paste0(station_name, "_All_Plots_Combined.pdf"))
  pdf_combine(input = station_plot_paths, output = combined_pdf_path)
  cat("\nAll plots for", station_name, "combined into:", combined_pdf_path, "\n")
  
  # ---- summary row ----------------------------------------------------------
  
  summary_row <- tibble(
    Station = station_name,
    Number_of_Weeks = nrow(data),
    Kmeans_BetweenSS_Ratio = betweenss_ratio,
    Kmeans_Avg_Silhouette = avg_sil,
    DBSCAN_eps = db_result$eps,
    DBSCAN_minPts = db_result$minPts,
    DBSCAN_Noise = sum(db_model$cluster == 0),
    DBSCAN_Noise_Percent = 100 * mean(db_model$cluster == 0),
    Ensemble_Avg_Silhouette = ensemble_avg_sil,
    SVM_CV_Accuracy = svm_result$accuracy,
    SVM_CV_Balanced_Accuracy = svm_result$balanced_accuracy,
    RF_CV_Accuracy = rf_result$accuracy,
    RF_CV_Balanced_Accuracy = rf_result$balanced_accuracy,
    NB_CV_Accuracy = nb_result$accuracy,
    NB_CV_Balanced_Accuracy = nb_result$balanced_accuracy
  )
  
  list(data = data, correlation = corr_result$matrix, kmeans = km_model,
       kmeans_table = kmeans_table, dbscan = db_model, dbscan_table = dbscan_table,
       ensemble_table = ensemble_table, svm_table = svm_table, rf_table = rf_table,
       nb_table = nb_table, feature_comparison = feature_comparison,
       clustering_comparison = clustering_comparison,
       combined_tables = combined_tables, combined_pdf_path = combined_pdf_path,
       summary = summary_row)
}

# =============================================================================
# Run both stations
# =============================================================================

WIZ_result  <- analyse_station("WIZ",  station_files["WIZ"])
WSRZ_result <- analyse_station("WSRZ", station_files["WSRZ"])

# =============================================================================
# Compare WIZ and WSRZ
# =============================================================================

cat("\n\n==================== WIZ vs WSRZ ====================\n")

comparison_summary <- bind_rows(WIZ_result$summary, WSRZ_result$summary)
print(comparison_summary)
write_csv(comparison_summary, file.path(CLUSTER_DIR, "WIZ_WSRZ_Model_Comparison.csv"))

weekly_comparison <- full_join(
  WIZ_result$data %>%
    select(Week_Start, Week_End,
           WIZ_True_Class = Model_Label, WIZ_Kmeans = Kmeans_Cluster,
           WIZ_DBSCAN = DBSCAN_Cluster, WIZ_Ensemble = Ensemble_Cluster,
           WIZ_SVM = SVM_Prediction, WIZ_RF = RF_Prediction, WIZ_NB = NB_Prediction),
  WSRZ_result$data %>%
    select(Week_Start, Week_End,
           WSRZ_True_Class = Model_Label, WSRZ_Kmeans = Kmeans_Cluster,
           WSRZ_DBSCAN = DBSCAN_Cluster, WSRZ_Ensemble = Ensemble_Cluster,
           WSRZ_SVM = SVM_Prediction, WSRZ_RF = RF_Prediction, WSRZ_NB = NB_Prediction),
  by = c("Week_Start", "Week_End")
)

print(weekly_comparison)
write_csv(weekly_comparison, file.path(CLUSTER_DIR, "WIZ_WSRZ_Weekly_Cluster_Comparison.csv"))

true_class_comparison <- bind_rows(
  WIZ_result$data  %>% count(Model_Label, name = "Number_of_Weeks") %>% mutate(Station = "WIZ"),
  WSRZ_result$data %>% count(Model_Label, name = "Number_of_Weeks") %>% mutate(Station = "WSRZ")
)

print(true_class_comparison)
write_csv(true_class_comparison, file.path(CLUSTER_DIR, "WIZ_WSRZ_True_Class_Distribution.csv"))

# ---- combine everything across BOTH stations into single files ------------

all_stations_feature_comparison <- bind_rows(
  WIZ_result$feature_comparison$results,
  WSRZ_result$feature_comparison$results
)
write_csv(all_stations_feature_comparison,
          file.path(CLUSTER_DIR, "WIZ_WSRZ_FeatureSet_Balance_Accuracy_Combined.csv"))

all_stations_clustering_comparison <- bind_rows(
  WIZ_result$clustering_comparison$metrics,
  WSRZ_result$clustering_comparison$metrics
)
write_csv(all_stations_clustering_comparison,
          file.path(CLUSTER_DIR, "WIZ_WSRZ_Clustering_FeatureSet_Comparison_Combined.csv"))

all_stations_confusion_tables <- bind_rows(
  WIZ_result$combined_tables$combined,
  WSRZ_result$combined_tables$combined
)
write_csv(all_stations_confusion_tables,
          file.path(CLUSTER_DIR, "WIZ_WSRZ_All_Confusion_Tables_Combined.csv"))

# one PDF with every plot from both stations
grand_pdf_path <- file.path(CLUSTER_DIR, "WIZ_WSRZ_All_Plots_Combined.pdf")
pdf_combine(input = c(WIZ_result$combined_pdf_path, WSRZ_result$combined_pdf_path),
            output = grand_pdf_path)

cat("\nAll results combined into:\n")
cat(" -", file.path(CLUSTER_DIR, "WIZ_WSRZ_FeatureSet_Balance_Accuracy_Combined.csv"), "\n")
cat(" -", file.path(CLUSTER_DIR, "WIZ_WSRZ_Clustering_FeatureSet_Comparison_Combined.csv"), "\n")
cat(" -", file.path(CLUSTER_DIR, "WIZ_WSRZ_All_Confusion_Tables_Combined.csv"), "\n")
cat(" -", grand_pdf_path, "\n")

cat("\nAnalysis complete.\n")