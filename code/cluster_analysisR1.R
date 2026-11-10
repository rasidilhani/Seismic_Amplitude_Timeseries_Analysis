# =============================================================================
# WIZ and WSRZ
# 13-feature correlation -> 7 selected features ->
# K-means -> DBSCAN -> SVM -> Station comparison
# =============================================================================

library(tidyverse)
library(corrplot)
library(dbscan)
library(e1071)
library(caret)

set.seed(123456789, kind = "Mersenne-Twister")


# =============================================================================
# 1. Paths
# =============================================================================

DATA_DIR <- paste0(
  "C:/Users/UserA1/Documents/GitHub/",
  "Seismic_Amplitude_Timeseries_Analysis/data"
)

RESULTS_DIR <- paste0(
  "C:/Users/UserA1/Documents/GitHub/",
  "Seismic_Amplitude_Timeseries_Analysis/results"
)

CLUSTER_DIR <- file.path(
  RESULTS_DIR,
  "cluster_analysis_R1"
)

dir.create(
  CLUSTER_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


VAL_FILE <- file.path(
  DATA_DIR,
  "VAL_WhakaariWhiteIsland.csv"
)


station_files <- c(
  
  WIZ = file.path(
    DATA_DIR,
    "weekly_HC_RSAM_WIZ_20140701_to_20191231.csv"
  ),
  
  WSRZ = file.path(
    DATA_DIR,
    "weekly_HC_RSAM_WSRZ_20140701_to_20191231.csv"
  )
)


# =============================================================================
# 2. Analysis period
# =============================================================================

ANALYSIS_START <- as.POSIXct(
  "2014-07-01 00:00:00",
  tz = "Pacific/Auckland"
)

ANALYSIS_END <- as.POSIXct(
  "2019-12-31 23:59:59",
  tz = "Pacific/Auckland"
)

VAL_SYSTEM_VERSION <- 3


# =============================================================================
# 3. ORIGINAL 13 FEATURES
# =============================================================================
#
# These are used FIRST for the correlation analysis.
#
# =============================================================================

all_features <- c(
  
  "H_Shannon",
  "C_Shannon",
  
  "Fisher_Info",
  "C_Fisher",
  
  "H_Renyi",
  "C_Renyi",
  
  "H_Tsallis",
  "C_Tsallis",
  
  "Var_H_Shannon",
  "Var_C_Shannon",
  
  "Var_H_Renyi",
  "Var_H_Tsallis",
  
  "Disequilibrium"
)


cat("\nOriginal feature set:\n")
cat("Number of features =", length(all_features), "\n\n")

print(all_features)


# =============================================================================
# 4. SELECTED 7 FEATURES
# =============================================================================
#
# These seven features are used for:
#
#   K-means
#   DBSCAN
#   SVM
#
# =============================================================================

feature_cols <- c(
  
  "Fisher_Info",
  "C_Fisher",
  
  "Var_H_Shannon",
  "Var_C_Shannon",
  
  "Var_H_Tsallis",
  "Var_H_Renyi",
  
  "Disequilibrium"
)


cat("\nSelected features for cluster/classification analysis:\n")
cat("Number of selected features =", length(feature_cols), "\n\n")

print(feature_cols)


# =============================================================================
# 5. Three VAL activity classes
# =============================================================================

activity_levels <- c(
  
  "Background unrest",
  
  "Moderate to heightened unrest",
  
  "Minor to moderate eruption"
  
)


model_label_levels <- c(
  
  "1_Background_unrest",
  
  "2_Moderate_to_heightened_unrest",
  
  "3_Minor_to_moderate_eruption"
  
)


# =============================================================================
# 6. Read VAL data
# =============================================================================

val_df <- read_csv(
  VAL_FILE,
  show_col_types = FALSE
) %>%
  
  mutate(
    
    PeriodStart_Local = as.POSIXct(
      
      stringr::str_remove(
        `PeriodStart LocalTime`,
        " NZ[SD]T$"
      ),
      
      format = "%Y-%m-%d %H:%M:%S",
      tz = "Pacific/Auckland"
      
    ),
    
    PeriodEnd_Local = as.POSIXct(
      
      stringr::str_remove(
        `PeriodEnd LocalTime`,
        " NZ[SD]T$"
      ),
      
      format = "%Y-%m-%d %H:%M:%S",
      tz = "Pacific/Auckland"
      
    )
    
  )


# -----------------------------------------------------------------------------
# Replace open-ended VAL periods
# -----------------------------------------------------------------------------

val_df$PeriodEnd_Local[
  
  is.na(val_df$PeriodEnd_Local) |
    
    stringr::str_starts(
      val_df$`PeriodEnd LocalTime`,
      "9999"
    )
  
] <- ANALYSIS_END


# -----------------------------------------------------------------------------
# Keep only VAL System Version 3
# -----------------------------------------------------------------------------

val_df <- val_df %>%
  
  filter(
    
    `VAL System Version` == VAL_SYSTEM_VERSION,
    
    PeriodEnd_Local > ANALYSIS_START,
    
    PeriodStart_Local < ANALYSIS_END
    
  )


# =============================================================================
# 7. Convert VAL descriptions into THREE activity classes
# =============================================================================

val_df <- val_df %>%
  
  mutate(
    
    Model_Label_Num = case_when(
      
      Description ==
        "Minor volcanic unrest" ~ 1L,
      
      Description ==
        "Moderate to heightened volcanic unrest" ~ 2L,
      
      Description %in% c(
        "Minor volcanic eruption",
        "Moderate volcanic eruption"
      ) ~ 3L,
      
      TRUE ~ NA_integer_
      
    ),
    
    Activity_3Class = case_when(
      
      Model_Label_Num == 1 ~
        "Background unrest",
      
      Model_Label_Num == 2 ~
        "Moderate to heightened unrest",
      
      Model_Label_Num == 3 ~
        "Minor to moderate eruption",
      
      TRUE ~ NA_character_
      
    ),
    
    Model_Label = case_when(
      
      Model_Label_Num == 1 ~
        "1_Background_unrest",
      
      Model_Label_Num == 2 ~
        "2_Moderate_to_heightened_unrest",
      
      Model_Label_Num == 3 ~
        "3_Minor_to_moderate_eruption",
      
      TRUE ~ NA_character_
      
    )
    
  ) %>%
  
  filter(
    !is.na(Model_Label_Num)
  )


# =============================================================================
# 8. Function for assigning VAL class to a weekly window
# =============================================================================
#
# If a week overlaps more than one VAL period,
# use the highest activity class.
#
# =============================================================================

assign_val_class <- function(
    week_start,
    week_end
) {
  
  overlap <- val_df %>%
    
    filter(
      
      PeriodStart_Local < week_end,
      
      PeriodEnd_Local > week_start
      
    )
  
  
  if (nrow(overlap) == 0) {
    
    return(NA_integer_)
    
  }
  
  
  max(
    overlap$Model_Label_Num,
    na.rm = TRUE
  )
  
}


# =============================================================================
# 9. Function for 13-feature correlogram
# =============================================================================

make_correlogram <- function(
    data,
    station_name
) {
  
  # ---------------------------------------------------------------------------
  # Check that all 13 features exist
  # ---------------------------------------------------------------------------
  
  missing_features <- setdiff(
    all_features,
    names(data)
  )
  
  
  if (length(missing_features) > 0) {
    
    stop(
      paste(
        station_name,
        "is missing these features:",
        paste(missing_features, collapse = ", ")
      )
    )
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Spearman correlation
  # ---------------------------------------------------------------------------
  
  cor_matrix <- cor(
    
    data[, all_features],
    
    use = "pairwise.complete.obs",
    
    method = "spearman"
    
  )
  
  
  # ---------------------------------------------------------------------------
  # Display correlogram
  # ---------------------------------------------------------------------------
  
  corrplot(
    
    cor_matrix,
    
    method = "color",
    
    type = "upper",
    
    order = "hclust",
    
    tl.col = "black",
    
    tl.srt = 45,
    
    tl.cex = 0.75,
    
    addCoef.col = "black",
    
    number.cex = 0.45,
    
    diag = FALSE,
    
    title = paste(
      station_name,
      "- Correlation of 13 Original Features"
    ),
    
    mar = c(0, 0, 2, 0)
    
  )
  
  
  # ---------------------------------------------------------------------------
  # Save correlogram
  # ---------------------------------------------------------------------------
  
  pdf(
    
    file.path(
      CLUSTER_DIR,
      paste0(
        station_name,
        "_13_Feature_Correlogram.pdf"
      )
    ),
    
    width = 11,
    height = 9
    
  )
  
  
  corrplot(
    
    cor_matrix,
    
    method = "color",
    
    type = "upper",
    
    order = "hclust",
    
    tl.col = "black",
    
    tl.srt = 45,
    
    tl.cex = 0.75,
    
    addCoef.col = "black",
    
    number.cex = 0.45,
    
    diag = FALSE,
    
    title = paste(
      station_name,
      "- Correlation of 13 Original Features"
    ),
    
    mar = c(0, 0, 2, 0)
    
  )
  
  
  dev.off()
  
  
  # Save numerical correlation matrix
  
  write.csv(
    
    cor_matrix,
    
    file.path(
      CLUSTER_DIR,
      paste0(
        station_name,
        "_13_Feature_Correlation_Matrix.csv"
      )
    )
    
  )
  
  
  return(cor_matrix)
  
}


# =============================================================================
# 10. Simple DBSCAN search
# =============================================================================
#
# DBSCAN does NOT naturally force exactly 3 clusters.
#
# Therefore:
#   search eps and minPts
#   keep solutions giving exactly 3 non-noise clusters
#   choose the solution with the smallest amount of noise
#
# Cluster 0 = noise.
#
# =============================================================================

find_three_cluster_dbscan <- function(X) {
  
  search_grid <- expand_grid(
    
    eps = seq(
      0.05,
      5.00,
      by = 0.05
    ),
    
    minPts = 3:20
    
  )
  
  
  search_grid$n_clusters <- NA_integer_
  search_grid$noise <- NA_integer_
  
  
  for (i in seq_len(nrow(search_grid))) {
    
    fit <- dbscan::dbscan(
      
      X,
      
      eps = search_grid$eps[i],
      
      minPts = search_grid$minPts[i]
      
    )
    
    
    clusters <- fit$cluster
    
    
    search_grid$n_clusters[i] <-
      
      length(
        setdiff(
          unique(clusters),
          0
        )
      )
    
    
    search_grid$noise[i] <-
      
      sum(
        clusters == 0
      )
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Keep solutions with exactly 3 clusters
  # ---------------------------------------------------------------------------
  
  possible <- search_grid %>%
    
    filter(
      n_clusters == 3
    )
  
  
  if (nrow(possible) == 0) {
    
    stop(
      "No DBSCAN solution with exactly 3 clusters was found."
    )
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Choose solution with least noise
  # ---------------------------------------------------------------------------
  
  best <- possible %>%
    
    arrange(
      noise,
      abs(eps - 1.1),
      minPts
    ) %>%
    
    slice(1)
  
  
  final_fit <- dbscan::dbscan(
    
    X,
    
    eps = best$eps,
    
    minPts = best$minPts
    
  )
  
  
  list(
    
    model = final_fit,
    
    eps = best$eps,
    
    minPts = best$minPts,
    
    search = search_grid
    
  )
  
}


# =============================================================================
# 11. Simple 5-fold SVM cross-validation
# =============================================================================
#
# K-means and DBSCAN are UNSUPERVISED.
#
# SVM is SUPERVISED because it uses the true VAL classes.
#
# =============================================================================

run_svm_cv <- function(
    X,
    y
) {
  
  y <- droplevels(y)
  
  
  # ---------------------------------------------------------------------------
  # Number of folds
  # ---------------------------------------------------------------------------
  
  smallest_class <- min(
    table(y)
  )
  
  
  k <- min(
    5,
    smallest_class
  )
  
  
  if (k < 2) {
    
    stop(
      "Not enough observations for SVM cross-validation."
    )
    
  }
  
  
  folds <- caret::createFolds(
    
    y,
    
    k = k,
    
    list = TRUE,
    
    returnTrain = FALSE
    
  )
  
  
  svm_prediction <- factor(
    
    rep(
      NA_character_,
      length(y)
    ),
    
    levels = levels(y)
    
  )
  
  
  # ---------------------------------------------------------------------------
  # Cross-validation
  # ---------------------------------------------------------------------------
  
  for (i in seq_along(folds)) {
    
    test_index <- folds[[i]]
    
    train_index <- setdiff(
      seq_along(y),
      test_index
    )
    
    
    y_train <- droplevels(
      y[train_index]
    )
    
    
    # -------------------------------------------------------------------------
    # Class weights
    # -------------------------------------------------------------------------
    
    class_count <- table(
      y_train
    )
    
    
    class_weights <-
      
      sum(class_count) /
      
      (
        length(class_count) *
          class_count
      )
    
    
    class_weights <- as.numeric(
      class_weights
    )
    
    
    names(class_weights) <-
      names(class_count)
    
    
    # -------------------------------------------------------------------------
    # SVM
    # -------------------------------------------------------------------------
    
    svm_model <- svm(
      
      x = X[
        train_index,
        ,
        drop = FALSE
      ],
      
      y = y_train,
      
      kernel = "radial",
      
      cost = 1,
      
      gamma = 1 / ncol(X),
      
      class.weights = class_weights,
      
      scale = FALSE
      
    )
    
    
    svm_prediction[test_index] <-
      
      as.character(
        
        predict(
          
          svm_model,
          
          X[
            test_index,
            ,
            drop = FALSE
          ]
          
        )
        
      )
    
  }
  
  
  accuracy <- mean(
    svm_prediction == y
  )
  
  
  list(
    
    prediction = svm_prediction,
    
    accuracy = accuracy,
    
    folds = k
    
  )
  
}


# =============================================================================
# 12. Function for analysing ONE station
# =============================================================================

analyse_station <- function(
    station_name,
    file_name
) {
  
  cat("\n\n")
  cat("============================================================\n")
  cat("STATION:", station_name, "\n")
  cat("============================================================\n")
  
  
  # ===========================================================================
  # A. Read weekly data
  # ===========================================================================
  
  data <- read_csv(
    file_name,
    show_col_types = FALSE
  ) %>%
    
    mutate(
      
      Week_Start = as.POSIXct(
        Week_Start,
        tz = "Pacific/Auckland"
      ),
      
      Week_End = as.POSIXct(
        Week_End,
        tz = "Pacific/Auckland"
      )
      
    ) %>%
    
    filter(
      
      Week_Start >= ANALYSIS_START,
      
      Week_Start <= ANALYSIS_END
      
    )
  
  
  cat(
    "\nNumber of weekly observations:",
    nrow(data),
    "\n"
  )
  
  
  # ===========================================================================
  # B. CORRELOGRAM USING ALL 13 FEATURES
  # ===========================================================================
  
  cat(
    "\nCreating 13-feature correlogram...\n"
  )
  
  
  correlation_matrix <- make_correlogram(
    
    data,
    
    station_name
    
  )
  
  
  # ===========================================================================
  # C. Assign true VAL classes
  # ===========================================================================
  
  data$Model_Label_Num <- mapply(
    
    assign_val_class,
    
    data$Week_Start,
    
    data$Week_End
    
  )
  
  
  data <- data %>%
    
    mutate(
      
      Activity_3Class = case_when(
        
        Model_Label_Num == 1 ~
          "Background unrest",
        
        Model_Label_Num == 2 ~
          "Moderate to heightened unrest",
        
        Model_Label_Num == 3 ~
          "Minor to moderate eruption",
        
        TRUE ~ NA_character_
        
      ),
      
      Model_Label = case_when(
        
        Model_Label_Num == 1 ~
          "1_Background_unrest",
        
        Model_Label_Num == 2 ~
          "2_Moderate_to_heightened_unrest",
        
        Model_Label_Num == 3 ~
          "3_Minor_to_moderate_eruption",
        
        TRUE ~ NA_character_
        
      ),
      
      Model_Label = factor(
        
        Model_Label,
        
        levels = model_label_levels
        
      )
      
    ) %>%
    
    filter(
      
      !is.na(Model_Label),
      
      if_all(
        all_of(feature_cols),
        ~ is.finite(.x)
      )
      
    )
  
  
  cat("\nTrue VAL classes:\n")
  
  print(
    table(data$Model_Label)
  )
  
  
  # ===========================================================================
  # D. Select the 7 features
  # ===========================================================================
  
  cluster_data <- data %>%
    
    select(
      all_of(feature_cols)
    )
  
  
  # ===========================================================================
  # E. Scale
  # ===========================================================================
  
  X <- scale(
    cluster_data
  )
  
  
  # ===========================================================================
  # F. K-MEANS
  # ===========================================================================
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("K-MEANS\n")
  cat("------------------------------------------------------------\n")
  
  
  set.seed(123)
  
  
  km_model <- kmeans(
    
    X,
    
    centers = 3,
    
    nstart = 50
    
  )
  
  
  data$Kmeans_Cluster <-
    
    factor(
      km_model$cluster
    )
  
  
  cat("\nNumber of observations in each K-means cluster:\n")
  
  print(
    table(
      data$Kmeans_Cluster
    )
  )
  
  
  cat("\nK-means versus true VAL class:\n")
  
  kmeans_table <- table(
    
    True_Class =
      data$Model_Label,
    
    Kmeans_Cluster =
      data$Kmeans_Cluster
    
  )
  
  
  print(
    kmeans_table
  )
  
  
  cat("\nK-means centres:\n")
  
  print(
    km_model$centers
  )
  
  
  # ===========================================================================
  # G. DBSCAN
  # ===========================================================================
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("DBSCAN\n")
  cat("------------------------------------------------------------\n")
  
  
  db_result <- find_three_cluster_dbscan(
    X
  )
  
  
  db_model <- db_result$model
  
  
  cat(
    "\nSelected eps =",
    db_result$eps,
    "\n"
  )
  
  
  cat(
    "Selected minPts =",
    db_result$minPts,
    "\n"
  )
  
  
  data$DBSCAN_Cluster <- factor(
    
    db_model$cluster,
    
    levels = 0:3,
    
    labels = c(
      "Noise",
      "Cluster 1",
      "Cluster 2",
      "Cluster 3"
    )
    
  )
  
  
  cat("\nDBSCAN cluster sizes:\n")
  
  print(
    table(
      data$DBSCAN_Cluster
    )
  )
  
  
  cat("\nDBSCAN versus true VAL class:\n")
  
  dbscan_table <- table(
    
    True_Class =
      data$Model_Label,
    
    DBSCAN_Cluster =
      data$DBSCAN_Cluster
    
  )
  
  
  print(
    dbscan_table
  )
  
  
  # ===========================================================================
  # H. SVM
  # ===========================================================================
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("SVM - THREE CLASS\n")
  cat("------------------------------------------------------------\n")
  
  
  svm_result <- run_svm_cv(
    
    X,
    
    data$Model_Label
    
  )
  
  
  data$SVM_Prediction <-
    svm_result$prediction
  
  
  cat(
    "\nNumber of CV folds:",
    svm_result$folds,
    "\n"
  )
  
  
  cat(
    "\nSVM cross-validated accuracy:",
    round(
      svm_result$accuracy,
      4
    ),
    "\n"
  )
  
  
  cat("\nSVM confusion matrix:\n")
  
  
  svm_table <- table(
    
    True_Class =
      data$Model_Label,
    
    SVM_Prediction =
      data$SVM_Prediction
    
  )
  
  
  print(
    svm_table
  )
  
  
  # ===========================================================================
  # I. Save result tables
  # ===========================================================================
  
  write.csv(
    
    as.data.frame.matrix(
      kmeans_table
    ),
    
    file.path(
      CLUSTER_DIR,
      paste0(
        station_name,
        "_Kmeans_vs_VAL.csv"
      )
    )
    
  )
  
  
  write.csv(
    
    as.data.frame.matrix(
      dbscan_table
    ),
    
    file.path(
      CLUSTER_DIR,
      paste0(
        station_name,
        "_DBSCAN_vs_VAL.csv"
      )
    )
    
  )
  
  
  write.csv(
    
    as.data.frame.matrix(
      svm_table
    ),
    
    file.path(
      CLUSTER_DIR,
      paste0(
        station_name,
        "_SVM_vs_VAL.csv"
      )
    )
    
  )
  
  
  write_csv(
    
    data,
    
    file.path(
      CLUSTER_DIR,
      paste0(
        station_name,
        "_Weekly_Cluster_Results.csv"
      )
    )
    
  )
  
  
  # ===========================================================================
  # J. Simple station summary
  # ===========================================================================
  
  station_summary <- tibble(
    
    Station =
      station_name,
    
    Number_of_Weeks =
      nrow(data),
    
    Kmeans_BetweenSS_Ratio =
      km_model$betweenss /
      km_model$totss,
    
    DBSCAN_eps =
      db_result$eps,
    
    DBSCAN_minPts =
      db_result$minPts,
    
    DBSCAN_Noise =
      sum(
        db_model$cluster == 0
      ),
    
    DBSCAN_Noise_Percent =
      100 *
      mean(
        db_model$cluster == 0
      ),
    
    SVM_CV_Accuracy =
      svm_result$accuracy
    
  )
  
  
  return(
    
    list(
      
      data = data,
      
      correlation =
        correlation_matrix,
      
      kmeans =
        km_model,
      
      kmeans_table =
        kmeans_table,
      
      dbscan =
        db_model,
      
      dbscan_table =
        dbscan_table,
      
      svm_table =
        svm_table,
      
      summary =
        station_summary
      
    )
    
  )
  
}


# =============================================================================
# 13. RUN WIZ
# =============================================================================

WIZ_result <- analyse_station(
  
  station_name = "WIZ",
  
  file_name =
    station_files["WIZ"]
  
)


# =============================================================================
# 14. RUN WSRZ
# =============================================================================

WSRZ_result <- analyse_station(
  
  station_name = "WSRZ",
  
  file_name =
    station_files["WSRZ"]
  
)


# =============================================================================
# 15. COMPARE WIZ AND WSRZ
# =============================================================================

cat("\n\n")
cat("============================================================\n")
cat("WIZ AND WSRZ COMPARISON\n")
cat("============================================================\n")


comparison_summary <- bind_rows(
  
  WIZ_result$summary,
  
  WSRZ_result$summary
  
)


print(
  comparison_summary
)


write_csv(
  
  comparison_summary,
  
  file.path(
    CLUSTER_DIR,
    "WIZ_WSRZ_Model_Comparison.csv"
  )
  
)


# =============================================================================
# 16. Compare weekly classifications side-by-side
# =============================================================================

WIZ_compare <- WIZ_result$data %>%
  
  select(
    
    Week_Start,
    Week_End,
    
    WIZ_True_Class =
      Model_Label,
    
    WIZ_Kmeans =
      Kmeans_Cluster,
    
    WIZ_DBSCAN =
      DBSCAN_Cluster,
    
    WIZ_SVM =
      SVM_Prediction
    
  )


WSRZ_compare <- WSRZ_result$data %>%
  
  select(
    
    Week_Start,
    Week_End,
    
    WSRZ_True_Class =
      Model_Label,
    
    WSRZ_Kmeans =
      Kmeans_Cluster,
    
    WSRZ_DBSCAN =
      DBSCAN_Cluster,
    
    WSRZ_SVM =
      SVM_Prediction
    
  )


weekly_comparison <- full_join(
  
  WIZ_compare,
  
  WSRZ_compare,
  
  by = c(
    "Week_Start",
    "Week_End"
  )
  
)


print(
  weekly_comparison
)


write_csv(
  
  weekly_comparison,
  
  file.path(
    CLUSTER_DIR,
    "WIZ_WSRZ_Weekly_Cluster_Comparison.csv"
  )
  
)


# =============================================================================
# 17. Compare true class distributions
# =============================================================================

true_class_comparison <- bind_rows(
  
  WIZ_result$data %>%
    
    count(
      Model_Label,
      name = "Number_of_Weeks"
    ) %>%
    
    mutate(
      Station = "WIZ"
    ),
  
  
  WSRZ_result$data %>%
    
    count(
      Model_Label,
      name = "Number_of_Weeks"
    ) %>%
    
    mutate(
      Station = "WSRZ"
    )
  
)


print(
  true_class_comparison
)


write_csv(
  
  true_class_comparison,
  
  file.path(
    CLUSTER_DIR,
    "WIZ_WSRZ_True_Class_Distribution.csv"
  )
  
)


cat("\nAnalysis complete.\n")

