# ============================================================
# SVM CLASSIFICATION - 7 FEATURES
# Compare SVM kernels using 5-fold cross-validation
# PCA and easylabel analysis for Radial SVM
# ============================================================


# ============================================================
# 1. CLEAR ENVIRONMENT
# ============================================================

rm(list = ls())


# ============================================================
# 2. LOAD PACKAGES
# ============================================================

library(readr)
library(here)
library(dplyr)
library(e1071)
library(caret)
library(writexl)
library(ggplot2)
library(easylabel)


# ============================================================
# 3. FILE LOCATIONS
# ============================================================

DATA_DIR <- here("data")

OUT_DIR <- here(
  "results",
  "Weekly_OP_Analysis_2007to2019",
  "SVM_classification"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. LOAD WEEKLY ORDINAL-PATTERN FEATURE DATA
# ============================================================

df <- read_csv(
  here(
    "results",
    "Weekly_OP_Analysis_2007to2019",
    "weekly_entropy_complexity_rsam_avg_2007to2019.csv"
  )
)


# ============================================================
# 5. SELECT THE 7 FEATURES
# ============================================================

features_7 <- c(
  "Fisher_Info",
  "C_Fisher",
  "Var_H_Shannon",
  "Var_C_Shannon",
  "Var_H_Tsallis",
  "Var_H_Renyi",
  "Disequilibrium"
)


# Keep only rows with complete values for all 7 features
svm_data <- df[
  complete.cases(df[, features_7]),
]


# ============================================================
# 6. LOAD VAL GROUND-TRUTH DATA
# ============================================================

val_df <- read_csv(
  file.path(
    DATA_DIR,
    "VAL_WhakaariWhiteIsland.csv"
  )
)


# Convert VAL dates to UTC
val_df <- val_df %>%
  mutate(
    
    PeriodStart_UTC = as.POSIXct(
      `PeriodStart UTC`,
      format = "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    ),
    
    PeriodEnd_UTC = as.POSIXct(
      `PeriodEnd UTC`,
      format = "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    )
  )


# ============================================================
# 7. DEAL WITH OPEN-ENDED VAL PERIODS
# ============================================================

open_ended <-
  is.na(val_df$PeriodEnd_UTC) |
  format(val_df$`PeriodEnd UTC`, "%Y") == "9999"

val_df$PeriodEnd_UTC[open_ended] <- Sys.time()


# ============================================================
# 8. CREATE BINARY CLASS LABELS
#
# VAL = 1  --> Not unrest
# VAL > 1  --> Unrest
# ============================================================

val_df$true_label <- ifelse(
  val_df$VAL == 1,
  "Not unrest",
  "Unrest"
)


# ============================================================
# 9. FUNCTION TO ASSIGN VAL LABEL TO EACH WEEK
# ============================================================

assign_true_label <- function(week_start, week_end) {
  
  overlap_start <- pmax(
    val_df$PeriodStart_UTC,
    week_start
  )
  
  overlap_end <- pmin(
    val_df$PeriodEnd_UTC,
    week_end
  )
  
  overlap_seconds <- as.numeric(
    overlap_end - overlap_start,
    units = "secs"
  )
  
  overlap_seconds[overlap_seconds < 0] <- 0
  
  
  # If there is no overlap
  if (all(overlap_seconds == 0)) {
    return(NA_character_)
  }
  
  
  # Use the VAL state with the largest overlap
  val_df$true_label[
    which.max(overlap_seconds)
  ]
}


# ============================================================
# 10. ASSIGN TRUE LABELS TO ALL WEEKS
# ============================================================

svm_data$True_Label <- mapply(
  assign_true_label,
  svm_data$Week_Start,
  svm_data$Week_End
)


# Remove weeks without a VAL label
svm_data <- svm_data %>%
  filter(!is.na(True_Label))


# Convert label to factor
svm_data$True_Label <- factor(
  svm_data$True_Label,
  levels = c(
    "Not unrest",
    "Unrest"
  )
)


# Add a unique ID to every week
svm_data$Row_ID <- seq_len(
  nrow(svm_data)
)


# Check the number of observations in each class
cat("\nNumber of observations in each class:\n")

print(
  table(svm_data$True_Label)
)


# ============================================================
# 11. CREATE 5-FOLD CROSS-VALIDATION
# ============================================================

set.seed(123)

folds <- createFolds(
  svm_data$True_Label,
  k = 5,
  list = TRUE
)


# ============================================================
# 12. FUNCTION TO RUN ONE SVM KERNEL
# ============================================================

run_svm <- function(kernel_name) {
  
  all_predictions <- data.frame()
  
  
  # Run the five folds
  for (i in 1:5) {
    
    
    # --------------------------------------------------------
    # Training and test observations
    # --------------------------------------------------------
    
    test_index <- folds[[i]]
    
    train_data <- svm_data[-test_index, ]
    
    test_data <- svm_data[test_index, ]
    
    
    # --------------------------------------------------------
    # Predictor variables
    # --------------------------------------------------------
    
    train_x <- train_data[, features_7]
    
    test_x <- test_data[, features_7]
    
    
    # --------------------------------------------------------
    # Class labels
    # --------------------------------------------------------
    
    train_y <- train_data$True_Label
    
    test_y <- test_data$True_Label
    
    
    # --------------------------------------------------------
    # Scale training data
    # --------------------------------------------------------
    
    train_x <- scale(train_x)
    
    
    # Scale test data using training mean and SD
    test_x <- scale(
      test_x,
      center = attr(train_x, "scaled:center"),
      scale = attr(train_x, "scaled:scale")
    )
    
    
    # --------------------------------------------------------
    # Calculate class weights
    # --------------------------------------------------------
    
    unrest_weight <-
      sum(train_y == "Not unrest") /
      sum(train_y == "Unrest")
    
    
    weights <- c(
      "Not unrest" = 1,
      "Unrest" = unrest_weight
    )
    
    
    # --------------------------------------------------------
    # Train SVM
    # --------------------------------------------------------
    
    model <- svm(
      x = train_x,
      y = train_y,
      kernel = kernel_name,
      class.weights = weights
    )
    
    
    # --------------------------------------------------------
    # Predict test observations
    # --------------------------------------------------------
    
    prediction <- predict(
      model,
      test_x
    )
    
    
    # --------------------------------------------------------
    # Save predictions
    # --------------------------------------------------------
    
    fold_results <- data.frame(
      
      Row_ID = test_data$Row_ID,
      
      Week_Start = test_data$Week_Start,
      
      Week_End = test_data$Week_End,
      
      True_Label = test_y,
      
      Predicted_Label = prediction
    )
    
    
    all_predictions <- bind_rows(
      all_predictions,
      fold_results
    )
  }
  
  
  # ==========================================================
  # CLASSIFY EACH RESULT AS TP, TN, FP OR FN
  # ==========================================================
  
  all_predictions <- all_predictions %>%
    
    mutate(
      
      Match_Status = case_when(
        
        True_Label == "Unrest" &
          Predicted_Label == "Unrest" ~ "TP",
        
        True_Label == "Not unrest" &
          Predicted_Label == "Not unrest" ~ "TN",
        
        True_Label == "Not unrest" &
          Predicted_Label == "Unrest" ~ "FP",
        
        True_Label == "Unrest" &
          Predicted_Label == "Not unrest" ~ "FN"
      )
    )
  
  
  # ==========================================================
  # CONFUSION MATRIX
  # ==========================================================
  
  cm <- confusionMatrix(
    all_predictions$Predicted_Label,
    all_predictions$True_Label,
    positive = "Unrest"
  )
  
  
  # Count TP, TN, FP and FN
  counts <- table(
    factor(
      all_predictions$Match_Status,
      levels = c(
        "TP",
        "TN",
        "FP",
        "FN"
      )
    )
  )
  
  
  # ==========================================================
  # PERFORMANCE SUMMARY
  # ==========================================================
  
  results <- data.frame(
    
    Kernel = kernel_name,
    
    TP = as.numeric(
      counts["TP"]
    ),
    
    TN = as.numeric(
      counts["TN"]
    ),
    
    FP = as.numeric(
      counts["FP"]
    ),
    
    FN = as.numeric(
      counts["FN"]
    ),
    
    Accuracy =
      as.numeric(
        cm$overall["Accuracy"]
      ),
    
    Precision =
      as.numeric(
        cm$byClass["Pos Pred Value"]
      ),
    
    Recall =
      as.numeric(
        cm$byClass["Sensitivity"]
      ),
    
    Specificity =
      as.numeric(
        cm$byClass["Specificity"]
      ),
    
    F1 =
      as.numeric(
        cm$byClass["F1"]
      )
  )
  
  
  # Return results
  return(
    list(
      
      summary = results,
      
      predictions = all_predictions,
      
      confusion_matrix = cm
    )
  )
}


# ============================================================
# 13. RUN THE FOUR SVM KERNELS
# ============================================================

svm_linear <- run_svm(
  "linear"
)

svm_radial <- run_svm(
  "radial"
)

svm_polynomial <- run_svm(
  "polynomial"
)

svm_sigmoid <- run_svm(
  "sigmoid"
)


# ============================================================
# 14. COMPARE THE FOUR KERNELS
# ============================================================

kernel_comparison <- bind_rows(
  
  svm_linear$summary,
  
  svm_radial$summary,
  
  svm_polynomial$summary,
  
  svm_sigmoid$summary
)


cat("\nSVM Kernel Comparison:\n")

print(
  kernel_comparison
)


# ============================================================
# 15. SAVE KERNEL COMPARISON
# ============================================================

write_xlsx(
  
  kernel_comparison,
  
  file.path(
    OUT_DIR,
    "SVM_kernel_comparison_7features.xlsx"
  )
)


# ============================================================
# 16. CHECK RADIAL SVM RESULTS
# ============================================================

cat("\nRadial SVM Summary:\n")

print(
  svm_radial$summary
)


cat("\nRadial SVM TP, TN, FP and FN:\n")

print(
  table(
    svm_radial$predictions$Match_Status
  )
)


cat("\nRadial SVM Confusion Matrix:\n")

print(
  svm_radial$confusion_matrix
)


# ============================================================
# 17. PCA USING THE 7 FEATURES
# ============================================================

pca_result <- prcomp(
  scale(
    svm_data[, features_7]
  )
)


# Create PCA data frame
pca_data <- data.frame(
  
  Row_ID = svm_data$Row_ID,
  
  PC1 = pca_result$x[, 1],
  
  PC2 = pca_result$x[, 2],
  
  Week_Start = svm_data$Week_Start
)


# ============================================================
# 18. ADD RADIAL SVM RESULTS TO PCA DATA
# ============================================================

pca_radial <- pca_data %>%
  
  left_join(
    
    svm_radial$predictions %>%
      
      select(
        Row_ID,
        True_Label,
        Predicted_Label,
        Match_Status
      ),
    
    by = "Row_ID"
  )


# Check that PCA and SVM counts agree
cat("\nPCA classification counts:\n")

print(
  table(
    pca_radial$Match_Status
  )
)


# ============================================================
# 19. PCA DIAGNOSTIC PLOT
# ============================================================

p_radial <- ggplot(
  pca_radial,
  aes(
    x = PC1,
    y = PC2,
    color = Match_Status
  )
) +
  
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  
  labs(
    title = "Radial SVM - 7 Features",
    x = "PC1",
    y = "PC2",
    color = "Classification"
  ) +
  
  theme_bw()


print(
  p_radial
)


# ============================================================
# 20. EASYLABEL - TRUE POSITIVES
# ============================================================

easylabel(
  
  pca_radial %>%
    filter(
      Match_Status == "TP"
    ),
  
  x = "PC1",
  
  y = "PC2",
  
  col = "Match_Status",
  
  labs = "Week_Start",
  
  size = 8
)


# ============================================================
# 21. EASYLABEL - FALSE NEGATIVES
# ============================================================

easylabel(
  
  pca_radial %>%
    filter(
      Match_Status == "FN"
    ),
  
  x = "PC1",
  
  y = "PC2",
  
  col = "Match_Status",
  
  labs = "Week_Start",
  
  size = 8
)


# ============================================================
# 22. EASYLABEL - FALSE POSITIVES
# ============================================================

easylabel(
  
  pca_radial %>%
    filter(
      Match_Status == "FP"
    ),
  
  x = "PC1",
  
  y = "PC2",
  
  col = "Match_Status",
  
  labs = "Week_Start",
  
  size = 8
)
