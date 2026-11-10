# Weekly ordinal-pattern analysis validated against Whakaari VAL records
# Analysis period: VAL system version 3, 2014-07-01 onward
#
# CHANGE FROM ORIGINAL: Sections 5-7 replaced a single stratified 70/30
# train/test split with repeated stratified 5-fold CV (5 folds x 20 repeats
# for multinomial logistic regression; 5 folds x 10 repeats for the SVM grid search, since it searches
# 40 cost/gamma combinations and runtime scales with repeats x combinations
# -- raise SVM_N_REPEATS back to 20 if runtime allows). This avoids the
# single-split instability problem: with a rare eruption class, one random
# 30% test split can land very few (or unrepresentatively many) eruption
# cases, making a single reported confusion matrix unreliable.
#
# The VAL labeling logic in Section 3 is restricted to VAL System Version 3.
# Model labels are:
#   1 = background unrest
#   2 = moderate to heightened volcanic unrest
#   3 = minor to moderate volcanic eruption

library(tidyverse)
library(here)
library(dbscan)
library(e1071)
library(caret)
library(ggalluvial)
library(writexl)
library(easylabel)
library(nnet)

set.seed(123456789, kind = "Mersenne-Twister")

# -----------------------------------------------------------------------------
# 1. Paths and settings
# -----------------------------------------------------------------------------
DATA_DIR <- here("data")
RESULTS_DIR <- here("results", "Weekly_OP_Analysis_2007to2019")
OUT_DIR <- file.path(RESULTS_DIR, "Cluster Analysis_System3_2014to2019")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OP_FILE <- file.path(
  RESULTS_DIR,
  "weekly_entropy_complexity_rsam_avg_2007to2019.csv"
)
VAL_FILE <- file.path(DATA_DIR, "VAL_WhakaariWhiteIsland.csv")

ANALYSIS_START_UTC <- as.POSIXct("2014-07-01 00:00:00", tz = "UTC")
VAL_SYSTEM_VERSION <- 3
FEATURE_CORRELATION_CUTOFF <- 0.90

features <- c(
  "Fisher_Info", "C_Fisher",
  "Var_H_Shannon", "Var_C_Shannon",
  "Var_H_Tsallis", "Var_H_Renyi",
  "Disequilibrium"
)

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

activity_to_model_level <- function(activity) {
  case_when(
    activity == "Background unrest" ~ 1L,
    activity == "Moderate to heightened unrest" ~ 2L,
    activity == "Minor to moderate eruption" ~ 3L,
    TRUE ~ NA_integer_
  )
}

model_level_to_activity <- function(level) {
  case_when(
    level == 1L ~ "Background unrest",
    level == 2L ~ "Moderate to heightened unrest",
    level == 3L ~ "Minor to moderate eruption",
    TRUE ~ NA_character_
  )
}

model_level_to_label <- function(level) {
  case_when(
    level == 1L ~ "1_Background_unrest",
    level == 2L ~ "2_Moderate_to_heightened_unrest",
    level == 3L ~ "3_Minor_to_moderate_eruption",
    TRUE ~ NA_character_
  )
}

EPS_VALUE <- 1.1
MINPTS_MULTIPLIER <- 2

# Repeated stratified CV settings (replaces the single 70/30 split)
N_FOLDS <- 5
MULTINOM_N_REPEATS <- 20
SVM_N_REPEATS <- 10  # lower because the cost/gamma grid is large;
# raise to 20 if runtime allows

# Set to TRUE when running interactively in RStudio. Each easylabel window
# opens in turn; close the current window to continue to the next PCA plot.
RUN_EASYLABEL <- interactive()

serif_theme <- theme_bw(base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# -----------------------------------------------------------------------------
# 2. Read weekly ordinal-pattern features
# -----------------------------------------------------------------------------
weekly_df <- read_csv(OP_FILE, show_col_types = FALSE) %>%
  mutate(
    Week_Start = as.POSIXct(Week_Start, tz = "UTC"),
    Week_End = as.POSIXct(Week_End, tz = "UTC")
  ) %>%
  filter(Week_Start >= ANALYSIS_START_UTC) %>%
  filter(if_all(all_of(features), is.finite))

if (anyNA(weekly_df$Week_Start) || anyNA(weekly_df$Week_End)) {
  stop("Week_Start or Week_End could not be converted to UTC date-times.")
}

# -----------------------------------------------------------------------------
# 3. Read VAL system 3 records and define model labels
# -----------------------------------------------------------------------------
val_df <- read_csv(VAL_FILE, show_col_types = FALSE) %>%
  mutate(
    PeriodStart_UTC = as.POSIXct(
      `PeriodStart UTC`, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"
    ),
    PeriodEnd_UTC = as.POSIXct(
      `PeriodEnd UTC`, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"
    )
  )

val_df$PeriodEnd_UTC[
  is.na(val_df$PeriodEnd_UTC) |
    stringr::str_starts(val_df$`PeriodEnd UTC`, "9999")
] <- max(weekly_df$Week_End)

val_df <- val_df %>%
  filter(
    `VAL System Version` == VAL_SYSTEM_VERSION,
    PeriodEnd_UTC > min(weekly_df$Week_Start),
    PeriodStart_UTC < max(weekly_df$Week_End)
  ) %>%
  mutate(
    Activity_3Class = case_when(
      Description == "Minor volcanic unrest" ~ "Background unrest",
      Description == "Moderate to heightened volcanic unrest" ~
        "Moderate to heightened unrest",
      Description %in% c("Minor volcanic eruption", "Moderate volcanic eruption") ~
        "Minor to moderate eruption",
      TRUE ~ NA_character_
    ),
    Model_Label_Num = activity_to_model_level(Activity_3Class),
    Model_Label = factor(
      model_level_to_label(Model_Label_Num),
      levels = model_label_levels
    )
  )

# -----------------------------------------------------------------------------
# 3a. Descriptive Sankey plots for the VAL system 3 record
# -----------------------------------------------------------------------------
val_colors <- c(
  "VAL 1" = "#4DAF4A",
  "VAL 2" = "#FFB000",
  "VAL 3" = "#F07F2F",
  "VAL 4" = "#C43C55"
)

activity_colors <- c(
  "Background unrest" = "#4DAF4A",
  "Moderate to heightened unrest" = "#FFB000",
  "Minor to moderate eruption" = "#D73027"
)

model_label_colors <- c(
  "1_Background_unrest" = "#4DAF4A",
  "2_Moderate_to_heightened_unrest" = "#FFB000",
  "3_Minor_to_moderate_eruption" = "#D73027"
)

# Plot 1: numerical VAL -> system version -> official description.
val_system_sankey <- val_df %>%
  transmute(
    VAL_Level = paste("VAL", VAL),
    VAL_System = paste("System version", `VAL System Version`),
    Description = stringr::str_wrap(Description, width = 25)
  ) %>%
  count(VAL_Level, VAL_System, Description, name = "Number_of_periods")

p_val_system_sankey <- ggplot(
  val_system_sankey,
  aes(
    axis1 = VAL_Level,
    axis2 = VAL_System,
    axis3 = Description,
    y = Number_of_periods
  )
) +
  geom_alluvium(aes(fill = VAL_Level), width = 0.15, alpha = 0.78) +
  geom_stratum(width = 0.15, fill = "grey96", color = "grey35") +
  geom_text(
    stat = "stratum", aes(label = after_stat(stratum)),
    family = "serif", size = 2.7, lineheight = 0.9
  ) +
  scale_x_discrete(
    limits = c("Volcanic Alert Level", "VAL system", "Activity description"),
    expand = c(0.12, 0.12)
  ) +
  scale_fill_manual(values = val_colors) +
  labs(
    title = "VAL system 3 activity descriptions",
    x = NULL, y = "Number of recorded VAL periods", fill = NULL
  ) +
  serif_theme

print(p_val_system_sankey)
ggsave(
  file.path(OUT_DIR, "VAL_System_Description_Sankey.pdf"),
  p_val_system_sankey, width = 10, height = 5.5
)

# Plot 2: numerical VAL -> the three requested model classes.
val_activity_sankey <- val_df %>%
  transmute(
    VAL_Level = paste("VAL", VAL),
    Activity_3Class
  ) %>%
  filter(!is.na(Activity_3Class)) %>%
  count(VAL_Level, Activity_3Class, name = "Number_of_periods")

p_val_activity_sankey <- ggplot(
  val_activity_sankey,
  aes(axis1 = VAL_Level, axis2 = Activity_3Class, y = Number_of_periods)
) +
  geom_alluvium(aes(fill = VAL_Level), width = 0.18, alpha = 0.78) +
  geom_stratum(width = 0.18, fill = "grey96", color = "grey35") +
  geom_text(
    stat = "stratum", aes(label = after_stat(stratum)),
    family = "serif", size = 3
  ) +
  scale_x_discrete(
    limits = c("Volcanic Alert Level", "Activity class"),
    expand = c(0.15, 0.15)
  ) +
  scale_fill_manual(values = val_colors) +
  labs(
    title = "VAL system 3 model label classification",
    x = NULL, y = "Number of recorded VAL periods", fill = NULL
  ) +
  serif_theme

print(p_val_activity_sankey)
ggsave(
  file.path(OUT_DIR, "VAL_Activity_Class_Sankey.pdf"),
  p_val_activity_sankey, width = 8, height = 5
)

assign_val_label <- function(week_start, week_end, val_table) {
  overlap <- val_table$PeriodStart_UTC < week_end &
    val_table$PeriodEnd_UTC > week_start
  
  if (!any(overlap, na.rm = TRUE)) {
    return(c(
      Activity_3Class = NA_character_,
      Model_Label_Num = NA_character_,
      Model_Label = NA_character_,
      VAL_Levels = NA_character_,
      VAL_Descriptions = NA_character_
    ))
  }
  
  descriptions <- unique(val_table$Description[overlap])
  val_levels <- unique(val_table$VAL[overlap])
  model_levels <- unique(na.omit(val_table$Model_Label_Num[overlap]))
  model_level <- if (length(model_levels) > 0) max(model_levels) else NA_integer_
  
  c(
    Activity_3Class = model_level_to_activity(model_level),
    Model_Label_Num = as.character(model_level),
    Model_Label = model_level_to_label(model_level),
    VAL_Levels = paste(val_levels, collapse = "; "),
    VAL_Descriptions = paste(descriptions, collapse = "; ")
  )
}

assigned <- t(mapply(
  assign_val_label,
  weekly_df$Week_Start,
  weekly_df$Week_End,
  MoreArgs = list(val_table = val_df)
))

weekly_df$Activity_3Class <- assigned[, "Activity_3Class"]
weekly_df$Model_Label_Num <- as.integer(assigned[, "Model_Label_Num"])
weekly_df$Model_Label <- assigned[, "Model_Label"]
weekly_df$VAL_Levels <- assigned[, "VAL_Levels"]
weekly_df$VAL_Descriptions <- assigned[, "VAL_Descriptions"]

weekly_df <- weekly_df %>%
  mutate(
    Activity_3Class = factor(
      Activity_3Class,
      levels = activity_levels
    ),
    Model_Label = factor(
      Model_Label,
      levels = model_label_levels
    )
  ) %>%
  filter(!is.na(Activity_3Class), !is.na(Model_Label), !is.na(Model_Label_Num))

cat("VAL system 3 model labels by weekly window:\n")
print(table(weekly_df$Model_Label, useNA = "ifany"))
cat("\nActivity labels by weekly window:\n")
print(table(weekly_df$Activity_3Class, useNA = "ifany"))
cat("\nWeekly data range used for modelling:\n")
print(range(weekly_df$Week_Start))

feature_cor <- cor(weekly_df[, features], use = "pairwise.complete.obs")
feature_cor_pairs <- purrr::map_dfr(
  combn(features, 2, simplify = FALSE),
  ~ tibble(
    Feature_1 = .x[1],
    Feature_2 = .x[2],
    Correlation = unname(feature_cor[.x[1], .x[2]]),
    Abs_Correlation = abs(Correlation)
  )
) %>%
  arrange(desc(Abs_Correlation))

high_feature_correlations <- feature_cor_pairs %>%
  filter(Abs_Correlation >= FEATURE_CORRELATION_CUTOFF)

cat(sprintf(
  "\nFeature pairs with |correlation| >= %.2f:\n",
  FEATURE_CORRELATION_CUTOFF
))
print(high_feature_correlations)

# -----------------------------------------------------------------------------
# 4. DBSCAN: unsupervised anomaly detection (unchanged -- no split involved)
# -----------------------------------------------------------------------------
x_all <- scale(weekly_df[, features])
minPts_value <- MINPTS_MULTIPLIER * length(features)

pdf(file.path(OUT_DIR, "DBSCAN_kNN_distance_plot.pdf"), width = 7, height = 5)
kNNdistplot(x_all, k = minPts_value)
abline(h = EPS_VALUE, col = "red", lty = 2)
title("k-NN distance plot")
dev.off()

db_fit <- dbscan(x_all, eps = EPS_VALUE, minPts = minPts_value)

weekly_df <- weekly_df %>%
  mutate(
    DBSCAN_Cluster = db_fit$cluster,
    DBSCAN_State = factor(
      if_else(DBSCAN_Cluster == 0, "Anomaly", "Clustered"),
      levels = c("Anomaly", "Clustered")
    )
  )

cat("\nDBSCAN cluster counts:\n")
print(table(weekly_df$DBSCAN_Cluster))
cat("\nDBSCAN state versus three-class VAL model label:\n")
print(table(weekly_df$Model_Label, weekly_df$DBSCAN_State))

db_obs <- factor(
  if_else(weekly_df$Model_Label_Num == 3L, "Class_3_eruption", "Class_1_or_2_non_eruption"),
  levels = c("Class_3_eruption", "Class_1_or_2_non_eruption")
)
db_pred <- factor(
  if_else(weekly_df$DBSCAN_State == "Anomaly", "Class_3_eruption", "Class_1_or_2_non_eruption"),
  levels = c("Class_3_eruption", "Class_1_or_2_non_eruption")
)
db_cm <- confusionMatrix(db_pred, db_obs, positive = "Class_3_eruption", mode = "everything")

dbscan_unrecorded_candidates <- weekly_df %>%
  filter(DBSCAN_State == "Anomaly", Model_Label_Num %in% c(1L, 2L)) %>%
  select(
    any_of(c("Week_Index", "Variable")),
    Week_Start, Week_End, Model_Label_Num, Model_Label, Activity_3Class,
    VAL_Levels, VAL_Descriptions, DBSCAN_Cluster, DBSCAN_State,
    all_of(features)
  )

cat("\nDBSCAN anomaly candidates not recorded as VAL eruption weeks:\n")
print(nrow(dbscan_unrecorded_candidates))

pca_fit <- prcomp(x_all)

point_labels <- if ("Week_Index" %in% names(weekly_df)) {
  paste0(
    "Week ", weekly_df$Week_Index,
    " | ", format(weekly_df$Week_Start, "%Y-%m-%d")
  )
} else {
  paste0("Week starting ", format(weekly_df$Week_Start, "%Y-%m-%d"))
}

pca_df <- weekly_df %>%
  mutate(
    PC1 = pca_fit$x[, 1],
    PC2 = pca_fit$x[, 2],
    Point_Label = point_labels
  )

p_dbscan <- ggplot(
  pca_df,
  aes(PC1, PC2, color = Activity_3Class, shape = DBSCAN_State)
) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_color_manual(values = activity_colors) +
  labs(
    title = "DBSCAN anomalies and VAL activity classes",
    color = "VAL class",
    shape = "DBSCAN state"
  ) +
  serif_theme

print(p_dbscan)
ggsave(
  file.path(OUT_DIR, "DBSCAN_VAL_PCA_plot.pdf"),
  p_dbscan, width = 7, height = 5
)

if (RUN_EASYLABEL) {
  easylabel(
    pca_df,
    x = "PC1",
    y = "PC2",
    col = "DBSCAN_State",
    labs = "Point_Label",
    size = 7
  )
}

sankey_df <- weekly_df %>%
  count(Activity_3Class, DBSCAN_State, name = "Weeks")

p_sankey <- ggplot(
  sankey_df,
  aes(axis1 = Activity_3Class, axis2 = DBSCAN_State, y = Weeks)
) +
  geom_alluvium(aes(fill = Activity_3Class), width = 0.2, alpha = 0.75) +
  geom_stratum(width = 0.2, fill = "grey95", color = "grey30") +
  geom_text(
    stat = "stratum", aes(label = after_stat(stratum)),
    family = "serif", size = 3
  ) +
  scale_x_discrete(
    limits = c("VAL activity class", "DBSCAN result"),
    expand = c(0.15, 0.15)
  ) +
  scale_fill_manual(values = activity_colors) +
  labs(
    title = "VAL activity class and DBSCAN anomaly status",
    x = NULL, y = "Number of weekly windows", fill = "VAL class"
  ) +
  serif_theme
print(p_sankey)
ggsave(
  file.path(OUT_DIR, "VAL_DBSCAN_Sankey.pdf"),
  p_sankey, width = 7, height = 5
)

# -----------------------------------------------------------------------------
# 5. Repeated stratified CV setup for supervised multiclass models
#    Class 3 is rare after 2014-07-01, so the effective fold count is capped by
#    the smallest class count.
# -----------------------------------------------------------------------------
all_x <- weekly_df[, features]
all_y <- factor(weekly_df$Model_Label, levels = model_label_levels)

all_counts <- table(all_y)
cat("\nFull-dataset class counts for CV:\n")
print(all_counts)
if (any(all_counts < 2)) {
  stop("A supervised model cannot be fitted: one class has fewer than two observations.")
}

CV_N_FOLDS <- min(N_FOLDS, as.integer(min(all_counts)))
if (CV_N_FOLDS < N_FOLDS) {
  cat(sprintf(
    "\nUsing %d folds instead of %d because the smallest class has %d weekly windows.\n",
    CV_N_FOLDS, N_FOLDS, min(all_counts)
  ))
}

class_weights <- setNames(
  length(all_y) / (length(model_label_levels) * as.numeric(all_counts)),
  names(all_counts)
)

weights_for_labels <- function(y) {
  unname(class_weights[as.character(y)])
}

metric_values <- function(cm, metric) {
  if (is.null(dim(cm$byClass))) {
    return(as.numeric(cm$byClass[metric]))
  }
  as.numeric(cm$byClass[, metric])
}

metric_mean <- function(cm, metric, na_to_zero = FALSE) {
  values <- metric_values(cm, metric)
  if (na_to_zero) {
    values[is.na(values)] <- 0
  }
  mean(values, na.rm = TRUE)
}

metrics_from_confusion <- function(cm) {
  data.frame(
    Accuracy = unname(cm$overall["Accuracy"]),
    Kappa = unname(cm$overall["Kappa"]),
    Mean_Balanced_Accuracy = metric_mean(cm, "Balanced Accuracy"),
    Mean_F1 = metric_mean(cm, "F1", na_to_zero = TRUE),
    Mean_Sensitivity = metric_mean(cm, "Sensitivity"),
    Mean_Specificity = metric_mean(cm, "Specificity")
  )
}

# -----------------------------------------------------------------------------
# 6. Weighted RBF SVM: cost/gamma grid searched via repeated stratified CV
# -----------------------------------------------------------------------------
svm_cost_grid  <- c(0.1, 1, 10, 100, 1000)
svm_gamma_grid <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5)
svm_grid <- expand.grid(cost = svm_cost_grid, gamma = svm_gamma_grid)

evaluate_svm <- function(x, y, cost, gamma, class_weights, n_folds, n_repeats) {
  acc_rep <- numeric(n_repeats); kappa_rep <- numeric(n_repeats)
  bal_acc_rep <- numeric(n_repeats); f1_rep <- numeric(n_repeats)
  sens_rep <- numeric(n_repeats); spec_rep <- numeric(n_repeats)
  
  for (r in seq_len(n_repeats)) {
    set.seed(1000 + r)
    folds <- createFolds(y, k = n_folds, list = TRUE)
    fa <- rep(NA_real_, n_folds); fk <- rep(NA_real_, n_folds)
    fb <- rep(NA_real_, n_folds); ff <- rep(NA_real_, n_folds)
    fs <- rep(NA_real_, n_folds); fp <- rep(NA_real_, n_folds)
    
    for (f in seq_len(n_folds)) {
      test_idx <- folds[[f]]
      train_data <- data.frame(Model_Label = y[-test_idx], x[-test_idx, , drop = FALSE])
      test_data  <- x[test_idx, , drop = FALSE]
      test_truth <- y[test_idx]
      
      fit <- svm(Model_Label ~ ., data = train_data, kernel = "radial",
                 cost = cost, gamma = gamma, class.weights = class_weights)
      pred <- factor(predict(fit, newdata = test_data), levels = model_label_levels)
      
      cm <- confusionMatrix(pred, test_truth, mode = "everything")
      fold_metrics <- metrics_from_confusion(cm)
      fa[f] <- fold_metrics$Accuracy
      fk[f] <- fold_metrics$Kappa
      fb[f] <- fold_metrics$Mean_Balanced_Accuracy
      ff[f] <- fold_metrics$Mean_F1
      fs[f] <- fold_metrics$Mean_Sensitivity
      fp[f] <- fold_metrics$Mean_Specificity
    }
    acc_rep[r] <- mean(fa, na.rm = TRUE); kappa_rep[r] <- mean(fk, na.rm = TRUE)
    bal_acc_rep[r] <- mean(fb, na.rm = TRUE); f1_rep[r] <- mean(ff, na.rm = TRUE)
    sens_rep[r] <- mean(fs, na.rm = TRUE); spec_rep[r] <- mean(fp, na.rm = TRUE)
  }
  
  data.frame(
    cost = cost, gamma = gamma,
    Mean_Accuracy = mean(acc_rep), Mean_Kappa = mean(kappa_rep),
    Mean_Balanced_Accuracy = mean(bal_acc_rep), SD_Balanced_Accuracy = sd(bal_acc_rep),
    Mean_F1 = mean(f1_rep), Mean_Sensitivity = mean(sens_rep), Mean_Specificity = mean(spec_rep)
  )
}

cat(sprintf("\nSearching SVM cost/gamma grid (%d combinations) via repeated stratified CV...\n",
            nrow(svm_grid)))
svm_grid_results <- do.call(rbind, lapply(seq_len(nrow(svm_grid)), function(i) {
  evaluate_svm(all_x, all_y, svm_grid$cost[i], svm_grid$gamma[i],
               class_weights, CV_N_FOLDS, SVM_N_REPEATS)
}))
svm_grid_results <- svm_grid_results %>% arrange(desc(Mean_Balanced_Accuracy))
best_cost  <- svm_grid_results$cost[1]
best_gamma <- svm_grid_results$gamma[1]
cat(sprintf("Selected SVM cost=%.3f, gamma=%.3f (mean Balanced Accuracy = %.3f)\n",
            best_cost, best_gamma, svm_grid_results$Mean_Balanced_Accuracy[1]))
print(head(svm_grid_results, 10))

# ---- Out-of-fold SVM predictions for the selected cost/gamma ----
set.seed(123456789, kind = "Mersenne-Twister")
cv_folds_best <- createFolds(all_y, k = CV_N_FOLDS, list = TRUE)
svm_oof_pred <- factor(rep(NA_character_, nrow(weekly_df)), levels = model_label_levels)

for (f in seq_len(CV_N_FOLDS)) {
  test_idx <- cv_folds_best[[f]]
  train_data <- data.frame(Model_Label = all_y[-test_idx], all_x[-test_idx, , drop = FALSE])
  test_data  <- all_x[test_idx, , drop = FALSE]
  
  fit <- svm(Model_Label ~ ., data = train_data, kernel = "radial",
             cost = best_cost, gamma = best_gamma, class.weights = class_weights)
  svm_oof_pred[test_idx] <- predict(fit, newdata = test_data)
}

svm_cm <- confusionMatrix(svm_oof_pred, all_y, mode = "everything")
svm_oof_metrics <- metrics_from_confusion(svm_cm)
cat(sprintf("\nSVM confusion matrix (out-of-fold, single %d-fold CV pass):\n", CV_N_FOLDS))
print(svm_cm)

# -----------------------------------------------------------------------------
# 7. Weighted multinomial logistic regression via repeated stratified CV
# -----------------------------------------------------------------------------
evaluate_multinom <- function(x, y, class_weights, n_folds, n_repeats) {
  acc_rep <- numeric(n_repeats); kappa_rep <- numeric(n_repeats)
  bal_acc_rep <- numeric(n_repeats); f1_rep <- numeric(n_repeats)
  sens_rep <- numeric(n_repeats); spec_rep <- numeric(n_repeats)
  
  for (r in seq_len(n_repeats)) {
    set.seed(1000 + r)
    folds <- createFolds(y, k = n_folds, list = TRUE)
    fa <- rep(NA_real_, n_folds); fk <- rep(NA_real_, n_folds)
    fb <- rep(NA_real_, n_folds); ff <- rep(NA_real_, n_folds)
    fs <- rep(NA_real_, n_folds); fp <- rep(NA_real_, n_folds)
    
    for (f in seq_len(n_folds)) {
      test_idx <- folds[[f]]
      train_x <- x[-test_idx, , drop = FALSE]; train_y <- y[-test_idx]
      test_x  <- x[test_idx, , drop = FALSE];  test_truth  <- y[test_idx]
      
      train_data <- data.frame(Model_Label = train_y, train_x)
      weights_train <- weights_for_labels(train_y)
      
      fit <- multinom(
        Model_Label ~ .,
        data = train_data,
        weights = weights_train,
        trace = FALSE,
        maxit = 1000
      )
      pred <- factor(predict(fit, newdata = data.frame(test_x)), levels = model_label_levels)
      
      cm <- confusionMatrix(pred, test_truth, mode = "everything")
      fold_metrics <- metrics_from_confusion(cm)
      fa[f] <- fold_metrics$Accuracy
      fk[f] <- fold_metrics$Kappa
      fb[f] <- fold_metrics$Mean_Balanced_Accuracy
      ff[f] <- fold_metrics$Mean_F1
      fs[f] <- fold_metrics$Mean_Sensitivity
      fp[f] <- fold_metrics$Mean_Specificity
    }
    acc_rep[r] <- mean(fa, na.rm = TRUE); kappa_rep[r] <- mean(fk, na.rm = TRUE)
    bal_acc_rep[r] <- mean(fb, na.rm = TRUE); f1_rep[r] <- mean(ff, na.rm = TRUE)
    sens_rep[r] <- mean(fs, na.rm = TRUE); spec_rep[r] <- mean(fp, na.rm = TRUE)
  }
  
  data.frame(
    Mean_Accuracy = mean(acc_rep), Mean_Kappa = mean(kappa_rep),
    Mean_Balanced_Accuracy = mean(bal_acc_rep), SD_Balanced_Accuracy = sd(bal_acc_rep),
    Mean_F1 = mean(f1_rep), Mean_Sensitivity = mean(sens_rep), Mean_Specificity = mean(spec_rep)
  )
}

cat("\nEvaluating weighted multinomial logistic regression via repeated stratified CV...\n")
multinom_results <- evaluate_multinom(
  all_x, all_y, class_weights, CV_N_FOLDS, MULTINOM_N_REPEATS
)
cat(sprintf(
  "Weighted multinomial logistic regression -> BalAcc=%.3f  Sens=%.3f  Spec=%.3f  F1=%.3f\n",
  multinom_results$Mean_Balanced_Accuracy, multinom_results$Mean_Sensitivity,
  multinom_results$Mean_Specificity, multinom_results$Mean_F1
))

# ---- Out-of-fold multinomial predictions (same folds as SVM, for comparability) ----
multinom_oof_pred <- factor(rep(NA_character_, nrow(weekly_df)), levels = model_label_levels)

for (f in seq_len(CV_N_FOLDS)) {
  test_idx <- cv_folds_best[[f]]
  train_x <- all_x[-test_idx, , drop = FALSE]; train_y <- all_y[-test_idx]
  test_x  <- all_x[test_idx, , drop = FALSE]
  
  train_data <- data.frame(Model_Label = train_y, train_x)
  weights_train <- weights_for_labels(train_y)
  
  fit <- multinom(
    Model_Label ~ .,
    data = train_data,
    weights = weights_train,
    trace = FALSE,
    maxit = 1000
  )
  multinom_oof_pred[test_idx] <- predict(fit, newdata = data.frame(test_x))
}

multinom_cm <- confusionMatrix(multinom_oof_pred, all_y, mode = "everything")
multinom_oof_metrics <- metrics_from_confusion(multinom_cm)
cat(sprintf(
  "\nWeighted multinomial logistic regression confusion matrix (out-of-fold, single %d-fold CV pass):\n",
  CV_N_FOLDS
))
print(multinom_cm)

# -----------------------------------------------------------------------------
# 7a. PCA displays of out-of-fold SVM and multinomial predictions
# -----------------------------------------------------------------------------
# PCA is calculated only once from the seven standardized features. The model
# predictions are then displayed on the same two-dimensional projection.
# The classifiers themselves still use all seven original feature dimensions.
model_pca_df <- pca_df %>%
  mutate(
    Ground_Truth = factor(all_y, levels = model_label_levels),
    SVM_Prediction = factor(svm_oof_pred, levels = model_label_levels),
    Multinom_Prediction = factor(multinom_oof_pred, levels = model_label_levels)
  )

truth_shapes <- c(
  "1_Background_unrest" = 16,
  "2_Moderate_to_heightened_unrest" = 17,
  "3_Minor_to_moderate_eruption" = 15
)

p_svm_pca <- ggplot(
  model_pca_df,
  aes(PC1, PC2, color = SVM_Prediction, shape = Ground_Truth)
) +
  geom_point(size = 1.8, alpha = 0.72) +
  scale_color_manual(values = model_label_colors) +
  scale_shape_manual(values = truth_shapes) +
  labs(
    title = "Weighted RBF SVM three-class predictions in PCA space",
    x = "PC1", y = "PC2",
    color = "Out-of-fold prediction",
    shape = "VAL ground truth"
  ) +
  serif_theme

print(p_svm_pca)
ggsave(
  file.path(OUT_DIR, "SVM_OOF_Predictions_PCA.pdf"),
  p_svm_pca, width = 7, height = 5
)

if (RUN_EASYLABEL) {
  easylabel(
    model_pca_df,
    x = "PC1",
    y = "PC2",
    col = "SVM_Prediction",
    labs = "Point_Label",
    size = 7
  )
}

p_multinom_pca <- ggplot(
  model_pca_df,
  aes(PC1, PC2, color = Multinom_Prediction, shape = Ground_Truth)
) +
  geom_point(size = 1.8, alpha = 0.72) +
  scale_color_manual(values = model_label_colors) +
  scale_shape_manual(values = truth_shapes) +
  labs(
    title = "Weighted multinomial logistic predictions in PCA space",
    x = "PC1", y = "PC2",
    color = "Out-of-fold prediction",
    shape = "VAL ground truth"
  ) +
  serif_theme

print(p_multinom_pca)
ggsave(
  file.path(OUT_DIR, "Multinom_OOF_Predictions_PCA.pdf"),
  p_multinom_pca, width = 7, height = 5
)

if (RUN_EASYLABEL) {
  easylabel(
    model_pca_df,
    x = "PC1",
    y = "PC2",
    col = "Multinom_Prediction",
    labs = "Point_Label",
    size = 7
  )
}

# Full-dataset refit, kept for coefficient interpretation only (not for the
# performance metrics above, which come from CV).
full_weights <- weights_for_labels(all_y)
multinom_fit_full <- multinom(
  Model_Label ~ .,
  data = data.frame(Model_Label = all_y, all_x),
  weights = full_weights,
  trace = FALSE,
  maxit = 1000
)

# -----------------------------------------------------------------------------
# 8. Collect performance measures and save all results
# -----------------------------------------------------------------------------
db_metrics <- metrics_from_confusion(db_cm)
comparison <- data.frame(
  Model = c(
    "DBSCAN anomaly flag",
    "Weighted RBF SVM (CV)",
    "Weighted multinomial logistic regression (CV)"
  ),
  Target = c(
    "Class 3 eruption versus class 1/2",
    "Three-class VAL model label",
    "Three-class VAL model label"
  ),
  Accuracy = c(
    db_metrics$Accuracy,
    svm_grid_results$Mean_Accuracy[1],
    multinom_results$Mean_Accuracy
  ),
  Kappa = c(
    db_metrics$Kappa,
    svm_grid_results$Mean_Kappa[1],
    multinom_results$Mean_Kappa
  ),
  Mean_Balanced_Accuracy = c(
    db_metrics$Mean_Balanced_Accuracy,
    svm_grid_results$Mean_Balanced_Accuracy[1],
    multinom_results$Mean_Balanced_Accuracy
  ),
  Mean_Sensitivity = c(
    db_metrics$Mean_Sensitivity,
    svm_grid_results$Mean_Sensitivity[1],
    multinom_results$Mean_Sensitivity
  ),
  Mean_Specificity = c(
    db_metrics$Mean_Specificity,
    svm_grid_results$Mean_Specificity[1],
    multinom_results$Mean_Specificity
  ),
  Mean_F1 = c(
    db_metrics$Mean_F1,
    svm_grid_results$Mean_F1[1],
    multinom_results$Mean_F1
  )
)
cat("\n==== Model comparison (DBSCAN: full-data; SVM/multinomial: repeated-CV mean) ====\n")
print(comparison, row.names = FALSE)

confusion_to_df <- function(cm) {
  as.data.frame.matrix(cm$table) %>%
    rownames_to_column("Prediction")
}

write_xlsx(
  list(
    Weekly_Labels_DBSCAN   = weekly_df,
    Model_Comparison       = comparison,
    DBSCAN_Confusion       = confusion_to_df(db_cm),
    SVM_Confusion_OOF      = confusion_to_df(svm_cm),
    Multinom_Confusion_OOF = confusion_to_df(multinom_cm),
    SVM_Grid_Results       = svm_grid_results,
    Feature_Correlation    = as.data.frame(feature_cor) %>% rownames_to_column("Feature"),
    High_Correlation_Pairs = high_feature_correlations,
    DBSCAN_Unrecorded      = dbscan_unrecorded_candidates,
    VAL_System_Sankey      = val_system_sankey,
    VAL_Activity_Sankey    = val_activity_sankey,
    VAL_DBSCAN_Sankey      = sankey_df
  ),
  path = file.path(OUT_DIR, "OP_VAL_System3_DBSCAN_SVM_Multinom_results.xlsx")
)

capture.output(
  summary(multinom_fit_full),
  file = file.path(OUT_DIR, "weighted_multinomial_logistic_summary.txt")
)

cat("\nAll outputs saved in:\n", OUT_DIR, "\n")
