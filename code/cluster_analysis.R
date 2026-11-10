# Weekly ordinal-pattern analysis validated against Whakaari VAL records
# Analysis period: VAL system version 3, 2014-07-01 to 2019-12-31
# Run separately for TWO stations: WIZ and WSRZ (same method, same VAL
# ground truth, independent weekly feature files)
# Model labels are:
#   1 = background unrest
#   2 = moderate to heightened volcanic unrest
#   3 = minor to moderate volcanic eruption
#
# NOTE: multinomial logistic regression has been removed from this version.
# Only DBSCAN (unsupervised anomaly detection) and weighted RBF SVM
# (supervised three-class classification) are evaluated, per the decision
# that logistic regression was not competitive with the SVM.

library(tidyverse)
library(dbscan)
library(e1071)
library(caret)
library(ggalluvial)
library(writexl)
library(easylabel)

set.seed(123456789, kind = "Mersenne-Twister")

# -----------------------------------------------------------------------------
# 1. Paths and settings
# -----------------------------------------------------------------------------
DATA_DIR    <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"
RESULTS_DIR <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"

CLUSTER_DIR <- file.path(RESULTS_DIR, "cluster_analysis")
dir.create(CLUSTER_DIR, recursive = TRUE, showWarnings = FALSE)

VAL_FILE <- file.path(DATA_DIR, "VAL_WhakaariWhiteIsland.csv")

# One weekly OP-feature file per station -- adjust these two paths if your
# actual filenames differ.
station_files <- c(
  WIZ  = file.path(DATA_DIR, "weekly_HC_RSAM_WIZ_20140701_to_20191231.csv"),
  WSRZ = file.path(DATA_DIR, "weekly_HC_RSAM_WSRZ_20140701_to_20191231.csv")
)

ANALYSIS_START_LOCAL <- as.POSIXct("2014-07-01 00:00:00", tz = "UTC")  # NZ local; tz="UTC" is just a storage tag, not a conversion
ANALYSIS_END_LOCAL   <- as.POSIXct("2019-12-31 23:59:59", tz = "UTC")  # NZ local; tz="UTC" is just a storage tag, not a conversion
VAL_SYSTEM_VERSION <- 3
FEATURE_CORRELATION_CUTOFF <- 0.90

features <- c(
  "Fisher_Info", "C_Fisher",
  "Var_H_Shannon", "Var_C_Shannon",
  "Var_H_Tsallis", "Var_H_Renyi",
  "Disequilibrium"
)

# The 7 features above are used throughout this script -- correlation
# analysis, DBSCAN, the SVM, and the pairwise correlation plot below:
#   Fisher_Info    -- Fisher information
#   C_Fisher       -- Fisher-based statistical complexity
#   Var_H_Shannon  -- variance of Shannon entropy
#   Var_C_Shannon  -- variance of Shannon (Jensen-Shannon) complexity
#   Var_H_Tsallis  -- variance of Tsallis entropy
#   Var_H_Renyi    -- variance of Renyi entropy
#   Disequilibrium -- ordinal-pattern disequilibrium
cat(sprintf("Features used in this analysis (n = %d):\n", length(features)))
cat(paste(" -", features, collapse = "\n"), "\n\n")

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

# Repeated stratified CV settings
N_FOLDS <- 5
SVM_N_REPEATS <- 10  # lower because the cost/gamma grid is large;
# raise to 20 if runtime allows

# Set to TRUE when running interactively in RStudio. Each easylabel window
# opens in turn; close the current window to continue to the next PCA plot.
RUN_EASYLABEL <- interactive()

serif_theme <- theme_bw(base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    legend.position = "bottom",
    legend.box = "vertical",      # stack separate legends (e.g. color + shape)
    legend.margin = margin(t = 0, b = 0),
    legend.box.spacing = unit(0.15, "cm"),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank()
  )

# Reusable guide spec: wraps any color/fill/shape legend onto two rows instead
# of one, so long legends (e.g. VAL descriptions, 3-4 category labels) don't
# get clipped at the plot edge. Add `+ legend_two_rows` to every plot -- if a
# given aesthetic isn't used in that plot, ggplot silently ignores that part
# of the guides() call.
legend_two_rows <- guides(
  color = guide_legend(nrow = 2, byrow = TRUE),
  fill  = guide_legend(nrow = 2, byrow = TRUE),
  shape = guide_legend(nrow = 2, byrow = TRUE)
)

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

truth_shapes <- c(
  "1_Background_unrest" = 16,
  "2_Moderate_to_heightened_unrest" = 17,
  "3_Minor_to_moderate_eruption" = 15
)

# -----------------------------------------------------------------------------
# 2. Read VAL system 3 records ONCE and define model labels
#    (doesn't depend on station -- shared across both runs)
# -----------------------------------------------------------------------------
# Using the "PeriodStart LocalTime"/"PeriodEnd LocalTime" columns (F & G in
# the source file) rather than the UTC ones, since the weekly OP-feature
# files' Week_Start/Week_End are themselves NZ local time -- comparing
# local-to-local avoids any UTC conversion (and the DST-transition
# arithmetic that comes with it) entirely. Each value ends with a literal
# " NZST" or " NZDT" suffix (e.g. "2016-04-28 11:30:00 NZST"); that suffix
# is just a label confirming which one was in effect and isn't needed for
# comparison, so it's stripped before parsing.
val_df <- read_csv(VAL_FILE, show_col_types = FALSE) %>%
  mutate(
    PeriodStart_Local = as.POSIXct(
      stringr::str_remove(`PeriodStart LocalTime`, " NZ[SD]T$"),
      format = "%Y-%m-%d %H:%M:%S", tz = "UTC"  # tz="UTC" is just a storage tag; values are NZ local
    ),
    PeriodEnd_Local = as.POSIXct(
      stringr::str_remove(`PeriodEnd LocalTime`, " NZ[SD]T$"),
      format = "%Y-%m-%d %H:%M:%S", tz = "UTC"
    )
  )

# Open-ended "9999" sentinel periods are capped at ANALYSIS_END_LOCAL, since
# that's the last date we'll ever need VAL coverage for.
val_df$PeriodEnd_Local[
  is.na(val_df$PeriodEnd_Local) |
    stringr::str_starts(val_df$`PeriodEnd LocalTime`, "9999")
] <- ANALYSIS_END_LOCAL

val_df <- val_df %>%
  filter(
    `VAL System Version` == VAL_SYSTEM_VERSION,
    PeriodEnd_Local > ANALYSIS_START_LOCAL,
    PeriodStart_Local < ANALYSIS_END_LOCAL
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
# 2a. Descriptive Sankey plots for the VAL system 3 record (shared, saved
#     once at the top level of CLUSTER_DIR -- not station-specific)
# -----------------------------------------------------------------------------
val_system_sankey <- val_df %>%
  transmute(
    VAL_Level = paste("VAL", VAL),
    Description = stringr::str_wrap(Description, width = 25)
  ) %>%
  count(VAL_Level, Description, name = "Number_of_periods")

p_val_system_sankey <- ggplot(
  val_system_sankey,
  aes(
    axis1 = VAL_Level,
    axis2 = Description,
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
    limits = c("Volcanic Alert Level", "Activity description"),
    expand = c(0.12, 0.12)
  ) +
  scale_fill_manual(values = val_colors) +
  labs(
    title = "VAL Level and Description",
    x = NULL, y = "Number of recorded VAL periods", fill = NULL
  ) +
  serif_theme +
  legend_two_rows

print(p_val_system_sankey)
ggsave(
  file.path(CLUSTER_DIR, "VAL_System_Description_Sankey.pdf"),
  p_val_system_sankey, width = 10, height = 5.5
)

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
    title = "VAL Level and Model Class",
    x = NULL, y = "Number of recorded VAL periods", fill = NULL
  ) +
  serif_theme +
  legend_two_rows

print(p_val_activity_sankey)
ggsave(
  file.path(CLUSTER_DIR, "VAL_Activity_Class_Sankey.pdf"),
  p_val_activity_sankey, width = 9, height = 5.5
)

assign_val_label <- function(week_start, week_end, val_table) {
  overlap <- val_table$PeriodStart_Local < week_end &
    val_table$PeriodEnd_Local > week_start
  
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

confusion_to_df <- function(cm) {
  as.data.frame.matrix(cm$table) %>%
    rownames_to_column("Prediction")
}

# NOTE: an earlier version of this script converted Week_Start/Week_End to
# UTC using a hand-rolled NZ daylight-saving calendar rule, to compare
# against the VAL file's UTC columns. That's no longer needed -- the VAL
# file's "PeriodStart LocalTime"/"PeriodEnd LocalTime" columns are used
# instead (see val_df above), so both sides of the comparison are already
# in NZ local time and no conversion happens at all.

# Pairwise feature correlation heatmap. Always uses Spearman correlation
# (rank-based, robust to the non-normal/bounded distributions typical of
# entropy and complexity measures), regardless of what method the
# High_Correlation_Pairs table elsewhere in this script uses.
plot_correlation_heatmap <- function(d, station_name, feature_cols, method = "spearman") {
  cor_mat <- cor(d[, feature_cols], method = method, use = "pairwise.complete.obs")
  
  cor_df <- as.data.frame(cor_mat) %>%
    rownames_to_column("Feature_1") %>%
    pivot_longer(-Feature_1, names_to = "Feature_2", values_to = "Correlation") %>%
    mutate(
      Feature_1 = factor(Feature_1, levels = feature_cols),
      Feature_2 = factor(Feature_2, levels = feature_cols)
    )
  
  ggplot(cor_df, aes(x = Feature_2, y = Feature_1, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Correlation)), family = "serif", size = 3) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho"
    ) +
    labs(title = paste0(station_name, ": Feature Correlation")) +
    theme_bw(base_family = "serif") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

# =============================================================================
# Per-station analysis: DBSCAN + weighted RBF SVM (repeated CV), against the
# shared VAL ground truth built above.
# =============================================================================
run_station_cluster_analysis <- function(station_name, op_file, val_df, out_dir) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("\n\n########## %s ##########\n", station_name))
  
  # ---- Read weekly ordinal-pattern features for this station ----
  # Week_Start/Week_End in these files are NZ local time (NZST/NZDT
  # depending on the date). val_df above now uses VAL's own local-time
  # columns ("PeriodStart LocalTime"/"PeriodEnd LocalTime") rather than
  # its UTC columns, so both sides of every comparison below are already
  # in NZ local time -- no timezone conversion is needed at all.
  weekly_df <- read_csv(op_file, show_col_types = FALSE) %>%
    mutate(
      Week_Start = as.POSIXct(Week_Start, tz = "UTC"),  # tz="UTC" is just a storage tag; values are NZ local
      Week_End   = as.POSIXct(Week_End,   tz = "UTC")
    ) %>%
    filter(Week_Start >= ANALYSIS_START_LOCAL, Week_Start <= ANALYSIS_END_LOCAL) %>%
    filter(if_all(all_of(features), is.finite))
  
  if (anyNA(weekly_df$Week_Start) || anyNA(weekly_df$Week_End)) {
    stop(sprintf("[%s] Week_Start or Week_End could not be converted to UTC date-times.", station_name))
  }
  
  # ---- Assign VAL model labels to each weekly window ----
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
      Activity_3Class = factor(Activity_3Class, levels = activity_levels),
      Model_Label = factor(Model_Label, levels = model_label_levels)
    ) %>%
    filter(!is.na(Activity_3Class), !is.na(Model_Label), !is.na(Model_Label_Num))
  
  cat(sprintf("[%s] VAL system 3 model labels by weekly window:\n", station_name))
  print(table(weekly_df$Model_Label, useNA = "ifany"))
  cat(sprintf("\n[%s] Activity labels by weekly window:\n", station_name))
  print(table(weekly_df$Activity_3Class, useNA = "ifany"))
  cat(sprintf("\n[%s] Weekly data range used for modelling:\n", station_name))
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
    "\n[%s] Feature pairs with |correlation| >= %.2f:\n",
    station_name, FEATURE_CORRELATION_CUTOFF
  ))
  print(high_feature_correlations)
  
  # ---------------------------------------------------------------------------
  # Pairwise feature correlation heatmap (Spearman). Note this is a separate
  # computation from feature_cor above, which is Pearson (used only for the
  # >= 0.90 cutoff table and the Excel export) -- tell me if you'd rather
  # unify both on Spearman for consistency.
  # ---------------------------------------------------------------------------
  p_corr <- plot_correlation_heatmap(weekly_df, station_name, features, method = "spearman")
  print(p_corr)
  ggsave(
    file.path(out_dir, sprintf("Feature_Correlation_Spearman_%s.pdf", station_name)),
    p_corr, width = 7, height = 6
  )
  
  # ---------------------------------------------------------------------------
  # DBSCAN: unsupervised anomaly detection (full data, no CV split)
  # ---------------------------------------------------------------------------
  x_all <- scale(weekly_df[, features])
  minPts_value <- MINPTS_MULTIPLIER * length(features)
  
  pdf(file.path(out_dir, sprintf("DBSCAN_kNN_distance_plot_%s.pdf", station_name)), width = 7, height = 5)
  kNNdistplot(x_all, k = minPts_value)
  abline(h = EPS_VALUE, col = "red", lty = 2)
  title(sprintf("%s: k-NN distance plot", station_name))
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
  
  cat(sprintf("\n[%s] DBSCAN cluster counts:\n", station_name))
  print(table(weekly_df$DBSCAN_Cluster))
  cat(sprintf("\n[%s] DBSCAN state versus three-class VAL model label:\n", station_name))
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
  
  cat(sprintf("\n[%s] DBSCAN anomaly candidates not recorded as VAL eruption weeks: %d\n",
              station_name, nrow(dbscan_unrecorded_candidates)))
  
  pca_fit <- prcomp(x_all)
  
  # Percentage of total variance explained by each PC, for axis labels
  # (e.g. "PC1 (42.3% variance)") -- computed once per station and reused
  # on every PCA scatter plot below.
  pca_var_pct <- (pca_fit$sdev^2 / sum(pca_fit$sdev^2)) * 100
  pc1_axis_label <- sprintf("PC1 (%.1f%%)", pca_var_pct[1])
  pc2_axis_label <- sprintf("PC2 (%.1f%%)", pca_var_pct[2])
  cat(sprintf(
    "\n[%s] PCA variance explained: PC1 = %.1f%%, PC2 = %.1f%% (cumulative %.1f%%)\n",
    station_name, pca_var_pct[1], pca_var_pct[2], pca_var_pct[1] + pca_var_pct[2]
  ))
  
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
      title = paste0(station_name, ": DBSCAN"),
      x = pc1_axis_label,
      y = pc2_axis_label,
      color = "VAL class",
      shape = "DBSCAN state"
    ) +
    serif_theme +
    legend_two_rows
  
  print(p_dbscan)
  ggsave(
    file.path(out_dir, sprintf("DBSCAN_VAL_PCA_plot_%s.pdf", station_name)),
    p_dbscan, width = 8.5, height = 5.5
  )
  
  if (RUN_EASYLABEL) {
    easylabel(pca_df, x = "PC1", y = "PC2", col = "DBSCAN_State", labs = "Point_Label", size = 7)
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
      title = paste0(station_name, ": VAL vs. DBSCAN"),
      x = NULL, y = "Number of weekly windows", fill = "VAL class"
    ) +
    serif_theme +
    legend_two_rows
  print(p_sankey)
  ggsave(
    file.path(out_dir, sprintf("VAL_DBSCAN_Sankey_%s.pdf", station_name)),
    p_sankey, width = 8, height = 5.5
  )
  
  # ---------------------------------------------------------------------------
  # Repeated stratified CV setup for the supervised multiclass model
  # ---------------------------------------------------------------------------
  all_x <- weekly_df[, features]
  all_y <- factor(weekly_df$Model_Label, levels = model_label_levels)
  
  all_counts <- table(all_y)
  cat(sprintf("\n[%s] Full-dataset class counts for CV:\n", station_name))
  print(all_counts)
  if (any(all_counts < 2)) {
    stop(sprintf("[%s] A supervised model cannot be fitted: one class has fewer than two observations.", station_name))
  }
  
  CV_N_FOLDS <- min(N_FOLDS, as.integer(min(all_counts)))
  if (CV_N_FOLDS < N_FOLDS) {
    cat(sprintf(
      "\n[%s] Using %d folds instead of %d because the smallest class has %d weekly windows.\n",
      station_name, CV_N_FOLDS, N_FOLDS, min(all_counts)
    ))
  }
  
  class_weights <- setNames(
    length(all_y) / (length(model_label_levels) * as.numeric(all_counts)),
    names(all_counts)
  )
  
  weights_for_labels <- function(y) {
    unname(class_weights[as.character(y)])
  }
  
  # ---------------------------------------------------------------------------
  # Weighted RBF SVM: cost/gamma grid searched via repeated stratified CV
  # ---------------------------------------------------------------------------
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
  
  cat(sprintf("\n[%s] Searching SVM cost/gamma grid (%d combinations) via repeated stratified CV...\n",
              station_name, nrow(svm_grid)))
  svm_grid_results <- do.call(rbind, lapply(seq_len(nrow(svm_grid)), function(i) {
    evaluate_svm(all_x, all_y, svm_grid$cost[i], svm_grid$gamma[i],
                 class_weights, CV_N_FOLDS, SVM_N_REPEATS)
  }))
  svm_grid_results <- svm_grid_results %>% arrange(desc(Mean_Balanced_Accuracy))
  best_cost  <- svm_grid_results$cost[1]
  best_gamma <- svm_grid_results$gamma[1]
  cat(sprintf("[%s] Selected SVM cost=%.3f, gamma=%.3f (mean Balanced Accuracy = %.3f)\n",
              station_name, best_cost, best_gamma, svm_grid_results$Mean_Balanced_Accuracy[1]))
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
  cat(sprintf("\n[%s] SVM confusion matrix (out-of-fold, single %d-fold CV pass):\n", station_name, CV_N_FOLDS))
  print(svm_cm)
  
  # ---------------------------------------------------------------------------
  # PCA display of out-of-fold SVM predictions
  # ---------------------------------------------------------------------------
  model_pca_df <- pca_df %>%
    mutate(
      Ground_Truth = factor(all_y, levels = model_label_levels),
      SVM_Prediction = factor(svm_oof_pred, levels = model_label_levels)
    )
  
  p_svm_pca <- ggplot(
    model_pca_df,
    aes(PC1, PC2, color = SVM_Prediction, shape = Ground_Truth)
  ) +
    geom_point(size = 1.8, alpha = 0.72) +
    scale_color_manual(values = model_label_colors) +
    scale_shape_manual(values = truth_shapes) +
    labs(
      title = paste0(station_name, ": SVM"),
      x = pc1_axis_label, y = pc2_axis_label,
      color = "SVM prediction",
      shape = "VAL truth"
    ) +
    serif_theme +
    legend_two_rows
  
  print(p_svm_pca)
  ggsave(
    file.path(out_dir, sprintf("SVM_OOF_Predictions_PCA_%s.pdf", station_name)),
    p_svm_pca, width = 8.5, height = 5.5
  )
  
  if (RUN_EASYLABEL) {
    easylabel(model_pca_df, x = "PC1", y = "PC2", col = "SVM_Prediction", labs = "Point_Label", size = 7)
  }
  
  # ---------------------------------------------------------------------------
  # Collect performance measures and save all results for this station
  # ---------------------------------------------------------------------------
  db_metrics <- metrics_from_confusion(db_cm)
  comparison <- data.frame(
    Model = c(
      "DBSCAN anomaly flag",
      "Weighted RBF SVM (CV)"
    ),
    Target = c(
      "Class 3 eruption versus class 1/2",
      "Three-class VAL model label"
    ),
    Accuracy = c(db_metrics$Accuracy, svm_grid_results$Mean_Accuracy[1]),
    Kappa = c(db_metrics$Kappa, svm_grid_results$Mean_Kappa[1]),
    Mean_Balanced_Accuracy = c(db_metrics$Mean_Balanced_Accuracy, svm_grid_results$Mean_Balanced_Accuracy[1]),
    Mean_Sensitivity = c(db_metrics$Mean_Sensitivity, svm_grid_results$Mean_Sensitivity[1]),
    Mean_Specificity = c(db_metrics$Mean_Specificity, svm_grid_results$Mean_Specificity[1]),
    Mean_F1 = c(db_metrics$Mean_F1, svm_grid_results$Mean_F1[1])
  )
  cat(sprintf("\n==== [%s] Model comparison (DBSCAN: full-data; SVM: repeated-CV mean) ====\n", station_name))
  print(comparison, row.names = FALSE)
  
  write_xlsx(
    list(
      Weekly_Labels_DBSCAN   = weekly_df,
      Model_Comparison       = comparison,
      DBSCAN_Confusion       = confusion_to_df(db_cm),
      SVM_Confusion_OOF      = confusion_to_df(svm_cm),
      SVM_Grid_Results       = svm_grid_results,
      Feature_Correlation    = as.data.frame(feature_cor) %>% rownames_to_column("Feature"),
      High_Correlation_Pairs = high_feature_correlations,
      DBSCAN_Unrecorded      = dbscan_unrecorded_candidates,
      VAL_DBSCAN_Sankey      = sankey_df
    ),
    path = file.path(out_dir, sprintf("OP_VAL_System3_DBSCAN_SVM_results_%s.xlsx", station_name))
  )
  
  cat(sprintf("\n[%s] outputs saved in:\n  %s\n", station_name, out_dir))
  
  invisible(list(
    weekly_df = weekly_df, comparison = comparison,
    db_cm = db_cm, svm_cm = svm_cm
  ))
}

# -----------------------------------------------------------------------------
# Run for both stations
# -----------------------------------------------------------------------------
station_results <- list()
for (station_name in names(station_files)) {
  out_dir <- file.path(CLUSTER_DIR, station_name)
  station_results[[station_name]] <- run_station_cluster_analysis(
    station_name = station_name,
    op_file      = station_files[[station_name]],
    val_df       = val_df,
    out_dir      = out_dir
  )
}

cat("\nAll stations processed. Results under:\n", CLUSTER_DIR, "\n")


##############################################################################

#Checking the number of rows that used for this analysis
for (st in names(station_files)) {
  raw <- read_csv(station_files[st], show_col_types = FALSE)
  cat(sprintf("[%s] raw row count in CSV file: %d\n", st, nrow(raw)))
  
  n_bad_feature_rows <- raw %>%
    filter(!if_all(all_of(features), is.finite)) %>%
    nrow()
  cat(sprintf("[%s] rows that would be DROPPED by if_all(...is.finite) due to at least one NA feature: %d\n",
              st, n_bad_feature_rows))
  
  if (n_bad_feature_rows > 0) {
    cat(sprintf("[%s] which weeks, and which feature is NA:\n", st))
    raw %>%
      filter(!if_all(all_of(features), is.finite)) %>%
      select(Week_Start, Week_End, all_of(features)) %>%
      print(n = Inf)
  }
}
