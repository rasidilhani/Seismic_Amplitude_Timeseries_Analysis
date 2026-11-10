# =============================================================================
# Weekly ordinal-pattern analysis validated against Whakaari VAL records
# Analysis period: VAL system version 3, 2014-07-01 to 2019-12-31
#
# Run separately for TWO stations:
#   WIZ
#   WSRZ
#
# Models:
#   1 = background unrest
#   2 = moderate to heightened volcanic unrest
#   3 = minor to moderate volcanic eruption
#
# DBSCAN:
#   Unsupervised THREE-cluster analysis
#   Cluster 0 = Noise
#   Cluster 1 = DBSCAN Cluster 1
#   Cluster 2 = DBSCAN Cluster 2
#   Cluster 3 = DBSCAN Cluster 3
#
# SVM:
#   Supervised three-class weighted RBF SVM
# =============================================================================


library(tidyverse)
library(dbscan)
library(e1071)
library(caret)
library(ggalluvial)
library(writexl)
library(easylabel)

set.seed(123456789, kind = "Mersenne-Twister")


# =============================================================================
# 1. Paths and settings
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
  "cluster_analysis"
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


# Weekly OP-feature file for each station
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


ANALYSIS_START_LOCAL <- as.POSIXct(
  "2014-07-01 00:00:00",
  tz = "Pacific/Auckland"
)

ANALYSIS_END_LOCAL <- as.POSIXct(
  "2019-12-31 23:59:59",
  tz = "Pacific/Auckland"
)

VAL_SYSTEM_VERSION <- 3

FEATURE_CORRELATION_CUTOFF <- 0.90


# =============================================================================
# Features
# =============================================================================

features <- c(
  
  "Fisher_Info",
  "C_Fisher",
  
  "Var_H_Shannon",
  "Var_C_Shannon",
  
  "Var_H_Tsallis",
  "Var_H_Renyi",
  
  "Disequilibrium"
  
)


cat(
  sprintf(
    "Features used in this analysis (n = %d):\n",
    length(features)
  )
)

cat(
  paste(
    " -",
    features,
    collapse = "\n"
  ),
  "\n\n"
)


# =============================================================================
# Three VAL activity classes
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


activity_to_model_level <- function(activity) {
  
  case_when(
    
    activity == "Background unrest" ~ 1L,
    
    activity ==
      "Moderate to heightened unrest" ~ 2L,
    
    activity ==
      "Minor to moderate eruption" ~ 3L,
    
    TRUE ~ NA_integer_
    
  )
  
}


model_level_to_activity <- function(level) {
  
  case_when(
    
    level == 1L ~
      "Background unrest",
    
    level == 2L ~
      "Moderate to heightened unrest",
    
    level == 3L ~
      "Minor to moderate eruption",
    
    TRUE ~ NA_character_
    
  )
  
}


model_level_to_label <- function(level) {
  
  case_when(
    
    level == 1L ~
      "1_Background_unrest",
    
    level == 2L ~
      "2_Moderate_to_heightened_unrest",
    
    level == 3L ~
      "3_Minor_to_moderate_eruption",
    
    TRUE ~ NA_character_
    
  )
  
}


# =============================================================================
# DBSCAN settings
# =============================================================================

# Your original eps value.
# It is retained as the preferred value.
# The search below finds eps values producing exactly THREE clusters.
EPS_VALUE <- 1.1

MINPTS_MULTIPLIER <- 2

# Search ranges for a three-cluster DBSCAN solution
DBSCAN_EPS_GRID <- seq(
  0.05,
  5.00,
  by = 0.05
)

# Search several minPts values as well.
# This is needed because DBSCAN does not guarantee a fixed number of clusters.
DBSCAN_MINPTS_GRID <- 3:20


# =============================================================================
# Repeated stratified CV settings
# =============================================================================

N_FOLDS <- 5

SVM_N_REPEATS <- 10


# =============================================================================
# Easylabel
# =============================================================================

RUN_EASYLABEL <- interactive()


# =============================================================================
# Plot theme
# =============================================================================

serif_theme <- theme_bw(
  base_family = "serif"
) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 11
    ),
    
    legend.position = "bottom",
    
    legend.box = "vertical",
    
    legend.margin = margin(
      t = 0,
      b = 0
    ),
    
    legend.box.spacing = unit(
      0.15,
      "cm"
    ),
    
    legend.key.size = unit(
      0.4,
      "cm"
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    legend.title = element_text(
      size = 9,
      face = "bold"
    ),
    
    panel.grid.minor =
      element_blank()
    
  )


legend_two_rows <- guides(
  
  color = guide_legend(
    nrow = 2,
    byrow = TRUE
  ),
  
  fill = guide_legend(
    nrow = 2,
    byrow = TRUE
  ),
  
  shape = guide_legend(
    nrow = 2,
    byrow = TRUE
  )
  
)


# =============================================================================
# Colours
# =============================================================================

val_colors <- c(
  
  "VAL 1" = "#4DAF4A",
  
  "VAL 2" = "#FFB000",
  
  "VAL 3" = "#F07F2F",
  
  "VAL 4" = "#C43C55"
  
)


activity_colors <- c(
  
  "Background unrest" =
    "#4DAF4A",
  
  "Moderate to heightened unrest" =
    "#FFB000",
  
  "Minor to moderate eruption" =
    "#D73027"
  
)


model_label_colors <- c(
  
  "1_Background_unrest" =
    "#4DAF4A",
  
  "2_Moderate_to_heightened_unrest" =
    "#FFB000",
  
  "3_Minor_to_moderate_eruption" =
    "#D73027"
  
)


truth_shapes <- c(
  
  "1_Background_unrest" = 16,
  
  "2_Moderate_to_heightened_unrest" = 17,
  
  "3_Minor_to_moderate_eruption" = 15
  
)


# =============================================================================
# 2. Read VAL system 3 records
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
# Cap open-ended records
# -----------------------------------------------------------------------------

val_df$PeriodEnd_Local[
  
  is.na(
    val_df$PeriodEnd_Local
  ) |
    
    stringr::str_starts(
      val_df$`PeriodEnd LocalTime`,
      "9999"
    )
  
] <- ANALYSIS_END_LOCAL


# -----------------------------------------------------------------------------
# Keep VAL system version 3
# -----------------------------------------------------------------------------

val_df <- val_df %>%
  
  filter(
    
    `VAL System Version` ==
      VAL_SYSTEM_VERSION,
    
    PeriodEnd_Local >
      ANALYSIS_START_LOCAL,
    
    PeriodStart_Local <
      ANALYSIS_END_LOCAL
    
  ) %>%
  
  mutate(
    
    Activity_3Class = case_when(
      
      Description ==
        "Minor volcanic unrest" ~
        
        "Background unrest",
      
      Description ==
        "Moderate to heightened volcanic unrest" ~
        
        "Moderate to heightened unrest",
      
      Description %in% c(
        
        "Minor volcanic eruption",
        
        "Moderate volcanic eruption"
        
      ) ~
        
        "Minor to moderate eruption",
      
      TRUE ~ NA_character_
      
    ),
    
    Model_Label_Num =
      activity_to_model_level(
        Activity_3Class
      ),
    
    Model_Label = factor(
      
      model_level_to_label(
        Model_Label_Num
      ),
      
      levels =
        model_label_levels
      
    )
    
  )


# =============================================================================
# 2a. Sankey: VAL and description
# =============================================================================

val_system_sankey <- val_df %>%
  
  transmute(
    
    VAL_Level =
      paste(
        "VAL",
        VAL
      ),
    
    Description =
      stringr::str_wrap(
        Description,
        width = 25
      )
    
  ) %>%
  
  count(
    
    VAL_Level,
    
    Description,
    
    name =
      "Number_of_periods"
    
  )


p_val_system_sankey <- ggplot(
  
  val_system_sankey,
  
  aes(
    
    axis1 = VAL_Level,
    
    axis2 = Description,
    
    y = Number_of_periods
    
  )
  
) +
  
  geom_alluvium(
    
    aes(
      fill = VAL_Level
    ),
    
    width = 0.15,
    
    alpha = 0.78
    
  ) +
  
  geom_stratum(
    
    width = 0.15,
    
    fill = "grey96",
    
    color = "grey35"
    
  ) +
  
  geom_text(
    
    stat = "stratum",
    
    aes(
      label =
        after_stat(stratum)
    ),
    
    family = "serif",
    
    size = 2.7,
    
    lineheight = 0.9
    
  ) +
  
  scale_x_discrete(
    
    limits = c(
      
      "Volcanic Alert Level",
      
      "Activity description"
      
    ),
    
    expand =
      c(0.12, 0.12)
    
  ) +
  
  scale_fill_manual(
    values = val_colors
  ) +
  
  labs(
    
    title =
      "VAL Level and Description",
    
    x = NULL,
    
    y =
      "Number of recorded VAL periods",
    
    fill = NULL
    
  ) +
  
  serif_theme +
  
  legend_two_rows


print(
  p_val_system_sankey
)


ggsave(
  
  file.path(
    
    CLUSTER_DIR,
    
    "VAL_System_Description_Sankey.pdf"
    
  ),
  
  p_val_system_sankey,
  
  width = 10,
  
  height = 5.5
  
)


# =============================================================================
# Sankey: VAL and three model classes
# =============================================================================

val_activity_sankey <- val_df %>%
  
  transmute(
    
    VAL_Level =
      paste(
        "VAL",
        VAL
      ),
    
    Activity_3Class
    
  ) %>%
  
  filter(
    !is.na(
      Activity_3Class
    )
  ) %>%
  
  count(
    
    VAL_Level,
    
    Activity_3Class,
    
    name =
      "Number_of_periods"
    
  )


p_val_activity_sankey <- ggplot(
  
  val_activity_sankey,
  
  aes(
    
    axis1 =
      VAL_Level,
    
    axis2 =
      Activity_3Class,
    
    y =
      Number_of_periods
    
  )
  
) +
  
  geom_alluvium(
    
    aes(
      fill =
        VAL_Level
    ),
    
    width = 0.18,
    
    alpha = 0.78
    
  ) +
  
  geom_stratum(
    
    width = 0.18,
    
    fill = "grey96",
    
    color = "grey35"
    
  ) +
  
  geom_text(
    
    stat = "stratum",
    
    aes(
      label =
        after_stat(stratum)
    ),
    
    family = "serif",
    
    size = 3
    
  ) +
  
  scale_x_discrete(
    
    limits = c(
      
      "Volcanic Alert Level",
      
      "Activity class"
      
    ),
    
    expand =
      c(0.15, 0.15)
    
  ) +
  
  scale_fill_manual(
    values =
      val_colors
  ) +
  
  labs(
    
    title =
      "VAL Level and Model Class",
    
    x = NULL,
    
    y =
      "Number of recorded VAL periods",
    
    fill = NULL
    
  ) +
  
  serif_theme +
  
  legend_two_rows


print(
  p_val_activity_sankey
)


ggsave(
  
  file.path(
    
    CLUSTER_DIR,
    
    "VAL_Activity_Class_Sankey.pdf"
    
  ),
  
  p_val_activity_sankey,
  
  width = 9,
  
  height = 5.5
  
)


# =============================================================================
# Assign VAL label to weekly window
# =============================================================================

assign_val_label <- function(
    week_start,
    week_end,
    val_table
) {
  
  overlap <-
    
    val_table$PeriodStart_Local <
    week_end &
    
    val_table$PeriodEnd_Local >
    week_start
  
  
  if (
    !any(
      overlap,
      na.rm = TRUE
    )
  ) {
    
    return(
      
      c(
        
        Activity_3Class =
          NA_character_,
        
        Model_Label_Num =
          NA_character_,
        
        Model_Label =
          NA_character_,
        
        VAL_Levels =
          NA_character_,
        
        VAL_Descriptions =
          NA_character_
        
      )
      
    )
    
  }
  
  
  descriptions <-
    
    unique(
      val_table$Description[
        overlap
      ]
    )
  
  
  val_levels <-
    
    unique(
      val_table$VAL[
        overlap
      ]
    )
  
  
  model_levels <-
    
    unique(
      
      na.omit(
        
        val_table$Model_Label_Num[
          overlap
        ]
        
      )
      
    )
  
  
  model_level <-
    
    if (
      length(
        model_levels
      ) > 0
    ) {
      
      max(
        model_levels
      )
      
    } else {
      
      NA_integer_
      
    }
  
  
  c(
    
    Activity_3Class =
      model_level_to_activity(
        model_level
      ),
    
    Model_Label_Num =
      as.character(
        model_level
      ),
    
    Model_Label =
      model_level_to_label(
        model_level
      ),
    
    VAL_Levels =
      paste(
        val_levels,
        collapse = "; "
      ),
    
    VAL_Descriptions =
      paste(
        descriptions,
        collapse = "; "
      )
    
  )
  
}


# =============================================================================
# Confusion matrix functions
# =============================================================================

metric_values <- function(
    cm,
    metric
) {
  
  if (
    is.null(
      dim(
        cm$byClass
      )
    )
  ) {
    
    return(
      as.numeric(
        cm$byClass[
          metric
        ]
      )
    )
    
  }
  
  
  as.numeric(
    cm$byClass[
      ,
      metric
    ]
  )
  
}


metric_mean <- function(
    cm,
    metric,
    na_to_zero = FALSE
) {
  
  values <-
    metric_values(
      cm,
      metric
    )
  
  
  if (
    na_to_zero
  ) {
    
    values[
      is.na(values)
    ] <- 0
    
  }
  
  
  mean(
    values,
    na.rm = TRUE
  )
  
}


metrics_from_confusion <- function(
    cm
) {
  
  data.frame(
    
    Accuracy =
      unname(
        cm$overall[
          "Accuracy"
        ]
      ),
    
    Kappa =
      unname(
        cm$overall[
          "Kappa"
        ]
      ),
    
    Mean_Balanced_Accuracy =
      metric_mean(
        cm,
        "Balanced Accuracy"
      ),
    
    Mean_F1 =
      metric_mean(
        cm,
        "F1",
        na_to_zero = TRUE
      ),
    
    Mean_Sensitivity =
      metric_mean(
        cm,
        "Sensitivity"
      ),
    
    Mean_Specificity =
      metric_mean(
        cm,
        "Specificity"
      )
    
  )
  
}


confusion_to_df <- function(
    cm
) {
  
  as.data.frame.matrix(
    cm$table
  ) %>%
    
    rownames_to_column(
      "Prediction"
    )
  
}


# =============================================================================
# Correlation heatmap
# =============================================================================

plot_correlation_heatmap <- function(
    d,
    station_name,
    feature_cols,
    method = "spearman"
) {
  
  cor_mat <- cor(
    
    d[
      ,
      feature_cols
    ],
    
    method = method,
    
    use =
      "pairwise.complete.obs"
    
  )
  
  
  cor_df <-
    
    as.data.frame(
      cor_mat
    ) %>%
    
    rownames_to_column(
      "Feature_1"
    ) %>%
    
    pivot_longer(
      
      -Feature_1,
      
      names_to =
        "Feature_2",
      
      values_to =
        "Correlation"
      
    ) %>%
    
    mutate(
      
      Feature_1 =
        factor(
          Feature_1,
          levels = feature_cols
        ),
      
      Feature_2 =
        factor(
          Feature_2,
          levels = feature_cols
        )
      
    )
  
  
  ggplot(
    
    cor_df,
    
    aes(
      
      x = Feature_2,
      
      y = Feature_1,
      
      fill =
        Correlation
      
    )
    
  ) +
    
    geom_tile(
      color = "white"
    ) +
    
    geom_text(
      
      aes(
        label =
          sprintf(
            "%.2f",
            Correlation
          )
      ),
      
      family =
        "serif",
      
      size = 3
      
    ) +
    
    scale_fill_gradient2(
      
      low =
        "#2166AC",
      
      mid =
        "white",
      
      high =
        "#B2182B",
      
      midpoint = 0,
      
      limits =
        c(-1, 1),
      
      name =
        "Spearman\nrho"
      
    ) +
    
    labs(
      
      title =
        paste0(
          station_name,
          ": Feature Correlation"
        )
      
    ) +
    
    theme_bw(
      base_family =
        "serif"
    ) +
    
    theme(
      
      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5,
          size = 11
        ),
      
      axis.title =
        element_blank(),
      
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1
        ),
      
      panel.grid =
        element_blank(),
      
      legend.position =
        "right"
      
    )
  
}


# =============================================================================
# 3. Run analysis separately for each station
# =============================================================================

run_station_cluster_analysis <- function(
    station_name,
    op_file,
    val_df,
    out_dir
) {
  
  
  dir.create(
    
    out_dir,
    
    recursive = TRUE,
    
    showWarnings = FALSE
    
  )
  
  
  cat(
    sprintf(
      "\n\n########## %s ##########\n",
      station_name
    )
  )
  
  
  # ===========================================================================
  # Read weekly ordinal-pattern features
  # ===========================================================================
  
  weekly_df <- read_csv(
    op_file,
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
      
      Week_Start >=
        ANALYSIS_START_LOCAL,
      
      Week_Start <=
        ANALYSIS_END_LOCAL
      
    ) %>%
    
    filter(
      
      if_all(
        all_of(features),
        is.finite
      )
      
    )
  
  
  if (
    
    anyNA(
      weekly_df$Week_Start
    ) ||
    
    anyNA(
      weekly_df$Week_End
    )
    
  ) {
    
    stop(
      
      sprintf(
        
        "[%s] Week_Start or Week_End could not be converted to date-times.",
        
        station_name
        
      )
      
    )
    
  }
  
  
  # ===========================================================================
  # Assign VAL labels
  # ===========================================================================
  
  assigned <- t(
    
    mapply(
      
      assign_val_label,
      
      weekly_df$Week_Start,
      
      weekly_df$Week_End,
      
      MoreArgs =
        list(
          val_table =
            val_df
        )
      
    )
    
  )
  
  
  weekly_df$Activity_3Class <-
    assigned[
      ,
      "Activity_3Class"
    ]
  
  
  weekly_df$Model_Label_Num <-
    as.integer(
      assigned[
        ,
        "Model_Label_Num"
      ]
    )
  
  
  weekly_df$Model_Label <-
    assigned[
      ,
      "Model_Label"
    ]
  
  
  weekly_df$VAL_Levels <-
    assigned[
      ,
      "VAL_Levels"
    ]
  
  
  weekly_df$VAL_Descriptions <-
    assigned[
      ,
      "VAL_Descriptions"
    ]
  
  
  weekly_df <- weekly_df %>%
    
    mutate(
      
      Activity_3Class =
        factor(
          
          Activity_3Class,
          
          levels =
            activity_levels
          
        ),
      
      Model_Label =
        factor(
          
          Model_Label,
          
          levels =
            model_label_levels
          
        )
      
    ) %>%
    
    filter(
      
      !is.na(
        Activity_3Class
      ),
      
      !is.na(
        Model_Label
      ),
      
      !is.na(
        Model_Label_Num
      )
      
    )
  
  
  cat(
    
    sprintf(
      
      "[%s] VAL system 3 model labels by weekly window:\n",
      
      station_name
      
    )
    
  )
  
  
  print(
    
    table(
      
      weekly_df$Model_Label,
      
      useNA =
        "ifany"
      
    )
    
  )
  
  
  cat(
    
    sprintf(
      
      "\n[%s] Activity labels by weekly window:\n",
      
      station_name
      
    )
    
  )
  
  
  print(
    
    table(
      
      weekly_df$Activity_3Class,
      
      useNA =
        "ifany"
      
    )
    
  )
  
  
  cat(
    
    sprintf(
      
      "\n[%s] Weekly data range used for modelling:\n",
      
      station_name
      
    )
    
  )
  
  
  print(
    
    range(
      weekly_df$Week_Start
    )
    
  )
  
  
  # ===========================================================================
  # Feature correlations
  # ===========================================================================
  
  feature_cor <- cor(
    
    weekly_df[
      ,
      features
    ],
    
    use =
      "pairwise.complete.obs"
    
  )
  
  
  feature_cor_pairs <-
    
    purrr::map_dfr(
      
      combn(
        features,
        2,
        simplify = FALSE
      ),
      
      ~ tibble(
        
        Feature_1 =
          .x[1],
        
        Feature_2 =
          .x[2],
        
        Correlation =
          unname(
            
            feature_cor[
              .x[1],
              .x[2]
            ]
            
          ),
        
        Abs_Correlation =
          abs(
            Correlation
          )
        
      )
      
    ) %>%
    
    arrange(
      desc(
        Abs_Correlation
      )
    )
  
  
  high_feature_correlations <-
    
    feature_cor_pairs %>%
    
    filter(
      
      Abs_Correlation >=
        FEATURE_CORRELATION_CUTOFF
      
    )
  
  
  cat(
    
    sprintf(
      
      "\n[%s] Feature pairs with |correlation| >= %.2f:\n",
      
      station_name,
      
      FEATURE_CORRELATION_CUTOFF
      
    )
    
  )
  
  
  print(
    high_feature_correlations
  )
  
  
  # ===========================================================================
  # Spearman correlation heatmap
  # ===========================================================================
  
  p_corr <-
    
    plot_correlation_heatmap(
      
      weekly_df,
      
      station_name,
      
      features,
      
      method = "spearman"
      
    )
  
  
  print(
    p_corr
  )
  
  
  ggsave(
    
    file.path(
      
      out_dir,
      
      sprintf(
        
        "Feature_Correlation_Spearman_%s.pdf",
        
        station_name
        
      )
      
    ),
    
    p_corr,
    
    width = 7,
    
    height = 6
    
  )
  
  
  # ===========================================================================
  # DBSCAN THREE-CLUSTER ANALYSIS
  # ===========================================================================
  
  x_all <- scale(
    
    weekly_df[
      ,
      features
    ]
    
  )
  
  
  # ---------------------------------------------------------------------------
  # Search both eps and minPts for exactly THREE non-noise clusters
  # ---------------------------------------------------------------------------
  
  dbscan_search <- data.frame()
  
  for (m in DBSCAN_MINPTS_GRID) {
    for (e in DBSCAN_EPS_GRID) {
      
      temp_fit <- dbscan(
        x_all,
        eps = e,
        minPts = m
      )
      
      temp_cluster <- temp_fit$cluster
      
      n_clusters <- length(
        unique(temp_cluster[temp_cluster != 0])
      )
      
      noise_n <- sum(temp_cluster == 0)
      noise_prop <- mean(temp_cluster == 0)
      
      dbscan_search <- rbind(
        dbscan_search,
        data.frame(
          eps = e,
          minPts = m,
          n_clusters = n_clusters,
          noise_n = noise_n,
          noise_prop = noise_prop
        )
      )
    }
  }
  
  # Keep only parameter combinations producing exactly three clusters
  three_cluster_search <- dbscan_search %>%
    filter(n_clusters == 3)
  
  if (nrow(three_cluster_search) == 0) {
    stop(
      sprintf(
        paste0(
          "[%s] No combination of eps and minPts in the current search grids ",
          "produced exactly three DBSCAN clusters."
        ),
        station_name
      )
    )
  }
  
  # Choose the three-cluster solution with the smallest noise proportion.
  # Ties are resolved by preferring eps closest to 1.1 and minPts closest
  # to the original rule 2 * number of features.
  original_minPts <- MINPTS_MULTIPLIER * length(features)
  
  best_dbscan <- three_cluster_search %>%
    mutate(
      Distance_From_Original_Eps = abs(eps - EPS_VALUE),
      Distance_From_Original_MinPts = abs(minPts - original_minPts)
    ) %>%
    arrange(
      noise_prop,
      Distance_From_Original_Eps,
      Distance_From_Original_MinPts
    ) %>%
    slice(1)
  
  best_eps <- best_dbscan$eps
  best_minPts <- best_dbscan$minPts
  
  cat(
    sprintf(
      "\n[%s] DBSCAN selected eps = %.2f, minPts = %d\n",
      station_name,
      best_eps,
      best_minPts
    )
  )
  
  cat(
    sprintf(
      "[%s] Number of noise observations = %d (%.2f%%)\n",
      station_name,
      best_dbscan$noise_n,
      100 * best_dbscan$noise_prop
    )
  )
  
  # ---------------------------------------------------------------------------
  # k-NN distance plot
  # ---------------------------------------------------------------------------
  
  pdf(
    
    file.path(
      
      out_dir,
      
      sprintf(
        
        "DBSCAN_kNN_distance_plot_%s.pdf",
        
        station_name
        
      )
      
    ),
    
    width = 7,
    
    height = 5
    
  )
  
  
  kNNdistplot(
    
    x_all,
    
    k =
      best_minPts
    
  )
  
  
  abline(
    
    h = best_eps,
    
    col = "red",
    
    lty = 2
    
  )
  
  
  title(
    
    sprintf(
      
      "%s: k-NN distance plot",
      
      station_name
      
    )
    
  )
  
  
  dev.off()
  
  
  # ---------------------------------------------------------------------------
  # Final three-cluster DBSCAN
  # ---------------------------------------------------------------------------
  
  db_fit <- dbscan(
    
    x_all,
    
    eps =
      best_eps,
    
    minPts =
      best_minPts
    
  )
  
  
  weekly_df$DBSCAN_Cluster <-
    db_fit$cluster
  
  
  weekly_df <- weekly_df %>%
    
    mutate(
      
      DBSCAN_Result = case_when(
        
        DBSCAN_Cluster == 0 ~
          "Noise",
        
        DBSCAN_Cluster == 1 ~
          "Cluster 1",
        
        DBSCAN_Cluster == 2 ~
          "Cluster 2",
        
        DBSCAN_Cluster == 3 ~
          "Cluster 3",
        
        TRUE ~
          "Other"
        
      ),
      
      DBSCAN_Result =
        factor(
          
          DBSCAN_Result,
          
          levels = c(
            
            "Noise",
            
            "Cluster 1",
            
            "Cluster 2",
            
            "Cluster 3"
            
          )
          
        )
      
    )
  
  
  # ---------------------------------------------------------------------------
  # DBSCAN counts
  # ---------------------------------------------------------------------------
  
  cat(
    
    sprintf(
      
      "\n[%s] DBSCAN cluster counts:\n",
      
      station_name
      
    )
    
  )
  
  
  print(
    
    table(
      weekly_df$DBSCAN_Cluster
    )
    
  )
  
  
  # ---------------------------------------------------------------------------
  # Three VAL classes versus DBSCAN clusters
  # ---------------------------------------------------------------------------
  
  dbscan_table <- table(
    
    VAL_Class =
      weekly_df$Model_Label,
    
    DBSCAN_Cluster =
      weekly_df$DBSCAN_Result
    
  )
  
  
  cat(
    
    sprintf(
      
      "\n[%s] Three VAL classes versus DBSCAN clusters:\n",
      
      station_name
      
    )
    
  )
  
  
  print(
    dbscan_table
  )
  
  
  # ---------------------------------------------------------------------------
  # Noise observations
  # ---------------------------------------------------------------------------
  
  dbscan_noise <-
    
    weekly_df %>%
    
    filter(
      DBSCAN_Cluster == 0
    ) %>%
    
    select(
      
      any_of(
        c(
          "Week_Index",
          "Variable"
        )
      ),
      
      Week_Start,
      
      Week_End,
      
      Model_Label_Num,
      
      Model_Label,
      
      Activity_3Class,
      
      VAL_Levels,
      
      VAL_Descriptions,
      
      DBSCAN_Cluster,
      
      DBSCAN_Result,
      
      all_of(features)
      
    )
  
  
  cat(
    
    sprintf(
      
      "\n[%s] DBSCAN noise observations: %d\n",
      
      station_name,
      
      nrow(
        dbscan_noise
      )
      
    )
    
  )
  
  
  # ===========================================================================
  # PCA
  # ===========================================================================
  
  pca_fit <-
    prcomp(
      x_all
    )
  
  
  pca_var_pct <-
    
    (
      pca_fit$sdev^2 /
        sum(
          pca_fit$sdev^2
        )
    ) * 100
  
  
  pc1_axis_label <-
    
    sprintf(
      
      "PC1 (%.1f%%)",
      
      pca_var_pct[1]
      
    )
  
  
  pc2_axis_label <-
    
    sprintf(
      
      "PC2 (%.1f%%)",
      
      pca_var_pct[2]
      
    )
  
  
  cat(
    
    sprintf(
      
      paste0(
        
        "\n[%s] PCA variance explained: ",
        "PC1 = %.1f%%, ",
        "PC2 = %.1f%% ",
        "(cumulative %.1f%%)\n"
        
      ),
      
      station_name,
      
      pca_var_pct[1],
      
      pca_var_pct[2],
      
      pca_var_pct[1] +
        pca_var_pct[2]
      
    )
    
  )
  
  
  # ---------------------------------------------------------------------------
  # Labels for easylabel
  # ---------------------------------------------------------------------------
  
  point_labels <-
    
    if (
      "Week_Index" %in%
      names(
        weekly_df
      )
    ) {
      
      paste0(
        
        "Week ",
        
        weekly_df$Week_Index,
        
        " | ",
        
        format(
          
          weekly_df$Week_Start,
          
          "%Y-%m-%d"
          
        )
        
      )
      
    } else {
      
      paste0(
        
        "Week starting ",
        
        format(
          
          weekly_df$Week_Start,
          
          "%Y-%m-%d"
          
        )
        
      )
      
    }
  
  
  pca_df <- weekly_df %>%
    
    mutate(
      
      PC1 =
        pca_fit$x[
          ,
          1
        ],
      
      PC2 =
        pca_fit$x[
          ,
          2
        ],
      
      Point_Label =
        point_labels
      
    )
  
  
  # ===========================================================================
  # PCA plot: VAL classes and DBSCAN clusters
  # ===========================================================================
  
  p_dbscan <- ggplot(
    
    pca_df,
    
    aes(
      
      PC1,
      
      PC2,
      
      color =
        Activity_3Class,
      
      shape =
        DBSCAN_Result
      
    )
    
  ) +
    
    geom_point(
      
      size = 1.8,
      
      alpha = 0.75
      
    ) +
    
    scale_color_manual(
      
      values =
        activity_colors
      
    ) +
    
    labs(
      
      title =
        paste0(
          station_name,
          ": DBSCAN Three-Cluster Solution"
        ),
      
      x =
        pc1_axis_label,
      
      y =
        pc2_axis_label,
      
      color =
        "VAL class",
      
      shape =
        "DBSCAN result"
      
    ) +
    
    serif_theme +
    
    legend_two_rows
  
  
  print(
    p_dbscan
  )
  
  
  ggsave(
    
    file.path(
      
      out_dir,
      
      sprintf(
        
        "DBSCAN_VAL_PCA_plot_%s.pdf",
        
        station_name
        
      )
      
    ),
    
    p_dbscan,
    
    width = 8.5,
    
    height = 5.5
    
  )
  
  
  # ---------------------------------------------------------------------------
  # Easylabel DBSCAN plot
  # ---------------------------------------------------------------------------
  
  if (
    RUN_EASYLABEL
  ) {
    
    easylabel(
      
      pca_df,
      
      x = "PC1",
      
      y = "PC2",
      
      col =
        "DBSCAN_Result",
      
      labs =
        "Point_Label",
      
      size = 7
      
    )
    
  }
  
  
  # ===========================================================================
  # Sankey: VAL class vs DBSCAN cluster
  # ===========================================================================
  
  sankey_df <-
    
    weekly_df %>%
    
    count(
      
      Activity_3Class,
      
      DBSCAN_Result,
      
      name =
        "Weeks"
      
    )
  
  
  p_sankey <- ggplot(
    
    sankey_df,
    
    aes(
      
      axis1 =
        Activity_3Class,
      
      axis2 =
        DBSCAN_Result,
      
      y =
        Weeks
      
    )
    
  ) +
    
    geom_alluvium(
      
      aes(
        fill =
          Activity_3Class
      ),
      
      width = 0.2,
      
      alpha = 0.75
      
    ) +
    
    geom_stratum(
      
      width = 0.2,
      
      fill = "grey95",
      
      color = "grey30"
      
    ) +
    
    geom_text(
      
      stat = "stratum",
      
      aes(
        label =
          after_stat(stratum)
      ),
      
      family =
        "serif",
      
      size = 3
      
    ) +
    
    scale_x_discrete(
      
      limits = c(
        
        "VAL activity class",
        
        "DBSCAN result"
        
      ),
      
      expand =
        c(0.15, 0.15)
      
    ) +
    
    scale_fill_manual(
      
      values =
        activity_colors
      
    ) +
    
    labs(
      
      title =
        paste0(
          
          station_name,
          
          ": VAL vs. DBSCAN"
          
        ),
      
      x = NULL,
      
      y =
        "Number of weekly windows",
      
      fill =
        "VAL class"
      
    ) +
    
    serif_theme +
    
    legend_two_rows
  
  
  print(
    p_sankey
  )
  
  
  ggsave(
    
    file.path(
      
      out_dir,
      
      sprintf(
        
        "VAL_DBSCAN_Sankey_%s.pdf",
        
        station_name
        
      )
      
    ),
    
    p_sankey,
    
    width = 8,
    
    height = 5.5
    
  )
  
  
  # ===========================================================================
  # Supervised three-class SVM
  # ===========================================================================
  
  all_x <-
    
    weekly_df[
      ,
      features
    ]
  
  
  all_y <- factor(
    
    weekly_df$Model_Label,
    
    levels =
      model_label_levels
    
  )
  
  
  all_counts <-
    table(
      all_y
    )
  
  
  cat(
    
    sprintf(
      
      "\n[%s] Full-dataset class counts for CV:\n",
      
      station_name
      
    )
    
  )
  
  
  print(
    all_counts
  )
  
  
  if (
    any(
      all_counts < 2
    )
  ) {
    
    stop(
      
      sprintf(
        
        paste0(
          
          "[%s] A supervised model cannot be fitted: ",
          "one class has fewer than two observations."
          
        ),
        
        station_name
        
      )
      
    )
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Number of CV folds
  # ---------------------------------------------------------------------------
  
  CV_N_FOLDS <-
    
    min(
      
      N_FOLDS,
      
      as.integer(
        min(
          all_counts
        )
      )
      
    )
  
  
  if (
    CV_N_FOLDS <
    N_FOLDS
  ) {
    
    cat(
      
      sprintf(
        
        paste0(
          
          "\n[%s] Using %d folds instead of %d ",
          "because the smallest class has %d weekly windows.\n"
          
        ),
        
        station_name,
        
        CV_N_FOLDS,
        
        N_FOLDS,
        
        min(
          all_counts
        )
        
      )
      
    )
    
  }
  
  
  # ===========================================================================
  # Class weights
  # ===========================================================================
  
  class_weights <- setNames(
    
    length(
      all_y
    ) /
      
      (
        length(
          model_label_levels
        ) *
          
          as.numeric(
            all_counts
          )
      ),
    
    names(
      all_counts
    )
    
  )
  
  
  weights_for_labels <- function(
    y
  ) {
    
    unname(
      
      class_weights[
        as.character(
          y
        )
      ]
      
    )
    
  }
  
  
  # ===========================================================================
  # Weighted RBF SVM
  # ===========================================================================
  
  svm_cost_grid <- c(
    
    0.1,
    
    1,
    
    10,
    
    100,
    
    1000
    
  )
  
  
  svm_gamma_grid <- c(
    
    0.001,
    
    0.005,
    
    0.01,
    
    0.05,
    
    0.1,
    
    0.5,
    
    1,
    
    5
    
  )
  
  
  svm_grid <- expand.grid(
    
    cost =
      svm_cost_grid,
    
    gamma =
      svm_gamma_grid
    
  )
  
  
  # ===========================================================================
  # SVM evaluation function
  # ===========================================================================
  
  evaluate_svm <- function(
    
    x,
    
    y,
    
    cost,
    
    gamma,
    
    class_weights,
    
    n_folds,
    
    n_repeats
    
  ) {
    
    
    acc_rep <-
      numeric(
        n_repeats
      )
    
    kappa_rep <-
      numeric(
        n_repeats
      )
    
    bal_acc_rep <-
      numeric(
        n_repeats
      )
    
    f1_rep <-
      numeric(
        n_repeats
      )
    
    sens_rep <-
      numeric(
        n_repeats
      )
    
    spec_rep <-
      numeric(
        n_repeats
      )
    
    
    for (
      r in seq_len(
        n_repeats
      )
    ) {
      
      
      set.seed(
        1000 + r
      )
      
      
      folds <- createFolds(
        
        y,
        
        k =
          n_folds,
        
        list = TRUE
        
      )
      
      
      fa <-
        rep(
          NA_real_,
          n_folds
        )
      
      fk <-
        rep(
          NA_real_,
          n_folds
        )
      
      fb <-
        rep(
          NA_real_,
          n_folds
        )
      
      ff <-
        rep(
          NA_real_,
          n_folds
        )
      
      fs <-
        rep(
          NA_real_,
          n_folds
        )
      
      fp <-
        rep(
          NA_real_,
          n_folds
        )
      
      
      for (
        f in seq_len(
          n_folds
        )
      ) {
        
        
        test_idx <-
          folds[[f]]
        
        
        train_data <-
          data.frame(
            
            Model_Label =
              y[
                -test_idx
              ],
            
            x[
              -test_idx,
              ,
              drop = FALSE
            ]
            
          )
        
        
        test_data <-
          
          x[
            test_idx,
            ,
            drop = FALSE
          ]
        
        
        test_truth <-
          y[
            test_idx
          ]
        
        
        fit <- svm(
          
          Model_Label ~ .,
          
          data =
            train_data,
          
          kernel =
            "radial",
          
          cost =
            cost,
          
          gamma =
            gamma,
          
          class.weights =
            class_weights
          
        )
        
        
        pred <- factor(
          
          predict(
            
            fit,
            
            newdata =
              test_data
            
          ),
          
          levels =
            model_label_levels
          
        )
        
        
        cm <- confusionMatrix(
          
          pred,
          
          test_truth,
          
          mode =
            "everything"
          
        )
        
        
        fold_metrics <-
          
          metrics_from_confusion(
            cm
          )
        
        
        fa[f] <-
          fold_metrics$Accuracy
        
        fk[f] <-
          fold_metrics$Kappa
        
        fb[f] <-
          fold_metrics$Mean_Balanced_Accuracy
        
        ff[f] <-
          fold_metrics$Mean_F1
        
        fs[f] <-
          fold_metrics$Mean_Sensitivity
        
        fp[f] <-
          fold_metrics$Mean_Specificity
        
      }
      
      
      acc_rep[r] <-
        mean(
          fa,
          na.rm = TRUE
        )
      
      kappa_rep[r] <-
        mean(
          fk,
          na.rm = TRUE
        )
      
      bal_acc_rep[r] <-
        mean(
          fb,
          na.rm = TRUE
        )
      
      f1_rep[r] <-
        mean(
          ff,
          na.rm = TRUE
        )
      
      sens_rep[r] <-
        mean(
          fs,
          na.rm = TRUE
        )
      
      spec_rep[r] <-
        mean(
          fp,
          na.rm = TRUE
        )
      
    }
    
    
    data.frame(
      
      cost =
        cost,
      
      gamma =
        gamma,
      
      Mean_Accuracy =
        mean(
          acc_rep
        ),
      
      Mean_Kappa =
        mean(
          kappa_rep
        ),
      
      Mean_Balanced_Accuracy =
        mean(
          bal_acc_rep
        ),
      
      SD_Balanced_Accuracy =
        sd(
          bal_acc_rep
        ),
      
      Mean_F1 =
        mean(
          f1_rep
        ),
      
      Mean_Sensitivity =
        mean(
          sens_rep
        ),
      
      Mean_Specificity =
        mean(
          spec_rep
        )
      
    )
    
  }
  
  
  # ===========================================================================
  # Search SVM cost/gamma grid
  # ===========================================================================
  
  cat(
    
    sprintf(
      
      paste0(
        
        "\n[%s] Searching SVM cost/gamma grid ",
        "(%d combinations) via repeated stratified CV...\n"
        
      ),
      
      station_name,
      
      nrow(
        svm_grid
      )
      
    )
    
  )
  
  
  svm_grid_results <- do.call(
    
    rbind,
    
    lapply(
      
      seq_len(
        nrow(
          svm_grid
        )
      ),
      
      function(
    i
      ) {
        
        evaluate_svm(
          
          all_x,
          
          all_y,
          
          svm_grid$cost[i],
          
          svm_grid$gamma[i],
          
          class_weights,
          
          CV_N_FOLDS,
          
          SVM_N_REPEATS
          
        )
        
      }
    
    )
    
  )
  
  
  svm_grid_results <-
    
    svm_grid_results %>%
    
    arrange(
      desc(
        Mean_Balanced_Accuracy
      )
    )
  
  
  best_cost <-
    
    svm_grid_results$cost[
      1
    ]
  
  
  best_gamma <-
    
    svm_grid_results$gamma[
      1
    ]
  
  
  cat(
    
    sprintf(
      
      paste0(
        
        "[%s] Selected SVM cost=%.3f, ",
        "gamma=%.3f ",
        "(mean Balanced Accuracy = %.3f)\n"
        
      ),
      
      station_name,
      
      best_cost,
      
      best_gamma,
      
      svm_grid_results$
        Mean_Balanced_Accuracy[1]
      
    )
    
  )
  
  
  print(
    
    head(
      svm_grid_results,
      10
    )
    
  )
  
  
  # ===========================================================================
  # Out-of-fold SVM predictions
  # ===========================================================================
  
  set.seed(
    123456789,
    kind =
      "Mersenne-Twister"
  )
  
  
  cv_folds_best <- createFolds(
    
    all_y,
    
    k =
      CV_N_FOLDS,
    
    list = TRUE
    
  )
  
  
  svm_oof_pred <- factor(
    
    rep(
      
      NA_character_,
      
      nrow(
        weekly_df
      )
      
    ),
    
    levels =
      model_label_levels
    
  )
  
  
  for (
    f in seq_len(
      CV_N_FOLDS
    )
  ) {
    
    
    test_idx <-
      
      cv_folds_best[[f]]
    
    
    train_data <- data.frame(
      
      Model_Label =
        all_y[
          -test_idx
        ],
      
      all_x[
        -test_idx,
        ,
        drop = FALSE
      ]
      
    )
    
    
    test_data <-
      
      all_x[
        test_idx,
        ,
        drop = FALSE
      ]
    
    
    fit <- svm(
      
      Model_Label ~ .,
      
      data =
        train_data,
      
      kernel =
        "radial",
      
      cost =
        best_cost,
      
      gamma =
        best_gamma,
      
      class.weights =
        class_weights
      
    )
    
    
    svm_oof_pred[
      test_idx
    ] <-
      
      predict(
        
        fit,
        
        newdata =
          test_data
        
      )
    
  }
  
  
  svm_cm <- confusionMatrix(
    
    svm_oof_pred,
    
    all_y,
    
    mode =
      "everything"
    
  )
  
  
  cat(
    
    sprintf(
      
      paste0(
        
        "\n[%s] SVM confusion matrix ",
        "(out-of-fold, single %d-fold CV pass):\n"
        
      ),
      
      station_name,
      
      CV_N_FOLDS
      
    )
    
  )
  
  
  print(
    svm_cm
  )
  
  
  # ===========================================================================
  # PCA display of SVM predictions
  # ===========================================================================
  
  model_pca_df <- pca_df %>%
    
    mutate(
      
      Ground_Truth =
        factor(
          
          all_y,
          
          levels =
            model_label_levels
          
        ),
      
      SVM_Prediction =
        factor(
          
          svm_oof_pred,
          
          levels =
            model_label_levels
          
        )
      
    )
  
  
  p_svm_pca <- ggplot(
    
    model_pca_df,
    
    aes(
      
      PC1,
      
      PC2,
      
      color =
        SVM_Prediction,
      
      shape =
        Ground_Truth
      
    )
    
  ) +
    
    geom_point(
      
      size = 1.8,
      
      alpha = 0.72
      
    ) +
    
    scale_color_manual(
      
      values =
        model_label_colors
      
    ) +
    
    scale_shape_manual(
      
      values =
        truth_shapes
      
    ) +
    
    labs(
      
      title =
        paste0(
          station_name,
          ": SVM"
        ),
      
      x =
        pc1_axis_label,
      
      y =
        pc2_axis_label,
      
      color =
        "SVM prediction",
      
      shape =
        "VAL truth"
      
    ) +
    
    serif_theme +
    
    legend_two_rows
  
  
  print(
    p_svm_pca
  )
  
  
  ggsave(
    
    file.path(
      
      out_dir,
      
      sprintf(
        
        "SVM_OOF_Predictions_PCA_%s.pdf",
        
        station_name
        
      )
      
    ),
    
    p_svm_pca,
    
    width = 8.5,
    
    height = 5.5
    
  )
  
  
  if (
    RUN_EASYLABEL
  ) {
    
    easylabel(
      
      model_pca_df,
      
      x = "PC1",
      
      y = "PC2",
      
      col =
        "SVM_Prediction",
      
      labs =
        "Point_Label",
      
      size = 7
      
    )
    
  }
  
  
  # ===========================================================================
  # SVM performance measures
  #
  # DBSCAN is unsupervised.
  # Therefore, Cluster 1 is NOT automatically VAL class 1,
  # Cluster 2 is NOT automatically VAL class 2, etc.
  #
  # The DBSCAN three-class relationship is reported using
  # the contingency table and Sankey diagram instead.
  # ===========================================================================
  
  comparison <- data.frame(
    
    Model =
      "Weighted RBF SVM (CV)",
    
    Target =
      "Three-class VAL model label",
    
    Accuracy =
      svm_grid_results$
      Mean_Accuracy[1],
    
    Kappa =
      svm_grid_results$
      Mean_Kappa[1],
    
    Mean_Balanced_Accuracy =
      svm_grid_results$
      Mean_Balanced_Accuracy[1],
    
    Mean_Sensitivity =
      svm_grid_results$
      Mean_Sensitivity[1],
    
    Mean_Specificity =
      svm_grid_results$
      Mean_Specificity[1],
    
    Mean_F1 =
      svm_grid_results$
      Mean_F1[1]
    
  )
  
  
  cat(
    
    sprintf(
      
      "\n==== [%s] SVM model performance ====\n",
      
      station_name
      
    )
    
  )
  
  
  print(
    
    comparison,
    
    row.names = FALSE
    
  )
  
  
  # ===========================================================================
  # DBSCAN summary table
  # ===========================================================================
  
  dbscan_summary <- data.frame(
    
    Station =
      station_name,
    
    EPS =
      best_eps,
    
    MinPts =
      best_minPts,
    
    Number_of_Clusters = 3,
    
    Noise_Count =
      sum(
        weekly_df$DBSCAN_Cluster == 0
      ),
    
    Noise_Percentage =
      mean(
        weekly_df$DBSCAN_Cluster == 0
      ) * 100,
    
    Cluster_1_Count =
      sum(
        weekly_df$DBSCAN_Cluster == 1
      ),
    
    Cluster_2_Count =
      sum(
        weekly_df$DBSCAN_Cluster == 2
      ),
    
    Cluster_3_Count =
      sum(
        weekly_df$DBSCAN_Cluster == 3
      )
    
  )
  
  
  cat(
    
    sprintf(
      
      "\n[%s] DBSCAN three-cluster summary:\n",
      
      station_name
      
    )
    
  )
  
  
  print(
    dbscan_summary
  )
  
  
  # ===========================================================================
  # Save Excel outputs
  # ===========================================================================
  
  write_xlsx(
    
    list(
      
      Weekly_Labels_DBSCAN =
        weekly_df,
      
      Model_Comparison =
        comparison,
      
      DBSCAN_3Class_Table =
        as.data.frame(
          dbscan_table
        ),
      
      DBSCAN_Summary =
        dbscan_summary,
      
      DBSCAN_Parameter_Search =
        dbscan_search,
      
      DBSCAN_Noise =
        dbscan_noise,
      
      SVM_Confusion_OOF =
        confusion_to_df(
          svm_cm
        ),
      
      SVM_Grid_Results =
        svm_grid_results,
      
      Feature_Correlation =
        as.data.frame(
          feature_cor
        ) %>%
        
        rownames_to_column(
          "Feature"
        ),
      
      High_Correlation_Pairs =
        high_feature_correlations,
      
      VAL_DBSCAN_Sankey =
        sankey_df
      
    ),
    
    path = file.path(
      
      out_dir,
      
      sprintf(
        
        "OP_VAL_System3_DBSCAN_SVM_results_%s.xlsx",
        
        station_name
        
      )
      
    )
    
  )
  
  
  cat(
    
    sprintf(
      
      "\n[%s] outputs saved in:\n  %s\n",
      
      station_name,
      
      out_dir
      
    )
    
  )
  
  
  invisible(
    
    list(
      
      weekly_df =
        weekly_df,
      
      comparison =
        comparison,
      
      dbscan_table =
        dbscan_table,
      
      dbscan_summary =
        dbscan_summary,
      
      dbscan_search =
        dbscan_search,
      
      svm_cm =
        svm_cm
      
    )
    
  )
  
}


# =============================================================================
# 4. Run for both stations
# =============================================================================

station_results <- list()


for (
  station_name in names(
    station_files
  )
) {
  
  
  out_dir <- file.path(
    
    CLUSTER_DIR,
    
    station_name
    
  )
  
  
  station_results[[station_name]] <-
    
    run_station_cluster_analysis(
      
      station_name =
        station_name,
      
      op_file =
        station_files[[station_name]],
      
      val_df =
        val_df,
      
      out_dir =
        out_dir
      
    )
  
}


cat(
  
  "\nAll stations processed. Results under:\n",
  
  CLUSTER_DIR,
  
  "\n"
  
)


# =============================================================================
# 5. Check number of rows used in analysis
# =============================================================================

for (
  st in names(
    station_files
  )
) {
  
  
  raw <- read_csv(
    
    station_files[
      st
    ],
    
    show_col_types = FALSE
    
  )
  
  
  cat(
    
    sprintf(
      
      "[%s] raw row count in CSV file: %d\n",
      
      st,
      
      nrow(
        raw
      )
      
    )
    
  )
  
  
  n_bad_feature_rows <-
    
    raw %>%
    
    filter(
      
      !if_all(
        
        all_of(
          features
        ),
        
        is.finite
        
      )
      
    ) %>%
    
    nrow()
  
  
  cat(
    
    sprintf(
      
      paste0(
        
        "[%s] rows that would be DROPPED by ",
        "if_all(...is.finite) due to at least ",
        "one NA feature: %d\n"
        
      ),
      
      st,
      
      n_bad_feature_rows
      
    )
    
  )
  
  
  if (
    n_bad_feature_rows > 0
  ) {
    
    
    cat(
      
      sprintf(
        
        "[%s] which weeks, and which feature is NA:\n",
        
        st
        
      )
      
    )
    
    
    raw %>%
      
      filter(
        
        !if_all(
          
          all_of(
            features
          ),
          
          is.finite
          
        )
        
      ) %>%
      
      select(
        
        Week_Start,
        
        Week_End,
        
        all_of(
          features
        )
        
      ) %>%
      
      print(
        n = Inf
      )
    
  }
  
}