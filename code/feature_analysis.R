## ============================================================
## WIZ High vs Normal Activity Classification
## Methods: 1) Logistic Regression / LDA (against real RSAM/Displacement
##             ground truth, built via a data-driven percentile threshold)
##          2) Gaussian Mixture Model (soft clustering)
##          5) DBSCAN (density-based clustering)
## Run separately for Displacement and RSAM -- each uses its OWN raw
## column (displacement_avg_m / rsam_avg) to define its own ground truth.
## ============================================================

library(here)
library(readxl)
library(readr)
library(lubridate)
library(writexl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mclust)     # Method 2: GMM (Gaussian Mixture Models)
library(dbscan)     # Method 5: DBSCAN
library(MASS)       # Method 1: LDA (alternative to logistic regression)
library(patchwork)

set.seed(123)

## ------------------------------------------------------------
## 0. EDIT THIS: your 13 feature columns
## ------------------------------------------------------------
feature_cols <- c(
  "H_Shannon", "C_Shannon",
  "H_Renyi",   "C_Renyi",
  "H_Tsallis", "C_Tsallis",
  "H_Fisher",  "C_Fisher",
  "Disequilibrium",
  "Var_H_Shannon", "Var_C_Shannon",
  "Var_H_Renyi",   "Var_H_Tsallis"
)
stopifnot(length(feature_cols) == 13)

## ------------------------------------------------------------
## 1. Paths
## ------------------------------------------------------------
input_path <- here("results", "WIZ_OrdinalPatterns_D5_Monthly_Quarterly.xlsx")
out_dir    <- here("results", "Stepwise_RSAM_Yearly_Monthly_Clustering")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## --- Raw annual eruption/amplitude files ---
## Columns confirmed: day, unix_timestamp, datetime_nz, displacement_avg_m, rsam_avg
eruption_dir         <- here("data")
eruption_file_pattern <- "^WIZ_NZ_[0-9]{4}\\.csv$"
eruption_date_col    <- "unix_timestamp"   # numeric epoch seconds -- reliable, unlike the datetime_nz string column

## Each variable defines "eruptive" from its OWN raw column, since you
## treat Displacement and RSAM separately throughout.
raw_value_col <- c(Displacement = "displacement_avg_m", RSAM = "rsam_avg")

## EDIT ME: percentile used as the data-driven eruptive threshold.
## Higher = stricter ("only the most extreme values count as eruptive").
## Check the diagnostic plot (Section 3b) and adjust if it doesn't match
## what you see visually in the monthly time series.
eruption_percentile <- 0.95

## ------------------------------------------------------------
## 2. Load ordinal pattern data
## ------------------------------------------------------------
df_all <- read_xlsx(input_path, sheet = 1) %>%
  filter(
    Period_Type == "Monthly",
    Year >= 2007,
    Year <= 2022
  ) %>%
  mutate(
    Period = factor(Period, levels = month.abb),
    Year   = factor(Year)
  )

## ------------------------------------------------------------
## 3a. Build ground-truth eruption labels per Year-Period (Monthly)
##     for ONE raw value column, using a percentile threshold.
## ------------------------------------------------------------
build_eruption_groundtruth <- function(dir_path, pattern, date_col,
                                       value_col, percentile = 0.95) {
  
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    message("No files matched '", pattern, "' in ", dir_path, ". Skipping ground truth for ", value_col)
    return(NULL)
  }
  
  raw_all <- lapply(files, function(f) {
    d <- read_csv(f, show_col_types = FALSE)
    ## Same approach as your seismic_all pipeline: numeric epoch -> UTC -> NZ local time
    utc_time <- as.POSIXct(d[[date_col]], origin = "1970-01-01", tz = "UTC")
    d$datetime <- with_tz(utc_time, tzone = "Pacific/Auckland")
    d %>%
      transmute(
        datetime = datetime,
        Year     = factor(format(datetime, "%Y")),
        Period   = factor(format(datetime, "%b"), levels = month.abb),
        value    = .data[[value_col]]
      )
  }) %>% bind_rows()
  
  threshold <- quantile(raw_all$value, probs = percentile, na.rm = TRUE)
  raw_all$is_eruptive_point <- as.integer(raw_all$value >= threshold)
  
  monthly_summary <- raw_all %>%
    group_by(Year, Period) %>%
    summarise(
      n_points      = n(),
      n_eruptive    = sum(is_eruptive_point, na.rm = TRUE),
      eruptive_frac = n_eruptive / n_points,
      known_high    = factor(if_else(n_eruptive > 0, "High", "Normal"),
                             levels = c("Normal", "High")),
      .groups = "drop"
    )
  
  list(monthly = monthly_summary, raw = raw_all, threshold = threshold)
}

## ------------------------------------------------------------
## 3b. Diagnostic plot: raw time series + threshold line, so you can
##     visually confirm the percentile matches what you'd flag by eye.
## ------------------------------------------------------------
plot_groundtruth_diagnostic <- function(gt, variable_name, value_col_name) {
  p <- ggplot(gt$raw, aes(x = datetime, y = value)) +
    geom_line(linewidth = 0.2, color = "grey40") +
    geom_point(data = filter(gt$raw, is_eruptive_point == 1),
               aes(x = datetime, y = value), color = "red", size = 0.6) +
    geom_hline(yintercept = gt$threshold, linetype = "dashed", color = "red") +
    labs(title = paste0(variable_name, ": ", value_col_name,
                        " with eruptive threshold (", scales::percent(eruption_percentile),
                        " percentile = ", round(gt$threshold, 4), ")"),
         x = "Date", y = value_col_name) +
    theme_minimal(base_size = 11)
  print(p)
  ggsave(file.path(out_dir, paste0(variable_name, "_groundtruth_diagnostic.pdf")),
         p, width = 12, height = 4, dpi = 300)
  p
}

## ------------------------------------------------------------
## 4. Helper: print THEN save a plot
## ------------------------------------------------------------
print_save <- function(p, filename, width = 6, height = 5) {
  print(p)
  ggsave(file.path(out_dir, filename), p, width = width, height = height, dpi = 300)
}

## ------------------------------------------------------------
## 5. Method 2 -- Gaussian Mixture Model (soft clustering)
## ------------------------------------------------------------
run_gmm <- function(df_var, feature_cols, variable_name, pca_df) {
  
  X_scaled <- scale(df_var %>% dplyr::select(all_of(feature_cols)))
  gmm <- Mclust(X_scaled, G = 2)
  
  df_var$cluster_gmm <- factor(gmm$classification)
  probs <- gmm$z
  high_col <- which.max(colMeans(probs))
  df_var$prob_high_gmm <- probs[, high_col]
  df_var$label_gmm <- factor(if_else(df_var$prob_high_gmm > 0.5, "High", "Normal"),
                             levels = c("Normal", "High"))
  
  pca_df$label_gmm <- df_var$label_gmm
  pca_df$prob_high_gmm <- df_var$prob_high_gmm
  
  p_gmm <- ggplot(pca_df, aes(PC1, PC2, color = label_gmm, alpha = prob_high_gmm)) +
    geom_point(size = 2) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    labs(title = paste0(variable_name, ": GMM soft clustering"),
         color = "GMM label", alpha = "P(High)") +
    theme_minimal(base_size = 11)
  
  print_save(p_gmm, paste0(variable_name, "_GMM_pca.pdf"))
  
  list(df = df_var, model = gmm, plot = p_gmm)
}

## ------------------------------------------------------------
## 6. Method 5 -- DBSCAN (density-based)
## ------------------------------------------------------------
run_dbscan <- function(df_var, feature_cols, variable_name, pca_df, minPts = NULL) {
  
  X_scaled <- scale(df_var %>% dplyr::select(all_of(feature_cols)))
  if (is.null(minPts)) minPts <- ncol(X_scaled) + 1
  
  kNNdistplot(X_scaled, k = minPts)
  title(main = paste0(variable_name, ": k-NN distance plot (choose eps at the knee)"),
        xlab = "Points sorted by distance", ylab = paste0(minPts, "-NN distance"))
  
  ## Also save this base-R plot to PDF (ggsave only works on ggplot objects)
  pdf(file.path(out_dir, paste0(variable_name, "_kNNdistplot.pdf")), width = 6, height = 5)
  kNNdistplot(X_scaled, k = minPts)
  title(main = paste0(variable_name, ": k-NN distance plot (choose eps at the knee)"),
        xlab = "Points sorted by distance", ylab = paste0(minPts, "-NN distance"))
  dev.off()
  
  dist_k <- sort(kNNdist(X_scaled, k = minPts))
  d1 <- diff(dist_k)
  eps_auto <- dist_k[which.max(d1)]
  message(variable_name, ": auto-selected eps = ", round(eps_auto, 3),
          " (inspect k-NN plot and adjust if needed)")
  
  db <- dbscan(X_scaled, eps = eps_auto, minPts = minPts)
  df_var$cluster_dbscan <- factor(db$cluster)
  
  clust_sizes <- table(db$cluster[db$cluster != 0])
  if (length(clust_sizes) >= 1) {
    high_clust <- names(clust_sizes)[which.min(clust_sizes)]
    df_var$label_dbscan <- factor(
      case_when(
        db$cluster == 0 ~ "Noise",
        db$cluster == as.integer(high_clust) ~ "High",
        TRUE ~ "Normal"
      ),
      levels = c("Normal", "High", "Noise")
    )
  } else {
    df_var$label_dbscan <- factor("Noise", levels = c("Normal", "High", "Noise"))
  }
  
  pca_df$label_dbscan <- df_var$label_dbscan
  
  p_db <- ggplot(pca_df, aes(PC1, PC2, color = label_dbscan)) +
    geom_point(size = 2) +
    labs(title = paste0(variable_name, ": DBSCAN (eps=", round(eps_auto, 2),
                        ", minPts=", minPts, ")"),
         color = "DBSCAN label") +
    theme_minimal(base_size = 11)
  
  print_save(p_db, paste0(variable_name, "_DBSCAN_pca.pdf"))
  
  list(df = df_var, model = db, eps = eps_auto, minPts = minPts, plot = p_db)
}

## ------------------------------------------------------------
## 7. Method 1 -- Logistic Regression / LDA against real ground truth
## ------------------------------------------------------------
run_supervised <- function(df_var, feature_cols, variable_name, pca_df) {
  
  if (!"known_high" %in% names(df_var) || all(is.na(df_var$known_high))) {
    message(variable_name, ": no ground-truth labels available -- skipping Method 1.")
    return(NULL)
  }
  
  X <- df_var %>% dplyr::select(all_of(feature_cols))
  X_scaled <- as.data.frame(scale(X))
  X_scaled$known_high <- df_var$known_high
  
  logit_fit <- glm(known_high ~ ., data = X_scaled, family = binomial)
  print(summary(logit_fit))
  
  df_var$prob_high_logit <- predict(logit_fit, type = "response")
  df_var$label_logit <- factor(if_else(df_var$prob_high_logit > 0.5, "High", "Normal"),
                               levels = c("Normal", "High"))
  
  lda_fit <- lda(known_high ~ ., data = X_scaled)
  lda_pred <- predict(lda_fit)
  df_var$label_lda <- lda_pred$class
  
  conf_logit <- table(Truth = df_var$known_high, Predicted = df_var$label_logit)
  conf_lda   <- table(Truth = df_var$known_high, Predicted = df_var$label_lda)
  message(variable_name, " logistic regression accuracy: ",
          round(100 * sum(diag(conf_logit)) / sum(conf_logit), 1), "%")
  message(variable_name, " LDA accuracy: ",
          round(100 * sum(diag(conf_lda)) / sum(conf_lda), 1), "%")
  
  pca_df$known_high      <- df_var$known_high
  pca_df$label_logit     <- df_var$label_logit
  pca_df$prob_high_logit <- df_var$prob_high_logit
  
  p_truth <- ggplot(pca_df, aes(PC1, PC2, color = known_high)) +
    geom_point(size = 2) +
    labs(title = paste0(variable_name, ": Known eruption ground truth"), color = "Truth") +
    theme_minimal(base_size = 11)
  
  p_logit <- ggplot(pca_df, aes(PC1, PC2, color = label_logit, alpha = prob_high_logit)) +
    geom_point(size = 2) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    labs(title = paste0(variable_name, ": Logistic regression prediction"),
         color = "Predicted", alpha = "P(High)") +
    theme_minimal(base_size = 11)
  
  print_save(p_truth, paste0(variable_name, "_ground_truth_pca.pdf"))
  print_save(p_logit, paste0(variable_name, "_logit_pca.pdf"))
  
  list(df = df_var, logit_model = logit_fit, lda_model = lda_fit,
       conf_logit = conf_logit, conf_lda = conf_lda,
       plots = list(truth = p_truth, logit = p_logit))
}

## ------------------------------------------------------------
## 8. Master pipeline per variable
## ------------------------------------------------------------
run_all_methods <- function(data, variable_name, feature_cols, out_dir,
                            raw_value_col, percentile) {
  
  message("\n===== ", variable_name, " =====")
  
  ## Build this variable's own ground truth from its own raw column
  gt <- build_eruption_groundtruth(eruption_dir, eruption_file_pattern, eruption_date_col,
                                   value_col = raw_value_col, percentile = percentile)
  
  df_var <- data %>%
    filter(Variable == variable_name) %>%
    drop_na(all_of(feature_cols))
  
  if (nrow(df_var) < 10) {
    warning("Too few rows for ", variable_name, "; skipping.")
    return(NULL)
  }
  
  if (!is.null(gt)) {
    plot_groundtruth_diagnostic(gt, variable_name, raw_value_col)
    df_var <- df_var %>% left_join(gt$monthly, by = c("Year", "Period"))
  }
  
  X_scaled <- scale(df_var %>% dplyr::select(all_of(feature_cols)))
  pca <- prcomp(X_scaled)
  pca_df <- as.data.frame(pca$x[, 1:2])
  var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)
  pca_df$Year   <- df_var$Year
  pca_df$Period <- df_var$Period
  
  gmm_res    <- run_gmm(df_var, feature_cols, variable_name, pca_df)
  df_var     <- gmm_res$df
  dbscan_res <- run_dbscan(df_var, feature_cols, variable_name, pca_df)
  df_var     <- dbscan_res$df
  supervised_res <- run_supervised(df_var, feature_cols, variable_name, pca_df)
  if (!is.null(supervised_res)) df_var <- supervised_res$df
  
  plots <- list(gmm_res$plot, dbscan_res$plot)
  if (!is.null(supervised_res)) plots <- c(supervised_res$plots, plots)
  combo <- wrap_plots(plots, ncol = 2)
  print_save(combo, paste0(variable_name, "_all_methods_combined.pdf"), width = 11, height = 8)
  
  write_xlsx(df_var, path = file.path(out_dir, paste0(variable_name, "_all_methods_results.xlsx")))
  
  list(df = df_var, gmm = gmm_res, dbscan = dbscan_res, supervised = supervised_res,
       ground_truth_threshold = if (!is.null(gt)) gt$threshold else NA,
       pca_var_explained = var_explained)
}

## ------------------------------------------------------------
## 9. Run for each variable separately (each with its OWN raw column)
## ------------------------------------------------------------
results_displacement <- run_all_methods(df_all, "Displacement", feature_cols, out_dir,
                                        raw_value_col = raw_value_col["Displacement"],
                                        percentile = eruption_percentile)

results_rsam <- run_all_methods(df_all, "RSAM", feature_cols, out_dir,
                                raw_value_col = raw_value_col["RSAM"],
                                percentile = eruption_percentile)

