# =====================================================================
# Event coincidence analysis (ECA) - weekly - grid search - ALL features
#   coincidences, circular-shift test, grid-corrected p,
#   Benjamini-Hochberg q-values across features
# All times in UTC.
# =====================================================================

library(readr)
library(lubridate)

# ---------------------------------------------------------------------
# Settings - change only these lines
# ---------------------------------------------------------------------
data_dir    <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"
results_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"

station   <- "WIZ"                     # "WIZ" or "WSRZ"
direction <- "upper"                   # "upper" = unusually high, "lower" = unusually low

all_features <- c(
  "H_Shannon", "C_Shannon", "Fisher_Info", "C_Fisher",
  "H_Renyi", "C_Renyi", "H_Tsallis", "C_Tsallis",
  "Var_H_Shannon", "Var_C_Shannon", "Var_H_Renyi", "Var_H_Tsallis",
  "Disequilibrium"
)

thr_grid <- c(0.80, 0.85, 0.90, 0.95)  # thresholds (percentiles)
dT_grid  <- c(1, 2, 3)                 # length of the look-back window (weeks)
tau_grid <- c(0, 1)                    # extra weeks skipped before the window

# eruption onsets (UTC, from the VAL table)
onset <- ymd_hms(c("2016-04-27 23:30:00",
                   "2016-09-13 03:15:00",
                   "2019-12-09 01:30:00"), tz = "UTC")
#onset <- ymd_hms(c("2012-08-06 05:30:00", 
#                   "2012-12-12 00:30:00", 
#                   "2013-02-22 23:00:00", 
#                   "2013-08-19 23:10:00",
#                  "2013-10-11 23:00:00",
#                  "2016-04-27 23:30:00",
#                  "2016-09-13 03:15:00",
#                  "2019-12-09 01:30:00"), tz = "UTC")

out_dir <- file.path(results_dir, paste0("ECA_results_", station, "_", direction))
dir.create(out_dir, showWarnings = FALSE)

# ---------------------------------------------------------------------
# Step 1: read the weekly results and number the weeks
# ---------------------------------------------------------------------
file_name <- paste0("HC_RSAM_", station, "_20140701_to_20191231_clean.csv")
#file_name <- paste0("HC_RSAM_", station, "_6daysoverlap_20070101_to_20191231_clean.csv")
feat <- read_csv(file.path(data_dir, file_name), show_col_types = FALSE)

week_end <- feat$Week_End
week     <- round(as.numeric(difftime(feat$Week_Start, min(feat$Week_Start),
                                      units = "days")) / 6) + 1   # weeks start every 6 days
n_weeks  <- max(week)
n_shift  <- n_weeks - 1                                           # all possible shifts

cat(station, "| rows:", nrow(feat), "| weeks in period:", n_weeks, "\n")

# ---------------------------------------------------------------------
# Step 2: eruption weeks (first week that ends after the onset)
# ---------------------------------------------------------------------
erupt_week <- rep(NA, length(onset))
for (j in 1:length(onset)) {
  erupt_week[j] <- week[which(week_end > onset[j])[1]]
}
print(data.frame(onset = onset, erupt_week = erupt_week))
N_B <- length(erupt_week)

# all combinations of threshold, dT and tau
grid  <- expand.grid(thr = thr_grid, dT = dT_grid, tau = tau_grid)
n_set <- nrow(grid)

# function: how many eruptions have a feature event in their window
count_hits <- function(ev_week, lo, hi) {
  hits <- 0
  for (j in 1:N_B) {
    gap <- erupt_week[j] - ev_week
    if (any(gap >= lo & gap <= hi)) hits <- hits + 1
  }
  hits
}

# ---------------------------------------------------------------------
# Step 3: loop over all features
# ---------------------------------------------------------------------
summary_all <- NULL

for (f in all_features) {
  
  cat("\n=====", f, "=====\n")
  x <- feat[[f]]
  
  settings <- grid
  settings$n_events <- NA
  settings$hits_B   <- NA     # eruptions with a feature event before them
  settings$hits_A   <- NA     # feature events followed by an eruption
  null_B  <- matrix(NA, n_shift, n_set)
  ev_list <- list()
  
  # -------------------------------------------------------------------
  # 3a: real coincidences and shifted (null) coincidences
  # -------------------------------------------------------------------
  for (k in 1:n_set) {
    
    # weeks above (or below) the threshold
    if (direction == "upper") {
      ev <- !is.na(x) & x > quantile(x, settings$thr[k], na.rm = TRUE)
    } else {
      ev <- !is.na(x) & x < quantile(x, 1 - settings$thr[k], na.rm = TRUE)
    }
    ev_week <- week[ev & !c(FALSE, head(ev, -1))]    # keep first week of each run
    ev_list[[k]] <- ev_week
    
    # window: lo ... hi weeks before the eruption week
    lo <- 1 + settings$tau[k]
    hi <- settings$dT[k] + settings$tau[k]
    
    # precursor hits (real)
    settings$hits_B[k] <- count_hits(ev_week, lo, hi)
    
    # trigger hits (real)
    kA <- 0
    for (i in seq_along(ev_week)) {
      gap <- erupt_week - ev_week[i]
      if (any(gap >= lo & gap <= hi)) kA <- kA + 1
    }
    settings$hits_A[k]   <- kA
    settings$n_events[k] <- length(ev_week)
    
    # precursor hits after circular shift by s weeks
    for (s in 1:n_shift) {
      ev_s <- (ev_week + s - 1) %% n_weeks + 1
      null_B[s, k] <- count_hits(ev_s, lo, hi)
    }
  }
  
  # -------------------------------------------------------------------
  # 3b: p-value of every setting (circular-shift test)
  # -------------------------------------------------------------------
  settings$p <- NA
  for (k in 1:n_set) {
    settings$p[k] <- (1 + sum(null_B[, k] >= settings$hits_B[k])) / (n_shift + 1)
  }
  settings$precursor_rate <- round(settings$hits_B / N_B, 2)
  settings$trigger_rate   <- round(settings$hits_A / settings$n_events, 3)
  
  # -------------------------------------------------------------------
  # 3c: p-value corrected for picking the best of the 24 settings
  # -------------------------------------------------------------------
  p_null <- matrix(NA, n_shift, n_set)
  for (k in 1:n_set) {
    for (s in 1:n_shift) {
      p_null[s, k] <- (sum(null_B[, k] >= null_B[s, k]) +
                         (settings$hits_B[k] >= null_B[s, k])) / (n_shift + 1)
    }
  }
  best_null <- apply(p_null, 1, min)          # best p of each shifted series
  b         <- which.min(settings$p)          # best real setting
  p_grid    <- (1 + sum(best_null <= settings$p[b])) / (n_shift + 1)
  
  cat("best: thr", settings$thr[b], "| dT", settings$dT[b], "| tau", settings$tau[b],
      "| hits", settings$hits_B[b], "of", N_B,
      "| p =", round(settings$p[b], 4), "| corrected =", round(p_grid, 4), "\n")
  
  # -------------------------------------------------------------------
  # 3d: save the grid table of this feature
  # -------------------------------------------------------------------
  settings$p <- round(settings$p, 4)
  write_csv(settings, file.path(out_dir, paste0("ECA_grid_", station, "_", f, ".csv")))
  
  # -------------------------------------------------------------------
  # 3e: figure of this feature (threshold and events of the best setting)
  # -------------------------------------------------------------------
  thr_q    <- ifelse(direction == "upper", settings$thr[b], 1 - settings$thr[b])
  thr_best <- quantile(x, thr_q, na.rm = TRUE)
  ev_best  <- ev_list[[b]]
  ok       <- !is.na(x)
  
  pdf(file.path(out_dir, paste0("ECA_", station, "_", f, ".pdf")), width = 10, height = 4)
  par(mar = c(4, 4, 3, 1))
  plot(week[ok], x[ok], type = "l", xlab = "week number", ylab = f,
       xlim = range(week),
       main = paste0(station, " - ", f, "  (threshold ", settings$thr[b],
                     ", dT = ", settings$dT[b], ", tau = ", settings$tau[b],
                     ", hits ", settings$hits_B[b], "/", N_B,
                     ", corrected p = ", round(p_grid, 3), ")"))
  abline(h = thr_best, col = "red", lty = 2)
  points(ev_best, x[match(ev_best, week)], col = "red", pch = 16)
  abline(v = erupt_week, col = "blue", lwd = 2)
  dev.off()
  
  # -------------------------------------------------------------------
  # 3f: one summary row for this feature
  # -------------------------------------------------------------------
  summary_all <- rbind(summary_all, data.frame(
    feature        = f,
    thr            = settings$thr[b],
    dT             = settings$dT[b],
    tau            = settings$tau[b],
    n_events       = settings$n_events[b],
    hits_B         = settings$hits_B[b],
    hits_A         = settings$hits_A[b],
    precursor_rate = settings$precursor_rate[b],
    trigger_rate   = settings$trigger_rate[b],
    best_p         = round(settings$p[b], 4),
    p_grid         = round(p_grid, 4)))
}

# ---------------------------------------------------------------------
# Step 4: compare the features - Benjamini-Hochberg over all features
# ---------------------------------------------------------------------
summary_all$q_BH <- round(p.adjust(summary_all$p_grid, method = "BH"), 4)
summary_all <- summary_all[order(summary_all$p_grid), ]
print(summary_all)

write_csv(summary_all, file.path(out_dir, paste0("ECA_summary_", station, "_all_features.csv")))

# ---------------------------------------------------------------------
# Step 5: summary figure - corrected p-value of every feature
# ---------------------------------------------------------------------
draw_summary <- function() {
  par(mar = c(4, 9, 3, 1))
  barplot(rev(-log10(summary_all$p_grid)), names.arg = rev(summary_all$feature),
          horiz = TRUE, las = 1, col = ifelse(rev(summary_all$q_BH) < 0.05, "red", "grey60"),
          xlab = "-log10(corrected p)",
          main = paste(station, "- ECA, best setting per feature"))
  abline(v = -log10(0.05), lty = 2)
}
draw_summary()
pdf(file.path(out_dir, paste0("ECA_summary_", station, "_all_features.pdf")), width = 8, height = 6)
draw_summary()
dev.off()

cat("\nAll tables and figures saved in:", out_dir, "\n")

#End of the code
##############################################################
#################################################################
# =====================================================================
# Event coincidence analysis (ECA) - weekly - grid search - ALL features
#   coincidences, circular-shift test, grid-corrected p,
#   Benjamini-Hochberg q-values across features
#   + all grid thresholds drawn on every figure
#   + combined multi-panel figure of all features
# All times in UTC.
# =====================================================================

library(readr)
library(lubridate)

# ---------------------------------------------------------------------
# Settings - change only these lines
# ---------------------------------------------------------------------
data_dir    <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"
results_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"

station   <- "WIZ"                     # "WIZ" or "WSRZ"
direction <- "upper"                   # default: "upper" = unusually high, "lower" = unusually low
lower_features <- c("H_Shannon", "H_Renyi", "H_Tsallis")   # these use "lower" (low entropy = unusual)

tick_step <- 25                        # x-axis tick every tick_step weeks
date_fmt  <- "%Y-%m-%d"                # format of the week-end date labels

all_features <- c(
  "H_Shannon", "C_Shannon", "Fisher_Info", "C_Fisher",
  "H_Renyi", "C_Renyi", "H_Tsallis", "C_Tsallis",
  "Var_H_Shannon", "Var_C_Shannon", "Var_H_Renyi", "Var_H_Tsallis",
  "Disequilibrium"
)

# direction of every feature
feat_dir <- setNames(ifelse(all_features %in% lower_features, "lower", direction), all_features)
dir_note <- if (length(lower_features) > 0) {
  paste0(paste(lower_features, collapse = ", "), ": lower; others: ", direction)
} else direction

thr_grid <- c(0.80, 0.85, 0.90, 0.95)  # thresholds (percentiles)
dT_grid  <- c(1, 2, 3)                 # length of the look-back window (weeks)
tau_grid <- c(0, 1)                    # extra weeks skipped before the window

# one colour per threshold (colour-blind safe); recycled if the grid is longer
thr_cols <- rep(c("#E69F00", "#009E73", "#0072B2", "#CC79A7", "#D55E00"),
                length.out = length(thr_grid))

# eruption onsets (UTC, from the VAL table)
onset <- ymd_hms(c("2016-04-27 23:30:00",
                   "2016-09-13 03:15:00",
                   "2019-12-09 01:30:00"), tz = "UTC")
#onset <- ymd_hms(c("2012-08-06 05:30:00", 
#                   "2012-12-12 00:30:00", 
#                   "2013-02-22 23:00:00", 
#                   "2013-08-19 23:10:00",
#                   "2013-10-11 23:00:00",
#                   "2016-04-27 23:30:00",
#                   "2016-09-13 03:15:00",
#                   "2019-12-09 01:30:00"), tz = "UTC")

dir_tag <- if (length(lower_features) > 0) "mixed" else direction
out_dir <- file.path(results_dir, paste0("ECA_results_", station, "_", dir_tag))
dir.create(out_dir, showWarnings = FALSE)

# ---------------------------------------------------------------------
# Step 1: read the weekly results and number the weeks
# ---------------------------------------------------------------------
file_name <- paste0("HC_RSAM_", station, "_20140701_to_20191231_clean.csv")
#file_name <- paste0("HC_RSAM_", station, "_6daysoverlap_20070101_to_20191231_clean.csv")
feat <- read_csv(file.path(data_dir, file_name), show_col_types = FALSE)

week_end <- feat$Week_End
week     <- round(as.numeric(difftime(feat$Week_Start, min(feat$Week_Start),
                                      units = "days")) / 6) + 1   # weeks start every 6 days
n_weeks  <- max(week)
n_shift  <- n_weeks - 1                                           # all possible shifts

cat(station, "| rows:", nrow(feat), "| weeks in period:", n_weeks, "\n")

# ---------------------------------------------------------------------
# Step 2: eruption weeks (first week that ends after the onset)
# ---------------------------------------------------------------------
erupt_week <- rep(NA, length(onset))
for (j in 1:length(onset)) {
  erupt_week[j] <- week[which(week_end > onset[j])[1]]
}
print(data.frame(onset = onset, erupt_week = erupt_week))
N_B <- length(erupt_week)

# all combinations of threshold, dT and tau
grid  <- expand.grid(thr = thr_grid, dT = dT_grid, tau = tau_grid)
n_set <- nrow(grid)

# function: how many eruptions have a feature event in their window
count_hits <- function(ev_week, lo, hi) {
  hits <- 0
  for (j in 1:N_B) {
    gap <- erupt_week[j] - ev_week
    if (any(gap >= lo & gap <= hi)) hits <- hits + 1
  }
  hits
}

# function: draw every threshold of the grid as a dashed line,
#           the best one thicker and its label in bold
draw_thresholds <- function(x, best_thr, cex_lab = 0.75, dir = direction) {
  for (i in seq_along(thr_grid)) {
    q       <- if (dir == "upper") thr_grid[i] else 1 - thr_grid[i]
    hv      <- quantile(x, q, na.rm = TRUE)
    is_best <- isTRUE(all.equal(thr_grid[i], best_thr))
    abline(h = hv, col = thr_cols[i], lty = 2, lwd = ifelse(is_best, 2.2, 1))
    text(par("usr")[2], hv, paste0("P", round(q * 100)), pos = 4, xpd = NA,
         col = thr_cols[i], cex = cex_lab, font = ifelse(is_best, 2, 1))
  }
}

# week-end date of any week number (interpolated if that week has no row)
week_end_num <- as.numeric(as.POSIXct(week_end, tz = "UTC"))
week_date <- function(wk) {
  as.POSIXct(approx(week, week_end_num, xout = wk, rule = 2, ties = mean)$y,
             origin = "1970-01-01", tz = "UTC")
}
tick_wk <- seq(tick_step, max(week), by = tick_step)

# function: x-axis with a tick every tick_step weeks, labelled by week-end date
draw_week_axis <- function(labels = TRUE, cex_lab = 0.75) {
  axis(1, at = tick_wk, labels = FALSE)
  if (labels) {
    y_in <- grconvertY(par("usr")[3], "user", "inches") - 0.12
    text(tick_wk, grconvertY(y_in, "inches", "user"),
         format(week_date(tick_wk), date_fmt),
         srt = 45, adj = c(1, 1), xpd = NA, cex = cex_lab)
  }
}

# ---------------------------------------------------------------------
# Step 3: loop over all features
# ---------------------------------------------------------------------
summary_all <- NULL
panel_data  <- list()     # kept for the combined figure (Step 6)

for (f in all_features) {
  
  dir_f <- feat_dir[[f]]
  cat("\n=====", f, "(", dir_f, ") =====\n")
  x <- feat[[f]]
  
  settings <- grid
  settings$n_events <- NA
  settings$hits_B   <- NA     # eruptions with a feature event before them
  settings$hits_A   <- NA     # feature events followed by an eruption
  null_B  <- matrix(NA, n_shift, n_set)
  ev_list <- list()
  
  # -------------------------------------------------------------------
  # 3a: real coincidences and shifted (null) coincidences
  # -------------------------------------------------------------------
  for (k in 1:n_set) {
    
    # weeks above (or below) the threshold
    if (dir_f == "upper") {
      ev <- !is.na(x) & x > quantile(x, settings$thr[k], na.rm = TRUE)
    } else {
      ev <- !is.na(x) & x < quantile(x, 1 - settings$thr[k], na.rm = TRUE)
    }
    ev_week <- week[ev & !c(FALSE, head(ev, -1))]    # keep first week of each run
    ev_list[[k]] <- ev_week
    
    # window: lo ... hi weeks before the eruption week
    lo <- 1 + settings$tau[k]
    hi <- settings$dT[k] + settings$tau[k]
    
    # precursor hits (real)
    settings$hits_B[k] <- count_hits(ev_week, lo, hi)
    
    # trigger hits (real)
    kA <- 0
    for (i in seq_along(ev_week)) {
      gap <- erupt_week - ev_week[i]
      if (any(gap >= lo & gap <= hi)) kA <- kA + 1
    }
    settings$hits_A[k]   <- kA
    settings$n_events[k] <- length(ev_week)
    
    # precursor hits after circular shift by s weeks
    for (s in 1:n_shift) {
      ev_s <- (ev_week + s - 1) %% n_weeks + 1
      null_B[s, k] <- count_hits(ev_s, lo, hi)
    }
  }
  
  # -------------------------------------------------------------------
  # 3b: p-value of every setting (circular-shift test)
  # -------------------------------------------------------------------
  settings$p <- NA
  for (k in 1:n_set) {
    settings$p[k] <- (1 + sum(null_B[, k] >= settings$hits_B[k])) / (n_shift + 1)
  }
  settings$precursor_rate <- round(settings$hits_B / N_B, 2)
  settings$trigger_rate   <- round(settings$hits_A / settings$n_events, 3)
  
  # -------------------------------------------------------------------
  # 3c: p-value corrected for picking the best of the 24 settings
  # -------------------------------------------------------------------
  p_null <- matrix(NA, n_shift, n_set)
  for (k in 1:n_set) {
    for (s in 1:n_shift) {
      p_null[s, k] <- (sum(null_B[, k] >= null_B[s, k]) +
                         (settings$hits_B[k] >= null_B[s, k])) / (n_shift + 1)
    }
  }
  best_null <- apply(p_null, 1, min)          # best p of each shifted series
  b         <- which.min(settings$p)          # best real setting
  p_grid    <- (1 + sum(best_null <= settings$p[b])) / (n_shift + 1)
  
  cat("best: thr", settings$thr[b], "| dT", settings$dT[b], "| tau", settings$tau[b],
      "| hits", settings$hits_B[b], "of", N_B,
      "| p =", round(settings$p[b], 4), "| corrected =", round(p_grid, 4), "\n")
  
  # -------------------------------------------------------------------
  # 3d: save the grid table of this feature
  # -------------------------------------------------------------------
  settings$p <- round(settings$p, 4)
  write_csv(settings, file.path(out_dir, paste0("ECA_grid_", station, "_", f, ".csv")))
  
  # -------------------------------------------------------------------
  # 3e: figure of this feature (all thresholds, events of the best setting)
  # -------------------------------------------------------------------
  ev_best <- ev_list[[b]]
  ok      <- !is.na(x)
  
  pdf(file.path(out_dir, paste0("ECA_", station, "_", f, ".pdf")), width = 10, height = 4.5)
  par(mar = c(6, 4, 3, 4))                    # room for date labels and threshold labels
  plot(week[ok], x[ok], type = "l", xlab = "", ylab = f, xaxt = "n",
       xlim = range(week),
       main = paste0(station, " - ", f, " [", dir_f, "]  (threshold ", settings$thr[b],
                     ", dT = ", settings$dT[b], ", tau = ", settings$tau[b],
                     ", hits ", settings$hits_B[b], "/", N_B,
                     ", corrected p = ", round(p_grid, 3), ")"))
  draw_week_axis()
  mtext("week end date (UTC)", side = 1, line = 4.5)
  draw_thresholds(x, settings$thr[b], dir = dir_f)
  points(ev_best, x[match(ev_best, week)], col = "red", pch = 16)
  abline(v = erupt_week, col = "blue", lwd = 2)
  dev.off()
  
  # keep what the combined figure needs
  panel_data[[f]] <- list(x = x, ev = ev_best, thr = settings$thr[b], dir = dir_f,
                          dT = settings$dT[b], tau = settings$tau[b],
                          hits = settings$hits_B[b], p_grid = p_grid)
  
  # -------------------------------------------------------------------
  # 3f: one summary row for this feature
  # -------------------------------------------------------------------
  summary_all <- rbind(summary_all, data.frame(
    feature        = f,
    direction      = dir_f,
    thr            = settings$thr[b],
    dT             = settings$dT[b],
    tau            = settings$tau[b],
    n_events       = settings$n_events[b],
    hits_B         = settings$hits_B[b],
    hits_A         = settings$hits_A[b],
    precursor_rate = settings$precursor_rate[b],
    trigger_rate   = settings$trigger_rate[b],
    best_p         = round(settings$p[b], 4),
    p_grid         = round(p_grid, 4)))
}

# ---------------------------------------------------------------------
# Step 4: compare the features - Benjamini-Hochberg over all features
# ---------------------------------------------------------------------
summary_all$q_BH <- round(p.adjust(summary_all$p_grid, method = "BH"), 4)
summary_all <- summary_all[order(summary_all$p_grid), ]
print(summary_all)

write_csv(summary_all, file.path(out_dir, paste0("ECA_summary_", station, "_all_features.csv")))

# ---------------------------------------------------------------------
# Step 5: summary figure - corrected p-value of every feature
# ---------------------------------------------------------------------
draw_summary <- function() {
  par(mar = c(4, 9, 3, 1))
  barplot(rev(-log10(summary_all$p_grid)), names.arg = rev(summary_all$feature),
          horiz = TRUE, las = 1, col = ifelse(rev(summary_all$q_BH) < 0.05, "red", "grey60"),
          xlab = "-log10(corrected p)",
          main = paste(station, "- ECA, best setting per feature"))
  abline(v = -log10(0.05), lty = 2)
}
draw_summary()
pdf(file.path(out_dir, paste0("ECA_summary_", station, "_all_features.pdf")), width = 8, height = 6)
draw_summary()
dev.off()

# ---------------------------------------------------------------------
# Step 6: combined figure - all features stacked, shared week axis
#   ncol = 1 : one column, all panels on one (tall) page
#   ncol = 2 : two columns, fits an A4 page
#   order_by_p = TRUE puts the most significant feature on top
# ---------------------------------------------------------------------
draw_panels <- function(ncol = 1, order_by_p = FALSE) {
  feats <- if (order_by_p) summary_all$feature else all_features
  n     <- length(feats)
  nrow  <- ceiling(n / ncol)
  
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  par(mfrow = c(nrow, ncol))
  par(mar = c(0.4, 4, 0.4, 2.6), oma = c(6, 0, 3.5, 0),
      mgp = c(2.5, 0.45, 0), tcl = -0.25, cex = 0.9)
  
  for (i in seq_len(n)) {
    f  <- feats[i]
    d  <- panel_data[[f]]
    ok <- !is.na(d$x)
    q  <- summary_all$q_BH[summary_all$feature == f]
    
    plot(week[ok], d$x[ok], type = "l", xlim = range(week),
         xaxt = "n", xlab = "", ylab = "", las = 1, cex.axis = 0.7)
    draw_week_axis(labels = (i > n - ncol), cex_lab = 0.6)   # dates only on bottom panels
    draw_thresholds(d$x, d$thr, cex_lab = 0.55, dir = d$dir)
    points(d$ev, d$x[match(d$ev, week)], col = "red", pch = 16, cex = 0.6)
    abline(v = erupt_week, col = "blue", lwd = 1.5)
    
    lab <- paste0(f, " [", d$dir, "]   thr ", d$thr, ", dT ", d$dT, ", tau ", d$tau,
                  ", hits ", d$hits, "/", N_B,
                  ", p_grid ", round(d$p_grid, 3), ", q ", round(q, 3))
    legend("topleft", legend = lab, bty = "o", box.col = NA,
           bg = adjustcolor("white", 0.75), cex = 0.7, inset = 0.005,
           text.font = ifelse(q < 0.05, 2, 1))           # bold if q < 0.05
  }
  
  mtext("week end date (UTC)", side = 1, outer = TRUE, line = 4.5)
  mtext(paste0(station, " - ECA, all features (", dir_note, ")"),
        side = 3, outer = TRUE, line = 2, font = 2)
  mtext("dashed: threshold percentiles (bold = best setting)   red: events of best setting   blue: eruption onsets",
        side = 3, outer = TRUE, line = 0.6, cex = 0.7)
}

# one column, all 13 panels on one page
pdf(file.path(out_dir, paste0("ECA_", station, "_all_features_1col.pdf")),
    width = 10, height = 1.5 * length(all_features) + 1.2)
draw_panels(ncol = 1)
dev.off()

# two columns on an A4 page
pdf(file.path(out_dir, paste0("ECA_", station, "_all_features_2col.pdf")),
    width = 8.27, height = 11.69)
draw_panels(ncol = 2)
dev.off()

# =====================================================================
# Step 7: composite features
#   7a  percentile transform of every feature (rolling or global)
#   7b  group the features by correlation, orient them inside a group
#   7c  ECA core (same grid + circular-shift test + grid correction)
#   7d  median and vote-count composite of every group with >= 2 members
#   7e  Benjamini-Hochberg over the composites, comparison table
#   7f  stacked figure of the composites
# Groups are defined from the correlation structure only, never from
# the ECA results, so the corrections above stay valid.
# =====================================================================

# --- settings for Step 7 ---------------------------------------------
roll_win      <- 26     # trailing window (weeks) for the rolling percentile
use_rolling   <- TRUE   # FALSE = global percentile ranks (comparable with Step 3)
r_min         <- 0.7    # features joined at |Spearman r| >= r_min form a group
orient_manual <- c()    # optional manual orientation, e.g. c(H_Shannon = -1)

comp_dir <- file.path(out_dir, "composites")
dir.create(comp_dir, showWarnings = FALSE)

# ---------------------------------------------------------------------
# 7a: percentile transform
#     rolling: share of the last roll_win weeks that are <= this week
#     (causal - uses only past and current weeks)
# ---------------------------------------------------------------------
roll_pct <- function(x, w, min_n = ceiling(w / 2)) {
  out <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    if (is.na(x[i])) next
    win <- x[max(1, i - w + 1):i]
    win <- win[!is.na(win)]
    if (length(win) >= min_n) out[i] <- mean(win <= x[i])
  }
  out
}
glob_pct <- function(x) {
  rank(x, na.last = "keep", ties.method = "average") / sum(!is.na(x))
}

P <- sapply(all_features, function(f) {
  if (use_rolling) roll_pct(feat[[f]], roll_win) else glob_pct(feat[[f]])
})
# flip the "lower" features so that "unusual" means a high percentile for every feature
for (f in all_features) if (feat_dir[[f]] == "lower") P[, f] <- 1 - P[, f]

# ---------------------------------------------------------------------
# 7b: correlation groups and orientation
# ---------------------------------------------------------------------
X       <- as.data.frame(feat[, all_features])
cor_mat <- cor(X, method = "spearman", use = "pairwise.complete.obs")
hc      <- hclust(as.dist(1 - abs(cor_mat)), method = "average")
grp     <- cutree(hc, h = 1 - r_min)
groups  <- split(names(grp), grp)          # members keep the all_features order

write_csv(data.frame(feature = rownames(cor_mat), round(cor_mat, 3), check.names = FALSE),
          file.path(comp_dir, paste0("Feature_spearman_", station, ".csv")))

pdf(file.path(comp_dir, paste0("Feature_dendrogram_", station, ".pdf")), width = 8, height = 5)
par(mar = c(2, 4, 3, 1))
plot(hc, xlab = "", sub = "", ylab = "1 - |Spearman r|",
     main = paste(station, "- feature groups (average linkage)"))
abline(h = 1 - r_min, lty = 2, col = "red")
rect.hclust(hc, h = 1 - r_min, border = "grey40")
dev.off()

# orientation: each member follows the sign of its correlation with the
# group's reference feature (its first member), so they do not cancel
orient <- setNames(rep(1, length(all_features)), all_features)
ref    <- setNames(rep(NA_character_, length(all_features)), all_features)
for (g in groups) {
  for (f in g) {
    # sign of the correlation, corrected for the features' own directions
    d_f   <- ifelse(feat_dir[[f]]    == "lower", -1, 1)
    d_ref <- ifelse(feat_dir[[g[1]]] == "lower", -1, 1)
    orient[f] <- ifelse(cor_mat[f, g[1]] * d_f * d_ref >= 0, 1, -1)
    ref[f]    <- g[1]
  }
}
orient[names(orient_manual)] <- orient_manual
if (any(orient == -1)) {
  cat("Note: flipped against their own direction inside their group:",
      names(orient)[orient == -1], "\n")
}

O <- P
for (f in all_features) if (orient[f] == -1) O[, f] <- 1 - P[, f]

group_table <- data.frame(
  feature     = all_features,
  direction   = feat_dir[all_features],
  group       = grp[all_features],
  reference   = ref[all_features],
  r_with_ref  = round(sapply(all_features, function(f) cor_mat[f, ref[f]]), 3),
  orientation = orient[all_features])
print(group_table, row.names = FALSE)
write_csv(group_table, file.path(comp_dir, paste0("Feature_groups_", station, ".csv")))

# ---------------------------------------------------------------------
# 7c: ECA core - identical to Steps 3a-3c, but the events come from
#     ev_fun(thr), which returns TRUE/FALSE for every row
# ---------------------------------------------------------------------
eca_core <- function(ev_fun) {
  settings <- grid
  settings$n_events <- NA
  settings$hits_B   <- NA
  settings$hits_A   <- NA
  null_B  <- matrix(NA, n_shift, n_set)
  ev_list <- list()
  
  for (k in 1:n_set) {
    ev <- ev_fun(settings$thr[k])
    ev[is.na(ev)] <- FALSE
    ev_week <- week[ev & !c(FALSE, head(ev, -1))]    # first week of each run
    ev_list[[k]] <- ev_week
    
    lo <- 1 + settings$tau[k]
    hi <- settings$dT[k] + settings$tau[k]
    
    settings$hits_B[k] <- count_hits(ev_week, lo, hi)
    kA <- 0
    for (i in seq_along(ev_week)) {
      gap <- erupt_week - ev_week[i]
      if (any(gap >= lo & gap <= hi)) kA <- kA + 1
    }
    settings$hits_A[k]   <- kA
    settings$n_events[k] <- length(ev_week)
    
    for (s in 1:n_shift) {
      ev_s <- (ev_week + s - 1) %% n_weeks + 1
      null_B[s, k] <- count_hits(ev_s, lo, hi)
    }
  }
  
  settings$p <- NA
  for (k in 1:n_set) {
    settings$p[k] <- (1 + sum(null_B[, k] >= settings$hits_B[k])) / (n_shift + 1)
  }
  settings$precursor_rate <- round(settings$hits_B / N_B, 2)
  settings$trigger_rate   <- round(settings$hits_A / settings$n_events, 3)
  
  p_null <- matrix(NA, n_shift, n_set)
  for (k in 1:n_set) {
    for (s in 1:n_shift) {
      p_null[s, k] <- (sum(null_B[, k] >= null_B[s, k]) +
                         (settings$hits_B[k] >= null_B[s, k])) / (n_shift + 1)
    }
  }
  best_null <- apply(p_null, 1, min)
  b         <- which.min(settings$p)
  p_grid    <- (1 + sum(best_null <= settings$p[b])) / (n_shift + 1)
  
  list(settings = settings, b = b, p_grid = p_grid, ev_list = ev_list)
}

# ---------------------------------------------------------------------
# 7d: composites of every group with >= 2 members
#     median : median of the oriented percentiles, thresholded at its own
#              grid percentiles (as in Step 3)
#     vote   : number of members above the grid threshold;
#              event = strict majority of the members
# ---------------------------------------------------------------------
comp_info    <- list()
summary_comp <- NULL

for (k_g in seq_along(groups)) {
  g <- groups[[k_g]]
  if (length(g) < 2) {
    cat("Group", k_g, "has one member (", g, ") - already covered by Step 3\n")
    next
  }
  m         <- length(g)
  Og        <- O[, g, drop = FALSE]
  med       <- apply(Og, 1, median)        # NA unless every member is available
  min_votes <- floor(m / 2) + 1
  
  for (type in c("median", "vote")) {
    name <- paste0("G", k_g, "_", type)
    cat("\n=====", name, ":", paste(g, collapse = ", "), "=====\n")
    
    if (type == "median") {
      res <- eca_core(function(thr) !is.na(med) & med > quantile(med, thr, na.rm = TRUE))
    } else {
      res <- eca_core(function(thr) {
        v <- rowSums(Og > thr)
        !is.na(v) & v >= min_votes
      })
    }
    s <- res$settings
    b <- res$b
    cat("best: thr", s$thr[b], "| dT", s$dT[b], "| tau", s$tau[b],
        "| hits", s$hits_B[b], "of", N_B,
        "| p =", round(s$p[b], 4), "| corrected =", round(res$p_grid, 4), "\n")
    
    s_out   <- s
    s_out$p <- round(s_out$p, 4)
    write_csv(s_out, file.path(comp_dir, paste0("ECA_grid_", station, "_", name, ".csv")))
    
    comp_info[[name]] <- list(type = type, members = g, Og = Og, med = med,
                              min_votes = min_votes, res = res)
    
    summary_comp <- rbind(summary_comp, data.frame(
      composite      = name,
      type           = type,
      members        = paste(g, collapse = "+"),
      min_votes      = ifelse(type == "vote", min_votes, NA),
      thr            = s$thr[b],
      dT             = s$dT[b],
      tau            = s$tau[b],
      n_events       = s$n_events[b],
      hits_B         = s$hits_B[b],
      hits_A         = s$hits_A[b],
      precursor_rate = s$precursor_rate[b],
      trigger_rate   = s$trigger_rate[b],
      best_p         = round(s$p[b], 4),
      p_grid         = round(res$p_grid, 4)))
  }
}

# ---------------------------------------------------------------------
# 7e: Benjamini-Hochberg over the composites (their own family)
#     + one table next to the single features of Step 3
# ---------------------------------------------------------------------
if (!is.null(summary_comp)) {
  summary_comp$q_BH <- round(p.adjust(summary_comp$p_grid, method = "BH"), 4)
  summary_comp <- summary_comp[order(summary_comp$p_grid), ]
  print(summary_comp, row.names = FALSE)
  write_csv(summary_comp, file.path(comp_dir, paste0("ECA_summary_", station, "_composites.csv")))
  
  comparison <- rbind(
    data.frame(name = summary_all$feature, set = "single (Step 3)",
               thr = summary_all$thr, dT = summary_all$dT, tau = summary_all$tau,
               hits_B = summary_all$hits_B, n_events = summary_all$n_events,
               p_grid = summary_all$p_grid, q_BH_within_set = summary_all$q_BH),
    data.frame(name = paste0(summary_comp$composite, " (", summary_comp$members, ")"),
               set = paste0("composite (", ifelse(use_rolling,
                                                  paste0("rolling ", roll_win, " wk"),
                                                  "global"), ")"),
               thr = summary_comp$thr, dT = summary_comp$dT, tau = summary_comp$tau,
               hits_B = summary_comp$hits_B, n_events = summary_comp$n_events,
               p_grid = summary_comp$p_grid, q_BH_within_set = summary_comp$q_BH))
  write_csv(comparison, file.path(comp_dir, paste0("ECA_comparison_", station, ".csv")))
}

# ---------------------------------------------------------------------
# 7f: stacked figure of the composites
# ---------------------------------------------------------------------
draw_comp_panels <- function() {
  nm <- names(comp_info)
  n  <- length(nm)
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  par(mfrow = c(n, 1))
  par(mar = c(0.4, 4, 0.4, 2.6), oma = c(6, 0, 3.5, 0),
      mgp = c(2.5, 0.45, 0), tcl = -0.25, cex = 0.9)
  
  for (i in seq_len(n)) {
    ci <- comp_info[[nm[i]]]
    s  <- ci$res$settings
    b  <- ci$res$b
    ev <- ci$res$ev_list[[b]]
    q  <- summary_comp$q_BH[summary_comp$composite == nm[i]]
    
    if (ci$type == "median") {
      y  <- ci$med
      ok <- !is.na(y)
      plot(week[ok], y[ok], type = "l", xlim = range(week), ylim = c(0, 1),
           xaxt = "n", xlab = "", ylab = "", las = 1, cex.axis = 0.7)
      draw_thresholds(y, s$thr[b], cex_lab = 0.55, dir = "upper")
      lab_rule <- "median of oriented percentiles"
    } else {
      y  <- rowSums(ci$Og > s$thr[b])        # votes at the best threshold
      ok <- !is.na(y)
      plot(week[ok], y[ok], type = "s", xlim = range(week),
           ylim = c(0, length(ci$members)),
           xaxt = "n", xlab = "", ylab = "", las = 1, cex.axis = 0.7)
      col_b <- thr_cols[match(s$thr[b], thr_grid)]
      abline(h = ci$min_votes - 0.5, lty = 2, lwd = 2, col = col_b)
      text(par("usr")[2], ci$min_votes - 0.5, paste0(">=", ci$min_votes),
           pos = 4, xpd = NA, col = col_b, cex = 0.55, font = 2)
      lab_rule <- paste0("votes: members above P", round(s$thr[b] * 100))
    }
    points(ev, y[match(ev, week)], col = "red", pch = 16, cex = 0.6)
    abline(v = erupt_week, col = "blue", lwd = 1.5)
    draw_week_axis(labels = (i == n), cex_lab = 0.6)
    
    lab <- paste0(nm[i], " [", paste(ci$members, collapse = ", "), "]  ", lab_rule,
                  "   thr ", s$thr[b], ", dT ", s$dT[b], ", tau ", s$tau[b],
                  ", hits ", s$hits_B[b], "/", N_B,
                  ", p_grid ", round(ci$res$p_grid, 3), ", q ", round(q, 3))
    legend("topleft", legend = lab, bty = "o", box.col = NA,
           bg = adjustcolor("white", 0.75), cex = 0.6, inset = 0.005,
           text.font = ifelse(q < 0.05, 2, 1))
  }
  
  mtext("week end date (UTC)", side = 1, outer = TRUE, line = 4.5)
  mtext(paste0(station, " - ECA, composite features (", dir_note, "; ",
               ifelse(use_rolling, paste0("rolling ", roll_win, "-week percentiles"),
                      "global percentiles"), ")"),
        side = 3, outer = TRUE, line = 2, font = 2)
  mtext("dashed: thresholds (bold = best setting)   red: events of best setting   blue: eruption onsets",
        side = 3, outer = TRUE, line = 0.6, cex = 0.7)
}

if (length(comp_info) > 0) {
  pdf(file.path(comp_dir, paste0("ECA_", station, "_composites.pdf")),
      width = 10, height = 1.6 * length(comp_info) + 1.2)
  draw_comp_panels()
  dev.off()
}

cat("\nAll tables and figures saved in:", out_dir, "\n")