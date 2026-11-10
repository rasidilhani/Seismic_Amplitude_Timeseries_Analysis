# WIZ_OrdinalPatterns_D5_Hourly_Every3Days_2019.R
#
# Same pipeline as WIZ_OrdinalPatterns_D5_Hourly_Last6Months_2019.R (same two
# packages, same compute_hour() math), but instead of processing EVERY day
# from Jul 1 - Dec 31, it processes only every 3rd day starting 3 July:
#   Jul 3, 6, 9, 12, 15, ... , Dec 9, Dec 12, ... through Dec 2019.
#
# Start date is 3 July (not 1 July) specifically so that 9 December falls
# exactly on the 3-day grid: 9 Dec is day-of-year 343, and 343 - 184 = 159,
# which divides evenly by 3 (53 steps). Starting from 1 or 2 July does NOT
# land on 9 Dec -- the nearest grid points from those starts are 8 Dec / 11 Dec.
# This script PRINTS the full resulting date list at the start so you can
# double-check it -- change `start_doy` / `step_days` below if needed.
#
# Results are saved incrementally: as soon as one day's hourly results are
# computed, they are appended to the single output CSV immediately -- so
# nothing is lost if the run is interrupted partway through.

rm(list = ls())

library(tidyverse)
library(lubridate)
library(StatOrdPattHxC)
library(ordinalpatterns)

# ── Parameters (identical to previous scripts) ─────────────────────────────
D       <- 5
tau     <- 1
z_alpha <- qnorm(0.975)
BETA    <- 1.5
var_col <- "rsam"

step_days <- 3   # <-- change this if you want a different sampling interval
start_doy <- 184 # <-- 3 July 2019; chosen so 9 Dec (doy 343) lands exactly on the grid

# ── Folders ─────────────────────────────────────────────────────────────────
input_dir  <- "E:/Seismic_Amplitude_Timeseries_Analysis/2019 secods RSAM"
output_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Build the day-of-year sequence: every `step_days` days, from start_doy ──
# 2019 is not a leap year: Jan(31)+Feb(28)+Mar(31)+Apr(30)+May(31)+Jun(30)=181
# -> day 184 = 3 July 2019 (start), day 365 = 31 Dec 2019 (end of range)
# start_doy = 184 is deliberate: it's the only nearby start date that puts
# 9 December (doy 343) exactly on the 3-day grid.
selected_doy <- seq(start_doy, 365, by = step_days)

# Print the resulting calendar dates for a sanity check against what you intended
selected_dates <- as.Date(selected_doy - 1, origin = "2019-01-01")
cat("Selected days (every", step_days, "days, starting 3 July 2019):\n")
print(data.frame(doy = selected_doy, date = selected_dates))
cat("\nTotal days selected:", length(selected_doy), "\n\n")

# ── Match selected days to files on the external drive ─────────────────────
all_files <- list.files(input_dir, pattern = "^2019\\.\\d{3}\\.WIZ\\.RSAM\\.csv$",
                        full.names = TRUE)

file_info <- tibble::tibble(path = all_files) %>%
  dplyr::mutate(
    fname = basename(path),
    doy   = as.integer(stringr::str_extract(fname, "(?<=^2019\\.)\\d{3}"))
  ) %>%
  dplyr::filter(doy %in% selected_doy) %>%
  dplyr::arrange(doy)

cat("Files found on drive matching selected days:", nrow(file_info), "of", length(selected_doy), "expected\n")
missing_doy <- setdiff(selected_doy, file_info$doy)
if (length(missing_doy) > 0) {
  cat("WARNING: no file found for these day(s)-of-year:", paste(missing_doy, collapse = ", "), "\n")
}
cat("\n")
stopifnot(nrow(file_info) > 0)

# ── Per-hour feature computation (unchanged from previous scripts) ─────────
compute_hour <- function(hour, d) {
  
  sub_df <- d %>% dplyr::filter(hour_bin == hour)
  series <- sub_df[[var_col]]
  series <- series[is.finite(series)]
  n_i    <- length(series)
  
  n_pats <- factorial(D)
  
  if (n_i < D) {
    warning(sprintf("Hour %s has only %d points (< D = %d) -- skipping.",
                    format(hour, tz = "Pacific/Auckland"), n_i, D))
    return(NULL)
  }
  
  n_eff <- n_i - D + 1
  
  xpd   <- ordinalpatterns::op_pd(series, D = D, tau = tau)
  probs <- xpd$probabilities
  
  if (length(probs) != n_pats) {
    warning(sprintf(
      "Hour %s: op_pd() returned %d probabilities, expected %d -- investigate.",
      format(hour, tz = "Pacific/Auckland"), length(probs), n_pats
    ))
  }
  
  Hs <- ordinalpatterns::permutation_entropy(probs, normalized = TRUE)
  Cs <- ordinalpatterns::statistical_complexity(probs, entropy = Hs, normalized = TRUE)
  Fi <- ordinalpatterns::fisher_information(probs)
  
  Hr <- StatOrdPattHxC::HRenyi(probs, beta = BETA)
  Ht <- StatOrdPattHxC::HTsallis(probs, beta = BETA)
  
  JS <- ordinalpatterns::jsd(probs) / log(2)
  
  Cr <- JS * Hr
  Ct <- JS * Ht
  Cf <- JS * Fi
  
  Var_Hs <- suppressWarnings(StatOrdPattHxC::sigma2q(series, emb = D, ent = "S"))
  Var_Hr <- suppressWarnings(StatOrdPattHxC::sigma2q(series, emb = D, ent = "R", beta = BETA))
  Var_Ht <- suppressWarnings(StatOrdPattHxC::sigma2q(series, emb = D, ent = "T", beta = BETA))
  Var_Hf <- suppressWarnings(StatOrdPattHxC::sigma2q(series, emb = D, ent = "F"))
  
  Var_HI <- suppressWarnings(StatOrdPattHxC::asymptoticVarHShannonMultinomial(probs, n_eff))
  Var_CI <- suppressWarnings(StatOrdPattHxC::varC(probs, n_eff))
  
  a_ratio <- ifelse(Var_HI > 0, Var_Hs / Var_HI, NA)
  Var_Cs  <- a_ratio * Var_CI
  
  semi <- function(v) ifelse(!is.finite(v) | v <= 0, NA, sqrt(v) / sqrt(n_eff) * z_alpha)
  
  tibble::tibble(
    Hour_Start = hour, N_points = n_i, N_eff = n_eff,
    
    H_Shannon = Hs, C_Shannon = Cs, 
    Fisher_Info = Fi, C_Fisher = Cf,
    H_Renyi   = Hr, C_Renyi   = Cr,
    H_Tsallis = Ht, C_Tsallis = Ct,
    
    Var_H_Shannon = Var_Hs, Var_C_Shannon = Var_Cs,
    Var_H_Renyi   = Var_Hr,
    Var_H_Tsallis = Var_Ht,
    Var_H_Fisher  = Var_Hf,
    
    Semi_H_Shannon = semi(Var_Hs), Semi_C_Shannon = semi(Var_Cs),
    Semi_H_Renyi   = semi(Var_Hr),
    Semi_H_Tsallis = semi(Var_Ht),
    Semi_H_Fisher  = semi(Var_Hf)
  )
}

# ── Per-day processing: read, average 100Hz -> 1Hz, bin hourly, compute ────
process_day <- function(path, fname) {
  
  d <- readr::read_csv(path, show_col_types = FALSE)
  
  if (!all(c("unix_timestamp", var_col) %in% names(d))) {
    warning(sprintf("%s: missing required columns -- skipping file.", fname))
    return(NULL)
  }
  
  d <- d %>%
    dplyr::group_by(unix_timestamp) %>%
    dplyr::summarise(rsam = mean(.data[[var_col]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(unix_timestamp)
  
  utc_time      <- as.POSIXct(d$unix_timestamp, origin = "1970-01-01", tz = "UTC")
  d$datetime_nz <- lubridate::with_tz(utc_time, tzone = "Pacific/Auckland")
  d$hour_bin    <- lubridate::floor_date(d$datetime_nz, unit = "hour")
  
  day_hours <- sort(unique(d$hour_bin))
  
  purrr::map_dfr(day_hours, ~ compute_hour(.x, d))
}

# ── Output file: single CSV, unique to the second, written incrementally ───
timestamp_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_file <- file.path(output_dir,
                      paste0("WIZ_OP_Hourly_Every_", step_days,
                             "Days_2019_", timestamp_tag, ".csv"))

cat("Output file (created after first day completes):\n", out_file, "\n\n")

# ── Run across selected days, appending to the CSV as each day finishes ────
n_files <- nrow(file_info)
t_start <- Sys.time()
file_created <- FALSE

for (i in seq_len(n_files)) {
  fname <- file_info$fname[i]
  path  <- file_info$path[i]
  
  cat(sprintf("[%d/%d] Processing %s ... ", i, n_files, fname))
  
  t_day <- system.time({
    day_result <- process_day(path, fname)
  })
  
  n_hours <- ifelse(is.null(day_result), 0L, nrow(day_result))
  cat(sprintf("%.1f sec, %d hours\n", t_day["elapsed"], n_hours))
  
  if (!is.null(day_result) && nrow(day_result) > 0) {
    # First write creates the file with headers; every write after that
    # appends without repeating the header -- so results are saved to disk
    # the moment each day finishes, not held in memory until the very end.
    readr::write_csv(day_result, out_file, append = file_created, col_names = !file_created)
    file_created <- TRUE
  }
}

t_end <- Sys.time()
cat(sprintf("\nAll selected days processed in %.1f minutes.\n",
            as.numeric(difftime(t_end, t_start, units = "mins"))))
cat("Final combined results saved to:\n", out_file, "\n")


############################################
# WIZ_OrdinalPatterns_D5_Hourly_6July2019.R
#
# Single-day version of WIZ_OrdinalPatterns_D5_Hourly_Every3Days_2019.R.
# Runs the SAME pipeline (same packages, same compute_hour() math) but
# restricted to ONLY 6 July 2019 (day-of-year 187) instead of the full
# every-3-days range from July to December.
#
# Results are still saved incrementally (though with only one day, this
# just means the CSV is written as soon as that day's hourly results
# are computed).

rm(list = ls())

library(tidyverse)
library(lubridate)
library(StatOrdPattHxC)
library(ordinalpatterns)

# ── Parameters (identical to previous scripts) ─────────────────────────────
D       <- 5
tau     <- 1
z_alpha <- qnorm(0.975)
BETA    <- 1.5
var_col <- "rsam"

# ── Folders ─────────────────────────────────────────────────────────────────
input_dir  <- "E:/Seismic_Amplitude_Timeseries_Analysis/2019 secods RSAM"
output_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Single day selection: 6 July 2019 = day-of-year 187 ─────────────────────
# 2019 is not a leap year: Jan(31)+Feb(28)+Mar(31)+Apr(30)+May(31)+Jun(30)=181
# -> day 187 = 6 July 2019
selected_doy <- 187

# Print the resulting calendar date for a sanity check against what you intended
selected_dates <- as.Date(selected_doy - 1, origin = "2019-01-01")
cat("Selected day:\n")
print(data.frame(doy = selected_doy, date = selected_dates))
cat("\n")

# ── Match selected day to the file on the external drive ───────────────────
all_files <- list.files(input_dir, pattern = "^2019\\.\\d{3}\\.WIZ\\.RSAM\\.csv$",
                        full.names = TRUE)

file_info <- tibble::tibble(path = all_files) %>%
  dplyr::mutate(
    fname = basename(path),
    doy   = as.integer(stringr::str_extract(fname, "(?<=^2019\\.)\\d{3}"))
  ) %>%
  dplyr::filter(doy %in% selected_doy) %>%
  dplyr::arrange(doy)

cat("Files found on drive matching selected day:", nrow(file_info), "of", length(selected_doy), "expected\n")
missing_doy <- setdiff(selected_doy, file_info$doy)
if (length(missing_doy) > 0) {
  cat("WARNING: no file found for this day-of-year:", paste(missing_doy, collapse = ", "), "\n")
}
cat("\n")
stopifnot(nrow(file_info) > 0)

# ── Per-hour feature computation (unchanged from previous scripts) ─────────
compute_hour <- function(hour, d) {
  
  sub_df <- d %>% dplyr::filter(hour_bin == hour)
  series <- sub_df[[var_col]]
  series <- series[is.finite(series)]
  n_i    <- length(series)
  
  n_pats <- factorial(D)
  
  if (n_i < D) {
    warning(sprintf("Hour %s has only %d points (< D = %d) -- skipping.",
                    format(hour, tz = "Pacific/Auckland"), n_i, D))
    return(NULL)
  }
  
  n_eff <- n_i - D + 1
  
  xpd   <- ordinalpatterns::op_pd(series, D = D, tau = tau)
  probs <- xpd$probabilities
  
  if (length(probs) != n_pats) {
    warning(sprintf(
      "Hour %s: op_pd() returned %d probabilities, expected %d -- investigate.",
      format(hour, tz = "Pacific/Auckland"), length(probs), n_pats
    ))
  }
  
  Hs <- ordinalpatterns::permutation_entropy(probs, normalized = TRUE)
  Cs <- ordinalpatterns::statistical_complexity(probs, entropy = Hs, normalized = TRUE)
  Fi <- ordinalpatterns::fisher_information(probs)
  
  Hr <- StatOrdPattHxC::HRenyi(probs, beta = BETA)
  Ht <- StatOrdPattHxC::HTsallis(probs, beta = BETA)
  
  JS <- ordinalpatterns::jsd(probs) / log(2)
  
  Cr <- JS * Hr
  Ct <- JS * Ht
  
  Var_Hs <- suppressWarnings(StatOrdPattHxC::sigma2q(series, emb = D, ent = "S"))
  Var_Hr <- suppressWarnings(StatOrdPattHxC::sigma2q(series, emb = D, ent = "R", beta = BETA))
  Var_Ht <- suppressWarnings(StatOrdPattHxC::sigma2q(series, emb = D, ent = "T", beta = BETA))
  Var_Hf <- suppressWarnings(StatOrdPattHxC::sigma2q(series, emb = D, ent = "F"))
  
  Var_HI <- suppressWarnings(StatOrdPattHxC::asymptoticVarHShannonMultinomial(probs, n_eff))
  Var_CI <- suppressWarnings(StatOrdPattHxC::varC(probs, n_eff))
  
  a_ratio <- ifelse(Var_HI > 0, Var_Hs / Var_HI, NA)
  Var_Cs  <- a_ratio * Var_CI
  
  semi <- function(v) ifelse(!is.finite(v) | v <= 0, NA, sqrt(v) / sqrt(n_eff) * z_alpha)
  
  tibble::tibble(
    Hour_Start = hour, N_points = n_i, N_eff = n_eff,
    
    H_Shannon = Hs, C_Shannon = Cs, Fisher_Info = Fi,
    H_Renyi   = Hr, C_Renyi   = Cr,
    H_Tsallis = Ht, C_Tsallis = Ct,
    
    Var_H_Shannon = Var_Hs, Var_C_Shannon = Var_Cs,
    Var_H_Renyi   = Var_Hr,
    Var_H_Tsallis = Var_Ht,
    Var_H_Fisher  = Var_Hf,
    
    Semi_H_Shannon = semi(Var_Hs), Semi_C_Shannon = semi(Var_Cs),
    Semi_H_Renyi   = semi(Var_Hr),
    Semi_H_Tsallis = semi(Var_Ht),
    Semi_H_Fisher  = semi(Var_Hf)
  )
}

# ── Per-day processing: read, average 100Hz -> 1Hz, bin hourly, compute ────
process_day <- function(path, fname) {
  
  d <- readr::read_csv(path, show_col_types = FALSE)
  
  if (!all(c("unix_timestamp", var_col) %in% names(d))) {
    warning(sprintf("%s: missing required columns -- skipping file.", fname))
    return(NULL)
  }
  
  d <- d %>%
    dplyr::group_by(unix_timestamp) %>%
    dplyr::summarise(rsam = mean(.data[[var_col]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(unix_timestamp)
  
  utc_time      <- as.POSIXct(d$unix_timestamp, origin = "1970-01-01", tz = "UTC")
  d$datetime_nz <- lubridate::with_tz(utc_time, tzone = "Pacific/Auckland")
  d$hour_bin    <- lubridate::floor_date(d$datetime_nz, unit = "hour")
  
  day_hours <- sort(unique(d$hour_bin))
  
  purrr::map_dfr(day_hours, ~ compute_hour(.x, d))
}

# ── Output file: single CSV, unique to the second, written incrementally ───
timestamp_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_file <- file.path(output_dir,
                      paste0("WIZ_OrdinalPatterns_D5_Hourly_6July2019_", timestamp_tag, ".csv"))

cat("Output file (created after day completes):\n", out_file, "\n\n")

# ── Run for the single selected day, appending to the CSV as it finishes ───
n_files <- nrow(file_info)
t_start <- Sys.time()
file_created <- FALSE

for (i in seq_len(n_files)) {
  fname <- file_info$fname[i]
  path  <- file_info$path[i]
  
  cat(sprintf("[%d/%d] Processing %s ... ", i, n_files, fname))
  
  t_day <- system.time({
    day_result <- process_day(path, fname)
  })
  
  n_hours <- ifelse(is.null(day_result), 0L, nrow(day_result))
  cat(sprintf("%.1f sec, %d hours\n", t_day["elapsed"], n_hours))
  
  if (!is.null(day_result) && nrow(day_result) > 0) {
    # First write creates the file with headers; every write after that
    # appends without repeating the header -- so results are saved to disk
    # the moment each day finishes, not held in memory until the very end.
    readr::write_csv(day_result, out_file, append = file_created, col_names = !file_created)
    file_created <- TRUE
  }
}

t_end <- Sys.time()
cat(sprintf("\nDay processed in %.1f minutes.\n",
            as.numeric(difftime(t_end, t_start, units = "mins"))))
cat("Final results saved to:\n", out_file, "\n")