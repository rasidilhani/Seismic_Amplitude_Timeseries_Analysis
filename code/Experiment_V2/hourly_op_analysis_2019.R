#############################################
library(ordinalpatterns)

P <- c(0.5, 0.3, 0.2, 0, 0)
cat("jsd(P, P)          =", jsd(P, P), "  (expect ~0)\n")

A <- c(1, 0)
B <- c(0, 1)
cat("jsd(A, B)          =", jsd(A, B), "  (1 = normalized, 0.693 = unnormalized)\n")

p_uniform_5 <- rep(1/5, 5)
cat("jsd(uniform-5)     =", jsd(p_uniform_5), "  (expect ~0 if Q=NULL is 1/length(P))\n")

Jensen_Shannon_manual <- function(p, q) {
  m  <- 0.5 * (p + q)
  js <- 0.5 * sum(ifelse(p == 0, 0, p * log((p + 1e-12) / (m + 1e-12)))) +
    0.5 * sum(ifelse(q == 0, 0, q * log((q + 1e-12) / (m + 1e-12))))
  js / log(2)
}

set.seed(1)
x    <- rnorm(2000)
D    <- 5
xpd  <- op_pd(x, D = D, tau = 1)
prob <- xpd$probabilities

cat("length(prob)       =", length(prob), " (expect", factorial(D), "if full D! support)\n")

Pe <- rep(1 / length(prob), length(prob))
cat("manual JS (norm.)  =", Jensen_Shannon_manual(prob, Pe), "\n")
cat("jsd(prob)          =", jsd(prob), "\n")
cat("jsd(prob)/log(2)   =", jsd(prob) / log(2), "\n")

###########################################################
# test_one_hour_2019.R
#
# Sanity-check run: picks ONE hour out of the 2019 per-second RSAM data and
# runs the exact same pipeline as WIZ_OrdinalPatterns_D5_Hourly_2019.R,
# printing every intermediate value AND timing so you can eyeball it before
# committing to the full-year run.
#
# CHANGE FROM PREVIOUS VERSION:
#   Raw file is 100 Hz (100 rows share each unix_timestamp). We now collapse
#   to one averaged value per second BEFORE any ordinal-pattern computation,
#   so each hour has 3600 points instead of 360000. This is both physically
#   correct (RSAM should be one amplitude value per second) and removes what
#   was almost certainly the runtime bottleneck: op_pd() + four sigma2q()
#   calls were each processing 360k raw points per hour.

rm(list = ls())

library(tidyverse)
library(lubridate)
library(StatOrdPattHxC)
library(ordinalpatterns)
library(here)

# ── Parameters (identical to the main script) ─────────────────────────────
D       <- 5
tau     <- 1
z_alpha <- qnorm(0.975)
BETA    <- 1.5

var_col <- "rsam"

# ── File to check ──────────────────────────────────────────────────────────
data_path <- here("data", "2019.113.WIZ.RSAM.csv")

cat("Reading:", data_path, "\n")
d <- readr::read_csv(data_path, show_col_types = FALSE)

stopifnot(all(c("unix_timestamp", var_col) %in% names(d)))

cat("Rows before averaging:", nrow(d), "\n")

# ── Collapse 100 Hz raw samples to 1 Hz averages ───────────────────────────
# Check block size is consistent (i.e. no dropped samples within a second)
block_sizes <- d %>% dplyr::count(unix_timestamp) %>% dplyr::pull(n)
cat("Samples-per-second summary:\n")
print(summary(block_sizes))
if (length(unique(block_sizes)) > 1) {
  cat("WARNING: inconsistent samples-per-second -- some seconds have dropped/extra samples.\n")
}

d <- d %>%
  dplyr::group_by(unix_timestamp) %>%
  dplyr::summarise(rsam = mean(.data[[var_col]], na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(unix_timestamp)

var_col <- "rsam"   # unchanged, just re-affirming after summarise
cat("Rows after averaging to 1 Hz:", nrow(d), "\n\n")

# ── Timestamp handling: unix -> UTC -> NZ (same pattern as every other script) ─
utc_time      <- as.POSIXct(d$unix_timestamp, origin = "1970-01-01", tz = "UTC")
d$datetime_nz <- lubridate::with_tz(utc_time, tzone = "Pacific/Auckland")
d$year        <- lubridate::year(d$datetime_nz)
d$hour_bin    <- lubridate::floor_date(d$datetime_nz, unit = "hour")

cat("Total rows in file:", nrow(d), "\n")
cat("Date range (NZ):", format(min(d$datetime_nz)), "to", format(max(d$datetime_nz)), "\n")

all_hours <- sort(unique(d$hour_bin))
cat("Number of distinct hour bins found:", length(all_hours), "\n\n")

# ── Per-hour feature computation (same math as the main script) ───────────
compute_hour <- function(hour) {
  
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
  cat("N points this hour (post-averaging):", n_i, " | N_eff:", n_eff, "\n")
  
  # ---- ordinalpatterns: Shannon entropy, Shannon complexity, Fisher info ----
  t_opd <- system.time({
    xpd   <- ordinalpatterns::op_pd(series, D = D, tau = tau)
    probs <- xpd$probabilities
  })
  cat("op_pd():", t_opd["elapsed"], "sec\n")
  
  if (length(probs) != n_pats) {
    warning(sprintf(
      "Hour %s: op_pd() returned %d probabilities, expected %d -- investigate.",
      format(hour, tz = "Pacific/Auckland"), length(probs), n_pats
    ))
  }
  
  Hs <- ordinalpatterns::permutation_entropy(probs, normalized = TRUE)
  Cs <- ordinalpatterns::statistical_complexity(probs, entropy = Hs, normalized = TRUE)
  Fi <- ordinalpatterns::fisher_information(probs)
  
  # ---- StatOrdPattHxC: Renyi, Tsallis, variances ----
  Hr <- StatOrdPattHxC::HRenyi(probs, beta = BETA)
  Ht <- StatOrdPattHxC::HTsallis(probs, beta = BETA)
  
  JS <- ordinalpatterns::jsd(probs) / log(2)   # natural-log based, normalised to log2
  
  Cr <- JS * Hr
  Ct <- JS * Ht
  Cf <- JS * Fi
  
  t_sigma <- system.time({
    Var_Hs <- suppressWarnings(StatOrdPattHxC::sigma2q(series, emb = D, ent = "S"))
    Var_Hr <- suppressWarnings(StatOrdPattHxC::sigma2q(series, emb = D, ent = "R", beta = BETA))
    Var_Ht <- suppressWarnings(StatOrdPattHxC::sigma2q(series, emb = D, ent = "T", beta = BETA))
    Var_Hf <- suppressWarnings(StatOrdPattHxC::sigma2q(series, emb = D, ent = "F"))
  })
  cat("4x sigma2q():", t_sigma["elapsed"], "sec  (", t_sigma["elapsed"] / 4, "sec/call )\n")
  
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

test_hour <- all_hours[1]
cat("Testing hour:", format(test_hour, tz = "Pacific/Auckland"), "\n")

sub_df <- d %>% dplyr::filter(hour_bin == test_hour)
cat("N points in this hour:", nrow(sub_df), "\n")
cat("First few RSAM values:", head(sub_df[[var_col]], 5), "\n\n")

t_total <- system.time({
  result <- compute_hour(test_hour)
})
cat("\nTOTAL TIME for one hour:", t_total["elapsed"], "sec\n\n")

cat("── Result for this hour ──────────────────────────────\n")
print(as.data.frame(t(result)))

cat("\n── Sanity checks ─────────────────────────────────────\n")
cat("H_Shannon in [0,1]?      ", result$H_Shannon >= 0 && result$H_Shannon <= 1, " (value:", result$H_Shannon, ")\n")
cat("C_Shannon in [0,1]?      ", result$C_Shannon >= 0 && result$C_Shannon <= 1, " (value:", result$C_Shannon, ")\n")
cat("Var_H_Shannon positive?  ", is.na(result$Var_H_Shannon) || result$Var_H_Shannon > 0, " (value:", result$Var_H_Shannon, ")\n")

###########################################################
# WIZ_OrdinalPatterns_D5_Hourly_Last6Months_2019.R
#
# Reads daily 100 Hz RSAM files (2019.001.WIZ.RSAM.csv ... 2019.365.WIZ.RSAM.csv)
# from an external drive, keeps only the LAST 6 MONTHS of 2019 (Jul-Dec, i.e.
# day-of-year >= 182), averages each second's 100 raw samples down to 1 Hz,
# computes hourly ordinal-pattern Shannon/Renyi/Tsallis entropy, complexity,
# Fisher information and variances (same math as test_one_hour_2019.R), and
# saves ALL hours from ALL selected days into a single combined CSV.
#
# Input folder (external drive):
#   E:/Seismic_Amplitude_Timeseries_Analysis/2019 secods RSAM
# Output folder (local):
#   C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results

rm(list = ls())

library(tidyverse)
library(lubridate)
library(StatOrdPattHxC)
library(ordinalpatterns)

# ── Parameters (identical to test_one_hour_2019.R) ─────────────────────────
D       <- 5
tau     <- 1
z_alpha <- qnorm(0.975)
BETA    <- 1.5
var_col <- "rsam"

# ── Folders ─────────────────────────────────────────────────────────────────
input_dir  <- "E:/Seismic_Amplitude_Timeseries_Analysis/2019 secods RSAM"
output_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Find daily files and keep only the last 6 months (day-of-year >= 182) ──
# Filenames: 2019.DDD.WIZ.RSAM.csv  (DDD = 001-365, 1-indexed day of year)
# 2019 is not a leap year: Jan(31)+Feb(28)+Mar(31)+Apr(30)+May(31)+Jun(30) = 181
# -> day 182 = 1 July 2019, so DDD >= 182 gives Jul-Dec inclusive.
all_files <- list.files(input_dir, pattern = "^2019\\.\\d{3}\\.WIZ\\.RSAM\\.csv$",
                        full.names = TRUE)

file_info <- tibble::tibble(path = all_files) %>%
  dplyr::mutate(
    fname = basename(path),
    doy   = as.integer(stringr::str_extract(fname, "(?<=^2019\\.)\\d{3}"))
  ) %>%
  dplyr::filter(doy >= 182) %>%
  dplyr::arrange(doy)

cat("Total daily files found on drive:", length(all_files), "\n")
cat("Files selected (last 6 months, DDD >= 182):", nrow(file_info), "\n")
if (nrow(file_info) > 0) {
  cat("First file:", file_info$fname[1], " | Last file:", file_info$fname[nrow(file_info)], "\n\n")
}
stopifnot(nrow(file_info) > 0)

# ── Per-hour feature computation (unchanged from test_one_hour_2019.R) ─────
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
  
  # Collapse 100 Hz raw samples to 1 Hz averages (same approach validated in
  # test_one_hour_2019.R -- this is what brought runtime down to something
  # feasible for a full 6-month run).
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

# ── Run across all selected days, with progress + timing ───────────────────
n_files <- nrow(file_info)
all_results <- vector("list", n_files)

t_start <- Sys.time()

for (i in seq_len(n_files)) {
  fname <- file_info$fname[i]
  path  <- file_info$path[i]
  
  cat(sprintf("[%d/%d] Processing %s ... ", i, n_files, fname))
  t_day <- system.time({
    all_results[[i]] <- process_day(path, fname)
  })
  cat(sprintf("%.1f sec, %d hours\n", t_day["elapsed"],
              ifelse(is.null(all_results[[i]]), 0L, nrow(all_results[[i]]))))
}

t_end <- Sys.time()
cat(sprintf("\nAll files processed in %.1f minutes.\n",
            as.numeric(difftime(t_end, t_start, units = "mins"))))

# ── Combine and save ─────────────────────────────────────────────────────
results <- dplyr::bind_rows(all_results) %>%
  dplyr::arrange(Hour_Start)

cat("Total hourly rows in combined result:", nrow(results), "\n")

# Unique filename down to the second, so re-runs never silently overwrite
# a previous result file.
timestamp_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_file <- file.path(output_dir,
                      paste0("WIZ_OrdinalPatterns_D5_Hourly_Last6Months_2019_",
                             timestamp_tag, ".csv"))

readr::write_csv(results, out_file)

cat("Saved combined results to:\n", out_file, "\n")