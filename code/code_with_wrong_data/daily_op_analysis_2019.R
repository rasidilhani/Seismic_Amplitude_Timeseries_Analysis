## ============================================================================
## Daily Ordinal-Pattern Entropy & Complexity Analysis
## Shannon, Renyi, Tsallis, Fisher measures + associated statistical complexity
## ============================================================================
##
## ASSUMPTIONS MADE (edit the CONFIG block below if any of these are wrong):
##  1. unix_timestamp is standard UTC time (as you specified) -> converted
##     with tz = "UTC", NOT Pacific/Auckland.
##  2. Windows are now DAILY (calendar-day, midnight-to-midnight in UTC),
##     with a 1-hour overlap pulled in from the END of the previous day.
##     At 10-min sampling that overlap is 6 records, so a full day's
##     window has 144 (day) + 6 (overlap) = 150 points. The very first
##     day in the data has no prior data to borrow from, so it will
##     naturally have fewer points (that's expected, not an error).
##  3. Only rows with year(datetime) == 2019 are used.
##  4. script loops over rsam_avg and writes result file
##  5. D (embedding dimension) = 3, tau (delay) = 1, BETA (Renyi/Tsallis order) = 1.5
##     -- adjust in CONFIG for your data.
## ============================================================================


## ---- Packages --------------------------------------------------------------
library(dplyr)
library(tibble)
library(purrr)
library(lubridate)
library(readr)
library(ordinalpatterns)
library(StatOrdPattHxC)

## ---- CONFIG ------------------------------------------------------------
DATA_DIR    <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"
RESULTS_DIR <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results/Daily_OP_Analysis_2019"

# Base name of your data file inside DATA_DIR (extension auto-detected below)
DATA_FILE_BASENAME <- "WIZ_NZ_2019"

VAR_COLS    <- c("rsam_avg")  # variables to analyse

D           <- 3     # embedding dimension for ordinal patterns
TAU         <- 1      # embedding delay
BETA        <- 1.5      # order parameter for Renyi / Tsallis entropy
Z_ALPHA     <- 1.96   # z-value for ~95% semi-CI

WINDOW_DAYS   <- 1     # length of each daily window (calendar days)
OVERLAP_HOURS <- 1     # hours borrowed from the END of the previous day
SAMPLE_MIN    <- 10    # sampling interval in minutes (used only for gap-checking messages)

YEAR_FILTER <- 2019

# Subfolder (under RESULTS_DIR) that results are saved into
OUTPUT_SUBDIR <- "Weekly_OP_Analysis_2019"
RESULTS_OUT_DIR <- file.path(RESULTS_DIR, OUTPUT_SUBDIR)

if (!dir.exists(RESULTS_OUT_DIR)) dir.create(RESULTS_OUT_DIR, recursive = TRUE)

## ---- Load & prepare data ----------------------------------------------
matched_files <- list.files(
  DATA_DIR,
  pattern = paste0("^", DATA_FILE_BASENAME, "\\.[A-Za-z0-9]+$"),
  full.names = TRUE
)
if (length(matched_files) == 0) {
  stop("No file found in ", DATA_DIR, " matching basename '", DATA_FILE_BASENAME, "'. ",
       "Check the filename/extension and DATA_DIR.")
}
if (length(matched_files) > 1) {
  message("Multiple files matched '", DATA_FILE_BASENAME, "', using the first: ", matched_files[1])
}
data_path <- matched_files[1]
message("Reading data from: ", data_path)

raw <- readr::read_csv(data_path, col_types = readr::cols(.default = readr::col_guess()))

# If comma-delimited reading didn't work, try tab, then whitespace.
if (ncol(raw) < 4) {
  message("Comma-delimited read produced ", ncol(raw),
          " column(s) -- retrying as tab-delimited.")
  raw <- readr::read_delim(data_path, delim = "\t",
                           col_types = readr::cols(.default = readr::col_guess()))
}
if (ncol(raw) < 4) {
  message("Tab-delimited read produced ", ncol(raw),
          " column(s) -- retrying as whitespace-delimited.")
  raw <- readr::read_table(data_path, col_types = readr::cols(.default = readr::col_guess()))
}

required_cols <- c("day", "unix_timestamp", "displacement_avg_m", "rsam_avg")
missing_cols  <- setdiff(required_cols, names(raw))

if (length(missing_cols) > 0) {
  stop(
    "Could not find expected column(s): ", paste(missing_cols, collapse = ", "),
    ".\nColumns actually found in the file: ", paste(names(raw), collapse = ", "),
    "\nOpen the file and check its delimiter/header row, then adjust the ",
    "read_delim()/read_table() call above to match."
  )
}

raw <- raw %>%
  dplyr::mutate(
    day                = as.numeric(day),
    unix_timestamp      = as.numeric(unix_timestamp),
    displacement_avg_m = as.numeric(displacement_avg_m),
    rsam_avg           = as.numeric(rsam_avg)
  )

d <- raw %>%
  dplyr::mutate(
    datetime = as.POSIXct(unix_timestamp, origin = "1970-01-01", tz = "UTC")
  ) %>%
  dplyr::filter(lubridate::year(datetime) == YEAR_FILTER) %>%
  dplyr::arrange(datetime)

if (nrow(d) == 0) stop("No rows left after filtering for year ", YEAR_FILTER, ".")

## ---- Build daily window definitions (with 1-hour overlap from the previous day) ----
t_min <- lubridate::floor_date(min(d$datetime), unit = "day")
t_max <- max(d$datetime)

# One window per calendar day, stepping by WINDOW_DAYS (= 1 day)
day_starts <- seq(from = t_min, to = t_max, by = paste(WINDOW_DAYS, "days"))
# drop a trailing day that would start on/after the last timestamp (no data)
day_starts <- day_starts[day_starts < t_max]

day_defs <- tibble::tibble(
  day_index    = seq_along(day_starts),
  day_start    = day_starts,                                             # calendar-day boundary (midnight)
  day_end      = day_starts + lubridate::days(WINDOW_DAYS),               # exclusive upper bound
  window_start = day_starts - lubridate::hours(OVERLAP_HOURS),            # pulls in last hour of previous day
  window_end   = day_starts + lubridate::days(WINDOW_DAYS)                # window end == day end (no forward overlap)
)

## ---- Per-day feature computation ---------------------------------------
compute_day <- function(day_index, day_start, day_end, window_start, window_end, d, var_col) {
  
  sub_df <- d %>%
    dplyr::filter(datetime >= window_start, datetime < window_end)
  
  series <- sub_df[[var_col]]
  series <- series[is.finite(series)]
  n_i    <- length(series)
  
  n_pats <- factorial(D)
  
  if (n_i < D) {
    warning(sprintf(
      "Day %d (%s to %s, var = %s) has only %d points (< D = %d) -- skipping.",
      day_index, format(day_start, tz = "UTC"), format(day_end, tz = "UTC"),
      var_col, n_i, D
    ))
    return(NULL)
  }
  
  n_eff <- n_i - D + 1
  
  # Sanity check: at SAMPLE_MIN-min sampling, a full day (WINDOW_DAYS) plus
  # the OVERLAP_HOURS borrowed from the previous day should have
  # (WINDOW_DAYS*24*60 + OVERLAP_HOURS*60) / SAMPLE_MIN points. The very
  # first day has no prior data to borrow, so it will fall short of this
  # by design -- not an error.
  expected_full_day <- (WINDOW_DAYS * 24 * 60 + OVERLAP_HOURS * 60) / SAMPLE_MIN
  if (n_i < expected_full_day) {
    message(sprintf(
      "  Day %d (%s to %s, var = %s): %d / %d expected points (gaps present, or first day with no prior overlap).",
      day_index, format(day_start, tz = "UTC"), format(day_end, tz = "UTC"),
      var_col, n_i, expected_full_day
    ))
  }
  
  xpd   <- ordinalpatterns::op_pd(series, D = D, tau = TAU)
  probs <- xpd$probabilities
  
  if (length(probs) != n_pats) {
    warning(sprintf(
      "Day %d (var = %s): op_pd() returned %d probabilities, expected %d -- investigate.",
      day_index, var_col, length(probs), n_pats
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
  
  semi <- function(v) ifelse(!is.finite(v) | v <= 0, NA, sqrt(v) / sqrt(n_eff) * Z_ALPHA)
  
  tibble::tibble(
    Variable      = var_col,
    Day_Index     = day_index,
    Day_Start     = day_start,
    Day_End       = day_end,
    Window_Start  = window_start,
    Window_End    = window_end,
    N_points      = n_i,
    N_eff         = n_eff,
    
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

## ---- Run for each variable and save (incremental, day-by-day) ---------
for (var_col in VAR_COLS) {
  
  message("Computing daily entropy/complexity for: ", var_col)
  
  out_file <- file.path(RESULTS_OUT_DIR, paste0("daily_entropy_complexity_", var_col, ".csv"))
  
  # Start each run with a fresh file so re-runs don't append to stale results
  if (file.exists(out_file)) {
    file.remove(out_file)
  }
  
  n_days_written <- 0
  
  for (i in seq_len(nrow(day_defs))) {
    
    di <- day_defs$day_index[i]
    ds <- day_defs$day_start[i]
    de <- day_defs$day_end[i]
    ws <- day_defs$window_start[i]
    we <- day_defs$window_end[i]
    
    day_result <- compute_day(di, ds, de, ws, we, d = d, var_col = var_col)
    
    if (!is.null(day_result)) {
      readr::write_csv(
        day_result,
        out_file,
        append    = file.exists(out_file),
        col_names = !file.exists(out_file)
      )
      n_days_written <- n_days_written + 1
      message(sprintf("  Day %d done and saved -> %s", di, out_file))
    }
  }
  
  message("  -> finished ", var_col, ": ", n_days_written, " days saved to ", out_file)
}

message("Done.")

############################################################################################
## ============================================================================
## Daily Ordinal-Pattern Entropy & Complexity Analysis
## Shannon, Renyi, Tsallis, Fisher measures + associated statistical complexity
## ============================================================================
##
## ASSUMPTIONS MADE (edit the CONFIG block below if any of these are wrong):
##  1. unix_timestamp is standard UTC time (as you specified) -> converted
##     with tz = "UTC", NOT Pacific/Auckland.
##  2. Windows are now DAILY (calendar-day, midnight-to-midnight in UTC),
##     with a 1-hour overlap pulled in from the END of the previous day.
##     At 10-min sampling that overlap is 6 records, so a full day's
##     window has 144 (day) + 6 (overlap) = 150 points. The very first
##     day in the data has no prior data to borrow from, so it will
##     naturally have fewer points (that's expected, not an error).
##  3. Only rows with year(datetime) == 2019 are used.
##  4. script loops over BOTH (Displacement_avg_m and rsam_avg) and writes two results file
##     per variable.
##  5. D (embedding dimension) = 3, tau (delay) = 1, BETA (Renyi/Tsallis order) = 1.5
##     -- adjust in CONFIG for your data.
## ============================================================================


## ---- Packages --------------------------------------------------------------
## ============================================================================
## Daily Ordinal-Pattern Entropy & Complexity Analysis
## Shannon, Renyi, Tsallis, Fisher measures + associated statistical complexity
## ============================================================================
##
## ASSUMPTIONS MADE (edit the CONFIG block below if any of these are wrong):
##  1. unix_timestamp is standard UTC time (as you specified) -> converted
##     with tz = "UTC", NOT Pacific/Auckland.
##  2. Windows are now DAILY (calendar-day, midnight-to-midnight in UTC),
##     with a 1-hour overlap pulled in from the END of the previous day.
##     At 10-min sampling that overlap is 6 records, so a full day's
##     window has 144 (day) + 6 (overlap) = 150 points. The very first
##     day in the data has no prior data to borrow from, so it will
##     naturally have fewer points (that's expected, not an error).
##  3. Rows from YEAR_START through YEAR_END (inclusive, currently 2007-2019)
##     are all loaded together and processed as ONE continuous run -- daily
##     windows span the whole 13-year period, not per-year, and all results
##     for all years are written to a single output CSV per variable.
##  4. VAR_COLS currently only includes rsam_avg (displacement dropped).
##  5. D (embedding dimension) = 3, tau (delay) = 1, BETA (Renyi/Tsallis order) = 1.5
##     -- adjust in CONFIG for your data.
## ============================================================================


## ---- Packages --------------------------------------------------------------
## ============================================================================
## Daily Ordinal-Pattern Entropy & Complexity Analysis
## Shannon, Renyi, Tsallis, Fisher measures + associated statistical complexity
## ============================================================================
##
## ASSUMPTIONS MADE (edit the CONFIG block below if any of these are wrong):
##  1. unix_timestamp is standard UTC time (as you specified) -> converted
##     with tz = "UTC", NOT Pacific/Auckland.
##  2. Windows are now DAILY (calendar-day, midnight-to-midnight in UTC),
##     with a 1-hour overlap pulled in from the END of the previous day.
##     At 10-min sampling that overlap is 6 records, so a full day's
##     window has 144 (day) + 6 (overlap) = 150 points. The very first
##     day in the data has no prior data to borrow from, so it will
##     naturally have fewer points (that's expected, not an error).
##  3. Rows from YEAR_START through YEAR_END (inclusive, currently 2007-2019)
##     are all loaded together and processed as ONE continuous run -- daily
##     windows span the whole 13-year period, not per-year, and all results
##     for all years are written to a single output CSV per variable.
##  4. VAR_COLS currently only includes rsam_avg (displacement dropped).
##  5. D (embedding dimension) = 3, tau (delay) = 1, BETA (Renyi/Tsallis order) = 1.5
##     -- adjust in CONFIG for your data.
## ============================================================================


## ---- Packages --------------------------------------------------------------
## ============================================================================
## Daily Ordinal-Pattern Entropy & Complexity Analysis
## Shannon, Renyi, Tsallis, Fisher measures + associated statistical complexity
## ============================================================================
##
## ASSUMPTIONS MADE (edit the CONFIG block below if any of these are wrong):
##  1. unix_timestamp is standard UTC time (as you specified) -> converted
##     with tz = "UTC", NOT Pacific/Auckland.
##  2. Windows are now DAILY (calendar-day, midnight-to-midnight in UTC),
##     with a 1-hour overlap pulled in from the END of the previous day.
##     At 10-min sampling that overlap is 6 records, so a full day's
##     window has 144 (day) + 6 (overlap) = 150 points. The very first
##     day in the data has no prior data to borrow from, so it will
##     naturally have fewer points (that's expected, not an error).
##  3. Rows from YEAR_START through YEAR_END (inclusive, currently 2007-2019)
##     are all loaded together and processed as ONE continuous run -- daily
##     windows span the whole 13-year period, not per-year, and all results
##     for all years are written to a single output CSV per variable.
##  4. VAR_COLS currently only includes rsam_avg (displacement dropped).
##  5. D (embedding dimension) = 3, tau (delay) = 1, BETA (Renyi/Tsallis order) = 1.5
##     -- adjust in CONFIG for your data.
## ============================================================================


## ---- Packages --------------------------------------------------------------
library(dplyr)
library(tibble)
library(purrr)
library(lubridate)
library(readr)
library(ordinalpatterns)
library(StatOrdPattHxC)

## ---- CONFIG ------------------------------------------------------------
DATA_DIR    <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"
RESULTS_DIR <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results/Daily_OP_Analysis_2019"

# Files are one-per-year: WIZ_NZ_<year>.csv (extension auto-detected below).
# All years from YEAR_START to YEAR_END are read and combined into one dataset.
DATA_FILE_PREFIX <- "WIZ_NZ_"

VAR_COLS    <- c("rsam_avg")  # variables to analyse

D           <- 3     # embedding dimension for ordinal patterns
TAU         <- 1      # embedding delay
BETA        <- 1.5      # order parameter for Renyi / Tsallis entropy
Z_ALPHA     <- 1.96   # z-value for ~95% semi-CI

WINDOW_DAYS   <- 1     # length of each daily window (calendar days)
OVERLAP_HOURS <- 1     # hours borrowed from the END of the previous day
SAMPLE_MIN    <- 10    # sampling interval in minutes (used only for gap-checking messages)

YEAR_START  <- 2007    # first year to include (inclusive)
YEAR_END    <- 2019    # last year to include (inclusive)

# Subfolder (under RESULTS_DIR) that results are saved into
OUTPUT_SUBDIR <- paste0("Daily_OP_Analysis_", YEAR_START, "_", YEAR_END)
RESULTS_OUT_DIR <- file.path(RESULTS_DIR, OUTPUT_SUBDIR)

if (!dir.exists(RESULTS_OUT_DIR)) dir.create(RESULTS_OUT_DIR, recursive = TRUE)

## ---- Load & prepare data (one file per year, combined) ------------------
read_one_file <- function(data_path) {
  message("Reading data from: ", data_path)
  
  raw <- readr::read_csv(data_path, col_types = readr::cols(.default = readr::col_guess()))
  
  # If comma-delimited reading didn't work, try tab, then whitespace.
  if (ncol(raw) < 4) {
    message("  Comma-delimited read produced ", ncol(raw),
            " column(s) -- retrying as tab-delimited.")
    raw <- readr::read_delim(data_path, delim = "\t",
                             col_types = readr::cols(.default = readr::col_guess()))
  }
  if (ncol(raw) < 4) {
    message("  Tab-delimited read produced ", ncol(raw),
            " column(s) -- retrying as whitespace-delimited.")
    raw <- readr::read_table(data_path, col_types = readr::cols(.default = readr::col_guess()))
  }
  
  required_cols <- c("day", "unix_timestamp", "displacement_avg_m", "rsam_avg")
  missing_cols  <- setdiff(required_cols, names(raw))
  
  if (length(missing_cols) > 0) {
    stop(
      "Could not find expected column(s) in ", data_path, ": ", paste(missing_cols, collapse = ", "),
      ".\nColumns actually found in the file: ", paste(names(raw), collapse = ", "),
      "\nOpen the file and check its delimiter/header row, then adjust the ",
      "read_delim()/read_table() call above to match."
    )
  }
  
  raw %>%
    dplyr::mutate(
      day                = as.numeric(day),
      unix_timestamp      = as.numeric(unix_timestamp),
      displacement_avg_m = as.numeric(displacement_avg_m),
      rsam_avg           = as.numeric(rsam_avg)
    ) %>%
    # Keep only the columns the analysis actually uses. Other columns
    # (e.g. datetime_nz) are dropped here because their type can vary
    # from file to file (e.g. parsed as text in one year's file, as a
    # datetime in another), which breaks bind_rows() when combining years.
    dplyr::select(dplyr::all_of(required_cols))
}

years_needed <- YEAR_START:YEAR_END
raw_list <- vector("list", length(years_needed))

for (i in seq_along(years_needed)) {
  yr <- years_needed[i]
  basename_yr <- paste0(DATA_FILE_PREFIX, yr)
  
  matched_files <- list.files(
    DATA_DIR,
    pattern = paste0("^", basename_yr, "\\.[A-Za-z0-9]+$"),
    full.names = TRUE
  )
  if (length(matched_files) == 0) {
    stop("No file found in ", DATA_DIR, " matching basename '", basename_yr, "'. ",
         "Check the filename/extension and DATA_DIR.")
  }
  if (length(matched_files) > 1) {
    message("Multiple files matched '", basename_yr, "', using the first: ", matched_files[1])
  }
  
  raw_list[[i]] <- read_one_file(matched_files[1])
}

raw <- dplyr::bind_rows(raw_list)

d <- raw %>%
  dplyr::mutate(
    datetime = as.POSIXct(unix_timestamp, origin = "1970-01-01", tz = "UTC")
  ) %>%
  dplyr::filter(lubridate::year(datetime) >= YEAR_START, lubridate::year(datetime) <= YEAR_END) %>%
  dplyr::arrange(datetime)

if (nrow(d) == 0) stop("No rows left after filtering for years ", YEAR_START, "-", YEAR_END, ".")

## ---- Build daily window definitions (with 1-hour overlap from the previous day) ----
t_min <- lubridate::floor_date(min(d$datetime), unit = "day")
t_max <- max(d$datetime)

# One window per calendar day, stepping by WINDOW_DAYS (= 1 day)
day_starts <- seq(from = t_min, to = t_max, by = paste(WINDOW_DAYS, "days"))
# drop a trailing day that would start on/after the last timestamp (no data)
day_starts <- day_starts[day_starts < t_max]

day_defs <- tibble::tibble(
  day_index    = seq_along(day_starts),
  day_start    = day_starts,                                             # calendar-day boundary (midnight)
  day_end      = day_starts + lubridate::days(WINDOW_DAYS),               # exclusive upper bound
  window_start = day_starts - lubridate::hours(OVERLAP_HOURS),            # pulls in last hour of previous day
  window_end   = day_starts + lubridate::days(WINDOW_DAYS)                # window end == day end (no forward overlap)
)

## ---- Per-day feature computation ---------------------------------------
compute_day <- function(day_index, day_start, day_end, window_start, window_end, d, var_col) {
  
  sub_df <- d %>%
    dplyr::filter(datetime >= window_start, datetime < window_end)
  
  series <- sub_df[[var_col]]
  series <- series[is.finite(series)]
  n_i    <- length(series)
  
  n_pats <- factorial(D)
  
  if (n_i < D) {
    warning(sprintf(
      "Day %d (%s to %s, var = %s) has only %d points (< D = %d) -- skipping.",
      day_index, format(day_start, tz = "UTC"), format(day_end, tz = "UTC"),
      var_col, n_i, D
    ))
    return(NULL)
  }
  
  n_eff <- n_i - D + 1
  
  # Sanity check: at SAMPLE_MIN-min sampling, a full day (WINDOW_DAYS) plus
  # the OVERLAP_HOURS borrowed from the previous day should have
  # (WINDOW_DAYS*24*60 + OVERLAP_HOURS*60) / SAMPLE_MIN points. The very
  # first day has no prior data to borrow, so it will fall short of this
  # by design -- not an error.
  expected_full_day <- (WINDOW_DAYS * 24 * 60 + OVERLAP_HOURS * 60) / SAMPLE_MIN
  if (n_i < expected_full_day) {
    message(sprintf(
      "  Day %d (%s to %s, var = %s): %d / %d expected points (gaps present, or first day with no prior overlap).",
      day_index, format(day_start, tz = "UTC"), format(day_end, tz = "UTC"),
      var_col, n_i, expected_full_day
    ))
  }
  
  xpd   <- ordinalpatterns::op_pd(series, D = D, tau = TAU)
  probs <- xpd$probabilities
  
  if (length(probs) != n_pats) {
    warning(sprintf(
      "Day %d (var = %s): op_pd() returned %d probabilities, expected %d -- investigate.",
      day_index, var_col, length(probs), n_pats
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
  
  semi <- function(v) ifelse(!is.finite(v) | v <= 0, NA, sqrt(v) / sqrt(n_eff) * Z_ALPHA)
  
  tibble::tibble(
    Variable      = var_col,
    Day_Index     = day_index,
    Day_Start     = day_start,
    Day_End       = day_end,
    Window_Start  = window_start,
    Window_End    = window_end,
    N_points      = n_i,
    N_eff         = n_eff,
    
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

## ---- Run for each variable and save (incremental, day-by-day) ---------
for (var_col in VAR_COLS) {
  
  message("Computing daily entropy/complexity for: ", var_col)
  
  out_file <- file.path(RESULTS_OUT_DIR, paste0("daily_entropy_complexity_", var_col, ".csv"))
  
  # Start each run with a fresh file so re-runs don't append to stale results
  if (file.exists(out_file)) {
    file.remove(out_file)
  }
  
  n_days_written <- 0
  
  for (i in seq_len(nrow(day_defs))) {
    
    di <- day_defs$day_index[i]
    ds <- day_defs$day_start[i]
    de <- day_defs$day_end[i]
    ws <- day_defs$window_start[i]
    we <- day_defs$window_end[i]
    
    day_result <- compute_day(di, ds, de, ws, we, d = d, var_col = var_col)
    
    if (!is.null(day_result)) {
      readr::write_csv(
        day_result,
        out_file,
        append    = file.exists(out_file),
        col_names = !file.exists(out_file)
      )
      n_days_written <- n_days_written + 1
      message(sprintf("  Day %d done and saved -> %s", di, out_file))
    }
  }
  
  message("  -> finished ", var_col, ": ", n_days_written, " days saved to ", out_file)
}

message("Done.")