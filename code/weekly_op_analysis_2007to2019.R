## ============================================================================
## Weekly Ordinal-Pattern Entropy & Complexity Analysis -- RSAM only, 2007-2019
## Shannon, Renyi, Tsallis, Fisher measures + associated statistical complexity
## Continuous timeline: weekly windows are built ONCE across the full
## 2007-2019 span, so a window CAN straddle a Dec 31 -> Jan 1 file boundary.
## ============================================================================
##
## ASSUMPTIONS MADE (edit the CONFIG block below if any of these are wrong):
##  1. unix_timestamp is standard UTC time (as you specified) -> converted
##     with tz = "UTC", NOT Pacific/Auckland.
##  2. Weekly windows are 7 days long, and each new window starts 6 days
##     after the previous one -> a 1-day overlap between consecutive weeks.
##     This reproduces your example: week1 = day 1-7, week2 = day 7-13,
##     week3 = day 13-19. Change WINDOW_DAYS / STEP_DAYS below if you meant
##     something else.
##  3. One input file per year: WIZ_NZ_2007, WIZ_NZ_2008, ..., WIZ_NZ_2019.
##     Each file is read and filtered to rows whose datetime year matches
##     its OWN year (guards against stray rows from adjacent years within a
##     single file), then ALL years are row-bound into one continuous,
##     time-sorted dataframe BEFORE any windowing happens. This means a
##     weekly window is free to span from late December of one year into
##     January of the next, using real data from both files -- there is no
##     artificial per-year restart and no dead zone at year boundaries.
##  4. If a whole year's file is missing, that's a genuine data gap: the
##     continuous timeline simply has a hole there. Windows that fall
##     entirely inside the hole are skipped (too few points, per the
##     existing n_i < D guard); windows that partially overlap real data
##     on either side are computed on however many points actually exist
##     in the window (never padded/interpolated) -- this is the "if there's
##     no data, that's OK, missing data is fine" behavior you asked for.
##  5. ONLY rsam_avg is analysed. 
##  6. D (embedding dimension) = 5, tau (delay) = 1 as you specified. BETA
##     (Renyi/Tsallis order) kept at 1.5. WINDOW_DAYS = 7, STEP_DAYS = 6
##     (7-day overlapping windows), unchanged from before.
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
RESULTS_DIR <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results/Weekly_OP_Analysis_2008to2019"

# Years to process -> files WIZ_NZ_2007 ... WIZ_NZ_2019 inside DATA_DIR
YEARS <- 2008:2019

VAR_COLS    <- c("rsam_avg")   # RSAM only

D           <- 5      # embedding dimension for ordinal patterns
TAU         <- 1      # embedding delay
BETA        <- 1.5    # order parameter for Renyi / Tsallis entropy
Z_ALPHA     <- 1.96   # z-value for ~95% semi-CI

WINDOW_DAYS <- 7       # length of each weekly window (days)
STEP_DAYS   <- 6       # advance between consecutive window starts (days)
# STEP_DAYS < WINDOW_DAYS => overlapping weeks

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

## ---- Helper: load & prepare one year's data -----------------------------
load_year_data <- function(year, data_dir) {
  
  basename_i <- paste0("WIZ_NZ_", year)
  
  matched_files <- list.files(
    data_dir,
    pattern = paste0("^", basename_i, "\\.[A-Za-z0-9]+$"),
    full.names = TRUE
  )
  if (length(matched_files) == 0) {
    warning("No file found in ", data_dir, " matching basename '", basename_i,
            "' -- skipping year ", year, " (leaves a real gap in the timeline).")
    return(NULL)
  }
  if (length(matched_files) > 1) {
    message("Multiple files matched '", basename_i, "', using the first: ", matched_files[1])
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
  
  required_cols <- c("day", "unix_timestamp", "rsam_avg")
  missing_cols  <- setdiff(required_cols, names(raw))
  
  if (length(missing_cols) > 0) {
    warning(
      "File for year ", year, " is missing expected column(s): ",
      paste(missing_cols, collapse = ", "),
      ". Columns actually found: ", paste(names(raw), collapse = ", "),
      " -- skipping this year."
    )
    return(NULL)
  }
  
  raw <- raw %>%
    dplyr::mutate(
      day            = as.numeric(day),
      unix_timestamp = as.numeric(unix_timestamp),
      rsam_avg       = as.numeric(rsam_avg)
    )
  
  d <- raw %>%
    dplyr::mutate(
      datetime = as.POSIXct(unix_timestamp, origin = "1970-01-01", tz = "UTC")
    ) %>%
    dplyr::filter(lubridate::year(datetime) == year) %>%
    dplyr::select(datetime, rsam_avg)
  
  if (nrow(d) == 0) {
    warning("No rows left after filtering file for year ", year, " to that year -- skipping.")
    return(NULL)
  }
  
  d
}

## ---- Load ALL years and concatenate into one continuous timeline --------
message("Loading and concatenating years ", min(YEARS), "-", max(YEARS), " ...")

d_list <- purrr::map(YEARS, load_year_data, data_dir = DATA_DIR)
d_list <- purrr::compact(d_list)  # drop NULLs (missing/bad years)

if (length(d_list) == 0) {
  stop("No usable data loaded for any year in ", min(YEARS), "-", max(YEARS), ".")
}

d_all <- dplyr::bind_rows(d_list) %>%
  dplyr::arrange(datetime) %>%
  dplyr::distinct(datetime, .keep_all = TRUE)  # guard against duplicate timestamps at file seams

message("Continuous timeline built: ", nrow(d_all), " rows from ",
        format(min(d_all$datetime), tz = "UTC"), " to ",
        format(max(d_all$datetime), tz = "UTC"))

## ---- Build overlapping weekly window definitions across the WHOLE span --
build_week_defs <- function(d) {
  
  t_min <- lubridate::floor_date(min(d$datetime), unit = "day")
  t_max <- max(d$datetime)
  
  week_starts <- seq(from = t_min, to = t_max, by = paste(STEP_DAYS, "days"))
  # drop windows whose start would leave no data at all
  week_starts <- week_starts[week_starts < t_max]
  
  tibble::tibble(
    week_index = seq_along(week_starts),
    week_start = week_starts,
    week_end   = week_starts + lubridate::days(WINDOW_DAYS)  # exclusive upper bound
  )
}

week_defs <- build_week_defs(d_all)
message("Total weekly windows across full span: ", nrow(week_defs))

## ---- Per-week feature computation ---------------------------------------
compute_week <- function(week_index, week_start, week_end, d, var_col) {
  
  sub_df <- d %>%
    dplyr::filter(datetime >= week_start, datetime < week_end)
  
  series <- sub_df[[var_col]]
  series <- series[is.finite(series)]
  n_i    <- length(series)
  
  n_pats <- factorial(D)
  
  if (n_i < D) {
    warning(sprintf(
      "Week %d (%s to %s, var = %s) has only %d points (< D = %d) -- skipping.",
      week_index, format(week_start, tz = "UTC"), format(week_end, tz = "UTC"),
      var_col, n_i, D
    ))
    return(NULL)
  }
  
  n_eff <- n_i - D + 1
  
  # Sanity check: at 10-min sampling a full 7-day window should have 1008
  # points. Gaps/missing data (including real missing years) will bring n_i
  # below that -- not an error, just worth knowing about.
  expected_full_week <- (WINDOW_DAYS * 24 * 60) / 10
  if (n_i < expected_full_week) {
    message(sprintf(
      "  Week %d (%s to %s, var = %s): %d / %d expected points (gaps present).",
      week_index, format(week_start, tz = "UTC"), format(week_end, tz = "UTC"),
      var_col, n_i, expected_full_week
    ))
  }
  
  xpd   <- ordinalpatterns::op_pd(series, D = D, tau = TAU)
  probs <- xpd$probabilities
  
  if (length(probs) != n_pats) {
    warning(sprintf(
      "Week %d (var = %s): op_pd() returned %d probabilities, expected %d -- investigate.",
      week_index, var_col, length(probs), n_pats
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
    Variable    = var_col,
    Week_Index  = week_index,
    Week_Start  = week_start,
    Week_End    = week_end,
    N_points    = n_i,
    N_eff       = n_eff,
    
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

## ---- Run for RSAM across the full continuous 2007-2019 span, save -------
for (var_col in VAR_COLS) {
  
  message("Computing weekly entropy/complexity for: ", var_col)
  
  out_file <- file.path(RESULTS_DIR, paste0("weekly_entropy_complexity_", var_col, "_2007to2019.csv"))
  
  # Start each run with a fresh file so re-runs don't append to stale results
  if (file.exists(out_file)) {
    file.remove(out_file)
  }
  
  n_weeks_written <- 0
  
  for (i in seq_len(nrow(week_defs))) {
    
    wi <- week_defs$week_index[i]
    ws <- week_defs$week_start[i]
    we <- week_defs$week_end[i]
    
    week_result <- compute_week(wi, ws, we, d = d_all, var_col = var_col)
    
    if (!is.null(week_result)) {
      week_result <- week_result %>%
        tibble::add_column(Year = lubridate::year(ws), .after = "Variable")
      
      readr::write_csv(
        week_result,
        out_file,
        append    = file.exists(out_file),
        col_names = !file.exists(out_file)
      )
      n_weeks_written <- n_weeks_written + 1
      message(sprintf("  Week %d done and saved -> %s", wi, out_file))
    }
  }
  
  message("  -> finished ", var_col, ": ", n_weeks_written, " weeks saved to ", out_file)
}

message("Done.")