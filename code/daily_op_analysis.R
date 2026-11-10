## ============================================================================
## Daily Ordinal-Pattern Entropy & Complexity Analysis -- RSAM only
## Analysis window: 2014-07-01 to 2019-12-31 (continuous timeline)
## Shannon, Renyi, Tsallis, Fisher measures + associated statistical complexity
## ============================================================================
##
## HOW THIS DIFFERS FROM THE WEEKLY SCRIPT:
##  1. D = 3 (embedding dimension), not 5, per your request.
##  2. Windows are DAILY, not weekly, with a 1-HOUR overlap (not 1-day) into
##     the previous day:
##       - Day 1 (2014-07-01 itself): [day_start, day_start + 1 day) ->
##         144 points (10-min sampling), NO overlap, since there's nothing
##         before START_DATE to overlap with.
##       - Every subsequent day i: [day_start(i) - 1 hour, day_start(i) + 1
##         day) -> 25 hours = 150 points: the last 6 points (1 hour) of the
##         PREVIOUS day, plus all 144 points of day i itself.
##     This matches your spec exactly: day 1 has 144 observations, every
##     other day has 150 (144 + 6 overlapping with the tail of the day
##     before it).
##  3. Data loading, the 5th-column datetime auto-detection (full-timestamp
##     vs. bare minute-mark), and the date-range restriction are UNCHANGED
##     from the validated weekly script -- copied verbatim.
##  4. Only RSAM is analysed. Default STATION_PREFIX is "WIZ_NZ" per your
##     request to start with WIZ; change the one CONFIG line to switch
##     stations later (same as the weekly script).
## ============================================================================


## ---- Make warnings appear immediately (not buffered to the end) ----------
options(warn = 1)

## ---- Packages --------------------------------------------------------------
library(dplyr)
library(tibble)
library(purrr)
library(lubridate)
library(readr)
library(stringr)
library(ordinalpatterns)
library(StatOrdPattHxC)

## ---- CONFIG ------------------------------------------------------------
DATA_DIR    <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"
RESULTS_DIR <- DATA_DIR   # results saved into the SAME folder as the CSVs

# Which station's annual files to use, e.g. files named
# "<STATION_PREFIX>_<year>_Annual_RSAM.csv". Change this one line to switch
# stations (e.g. "WIZ_NZ" <-> "WSRZ_NZ").
STATION_PREFIX <- "WIZ_NZ"

START_DATE <- as.Date("2014-07-01")
END_DATE   <- as.Date("2019-12-31")

# Year-files needed to cover [START_DATE, END_DATE]
YEARS <- lubridate::year(START_DATE):lubridate::year(END_DATE)

VAR_COLS    <- c("RSAM")   # RSAM only

D           <- 3      # embedding dimension for ordinal patterns (DAILY: 3)
TAU         <- 1      # embedding delay
BETA        <- 1.5    # order parameter for Renyi / Tsallis entropy
Z_ALPHA     <- 1.96   # z-value for ~95% semi-CI

OVERLAP_HOURS <- 1     # hours of overlap pulled from the tail of the
# PREVIOUS day into each day (except day 1)

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

## ---- Helper: parse Date column robustly (used only for the minute-mark
## fallback case below -- WIZ-style full timestamps don't need this at all)
parse_flexible_date <- function(x) {
  x_chr <- as.character(x)
  out <- suppressWarnings(lubridate::mdy(x_chr))       # matches "1/1/2014"
  if (all(is.na(out))) out <- suppressWarnings(lubridate::ymd(x_chr))  # "2014-01-01"
  if (all(is.na(out))) out <- suppressWarnings(lubridate::dmy(x_chr))  # last resort
  out
}

## ---- Helper: auto-detect and parse the 5th column into a real datetime --
## Two known formats, both confirmed via diagnostics on the real files:
##   (a) WIZ-style: already a full date+time string, e.g. "2014-01-01 00:09:59"
##   (b) WSRZ-style: a bare minute-mark with no date, e.g. "10:00:00",
##       "20:00:00", ..., "00:00:00" -- despite the HH:MM:SS *look*, the
##       leading number is actually just the MINUTE within the hour
##       (10/20/30/40/50/00 repeating); the hour has to be reconstructed by
##       watching that value roll over (e.g. 50 -> 00 = next hour), combined
##       with the Date column.
## Both `date_chr` and `col5_chr` must be raw, unparsed TEXT (see the
## col_types = cols(.default = col_character()) read below) so readr's own
## guessers never get a chance to mis-parse either one first.
parse_col5_to_datetime <- function(date_chr, col5_chr) {
  
  # --- Try (a): does col5 parse as a full date+time on its own? ---
  full_dt <- suppressWarnings(lubridate::parse_date_time(
    col5_chr, orders = c("ymd HMS", "ymd HM", "mdy HMS", "mdy HM"),
    tz = "UTC", quiet = TRUE
  ))
  frac_full <- mean(!is.na(full_dt))
  
  if (frac_full > 0.5) {
    if (frac_full < 1) {
      warning(sum(is.na(full_dt)), " row(s) failed full-timestamp parsing ",
              "out of ", length(full_dt), " -- these rows will be dropped.")
    }
    return(as.POSIXct(full_dt, tz = "UTC"))
  }
  
  # --- Fall back to (b): bare minute-mark + roll-over reconstruction ---
  date_parsed <- parse_flexible_date(date_chr)
  if (all(is.na(date_parsed))) {
    warning("Could not parse the Date column in either M/D/Y or Y-M-D form ",
            "-- cannot reconstruct datetimes for this file.")
    return(as.POSIXct(rep(NA, length(col5_chr)), tz = "UTC"))
  }
  
  leading_num  <- stringr::str_extract(col5_chr, "^\\-?\\d+(\\.\\d+)?")
  mins_in_hour <- suppressWarnings(as.numeric(leading_num))
  
  bad <- is.na(mins_in_hour) | mins_in_hour < 0 | mins_in_hour >= 60
  if (any(bad)) {
    warning(sum(bad), " row(s) had an unparseable/out-of-range minute-mark ",
            "value (examples: ", paste(utils::head(unique(col5_chr[bad]), 3), collapse = ", "),
            ") -- these rows will be dropped.")
    mins_in_hour[bad] <- NA_real_
  }
  
  secs_in_hour <- mins_in_hour * 60
  
  tibble::tibble(date_parsed, secs_in_hour) %>%
    dplyr::mutate(.row = dplyr::row_number()) %>%
    dplyr::group_by(date_parsed) %>%
    dplyr::arrange(.row, .by_group = TRUE) %>%
    dplyr::mutate(
      dropped_val = dplyr::coalesce(secs_in_hour < dplyr::lag(secs_in_hour), FALSE),
      hour_counter = cumsum(dropped_val),
      total_secs_of_day = hour_counter * 3600 + secs_in_hour
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.row) %>%
    dplyr::mutate(
      datetime = as.POSIXct(date_parsed, tz = "UTC") + total_secs_of_day
    ) %>%
    dplyr::pull(datetime)
}

## ---- Helper: load & prepare one year's data -----------------------------
load_year_data <- function(year, data_dir, station_prefix) {
  
  basename_i <- paste0(station_prefix, "_", year, "_Annual_RSAM")
  
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
  
  # Force EVERY column to raw text on read -- see weekly script for why
  # (readr's auto-guessing silently mis-parses both Date and the 5th column
  # in format-dependent ways otherwise).
  raw <- readr::read_csv(data_path,
                         col_types = readr::cols(.default = readr::col_character()),
                         show_col_types = FALSE)
  
  if (ncol(raw) < 5) {
    message("Comma-delimited read produced ", ncol(raw),
            " column(s) -- retrying as tab-delimited.")
    raw <- readr::read_delim(data_path, delim = "\t",
                             col_types = readr::cols(.default = readr::col_character()),
                             show_col_types = FALSE)
  }
  if (ncol(raw) < 5) {
    message("Tab-delimited read produced ", ncol(raw),
            " column(s) -- retrying as whitespace-delimited.")
    raw <- readr::read_table(data_path,
                             col_types = readr::cols(.default = readr::col_character()),
                             show_col_types = FALSE)
  }
  
  if (ncol(raw) < 5) {
    warning("File for year ", year, " has only ", ncol(raw),
            " column(s) -- expected at least 5 (Date, Year, DOY, Station, ",
            "timestamp, ...). Skipping.")
    return(NULL)
  }
  
  required_cols <- c("Date", "RSAM")
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
  
  col5_chr <- raw[[5]]  # by position -- header name/blank-ness varies by file
  
  d <- raw %>%
    dplyr::mutate(
      datetime = parse_col5_to_datetime(Date, col5_chr),
      RSAM     = as.numeric(RSAM)
    ) %>%
    dplyr::filter(!is.na(datetime)) %>%
    dplyr::filter(lubridate::year(datetime) == year) %>%
    dplyr::select(datetime, RSAM) %>%
    dplyr::arrange(datetime)
  
  if (nrow(d) == 0) {
    warning("No rows left after filtering file for year ", year, " to that year -- skipping.")
    return(NULL)
  }
  
  d
}

## ---- Diagnostics: show what's actually in DATA_DIR before loading -------
message("---- DIAGNOSTIC: files present in DATA_DIR ----")
all_files_in_dir <- list.files(DATA_DIR)
if (length(all_files_in_dir) == 0) {
  message("  DATA_DIR (", DATA_DIR, ") is empty or does not exist / is not reachable from R.")
} else {
  print(all_files_in_dir)
}
for (yr in YEARS) {
  expected <- paste0(STATION_PREFIX, "_", yr, "_Annual_RSAM")
  hits <- list.files(DATA_DIR, pattern = paste0("^", expected, "\\.[A-Za-z0-9]+$"), full.names = FALSE)
  if (length(hits) == 0) {
    message("  Year ", yr, ": NO file matched pattern '", expected, ".*' -- check filename/case.")
  } else {
    message("  Year ", yr, ": matched -> ", paste(hits, collapse = ", "))
  }
}
message("---- END DIAGNOSTIC ----")

## ---- Load ALL years and concatenate into one continuous timeline --------
message("Loading and concatenating years ", min(YEARS), "-", max(YEARS),
        " for station ", STATION_PREFIX, " ...")

d_list <- purrr::map(YEARS, load_year_data, data_dir = DATA_DIR, station_prefix = STATION_PREFIX)
d_list <- purrr::compact(d_list)  # drop NULLs (missing/bad years)

if (length(d_list) == 0) {
  stop("No usable data loaded for any year in ", min(YEARS), "-", max(YEARS), ".")
}

d_all <- dplyr::bind_rows(d_list) %>%
  dplyr::arrange(datetime) %>%
  dplyr::distinct(datetime, .keep_all = TRUE)  # guard against duplicate timestamps at file seams

## ---- Restrict to the requested analysis window ---------------------------
range_start <- as.POSIXct(START_DATE, tz = "UTC")
range_end   <- as.POSIXct(END_DATE, tz = "UTC") + lubridate::days(1)  # exclusive upper bound

d_all <- d_all %>%
  dplyr::filter(datetime >= range_start, datetime < range_end)

if (nrow(d_all) == 0) {
  stop("No data left after restricting to ", START_DATE, " -- ", END_DATE, ".")
}

message("Continuous timeline built (restricted to ", START_DATE, " -- ", END_DATE, "): ",
        nrow(d_all), " rows from ",
        format(min(d_all$datetime), tz = "UTC"), " to ",
        format(max(d_all$datetime), tz = "UTC"))

## ---- Build daily window definitions, 1-hour overlap (except day 1) ------
## Day 1: [START_DATE, START_DATE + 1 day)                    -> 144 points
## Day i>1: [day_start(i) - 1 hour, day_start(i) + 1 day)      -> 150 points
build_day_defs <- function(start_date, end_date, overlap_hours) {
  
  day_starts <- seq(as.POSIXct(start_date, tz = "UTC"),
                    as.POSIXct(end_date,   tz = "UTC"),
                    by = "day")
  
  n <- length(day_starts)
  
  window_start <- day_starts
  if (n > 1) {
    window_start[2:n] <- day_starts[2:n] - lubridate::hours(overlap_hours)
  }
  window_end <- day_starts + lubridate::days(1)
  
  tibble::tibble(
    day_index    = seq_len(n),
    nominal_date = as.Date(day_starts),
    window_start = window_start,
    window_end   = window_end
  )
}

day_defs <- build_day_defs(START_DATE, END_DATE, OVERLAP_HOURS)
message("Total daily windows across full span: ", nrow(day_defs))

## ---- Per-day feature computation -----------------------------------------
compute_day <- function(day_index, nominal_date, window_start, window_end, d, var_col) {
  
  sub_df <- d %>%
    dplyr::filter(datetime >= window_start, datetime < window_end)
  
  series <- sub_df[[var_col]]
  series <- series[is.finite(series)]
  n_i    <- length(series)
  
  n_pats <- factorial(D)
  
  if (n_i < D) {
    warning(sprintf(
      "Day %d (%s, %s to %s, var = %s) has only %d points (< D = %d) -- skipping.",
      day_index, nominal_date, format(window_start, tz = "UTC"), format(window_end, tz = "UTC"),
      var_col, n_i, D
    ))
    return(NULL)
  }
  
  n_eff <- n_i - D + 1
  
  # Expected point count: 144 for day 1 (no overlap), 150 for every day
  # after that (144 + 6 points/hour * OVERLAP_HOURS). Gaps in the raw data
  # will bring n_i below this -- not an error, just worth knowing about.
  expected_full_day <- if (day_index == 1) {
    (24 * 60) / 10
  } else {
    (24 * 60) / 10 + (OVERLAP_HOURS * 60) / 10
  }
  if (n_i < expected_full_day) {
    message(sprintf(
      "  Day %d (%s, var = %s): %d / %d expected points (gaps present).",
      day_index, nominal_date, var_col, n_i, expected_full_day
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
    Nominal_Date  = nominal_date,
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

## ---- Run for RSAM across the restricted span, save -----------------------
date_tag <- paste0(format(START_DATE, "%Y%m%d"), "_to_", format(END_DATE, "%Y%m%d"))

for (var_col in VAR_COLS) {
  
  message("Computing daily entropy/complexity for: ", var_col)
  
  out_file <- file.path(RESULTS_DIR, paste0("daily_entropy_complexity_", STATION_PREFIX, "_", var_col, "_", date_tag, ".csv"))
  
  # Start each run with a fresh file so re-runs don't append to stale results
  if (file.exists(out_file)) {
    file.remove(out_file)
  }
  
  n_days_written <- 0
  
  for (i in seq_len(nrow(day_defs))) {
    
    di <- day_defs$day_index[i]
    nd <- day_defs$nominal_date[i]
    ws <- day_defs$window_start[i]
    we <- day_defs$window_end[i]
    
    day_result <- compute_day(di, nd, ws, we, d = d_all, var_col = var_col)
    
    if (!is.null(day_result)) {
      day_result <- day_result %>%
        tibble::add_column(Year = lubridate::year(nd), .after = "Variable")
      
      readr::write_csv(
        day_result,
        out_file,
        append    = file.exists(out_file),
        col_names = !file.exists(out_file)
      )
      n_days_written <- n_days_written + 1
      message(sprintf("  Day %d (%s) done and saved -> %s", di, nd, out_file))
    }
  }
  
  message("  -> finished ", var_col, ": ", n_days_written, " days saved to ", out_file)
}

message("Done.")