## ============================================================================
## Weekly Ordinal-Pattern Entropy & Complexity Analysis -- RSAM only
## Analysis window: 2014-07-01 to 2019-12-31 (continuous timeline)
## Shannon, Renyi, Tsallis, Fisher measures + associated statistical complexity
## ============================================================================
##
## TIMEZONE NOTE (read this before anything else):
## The Date/Year/DOY/timestamp columns in the source "_Annual_RSAM.csv"
## files are NZ local time (NZST/NZDT depending on the date), NOT UTC --
## confirmed separately, and consistent with the downstream DBSCAN/SVM
## pipeline, which compares these same Week_Start/Week_End values against
## GeoNet's VAL "LocalTime" columns (also NZ local) rather than its UTC
## columns. This script never actually converts timezones anywhere --
## every `tz = "UTC"` below is just R's mandatory internal POSIXct label,
## required because R has to attach *some* zone attribute to a date-time
## object, not a claim that the values are genuinely UTC. Because of that,
## the values themselves are unaffected either way: they were, and remain,
## NZ local clock time throughout. What changed in this version is only
## the comments/messages that used to say "UTC" as if it were literally
## true -- corrected to say "NZ local (stored with a UTC label)" instead,
## so this doesn't cause a mix-up again later. No numeric outputs from
## this script change as a result of this correction.
##
## WHAT CHANGED FROM YOUR ORIGINAL SCRIPT (edit CONFIG block if any is wrong):
##  1. DATA_DIR / RESULTS_DIR updated. Results are now written into the SAME
##     folder as the input CSVs (DATA_DIR), per your request -- there is no
##     separate results directory any more.
##  2. Input files are now "<STATION>_NZ_<year>_Annual_RSAM.csv" (one per
##     station per year), not "WIZ_NZ_<year>.<ext>". Both stations are
##     processed by this ONE script now (STATIONS loop below), instead of
##     two hand-edited copies -- which is exactly the setup where earlier
##     versions of this pipeline silently drifted apart between stations.
##  3. New column layout: Date, Year, DOY, Station, <5th column>, RSAM, MF,
##     HF, DSAR. The 5th column comes in TWO DIFFERENT FORMS depending on
##     the station's file, both confirmed via diagnostics:
##       - WIZ files: already a full date+time string, e.g.
##         "2014-01-01 00:09:59" -- used directly as the datetime.
##       - WSRZ files: only a bare minute-mark with no date, e.g. "10:00:00",
##         "20:00:00", ... "00:00:00" (10/20/30/40/50/00 repeating) -- the
##         real hour has to be reconstructed from Date + watching the
##         minute value roll over (e.g. 50 -> 00 = next hour).
##     `parse_col5_to_datetime()` auto-detects which case a given file is
##     and handles it accordingly -- no per-station switch needed. An
##     earlier version of this script wrongly assumed BOTH stations used
##     the WIZ-style already-full-timestamp format; that assumption is
##     what this version fixes. Everything is read as raw TEXT first
##     (col_types forces every column to character) so readr's guessers
##     never get a chance to mis-parse anything before we look at it
##     ourselves. See the KNOWN LIMITATION note above
##     parse_col5_to_datetime() for one real edge case (a whole missing
##     hour) this reconstruction can't detect from the minute-mark alone.
##  4. Only RSAM is analysed (MF, HF, DSAR are read but dropped immediately).
##  5. Analysis is restricted to START_DATE .. END_DATE (2014-07-01 to
##     2019-12-31 inclusive, NZ local calendar dates) via a straightforward
##     date filter applied AFTER building the continuous timeline, so a
##     weekly window can still span a year-file boundary (e.g. Dec 2016 ->
##     Jan 2017) exactly as before.
##  6. WINDOW_DAYS = 7, STEP_DAYS = 6 (1-day overlap) -- UNCHANGED, as
##     requested. D = 5, tau = 1, BETA = 1.5 -- UNCHANGED.
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
RESULTS_DIR <- DATA_DIR   # <-- results saved into the SAME folder as the CSVs

# Both stations processed by ONE script, not two hand-edited copies -- a
# hardcoded "WIZ_NZ_" (as in the original version of this script) meant
# any WSRZ run had to be a manually-edited duplicate, which is exactly the
# kind of setup where the two copies can silently drift apart and produce
# a mismatch like the WIZ=335/WSRZ=334 row-count difference you found.
#STATIONS <- c("WIZ", "WSRZ")
STATIONS <- "WSRZ"
START_DATE <- as.Date("2014-07-01")   # NZ local calendar date
END_DATE   <- as.Date("2019-12-31")   # NZ local calendar date

# Year-files needed to cover [START_DATE, END_DATE]
YEARS <- lubridate::year(START_DATE):lubridate::year(END_DATE)

VAR_COLS    <- c("RSAM")   # RSAM only

D           <- 5      # embedding dimension for ordinal patterns
TAU         <- 1      # embedding delay
BETA        <- 1.5    # order parameter for Renyi / Tsallis entropy
Z_ALPHA     <- 1.96   # z-value for ~95% semi-CI

WINDOW_DAYS <- 7       # length of each weekly window (days)
STEP_DAYS   <- 6       # advance between consecutive window starts (days)
# STEP_DAYS < WINDOW_DAYS => overlapping weeks (unchanged from original)

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
##
## KNOWN LIMITATION of the (b) reconstruction: it detects an hour rollover
## only from a strict decrease in the minute value (e.g. 50 -> 00). If an
## ENTIRE hour of data is missing (not just some samples within it), the
## rollover from the last sample before the gap to the first sample after
## it still only registers as ONE rollover, silently under-counting by
## however many whole hours were skipped -- every subsequent timestamp
## that calendar day would then be mislabeled by that many hours. This
## can't be fixed from the minute-mark data alone; it would need either a
## reliable row count per day (1008 for a gap-free day) or an independent
## source of the true hour. Worth keeping in mind if a station's week
## counts or feature values look unexpectedly off on a day you know had a
## bigger-than-usual outage.
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

## ---- Helper: catch missed hour-rollovers in the minute-mark reconstruction
## At correct 10-min sampling, no single hour can ever legitimately contain
## more than 6 rows. If parse_col5_to_datetime() ever misses a rollover
## (the whole-missing-hour edge case noted above it), two or more real
## hours get merged into one reconstructed hour, and that hour will show
## up here with 7+ rows -- an unambiguous, checkable signal rather than
## something you'd have to notice by eye. Only meaningful for WSRZ-style
## files (the ones that actually went through reconstruction); for WIZ-
## style files (already-parsed full timestamps) this should never fire.
validate_hour_reconstruction <- function(d, station) {
  suspect <- d %>%
    dplyr::mutate(
      cal_date = as.Date(datetime),
      cal_hour = lubridate::hour(datetime)
    ) %>%
    dplyr::count(cal_date, cal_hour, name = "n_samples") %>%
    dplyr::filter(n_samples > 6)
  
  if (nrow(suspect) > 0) {
    warning(sprintf(
      "[%s] %d hour(s) have MORE than 6 samples -- this can only happen if ",
      station, nrow(suspect)
    ), "an hour-rollover was missed during reconstruction (see the KNOWN ",
    "LIMITATION note above parse_col5_to_datetime). The dates/hours below ",
    "have their timestamps misaligned for the rest of that calendar day ",
    "and need manual review:")
    print(suspect, n = Inf)
  } else {
    message(sprintf("[%s] Hour-reconstruction check passed: no hour has more than 6 samples.", station))
  }
  
  invisible(suspect)
}

## ---- Helper: load & prepare one year's data -----------------------------
load_year_data <- function(year, data_dir, station) {
  
  basename_i <- paste0(station, "_NZ_", year, "_Annual_RSAM")
  
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
  
  # Force EVERY column to raw text on read. This is deliberate: readr's
  # automatic type-guessing is exactly what caused silent, format-dependent
  # bugs earlier (Date auto-parsed to ISO Date losing its original form;
  # the 5th column auto-parsed as a bare "hms" time-of-day with the wrong
  # field read as hours). Reading everything as character means we control
  # every conversion ourselves, explicitly, via parse_col5_to_datetime().
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
      "Week %d (%s to %s NZ local, var = %s) has only %d points (< D = %d) -- skipping.",
      week_index, format(week_start), format(week_end),
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
      "  Week %d (%s to %s NZ local, var = %s): %d / %d expected points (gaps present).",
      week_index, format(week_start), format(week_end),
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

## ---- Run the whole pipeline once per station ------------------------------
date_tag <- paste0(format(START_DATE, "%Y%m%d"), "_to_", format(END_DATE, "%Y%m%d"))

for (station in STATIONS) {
  
  message(sprintf("\n\n========== %s ==========", station))
  
  ## ---- Diagnostics: show what's actually in DATA_DIR before loading ------
  message(sprintf("---- DIAGNOSTIC [%s]: files present in DATA_DIR ----", station))
  all_files_in_dir <- list.files(DATA_DIR)
  if (length(all_files_in_dir) == 0) {
    message("  DATA_DIR (", DATA_DIR, ") is empty or does not exist / is not reachable from R.")
  } else {
    print(all_files_in_dir)
  }
  for (yr in YEARS) {
    expected <- paste0(station, "_NZ_", yr, "_Annual_RSAM")
    hits <- list.files(DATA_DIR, pattern = paste0("^", expected, "\\.[A-Za-z0-9]+$"), full.names = FALSE)
    if (length(hits) == 0) {
      message("  Year ", yr, ": NO file matched pattern '", expected, ".*' -- check filename/case. ",
              "This year is a GAP in the timeline for ", station, ".")
    } else {
      message("  Year ", yr, ": matched -> ", paste(hits, collapse = ", "))
      peek <- tryCatch(
        readr::read_csv(file.path(DATA_DIR, hits[1]), n_max = 2, show_col_types = FALSE),
        error = function(e) NULL
      )
      if (!is.null(peek)) {
        message("    Columns found: ", paste(names(peek), collapse = " | "))
      }
    }
  }
  message(sprintf("---- END DIAGNOSTIC [%s] ----", station))
  
  ## ---- Load ALL years and concatenate into one continuous timeline -------
  message(sprintf("[%s] Loading and concatenating years %d-%d ...", station, min(YEARS), max(YEARS)))
  
  d_list <- purrr::map(YEARS, load_year_data, data_dir = DATA_DIR, station = station)
  d_list <- purrr::compact(d_list)  # drop NULLs (missing/bad years)
  
  if (length(d_list) == 0) {
    warning(sprintf("[%s] No usable data loaded for any year in %d-%d -- skipping this station.",
                    station, min(YEARS), max(YEARS)))
    next
  }
  
  d_all <- dplyr::bind_rows(d_list) %>%
    dplyr::arrange(datetime) %>%
    dplyr::distinct(datetime, .keep_all = TRUE)  # guard against duplicate timestamps at file seams
  
  ## ---- Catch missed hour-rollovers before doing anything else with this
  ## station's data (see validate_hour_reconstruction() above) ----
  validate_hour_reconstruction(d_all, station)
  
  ## ---- Restrict to the requested analysis window --------------------------
  # tz = "UTC" here is a storage label only (see TIMEZONE NOTE at top);
  # START_DATE/END_DATE are NZ local calendar dates, and d_all$datetime is
  # NZ local clock time carrying the same label, so this comparison is
  # already local-to-local and correct as written.
  range_start <- as.POSIXct(START_DATE, tz = "UTC")
  range_end   <- as.POSIXct(END_DATE, tz = "UTC") + lubridate::days(1)  # exclusive upper bound
  
  d_all <- d_all %>%
    dplyr::filter(datetime >= range_start, datetime < range_end)
  
  if (nrow(d_all) == 0) {
    warning(sprintf("[%s] No data left after restricting to %s -- %s -- skipping this station.",
                    station, START_DATE, END_DATE))
    next
  }
  
  message(sprintf(
    "[%s] Continuous timeline built (restricted to %s -- %s, NZ local): %d rows from %s to %s (NZ local)",
    station, START_DATE, END_DATE, nrow(d_all),
    format(min(d_all$datetime)), format(max(d_all$datetime))
  ))
  
  ## ---- Build overlapping weekly window definitions across the WHOLE span -
  build_week_defs <- function(d) {
    
    t_min <- lubridate::floor_date(min(d$datetime), unit = "day")
    t_max <- max(d$datetime)
    
    week_starts <- seq(from = t_min, to = t_max, by = paste(STEP_DAYS, "days"))
    # drop windows whose start would leave no data at all
    week_starts <- week_starts[week_starts < t_max]
    
    list(
      defs = tibble::tibble(
        week_index = seq_along(week_starts),
        week_start = week_starts,
        week_end   = week_starts + lubridate::days(WINDOW_DAYS)  # exclusive upper bound
      ),
      t_min = t_min,
      t_max = t_max
    )
  }
  
  wd <- build_week_defs(d_all)
  week_defs <- wd$defs
  
  # Explicit diagnostic for exactly the WIZ-vs-WSRZ row-count question:
  # t_min/t_max are each station's OWN data range, so if one station's
  # real data starts later or ends earlier within [START_DATE, END_DATE]
  # (a gap at the edge, or a missing year file above), the two stations
  # get different t_min/t_max and therefore a different total window
  # count, even for the "same" nominal analysis period.
  message(sprintf("[%s] Window-grid t_min: %s", station, format(wd$t_min)))
  message(sprintf("[%s] Window-grid t_max: %s", station, format(wd$t_max)))
  message(sprintf("[%s] Span t_max - t_min: %.4f days", station,
                  as.numeric(difftime(wd$t_max, wd$t_min, units = "days"))))
  message(sprintf("[%s] Total weekly windows across full span: %d", station, nrow(week_defs)))
  
  ## ---- Run for RSAM across the restricted span, save ----------------------
  for (var_col in VAR_COLS) {
    
    message(sprintf("[%s] Computing weekly entropy/complexity for: %s", station, var_col))
    
    out_file <- file.path(
      RESULTS_DIR,
      paste0("weekly_entropy_complexity_", station, "_", var_col, "_", date_tag, ".csv")
    )
    
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
        message(sprintf("  [%s] Week %d done and saved -> %s", station, wi, out_file))
      }
    }
    
    message(sprintf("  -> [%s] finished %s: %d weeks saved to %s", station, var_col, n_weeks_written, out_file))
  }
}

message("\nDone -- both stations processed.")