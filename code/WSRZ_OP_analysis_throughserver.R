## ============================================================================
## Weekly Ordinal-Pattern Entropy & Complexity Analysis -- RSAM only -- WSRZ
## Analysis window: 2014-07-01 to 2019-12-31 (continuous timeline)
## Shannon, Renyi, Tsallis, Fisher measures + associated statistical complexity
## ============================================================================
##
## WHAT CHANGED FROM YOUR PREVIOUS SCRIPT:
##  1. Input is no longer the manually-concatenated annual CSVs. Instead, for
##     every day in [START_DATE, END_DATE] this script rsync-syncs exactly
##     that day's file straight from the server into a small local CACHE_DIR,
##     then reads it from there. There is no permanent "downloaded dataset"
##     to keep track of -- CACHE_DIR is just a scratch mirror that gets
##     re-synced (and any changed days re-fetched) every time you run this.
##  2. Server layout (from your directory listing):
##       mlvolc/seismic_amplitude_timeseries_out/<year>/<year>.<DOY>/WIZ.NZ/<year>.<DOY>.WIZ.timeseries.csv
##     DOY is zero-padded to 3 digits.
##  3. New column layout per daily file: an UNNAMED first column that is
##     already a full timestamp with fractional seconds (e.g.
##     "2018-01-01 00:00:00.003131"), then RSAM, MF, HF, DSAR. No more
##     Date/Year/DOY/Station columns, and no minute-mark reconstruction --
##     confirmed via a real sample file. The first column is picked up BY
##     POSITION and renamed "datetime".
##  4. Only RSAM is analysed (MF, HF, DSAR are read but dropped immediately).
##  5. WINDOW_DAYS = 7, STEP_DAYS = 6 (1-day overlap) -- UNCHANGED. D = 5,
##     tau = 1, BETA = 1.5 -- UNCHANGED.
##
## BEFORE RUNNING:
##  - rsync must be available on this machine. On Windows this typically
##    means WSL, Git for Windows' bundled rsync, or cwRsync -- plain Windows
##    R does not ship rsync itself.
##  - Set the server password as an environment variable rather than pasting
##    it into this file, e.g. in the R console before sourcing this script:
##        Sys.setenv(RSYNC_PASSWORD = "the-password-you-were-given")
##    (or put that line in your .Renviron so you don't retype it each time).
## ============================================================================


## ---- Make warnings appear immediately (not buffered to the end) ----------
options(warn = 1)

## ---- Packages --------------------------------------------------------------
library(dplyr)
library(tibble)
library(purrr)
library(lubridate)
library(readr)
library(ordinalpatterns)
library(StatOrdPattHxC)

## ---- CONFIG ------------------------------------------------------------
STATION      <- "WSRZ"

START_DATE   <- as.Date("2014-07-01")
END_DATE     <- as.Date("2019-12-31")

RSYNC_HOST   <- "illslef1@core.geo.vuw.ac.nz"
RSYNC_MODULE <- "mlvolc/seismic_amplitude_timeseries_out"

## Local scratch mirror of just the days you need -- NOT a permanent copy of
## the dataset. Safe to delete any time; it will just be re-synced next run.
CACHE_DIR    <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/daily_cache"

## Where the final results CSV lands on YOUR computer.
RESULTS_DIR  <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"

VAR_COLS     <- c("RSAM")   # RSAM only

D            <- 5      # embedding dimension for ordinal patterns
TAU          <- 1      # embedding delay
BETA         <- 1.5    # order parameter for Renyi / Tsallis entropy
Z_ALPHA      <- 1.96   # z-value for ~95% semi-CI

WINDOW_DAYS  <- 7       # length of each weekly window (days)
STEP_DAYS    <- 6       # advance between consecutive window starts (days)
# STEP_DAYS < WINDOW_DAYS => overlapping weeks (unchanged from original)

if (!dir.exists(CACHE_DIR))   dir.create(CACHE_DIR, recursive = TRUE)
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

RSYNC_PASSWORD <- Sys.getenv("RSYNC_PASSWORD")
if (RSYNC_PASSWORD == "") {
  stop(
    "RSYNC_PASSWORD is not set. Before running this script, do:\n",
    "  Sys.setenv(RSYNC_PASSWORD = \"the-password-you-were-given\")\n",
    "then re-run."
  )
}

## rsync itself only exists inside your WSL/Ubuntu environment, not in
## Windows R directly, so every rsync call below is handed off to WSL with
## system2("wsl", ...). WSL needs Linux-style paths (e.g. "/mnt/c/Users/...")
## rather than Windows ones ("C:/Users/...") for anything it touches --
## R keeps using the normal Windows path everywhere else (e.g. when it reads
## the synced CSVs back with read_csv()). This helper does that conversion.
win_to_wsl_path <- function(path) {
  path  <- gsub("\\\\", "/", path)
  drive <- tolower(substr(path, 1, 1))
  rest  <- substr(path, 3, nchar(path))  # drop "C:"
  paste0("/mnt/", drive, rest)
}

## ---- Step 1: which daily files do we need for [START_DATE, END_DATE]? ----
all_dates <- seq(START_DATE, END_DATE, by = "day")

rel_paths <- sprintf(
  "%d/%d.%03d/%s.NZ/%d.%03d.%s.timeseries.csv",
  lubridate::year(all_dates), lubridate::year(all_dates), lubridate::yday(all_dates),
  STATION,
  lubridate::year(all_dates), lubridate::yday(all_dates), STATION
)

## Written inside CACHE_DIR (a real Windows folder that WSL can also see via
## /mnt/c/...) rather than tempdir(), so the WSL-side rsync process can find it.
## IMPORTANT: opened in binary ("wb") mode so R does NOT translate "\n" into
## Windows-style "\r\n" -- rsync's --files-from on the Linux/WSL side expects
## plain "\n"-terminated lines, and a stray trailing "\r" on every path would
## make every single filename fail to match (which is what silently produced
## the "0 files synced" result).
files_list_path <- file.path(CACHE_DIR, paste0(STATION, "_files_from.txt"))
files_con <- file(files_list_path, open = "wb")
writeLines(rel_paths, files_con, sep = "\n")
close(files_con)

## ---- Step 2: rsync-sync those files into CACHE_DIR, ONE FILE AT A TIME ---
## Testing against this specific server showed that a single batched
## --files-from request silently fails partway through: rsync correctly
## reports each file's size, then the transfer itself "vanishes" for every
## file. Requesting files one at a time -- exactly like running
## `rsync source/one/file dest/` by hand -- works reliably. So instead of one
## batched request, this builds a small bash script that loops over the file
## list and calls rsync separately for each file, with a modest number
## (NUM_PARALLEL) running at once so the whole sync doesn't take forever.
## Re-running only re-transfers files that actually changed on the server --
## rsync's own size/mtime comparison still applies to each individual call.
message("Syncing ", length(rel_paths), " daily file(s) for station ", STATION,
        " -- one rsync call per file, ", "this is slower than a single batch ",
        "sync but avoids a server-side issue with batched requests.")

NUM_PARALLEL <- 8

wsl_cache_dir       <- win_to_wsl_path(CACHE_DIR)
wsl_files_list_path <- win_to_wsl_path(files_list_path)

sync_src_base <- sprintf("rsync://%s/%s", RSYNC_HOST, RSYNC_MODULE)

loop_script_path     <- file.path(CACHE_DIR, paste0(STATION, "_sync_loop.sh"))
wsl_loop_script_path <- win_to_wsl_path(loop_script_path)

## Written as its own little bash script (rather than one long system2() call)
## so it can loop and run several rsync calls in parallel via xargs -- opened
## in binary ("wb") mode for the same CRLF reason as files_list_path above.
script_lines <- c(
  "#!/bin/bash",
  "sync_one() {",
  "  relpath=\"$1\"",
  "  destdir=\"$DEST_BASE/$(dirname \"$relpath\")\"",
  "  mkdir -p \"$destdir\"",
  "  if rsync -a \"$SRC_BASE/$relpath\" \"$destdir/\" >/dev/null 2>&1; then",
  "    echo \"OK   $relpath\"",
  "  else",
  "    echo \"MISS $relpath\"",
  "  fi",
  "}",
  "export -f sync_one",
  sprintf("export SRC_BASE=%s", shQuote(sync_src_base, type = "sh")),
  sprintf("export DEST_BASE=%s", shQuote(wsl_cache_dir, type = "sh")),
  sprintf("export RSYNC_PASSWORD=%s", shQuote(RSYNC_PASSWORD, type = "sh")),
  sprintf(
    "cat %s | xargs -P %d -I {} bash -c 'sync_one \"$@\"' _ {}",
    shQuote(wsl_files_list_path, type = "sh"), NUM_PARALLEL
  )
)

script_con <- file(loop_script_path, open = "wb")
writeLines(script_lines, script_con, sep = "\n")
close(script_con)

rsync_result <- tryCatch(
  system2("wsl", c("bash", wsl_loop_script_path), stdout = TRUE, stderr = TRUE),
  error = function(e) stop(
    "Could not run WSL -- is Windows Subsystem for Linux installed and is ",
    "'wsl' on your PATH? (", conditionMessage(e), ")"
  )
)

n_ok   <- sum(grepl("^OK", rsync_result))
n_miss <- sum(grepl("^MISS", rsync_result))
message("Sync finished: ", n_ok, " transferred, ", n_miss, " missing/failed.")
if (length(rsync_result) > 0 && length(rsync_result) <= 40) {
  cat(rsync_result, sep = "\n")  # print full detail only for short runs
}

## Fail fast with a clear message here rather than discovering much later,
## confusingly, that "no usable data loaded" -- if nothing actually landed in
## CACHE_DIR, something above (auth, connection, path) is wrong, and the
## output printed just above is the place to look.
n_synced <- length(list.files(CACHE_DIR, pattern = "\\.timeseries\\.csv$", recursive = TRUE))
if (n_synced == 0) {
  stop(
    "rsync did not download any files into CACHE_DIR (", CACHE_DIR, ").\n",
    "Scroll up to check the rsync output printed above for the actual error."
  )
}
message(n_synced, " daily file(s) now present in the local cache.")

## ---- Step 3: load & prepare one day's data --------------------------------
load_day_data <- function(date_i, station, cache_dir) {
  
  yr  <- lubridate::year(date_i)
  doy <- lubridate::yday(date_i)
  
  rel_path <- sprintf("%d/%d.%03d/%s.NZ/%d.%03d.%s.timeseries.csv",
                      yr, yr, doy, station, yr, doy, station)
  data_path <- file.path(cache_dir, rel_path)
  
  if (!file.exists(data_path)) {
    return(NULL)  # genuinely not on the server for this day -- a real gap
  }
  
  raw <- tryCatch(
    readr::read_csv(data_path, show_col_types = FALSE),
    error = function(e) {
      warning("Failed to read ", data_path, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(raw) || nrow(raw) == 0) return(NULL)
  
  ## The first column has no header in the file -- it's the datetime,
  ## picked up BY POSITION and renamed here.
  names(raw)[1] <- "datetime"
  
  if (!inherits(raw$datetime, "POSIXct")) {
    raw$datetime <- suppressWarnings(
      lubridate::ymd_hms(as.character(raw$datetime), tz = "UTC")
    )
  }
  
  if (!"RSAM" %in% names(raw)) {
    warning("File for ", date_i, " (", station, ") is missing an RSAM column -- skipping.")
    return(NULL)
  }
  
  raw %>%
    dplyr::mutate(RSAM = as.numeric(RSAM)) %>%
    dplyr::filter(!is.na(datetime), !is.na(RSAM)) %>%
    dplyr::select(datetime, RSAM) %>%
    dplyr::arrange(datetime)
}

## ---- Step 4: load ALL days and concatenate into one continuous timeline ---
message("Loading and concatenating ", length(all_dates), " day(s) for station ", STATION, " ...")

d_list <- purrr::map(all_dates, load_day_data, station = STATION, cache_dir = CACHE_DIR)
n_missing <- sum(purrr::map_lgl(d_list, is.null))
if (n_missing > 0) {
  message("  (", n_missing, " day(s) had no file on the server for ", STATION, " -- real gaps.)")
}
d_list <- purrr::compact(d_list)

if (length(d_list) == 0) {
  stop("No usable data loaded for ", STATION, " in ", START_DATE, " -- ", END_DATE, ".")
}

d_all <- dplyr::bind_rows(d_list) %>%
  dplyr::arrange(datetime) %>%
  dplyr::distinct(datetime, .keep_all = TRUE)  # guard against duplicate timestamps at file seams

message("Continuous timeline built: ", nrow(d_all), " rows from ",
        format(min(d_all$datetime), tz = "UTC"), " to ",
        format(max(d_all$datetime), tz = "UTC"))

## ---- Step 5: build overlapping weekly window definitions across the span --
build_week_defs <- function(d) {
  
  t_min <- lubridate::floor_date(min(d$datetime), unit = "day")
  t_max <- max(d$datetime)
  
  week_starts <- seq(from = t_min, to = t_max, by = paste(STEP_DAYS, "days"))
  week_starts <- week_starts[week_starts < t_max]
  
  tibble::tibble(
    week_index = seq_along(week_starts),
    week_start = week_starts,
    week_end   = week_starts + lubridate::days(WINDOW_DAYS)  # exclusive upper bound
  )
}

week_defs <- build_week_defs(d_all)
message("Total weekly windows across full span: ", nrow(week_defs))

## ---- Step 6: per-week feature computation (unchanged from before) --------
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

## ---- Step 7: run for RSAM across the span, save to YOUR computer ----------
date_tag <- paste0(format(START_DATE, "%Y%m%d"), "_to_", format(END_DATE, "%Y%m%d"))

for (var_col in VAR_COLS) {
  
  message("Computing weekly entropy/complexity for: ", var_col, " (", STATION, ")")
  
  out_file <- file.path(RESULTS_DIR, paste0("weekly_entropy_complexity_", STATION, "_", var_col, "_", date_tag, ".csv"))
  
  if (file.exists(out_file)) file.remove(out_file)
  
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