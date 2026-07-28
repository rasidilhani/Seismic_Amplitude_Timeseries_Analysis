## ============================================================================
## Weekly Ordinal-Pattern Entropy & Complexity Analysis -- ALL YEARS
## Shannon, Renyi, Tsallis, Fisher measures + associated statistical complexity
## ============================================================================
##
## CHANGES FROM THE 2019-ONLY VERSION:
##  - Loops over one input file per year: WIZ_NZ_2007, WIZ_NZ_2008, ...,
##    up to WIZ_NZ_2022 (edit YEARS below -- you said both "2007 to 2022"
##    and separately listed files ending at "...WIZ_NZ_2021"; I've defaulted
##    to 2007:2022 inclusive to match your stated range. If 2022 has no
##    file, that year is just skipped with a warning -- nothing breaks.)
##  - Each year's weekly windows are computed independently, since each
##    year lives in its own file (windows never span across year
##    boundaries).
##  - Still writes one incremental per-year, per-variable CSV exactly as
##    before (saved week-by-week, immediately) as a safety net.
##  - NEW: also builds ONE combined results workbook
##    (WIZ_weekly_entropy_complexity_<firstyear>_<lastyear>.xlsx) with two
##    sheets -- "rsam_avg" and "displacement_avg_m" -- each containing all
##    years stacked together (with a Year column). This workbook is
##    rewritten after every year finishes, so there's always an up-to-date
##    combined file on disk even if a later year fails.
## ============================================================================

library(dplyr)
library(tibble)
library(purrr)
library(lubridate)
library(readr)
library(ordinalpatterns)
library(StatOrdPattHxC)
library(openxlsx)   # install.packages("openxlsx") if not already installed

## ---- CONFIG ------------------------------------------------------------
DATA_DIR    <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"
RESULTS_DIR <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"

YEARS <- 2007:2022   # <-- EDIT if this range is wrong

VAR_COLS    <- c("displacement_avg_m", "rsam_avg")

D           <- 5      # embedding dimension for ordinal patterns
TAU         <- 1      # embedding delay
BETA        <- 1.5      # order parameter for Renyi / Tsallis entropy
Z_ALPHA     <- 1.96   # z-value for ~95% semi-CI

WINDOW_DAYS <- 7       # length of each weekly window (days)
STEP_DAYS   <- 6       # advance between consecutive window starts (days)

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

COMBINED_XLSX <- file.path(
  RESULTS_DIR,
  sprintf("WIZ_weekly_entropy_complexity_%d_%d.xlsx", min(YEARS), max(YEARS))
)

## ---- Helper: robust file reader (comma -> tab -> whitespace) --------------
read_year_file <- function(basename_year) {
  
  matched_files <- list.files(
    DATA_DIR,
    pattern = paste0("^", basename_year, "\\.[A-Za-z0-9]+$"),
    full.names = TRUE
  )
  
  if (length(matched_files) == 0) {
    warning("No file found for '", basename_year, "' in ", DATA_DIR, " -- skipping this year.")
    return(NULL)
  }
  if (length(matched_files) > 1) {
    message("Multiple files matched '", basename_year, "', using the first: ", matched_files[1])
  }
  
  data_path <- matched_files[1]
  message("Reading data from: ", data_path)
  
  raw <- readr::read_csv(data_path, col_types = readr::cols(.default = readr::col_guess()))
  
  if (ncol(raw) < 4) {
    message("Comma-delimited read produced ", ncol(raw), " column(s) -- retrying as tab-delimited.")
    raw <- readr::read_delim(data_path, delim = "\t",
                             col_types = readr::cols(.default = readr::col_guess()))
  }
  if (ncol(raw) < 4) {
    message("Tab-delimited read produced ", ncol(raw), " column(s) -- retrying as whitespace-delimited.")
    raw <- readr::read_table(data_path, col_types = readr::cols(.default = readr::col_guess()))
  }
  
  required_cols <- c("day", "unix_timestamp", "displacement_avg_m", "rsam_avg")
  missing_cols  <- setdiff(required_cols, names(raw))
  
  if (length(missing_cols) > 0) {
    warning(
      "File for '", basename_year, "' is missing column(s): ", paste(missing_cols, collapse = ", "),
      ". Columns found: ", paste(names(raw), collapse = ", "), " -- skipping this year."
    )
    return(NULL)
  }
  
  raw %>%
    dplyr::mutate(
      day                = as.numeric(day),
      unix_timestamp     = as.numeric(unix_timestamp),
      displacement_avg_m = as.numeric(displacement_avg_m),
      rsam_avg           = as.numeric(rsam_avg),
      datetime           = as.POSIXct(unix_timestamp, origin = "1970-01-01", tz = "UTC")
    ) %>%
    dplyr::arrange(datetime)
}

## ---- Per-week feature computation (same logic as before, + Year field) ---
compute_week <- function(week_index, week_start, week_end, d, var_col, year) {
  
  sub_df <- d %>% dplyr::filter(datetime >= week_start, datetime < week_end)
  series <- sub_df[[var_col]]
  series <- series[is.finite(series)]
  n_i    <- length(series)
  n_pats <- factorial(D)
  
  if (n_i < D) {
    warning(sprintf(
      "Year %d, Week %d (%s to %s, var = %s) has only %d points (< D = %d) -- skipping.",
      year, week_index, format(week_start, tz = "UTC"), format(week_end, tz = "UTC"),
      var_col, n_i, D
    ))
    return(NULL)
  }
  
  n_eff <- n_i - D + 1
  
  expected_full_week <- (WINDOW_DAYS * 24 * 60) / 10
  if (n_i < expected_full_week) {
    message(sprintf(
      "  Year %d, Week %d (%s to %s, var = %s): %d / %d expected points (gaps present).",
      year, week_index, format(week_start, tz = "UTC"), format(week_end, tz = "UTC"),
      var_col, n_i, expected_full_week
    ))
  }
  
  xpd   <- ordinalpatterns::op_pd(series, D = D, tau = TAU)
  probs <- xpd$probabilities
  
  if (length(probs) != n_pats) {
    warning(sprintf(
      "Year %d, Week %d (var = %s): op_pd() returned %d probabilities, expected %d -- investigate.",
      year, week_index, var_col, length(probs), n_pats
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
    Year        = year,
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

## ---- Main loop: every year x every variable, incremental + combined ------

combined_results <- setNames(vector("list", length(VAR_COLS)), VAR_COLS)
for (v in VAR_COLS) combined_results[[v]] <- list()

for (yr in YEARS) {
  
  message("=== Year ", yr, " ===")
  
  basename_year <- paste0("WIZ_NZ_", yr)
  d <- read_year_file(basename_year)
  
  if (is.null(d) || nrow(d) == 0) {
    message("  Skipping year ", yr, " (no usable data).")
    next
  }
  
  t_min <- lubridate::floor_date(min(d$datetime), unit = "day")
  t_max <- max(d$datetime)
  
  week_starts <- seq(from = t_min, to = t_max, by = paste(STEP_DAYS, "days"))
  week_starts <- week_starts[week_starts < t_max]
  
  week_defs <- tibble::tibble(
    week_index = seq_along(week_starts),
    week_start = week_starts,
    week_end   = week_starts + lubridate::days(WINDOW_DAYS)
  )
  
  for (var_col in VAR_COLS) {
    
    message("  Computing weekly entropy/complexity for: ", var_col)
    
    out_file <- file.path(RESULTS_DIR,
                          sprintf("weekly_entropy_complexity_%s_%d.csv", var_col, yr))
    if (file.exists(out_file)) file.remove(out_file)
    
    year_var_results <- list()
    
    for (i in seq_len(nrow(week_defs))) {
      
      wi <- week_defs$week_index[i]
      ws <- week_defs$week_start[i]
      we <- week_defs$week_end[i]
      
      week_result <- compute_week(wi, ws, we, d = d, var_col = var_col, year = yr)
      
      if (!is.null(week_result)) {
        readr::write_csv(
          week_result, out_file,
          append    = file.exists(out_file),
          col_names = !file.exists(out_file)
        )
        year_var_results[[length(year_var_results) + 1]] <- week_result
        message(sprintf("    Week %d done and saved -> %s", wi, out_file))
      }
    }
    
    if (length(year_var_results) > 0) {
      combined_results[[var_col]][[as.character(yr)]] <- dplyr::bind_rows(year_var_results)
    }
  }
  
  # Rebuild the combined multi-sheet workbook after every year, so there's
  # always an up-to-date combined file on disk even if a later year fails.
  sheets_out <- lapply(combined_results, function(yr_list) {
    if (length(yr_list) == 0) return(tibble::tibble())
    dplyr::bind_rows(yr_list)
  })
  openxlsx::write.xlsx(sheets_out, file = COMBINED_XLSX, overwrite = TRUE)
  message("  -> combined workbook updated: ", COMBINED_XLSX)
}

message("Done. Combined results saved to: ", COMBINED_XLSX)