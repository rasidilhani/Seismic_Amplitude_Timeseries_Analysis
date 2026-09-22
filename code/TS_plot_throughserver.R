## ============================================================================
## WIZ / WSRZ full-period and eruption-zoom RSAM time series plots
##
## Data comes straight from the server now, not the manually-concatenated
## annual CSVs. This reuses the exact sync approach proven in the OP-analysis
## scripts: one rsync call per daily file (a single batched --files-from
## request silently fails partway through on this server), run several at a
## time via a small bash/xargs loop script, with CRLF-safe file writing so
## WSL's rsync doesn't choke on Windows line endings.
##
## Server layout: mlvolc/seismic_amplitude_timeseries_out/<year>/<year>.<DOY>/
##                 <STATION>.NZ/<year>.<DOY>.<STATION>.timeseries.csv
## Daily file columns: an unnamed first column (full timestamp, fractional
## seconds), then RSAM, MF, HF, DSAR. Only RSAM is used here.
##
## BEFORE RUNNING:
##  - rsync must be reachable via WSL ('wsl' on PATH).
##  - Set the password as an environment variable, not in this file:
##        Sys.setenv(RSYNC_PASSWORD = "the-password-you-were-given")
## ============================================================================

options(warn = 1)

library(dplyr)
library(tibble)
library(purrr)
library(lubridate)
library(readr)
library(ggplot2)

## ---- CONFIG ----------------------------------------------------------------

STATIONS   <- c("WIZ", "WSRZ")
START_DATE <- as.Date("2014-07-01")
END_DATE   <- as.Date("2019-12-31")

RSYNC_HOST   <- "illslef1@core.geo.vuw.ac.nz"
RSYNC_MODULE <- "mlvolc/seismic_amplitude_timeseries_out"

## Scratch mirror of just the days needed - not a permanent copy. Safe to
## delete; it just gets re-synced (and any changed days re-fetched) next run.
CACHE_DIR   <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/daily_cache"
RESULTS_DIR <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"

TIME_AXIS_LABEL <- "UTC"   # everything here is UTC throughout - no local-time toggle

ANALYSIS_START <- as.POSIXct("2014-07-01 00:00:00", tz = "UTC")
ANALYSIS_END   <- as.POSIXct("2019-12-31 23:59:59", tz = "UTC")

station_colors <- c(WIZ = "#1B9E77", WSRZ = "#7570B3")

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

## rsync only exists inside WSL, not in Windows R directly, so every rsync
## call is handed off with system2("wsl", ...). WSL needs Linux-style paths
## ("/mnt/c/Users/...") for anything it touches; R keeps using normal Windows
## paths everywhere else (e.g. reading the synced CSVs back with read_csv()).
win_to_wsl_path <- function(path) {
  path  <- gsub("\\\\", "/", path)
  drive <- tolower(substr(path, 1, 1))
  rest  <- substr(path, 3, nchar(path))  # drop "C:"
  paste0("/mnt/", drive, rest)
}

## ---- Step 1: which daily files do we need, for BOTH stations? -------------

all_dates <- seq(START_DATE, END_DATE, by = "day")

rel_paths_for_station <- function(station, dates) {
  sprintf(
    "%d/%d.%03d/%s.NZ/%d.%03d.%s.timeseries.csv",
    lubridate::year(dates), lubridate::year(dates), lubridate::yday(dates),
    station,
    lubridate::year(dates), lubridate::yday(dates), station
  )
}

rel_paths <- unlist(lapply(STATIONS, rel_paths_for_station, dates = all_dates))

## Written inside CACHE_DIR (a real Windows folder WSL can also see via
## /mnt/c/...), opened in binary ("wb") mode so R does NOT translate "\n"
## into Windows "\r\n" - a stray trailing "\r" on every path makes every
## filename fail to match on the WSL/rsync side.
files_list_path <- file.path(CACHE_DIR, "timeseries_files_from.txt")
files_con <- file(files_list_path, open = "wb")
writeLines(rel_paths, files_con, sep = "\n")
close(files_con)

## ---- Step 2: rsync-sync those files into CACHE_DIR, ONE FILE AT A TIME ----
## A single batched --files-from request silently fails partway through on
## this server, so this loops over the file list and calls rsync separately
## for each file (like running `rsync source/one/file dest/` by hand),
## with NUM_PARALLEL running at once via xargs so the sync isn't too slow.
## Re-running only re-transfers files that actually changed on the server.

message("Syncing ", length(rel_paths), " daily file(s) across ", length(STATIONS),
        " station(s) - one rsync call per file.")

NUM_PARALLEL <- 8

wsl_cache_dir       <- win_to_wsl_path(CACHE_DIR)
wsl_files_list_path <- win_to_wsl_path(files_list_path)

sync_src_base <- sprintf("rsync://%s/%s", RSYNC_HOST, RSYNC_MODULE)

loop_script_path     <- file.path(CACHE_DIR, "timeseries_sync_loop.sh")
wsl_loop_script_path <- win_to_wsl_path(loop_script_path)

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
    "Could not run WSL - is Windows Subsystem for Linux installed and is ",
    "'wsl' on your PATH? (", conditionMessage(e), ")"
  )
)

n_ok   <- sum(grepl("^OK", rsync_result))
n_miss <- sum(grepl("^MISS", rsync_result))
message("Sync finished: ", n_ok, " transferred, ", n_miss, " missing/failed.")
if (length(rsync_result) > 0 && length(rsync_result) <= 40) {
  cat(rsync_result, sep = "\n")
}

n_synced <- length(list.files(CACHE_DIR, pattern = "\\.timeseries\\.csv$", recursive = TRUE))
if (n_synced == 0) {
  stop(
    "rsync did not download any files into CACHE_DIR (", CACHE_DIR, ").\n",
    "Scroll up to check the rsync output printed above for the actual error."
  )
}
message(n_synced, " daily file(s) now present in the local cache.")

## ---- Step 3: load one day's data for one station --------------------------

load_day_data <- function(date_i, station, cache_dir) {
  
  yr  <- lubridate::year(date_i)
  doy <- lubridate::yday(date_i)
  
  rel_path <- sprintf("%d/%d.%03d/%s.NZ/%d.%03d.%s.timeseries.csv",
                      yr, yr, doy, station, yr, doy, station)
  data_path <- file.path(cache_dir, rel_path)
  
  if (!file.exists(data_path)) {
    return(NULL)  # genuinely not on the server for this day - a real gap
  }
  
  raw <- tryCatch(
    readr::read_csv(data_path, show_col_types = FALSE),
    error = function(e) {
      warning("Failed to read ", data_path, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(raw) || nrow(raw) == 0) return(NULL)
  
  ## first column is unnamed in the file - it's the datetime, picked up BY
  ## POSITION and renamed here
  names(raw)[1] <- "datetime"
  
  if (!inherits(raw$datetime, "POSIXct")) {
    raw$datetime <- suppressWarnings(
      lubridate::ymd_hms(as.character(raw$datetime), tz = "UTC")
    )
  }
  
  if (!"RSAM" %in% names(raw)) {
    warning("File for ", date_i, " (", station, ") is missing an RSAM column - skipping.")
    return(NULL)
  }
  
  raw %>%
    dplyr::mutate(RSAM = as.numeric(RSAM)) %>%
    dplyr::filter(!is.na(datetime), !is.na(RSAM)) %>%
    dplyr::select(datetime, RSAM) %>%
    dplyr::arrange(datetime)
}

## ---- Step 4: load ALL days for one station into a continuous timeline -----

load_station <- function(station, dates, cache_dir) {
  
  message("Loading and concatenating ", length(dates), " day(s) for station ", station, " ...")
  
  d_list <- purrr::map(dates, load_day_data, station = station, cache_dir = cache_dir)
  n_missing <- sum(purrr::map_lgl(d_list, is.null))
  if (n_missing > 0) {
    message("  (", n_missing, " day(s) had no file on the server for ", station, " - real gaps.)")
  }
  d_list <- purrr::compact(d_list)
  
  if (length(d_list) == 0) {
    stop("No usable data loaded for ", station, " in ", START_DATE, " -- ", END_DATE, ".")
  }
  
  dplyr::bind_rows(d_list) %>%
    dplyr::arrange(datetime) %>%
    dplyr::distinct(datetime, .keep_all = TRUE) %>%   # guard duplicate timestamps at file seams
    dplyr::rename(utc_time = datetime) %>%
    dplyr::filter(utc_time >= ANALYSIS_START, utc_time <= ANALYSIS_END) %>%
    dplyr::mutate(station = station)
}

wiz  <- load_station("WIZ",  all_dates, CACHE_DIR)
wsrz <- load_station("WSRZ", all_dates, CACHE_DIR)

cat("WIZ range:  ", format(range(wiz$utc_time)),  "\n")
cat("WSRZ range: ", format(range(wsrz$utc_time)), "\n")

## ---- Step 5: eruption windows to highlight and zoom into -------------------
##
## Given as NZ local clock times, with offset_hours telling us how far ahead
## of UTC that clock time sits: 12 hours during NZST, 13 during NZDT. Since
## the RSAM data is UTC, every window is converted to UTC by subtracting
## its offset.

eruption_windows_raw <- tribble(
  ~label,                                       ~nz_start,              ~nz_end,                ~offset_hours,
  "28 Apr 2016 eruption",                       "2016-04-28 11:30:00",  "2016-04-28 18:44:59",   12,
  "13-15 Sep 2016 eruption",                    "2016-09-13 15:15:00",  "2016-09-15 12:14:59",   12,
  "9 Dec 2019 explosive eruption",              "2019-12-09 14:30:00",  "2019-12-09 16:24:59",   13,
  "9-12 Dec 2019 continued eruptive activity",  "2019-12-09 16:25:00",  "2019-12-12 10:19:59",   13
)

eruption_windows <- eruption_windows_raw %>%
  mutate(
    nz_start = as.POSIXct(nz_start, tz = "UTC"),   # tagged UTC only so R can do
    nz_end   = as.POSIXct(nz_end,   tz = "UTC"),   # arithmetic on it - the values
    # themselves are still the NZ local clock times written above
    start = nz_start - hours(offset_hours),
    end   = nz_end   - hours(offset_hours)
  ) %>%
  select(label, start, end)

## ---- Step 6: shared theme ---------------------------------------------------

theme_plain_serif <- theme_bw(base_size = 13, base_family = "serif") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(family = "serif"),
    axis.text  = element_text(family = "serif"),
    plot.title = element_text(family = "serif"),
    legend.position = "bottom"
  )

## ---- Step 7: full-period plots, one station each, no highlight -------------

plot_full_period <- function(d, station_name) {
  ggplot(d, aes(x = utc_time, y = RSAM)) +
    geom_line(color = station_colors[[station_name]], linewidth = 0.25) +
    labs(
      title = station_name,
      x = paste0("Date (", TIME_AXIS_LABEL, ")"), y = "RSAM (avg)"
    ) +
    theme_plain_serif
}

p_wiz_full  <- plot_full_period(wiz,  "WIZ")
p_wsrz_full <- plot_full_period(wsrz, "WSRZ")

## ---- Step 8: full-period plots with eruption windows highlighted -----------

plot_full_period_highlighted <- function(d, station_name, windows) {
  windows_valid <- windows %>% filter(!is.na(start), !is.na(end))
  
  p <- ggplot(d, aes(x = utc_time, y = RSAM))
  
  if (nrow(windows_valid) > 0) {
    p <- p + geom_rect(
      data = windows_valid, inherit.aes = FALSE,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      fill = "red", alpha = 0.25
    )
  }
  
  p +
    geom_line(color = station_colors[[station_name]], linewidth = 0.25) +
    labs(
      title = paste0(station_name, ": Eruptions Highlighted"),
      x = paste0("Date (", TIME_AXIS_LABEL, ")"), y = "RSAM (avg)"
    ) +
    theme_plain_serif
}

p_wiz_full_hl  <- plot_full_period_highlighted(wiz,  "WIZ",  eruption_windows)
p_wsrz_full_hl <- plot_full_period_highlighted(wsrz, "WSRZ", eruption_windows)

## ---- Display and save the four full-period plots ---------------------------

print(p_wiz_full)
print(p_wsrz_full)
print(p_wiz_full_hl)
print(p_wsrz_full_hl)

ggsave(file.path(RESULTS_DIR, "WIZ_2014to2019_full_rsam.pdf"),               p_wiz_full,     width = 12, height = 5)
ggsave(file.path(RESULTS_DIR, "WSRZ_2014to2019_full_rsam.pdf"),              p_wsrz_full,    width = 12, height = 5)
ggsave(file.path(RESULTS_DIR, "WIZ_2014to2019_full_rsam_highlighted.pdf"),   p_wiz_full_hl,  width = 12, height = 5)
ggsave(file.path(RESULTS_DIR, "WSRZ_2014to2019_full_rsam_highlighted.pdf"),  p_wsrz_full_hl, width = 12, height = 5)

## ---- Step 9: zoom plots, one per eruption window, both stations overlaid --

ZOOM_PAD_DAYS <- 4  # days of context shown before/after each eruption window

plot_zoom <- function(wiz_d, wsrz_d, label, win_start, win_end, pad_days = ZOOM_PAD_DAYS) {
  zoom_start <- win_start - days(pad_days)
  zoom_end   <- win_end   + days(pad_days)
  
  d_zoom <- bind_rows(
    wiz_d  %>% filter(utc_time >= zoom_start, utc_time <= zoom_end) %>% mutate(station = "WIZ"),
    wsrz_d %>% filter(utc_time >= zoom_start, utc_time <= zoom_end) %>% mutate(station = "WSRZ")
  )
  
  ggplot(d_zoom, aes(x = utc_time, y = RSAM, color = station)) +
    geom_rect(
      data = tibble(xmin = win_start, xmax = win_end),
      inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "red", alpha = 0.15
    ) +
    geom_line(linewidth = 0.4) +
    scale_color_manual(values = station_colors) +
    scale_x_datetime(date_labels = "%d %b\n%H:%M") +
    labs(
      title = label,
      x = paste0("Date/Time (", TIME_AXIS_LABEL, ")"), y = "RSAM (avg)", color = NULL
    ) +
    theme_plain_serif +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9))
}

eruption_windows_valid <- eruption_windows %>% filter(!is.na(start), !is.na(end))

if (nrow(eruption_windows_valid) < nrow(eruption_windows)) {
  warning(sprintf(
    "%d of %d eruption windows have no start/end date and were skipped -- fill in eruption_windows above.",
    nrow(eruption_windows) - nrow(eruption_windows_valid), nrow(eruption_windows)
  ))
}

zoom_plots <- pmap(
  eruption_windows_valid,
  function(label, start, end) plot_zoom(wiz, wsrz, label, start, end)
)
names(zoom_plots) <- eruption_windows_valid$label

walk(zoom_plots, print)

walk2(zoom_plots, names(zoom_plots), function(p, label) {
  fname <- paste0("zoom_", gsub("[^A-Za-z0-9]+", "_", label), ".pdf")
  ggsave(file.path(RESULTS_DIR, fname), p, width = 10, height = 5)
})

cat(sprintf("\nSaved %d full-period plots and %d zoom plots to %s\n",
            4, length(zoom_plots), RESULTS_DIR))