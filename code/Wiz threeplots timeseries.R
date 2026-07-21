# WIZ_ThreePlots_TimeSeries.R
#
# Produces three time series plots:
#   Plot 1: 2007-2022, 10-minute data, split into 4 panels (~4 years each)
#   Plot 2: 2019 only, raw seconds-resolution data, NOT averaged to 1-sec RSAM
#   Plot 3: December 2019 only, hourly-averaged RSAM
#
# ASSUMPTIONS (check these against your actual data before trusting the output):
#   - Plot 1 files are named WIZ_NZ_<year>.csv, one file per year, 2007-2022
#   - Column names for each data source are set via time_col / value_col_annual /
#     value_col_seconds below, based on the actual column headers you provided
#   - Plot 2's raw data is large enough that it must be downsampled for
#     PLOTTING ONLY (see PLOT2_DOWNSAMPLE below) -- otherwise R/ggplot cannot
#     render ~3 billion points for a full year of 100Hz data

rm(list = ls())

library(tidyverse)
library(lubridate)

# ── Configure column names here to match your actual CSV files ─────────────
time_col <- "unix_timestamp"   # same in both data sources

# Plot 1 (annual 10-min data) has columns:
#   day  unix_timestamp  datetime_nz  displacement_avg_m  rsam_avg
value_col_annual <- "rsam_avg"

# Plots 2 & 3 (seconds data) have columns:
#   unix_timestamp  rsam
value_col_seconds <- "rsam"

# =============================================================================
# PLOT 1: 2007-2022, 10-minute data, 4 panels (~4 years each)
# =============================================================================

data_dir_annual <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"
years <- 2007:2022

annual_files <- tibble::tibble(
  year = years,
  path = file.path(data_dir_annual, paste0("WIZ_NZ_", years, ".csv"))
) %>%
  dplyr::mutate(exists = file.exists(path))

cat("Plot 1 -- files found:", sum(annual_files$exists), "of", nrow(annual_files), "\n")
if (any(!annual_files$exists)) {
  cat("WARNING: missing files for years:",
      paste(annual_files$year[!annual_files$exists], collapse = ", "), "\n")
}

# Read and combine all available yearly files
read_annual_file <- function(path, year) {
  d <- readr::read_csv(path, show_col_types = FALSE)
  if (!all(c(time_col, value_col_annual) %in% names(d))) {
    warning(sprintf("%s: expected columns '%s'/'%s' not found -- check time_col/value_col_annual.",
                    basename(path), time_col, value_col_annual))
    return(NULL)
  }
  d %>%
    dplyr::transmute(
      datetime = as.POSIXct(.data[[time_col]], origin = "1970-01-01", tz = "UTC"),
      value    = .data[[value_col_annual]],
      year     = year
    )
}

annual_data <- annual_files %>%
  dplyr::filter(exists) %>%
  purrr::pmap_dfr(function(year, path, exists) read_annual_file(path, year))

cat("Rows in annual_data:", nrow(annual_data), "\n")
cat("Columns:", paste(names(annual_data), collapse = ", "), "\n")
stopifnot("year" %in% names(annual_data), nrow(annual_data) > 0)

# Split into 4 panels of ~4 years each
panel_breaks <- c(2007, 2011, 2015, 2019, 2023)  # panel boundaries (exclusive upper)
panel_labels <- c("2007-2010", "2011-2014", "2015-2018", "2019-2022")

annual_data <- annual_data %>%
  dplyr::mutate(
    panel = cut(year, breaks = panel_breaks, labels = panel_labels,
                right = FALSE, include.lowest = TRUE)
  )

p1 <- ggplot(annual_data, aes(x = datetime, y = value)) +
  geom_line(linewidth = 0.3, color = "black") +
  facet_wrap(~ panel, scales = "free_x", ncol = 1) +
  labs(x = "Year", y = "RSAM") +
  theme_classic(base_family = "serif", base_size = 12) +
  theme(strip.text = element_blank(),
        strip.background = element_blank())

print(p1)
ggsave(file.path(data_dir_annual, "Plot1_2007_2022_4panels.pdf"),
       p1, width = 10, height = 12, dpi = 150)

# =============================================================================
# PLOT 2: 2019 only, 100Hz averaged to 1-second RSAM, full year time series
# =============================================================================

data_dir_seconds <- "E:/Seismic_Amplitude_Timeseries_Analysis/2019 secods RSAM"

seconds_files <- list.files(data_dir_seconds,
                            pattern = "^2019\\.\\d{3}\\.WIZ\\.RSAM\\.csv$",
                            full.names = TRUE)
cat("\nPlot 2 -- files found:", length(seconds_files), "of 365 expected\n")

# Averages 100Hz samples sharing the same unix_timestamp (i.e. same second)
# down to a single 1-second RSAM value -- this is a real aggregation, not a
# plotting shortcut, so the result is exact 1-second resolution data.
average_day_to_1sec <- function(path) {
  d <- readr::read_csv(path, show_col_types = FALSE)
  if (!all(c(time_col, value_col_seconds) %in% names(d))) {
    warning(sprintf("%s: expected columns '%s'/'%s' not found.",
                    basename(path), time_col, value_col_seconds))
    return(NULL)
  }
  d %>%
    dplyr::group_by(.data[[time_col]]) %>%
    dplyr::summarise(rsam_1s = mean(.data[[value_col_seconds]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::transmute(
      datetime = as.POSIXct(.data[[time_col]], origin = "1970-01-01", tz = "UTC"),
      value    = rsam_1s
    )
}

cat("Averaging", length(seconds_files), "files from 100Hz to 1-second RSAM...\n")

rsam_2019 <- purrr::map_dfr(seconds_files, average_day_to_1sec)

# Even at 1-second resolution, a full year is ~31.5 million points -- still
# too many for ggplot/most graphics devices to render smoothly. This does NOT
# re-average the data (it's already exact 1-second RSAM); it only thins how
# many of those points get DRAWN. Set to 1 to plot every single point if your
# machine can handle it.
PLOT2_DRAW_EVERY_NTH <- 10   # <-- adjust: higher = faster/lighter plot to render

rsam_2019_for_plot <- rsam_2019[seq(1, nrow(rsam_2019), by = PLOT2_DRAW_EVERY_NTH), ]

p2 <- ggplot(rsam_2019_for_plot, aes(x = datetime, y = value)) +
  geom_line(linewidth = 0.15, color = "darkred") +
  labs(x = "Date", y = expression(paste("RSAM (", mu, "m ", s^-1, ")")))+
  theme_set(theme_classic(base_family = "serif", base_size = 12))

print(p2)
ggsave(file.path(data_dir_seconds, "Plot2_2019_1sec_RSAM.pdf"),
       p2, width = 12, height = 6, dpi = 150)

# =============================================================================
# PLOT 3: December 2019 only, hourly-averaged RSAM
# =============================================================================

# =============================================================================
# PLOT 3: December 2019 only, raw 100Hz time series (NOT averaged to hourly)
# =============================================================================

# December 2019: day-of-year 335 (1 Dec) to 365 (31 Dec)
dec_doy <- 335:365

dec_files <- tibble::tibble(path = seconds_files) %>%
  dplyr::mutate(
    fname = basename(path),
    doy   = as.integer(stringr::str_extract(fname, "(?<=^2019\\.)\\d{3}"))
  ) %>%
  dplyr::filter(doy %in% dec_doy) %>%
  dplyr::arrange(doy)

cat("\nPlot 3 -- December files found:", nrow(dec_files), "of", length(dec_doy), "expected\n")

read_raw_day <- function(path) {
  d <- readr::read_csv(path, show_col_types = FALSE)
  if (!all(c(time_col, value_col_seconds) %in% names(d))) {
    warning(sprintf("%s: expected columns '%s'/'%s' not found.",
                    basename(path), time_col, value_col_seconds))
    return(NULL)
  }
  tibble::tibble(
    datetime = as.POSIXct(d[[time_col]], origin = "1970-01-01", tz = "UTC"),
    value    = d[[value_col_seconds]]
  )
}

cat("Reading raw 100Hz data for December 2019...\n")
dec_raw <- purrr::map_dfr(dec_files$path, read_raw_day)

# December at 100Hz is ~267 million points -- too many to render directly.
# This thins how many points get DRAWN only; the read-in data itself is raw,
# un-averaged 100Hz.
PLOT3_DRAW_EVERY_NTH <- 50   # <-- adjust: higher = faster/lighter plot to render

dec_raw_for_plot <- dec_raw[seq(1, nrow(dec_raw), by = PLOT3_DRAW_EVERY_NTH), ]

p3 <- ggplot(dec_raw_for_plot, aes(x = datetime, y = value)) +
  geom_line(linewidth = 0.15, color = "darkgreen") +
  labs(x = "Date", y = expression(paste("RSAM (", mu, "m ", s^-1, ")"))) +
  theme_classic(base_family = "serif", base_size = 12)

print(p3)
ggsave(file.path(data_dir_seconds, "Plot3_Dec2019_raw100Hz.pdf"),
       p3, width = 10, height = 6, dpi = 150)

# =============================================================================
# PLOT 4: December 8-9, 2019 only, raw 100Hz time series (NOT averaged)
# =============================================================================

dec8_9_doy <- c(342, 343)   # 8 Dec, 9 Dec 2019

dec8_9_files <- tibble::tibble(path = seconds_files) %>%
  dplyr::mutate(
    fname = basename(path),
    doy   = as.integer(stringr::str_extract(fname, "(?<=^2019\\.)\\d{3}"))
  ) %>%
  dplyr::filter(doy %in% dec8_9_doy) %>%
  dplyr::arrange(doy)

cat("\nPlot 4 -- Dec 8-9 files found:", nrow(dec8_9_files), "of 2 expected\n")

cat("Reading raw 100Hz data for Dec 8-9, 2019...\n")
dec8_9_raw <- purrr::map_dfr(dec8_9_files$path, read_raw_day)

# 2 days at 100Hz is ~17 million points -- still large, so keep a light
# downsample for drawing (much lower than Plot 3's, since the range is small).
PLOT4_DRAW_EVERY_NTH <- 5   # <-- adjust as needed

dec8_9_raw_for_plot <- dec8_9_raw[seq(1, nrow(dec8_9_raw), by = PLOT4_DRAW_EVERY_NTH), ]

p4 <- ggplot(dec8_9_raw_for_plot, aes(x = datetime, y = value)) +
  geom_line(linewidth = 0.2, color = "darkgreen") +
  scale_x_datetime(date_breaks = "3 hours", date_labels = "%d %H:%M") +
  labs(x = "Date/Hour", y = expression(paste("RSAM (", mu, "m ", s^-1, ")"))) +
  theme_classic(base_family = "serif", base_size = 12)

print(p4)
ggsave(file.path(data_dir_seconds, "Plot4_Dec8_9_2019_raw100Hz.png"),
       p4, width = 10, height = 6, dpi = 150)

# =============================================================================
# PLOT 4: December 8-9, 2019 only, hourly-averaged RSAM
# =============================================================================

dec8_9_doy <- c(342, 343)   # 8 Dec, 9 Dec 2019

dec8_9_files <- tibble::tibble(path = seconds_files) %>%
  dplyr::mutate(
    fname = basename(path),
    doy   = as.integer(stringr::str_extract(fname, "(?<=^2019\\.)\\d{3}"))
  ) %>%
  dplyr::filter(doy %in% dec8_9_doy) %>%
  dplyr::arrange(doy)

cat("\nPlot 4 -- Dec 8-9 files found:", nrow(dec8_9_files), "of 2 expected\n")

cat("Reading raw 100Hz data for Dec 8-9, 2019...\n")
dec8_9_raw <- purrr::map_dfr(dec8_9_files$path, read_raw_day)

# 2 days at 100Hz is ~17 million points -- still large, so keep a light
# downsample for drawing (much lower than Plot 3's, since the range is small).
PLOT4_DRAW_EVERY_NTH <- 5   # <-- adjust as needed

dec8_9_raw_for_plot <- dec8_9_raw[seq(1, nrow(dec8_9_raw), by = PLOT4_DRAW_EVERY_NTH), ]

p4 <- ggplot(dec8_9_raw_for_plot, aes(x = datetime, y = value)) +
  geom_line(linewidth = 0.2, color = "darkgreen") +
  scale_x_datetime(date_breaks = "3 hours", date_labels = "%d %H:%M") +
  labs(x = "Date/Hour", y = expression(paste("RSAM (", mu, "m ", s^-1, ")"))) +
  theme_classic(base_family = "serif", base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(p4)
ggsave(file.path(data_dir_seconds, "Plot4_Dec8_9_2019_raw100Hz.png"),
       p4, width = 10, height = 6, dpi = 150)
cat("\nAll three plots generated and saved.\n")

##########################################################
# average 100Hz to 1Hz
# =============================================================================
# PLOT 3: December 2019 only, 100Hz averaged to 1-second RSAM
# =============================================================================

# December 2019: day-of-year 335 (1 Dec) to 365 (31 Dec)
dec_doy <- 335:365

dec_files <- tibble::tibble(path = seconds_files) %>%
  dplyr::mutate(
    fname = basename(path),
    doy   = as.integer(stringr::str_extract(fname, "(?<=^2019\\.)\\d{3}"))
  ) %>%
  dplyr::filter(doy %in% dec_doy) %>%
  dplyr::arrange(doy)

cat("\nPlot 3 -- December files found:", nrow(dec_files), "of", length(dec_doy), "expected\n")

cat("Averaging December 2019 from 100Hz to 1-second RSAM...\n")
dec_1sec <- purrr::map_dfr(dec_files$path, average_day_to_1sec)

# December at 1-second resolution is ~2.7 million points -- still thin for
# rendering only; the averaged data itself (dec_1sec) is untouched.
PLOT3_DRAW_EVERY_NTH <- 10   # <-- adjust: higher = faster/lighter plot to render

dec_1sec_for_plot <- dec_1sec[seq(1, nrow(dec_1sec), by = PLOT3_DRAW_EVERY_NTH), ]

p3 <- ggplot(dec_1sec_for_plot, aes(x = datetime, y = value)) +
  geom_line(linewidth = 0.15, color = "darkblue") +
  labs(x = "Date", y = expression(paste("RSAM (", mu, "m ", s^-1, ")"))) +
  theme_classic(base_family = "serif", base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(p3)
ggsave(file.path(data_dir_seconds, "Plot3_Dec2019_1secRSAM.png"),
       p3, width = 10, height = 6, dpi = 150)

# =============================================================================
# PLOT 4: December 8-9, 2019 only, 100Hz averaged to 1-second RSAM
# =============================================================================

dec8_9_doy <- c(342, 343)   # 8 Dec, 9 Dec 2019

dec8_9_files <- tibble::tibble(path = seconds_files) %>%
  dplyr::mutate(
    fname = basename(path),
    doy   = as.integer(stringr::str_extract(fname, "(?<=^2019\\.)\\d{3}"))
  ) %>%
  dplyr::filter(doy %in% dec8_9_doy) %>%
  dplyr::arrange(doy)

cat("\nPlot 4 -- Dec 8-9 files found:", nrow(dec8_9_files), "of 2 expected\n")

cat("Averaging Dec 8-9, 2019 from 100Hz to 1-second RSAM...\n")
dec8_9_1sec <- purrr::map_dfr(dec8_9_files$path, average_day_to_1sec)

# 2 days at 1-second resolution is ~172,800 points -- light enough to draw
# with little or no thinning.
PLOT4_DRAW_EVERY_NTH <- 1   # <-- adjust if still too heavy to render

dec8_9_1sec_for_plot <- dec8_9_1sec[seq(1, nrow(dec8_9_1sec), by = PLOT4_DRAW_EVERY_NTH), ]

p4 <- ggplot(dec8_9_1sec_for_plot, aes(x = datetime, y = value)) +
  geom_line(linewidth = 0.2, color = "darkblue") +
  scale_x_datetime(date_breaks = "3 hours", date_labels = "%d %H:%M") +
  labs(x = "Date/Hour", y = expression(paste("RSAM (", mu, "m ", s^-1, ")"))) +
  theme_classic(base_family = "serif", base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(p4)
ggsave(file.path(data_dir_seconds, "Plot4_Dec8_9_2019_1secRSAM.png"),
       p4, width = 10, height = 6, dpi = 150)