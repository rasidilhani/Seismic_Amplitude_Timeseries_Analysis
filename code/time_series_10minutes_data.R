setwd("C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis")

library(tidyverse)
library(lubridate)

# =============================================================================
# Settings
# =============================================================================
DATA_DIR    <- "data"
RESULTS_DIR <- "results"
dir.create(RESULTS_DIR, showWarnings = FALSE)

STATIONS <- c("WIZ", "WSRZ")
YEARS    <- 2014:2019   # one file per station per year is expected

SOURCE_DATA_TZ <- "NZ"   # <-- change to "UTC" if you confirm that instead

TIME_AXIS_LABEL <- if (SOURCE_DATA_TZ == "NZ") "NZ Local Time" else "UTC"

ANALYSIS_START <- as.POSIXct("2014-07-01 00:00:00", tz = "NZ")
ANALYSIS_END   <- as.POSIXct("2019-12-31 23:59:59", tz = "NZ")

station_colors <- c(WIZ = "#1B9E77", WSRZ = "#7570B3")

# =============================================================================
# 1. Read and combine all yearly files for one station
# =============================================================================
read_station_years <- function(station, years, data_dir) {
  purrr::map_dfr(years, function(yr) {
    f <- file.path(data_dir, sprintf("%s_NZ_%d_Annual_RSAM.csv", station, yr))
    if (!file.exists(f)) {
      warning(sprintf("File not found, skipping: %s", f))
      return(NULL)
    }
   
    d <- read_csv(f, show_col_types = FALSE, col_types = cols(.default = "c"))
    
      d <- tibble(
      Year = as.integer(d$Year),
      DOY  = as.integer(d$DOY),
      RSAM = as.numeric(d$RSAM)
    ) %>%
      group_by(Year, DOY) %>%
      mutate(
        sample_in_day = row_number(),
        utc_time = as.POSIXct(sprintf("%d-01-01 00:00:00", Year), tz = "NZ") +
          lubridate::days(DOY - 1) + lubridate::minutes(10 * sample_in_day)
      ) %>%
      ungroup() %>%
      select(utc_time, RSAM)
    
    d$year    <- yr
    d$station <- station
    d
  })
}

prepare_station <- function(d) {
  d %>%
    mutate(utc_time = as.POSIXct(utc_time, tz = "NZ")) %>%
    filter(utc_time >= ANALYSIS_START, utc_time <= ANALYSIS_END) %>%
    arrange(utc_time)
}

wiz  <- read_station_years("WIZ",  YEARS, DATA_DIR) %>% prepare_station()
wsrz <- read_station_years("WSRZ", YEARS, DATA_DIR) %>% prepare_station()

cat("WIZ range:  ", format(range(wiz$utc_time)),  "\n")
cat("WSRZ range: ", format(range(wsrz$utc_time)), "\n")

# =============================================================================
# 2. Eruption windows to highlight and zoom into
# =============================================================================

eruption_windows_raw <- tribble(
  ~label,                                       ~nz_start,              ~nz_end,                ~offset_hours,
  "28 Apr 2016 eruption",                       "2016-04-28 11:30:00",  "2016-04-28 18:44:59",   12,
  "13-15 Sep 2016 eruption",                    "2016-09-13 15:15:00",  "2016-09-15 12:14:59",   12,
  "9 Dec 2019 explosive eruption",              "2019-12-09 14:30:00",  "2019-12-09 16:24:59",   13,
  "9-12 Dec 2019 continued eruptive activity",  "2019-12-09 16:25:00",  "2019-12-12 10:19:59",   13
)

eruption_windows <- eruption_windows_raw %>%
  mutate(
    nz_start = as.POSIXct(nz_start, tz = "UTC"),   # tagged UTC only so R can
    nz_end   = as.POSIXct(nz_end,   tz = "UTC"),   # do arithmetic on it; the
    # values themselves are the
    # NZ local clock times above
    start = if (SOURCE_DATA_TZ == "NZ") nz_start else nz_start - hours(offset_hours),
    end   = if (SOURCE_DATA_TZ == "NZ") nz_end   else nz_end   - hours(offset_hours)
  ) %>%
  select(label, start, end)

# =============================================================================
# 3. Shared theme
# =============================================================================
theme_plain_serif <- theme_bw(base_size = 13, base_family = "serif") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(family = "serif"),
    axis.text  = element_text(family = "serif"),
    plot.title = element_text(family = "serif"),
    legend.position = "bottom"
  )

# =============================================================================
# 4. Plots A & B: full-period time series, one station each, no highlight
# =============================================================================
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

# =============================================================================
# 5. Plots C & D: full-period time series with eruption windows highlighted
# =============================================================================
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

# ---- Display the four full-period plots ----
print(p_wiz_full)
print(p_wsrz_full)
print(p_wiz_full_hl)
print(p_wsrz_full_hl)

# ---- Save the four full-period plots ----
ggsave(file.path(RESULTS_DIR, "WIZ_2014to2019_full_rsam.pdf"),               p_wiz_full,     width = 12, height = 5)
ggsave(file.path(RESULTS_DIR, "WSRZ_2014to2019_full_rsam.pdf"),              p_wsrz_full,    width = 12, height = 5)
ggsave(file.path(RESULTS_DIR, "WIZ_2014to2019_full_rsam_highlighted.pdf"),   p_wiz_full_hl,  width = 12, height = 5)
ggsave(file.path(RESULTS_DIR, "WSRZ_2014to2019_full_rsam_highlighted.pdf"),  p_wsrz_full_hl, width = 12, height = 5)

# =============================================================================
# 6. Zoom plots: one per eruption window, both stations overlaid
# =============================================================================
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