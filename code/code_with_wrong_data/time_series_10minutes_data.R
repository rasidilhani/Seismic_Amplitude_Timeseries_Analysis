setwd("C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis")

library(tidyverse)
library(lubridate)
library(Polychrome)

# ---- Read data ----
f <- "data/WIZ_NZ_2019.csv"
d <- read_csv(f, show_col_types = FALSE)
d$year <- 2019

# ---- Original timestamp (UTC), used directly for plotting ----
d$utc_time <- as.POSIXct(d$unix_timestamp, origin = "1970-01-01", tz = "UTC")

# ---- Define the three windows ----
# 1. Full year 2019
d_year <- d %>% filter(year(utc_time) == 2019)

# 2. December 2019
d_month <- d %>% filter(year(utc_time) == 2019, month(utc_time) == 12)

# 3. Custom window: 5 Dec to 13 Dec 2019
week_start <- as.Date("2019-12-05")
week_end   <- as.Date("2019-12-13")
d_week <- d %>% filter(as.Date(utc_time) >= week_start,
                       as.Date(utc_time) <= week_end)
cat("Week window:", format(week_start), "to", format(week_end), "\n")

# ---- Three distinct colors ----
plot_colors <- c(year = "#1B9E77", month = "#D95F02", week = "#7570B3")

# ---- Shared plain theme: no gridlines, serif font throughout ----
theme_plain_serif <- theme_bw(base_size = 13, base_family = "serif") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(family = "serif"),
    axis.text  = element_text(family = "serif"),
    plot.title = element_text(family = "serif")
  )

# ---- Plot 1: Full year 2019 ----
p_year <- ggplot(d_year, aes(x = utc_time, y = rsam_avg)) +
  geom_line(color = plot_colors["year"], linewidth = 0.3) +
  labs(title = "WIZ RSAM — Full Year 2019",
       x = "Date (UTC)", y = "RSAM (avg)") +
  theme_plain_serif

# ---- Plot 2: December 2019 ----
p_month <- ggplot(d_month, aes(x = utc_time, y = rsam_avg)) +
  geom_line(color = plot_colors["month"], linewidth = 0.4) +
  labs(title = "WIZ RSAM — December 2019",
       x = "Date (UTC)", y = "RSAM (avg)") +
  theme_plain_serif

# ---- Plot 3: Week containing 9 Dec 2019 ----
p_week <- ggplot(d_week, aes(x = utc_time, y = rsam_avg)) +
  geom_line(color = plot_colors["week"], linewidth = 0.6) +
  labs(title = paste0("WIZ RSAM — Week of ", format(week_start, "%d %b"),
                      " to ", format(week_end, "%d %b %Y"),
                      " (includes 9 Dec)"),
       x = "Date/Time (UTC)", y = "RSAM (avg)") +
  theme_plain_serif

# ---- Display in R session ----
print(p_year)
print(p_month)
print(p_week)

# ---- Save plots ----
dir.create("results", showWarnings = FALSE)

ggsave("results/WIZ_2019_full_year_rsam.pdf", p_year, width = 10, height = 5)
ggsave("results/WIZ_2019_december_rsam.pdf", p_month, width = 10, height = 5)
ggsave("results/WIZ_2019_week_of_dec9_rsam.pdf", p_week, width = 10, height = 5)

############################
# ---- Custom window: just 9 Dec 2019 ----
day_start <- as.POSIXct("2019-12-08 00:00:00", tz = "UTC")
day_end   <- as.POSIXct("2019-12-10 23:59:59", tz = "UTC")

d_day <- d %>% filter(utc_time >= day_start, utc_time <= day_end)

cat("Day window:", format(day_start), "to", format(day_end), "\n")

# ---- Plot: 8-10 Dec 2019, hourly ticks, labeled with date + time ----
p_day_hourly <- ggplot(d_day, aes(x = utc_time, y = rsam_avg)) +
  geom_line(color = plot_colors["week"], linewidth = 0.5) +
  scale_x_datetime(date_breaks = "3 hours", date_labels = "%d %b\n%H:%M") +
  labs(title = "WIZ RSAM — 8–10 December 2019 (hourly)",
       x = "Time (UTC)", y = "RSAM (avg)") +
  theme_plain_serif +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))

print(p_day_hourly)
ggsave("results/WIZ_2019_dec8to10_hourly_rsam.pdf", p_day_hourly, width = 12, height = 5)