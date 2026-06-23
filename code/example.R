setwd("C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis")

library(tidyverse)
library(lubridate)

csv_files <- list.files("data", pattern = "^WIZ_NZ_[0-9]{4}\\.csv$", full.names = TRUE)

seismic_all <- list()
for (f in csv_files) {
  yr <- as.integer(str_extract(f, "[0-9]{4}"))
  d <- read_csv(f, show_col_types = FALSE)
  d$year <- yr
  utc_time <- as.POSIXct(d$unix_timestamp, origin = "1970-01-01", tz = "UTC")
  d$datetime_nz <- with_tz(utc_time, tzone = "Pacific/Auckland")
  seismic_all[[f]] <- d
}
seismic_all <- bind_rows(seismic_all)

seismic_all$month <- month(seismic_all$datetime_nz, label = TRUE)
seismic_all$quarter <- paste0("Q", quarter(seismic_all$datetime_nz))

jan_2011 <- seismic_all %>% filter(year == 2011, month == "Jan")
nrow(jan_2011)
range(jan_2011$datetime_nz)

library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(purrr)
library(stringr)

files <- list.files("data", pattern = "WIZ_NZ_.*\\.csv$", full.names = TRUE)

seismic_all <- map_dfr(files, function(file) {
  
  year_value <- str_extract(basename(file), "[0-9]{4}")
  
  read_csv(file, show_col_types = FALSE) |>
    mutate(
      year = as.factor(year_value),
      datetime_nz = as.POSIXct(
        unix_timestamp,
        origin = "1970-01-01",
        tz = "Pacific/Auckland"
      ),
      doy = yday(datetime_nz)
    )
})

# 12 distinct colours for 2011--2022
year_colours <- c(
  "2011" = "#1b9e77",
  "2012" = "#d95f02",
  "2013" = "#7570b3",
  "2014" = "#e7298a",
  "2015" = "#66a61e",
  "2016" = "#e6ab02",
  "2017" = "#a6761d",
  "2018" = "#666666",
  "2019" = "#1f78b4",
  "2020" = "#b2df8a",
  "2021" = "#fb9a99",
  "2022" = "#cab2d6"
)

# Serif theme for report figures
report_theme <- theme_bw(base_family = "serif", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

# Annual Displacement for all years
ggplot(seismic_all,
       aes(x = doy,
           y = displacement_avg_m,
           colour = year,
           group = year)) +
  geom_line(alpha = 0.7, linewidth = 0.6) +
  scale_colour_manual(values = year_colours) +
  labs(
    title = "Annual Displacement Time Series (2011–2022)",
    x = "Day of Year",
    y = "Average Displacement (m)",
    colour = "Year"
  ) +
  report_theme

# Annual RSAM for all years
ggplot(seismic_all,
       aes(x = doy,
           y = rsam_avg,
           colour = year,
           group = year)) +
  geom_line(alpha = 0.7, linewidth = 0.6) +
  scale_colour_manual(values = year_colours) +
  labs(
    title = "Annual RSAM Time Series (2011–2022)",
    x = "Day of Year",
    y = "RSAM",
    colour = "Year"
  ) +
  report_theme

# Yearly RSAM using facets
ggplot(seismic_all,
       aes(x = doy,
           y = rsam_avg,
           colour = year,
           group = year)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~year, ncol = 3) +
  scale_colour_manual(values = year_colours) +
  labs(
    title = "RSAM Time Series by Year",
    x = "Day of Year",
    y = "RSAM",
    colour = "Year"
  ) +
  report_theme +
  theme(legend.position = "none")

# Monthly averages
seismic_monthly <- seismic_all |>
  mutate(month = month(datetime_nz, label = TRUE)) |>
  group_by(year, month) |>
  summarise(
    displacement = mean(displacement_avg_m, na.rm = TRUE),
    rsam = mean(rsam_avg, na.rm = TRUE),
    .groups = "drop"
  )

# Monthly mean displacement by year
ggplot(seismic_monthly,
       aes(x = month,
           y = displacement,
           colour = year,
           group = year)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(values = year_colours) +
  labs(
    title = "Monthly Mean Displacement by Year",
    x = "Month",
    y = "Mean Displacement (m)",
    colour = "Year"
  ) +
  report_theme

# Monthly mean rsam by year
ggplot(seismic_monthly,
       aes(x = month,
           y = rsam,
           colour = year,
           group = year)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(values = year_colours) +
  labs(
    title = "Monthly Mean RSAM by Year",
    x = "Month",
    y = "Mean RSAM (m)",
    colour = "Year"
  ) +
  report_theme
