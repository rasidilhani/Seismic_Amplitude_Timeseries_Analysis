library(readr)
library(dplyr)
library(ggplot2)

setwd("C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis")

WIZ1 <- read_csv("daily_cache/2016/2016.118/WIZ.NZ/2016.118.WIZ.timeseries.csv")
view(WIZ1)
#WIZ2 <- read_csv("daily_cache/2016/2016.119/WIZ.NZ/2016.119.WIZ.timeseries.csv")
#WIZ3 <- read_csv("daily_cache/2016/2016.120/WIZ.NZ/2016.120.WIZ.timeseries.csv")

#WIZ <- bind_rows(WIZ1, WIZ2, WIZ3) 
WIZ <- bind_rows(WIZ1, WIZ2) 
names(WIZ1)[1] <- "datetime"         

ggplot(WIZ1, aes(x = datetime, y = RSAM)) +
  geom_line() +
  labs(x = "Time (UTC)", y = "RSAM") +
  theme_bw()


## Simple zoom plot from the local daily_cache (no server / rsync needed)

library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)

CACHE_DIR <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/daily_cache"

## Read one day for one station from the cache
read_day <- function(date, station = "WIZ") {
  yr  <- year(date)
  doy <- yday(date)
  f <- file.path(CACHE_DIR, sprintf("%d/%d.%03d/%s.NZ/%d.%03d.%s.timeseries.csv",
                                    yr, yr, doy, station, yr, doy, station))
  if (!file.exists(f)) return(NULL)
  d <- read_csv(f, show_col_types = FALSE)
  names(d)[1] <- "datetime"          # first column is unnamed in the files
  d
}

## ---- 1. Choose the days to load -------------------------------------------
dates <- seq(as.Date("2016-04-27"), as.Date("2016-04-29"), by = "day")

WIZ <- bind_rows(lapply(dates, read_day, station = "WIZ")) %>%
  mutate(datetime = as_datetime(datetime, tz = "UTC")) %>%
  arrange(datetime)

## ---- 2. Eruption window (UTC) ---------------------------------------------
## 28 Apr 2016, 11:30-18:45 NZST  =  27 Apr 23:30 - 28 Apr 06:45 UTC
#eq_start <- as.POSIXct("2016-04-27 23:30:00", tz = "UTC")
#eq_end   <- as.POSIXct("2016-04-28 06:45:00", tz = "UTC")

## ---- 3. Filter to the zoom period (here: 2 days either side) --------------
#zoom <- WIZ %>%
#  dplyr::filter(datetime >= eq_start - days(2),
#                datetime <= eq_end   + days(2))

## ---- 4. Plot ----------------------------------------------------------------
ggplot(zoom, aes(x = datetime, y = RSAM)) +
  annotate("rect", xmin = eq_start, xmax = eq_end, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.2) +
  geom_line() +
  scale_x_datetime(date_labels = "%d %b\n%H:%M") +
  labs(title = "WIZ: 28 Apr 2016 eruption", x = "Time (UTC)", y = "RSAM") +
  theme_bw()
