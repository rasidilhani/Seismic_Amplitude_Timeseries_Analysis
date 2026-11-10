library(readr)
library(here)
library(dplyr)
library(ggplot2)
library(lubridate)     # needed for year(), yday(), ymd()

#CACHE_DIR <- here("daily_cache")    # relative to the project folder
CACHE_DIR <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/daily_cache"

## Read one day for one station from the cache
read_day <- function(date, station = "WIZ") {
  yr  <- year(date)
  doy <- yday(date)
  
  # e.g. daily_cache/2016/2016.118/WIZ.NZ/2016.118.WIZ.timeseries.csv
  f <- file.path(CACHE_DIR,
                 yr,
                 sprintf("%d.%03d", yr, doy),
                 paste0(station, ".NZ"),
                 sprintf("%d.%03d.%s.timeseries.csv", yr, doy, station))
  
  if (!file.exists(f)) return(NULL)      # day not in the cache
  
  d <- read_csv(f, show_col_types = FALSE)
  names(d)[1] <- "datetime"              # first column has no name in the files
  d
}

#d <- read_day(ymd("2016-04-27"))
#d <- read_day(ymd("2016-09-15"))
d <- read_day(ymd("2019-12-09"))
head(d)
nrow(d)

#month <- "2016-04"   
#month <- "2016-09"                         
month <- "2019-12"         

from <- ymd(paste0(month, "-01"))
to   <- from + months(1) - days(1)
all_dates <- seq(from, to, by = "day")

day_list <- list()
for (i in 1:length(all_dates)) {
  d <- read_day(all_dates[i])
  if (is.null(d)) {
    cat("missing:", format(all_dates[i]), "\n")
    next
  }
  day_list[[i]] <- d
}

x <- bind_rows(day_list)                    
x <- x %>% mutate(time_nz = with_tz(datetime, "UTC"))

nrow(x)                                    

#onset <- ymd_hms("2016-04-27 09:37:00", tz = "UTC")
#onset <- ymd_hms("2016-09-15 00:14:52", tz = "UTC")
onset <- ymd_hms("2019-12-09 01:11:00", tz = "UTC")

p <- ggplot(x, aes(x = datetime, y = RSAM)) +
  geom_line(linewidth = 0.3) +
  geom_vline(xintercept = onset, colour = "red", linewidth = 0.8) +
  scale_y_log10() +
  labs(x = "UTC time", y = "RSAM (log scale)",
       title = paste("WIZ RSAM -", month)) +
  theme_bw()

p                                         

ggsave(here("results", paste0("WIZ_RSAM_", month, ".pdf")), p, width = 10, height = 4)

###################################################################################
###############################################################################
# =====================================================================
# WIZ RSAM around one eruption (all times in UTC)
#   1. read whole months around the eruption
#   2. time series plot with eruption windows and weeks
#   3. ordinal-pattern histograms: 4 weeks before + eruption week
#   4. daily Shannon entropy and Fisher information (D = 3)
# =====================================================================

library(readr)
library(here)
library(dplyr)
library(ggplot2)
library(lubridate)
library(ordinalpatterns)
library(StatOrdPattHxC)

CACHE_DIR <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/daily_cache"

## Read one day for one station from the cache
read_day <- function(date, station = "WSRZ") {
  yr  <- year(date)
  doy <- yday(date)
  f <- file.path(CACHE_DIR,
                 yr,
                 sprintf("%d.%03d", yr, doy),
                 paste0(station, ".NZ"),
                 sprintf("%d.%03d.%s.timeseries.csv", yr, doy, station))
  if (!file.exists(f)) return(NULL)
  d <- read_csv(f, show_col_types = FALSE)
  names(d)[1] <- "datetime"
  d
}

# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------
D   <- 3              # embedding dimension
tau <- 1              # delay
n_weeks_before <- 4   # weeks before the eruption week

# eruptions from the VAL table (UTC)
eruptions <- data.frame(
  case  = c("2016-04", "2016-09", "2019-12", "2019-12"),
  start = ymd_hms(c("2016-04-27 23:30:00", "2016-09-13 03:15:00",
                    "2019-12-09 01:30:00", "2019-12-09 03:25:00"), tz = "UTC"),
  end   = ymd_hms(c("2016-04-28 06:44:59", "2016-09-15 00:14:59",
                    "2019-12-09 03:24:59", "2019-12-11 21:19:59"), tz = "UTC"),
  VAL   = c(3, 3, 4, 3)
)

#case <- "2016-04"
#case <- "2016-09"
case <- "2019-12"                                # <- choose the eruption

er <- eruptions[eruptions$case == case, ]        # one or two rows
eruption_day <- as.Date(min(er$start))           # day of the first onset

# ---------------------------------------------------------------------
# Period: whole months, from the month that contains "4 weeks before"
# to the end of the month that contains the eruption week
# ---------------------------------------------------------------------
from <- floor_date(eruption_day - days(7 * n_weeks_before), "month")
to   <- floor_date(eruption_day + days(6), "month") + months(1) - days(1)

# to choose the period yourself, use for example:
# from <- ymd("2019-11-01")
# to   <- ymd("2019-12-31")

cat("Reading from", format(from), "to", format(to), "\n")

# ---------------------------------------------------------------------
# Read the days
# ---------------------------------------------------------------------
all_dates <- seq(from, to, by = "day")

day_list <- list()
for (i in 1:length(all_dates)) {
  d <- read_day(all_dates[i])
  if (is.null(d)) {
    cat("missing:", format(all_dates[i]), "\n")
    next
  }
  day_list[[i]] <- d
}
x <- bind_rows(day_list)
x$day <- as.Date(x$datetime)                     # UTC day
nrow(x)

# ---------------------------------------------------------------------
# Weeks: eruption week = 7 days starting on the eruption day
#        week -1 = the 7 days before that, and so on
# ---------------------------------------------------------------------
week_start <- eruption_day - days(7 * (n_weeks_before:0))
week_lab   <- c(paste0("week -", n_weeks_before:1), "eruption week")

x$week <- NA
for (k in 1:length(week_start)) {
  in_week <- x$day >= week_start[k] & x$day < week_start[k] + 7
  x$week[in_week] <- week_lab[k]
}
x$week <- factor(x$week, levels = week_lab)

print(data.frame(week = week_lab, first_day = week_start, last_day = week_start + 6))

# ---------------------------------------------------------------------
# 1. Time series plot
#    red band = eruption window(s), dashed lines = week starts
# ---------------------------------------------------------------------
p_ts <- ggplot(x, aes(x = datetime, y = RSAM)) +
  annotate("rect", xmin = er$start, xmax = er$end, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.3) +
  geom_line(linewidth = 0.3) +
  geom_vline(xintercept = as.POSIXct(week_start[-length(week_start)]),
             linetype = 2, colour = "grey40") +                        # weeks -4 ... -1
  geom_vline(xintercept = as.POSIXct(week_start[length(week_start)]),
             linetype = 2, colour = "blue", linewidth = 0.8) +         # eruption week start
  scale_y_log10() +
  labs(x = "UTC time", y = "RSAM (log scale)",
       title = paste("WIZ RSAM -", case, "eruption")) +
  theme_bw()
p_ts

# ---------------------------------------------------------------------
# Ordinal-pattern measures with the packages
# returns NULL if the series is too short
# ---------------------------------------------------------------------
op_measures <- function(v) {
  series <- v[!is.na(v)]                          # op_pd needs a series without NA
  n_emb  <- length(series) - (D - 1) * tau
  if (n_emb < factorial(D) * 5) return(NULL)
  
  xpd   <- ordinalpatterns::op_pd(series, D = D, tau = tau)
  probs <- xpd$probabilities
  
  Hs <- ordinalpatterns::permutation_entropy(probs, normalized = TRUE)
  Cs <- ordinalpatterns::statistical_complexity(probs, entropy = Hs, normalized = TRUE)
  Fi <- ordinalpatterns::fisher_information(probs)
  JS <- ordinalpatterns::jsd(probs) / log(2)
  Cf <- JS * Fi
  
  list(probs = probs, n_emb = n_emb, H = Hs, C = Cs, F = Fi, C_F = Cf)
}

pattern_id <- paste0("P", 1:factorial(D))         # P1 ... P6 for D = 3

# ---------------------------------------------------------------------
# 2. Daily ordinal patterns, Shannon entropy and Fisher information
# ---------------------------------------------------------------------
days_all   <- sort(unique(x$day))
daily      <- NULL
daily_long <- NULL
for (i in 1:length(days_all)) {
  v <- x$RSAM[x$day == days_all[i]]
  m <- op_measures(v)
  if (is.null(m)) next                           # too few values that day
  daily      <- rbind(daily, data.frame(day = days_all[i], n_emb = m$n_emb,
                                        H = m$H, C = m$C, F = m$F, C_F = m$C_F))
  daily_long <- rbind(daily_long, data.frame(day = days_all[i], pattern = pattern_id,
                                             prob = m$probs))
}

er_days_start <- as.Date(er$start)
er_days_end   <- as.Date(er$end) + 1

p_H <- ggplot(daily, aes(x = day, y = H)) +
  annotate("rect", xmin = er_days_start, xmax = er_days_end, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.3) +
  geom_line() + geom_point(size = 1) +
  geom_vline(xintercept = week_start[-length(week_start)],
             linetype = 2, colour = "grey40") +
  geom_vline(xintercept = week_start[length(week_start)],
             linetype = 2, colour = "blue", linewidth = 0.8) +
  labs(x = "UTC day", y = "normalized Shannon entropy",
       title = paste("Daily Shannon entropy (D =", D, ") -", case)) +
  theme_bw()
p_H

p_C <- ggplot(daily, aes(x = day, y = C)) +
  annotate("rect", xmin = er_days_start, xmax = er_days_end, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.3) +
  geom_line() + geom_point(size = 1) +
  geom_vline(xintercept = week_start[-length(week_start)],
             linetype = 2, colour = "grey40") +
  geom_vline(xintercept = week_start[length(week_start)],
             linetype = 2, colour = "blue", linewidth = 0.8) +
  labs(x = "UTC day", y = "Complexity",
       title = paste("Daily Complexity (D =", D, ") -", case)) +
  theme_bw()
p_C


p_F <- ggplot(daily, aes(x = day, y = F)) +
  annotate("rect", xmin = er_days_start, xmax = er_days_end, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.3) +
  geom_line() + geom_point(size = 1) +
  geom_vline(xintercept = week_start[-length(week_start)],
             linetype = 2, colour = "grey40") +
  geom_vline(xintercept = week_start[length(week_start)],
             linetype = 2, colour = "blue", linewidth = 0.8) +
  labs(x = "UTC day", y = "Fisher information",
       title = paste("Daily Fisher information (D =", D, ") -", case)) +
  theme_bw()
p_F

# ---------------------------------------------------------------------
# 3. One OP histogram for the whole eruption month
#    all RSAM values of the month = one time series
# ---------------------------------------------------------------------
month_start <- floor_date(eruption_day, "month")
month_end   <- month_start + months(1) - days(1)

v_month <- x$RSAM[x$day >= month_start & x$day <= month_end]
m_month <- op_measures(v_month)

hist_month <- data.frame(pattern = factor(pattern_id, levels = pattern_id),
                         prob    = m_month$probs)
print(hist_month)

p_hist <- ggplot(hist_month, aes(x = pattern, y = prob)) +
  geom_col(fill = "steelblue") +
  labs(x = "ordinal pattern", y = "probability",
       title    = paste("Ordinal patterns (D =", D, ") -", format(month_start, "%B %Y")),
       subtitle = paste0("H = ", round(m_month$H, 3),
                         "   F = ", round(m_month$F, 3),
                         "   patterns = ", m_month$n_emb)) +
  theme_bw()
p_hist

# ---------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------
ggsave(here("results", paste0("WIZ_RSAM_ts_", case, ".pdf")),        p_ts,   width = 10, height = 4)
ggsave(here("results", paste0("WIZ_daily_H_", case, ".pdf")),        p_H,    width = 10, height = 4)
ggsave(here("results", paste0("WIZ_daily_F_", case, ".pdf")),        p_F,    width = 10, height = 4)
write_csv(daily,   here("results", paste0("WIZ_daily_OP_D3_",  case, ".csv")))