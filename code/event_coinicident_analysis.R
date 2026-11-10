# =====================================================================
# Event coincidence analysis (ECA) - weekly
# One row = one week. All times in UTC.
# =====================================================================

library(readr)
library(lubridate)

# ---------------------------------------------------------------------
# Settings - change only these lines
# ---------------------------------------------------------------------
data_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"

station   <- "WIZ"                 # "WIZ" or "WSRZ"
f         <- "Average"           # feature to test
threshold <- 0.80                  # top 10% of weeks = event
dT        <- 2                     # weeks before the eruption week we look
tau       <- 0                     # weeks we skip right before the eruption week

# eruption onsets (UTC, from the VAL table)
onset <- ymd_hms(c("2016-04-27 23:30:00",
                   "2016-09-13 03:15:00",
                   "2019-12-09 01:30:00",
                   "2019-12-09 03:25:00"), tz = "UTC")

# ---------------------------------------------------------------------
# Step 1: read the weekly results
# ---------------------------------------------------------------------
file_name <- paste0("HC_RSAM_", station, "_20140701_to_20191231_clean_Copy.csv")
feat <- read_csv(file.path(data_dir, file_name), show_col_types = FALSE)

n_weeks  <- nrow(feat)
week     <- 1:n_weeks              # week number = row number
week_end <- feat$Week_End

cat(station, "| weeks:", n_weeks, "\n")

# ---------------------------------------------------------------------
# Step 2: eruption weeks
# eruption week = first week that ends after the onset
# ---------------------------------------------------------------------
erupt_week <- rep(NA, length(onset))
for (j in 1:length(onset)) {
  erupt_week[j] <- which(week_end > onset[j])[1]
}
print(data.frame(onset = onset, erupt_week = erupt_week))

# ---------------------------------------------------------------------
# Step 3: feature events = weeks where the feature is unusually high
# ---------------------------------------------------------------------
x   <- feat[[f]]
thr <- quantile(x, threshold, na.rm = TRUE)

is_ev    <- !is.na(x) & x > thr                 # TRUE/FALSE for every week
first_ev <- is_ev & !c(FALSE, head(is_ev, -1))  # keep first week of each run
event_week <- week[first_ev]

cat(f, "| threshold =", round(thr, 4),
    "| weeks above =", sum(is_ev),
    "| events =", length(event_week), "\n")

plot(week, x, type = "l", xlab = "week number", ylab = f, main = paste(station, "-", f))
abline(h = thr, col = "red", lty = 2)
points(event_week, x[event_week], col = "red", pch = 16)
abline(v = erupt_week, col = "blue", lwd = 2)

# ---------------------------------------------------------------------
# Step 4: count coincidences
# an event "hits" an eruption if it is 1+tau ... dT+tau weeks earlier
# ---------------------------------------------------------------------
N_B <- length(erupt_week)   # eruptions
N_A <- length(event_week)   # feature events

# precursor: was each eruption preceded by a feature event?
hit_B        <- rep(FALSE, N_B)
weeks_before <- rep(NA, N_B)
for (j in 1:N_B) {
  gap <- erupt_week[j] - event_week
  hit_B[j] <- any(gap >= 1 + tau & gap <= dT + tau)
  if (any(gap >= 1)) weeks_before[j] <- min(gap[gap >= 1])
}

# trigger: was each feature event followed by an eruption?
hit_A <- rep(FALSE, N_A)
for (i in 1:N_A) {
  gap <- erupt_week - event_week[i]
  hit_A[i] <- any(gap >= 1 + tau & gap <= dT + tau)
}

cat("\ndT =", dT, "weeks | tau =", tau, "weeks\n")
cat("Precursor rate:", sum(hit_B), "of", N_B, "eruptions\n")
cat("Trigger rate:  ", sum(hit_A), "of", N_A, "events\n")
print(data.frame(eruption = format(onset, "%Y-%m-%d"),
                 hit = hit_B,
                 closest_event_weeks_before = weeks_before))

# ---------------------------------------------------------------------
# Step 5a: significance - Poisson formula
# ---------------------------------------------------------------------
p_one_B <- 1 - (1 - dT / n_weeks)^N_A
p_one_A <- 1 - (1 - dT / n_weeks)^N_B

p_prec_pois <- 1 - pbinom(sum(hit_B) - 1, N_B, p_one_B)
p_trig_pois <- 1 - pbinom(sum(hit_A) - 1, N_A, p_one_A)

cat("\nPrecursor p-value (Poisson):", round(p_prec_pois, 4), "\n")
cat("Trigger   p-value (Poisson):", round(p_trig_pois, 4), "\n")

# ---------------------------------------------------------------------
# Step 5b: significance - surrogates
# move all events forward by 1, 2, ..., n_weeks-1 weeks (wrap around)
# ---------------------------------------------------------------------
n_shift <- n_weeks - 1
null_B  <- rep(NA, n_shift)
null_A  <- rep(NA, n_shift)

for (s in 1:n_shift) {
  ev_s <- (event_week + s - 1) %% n_weeks + 1    # shifted event weeks
  
  kB <- 0
  for (j in 1:N_B) {
    gap <- erupt_week[j] - ev_s
    if (any(gap >= 1 + tau & gap <= dT + tau)) kB <- kB + 1
  }
  kA <- 0
  for (i in 1:N_A) {
    gap <- erupt_week - ev_s[i]
    if (any(gap >= 1 + tau & gap <= dT + tau)) kA <- kA + 1
  }
  null_B[s] <- kB
  null_A[s] <- kA
}

p_prec_surr <- (1 + sum(null_B >= sum(hit_B))) / (n_shift + 1)
p_trig_surr <- (1 + sum(null_A >= sum(hit_A))) / (n_shift + 1)

cat("\nPrecursor p-value (surrogate):", round(p_prec_surr, 4), "\n")
cat("Trigger   p-value (surrogate):", round(p_trig_surr, 4), "\n")
cat("How many eruptions the shifted events hit:\n")
print(table(null_B))


########################################################################
#########################################################################
# =====================================================================
# Event coincidence analysis (ECA) - weekly - grid search
# Tries every combination of threshold, dT and tau, and gives a
# p-value corrected for picking the best combination.
# All times in UTC.
# =====================================================================

library(readr)
library(lubridate)

# ---------------------------------------------------------------------
# Settings - change only these lines
# ---------------------------------------------------------------------
data_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"

station   <- "WIZ"                     # "WIZ" or "WSRZ"
#f         <- "Disequilibrium"           # feature to test
f         <- "Var_C_Shannon"
direction <- "upper"                   # "upper" = unusually high, "lower" = unusually low

thr_grid <- c(0.80, 0.85, 0.90, 0.95)  # thresholds (percentiles)
dT_grid  <- c(1, 2, 3)                 # weeks before the eruption week we look
tau_grid <- c(0, 1)                    # weeks we skip right before the eruption week

# eruption onsets (UTC, from the VAL table)
onset <- ymd_hms(c("2016-04-27 23:30:00",
                   "2016-09-13 03:15:00",
                   "2019-12-09 01:30:00"), tz = "UTC")

# ---------------------------------------------------------------------
# Step 1: read the weekly results
# ---------------------------------------------------------------------
file_name <- paste0("HC_RSAM_", station, "_20140701_to_20191231_clean.csv")
feat <- read_csv(file.path(data_dir, file_name), show_col_types = FALSE)

week_end <- feat$Week_End

# week number from the date (works even if some weeks are missing)
# weeks start every 6 days
week    <- round(as.numeric(difftime(feat$Week_Start, min(feat$Week_Start),
                                     units = "days")) / 6) + 1
n_weeks <- max(week)

cat(station, "| rows:", nrow(feat), "| weeks in period:", n_weeks, "\n")

# ---------------------------------------------------------------------
# Step 2: eruption weeks (first week that ends after the onset)
# ---------------------------------------------------------------------
erupt_week <- rep(NA, length(onset))
for (j in 1:length(onset)) {
  erupt_week[j] <- week[which(week_end > onset[j])[1]]
}
print(data.frame(onset = onset, erupt_week = erupt_week))

N_B     <- length(erupt_week)
n_shift <- n_weeks - 1
x       <- feat[[f]]

# ---------------------------------------------------------------------
# Step 3: all combinations of threshold, dT and tau
# ---------------------------------------------------------------------
settings <- expand.grid(thr = thr_grid, dT = dT_grid, tau = tau_grid)
n_set    <- nrow(settings)
cat("Number of settings:", n_set, "\n")

settings$n_events  <- NA
settings$hits_B    <- NA     # eruptions hit
settings$hits_A    <- NA     # events followed by an eruption
null_B <- matrix(NA, n_shift, n_set)   # hits for every shift and setting

# ---------------------------------------------------------------------
# Step 4: count coincidences for every setting (real + shifted events)
# ---------------------------------------------------------------------
for (k in 1:n_set) {
  
  # feature events for this threshold
  if (direction == "upper") {
    ev <- !is.na(x) & x > quantile(x, settings$thr[k], na.rm = TRUE)
  } else {
    ev <- !is.na(x) & x < quantile(x, 1 - settings$thr[k], na.rm = TRUE)
  }
  ev_week <- week[ev & !c(FALSE, head(ev, -1))]    # first week of each run
  
  lo <- 1 + settings$tau[k]                        # coincidence window
  hi <- settings$dT[k] + settings$tau[k]
  
  # real events: eruptions hit
  kB <- 0
  for (j in 1:N_B) {
    gap <- erupt_week[j] - ev_week
    if (any(gap >= lo & gap <= hi)) kB <- kB + 1
  }
  
  # real events: events followed by an eruption
  kA <- 0
  for (i in seq_along(ev_week)) {
    gap <- erupt_week - ev_week[i]
    if (any(gap >= lo & gap <= hi)) kA <- kA + 1
  }
  
  settings$n_events[k] <- length(ev_week)
  settings$hits_B[k]   <- kB
  settings$hits_A[k]   <- kA
  
  # shifted events
  for (s in 1:n_shift) {
    ev_s <- (ev_week + s - 1) %% n_weeks + 1
    kB_s <- 0
    for (j in 1:N_B) {
      gap <- erupt_week[j] - ev_s
      if (any(gap >= lo & gap <= hi)) kB_s <- kB_s + 1
    }
    null_B[s, k] <- kB_s
  }
}

# ---------------------------------------------------------------------
# Step 5: p-value for every setting
# ---------------------------------------------------------------------
settings$p <- NA
for (k in 1:n_set) {
  settings$p[k] <- (1 + sum(null_B[, k] >= settings$hits_B[k])) / (n_shift + 1)
}

settings$precursor_rate <- round(settings$hits_B / N_B, 2)
settings$trigger_rate   <- round(settings$hits_A / settings$n_events, 3)
settings$p              <- round(settings$p, 4)

print(settings[order(settings$p), ])

# ---------------------------------------------------------------------
# Step 6: p-value corrected for picking the best setting
# For every shift: what is the best p-value it would get over all
# settings? Compare our best p-value with those.
# ---------------------------------------------------------------------
p_null <- matrix(NA, n_shift, n_set)
for (k in 1:n_set) {
  for (s in 1:n_shift) {
    p_null[s, k] <- (sum(null_B[, k] >= null_B[s, k]) +          # shifts, incl. itself
                       (settings$hits_B[k] >= null_B[s, k])) /      # the real result
      (n_shift + 1)
  }
}
best_null <- apply(p_null, 1, min)
best_obs  <- min(settings$p)
p_grid    <- (1 + sum(best_null <= best_obs)) / (n_shift + 1)

cat("\nBest setting:\n")
print(settings[which.min(settings$p), ])
cat("Best p =", best_obs, "| corrected for", n_set, "settings:", round(p_grid, 4), "\n")

# ---------------------------------------------------------------------
# Step 7: plot the feature with the threshold of the best setting
# ---------------------------------------------------------------------
best     <- settings[which.min(settings$p), ]
thr_best <- quantile(x, ifelse(direction == "upper", best$thr, 1 - best$thr), na.rm = TRUE)
if (direction == "upper") ev <- !is.na(x) & x > thr_best else ev <- !is.na(x) & x < thr_best
ev_week_best <- week[ev & !c(FALSE, head(ev, -1))]

plot(week, x, type = "l", xlab = "week number", ylab = f,
     main = paste0(station, " - ", f, "  (threshold ", best$thr,
                   ", dT = ", best$dT, ", tau = ", best$tau, ")"))
abline(h = thr_best, col = "red", lty = 2)
points(ev_week_best, x[match(ev_week_best, week)], col = "red", pch = 16)
abline(v = erupt_week, col = "blue", lwd = 2)

# ---------------------------------------------------------------------
# Save the table
# ---------------------------------------------------------------------
write_csv(settings, file.path(data_dir, paste0("ECA_grid_", station, "_", f, ".csv")))

#############################################################
###############################################################
# Use the CoinCalc package
# =====================================================================
# Event coincidence analysis (ECA) - weekly - grid search
#   own code : coincidences, circular-shift test, grid-corrected p
#   CoinCalc : cross-check of the coincidence rates + Poisson and
#              waiting-time surrogate p-values for every setting
# All times in UTC.
# =====================================================================

library(readr)
library(lubridate)
library(CoinCalc)

# ---------------------------------------------------------------------
# Settings - change only these lines
# ---------------------------------------------------------------------
data_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"

station   <- "WIZ"                     # "WIZ" or "WSRZ"
f         <- "Fisher_Info"           # feature to test
direction <- "upper"                   # "upper" = unusually high, "lower" = unusually low

thr_grid <- c(0.80, 0.85, 0.90, 0.95)  # thresholds (percentiles)
dT_grid  <- c(1, 2, 3)                 # weeks before the eruption week we look
tau_grid <- c(0, 1)                    # weeks we skip right before the eruption week

n_reps_cc <- 1000                      # surrogates for the CoinCalc waiting-time test

# eruption onsets (UTC, from the VAL table)
onset <- ymd_hms(c("2016-04-27 23:30:00",
                   "2016-09-13 03:15:00",
                   "2019-12-09 01:30:00"), tz = "UTC")

# ---------------------------------------------------------------------
# Step 1: read the weekly results
# ---------------------------------------------------------------------
file_name <- paste0("HC_RSAM_", station, "_20140701_to_20191231_clean.csv")
feat <- read_csv(file.path(data_dir, file_name), show_col_types = FALSE)

week_end <- feat$Week_End

# week number from the date (weeks start every 6 days)
week    <- round(as.numeric(difftime(feat$Week_Start, min(feat$Week_Start),
                                     units = "days")) / 6) + 1
n_weeks <- max(week)

cat(station, "| rows:", nrow(feat), "| weeks in period:", n_weeks, "\n")

# ---------------------------------------------------------------------
# Step 2: eruption weeks (first week that ends after the onset)
# ---------------------------------------------------------------------
erupt_week <- rep(NA, length(onset))
for (j in 1:length(onset)) {
  erupt_week[j] <- week[which(week_end > onset[j])[1]]
}
print(data.frame(onset = onset, erupt_week = erupt_week))

N_B     <- length(erupt_week)
n_shift <- n_weeks - 1
x       <- feat[[f]]

# eruption series as 0/1 on the weekly grid (for CoinCalc)
B_bin <- rep(0, n_weeks)
B_bin[erupt_week] <- 1

# ---------------------------------------------------------------------
# Step 3: all combinations of threshold, dT and tau
# ---------------------------------------------------------------------
settings <- expand.grid(thr = thr_grid, dT = dT_grid, tau = tau_grid)
n_set    <- nrow(settings)
cat("Number of settings:", n_set, "\n")

settings$n_events <- NA
settings$hits_B   <- NA     # eruptions hit
settings$hits_A   <- NA     # events followed by an eruption
null_B  <- matrix(NA, n_shift, n_set)   # hits for every shift and setting
ev_list <- list()                       # event weeks for every setting

# ---------------------------------------------------------------------
# Step 4: count coincidences for every setting (real + shifted events)
# ---------------------------------------------------------------------
for (k in 1:n_set) {
  
  # feature events for this threshold
  if (direction == "upper") {
    ev <- !is.na(x) & x > quantile(x, settings$thr[k], na.rm = TRUE)
  } else {
    ev <- !is.na(x) & x < quantile(x, 1 - settings$thr[k], na.rm = TRUE)
  }
  ev_week <- week[ev & !c(FALSE, head(ev, -1))]    # first week of each run
  ev_list[[k]] <- ev_week
  
  lo <- 1 + settings$tau[k]                        # coincidence window
  hi <- settings$dT[k] + settings$tau[k]
  
  # real events: eruptions hit
  kB <- 0
  for (j in 1:N_B) {
    gap <- erupt_week[j] - ev_week
    if (any(gap >= lo & gap <= hi)) kB <- kB + 1
  }
  
  # real events: events followed by an eruption
  kA <- 0
  for (i in seq_along(ev_week)) {
    gap <- erupt_week - ev_week[i]
    if (any(gap >= lo & gap <= hi)) kA <- kA + 1
  }
  
  settings$n_events[k] <- length(ev_week)
  settings$hits_B[k]   <- kB
  settings$hits_A[k]   <- kA
  
  # shifted events
  for (s in 1:n_shift) {
    ev_s <- (ev_week + s - 1) %% n_weeks + 1
    kB_s <- 0
    for (j in 1:N_B) {
      gap <- erupt_week[j] - ev_s
      if (any(gap >= lo & gap <= hi)) kB_s <- kB_s + 1
    }
    null_B[s, k] <- kB_s
  }
}

# ---------------------------------------------------------------------
# Step 5: p-value for every setting (circular shift)
# ---------------------------------------------------------------------
settings$p <- NA
for (k in 1:n_set) {
  settings$p[k] <- (1 + sum(null_B[, k] >= settings$hits_B[k])) / (n_shift + 1)
}
settings$precursor_rate <- round(settings$hits_B / N_B, 2)
settings$trigger_rate   <- round(settings$hits_A / settings$n_events, 3)
settings$p              <- round(settings$p, 4)

# ---------------------------------------------------------------------
# Step 6: CoinCalc for every setting
# Our window: event 1+tau ... dT+tau weeks before the eruption week.
# CoinCalc window: tau_cc ... tau_cc + delT_cc  ->  tau_cc = 1 + tau,
#                                                   delT_cc = dT - 1
# ---------------------------------------------------------------------

# first a test call: look at the names of what CC.eca.ts returns
A_bin <- rep(0, n_weeks)
A_bin[ev_list[[1]]] <- 1
test <- CC.eca.ts(A_bin, B_bin, delT = settings$dT[1] - 1, tau = settings$tau[1] + 1,
                  sym = FALSE, sigtest = "poisson")
str(test)     # <- check the names here and change the four lines below if needed

settings$cc_prec_rate <- NA
settings$cc_trig_rate <- NA
settings$cc_p_poisson <- NA
settings$cc_p_wt      <- NA

for (k in 1:n_set) {
  A_bin <- rep(0, n_weeks)
  A_bin[ev_list[[k]]] <- 1
  
  r_pois <- CC.eca.ts(A_bin, B_bin, delT = settings$dT[k] - 1, tau = settings$tau[k] + 1,
                      sym = FALSE, sigtest = "poisson")
  r_wt   <- CC.eca.ts(A_bin, B_bin, delT = settings$dT[k] - 1, tau = settings$tau[k] + 1,
                      sym = FALSE, sigtest = "wt.surrogate", reps = n_reps_cc)
  
  settings$cc_prec_rate[k] <- r_pois[["precursor coincidence rate"]]
  settings$cc_trig_rate[k] <- r_pois[["trigger coincidence rate"]]
  settings$cc_p_poisson[k] <- r_pois[["p-value precursor"]]
  settings$cc_p_wt[k]      <- r_wt[["p-value precursor"]]
}

# check: our rates and CoinCalc rates should be the same
settings$same_prec <- abs(settings$precursor_rate - round(settings$cc_prec_rate, 2)) < 0.01
settings$same_trig <- abs(settings$trigger_rate   - round(settings$cc_trig_rate, 3)) < 0.001
cat("\nSettings where our precursor rate = CoinCalc:", sum(settings$same_prec), "of", n_set, "\n")
cat("Settings where our trigger rate   = CoinCalc:", sum(settings$same_trig), "of", n_set, "\n")

print(settings[order(settings$p), ])

# ---------------------------------------------------------------------
# Step 7: p-value corrected for picking the best setting
# ---------------------------------------------------------------------
p_null <- matrix(NA, n_shift, n_set)
for (k in 1:n_set) {
  for (s in 1:n_shift) {
    p_null[s, k] <- (sum(null_B[, k] >= null_B[s, k]) +          # shifts, incl. itself
                       (settings$hits_B[k] >= null_B[s, k])) /      # the real result
      (n_shift + 1)
  }
}
best_null <- apply(p_null, 1, min)
best_obs  <- min(settings$p)
p_grid    <- (1 + sum(best_null <= best_obs)) / (n_shift + 1)

cat("\nBest setting:\n")
print(settings[which.min(settings$p), ])
cat("Best p =", best_obs, "| corrected for", n_set, "settings:", round(p_grid, 4), "\n")

# ---------------------------------------------------------------------
# Step 8: plot the feature with the threshold of the best setting
# ---------------------------------------------------------------------
best     <- settings[which.min(settings$p), ]
thr_best <- quantile(x, ifelse(direction == "upper", best$thr, 1 - best$thr), na.rm = TRUE)
if (direction == "upper") {
  ev <- !is.na(x) & x > thr_best
} else {
  ev <- !is.na(x) & x < thr_best
}
ev_week_best <- week[ev & !c(FALSE, head(ev, -1))]

par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(week, x, type = "l", xlab = "week number", ylab = f,
     main = paste0(station, " - ", f, "  (threshold ", best$thr,
                   ", dT = ", best$dT, ", tau = ", best$tau, ")"))
abline(h = thr_best, col = "red", lty = 2)
points(ev_week_best, x[match(ev_week_best, week)], col = "red", pch = 16)
abline(v = erupt_week, col = "blue", lwd = 2)

# ---------------------------------------------------------------------
# Save the table
# ---------------------------------------------------------------------
write_csv(settings, file.path(data_dir, paste0("ECA_grid_CoinCalc_", station, "_", f, ".csv")))