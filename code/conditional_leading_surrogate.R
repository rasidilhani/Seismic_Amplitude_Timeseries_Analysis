# =====================================================================
# Weekly OP features vs eruptions  (all times UTC)
#
#   Part 1: conditional probability and lift
#   Part 2: lead-lag  (a) lead correlation  (b) superposed epoch
#   Part 3: surrogate test (shift the eruption labels in time)
#
# One row = one weekly window (7 days long, a new window every 6 days).
# =====================================================================

library(readr)
library(lubridate)

# ---------------------------------------------------------------------
# Settings - change only these lines
# ---------------------------------------------------------------------
data_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"
station  <- "WIZ"                                  # "WIZ" or "WSRZ"

features <- c("H_Shannon", "C_Shannon", "Fisher_Info", "C_Fisher",
              "H_Renyi", "C_Renyi", "H_Tsallis", "C_Tsallis",
              "Var_H_Shannon", "Var_C_Shannon", "Var_H_Renyi", "Var_H_Tsallis",
              "Disequilibrium")

n_bins    <- 5     # quintiles for the lift table
horizon   <- 2     # weeks: "an eruption starts within the next 2 weeks"
max_lead  <- 10    # weeks: longest lead time checked
quiet_gap <- 15    # weeks: more than 15 weeks from an eruption = quiet week
cor_cut   <- 0.9   # |r| above this = highly correlated features

# eruptions (UTC, VAL table). The two Dec 2019 windows are one episode.
eruptions <- data.frame(
  episode = c(1, 2, 3, 3),
  start   = ymd_hms(c("2016-04-27 23:30:00", "2016-09-13 03:15:00",
                      "2019-12-09 01:30:00", "2019-12-09 03:25:00"), tz = "UTC"),
  end     = ymd_hms(c("2016-04-28 06:44:59", "2016-09-15 00:14:59",
                      "2019-12-09 03:24:59", "2019-12-11 21:19:59"), tz = "UTC")
)
episode_names <- c("2016-04-27", "2016-09-13", "2019-12-09")
n_episodes    <- 3

out_dir <- file.path(data_dir, paste0("conditional_analysis_", station))
dir.create(out_dir, showWarnings = FALSE)

# ---------------------------------------------------------------------
# Step 1: read the weekly features and put them on a full weekly grid
# (a missing week becomes an empty row, so week numbers = real time)
# ---------------------------------------------------------------------
feat <- read_csv(file.path(data_dir, paste0("HC_RSAM_", station, "_20140701_to_20191231.csv")),
                 show_col_types = FALSE)

first_start <- min(feat$Week_Start)
row_week    <- round(as.numeric(difftime(feat$Week_Start, first_start, units = "days")) / 6) + 1

n          <- max(row_week)
week       <- 1:n
week_start <- first_start + days(6 * (week - 1))
week_end   <- week_start + days(7)

X <- matrix(NA, n, length(features))
colnames(X) <- features
X[row_week, ] <- as.matrix(feat[, features])

cat(station, "| rows in file:", nrow(feat), "| weeks in period:", n,
    "| missing weeks:", n - nrow(feat), "\n")

# ---------------------------------------------------------------------
# Step 2: labels for every week
#   during = 1 if the week overlaps an eruption
#   ahead  = 1 if an eruption starts 1 ... horizon weeks later
#   quiet  = TRUE if the week is far away from every eruption
# ---------------------------------------------------------------------
during      <- rep(0, n)
episode_now <- rep(NA, n)
for (e in 1:nrow(eruptions)) {
  hit <- week_start <= eruptions$end[e] & week_end >= eruptions$start[e]
  during[hit]      <- 1
  episode_now[hit] <- eruptions$episode[e]
}

# eruption week = first week that ends after the onset
onset_week <- rep(NA, n_episodes)
for (k in 1:n_episodes) {
  first_onset   <- min(eruptions$start[eruptions$episode == k])
  onset_week[k] <- which(week_end > first_onset)[1]
}
print(data.frame(eruption = episode_names, onset_week, week_start = week_start[onset_week]))

ahead        <- rep(0, n)
episode_next <- rep(NA, n)
dist_erupt   <- rep(NA, n)
for (w in 1:n) {
  gap <- onset_week - w                          # weeks from week w to each eruption
  if (during[w] == 0 & any(gap >= 1 & gap <= horizon)) {
    ahead[w]        <- 1
    episode_next[w] <- which(gap >= 1 & gap <= horizon)[1]
  }
  dist_erupt[w] <- min(abs(gap))
}
quiet <- during == 0 & dist_erupt > quiet_gap

cat("Weeks during an eruption:", sum(during),
    "| weeks before an eruption:", sum(ahead),
    "| quiet weeks:", sum(quiet), "\n")

# ---------------------------------------------------------------------
# Step 3: features that are not highly correlated
# keep a feature only if |r| < cor_cut with every feature already kept
# ---------------------------------------------------------------------
cm   <- cor(X, use = "pairwise.complete.obs")
keep <- features[1]
for (f in features[-1]) {
  if (all(abs(cm[f, keep]) < cor_cut)) keep <- c(keep, f)
}
cat("\nWeakly correlated features:\n"); print(keep)

# =====================================================================
# PART 1: P(eruption | feature bin) and lift
#   bins  = quintiles of the feature
#   lift  = P(eruption | bin) / P(eruption)
# =====================================================================
lift_table <- function(x, y, episode) {
  brk  <- unique(quantile(x, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  bin  <- cut(x, breaks = brk, include.lowest = TRUE, labels = FALSE)
  base <- mean(y[!is.na(bin)])
  out  <- NULL
  for (b in sort(unique(na.omit(bin)))) {
    i <- which(bin == b)
    out <- rbind(out, data.frame(
      bin        = b,
      from       = brk[b],
      to         = brk[b + 1],
      n_weeks    = length(i),
      n_event    = sum(y[i]),
      n_episodes = length(unique(na.omit(episode[i][y[i] == 1]))),
      p_event    = mean(y[i]),
      baseline   = base,
      lift       = mean(y[i]) / base))
  }
  out
}

lift_all <- NULL
for (f in features) {
  a <- lift_table(X[, f], during, episode_now)
  a$outcome <- "during"
  b <- lift_table(X[, f], ahead, episode_next)
  b$outcome <- "ahead"
  ab <- rbind(a, b)
  ab$feature <- f
  lift_all <- rbind(lift_all, ab)
}
write_csv(lift_all, file.path(out_dir, "part1_lift_table.csv"))

cat("\nHighest lifts with at least 2 different eruptions:\n")
top <- lift_all[lift_all$n_episodes >= 2, ]
print(head(top[order(-top$lift), ], 15))

# =====================================================================
# PART 2a: lead correlation
#   r(k) = cor( feature in week w , during in week w + k )
#   k = 0 same week (no warning), k = 1 one week earlier, ...
# =====================================================================
lead_r <- matrix(NA, max_lead + 1, length(features))
colnames(lead_r) <- features

for (f in features) {
  x <- X[, f]
  for (k in 0:max_lead) {
    lead_r[k + 1, f] <- cor(x[1:(n - k)], during[(1 + k):n], use = "complete.obs")
  }
}
write_csv(data.frame(lag_weeks = 0:max_lead, lead_r),
          file.path(out_dir, "part2_lead_correlation.csv"))

par(mfrow = c(4, 4), mar = c(3.5, 3.5, 2, 1), mgp = c(2, 0.6, 0))
for (f in features) {
  plot(0:max_lead, lead_r[, f], type = "b", pch = 16,
       xlab = "feature leads eruption by (weeks)", ylab = "r", main = f)
  abline(h = 0, col = "grey50")
}
dev.copy(pdf, file.path(out_dir, "part2_lead_correlation.pdf"), width = 11, height = 8.5)
dev.off()

# =====================================================================
# PART 2b: superposed epoch analysis
#   z = (feature - mean of quiet weeks) / sd of quiet weeks
#   0 = eruption week, -1 = the week before, ...
#   lead_weeks = weeks in a row with |z| > 2 just before the eruption
# =====================================================================
cols      <- c("red", "blue", "darkgreen")
epoch_res <- NULL

par(mfrow = c(4, 4), mar = c(3.5, 3.5, 2, 1), mgp = c(2, 0.6, 0))
for (f in features) {
  x <- X[, f]
  z <- (x - mean(x[quiet], na.rm = TRUE)) / sd(x[quiet], na.rm = TRUE)
  
  plot(NA, xlim = c(-max_lead, 1), ylim = c(-8, 8),
       xlab = "weeks relative to eruption week", ylab = "z (vs quiet)", main = f)
  abline(h = c(-2, 2), lty = 2, col = "grey40")
  abline(v = 0)
  
  for (k in 1:n_episodes) {
    w <- (onset_week[k] - max_lead):(onset_week[k] + 1)
    w <- w[w >= 1 & w <= n]
    lines(w - onset_week[k], z[w], type = "b", pch = 16, col = cols[k], lwd = 1.5)
    
    lead <- 0
    for (j in 1:max_lead) {
      zj <- z[onset_week[k] - j]
      if (!is.na(zj) & abs(zj) > 2) lead <- j else break
    }
    before <- (onset_week[k] - max_lead):(onset_week[k] - 1)
    epoch_res <- rbind(epoch_res, data.frame(
      feature       = f,
      eruption      = episode_names[k],
      z_week_before = round(z[onset_week[k] - 1], 2),
      max_abs_z     = round(max(abs(z[before]), na.rm = TRUE), 2),
      lead_weeks    = lead))
  }
}
plot.new()
legend("center", legend = episode_names, col = cols, lwd = 2, bty = "n", title = "Eruption")
dev.copy(pdf, file.path(out_dir, "part2_superposed_epoch.pdf"), width = 11, height = 8.5)
dev.off()

write_csv(epoch_res, file.path(out_dir, "part2_epoch_lead_weeks.csv"))
cat("\nFeatures with |z| > 2 right before an eruption:\n")
print(epoch_res[epoch_res$lead_weeks > 0, ])

# =====================================================================
# PART 3: surrogate test
#   shift the eruption labels by quiet_gap ... n - quiet_gap weeks
#   (wrap around at the end) and recompute. Every shift used once.
# =====================================================================
shifts  <- quiet_gap:(n - quiet_gap)
n_shift <- length(shifts)

max_lift <- function(bin, y) {
  ok <- !is.na(bin)
  max(tapply(y[ok], bin[ok], mean)) / mean(y[ok])
}

surr_res <- NULL
for (f in features) {
  x   <- X[, f]
  brk <- unique(quantile(x, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  bin <- cut(x, breaks = brk, include.lowest = TRUE, labels = FALSE)
  
  obs_during <- max_lift(bin, during)
  obs_ahead  <- max_lift(bin, ahead)
  obs_lead   <- max(abs(lead_r[2:(max_lead + 1), f]), na.rm = TRUE)   # lags 1 ... max_lead
  
  null_during <- rep(NA, n_shift)
  null_ahead  <- rep(NA, n_shift)
  null_lead   <- rep(NA, n_shift)
  
  for (s in 1:n_shift) {
    sh <- shifts[s]
    y1 <- c(during[(sh + 1):n], during[1:sh])
    y2 <- c(ahead[(sh + 1):n],  ahead[1:sh])
    
    null_during[s] <- max_lift(bin, y1)
    null_ahead[s]  <- max_lift(bin, y2)
    
    r <- rep(NA, max_lead)
    for (k in 1:max_lead) {
      r[k] <- cor(x[1:(n - k)], y1[(1 + k):n], use = "complete.obs")
    }
    null_lead[s] <- max(abs(r), na.rm = TRUE)
  }
  
  surr_res <- rbind(surr_res, data.frame(
    feature       = f,
    lift_during   = round(obs_during, 2),
    p_lift_during = (1 + sum(null_during >= obs_during)) / (n_shift + 1),
    lift_ahead    = round(obs_ahead, 2),
    p_lift_ahead  = (1 + sum(null_ahead >= obs_ahead)) / (n_shift + 1),
    max_lead_r    = round(obs_lead, 3),
    p_lead_r      = (1 + sum(null_lead >= obs_lead)) / (n_shift + 1)))
  cat("done:", f, "\n")
}

# correction for testing 13 features (Benjamini-Hochberg)
surr_res$q_lift_during <- p.adjust(surr_res$p_lift_during, method = "BH")
surr_res$q_lift_ahead  <- p.adjust(surr_res$p_lift_ahead,  method = "BH")
surr_res$q_lead_r      <- p.adjust(surr_res$p_lead_r,      method = "BH")

write_csv(surr_res, file.path(out_dir, "part3_surrogate_pvalues.csv"))
print(surr_res[order(surr_res$p_lead_r), ])

cat("\nAll results saved in:", out_dir, "\n")