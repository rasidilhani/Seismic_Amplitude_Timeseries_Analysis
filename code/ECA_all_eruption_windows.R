# =====================================================================
# Event coincidence analysis (ECA) - SLIDING 7-day windows - ALL features
#
# Input: one row per 7-day window; windows start one day apart
#        (6-day overlap), e.g. 2007-01-01 -> 2007-01-08, 2007-01-02 -> 2007-01-09, ...
#        Each row is time-stamped by its Week_End (the window has seen data
#        up to that moment only), so a row can be a precursor only if its
#        window ends before the eruption onset.
#
# Method:
#   - event days     : every row beyond the threshold (all of them count)
#   - episodes       : runs of consecutive rows beyond the threshold
#   - precursor hit  : an eruption with at least one event row whose window
#                      ended 7*tau ... 7*(tau+dT) days before the onset
#   - trigger rate   : share of episodes with an eruption following within
#                      that lead time
#   - significance   : circular shift of the whole event series by every
#                      possible number of rows (keeps the overlap/clustering),
#                      grid-corrected p, Benjamini-Hochberg across features
#   - Step 7         : composites of correlated features
# The step between rows is detected from the file, so the same script also
# works for 7-day windows every 6 days.
# All times in UTC.
#
# Eruptions (see the catalogue in the settings):
#   - the eruptions counted in the test are the major impulsive eruptions of
#     Dempsey et al. (2020); the October 2013 episode is one eruptive period
#   - an eruptive period has a start (used as the onset for the precursor
#     windows) and an end
#   - events that are not certainly eruptions (Sep-Dec 2012 lava eruption with
#     unknown onset, 2016-09 ash emission) are "neutral": the rows around them
#     are left out of the test, so they count neither as hits nor as false alarms
#   - windows that contain eruptive activity (syn-eruptive rows) are also left
#     out, so an eruption is not detected by its own signal
#   - 'scenario' switches between the main definition and two sensitivity runs
# =====================================================================

library(readr)
library(lubridate)

# ---------------------------------------------------------------------
# Settings - change only these lines
# ---------------------------------------------------------------------
data_dir    <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"
results_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"
file_name   <- "HC_RSAM_WIZ_6daysoverlap_20070101_to_20191231_clean.csv"    # <- your sliding-window results file

station   <- "WIZ"                     # "WIZ" or "WSRZ"
variable  <- "RSAM"                    # rows of the 'Variable' column to use
direction <- "upper"                   # default: "upper" = unusually high, "lower" = unusually low
lower_features <- c("H_Shannon", "H_Renyi", "H_Tsallis")   # low entropy = unusual

all_features <- c(
  "H_Shannon", "C_Shannon", "Fisher_Info", "C_Fisher",
  "H_Renyi", "C_Renyi", "H_Tsallis", "C_Tsallis",
  "Var_H_Shannon", "Var_C_Shannon", "Var_H_Renyi", "Var_H_Tsallis",
  "Disequilibrium"
)

thr_grid <- c(0.80, 0.85, 0.90, 0.95)  # thresholds (percentiles)
dT_grid  <- c(1, 2, 3)                 # length of the look-back window (weeks = 7 days)
tau_grid <- c(0, 1)                    # extra weeks skipped right before the onset

tick_by  <- "1 year"                   # x-axis tick spacing, e.g. "1 year", "6 months"
date_fmt <- "%Y-%m-%d"                 # format of the date labels

# one colour per threshold (colour-blind safe)
thr_cols <- rep(c("#E69F00", "#009E73", "#0072B2", "#CC79A7", "#D55E00"),
                length.out = length(thr_grid))

# ---------------------------------------------------------------------
# Eruption catalogue (all times UTC)
#   class "definite"     : major impulsive eruptions used by Dempsey et al. (2020)
#   class "probable"     : an eruption, but slow and over weeks, onset unknown
#   class "uncertain"    : not clear whether it was an eruption
#   class "not_eruption" : in the GeoNet VAL table, but not treated as an eruption
#   start      : onset used for the precursor windows
#   end        : end of the eruptive period (NA = a single eruption)
#   sub_events : further eruptions inside a period (drawn as thin lines only)
# ---------------------------------------------------------------------
eruptions <- data.frame(
  id    = c("2012-08", "2012-lava", "2013-geysers", "2013-08",
            "2013-10", "2016-04", "2016-09", "2019-12"),
  class = c("definite", "probable", "not_eruption", "definite",
            "definite", "definite", "uncertain", "definite"),
  start = c("2012-08-04 16:52:00", "2012-09-01 00:00:00", "2013-01-01 00:00:00",
            "2013-08-19 22:23:00", "2013-10-04 03:35:00", "2016-04-27 09:37:00",
            "2016-09-13 00:00:00", "2019-12-09 01:11:00"),
  end   = c(NA, "2012-12-31 23:59:00", "2013-04-10 23:59:00",
            NA, "2013-10-11 07:09:00", NA,
            NA, NA),
  sub_events = c("", "", "", "", "2013-10-08 02:05:00;2013-10-11 07:09:00", "", "", ""),
  note  = c("eruption",
            "slow lava eruption, some time Sep-Dec 2012",
            "many geysers, Jan - early Apr 2013",
            "steam and mud eruption",
            "eruptive period 3-11 Oct 2013 (several eruptions)",
            "eruption",
            "small ash emission, roughly 2016-09-13",
            "eruption"),
  stringsAsFactors = FALSE)

# Which classes count as eruptions in the test:
#   "eruption"   = counted (precursor windows, hits, p-values)
#   "neutral"    = left out of the test (rows around it are masked)
#   "background" = treated like any other time (only drawn in the figures)
scenario <- "definite"     # "definite" (main), "with_lava" or "with_all" (sensitivity)
scenario_roles <- list(
  definite  = c(definite = "eruption", probable = "neutral",  uncertain = "neutral",  not_eruption = "background"),
  with_lava = c(definite = "eruption", probable = "eruption", uncertain = "neutral",  not_eruption = "background"),
  with_all  = c(definite = "eruption", probable = "eruption", uncertain = "eruption", not_eruption = "background"))
mask_syn_eruptive <- TRUE  # leave out windows that contain eruptive activity
win_len_days      <- 7     # length of every feature window (days)

# raw RSAM for the top panel of the stacked figures
#   rsam_source = a FILE  : one cache file (.rds, .RData/.rda or .csv)
#               = a FOLDER: every .csv below it (sub-folders included) whose name contains
#                 the station, e.g. WIZ_NZ_2016_Annual_RSAM.csv or the server tree
#                 2016/2016.118/WIZ.NZ/2016.118.WIZ.timeseries.csv
# The time column (unix seconds, a date-time, or date-time text) and the RSAM column
# (rsam_column, else the first column whose name starts with "rsam") are found automatically.
rsam_source  <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/daily_cache"
rsam_column  <- "RSAM"                        # column with the RSAM values
rsam_tz      <- "UTC"                         # zone of text timestamps that carry none
# (columns with "nz" in the name are read as NZ time)
rsam_log     <- TRUE                          # log10 y-axis (RSAM spans orders of magnitude)
rsam_refresh <- FALSE                         # TRUE = re-read the source instead of the daily cache

# direction of every feature
feat_dir <- setNames(ifelse(all_features %in% lower_features, "lower", direction), all_features)
dir_note <- if (length(lower_features) > 0) {
  paste0(paste(lower_features, collapse = ", "), ": lower; others: ", direction)
} else direction

if (!scenario %in% names(scenario_roles)) stop("Unknown scenario: ", scenario)
out_dir <- file.path(results_dir, paste0("ECA_sliding_", station, "_", variable, "_", scenario))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------
# Step 1: read the windows and place them on a regular grid
# ---------------------------------------------------------------------
feat <- read_csv(file.path(data_dir, file_name), show_col_types = FALSE)
if ("Variable" %in% names(feat)) feat <- feat[feat$Variable == variable, ]

missing <- setdiff(all_features, names(feat))
if (length(missing) > 0) stop("Columns not in the file: ", paste(missing, collapse = ", "))

t_end <- as.POSIXct(feat$Week_End, tz = "UTC")
feat  <- feat[order(t_end), ]
t_end <- sort(t_end)

step_days <- median(as.numeric(diff(t_end), units = "days"))   # 1 for daily-sliding windows
t0        <- t_end[1]
gap_days  <- as.numeric(diff(t_end), units = "days")
if (any(gap_days <= 0)) {
  stop("Several rows have the same Week_End - check the 'Variable' filter and the file.")
}
# grid position of every row: one step per step_days, gaps become empty rows
idx <- c(1, 1 + cumsum(pmax(1, round(gap_days / step_days))))
n      <- max(idx)                                     # rows of the regular grid
t_grid <- t0 + (0:(n - 1)) * step_days * 86400         # window end of every grid row

X <- matrix(NA_real_, n, length(all_features), dimnames = list(NULL, all_features))
for (f in all_features) X[idx, f] <- feat[[f]]
n_shift <- n - 1

cat(station, variable, "| windows:", nrow(feat), "| step:", step_days, "day(s)",
    "| grid rows:", n, "| missing rows:", n - nrow(feat), "\n",
    "first window end:", format(t_grid[1]), "| last:", format(t_grid[n]), "\n")

# ---------------------------------------------------------------------
# Step 2: eruptions and their precursor windows
#   lead = onset - window end (days); a row is in the precursor window of
#   eruption j for setting (dT, tau) if  7*tau <= lead < 7*(tau + dT)
# ---------------------------------------------------------------------
# catalogue -> times and roles
eruptions$start <- ymd_hms(eruptions$start, tz = "UTC")
end_txt         <- eruptions$end
eruptions$end   <- eruptions$start
has_end         <- !is.na(end_txt) & end_txt != ""
eruptions$end[has_end] <- ymd_hms(end_txt[has_end], tz = "UTC")
if (any(eruptions$end < eruptions$start)) stop("An eruption ends before it starts - check the catalogue.")
eruptions$role  <- unname(scenario_roles[[scenario]][eruptions$class])
eruptions$in_record <- eruptions$start > t_grid[1] & eruptions$start <= t_grid[n]
eruptions <- eruptions[order(eruptions$start), ]

pos_of <- function(t) as.numeric(difftime(t, t0, units = "days")) / step_days + 1   # x position

# counted eruptions
ER <- eruptions[eruptions$role == "eruption", ]
if (any(!ER$in_record)) cat("Eruptions outside the record, dropped:", ER$id[!ER$in_record], "\n")
ER    <- ER[ER$in_record, ]
N_B   <- nrow(ER)
if (N_B == 0) stop("No eruption inside the record.")
onset     <- ER$start
onset_end <- ER$end
erupt_lab <- paste0("E", seq_len(N_B))               # E1, E2, ... in time order
eruptions$label <- ""
eruptions$label[match(ER$id, eruptions$id)] <- erupt_lab

erupt_pos     <- pos_of(onset)
erupt_end_pos <- pos_of(onset_end)
sub_t   <- ymd_hms(unlist(strsplit(ER$sub_events[ER$sub_events != ""], ";")), tz = "UTC", quiet = TRUE)
sub_pos <- pos_of(sub_t[!is.na(sub_t)])

# events that are not counted, drawn for context: N = neutral, B = background
OTH <- eruptions[eruptions$role != "eruption" & eruptions$end > t_grid[1] & eruptions$start <= t_grid[n], ]
if (nrow(OTH) > 0) {
  OTH$label <- ave(OTH$role, OTH$role, FUN = function(r) paste0(ifelse(r == "neutral", "N", "B"), seq_along(r)))
  eruptions$label[match(OTH$id, eruptions$id)] <- OTH$label
}

lead <- matrix(sapply(onset, function(o) as.numeric(difftime(o, t_grid, units = "days"))),
               nrow = n)                                                     # n x N_B

# rows left out of the test (mask): neither hits nor false alarms
#   syn-eruptive: the 7-day window (t - 7 d, t] touches the eruption or its period
#   neutral     : from the longest look-back before the start to the end of the window after it
max_look  <- 7 * (max(tau_grid) + max(dT_grid))      # days
touches   <- function(s, e) t_grid > s & t_grid < e + win_len_days * 86400
precursor_any <- rowSums(lead >= 0 & lead < max_look) > 0
mask_syn  <- rep(FALSE, n)
if (mask_syn_eruptive) for (j in seq_len(N_B)) mask_syn <- mask_syn | touches(onset[j], onset_end[j])
NEU       <- eruptions[eruptions$role == "neutral", ]
mask_neu  <- rep(FALSE, n)
for (j in seq_len(nrow(NEU))) {
  mask_neu <- mask_neu | (t_grid >= NEU$start[j] - max_look * 86400 &
                            t_grid <  NEU$end[j] + win_len_days * 86400)
}
mask_neu <- mask_neu & !precursor_any                # never hide a counted eruption's window
mask     <- mask_syn | mask_neu
cat("Scenario:", scenario, "|", N_B, "eruptions counted |",
    sum(mask_syn), "syn-eruptive rows and", sum(mask_neu & !mask_syn), "neutral rows left out\n")

grid  <- expand.grid(thr = thr_grid, dT = dT_grid, tau = tau_grid)
n_set <- nrow(grid)

win <- vector("list", n_set)
for (k in 1:n_set) {
  lo <- 7 * grid$tau[k]
  hi <- 7 * (grid$tau[k] + grid$dT[k])
  W  <- lead >= lo & lead < hi                         # n x N_B: row in window of eruption j
  win[[k]] <- list(W    = W,
                   a    = apply(W, 2, function(w) if (any(w)) min(which(w)) else NA),
                   len  = colSums(W),
                   trig = rowSums(W) > 0)              # row is in any precursor window
}
erupt_table <- data.frame(
  label                   = eruptions$label,
  id                      = eruptions$id,
  class                   = eruptions$class,
  role                    = eruptions$role,
  start_UTC               = format(eruptions$start, "%Y-%m-%d %H:%M"),
  end_UTC                 = ifelse(eruptions$end > eruptions$start,
                                   format(eruptions$end, "%Y-%m-%d %H:%M"), ""),
  start_NZ                = format(with_tz(eruptions$start, "Pacific/Auckland"), "%Y-%m-%d %H:%M %Z"),
  window_end_before_start = ifelse(eruptions$in_record,
                                   format(t_grid[pmax(1, floor(pos_of(eruptions$start)))], "%Y-%m-%d"), ""),
  note                    = eruptions$note)
print(erupt_table, row.names = FALSE)
write_csv(erupt_table, file.path(out_dir, paste0("Eruption_catalogue_", station, ".csv")))

# ---------------------------------------------------------------------
# Step 2b: raw RSAM for the top panel (daily median and daily maximum,
#          placed on the same time axis). The daily values are saved to
#          RSAM_daily_<station>.csv in the results folder and reused next time.
# ---------------------------------------------------------------------

# time column -> POSIXct in UTC (NULL if none found)
rsam_time <- function(df) {
  nm <- names(df)
  # 1. unix seconds (or milliseconds)
  u <- nm[grepl("unix", nm, ignore.case = TRUE) & sapply(df, is.numeric)]
  if (length(u) > 0) {
    x <- as.numeric(df[[u[1]]])
    if (median(x, na.rm = TRUE) > 1e11) x <- x / 1000
    return(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"))
  }
  # 2. a date-time column (an instant - its zone is already known)
  p <- nm[sapply(df, inherits, "POSIXct")]
  if (length(p) > 0) return(with_tz(df[[p[1]]], "UTC"))
  # 3. date-time text (the first text column that parses)
  for (cn in nm[sapply(df, is.character)]) {
    tz <- if (grepl("nz", cn, ignore.case = TRUE)) "Pacific/Auckland" else rsam_tz
    x  <- head(df[[cn]][!is.na(df[[cn]])], 100)
    ok <- suppressWarnings(!is.na(ymd_hms(x, tz = tz, quiet = TRUE)))
    if (length(x) > 0 && mean(ok) > 0.9) {
      return(with_tz(ymd_hms(df[[cn]], tz = tz, quiet = TRUE), "UTC"))
    }
  }
  NULL
}

# RSAM column name (NULL if none)
rsam_value_col <- function(df) {
  if (rsam_column %in% names(df)) return(rsam_column)
  r <- names(df)[grepl("^rsam", names(df), ignore.case = TRUE) & sapply(df, is.numeric)]
  if (length(r) > 0) r[1] else NULL
}

# one table -> daily median, maximum and count (NULL if it holds no RSAM)
rsam_to_daily <- function(df, label) {
  df <- as.data.frame(df)
  vc <- rsam_value_col(df)
  tt <- rsam_time(df)
  if (is.null(tt)) {
    # per-day files without a full timestamp: take the day from the name, e.g. 2007.001
    m <- regmatches(label, regexpr("(19|20)[0-9]{2}\\.[0-9]{3}", label))
    if (length(m) == 1) {
      d0 <- as.Date(sprintf("%s-01-01", substr(m, 1, 4))) + as.integer(substr(m, 6, 8)) - 1
      tt <- rep(as.POSIXct(d0, tz = "UTC") + 12 * 3600, nrow(df))
    }
  }
  if (is.null(vc) || is.null(tt)) {
    cols <- paste0(names(df), " <", sapply(df, function(z) class(z)[1]), ">")
    if (length(cols) > 8) cols <- c(head(cols, 8), "...")
    return(structure(list(), why = paste0(label, ": columns ", paste(cols, collapse = ", "))))
  }
  v  <- as.numeric(df[[vc]])
  ok <- !is.na(tt) & !is.na(v)
  if (rsam_log) ok <- ok & v > 0
  if (!any(ok)) return(NULL)
  day <- as.Date(tt[ok], tz = "UTC")
  data.frame(day = as.Date(names(tapply(v[ok], day, median))),
             med = as.numeric(tapply(v[ok], day, median)),
             max = as.numeric(tapply(v[ok], day, max)),
             n   = as.numeric(tapply(v[ok], day, length)))
}

# read the source and return one row per day with data
rsam_read_source <- function() {
  if (!file.exists(rsam_source)) {
    cat("RSAM source not found:", rsam_source, "- RSAM panel left empty\n")
    return(NULL)
  }
  
  if (dir.exists(rsam_source)) {
    files <- list.files(rsam_source, pattern = "\\.csv$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
    files <- files[grepl(station, basename(files), fixed = TRUE)]
    files <- files[normalizePath(files, mustWork = FALSE) !=
                     normalizePath(file.path(data_dir, file_name), mustWork = FALSE)]
    # keep files of the years in the record (when the year is in the file name)
    yr    <- suppressWarnings(as.integer(sub(".*?((19|20)[0-9]{2}).*", "\\1", basename(files))))
    yrs   <- as.integer(format(t_grid[1] - 7 * 86400, "%Y")):as.integer(format(t_grid[n], "%Y"))
    files <- files[is.na(yr) | yr %in% yrs]
    if (length(files) == 0) {
      cat("No CSV files with", station, "in their name under", rsam_source, "- RSAM panel left empty\n")
      return(NULL)
    }
    cat("Reading RSAM from", length(files), "file(s) under", rsam_source, "\n")
    parts <- vector("list", length(files))
    for (i in seq_along(files)) {
      parts[[i]] <- rsam_to_daily(read_csv(files[i], show_col_types = FALSE, progress = FALSE),
                                  basename(files[i]))
      if (i %% 500 == 0) cat("  ...", i, "of", length(files), "\n")
    }
  } else {
    ext <- tolower(tools::file_ext(rsam_source))
    obj <- switch(ext,
                  rds = readRDS(rsam_source),
                  rda = , rdata = {
                    e <- new.env()
                    load(rsam_source, envir = e)
                    tabs <- Filter(is.data.frame, mget(ls(e), envir = e))
                    if (length(tabs) == 0) stop("No data frame inside ", rsam_source)
                    tabs[[which.max(sapply(tabs, nrow))]]
                  },
                  csv = read_csv(rsam_source, show_col_types = FALSE),
                  stop("Unknown cache type '.", ext, "' - use .rds, .RData or .csv"))
    if (!is.data.frame(obj) && is.list(obj)) obj <- Filter(is.data.frame, obj)   # list of tables
    if (is.data.frame(obj)) obj <- list(obj)
    parts <- lapply(seq_along(obj), function(i) rsam_to_daily(obj[[i]], basename(rsam_source)))
  }
  
  skipped <- sapply(parts, function(p) !is.null(attr(p, "why")))
  if (any(skipped)) {
    cat(sum(skipped), "table(s) without a recognisable time or RSAM column were skipped, e.g.\n  ",
        attr(parts[[which(skipped)[1]]], "why"), "\n")
  }
  parts <- parts[!skipped & !sapply(parts, is.null)]
  if (length(parts) == 0) {
    cat("No RSAM values found - RSAM panel left empty\n")
    return(NULL)
  }
  d <- do.call(rbind, parts)
  # a day split over two files: keep the overall maximum and the median of the larger part
  d <- d[order(d$day, -d$n), ]
  mx <- tapply(d$max, d$day, max)
  d  <- d[!duplicated(d$day), ]
  d$max <- as.numeric(mx[as.character(d$day)])
  d[, c("day", "med", "max")]
}

read_rsam <- function() {
  cache <- file.path(out_dir, paste0("RSAM_daily_", station, ".csv"))
  if (!rsam_refresh && file.exists(cache)) {
    d <- read_csv(cache, show_col_types = FALSE, col_types = cols(day = col_date()))
    cat("RSAM: daily values read from", cache, "(set rsam_refresh <- TRUE to rebuild)\n")
  } else {
    d <- rsam_read_source()
    if (is.null(d)) return(NULL)
    write_csv(d, cache)
    cat("RSAM: daily values saved to", cache, "\n")
  }
  d <- as.data.frame(d)
  
  # one row per calendar day, NA on days without data (so lines break at gaps)
  all_days <- seq(min(d$day), max(d$day), by = "day")
  m <- match(all_days, d$day)
  d <- data.frame(day = all_days, med = d$med[m], max = d$max[m])
  d$pos <- as.numeric(difftime(as.POSIXct(d$day, tz = "UTC") + 12 * 3600, t0,
                               units = "days")) / step_days + 1
  cat("RSAM:", sum(!is.na(d$med)), "days with data,", format(min(d$day)), "to", format(max(d$day)), "\n")
  d
}
rsam_day <- read_rsam()

# ---------------------------------------------------------------------
# ECA core: events come from ev_fun(thr) = TRUE/FALSE for every grid row
# ---------------------------------------------------------------------
eca_core <- function(ev_fun) {
  settings <- grid
  settings$days_beyond <- NA     # rows beyond the threshold
  settings$n_events    <- NA     # episodes (runs of consecutive rows)
  settings$hits_B      <- NA     # eruptions with an event row in their window
  settings$hits_A      <- NA     # episodes followed by an eruption
  null_B  <- matrix(0L, n_shift, n_set)
  ev_list <- vector("list", n_set)
  s_vec   <- 1:n_shift
  
  for (k in 1:n_set) {
    ev <- ev_fun(settings$thr[k])
    ev[is.na(ev) | mask] <- FALSE                    # masked rows are never events
    ev_list[[k]] <- ev
    wk <- win[[k]]
    
    # precursor hits (real)
    settings$hits_B[k] <- sum(colSums(wk$W & ev) > 0)
    
    # episodes and trigger hits (real)
    start  <- ev & !c(FALSE, head(ev, -1))
    run_id <- cumsum(start)
    settings$days_beyond[k] <- sum(ev)
    settings$n_events[k]    <- sum(start)
    settings$hits_A[k]      <- if (any(ev)) sum(tapply(wk$trig[ev], run_id[ev], any)) else 0
    
    # precursor hits after circular shift by s rows (all shifts at once)
    cs <- c(0, cumsum(c(ev, ev)))
    for (j in 1:N_B) {
      if (is.na(wk$a[j])) next
      st <- ((wk$a[j] - 1 - s_vec) %% n) + 1
      null_B[, k] <- null_B[, k] + as.integer((cs[st + wk$len[j]] - cs[st]) > 0)
    }
  }
  
  # p-value of every setting
  settings$p <- (1 + colSums(null_B >= rep(settings$hits_B, each = n_shift))) / (n_shift + 1)
  settings$precursor_rate <- round(settings$hits_B / N_B, 2)
  settings$trigger_rate   <- round(settings$hits_A / settings$n_events, 3)
  
  # p-value corrected for picking the best of the settings
  p_null <- matrix(NA_real_, n_shift, n_set)
  for (k in 1:n_set) {
    ge <- sapply(0:N_B, function(h) sum(null_B[, k] >= h))   # ge[h + 1] = shifts with >= h hits
    p_null[, k] <- (ge[null_B[, k] + 1] +
                      (settings$hits_B[k] >= null_B[, k])) / (n_shift + 1)
  }
  best_null <- apply(p_null, 1, min)
  b         <- which.min(settings$p)
  p_grid    <- (1 + sum(best_null <= settings$p[b])) / (n_shift + 1)
  
  list(settings = settings, b = b, p_grid = p_grid, ev_list = ev_list)
}

# table of the episodes of one setting - to check every event by date
episode_table <- function(ev, y, k, dir) {
  start <- which(ev & !c(FALSE, head(ev, -1)))
  end   <- which(ev & !c(tail(ev, -1), FALSE))
  if (length(start) == 0) return(data.frame())
  W <- win[[k]]$W
  data.frame(
    first_window_end = format(t_grid[start], date_fmt),
    last_window_end  = format(t_grid[end], date_fmt),
    n_windows        = end - start + 1,
    extreme_value    = mapply(function(a, z) {
      if (dir == "upper") max(y[a:z], na.rm = TRUE) else min(y[a:z], na.rm = TRUE)
    }, start, end),
    precursor_of     = sapply(seq_along(start), function(i) {
      js <- which(colSums(W[start[i]:end[i], , drop = FALSE]) > 0)
      if (length(js) > 0) paste(format(onset[js], date_fmt), collapse = "; ") else ""
    }))
}

# ---------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------
tick_dates <- seq(floor_date(t_grid[1], "year"), t_grid[n], by = tick_by)
tick_dates <- tick_dates[tick_dates >= t_grid[1]]
tick_pos   <- as.numeric(difftime(tick_dates, t0, units = "days")) / step_days + 1

# x-axis: a tick at every tick_by, a date label on every 'every'-th tick
draw_axis <- function(labels = TRUE, cex_lab = 0.75, every = 1) {
  axis(1, at = tick_pos, labels = FALSE)
  if (labels) {
    sel  <- seq(1, length(tick_pos), by = every)
    y_in <- grconvertY(par("usr")[3], "user", "inches") - 0.12
    text(tick_pos[sel], grconvertY(y_in, "inches", "user"), format(tick_dates[sel], date_fmt),
         srt = 45, adj = c(1, 1), xpd = NA, cex = cex_lab)
  }
}

# every threshold of the grid, the best one thicker;
# labels are pushed apart when the lines are close together
draw_thresholds <- function(x, best_thr, cex_lab = 0.75, dir = direction) {
  q  <- if (dir == "upper") thr_grid else 1 - thr_grid
  hv <- sapply(q, function(qq) quantile(x, qq, na.rm = TRUE))
  is_best <- sapply(thr_grid, function(t) isTRUE(all.equal(t, best_thr)))
  abline(h = hv, col = thr_cols, lty = 2, lwd = ifelse(is_best, 2.2, 1))
  
  gap <- abs(diff(grconvertY(c(0, 0.13 * cex_lab / 0.75), "inches", "user")))
  ord <- order(hv)
  ly  <- hv[ord]
  for (i in seq_along(ly)[-1]) ly[i] <- max(ly[i], ly[i - 1] + gap)
  ly  <- ly - (mean(ly) - mean(hv[ord]))                 # keep the group centred
  lab_y <- numeric(length(hv))
  lab_y[ord] <- ly
  text(par("usr")[2], lab_y, paste0("P", round(q * 100)), pos = 4, xpd = NA,
       col = thr_cols, cex = cex_lab, font = ifelse(is_best, 2, 1))
}

# precursor windows of setting k (grey bands, outlined so short windows stay visible)
win_fill   <- adjustcolor("grey55", 0.35)
win_border <- "grey45"
draw_windows <- function(k) {
  usr <- par("usr")
  usr[3:4] <- if (par("ylog")) 10^usr[3:4] else usr[3:4]
  wk  <- win[[k]]
  for (j in 1:N_B) {
    if (is.na(wk$a[j])) next
    rect(wk$a[j] - 0.5, usr[3], wk$a[j] + wk$len[j] - 0.5, usr[4],
         col = win_fill, border = win_border, lwd = 0.6)
  }
}

# y value at a fraction of the plot height (works on log axes too)
npc_y <- function(frac) {
  u <- par("usr")[3:4]
  y <- u[1] + frac * diff(u)
  if (par("ylog")) 10^y else y
}

# everything that is labelled in the RSAM panel: counted eruptions, then the others
oth_pos  <- if (nrow(OTH) > 0) pos_of(OTH$start) else numeric(0)
oth_end  <- if (nrow(OTH) > 0) pos_of(OTH$end)   else numeric(0)
oth_col  <- ifelse(OTH$role == "neutral", "darkorange3", "grey45")
lab_pos  <- c(erupt_pos, oth_pos)
lab_txt  <- c(erupt_lab, OTH$label)
lab_col  <- c(rep("blue", N_B), oth_col)

# label row of every label: labels closer than min_gap (fraction of the time axis)
# go to the next row, so labels of clustered events do not overlap
label_tiers <- function(pos, min_gap = 0.02) {
  gap  <- min_gap * n
  tier <- integer(length(pos))
  last <- numeric(0)
  for (j in order(pos)) {
    t <- which(pos[j] - last > gap)[1]
    if (is.na(t)) {
      last <- c(last, -Inf)
      t <- length(last)
    }
    tier[j] <- t
    last[t] <- pos[j]
  }
  tier
}
lab_tier   <- label_tiers(lab_pos)
label_top  <- 0.95                                   # first label row (fraction of panel height)
label_step <- 0.10                                   # distance between label rows
data_top   <- label_top - label_step * max(lab_tier) - 0.02   # RSAM data stay below the labels

# y range of the panel in data units (works on log axes too)
panel_y <- function() {
  u <- par("usr")[3:4]
  if (par("ylog")) 10^u else u
}

# eruptions: blue line at the onset, blue band over an eruptive period,
# thin dotted blue lines at the further eruptions inside a period;
# neutral events: orange dashed (band if a period); background: grey dotted.
# With numbers = TRUE every line is labelled (E1, E2, ..., N1, B1, ...).
draw_eruptions <- function(lwd = 1.2, numbers = FALSE, cex_num = 0.6) {
  yy <- panel_y()
  if (length(oth_pos) > 0) {
    per <- oth_end > oth_pos + 0.5
    if (any(per)) rect(oth_pos[per], yy[1], oth_end[per], yy[2],
                       col = adjustcolor(oth_col[per], 0.15), border = NA)
    abline(v = oth_pos, col = oth_col, lwd = lwd * 0.8,
           lty = ifelse(OTH$role == "neutral", 2, 3))
  }
  per <- erupt_end_pos > erupt_pos + 0.5
  if (any(per)) rect(erupt_pos[per], yy[1], erupt_end_pos[per], yy[2],
                     col = adjustcolor("blue", 0.15), border = NA)
  if (length(sub_pos) > 0) abline(v = sub_pos, col = "blue", lwd = lwd * 0.6, lty = 3)
  abline(v = erupt_pos, col = "blue", lwd = lwd)
  if (!numbers) return(invisible())
  text(lab_pos, npc_y(label_top - label_step * (lab_tier - 1)), lab_txt,
       pos = 4, offset = 0.15, col = lab_col, cex = cex_num, font = 2)
}

# short text of the events, e.g. "E3 2013-10-04 to 10-11"
event_text <- function(lab, s, e) {
  paste0(lab, " ", format(s, "%Y-%m-%d"),
         ifelse(e > s, paste0(" to ", format(e, "%m-%d")), ""))
}

# panel title in the margin above the panel, so it never hides data
panel_label <- function(name, stats, two_lines = FALSE, bold_stats = FALSE, cex = 0.62) {
  if (two_lines) {
    mtext(name,  side = 3, adj = 0, line = 0.9,  cex = cex,        font = 2)
    mtext(stats, side = 3, adj = 0, line = 0.12, cex = cex * 0.92, font = ifelse(bold_stats, 2, 1))
  } else {
    mtext(name,  side = 3, adj = 0, line = 0.15, cex = cex,        font = 2)
    mtext(stats, side = 3, adj = 1, line = 0.15, cex = cex * 0.95, font = ifelse(bold_stats, 2, 1))
  }
}

# top panel of the stacked figures: the raw RSAM
draw_rsam_panel <- function() {
  if (is.null(rsam_day)) {
    plot(1, type = "n", xlim = c(1, n), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    text(n / 2, 0.5, paste("No RSAM data found in rsam_source:", rsam_source), cex = 0.7)
  } else {
    rng <- range(c(rsam_day$med, rsam_day$max), na.rm = TRUE)
    if (rsam_log) {
      lr   <- log10(rng)
      ylim <- 10^c(lr[1], lr[1] + diff(lr) / data_top)   # headroom for the E labels
    } else {
      lo   <- min(0, rng[1])
      ylim <- c(lo, lo + (rng[2] - lo) / data_top)
    }
    plot(rsam_day$pos, rsam_day$max, type = "l", col = "grey65", lwd = 0.5,
         xlim = c(1, n), ylim = ylim, log = ifelse(rsam_log, "y", ""),
         xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    lines(rsam_day$pos, rsam_day$med, lwd = 0.6)
    if (rsam_log) {
      k <- ceiling(log10(rng[1])):floor(log10(rng[2]))             # ticks only where the data are
      if (length(k) >= 2) {
        axis(2, at = 10^k, labels = parse(text = paste0("10^", k)), las = 1, cex.axis = 0.7)
      } else {
        axis(2, las = 1, cex.axis = 0.7)
      }
    } else {
      axis(2, las = 1, cex.axis = 0.7)
    }
  }
  draw_axis(labels = FALSE)
  draw_eruptions(numbers = TRUE)
  panel_label(paste0(station, " RSAM (10-minute values", ifelse(rsam_log, ", log scale", ""),
                     "): black = daily median, grey = daily maximum"),
              paste0("Eruptions counted (UTC, scenario '", scenario, "'): ",
                     paste(event_text(erupt_lab, onset, onset_end), collapse = ",  ")),
              two_lines = TRUE)
}

# colour key of the stacked figures (two lines in the outer top margin)
stacked_key <- function(what = "windows") {
  oth_txt <- if (nrow(OTH) > 0) {
    paste0("      Orange dashed: left out of the test, grey dotted: not treated as eruptions (",
           paste(event_text(OTH$label, OTH$start, OTH$end), collapse = ", "), ")")
  } else ""
  mtext(paste0("Blue: eruptions E1-E", N_B, " (line = onset, band = eruptive period)      ",
               "Grey bands: precursor window of the best setting ",
               "(dT weeks, ending tau weeks before each onset)"),
        side = 3, outer = TRUE, line = 2.0, cex = 0.62)
  mtext(substring(oth_txt, 7), side = 3, outer = TRUE, line = 1.15, cex = 0.62)
  mtext(paste0("Red: ", what, " beyond the best threshold (dark = inside a precursor window)      ",
               "Dashed: the four percentile thresholds (bold = best)"),
        side = 3, outer = TRUE, line = 0.3, cex = 0.62)
}

# rows beyond the threshold: small pale = outside a window, large red = inside (hit)
draw_beyond <- function(ev, y, k, cex = 1) {
  inw <- ev & win[[k]]$trig
  out <- which(ev & !inw)
  hit <- which(inw)
  points(out, y[out], pch = 16, col = adjustcolor("red", 0.45), cex = 0.45 * cex)
  points(hit, y[hit], pch = 16, col = "red3", cex = 0.9 * cex)
}

# legend in one row between the title and the plot, so it never covers data
event_legend <- function(cex = 0.7) {
  legend(x = mean(par("usr")[1:2]), y = par("usr")[4], xjust = 0.5, yjust = 0,
         ncol = 3, xpd = NA, bty = "n", cex = cex,
         legend = c("beyond the best threshold",
                    "beyond it inside a precursor window (hit)",
                    "precursor window (best dT, tau)",
                    "eruption onset (band = eruptive period)",
                    "left out of the test (not certainly an eruption)",
                    "not treated as an eruption",
                    "thresholds (bold = best)"),
         pch    = c(16, 16, 15, NA, NA, NA, NA),
         col    = c(adjustcolor("red", 0.45), "red3", win_fill, "blue", "darkorange3", "grey45", "grey40"),
         pt.cex = c(0.6, 1, 1.6, NA, NA, NA, NA),
         lty    = c(NA, NA, NA, 1, 2, 3, 2),
         lwd    = c(NA, NA, NA, 2, 1.6, 1.6, 1))
}

# ---------------------------------------------------------------------
# Step 3: loop over all features
# ---------------------------------------------------------------------
summary_all <- NULL
panel_data  <- list()

for (f in all_features) {
  
  dir_f <- feat_dir[[f]]
  cat("\n=====", f, "(", dir_f, ") =====\n")
  x  <- X[, f]
  xq <- x
  xq[mask] <- NA                  # thresholds come from the rows that take part in the test
  
  ev_fun <- function(thr) {
    if (dir_f == "upper") !is.na(xq) & xq > quantile(xq, thr, na.rm = TRUE)
    else                  !is.na(xq) & xq < quantile(xq, 1 - thr, na.rm = TRUE)
  }
  res      <- eca_core(ev_fun)
  settings <- res$settings
  b        <- res$b
  p_grid   <- res$p_grid
  ev_best  <- res$ev_list[[b]]
  
  cat("best: thr", settings$thr[b], "| dT", settings$dT[b], "| tau", settings$tau[b],
      "| hits", settings$hits_B[b], "of", N_B,
      "| p =", round(settings$p[b], 4), "| corrected =", round(p_grid, 4), "\n")
  
  # grid table and episode table of the best setting
  out   <- settings
  out$p <- round(out$p, 4)
  write_csv(out, file.path(out_dir, paste0("ECA_grid_", station, "_", f, ".csv")))
  write_csv(episode_table(ev_best, x, b, dir_f),
            file.path(out_dir, paste0("ECA_episodes_", station, "_", f, ".csv")))
  
  # figure of this feature
  pdf(file.path(out_dir, paste0("ECA_", station, "_", f, ".pdf")), width = 11, height = 4.5)
  par(mar = c(6, 4, 5.2, 4))
  plot(seq_len(n), x, type = "l", lwd = 0.6, xlab = "", ylab = f, xaxt = "n")
  title(main = paste0(station, " - ", f, " [", dir_f, "]  (threshold ", settings$thr[b],
                      ", dT = ", settings$dT[b], " wk, tau = ", settings$tau[b],
                      " wk, hits ", settings$hits_B[b], "/", N_B,
                      ", corrected p = ", round(p_grid, 3), ")"), line = 3.4)
  draw_windows(b)
  lines(seq_len(n), x, lwd = 0.6)
  draw_axis()
  mtext("window end date (UTC)", side = 1, line = 4.5)
  draw_thresholds(xq, settings$thr[b], dir = dir_f)
  draw_beyond(ev_best, x, b)
  draw_eruptions(lwd = 2)
  event_legend()
  dev.off()
  
  panel_data[[f]] <- list(x = x, xq = xq, ev = ev_best, k = b, thr = settings$thr[b], dir = dir_f,
                          dT = settings$dT[b], tau = settings$tau[b],
                          hits = settings$hits_B[b], p_grid = p_grid)
  
  summary_all <- rbind(summary_all, data.frame(
    feature        = f,
    direction      = dir_f,
    thr            = settings$thr[b],
    dT             = settings$dT[b],
    tau            = settings$tau[b],
    days_beyond    = settings$days_beyond[b],
    n_events       = settings$n_events[b],
    hits_B         = settings$hits_B[b],
    hits_A         = settings$hits_A[b],
    precursor_rate = settings$precursor_rate[b],
    trigger_rate   = settings$trigger_rate[b],
    best_p         = round(settings$p[b], 4),
    p_grid         = round(p_grid, 4)))
}

# ---------------------------------------------------------------------
# Step 4: Benjamini-Hochberg over all features
# ---------------------------------------------------------------------
summary_all$q_BH <- round(p.adjust(summary_all$p_grid, method = "BH"), 4)
summary_all <- summary_all[order(summary_all$p_grid), ]
print(summary_all, row.names = FALSE)
write_csv(summary_all, file.path(out_dir, paste0("ECA_summary_", station, "_all_features.csv")))

# ---------------------------------------------------------------------
# Step 5: summary figure - corrected p-value of every feature
# ---------------------------------------------------------------------
draw_summary <- function() {
  par(mar = c(4, 9, 3, 1))
  barplot(rev(-log10(summary_all$p_grid)), names.arg = rev(summary_all$feature),
          horiz = TRUE, las = 1, col = ifelse(rev(summary_all$q_BH) < 0.05, "red", "grey60"),
          xlab = "-log10(corrected p)",
          main = paste(station, "- ECA, best setting per feature"))
  abline(v = -log10(0.05), lty = 2)
}
pdf(file.path(out_dir, paste0("ECA_summary_", station, "_all_features.pdf")), width = 8, height = 6)
draw_summary()
dev.off()

# ---------------------------------------------------------------------
# Step 6: combined figure - all features stacked
#   ncol = 1 : one column on one tall page;  ncol = 2 : A4 page
# ---------------------------------------------------------------------
draw_panels <- function(ncol = 1, order_by_p = FALSE) {
  feats <- if (order_by_p) summary_all$feature else all_features
  nf    <- length(feats)
  nrow  <- ceiling(nf / ncol)
  old <- par(no.readonly = TRUE)
  on.exit({ layout(1); par(old) })
  
  # RSAM across the full width on top, then the features
  cells <- c(seq_len(nf) + 1, rep(0, nrow * ncol - nf))
  layout(rbind(rep(1, ncol), matrix(cells, nrow = nrow, ncol = ncol, byrow = TRUE)),
         heights = c(1.5, rep(1, nrow)))
  par(oma = c(6, 0, 5.0, 0), mgp = c(2.5, 0.45, 0), tcl = -0.25, cex = 0.9)
  
  par(mar = c(0.4, 4, 1.9, 2.6))
  draw_rsam_panel()
  
  two_lines <- ncol > 1
  par(mar = c(0.4, 4, ifelse(two_lines, 1.9, 1.15), 2.6))
  for (i in seq_len(nf)) {
    f <- feats[i]
    d <- panel_data[[f]]
    q <- summary_all$q_BH[summary_all$feature == f]
    plot(seq_len(n), d$x, type = "l", lwd = 0.5, xaxt = "n", xlab = "", ylab = "",
         las = 1, cex.axis = 0.7)
    draw_windows(d$k)
    lines(seq_len(n), d$x, lwd = 0.5)
    draw_axis(labels = (i > nf - ncol), cex_lab = 0.6, every = ncol)
    draw_thresholds(d$xq, d$thr, cex_lab = 0.55, dir = d$dir)
    draw_beyond(d$ev, d$x, d$k, cex = 0.6)
    draw_eruptions()
    panel_label(paste0(f, " [", d$dir, "]"),
                paste0("thr ", d$thr, ", dT ", d$dT, " wk, tau ", d$tau, " wk, hits ", d$hits,
                       "/", N_B, ", p_grid ", round(d$p_grid, 3), ", q ", round(q, 3)),
                two_lines = two_lines, bold_stats = (q < 0.05),
                cex = ifelse(two_lines, 0.55, 0.62))
  }
  mtext("date (UTC); every feature value is plotted at the end of its 7-day window",
        side = 1, outer = TRUE, line = 4.5, cex = 0.8)
  mtext(paste0(station, " - ECA, all features (", dir_note, ")"),
        side = 3, outer = TRUE, line = 3.1, font = 2)
  stacked_key("windows")
}

pdf(file.path(out_dir, paste0("ECA_", station, "_all_features_1col.pdf")),
    width = 11, height = 1.6 * (length(all_features) + 1.5) + 1.8)
draw_panels(ncol = 1)
dev.off()

pdf(file.path(out_dir, paste0("ECA_", station, "_all_features_2col.pdf")),
    width = 8.27, height = 11.69)
draw_panels(ncol = 2)
dev.off()

# =====================================================================
# Step 7: composite features
#   rolling percentile -> correlation groups -> median and vote composites
#   -> same ECA test -> BH over the composites
# Groups come from the correlation structure only, never from ECA results.
# =====================================================================
roll_win_days <- 182    # trailing window for the rolling percentile (days, ~6 months)
use_rolling   <- TRUE   # FALSE = global percentile ranks (comparable with Step 3)
r_min         <- 0.7    # features joined at |Spearman r| >= r_min form a group
orient_manual <- c()    # optional manual orientation, e.g. c(H_Shannon = -1)

comp_dir <- file.path(out_dir, "composites")
dir.create(comp_dir, showWarnings = FALSE)

# 7a: percentile transform (rolling = share of the trailing window <= this row)
roll_rows <- max(2, round(roll_win_days / step_days))
roll_pct <- function(x, w, min_n = ceiling(w / 2)) {
  out <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    if (is.na(x[i])) next
    v <- x[max(1, i - w + 1):i]
    v <- v[!is.na(v)]
    if (length(v) >= min_n) out[i] <- mean(v <= x[i])
  }
  out
}
glob_pct <- function(x) rank(x, na.last = "keep", ties.method = "average") / sum(!is.na(x))

P <- sapply(all_features, function(f) {
  if (use_rolling) roll_pct(X[, f], roll_rows) else glob_pct(X[, f])
})
for (f in all_features) if (feat_dir[[f]] == "lower") P[, f] <- 1 - P[, f]

# 7b: correlation groups and orientation
cor_mat <- cor(X, method = "spearman", use = "pairwise.complete.obs")
hc      <- hclust(as.dist(1 - abs(cor_mat)), method = "average")
grp     <- cutree(hc, h = 1 - r_min)
groups  <- split(names(grp), grp)

write_csv(data.frame(feature = rownames(cor_mat), round(cor_mat, 3), check.names = FALSE),
          file.path(comp_dir, paste0("Feature_spearman_", station, ".csv")))

pdf(file.path(comp_dir, paste0("Feature_dendrogram_", station, ".pdf")), width = 8, height = 5)
par(mar = c(2, 4, 3, 1))
plot(hc, xlab = "", sub = "", ylab = "1 - |Spearman r|",
     main = paste(station, "- feature groups (average linkage)"))
abline(h = 1 - r_min, lty = 2, col = "red")
rect.hclust(hc, h = 1 - r_min, border = "grey40")
dev.off()

orient <- setNames(rep(1, length(all_features)), all_features)
ref    <- setNames(rep(NA_character_, length(all_features)), all_features)
for (g in groups) {
  for (f in g) {
    d_f   <- ifelse(feat_dir[[f]]    == "lower", -1, 1)
    d_ref <- ifelse(feat_dir[[g[1]]] == "lower", -1, 1)
    orient[f] <- ifelse(cor_mat[f, g[1]] * d_f * d_ref >= 0, 1, -1)
    ref[f]    <- g[1]
  }
}
orient[names(orient_manual)] <- orient_manual
if (any(orient == -1)) {
  cat("Note: flipped against their own direction inside their group:",
      names(orient)[orient == -1], "\n")
}
O <- P
for (f in all_features) if (orient[f] == -1) O[, f] <- 1 - P[, f]

group_table <- data.frame(
  feature     = all_features,
  direction   = feat_dir[all_features],
  group       = grp[all_features],
  reference   = ref[all_features],
  r_with_ref  = round(sapply(all_features, function(f) cor_mat[f, ref[f]]), 3),
  orientation = orient[all_features])
print(group_table, row.names = FALSE)
write_csv(group_table, file.path(comp_dir, paste0("Feature_groups_", station, ".csv")))

# 7c: composites of every group with >= 2 members
comp_info    <- list()
summary_comp <- NULL

for (k_g in seq_along(groups)) {
  g <- groups[[k_g]]
  if (length(g) < 2) {
    cat("Group", k_g, "has one member (", g, ") - already covered by Step 3\n")
    next
  }
  m         <- length(g)
  Og        <- O[, g, drop = FALSE]
  med       <- apply(Og, 1, median)         # NA unless every member is available
  min_votes <- floor(m / 2) + 1
  
  for (type in c("median", "vote")) {
    name <- paste0("G", k_g, "_", type)
    cat("\n=====", name, ":", paste(g, collapse = ", "), "=====\n")
    if (type == "median") {
      medq <- med
      medq[mask] <- NA
      res <- eca_core(function(thr) !is.na(medq) & medq > quantile(medq, thr, na.rm = TRUE))
    } else {
      res <- eca_core(function(thr) {
        v <- rowSums(Og > thr)
        !is.na(v) & v >= min_votes
      })
    }
    s <- res$settings
    b <- res$b
    cat("best: thr", s$thr[b], "| dT", s$dT[b], "| tau", s$tau[b],
        "| hits", s$hits_B[b], "of", N_B,
        "| p =", round(s$p[b], 4), "| corrected =", round(res$p_grid, 4), "\n")
    
    s_out   <- s
    s_out$p <- round(s_out$p, 4)
    write_csv(s_out, file.path(comp_dir, paste0("ECA_grid_", station, "_", name, ".csv")))
    
    comp_info[[name]] <- list(type = type, members = g, Og = Og, med = med,
                              min_votes = min_votes, res = res)
    summary_comp <- rbind(summary_comp, data.frame(
      composite = name, type = type, members = paste(g, collapse = "+"),
      min_votes = ifelse(type == "vote", min_votes, NA),
      thr = s$thr[b], dT = s$dT[b], tau = s$tau[b],
      days_beyond = s$days_beyond[b], n_events = s$n_events[b],
      hits_B = s$hits_B[b], hits_A = s$hits_A[b],
      precursor_rate = s$precursor_rate[b], trigger_rate = s$trigger_rate[b],
      best_p = round(s$p[b], 4), p_grid = round(res$p_grid, 4)))
  }
}

# 7d: BH over the composites + comparison with the single features
if (!is.null(summary_comp)) {
  summary_comp$q_BH <- round(p.adjust(summary_comp$p_grid, method = "BH"), 4)
  summary_comp <- summary_comp[order(summary_comp$p_grid), ]
  print(summary_comp, row.names = FALSE)
  write_csv(summary_comp, file.path(comp_dir, paste0("ECA_summary_", station, "_composites.csv")))
  
  comparison <- rbind(
    data.frame(name = summary_all$feature, set = "single (Step 3)",
               thr = summary_all$thr, dT = summary_all$dT, tau = summary_all$tau,
               hits_B = summary_all$hits_B, n_events = summary_all$n_events,
               p_grid = summary_all$p_grid, q_BH_within_set = summary_all$q_BH),
    data.frame(name = paste0(summary_comp$composite, " (", summary_comp$members, ")"),
               set = paste0("composite (", ifelse(use_rolling,
                                                  paste0("rolling ", roll_win_days, " d"),
                                                  "global"), ")"),
               thr = summary_comp$thr, dT = summary_comp$dT, tau = summary_comp$tau,
               hits_B = summary_comp$hits_B, n_events = summary_comp$n_events,
               p_grid = summary_comp$p_grid, q_BH_within_set = summary_comp$q_BH))
  write_csv(comparison, file.path(comp_dir, paste0("ECA_comparison_", station, ".csv")))
}

# 7e: stacked figure of the composites
draw_comp_panels <- function() {
  nm  <- names(comp_info)
  nc  <- length(nm)
  old <- par(no.readonly = TRUE)
  on.exit({ layout(1); par(old) })
  layout(matrix(1:(nc + 1), ncol = 1), heights = c(1.5, rep(1, nc)))
  par(oma = c(6, 0, 5.0, 0), mgp = c(2.5, 0.45, 0), tcl = -0.25, cex = 0.9)
  
  par(mar = c(0.4, 4, 1.9, 2.6))
  draw_rsam_panel()
  
  par(mar = c(0.4, 4, 1.9, 2.6))
  for (i in seq_len(nc)) {
    ci <- comp_info[[nm[i]]]
    s  <- ci$res$settings
    b  <- ci$res$b
    ev <- ci$res$ev_list[[b]]
    q  <- summary_comp$q_BH[summary_comp$composite == nm[i]]
    
    if (ci$type == "median") {
      y <- ci$med
      plot(seq_len(n), y, type = "l", lwd = 0.5, ylim = c(0, 1),
           xaxt = "n", xlab = "", ylab = "", las = 1, cex.axis = 0.7)
      draw_windows(b)
      lines(seq_len(n), y, lwd = 0.5)
      yq <- y
      yq[mask] <- NA
      draw_thresholds(yq, s$thr[b], cex_lab = 0.55, dir = "upper")
      lab_rule <- "median of oriented percentiles"
    } else {
      y <- rowSums(ci$Og > s$thr[b])
      plot(seq_len(n), y, type = "s", lwd = 0.5, ylim = c(0, length(ci$members)),
           xaxt = "n", xlab = "", ylab = "", las = 1, cex.axis = 0.7)
      draw_windows(b)
      lines(seq_len(n), y, type = "s", lwd = 0.5)
      col_b <- thr_cols[match(s$thr[b], thr_grid)]
      abline(h = ci$min_votes - 0.5, lty = 2, lwd = 2, col = col_b)
      text(par("usr")[2], ci$min_votes - 0.5, paste0(">=", ci$min_votes),
           pos = 4, xpd = NA, col = col_b, cex = 0.55, font = 2)
      lab_rule <- paste0("votes: members above P", round(s$thr[b] * 100))
    }
    draw_beyond(ev, y, b, cex = 0.6)
    draw_eruptions()
    draw_axis(labels = (i == nc), cex_lab = 0.6)
    panel_label(paste0(nm[i], " [", paste(ci$members, collapse = ", "), "]"),
                paste0(lab_rule, "   thr ", s$thr[b], ", dT ", s$dT[b], " wk, tau ", s$tau[b],
                       " wk, hits ", s$hits_B[b], "/", N_B,
                       ", p_grid ", round(ci$res$p_grid, 3), ", q ", round(q, 3)),
                two_lines = TRUE, bold_stats = (q < 0.05))
  }
  mtext("date (UTC); every value is plotted at the end of its 7-day window",
        side = 1, outer = TRUE, line = 4.5, cex = 0.8)
  mtext(paste0(station, " - ECA, composite features (", dir_note, "; ",
               ifelse(use_rolling, paste0("rolling ", roll_win_days, "-day percentiles"),
                      "global percentiles"), ")"),
        side = 3, outer = TRUE, line = 3.1, font = 2)
  stacked_key("rows")
}

if (length(comp_info) > 0) {
  pdf(file.path(comp_dir, paste0("ECA_", station, "_composites.pdf")),
      width = 11, height = 1.7 * (length(comp_info) + 1) + 1.8)
  draw_comp_panels()
  dev.off()
}

cat("\nAll tables and figures saved in:", out_dir, "\n")