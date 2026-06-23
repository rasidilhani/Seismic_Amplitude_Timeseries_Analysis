library(tidyverse)
library(lubridate)
library(StatOrdPattHxC)
library(writexl)
library(here)

# ── Parameters ────────────────────────────────────────────────────────────
# Only D = 5 is computed 
D       <- 5
z_alpha <- qnorm(0.975)
BETA    <- 1.5

years_to_use <- 2011:2022

# Both Shiny-app variables get run through the ordinal pattern pipeline.
variables <- c(
  Displacement = "displacement_avg_m",
  RSAM         = "rsam_avg"
)

# ── Load & timestamp the seismic data ──────────────────────────────────────
# Identical logic to the Shiny app: build the timestamp in UTC first, then
# convert to NZ time, so the "1970-01-01" origin is anchored to UTC rather
# than being mis-read as already being in Pacific/Auckland.
csv_files <- list.files("data", pattern = "^WIZ_NZ_[0-9]{4}\\.csv$", full.names = TRUE)

seismic_all <- list()
for (f in csv_files) {
  d <- read_csv(f, show_col_types = FALSE)
  utc_time <- as.POSIXct(d$unix_timestamp, origin = "1970-01-01", tz = "UTC")
  d$datetime_nz <- with_tz(utc_time, tzone = "Pacific/Auckland")
  seismic_all[[f]] <- d
}
seismic_all <- bind_rows(seismic_all)

# Year / month / quarter are all derived from the NZ-local timestamp, not the
# filename — same reasoning as in the Shiny app: deriving "year" from the
# filename would put late-Dec UTC readings (which land on Jan 1 NZ time) into
# the wrong year/month bucket.
seismic_all$year    <- year(seismic_all$datetime_nz)
seismic_all$month   <- month(seismic_all$datetime_nz, label = TRUE)   # ordered Jan..Dec
seismic_all$quarter <- paste0("Q", quarter(seismic_all$datetime_nz))

seismic_all <- seismic_all %>% filter(year %in% years_to_use)

# ── Helper functions ──────────
Jensen_Shannon <- function(p, q) {
  m  <- 0.5 * (p + q)
  js <- 0.5 * sum(ifelse(p == 0, 0, p * log((p + 1e-12) / (m + 1e-12)))) +
    0.5 * sum(ifelse(q == 0, 0, q * log((q + 1e-12) / (m + 1e-12))))
  js / log(2)
}

Fisher_Ferri <- function(p) {
  total <- 0
  for (i in 1:(length(p) - 1)) {
    total <- total + (p[i + 1] - p[i])^2 / (p[i + 1] + p[i] + 1e-12)
  }
  0.5 * total
}

# Computes the full set of H/C metrics, their variances, and their
# semi-lengths for one numeric series at embedding dimension D.
compute_HC <- function(series, D, n_i, z_alpha) {
  prob  <- OPprob(series, emb = D)
  n_eff <- n_i - D + 1
  Pe    <- rep(1 / length(prob), length(prob))
  JS    <- Jensen_Shannon(prob, Pe)
  
  Hs <- HShannon(prob)
  Hr <- HRenyi(prob, beta = BETA)
  Ht <- HTsallis(prob, beta = BETA)
  Hf <- Fisher_Ferri(prob)
  
  Cs <- StatComplexity(prob)
  Cr <- JS * Hr
  Ct <- JS * Ht
  Cf <- JS * Hf
  
  Var_Hs <- suppressWarnings(sigma2q(series, emb = D, ent = "S"))
  Var_Hr <- suppressWarnings(sigma2q(series, emb = D, ent = "R", beta = BETA))
  Var_Ht <- suppressWarnings(sigma2q(series, emb = D, ent = "T", beta = BETA))
  Var_Hf <- suppressWarnings(sigma2q(series, emb = D, ent = "F"))
  
  Var_HI <- suppressWarnings(asymptoticVarHShannonMultinomial(prob, n_eff))
  Var_CI <- suppressWarnings(varC(prob, n_eff))
  
  a_ratio <- ifelse(Var_HI > 0, Var_Hs / Var_HI, NA)
  Var_Cs  <- a_ratio * Var_CI
  
  semi <- function(v) ifelse(!is.finite(v) | v <= 0, NA, sqrt(v) / sqrt(n_eff) * z_alpha)
  
  tibble(
    N_points = n_i, N_eff = n_eff,
    
    H_Shannon = Hs, C_Shannon = Cs,
    H_Renyi   = Hr, C_Renyi   = Cr,
    H_Tsallis = Ht, C_Tsallis = Ct,
    H_Fisher  = Hf, C_Fisher  = Cf,
    
    Var_H_Shannon = Var_Hs, Var_C_Shannon = Var_Cs,
    Var_H_Renyi   = Var_Hr,
    Var_H_Tsallis = Var_Ht,
    Var_H_Fisher  = Var_Hf,
    
    Semi_H_Shannon = semi(Var_Hs), Semi_C_Shannon = semi(Var_Cs),
    Semi_H_Renyi   = semi(Var_Hr),
    Semi_H_Tsallis = semi(Var_Ht),
    Semi_H_Fisher  = semi(Var_Hf)
  )
}

# Pulls one variable's series out of a (year + period)-filtered data frame,
# drops non-finite values, and guards against periods that are empty or too
# short to support D = 5 (OPprob needs at least D points to form one pattern).
run_one_period <- function(df, var_col, D, z_alpha) {
  series <- df[[var_col]]
  series <- series[is.finite(series)]
  n_i <- length(series)
  
  if (n_i < D) {
    return(tibble(
      N_points = n_i, N_eff = NA_real_,
      H_Shannon = NA, C_Shannon = NA,
      H_Renyi = NA, C_Renyi = NA,
      H_Tsallis = NA, C_Tsallis = NA,
      H_Fisher = NA, C_Fisher = NA,
      Var_H_Shannon = NA, Var_C_Shannon = NA,
      Var_H_Renyi = NA, Var_H_Tsallis = NA, Var_H_Fisher = NA,
      Semi_H_Shannon = NA, Semi_C_Shannon = NA,
      Semi_H_Renyi = NA, Semi_H_Tsallis = NA, Semi_H_Fisher = NA
    ))
  }
  
  compute_HC(series, D, n_i, z_alpha)
}

# ── Build the long-format results table ────────────────────────────────────
# One row per Year x Variable x Period_Type x Period (Jan..Dec or Q1..Q4).
all_results <- list()

for (yr in years_to_use) {
  year_df <- seismic_all %>% filter(year == yr)
  
  for (var_label in names(variables)) {
    var_col <- variables[[var_label]]
    
    # Monthly periods
    for (mo in levels(year_df$month)) {
      sub_df <- year_df %>% filter(month == mo)
      res <- run_one_period(sub_df, var_col, D, z_alpha)
      all_results[[length(all_results) + 1]] <- res %>%
        mutate(Year = yr, Variable = var_label, Period_Type = "Monthly",
               Period = as.character(mo), .before = 1)
    }
    
    # Quarterly periods
    for (q in c("Q1", "Q2", "Q3", "Q4")) {
      sub_df <- year_df %>% filter(quarter == q)
      res <- run_one_period(sub_df, var_col, D, z_alpha)
      all_results[[length(all_results) + 1]] <- res %>%
        mutate(Year = yr, Variable = var_label, Period_Type = "Quarterly",
               Period = q, .before = 1)
    }
  }
}

period_order <- c(month.abb, "Q1", "Q2", "Q3", "Q4")

ordinal_results <- bind_rows(all_results) %>%
  mutate(Period = factor(Period, levels = period_order)) %>%
  arrange(Variable, Year, Period_Type, Period) %>%
  mutate(Period = as.character(Period))

# ── Save everything as one combined, long-format workbook ─────────────────
dir.create(here("Results"), recursive = TRUE, showWarnings = FALSE)
out_path <- here("Results", "WIZ_OrdinalPatterns_D5_Monthly_Quarterly.xlsx")

write_xlsx(ordinal_results, out_path)

cat("Saved combined ordinal pattern results to:", out_path, "\n")
cat("Total rows:", nrow(ordinal_results), "\n")