# =====================================================================
# Replace negative variances with NA
# Reads the WIZ and WSRZ feature files, sets every negative value in a
# Var_* column to NA, and saves a new "_clean" file next to the original.
# The original files are not changed.
# =====================================================================

library(readr)

data_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"
stations <- c("WIZ", "WSRZ")

for (station in stations) {
  
  in_file  <- file.path(data_dir, paste0("HC_RSAM_", station, "_6daysoverlap_20070101_to_20191231.csv"))
  out_file <- file.path(data_dir, paste0("HC_RSAM_", station, "_6daysoverlap_20070101_to_20191231_clean.csv"))
  
  feat <- read_csv(in_file, show_col_types = FALSE)
  
  # all variance columns: Var_H_Shannon, Var_C_Shannon, Var_H_Renyi, ...
  var_cols <- grep("^Var_", names(feat), value = TRUE)
  
  cat("\n", station, "-", nrow(feat), "weeks\n")
  for (v in var_cols) {
    neg <- !is.na(feat[[v]]) & feat[[v]] < 0     # negative values in this column
    cat("  ", v, ": negative values replaced =", sum(neg), "\n")
    feat[[v]][neg] <- NA
  }
  
  write_csv(feat, out_file, na = "NA")
  cat("  saved:", out_file, "\n")
}
