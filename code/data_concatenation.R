# ============================================================
# Combine daily RSAM CSV files into annual files
# Years: 2007-2019
# Stations: WIZ.NZ and WSRZ.NZ
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 1. Main data directory
# ------------------------------------------------------------

base_dir <- "E:/data"

# Output directory
output_dir <- file.path(base_dir, "Annual_Data")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ------------------------------------------------------------
# 2. Function to combine one year for one station
# ------------------------------------------------------------

combine_year_station <- function(year, station) {
  
  year_dir <- file.path(base_dir, as.character(year))
  
  # Find folders such as 2007.001, 2007.002, ...
  day_folders <- list.dirs(
    year_dir,
    recursive = FALSE,
    full.names = TRUE
  )
  
  # Keep only folders matching YYYY.DDD
  day_folders <- day_folders[
    grepl(
      paste0("^", year, "\\.[0-9]{3}$"),
      basename(day_folders)
    )
  ]
  
  # Sort by day of year
  day_folders <- day_folders[
    order(as.integer(sub(".*\\.", "", basename(day_folders))))
  ]
  
  yearly_data <- list()
  
  for (day_folder in day_folders) {
    
    day_name <- basename(day_folder)
    
    # Example:
    # day_name = "2007.001"
    
    doy <- as.integer(sub(".*\\.", "", day_name))
    
    # Convert year + DOY to actual date
    date_value <- as.Date(
      doy - 1,
      origin = paste0(year, "-01-01")
    )
    
    station_folder <- file.path(
      day_folder,
      station
    )
    
    # Skip if station folder does not exist
    if (!dir.exists(station_folder)) {
      next
    }
    
    # Find CSV file
    csv_files <- list.files(
      station_folder,
      pattern = "\\.csv$",
      full.names = TRUE
    )
    
    # If your files do not actually end in .csv,
    # use this instead:
    #
    # csv_files <- list.files(
    #   station_folder,
    #   pattern = "timeseries",
    #   full.names = TRUE
    # )
    
    if (length(csv_files) == 0) {
      next
    }
    
    # Usually there should be only one file
    csv_file <- csv_files[1]
    
    # Read daily data
    daily_data <- read_csv(
      csv_file,
      col_types = cols(.default = col_character()),
      show_col_types = FALSE
    )
    
    # Add identifying information
    daily_data <- daily_data %>%
      mutate(
        Date = date_value,
        Year = year,
        DOY = doy,
        Station = station,
        .before = 1
      ) %>%
      mutate(
        across(
          any_of(c("RSAM", "MF", "HF", "DSAR")),
          as.numeric
        )
      )
    yearly_data[[day_name]] <- daily_data
    
    cat(
      "Read:",
      year,
      day_name,
      station,
      "-",
      nrow(daily_data),
      "rows\n"
    )
  }
  
  # ----------------------------------------------------------
  # Combine all daily files
  # ----------------------------------------------------------
  
  if (length(yearly_data) == 0) {
    
    cat(
      "No data found for",
      year,
      station,
      "\n"
    )
    
    return(NULL)
  }
  
  annual_data <- bind_rows(yearly_data)
  
  # ----------------------------------------------------------
  # Save annual CSV
  # ----------------------------------------------------------
  
  station_short <- sub("\\.NZ", "", station)
  
  output_file <- file.path(
    output_dir,
    paste0(
      station_short,
      "_NZ_",
      year,
      "_Annual_RSAM.csv"
    )
  )
  
  write_csv(
    annual_data,
    output_file
  )
  
  cat(
    "\nSaved:",
    output_file,
    "\nRows:",
    nrow(annual_data),
    "\n\n"
  )
  
  return(annual_data)
}


# ------------------------------------------------------------
# 3. Years and stations
# ------------------------------------------------------------

years <- 2007:2022

stations <- c(
  "WIZ.NZ",
  "WSRZ.NZ"
)


# ------------------------------------------------------------
# 4. Run for every year and station
# ------------------------------------------------------------

for (yr in years) {
  
  for (st in stations) {
    
    cat(
      "\n====================================\n",
      "Processing:",
      yr,
      st,
      "\n====================================\n"
    )
    
    combine_year_station(
      year = yr,
      station = st
    )
  }
}


# ============================================================
# END
# ============================================================