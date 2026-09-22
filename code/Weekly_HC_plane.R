## ============================================================================
## Weekly H-C (Entropy-Complexity) Plane Plots -- WIZ vs WSRZ, RSAM
## ONE FACETED figure per station (Shannon, Tsallis, Renyi, Fisher together
## as 4 panels) -- 2 stations x {noCI, withCI} = 4 PDFs total.
library(readr)
library(tidyverse)
library(lubridate)
library(StatOrdPattHxC)
library(scales)
library(writexl)  # install.packages("writexl") if not already installed
library(ggrepel)   # install.packages("ggrepel") if not already installed

## ---- CONFIG ---------------------------------------------------------------
D <- 5   # <-- MUST match the embedding dimension used to generate the CSVs

# Number of extreme points to label per measure
N_LABEL <- 4

data_dir    <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/data"
results_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"
out_dir     <- file.path(results_dir, "HC_Plane_Plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

date_tag <- "20140701_to_20191231"   # used only for OUTPUT plot filenames below

# Exact input filenames, as given.
input_files <- c(
  WIZ  = file.path(data_dir, "HC_RSAM_WIZ_20140701_to_20191231.csv"),
  WSRZ = file.path(data_dir, "HC_RSAM_WSRZ_20140701_to_20191231.csv")
)

missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0) {
  stop("These input files were not found:\n  ", paste(missing_inputs, collapse = "\n  "))
}

measure_levels <- c("Shannon", "Tsallis", "Renyi", "Fisher")

## ---- Theoretical Shannon H-C boundary (StatOrdPattHxC) --------------------
data("LinfLsup")

bounds <- LinfLsup %>%
  filter(as.integer(as.character(Dimension)) == D)

bound_group_col <- "Side"

# Crop the theoretical boundary curve to the H/C range actually spanned by
# the data (plus a small margin).
crop_bounds_to_data <- function(bounds, shannon_df, pad_frac = 0.05) {
  
  H_range <- range(shannon_df$H, na.rm = TRUE)
  C_range <- range(shannon_df$C, na.rm = TRUE)
  
  H_pad <- diff(H_range) * pad_frac
  C_pad <- diff(C_range) * pad_frac
  if (!is.finite(H_pad) || H_pad == 0) H_pad <- 0.01
  if (!is.finite(C_pad) || C_pad == 0) C_pad <- 0.01
  
  bounds %>%
    filter(
      H >= H_range[1] - H_pad, H <= H_range[2] + H_pad,
      C >= max(0, C_range[1] - C_pad), C <= C_range[2] + C_pad
    ) %>%
    mutate(Measure = factor("Shannon", levels = measure_levels))
}

add_shannon_boundary <- function(p, bounds_cropped) {
  p +
    geom_line(
      data = bounds_cropped,
      aes(x = H, y = C, group = .data[[bound_group_col]]),
      color = "grey35",
      linewidth = 0.5,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
}

## ---- Common theme ----------------------------------------------------------
hc_theme <- theme_classic(base_family = "serif", base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

## ---- Reshape one station's weekly CSV into long H-C-CI format -----------
build_long_df <- function(path) {
  
  # Week_Start / Week_End arrive as ISO 8601 with a literal T/Z, e.g.
  # "2014-07-01T00:00:00Z". Force them to raw text on read so readr's own
  # datetime guesser never gets a chance to parse it differently than
  # expected, then parse explicitly with lubridate::ymd_hms(), which
  # handles the T/Z form natively -- this way it works the same whether
  # the column shows up as text or an already-parsed timestamp.
  df <- read_csv(path,
                 col_types = readr::cols(Week_Start = readr::col_character(),
                                         Week_End   = readr::col_character(),
                                         .default   = readr::col_guess()),
                 show_col_types = FALSE)
  
  df <- df %>%
    mutate(
      Week_Start = as.Date(lubridate::ymd_hms(Week_Start, tz = "UTC", quiet = TRUE)),
      Week_End   = as.Date(lubridate::ymd_hms(Week_End,   tz = "UTC", quiet = TRUE))
    )
  
  df_shannon <- df %>%
    transmute(
      Week_Index, Week_Start, Week_End,
      Measure = "Shannon",
      H = H_Shannon, C = C_Shannon,
      Semi_H = Semi_H_Shannon, Semi_C = Semi_C_Shannon
    )
  
  df_tsallis <- df %>%
    transmute(
      Week_Index, Week_Start, Week_End,
      Measure = "Tsallis",
      H = H_Tsallis, C = C_Tsallis,
      Semi_H = Semi_H_Tsallis, Semi_C = NA_real_
    )
  
  df_renyi <- df %>%
    transmute(
      Week_Index, Week_Start, Week_End,
      Measure = "Renyi",
      H = H_Renyi, C = C_Renyi,
      Semi_H = Semi_H_Renyi, Semi_C = NA_real_
    )
  
  df_fisher <- df %>%
    transmute(
      Week_Index, Week_Start, Week_End,
      Measure = "Fisher",
      H = Fisher_Info, C = C_Fisher,
      Semi_H = Semi_H_Fisher, Semi_C = NA_real_
    )
  
  df_long <- bind_rows(df_shannon, df_tsallis, df_renyi, df_fisher) %>%
    mutate(Measure = factor(Measure, levels = measure_levels)) %>%
    filter(is.finite(H), is.finite(C)) %>%
    mutate(
      Semi_H = ifelse(is.finite(Semi_H) & Semi_H > 0, Semi_H, NA_real_),
      Semi_C = ifelse(is.finite(Semi_C) & Semi_C > 0, Semi_C, NA_real_)
    ) %>%
    group_by(Measure) %>%
    mutate(
      H_norm = scales::rescale(H, to = c(0, 1)),
      C_norm = scales::rescale(C, to = c(0, 1)),
      score        = H_norm - C_norm,   # low score = low H/high C corner (Shannon/Tsallis/Renyi)
      score_fisher = H_norm + C_norm    # high score = high H/high C corner (Fisher only)
    ) %>%
    ungroup()
  
  df_long
}

## ---- Pick the extreme points to label, per Measure facet ------------------
get_label_points <- function(df_long, n = N_LABEL) {
  
  non_fisher <- df_long %>%
    filter(Measure != "Fisher") %>%
    group_by(Measure) %>%
    slice_min(score, n = n, with_ties = FALSE) %>%
    ungroup()
  
  fisher <- df_long %>%
    filter(Measure == "Fisher") %>%
    group_by(Measure) %>%
    slice_max(score_fisher, n = n, with_ties = FALSE) %>%
    ungroup()
  
  bind_rows(non_fisher, fisher)
}

## ---- Plot builder: ALL FOUR measures faceted, ONE station -----------------
make_hc_plot <- function(df_long, bounds_cropped, label_points, station_label, with_ci) {
  
  p <- ggplot(df_long, aes(x = H, y = C, color = Week_Index)) +
    geom_point(size = 1.6, alpha = 0.8) +
    scale_color_viridis_c(name = "Week", option = "D") +
    facet_wrap(vars(Measure), ncol = 2, scales = "free") +
    labs(
      title = bquote(.(station_label) ~ italic(H) %*% italic(C) ~ "plane"),
      x = expression(italic(H)),
      y = expression(italic(C))
    ) +
    hc_theme
  
  if (with_ci) {
    # Horizontal CI (H direction) -- available for every measure.
    p <- p +
      geom_segment(
        data = df_long %>% filter(is.finite(Semi_H)),
        aes(x = pmax(0, H - Semi_H), xend = H + Semi_H, y = C, yend = C),
        linewidth = 0.3, alpha = 0.6
      )
    # Vertical CI (C direction) -- Shannon facet only (see ASSUMPTIONS above).
    p <- p +
      geom_segment(
        data = df_long %>% filter(Measure == "Shannon", is.finite(Semi_C)),
        aes(x = H, xend = H, y = pmax(0, C - Semi_C), yend = C + Semi_C),
        linewidth = 0.3, alpha = 0.6
      )
  }
  
  # Theoretical boundary curve -- Shannon facet only.
  p <- add_shannon_boundary(p, bounds_cropped)
  
  # Ring the exact point being labeled, so the ggrepel target is never
  # ambiguous regardless of where the label text gets nudged.
  p <- p +
    geom_point(
      data = label_points,
      aes(x = H, y = C),
      inherit.aes = FALSE,
      shape = 21, size = 3, stroke = 1,
      color = "black", fill = NA
    )
  
  p <- p +
    geom_text_repel(
      data = label_points,
      aes(x = H, y = C, label = paste0("week ", Week_Index)),
      inherit.aes = FALSE,
      color = "black",
      size = 3,
      fontface = "bold",
      segment.size = 0.3,
      segment.color = "grey30",
      min.segment.length = 0,   # always draw the leader line
      box.padding = 0.4,
      point.padding = 0.3,
      max.overlaps = Inf,
      seed = 42
    )
  
  p
}

## ---- Run for each station, save one faceted plot with and without CI -----
for (station_name in names(input_files)) {
  
  message("Building H-C plane plot for: ", station_name)
  
  df_long <- build_long_df(input_files[[station_name]])
  
  if (nrow(df_long) == 0) {
    warning("No usable weeks found for ", station_name, " -- skipping.")
    next
  }
  
  bounds_cropped <- crop_bounds_to_data(
    bounds,
    df_long %>% filter(Measure == "Shannon")
  )
  
  label_points <- get_label_points(df_long, n = N_LABEL)
  
  p_noci <- make_hc_plot(df_long, bounds_cropped, label_points, station_name, with_ci = FALSE)
  print(p_noci)
  ggsave(
    filename = file.path(out_dir, sprintf("HC_Plane_%s_D%d_noCI_%s.pdf", station_name, D, date_tag)),
    plot = p_noci, width = 20, height = 16, units = "cm", device = "pdf"
  )
  
  p_ci <- make_hc_plot(df_long, bounds_cropped, label_points, station_name, with_ci = TRUE)
  print(p_ci)
  ggsave(
    filename = file.path(out_dir, sprintf("HC_Plane_%s_D%d_withCI_%s.pdf", station_name, D, date_tag)),
    plot = p_ci, width = 20, height = 16, units = "cm", device = "pdf"
  )
  
  message("  -> saved noCI and withCI figures for ", station_name)
  
  ## ---- Export the labelled extreme points to Excel ----
  ## Two sheets from the SAME set of labelled points (this script only
  ## builds one faceted plot now, not a separate single-Shannon plot, so
  ## both sheets come from `label_points` above): one with all 4 measures
  ## (N_LABEL per measure), one filtered to just the Shannon rows for
  ## convenience.
  export_cols <- c("Week_Index", "Week_Start", "Week_End", "Measure",
                   "H", "C", "Semi_H", "Semi_C", "H_norm", "C_norm",
                   "score", "score_fisher")
  
  excel_out <- list(
    All_Measures = label_points %>%
      select(any_of(export_cols)) %>%
      arrange(Measure, Week_Index),
    Shannon_Only = label_points %>%
      filter(Measure == "Shannon") %>%
      select(any_of(export_cols)) %>%
      arrange(Week_Index)
  )
  
  excel_path <- file.path(out_dir, sprintf("HC_extreme_weeks_%s_D%d_%s.xlsx", station_name, D, date_tag))
  writexl::write_xlsx(excel_out, path = excel_path)
  message("  -> saved labelled extreme-week points to ", excel_path)
}

message("Done.")