## ============================================================================
## Weekly H-C (Entropy-Complexity) Plane Plots -- RSAM only, 2007-2019
##
## PART 1: One FACETED figure (4 panels: Shannon, Tsallis, Renyi, Fisher).
## PART 2: One STANDALONE figure for SHANNON ONLY (not faceted).
## Both parts are produced WITH and WITHOUT confidence intervals (4 PDFs
## total: facet-noCI, facet-withCI, shannon-noCI, shannon-withCI).
##
## FIXES APPLIED IN THIS VERSION (vs. the version you pasted):
##  - Removed a trailing comma in `input_files <- c(..., )` in the second
##    chunk, which was a syntax error that would stop the script.
##  - The second chunk previously looped over ALL FOUR measures
##    (Shannon, Tsallis, Renyi, Fisher), producing 8 PDFs. You asked for a
##    SINGLE plot for Shannon only, so that loop is now restricted to
##    "Shannon" -- producing exactly 2 PDFs (noCI / withCI).
##  - Removed leftover Displacement branch (var_label else-clause) since
##    input_files now only contains rsam_avg -- that branch could never
##    fire, but is now removed to avoid confusion.
##  - Renamed the single-plot output folder from "HC_Plane_Plots_single_2019"
##    to "HC_Plane_Plots_single_Shannon_2007to2019" to reflect the actual
##    date range and content.
##  - Updated stale comment ("3 most extreme points") to match N_LABEL <- 4.
##
## UPDATE 1 (kept from your version): Only weekly (overlapping) windows
## whose Week_Start falls on or before `cutoff_date` (9 Dec 2019, the
## Whakaari eruption date) are included. Any week starting after that date
## is excluded entirely. NOTE: with the full 2007-2019 dataset this only
## affects late-2019 weeks -- 2007-2018 weeks are unaffected. Confirm this
## is still what you want; remove the filter line in build_long_df() if not.
##
## UPDATE 2: The N_LABEL (=4) most extreme points in each facet/plot are
## labelled with their Week_Index number ("week N"), using ggrepel so
## labels don't overlap the points. For Shannon/Tsallis/Renyi this is the
## lowest H_norm - C_norm score (low-H / high-C corner). For Fisher this is
## REVERSED -- Fisher behaves oppositely to the other three measures, so
## its label corner is the highest H_norm + C_norm score instead.
##
## UPDATE 3: min.segment.length = 0 and a black ring marker on every
## labeled point, so the leader line is never silently omitted and it's
## always unambiguous which dot a label belongs to.
##
## UPDATE 4: Point colour is Week_Index (continuous viridis scale), showing
## temporal evolution rather than duplicating the H position on the x-axis.
## ============================================================================
##
## ASSUMPTIONS (edit CONFIG below if wrong):
##  1. D (embedding dimension) MUST match the D used to generate the weekly
##     CSVs, since it's used to pick the correct theoretical Shannon H-C
##     boundary curve from StatOrdPattHxC::LinfLsup.
##  2. CI layers: the weekly CSV only contains a variance/CI for H in all
##     four measures (Semi_H_*), but a variance/CI for C only for Shannon
##     (Semi_C_Shannon). So "with CI" plots show horizontal error bars (H
##     direction) on every measure, and vertical error bars (C direction)
##     on the Shannon plot only.
##  3. Week_Start / Week_End in the CSV are already correct UTC dates/times
##     -- cutoff filtering compares against them directly, forcing
##     tz = "UTC" whenever coerced to Date/POSIXct.
##  4. "Min entropy / max complexity" scoring: C rescaled to [0,1] per
##     measure (C_norm), score = H_norm - C_norm. The N_LABEL lowest-score
##     points are labelled (low-H/high-C corner) -- EXCEPT Fisher, which
##     uses score_fisher = H_norm + C_norm (largest values), since its
##     extreme corner is high-H/high-C, not low-H/high-C.
##  5. ONLY rsam_avg is analysed -- no Displacement file/column anywhere.
## ============================================================================

library(readr)
library(tidyverse)
library(lubridate)
library(here)
library(StatOrdPattHxC)
library(scales)
library(ggrepel)   # install.packages("ggrepel") if not already installed
library(writexl)   # install.packages("writexl") if not already installed -- for Excel export
library(dbscan)    # install.packages("dbscan") if not already installed -- cluster/noise detection
library(FNN)        # install.packages("FNN") if not already installed -- k-NN isolation distance

## ---- CONFIG ---------------------------------------------------------------
D <- 5   # <-- MUST match the embedding dimension used to generate the CSVs

# Only keep weekly (overlapping) windows starting on/before this date.
cutoff_date <- as.Date("2019-12-09", tz = "UTC")

# Number of extreme points to label per facet/plot (you asked for 10-15;
# 12 is a reasonable middle value -- change freely, e.g. to 10 or 15)
N_LABEL <- 12

## ---- Cluster-separation (outlier) detection settings ----
## DBSCAN flags points that don't belong to the main dense H-C trend as
## "noise" (cluster = 0), regardless of which direction they're separated
## in -- unlike the corner-based N_LABEL score above, this also catches
## points like low-H/low-C that aren't near either scored corner.
## minPts: a point needs at least this many neighbours within `eps` to
## count as part of a dense region. eps: neighbourhood radius, in the same
## units as normalized H_norm/C_norm (i.e. a fraction of the full data
## range, since both axes are rescaled to [0,1] before clustering).
## THESE TWO VALUES ARE DATA-DEPENDENT -- see the kNN-distance plot the
## script prints before clustering, and adjust eps to sit at the "elbow"
## of that plot if the default doesn't separate things the way you expect.
DBSCAN_MINPTS <- 5
DBSCAN_EPS    <- 0.03
K_NEIGHBORS   <- 5   # for the continuous kNN isolation-distance diagnostic

results_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results/Weekly_OP_Analysis_2008to2019"

out_dir_facet   <- file.path(results_dir, "HC_Plane_Plots_2008to2019")
out_dir_shannon <- file.path(results_dir, "HC_Plane_Plots_single_Shannon_2008to2019")
dir.create(out_dir_facet,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_shannon, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  rsam_avg = file.path(results_dir, "weekly_entropy_complexity_rsam_avg_2008to2019.csv")
)

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
    )
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

## ---- Reshape one variable's weekly CSV into long H-C-CI format -----------
build_long_df <- function(path) {
  
  df <- read_csv(path, show_col_types = FALSE)
  
  # Force Week_Start / Week_End to UTC (no silent local-tz shift).
  df <- df %>%
    mutate(
      Week_Start = as.Date(Week_Start, tz = "UTC"),
      Week_End   = as.Date(Week_End,   tz = "UTC")
    )
  
  # ---- Cutoff: drop any overlapping week starting AFTER 9 Dec 2019 ----
  df <- df %>% filter(Week_Start <= cutoff_date)
  
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
      score        = H_norm - C_norm,   # low score = low H/high C corner
      score_fisher = H_norm + C_norm    # high score = high H/high C corner (Fisher only)
    ) %>%
    ungroup()
  
  df_long
}

## ---- Pick the extreme points per facet to label ---------------------------
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

## ---- Pick the extreme points for ONE measure to label ---------------------
get_label_points_one <- function(df_measure, measure_name, n = N_LABEL) {
  if (measure_name == "Fisher") {
    df_measure %>% slice_max(score_fisher, n = n, with_ties = FALSE)
  } else {
    df_measure %>% slice_min(score, n = n, with_ties = FALSE)
  }
}

## ---- Detect points separated from the main H-C cluster (any direction) ---
## Works per Measure (since each measure has its own H_norm/C_norm scale).
## Adds:
##   knn_dist  -- mean distance to the K_NEIGHBORS nearest points in
##                (H_norm, C_norm) space. Large = isolated/separated point,
##                regardless of which corner/direction it sits in.
##   dbscan_cluster -- DBSCAN cluster id; 0 = "noise", i.e. not part of any
##                dense region -- this is the direct "separated from the
##                cluster" flag.
##   is_separated -- convenience boolean, TRUE when dbscan_cluster == 0.
add_separation_flags <- function(df_measure) {
  
  xy <- as.matrix(df_measure[, c("H_norm", "C_norm")])
  
  # kNN isolation distance (continuous diagnostic)
  k_use <- min(K_NEIGHBORS, nrow(xy) - 1)
  knn_res <- FNN::get.knn(xy, k = k_use)
  df_measure$knn_dist <- rowMeans(knn_res$nn.dist)
  
  # DBSCAN noise/cluster labeling (categorical "separated or not")
  db <- dbscan::dbscan(xy, eps = DBSCAN_EPS, minPts = DBSCAN_MINPTS)
  df_measure$dbscan_cluster <- db$cluster
  df_measure$is_separated   <- db$cluster == 0L
  
  df_measure
}

## ---- Diagnostic: kNN-distance plot to help choose DBSCAN_EPS -------------
## Sorted kth-nearest-neighbour distances -- look for the "elbow" (the
## point where the curve bends sharply upward) and set DBSCAN_EPS to that
## y-value. Run this once per measure/variable if the default eps doesn't
## separate points the way the plot suggests it should.
print_knn_distance_plot <- function(df_measure, measure_name) {
  xy <- as.matrix(df_measure[, c("H_norm", "C_norm")])
  k_use <- min(DBSCAN_MINPTS, nrow(xy) - 1)
  d_sorted <- sort(FNN::get.knn(xy, k = k_use)$nn.dist[, k_use])
  plot(d_sorted, type = "l",
       main = paste("kNN distance plot --", measure_name, "(pick eps at the elbow)"),
       xlab = "Points sorted by distance", ylab = paste0(DBSCAN_MINPTS, "-NN distance"))
  abline(h = DBSCAN_EPS, col = "red", lty = 2)
}



make_hc_plot_facet <- function(df_long, bounds_cropped, label_points, var_label, with_ci) {
  
  p <- ggplot(df_long, aes(x = H, y = C, color = Week_Index)) +
    geom_point(size = 1.6, alpha = 0.8) +
    scale_color_viridis_c(name = "Week", option = "D") +
    facet_wrap(vars(Measure), ncol = 2, scales = "free") +
    labs(
      title = bquote(.(var_label) ~ italic(H) %*% italic(C) ~ "plane"),
      x = expression(italic(H)),
      y = expression(italic(C))
    ) +
    hc_theme
  
  if (with_ci) {
    p <- p +
      geom_segment(
        data = df_long %>% filter(is.finite(Semi_H)),
        aes(x = pmax(0, H - Semi_H), xend = H + Semi_H, y = C, yend = C),
        linewidth = 0.3, alpha = 0.6
      ) +
      geom_segment(
        data = df_long %>% filter(Measure == "Shannon", is.finite(Semi_C)),
        aes(x = H, xend = H, y = pmax(0, C - Semi_C), yend = C + Semi_C),
        linewidth = 0.3, alpha = 0.6
      )
  }
  
  p <- add_shannon_boundary(
    p,
    bounds_cropped %>% mutate(Measure = factor("Shannon", levels = measure_levels))
  )
  
  # ---- Mark points DBSCAN flags as separated from the main cluster ----
  # Distinct marker (red X, no ring) from the corner-extreme points below.
  p <- p +
    geom_point(
      data = df_long %>% filter(is_separated),
      aes(x = H, y = C),
      inherit.aes = FALSE,
      shape = 4, size = 2.2, stroke = 0.9,
      color = "red"
    )
  
  # ---- Combine BOTH label sets into one repel layer -----------------------
  # Extreme (corner-score) points and DBSCAN-separated points are labelled
  # together in a single geom_text_repel call so ggrepel can push all the
  # labels apart from EACH OTHER, not just within their own set -- two
  # independent calls would let an extreme-label and a separated-label
  # collide since repulsion only works within one layer at a time.
  # color = I(...) uses the literal colour string directly (not mapped
  # through a scale), so it coexists with the continuous Week_Index colour
  # scale already used for the points.
  separated_not_extreme <- df_long %>%
    filter(is_separated) %>%
    anti_join(label_points, by = c("Week_Index", "Measure"))
  
  combined_labels <- bind_rows(
    label_points          %>% mutate(label_color = "black"),
    separated_not_extreme %>% mutate(label_color = "red")
  )
  
  p <- p +
    geom_point(
      data = label_points,
      aes(x = H, y = C),
      inherit.aes = FALSE,
      shape = 21, size = 3, stroke = 1,
      color = "black", fill = NA
    ) +
    geom_text_repel(
      data = combined_labels,
      aes(x = H, y = C, label = paste0("week ", Week_Index), color = I(label_color)),
      inherit.aes = FALSE,
      size = 2.3,
      fontface = "bold",
      segment.size = 0.25,
      segment.color = "grey30",
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.2,
      force = 2,          # push labels apart harder with more of them present
      max.iter = 5000,    # more iterations to find a clean layout at higher label counts
      max.overlaps = Inf,
      seed = 42
    )
  
  p
}

## ============================================================================
## PART 2: STANDALONE Shannon-only plot
## ============================================================================

make_hc_plot_shannon <- function(df_measure, bounds_cropped, label_points, var_label, with_ci) {
  
  p <- ggplot(df_measure, aes(x = H, y = C, color = Week_Index)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_viridis_c(name = "Week", option = "D") +
    labs(
      title = bquote(.(var_label) ~ "Shannon" ~ italic(H) %*% italic(C) ~ "plane"),
      x = expression(italic(H)),
      y = expression(italic(C))
    ) +
    hc_theme
  
  if (with_ci) {
    p <- p +
      geom_segment(
        data = df_measure %>% filter(is.finite(Semi_H)),
        aes(x = pmax(0, H - Semi_H), xend = H + Semi_H, y = C, yend = C),
        linewidth = 0.3, alpha = 0.6
      ) +
      geom_segment(
        data = df_measure %>% filter(is.finite(Semi_C)),
        aes(x = H, xend = H, y = pmax(0, C - Semi_C), yend = C + Semi_C),
        linewidth = 0.3, alpha = 0.6
      )
  }
  
  # ---- Mark points DBSCAN flags as separated from the main cluster ----
  p <- add_shannon_boundary(p, bounds_cropped)
  
  p <- p +
    geom_point(
      data = df_measure %>% filter(is_separated),
      aes(x = H, y = C),
      inherit.aes = FALSE,
      shape = 4, size = 2.5, stroke = 1,
      color = "red"
    )
  
  # ---- Combine extreme + separated labels into one repel layer ----
  separated_not_extreme <- df_measure %>%
    filter(is_separated) %>%
    anti_join(label_points, by = "Week_Index")
  
  combined_labels <- bind_rows(
    label_points          %>% mutate(label_color = "black"),
    separated_not_extreme %>% mutate(label_color = "red")
  )
  
  p <- p +
    geom_point(
      data = label_points,
      aes(x = H, y = C),
      inherit.aes = FALSE,
      shape = 21, size = 3.5, stroke = 1,
      color = "black", fill = NA
    ) +
    geom_text_repel(
      data = combined_labels,
      aes(x = H, y = C, label = paste0("week ", Week_Index), color = I(label_color)),
      inherit.aes = FALSE,
      size = 2.6,
      fontface = "bold",
      segment.size = 0.25,
      segment.color = "grey30",
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.2,
      force = 2,
      max.iter = 5000,
      max.overlaps = Inf,
      seed = 42
    )
  
  p
}

## ---- Run for each variable (RSAM only): facet plots + Shannon-only plots --
for (var_name in names(input_files)) {
  
  message("Building H-C plane plots for: ", var_name)
  
  df_long <- build_long_df(input_files[[var_name]])
  
  if (nrow(df_long) == 0) {
    warning("No weeks remain after cutoff filter for ", var_name, " -- skipping.")
    next
  }
  
  bounds_cropped <- crop_bounds_to_data(
    bounds,
    df_long %>% filter(Measure == "Shannon")
  )
  
  ## ---- Detect points separated from the main cluster, per measure ----
  df_long <- df_long %>%
    group_by(Measure) %>%
    group_modify(~ add_separation_flags(.x)) %>%
    ungroup()
  
  # Optional: uncomment to inspect the kNN-distance elbow plot per measure
  # and sanity-check DBSCAN_EPS before trusting the flags below.
  # for (m in measure_levels) {
  #   print_knn_distance_plot(df_long %>% filter(Measure == m), m)
  # }
  
  separated_points <- df_long %>% filter(is_separated)
  message(sprintf("  -> %d points flagged as separated from the main cluster (all measures, %s)",
                  nrow(separated_points), var_name))
  
  var_label <- "RSAM"
  
  ## ---- PART 1: faceted 4-panel plots (noCI / withCI) ----
  label_points_facet <- get_label_points(df_long, n = N_LABEL)
  
  p_facet_noci <- make_hc_plot_facet(df_long, bounds_cropped, label_points_facet, var_label, with_ci = FALSE)
  print(p_facet_noci)
  ggsave(
    filename = file.path(out_dir_facet, sprintf("HC_Plane_%s_D%d_noCI_upto9Dec_byWeek.pdf", var_name, D)),
    plot = p_facet_noci, width = 26, height = 20, units = "cm", device = "pdf"
  )
  
  p_facet_ci <- make_hc_plot_facet(df_long, bounds_cropped, label_points_facet, var_label, with_ci = TRUE)
  print(p_facet_ci)
  ggsave(
    filename = file.path(out_dir_facet, sprintf("HC_Plane_%s_D%d_withCI_upto9Dec_byWeek.pdf", var_name, D)),
    plot = p_facet_ci, width = 26, height = 20, units = "cm", device = "pdf"
  )
  
  message("  -> saved facet noCI/withCI figures to ", out_dir_facet)
  
  ## ---- PART 2: standalone Shannon-only plot (noCI / withCI) ----
  df_shannon <- df_long %>% filter(Measure == "Shannon")
  label_points_shannon <- get_label_points_one(df_shannon, "Shannon", n = N_LABEL)
  
  p_shannon_noci <- make_hc_plot_shannon(df_shannon, bounds_cropped, label_points_shannon, var_label, with_ci = FALSE)
  print(p_shannon_noci)
  ggsave(
    filename = file.path(out_dir_shannon, sprintf("HC_Plane_%s_Shannon_D%d_noCI_upto9Dec_byWeek.pdf", var_name, D)),
    plot = p_shannon_noci, width = 18, height = 15, units = "cm", device = "pdf"
  )
  
  p_shannon_ci <- make_hc_plot_shannon(df_shannon, bounds_cropped, label_points_shannon, var_label, with_ci = TRUE)
  print(p_shannon_ci)
  ggsave(
    filename = file.path(out_dir_shannon, sprintf("HC_Plane_%s_Shannon_D%d_withCI_upto9Dec_byWeek.pdf", var_name, D)),
    plot = p_shannon_ci, width = 18, height = 15, units = "cm", device = "pdf"
  )
  
  message("  -> saved Shannon-only noCI/withCI figures to ", out_dir_shannon)
  
  ## ---- Export the labelled extreme points and separated points to Excel ----
  ## Three sheets:
  ##  - Facet_Labels: corner-extreme points from the facet plot (N_LABEL per measure)
  ##  - Shannon_Labels: corner-extreme points from the standalone Shannon plot
  ##  - Separated_Points: ALL points DBSCAN flagged as separated from the
  ##    main cluster, any measure, any direction -- this is the one to use
  ##    for checking against eruption records, since it isn't limited to
  ##    the two corners the extreme-point score targets.
  export_cols <- c("Week_Index", "Week_Start", "Week_End", "Measure",
                   "H", "C", "Semi_H", "Semi_C", "H_norm", "C_norm",
                   "score", "score_fisher", "knn_dist", "dbscan_cluster")
  
  excel_out <- list(
    Facet_Labels      = label_points_facet   %>% select(any_of(export_cols)) %>% arrange(Measure, Week_Index),
    Shannon_Labels    = label_points_shannon %>% select(any_of(export_cols)) %>% arrange(Week_Index),
    Separated_Points  = separated_points     %>% select(any_of(export_cols)) %>% arrange(Measure, desc(knn_dist))
  )
  
  excel_path <- file.path(results_dir, sprintf("HC_extreme_weeks_%s_D%d.xlsx", var_name, D))
  writexl::write_xlsx(excel_out, path = excel_path)
  message("  -> saved labelled extreme-week points to ", excel_path)
}

message("Done.")
