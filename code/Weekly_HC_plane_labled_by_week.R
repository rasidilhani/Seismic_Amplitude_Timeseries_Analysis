## ============================================================================
## Weekly H-C (Entropy-Complexity) Plane Plots
## One figure per variable (rsam_avg, displacement_avg_m), each with 4 facets:
## Shannon, Tsallis, Renyi, Fisher. Produced with AND without confidence
## intervals.
##
## UPDATE 1: Only weekly (overlapping) windows whose Week_Start falls on or
## before `cutoff_date` are included -- i.e. the last window drawn is the
## final overlapping week that still touches 9 Dec 2019. Any week starting
## after that date is excluded entirely from the H-C planes.
##
## UPDATE 2: The 3 most extreme points in each facet are labelled with their
## Week_Index number ("week N"), using ggrepel so labels don't overlap the
## points. For Shannon/Tsallis/Renyi this is the lowest H_norm - C_norm
## score (low-H / high-C corner). For Fisher this is REVERSED -- Fisher
## behaves oppositely to the other three measures, so its label corner is
## the highest H_norm + C_norm score (high-H / high-C corner) instead.
##
## UPDATE 3 (bug fix): Added min.segment.length = 0 and a black ring marker
## on every labeled point. Previously, ggrepel could silently omit the
## leader line when it judged a label "close enough" to its point, and in
## crowded/compressed facets (Fisher's H_norm distribution can be skewed)
## that made a label look like it was pointing at the wrong neighbouring
## dot. The ring marker removes all ambiguity about which dot a label
## belongs to.
##
## UPDATE 4: Point colour changed from H (rescaled per measure) to Week_Index.
## Colouring by H was redundant with the x-axis (H IS the x-axis), so it
## only ever showed a left-to-right gradient that duplicated point position
## and added no new information. Colouring by week number instead shows
## temporal evolution -- e.g. whether early vs late weeks occupy different
## regions of the H-C plane, or whether the labelled extreme weeks cluster
## together in time. Uses a sequential viridis scale (appropriate for an
## ordered quantity like week number) instead of the old diverging
## red-yellow-blue scale (which was appropriate for "high vs low H" but
## not for a running index).
## ============================================================================
##
## ASSUMPTIONS (edit CONFIG below if wrong):
##  1. D (embedding dimension) here MUST match the D you used when you
##     generated the weekly CSVs, since it's only used to pick the correct
##     theoretical Shannon H-C boundary curve from StatOrdPattHxC::LinfLsup.
##  2. Point colour: now Week_Index, on a continuous viridis scale (see
##     UPDATE 4 above). H_norm / C_norm are still computed and used for
##     picking the extreme/labelled points -- they're just no longer used
##     for point colour.
##  3. CI layers: your weekly CSVs only contain a variance/CI for H in all
##     four measures (Semi_H_*), but a variance/CI for C only for Shannon
##     (Semi_C_Shannon). So the "with CI" plots show horizontal error bars
##     (H direction) on all four facets, and vertical error bars (C
##     direction) on the Shannon facet only -- there's simply no C variance
##     computed for Renyi/Tsallis/Fisher in the weekly script.
##  4. One combined figure per variable (all weeks together) -- not split by
##     month, since your weekly script doesn't produce a month field.
##  5. Week_Start / Week_End in the CSV are already correct UTC dates/times
##     (per your confirmation) -- so cutoff filtering below simply compares
##     against them directly, forcing tz = "UTC" whenever they are coerced
##     to Date/POSIXct, so no local-timezone shift is ever introduced here.
##  6. "Min entropy / max complexity" is judged the same way as in your
##     extreme-weeks Excel script: C is rescaled to [0,1] per measure too
##     (C_norm), and score = H_norm - C_norm. The 3 LOWEST-score points per
##     facet (per Measure) are labelled -- these are the low-H/high-C
##     corner points. Ties beyond 3 are broken arbitrarily by slice_min's
##     row order (with_ties = FALSE) so exactly 3 labels appear per facet.
##  7. EXCEPTION -- Fisher: its high-H corner coincides with high C, not low
##     C, so its label corner uses score_fisher = H_norm + C_norm (largest
##     values) instead of score (smallest values).
## ============================================================================

library(readr)
library(tidyverse)
library(lubridate)
library(here)
library(StatOrdPattHxC)
library(scales)
library(ggrepel)   # install.packages("ggrepel") if not already installed

## ---- CONFIG ---------------------------------------------------------------
D <- 5   # <-- MUST match the embedding dimension used to generate the CSVs

# Only keep weekly (overlapping) windows starting on/before this date.
# This keeps the last overlapping window that still touches 9 Dec 2019,
# and drops every week that starts after it.
cutoff_date <- as.Date("2019-12-09", tz = "UTC")

# Number of extreme points to label per facet
N_LABEL <- 4

results_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"
out_dir     <- file.path(results_dir, "HC_Plane_Plots_2019")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  rsam_avg           = file.path(results_dir, "weekly_entropy_complexity_rsam_avg_2019.csv"),
  displacement_avg_m = file.path(results_dir, "weekly_entropy_complexity_displacement_avg_m_2019.csv")
)

measure_levels <- c("Shannon", "Tsallis", "Renyi", "Fisher")

## ---- Theoretical Shannon H-C boundary (StatOrdPattHxC) --------------------
data("LinfLsup")

bounds <- LinfLsup %>%
  filter(as.integer(as.character(Dimension)) == D)

bound_group_col <- "Side"

# Crop the theoretical boundary curve to the H/C range actually spanned by
# the data (plus a small margin) -- shows only the relevant half of the
# upper/lower boundary curves near where the real data sits, instead of the
# full closed loop across the entire H domain.
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

## ---- Reshape one variable's weekly CSV into long H-C-CI format -----------
build_long_df <- function(path) {
  
  df <- read_csv(path, show_col_types = FALSE)
  
  # Force Week_Start / Week_End to be read as UTC (no silent local-tz shift),
  # whether they came in as character, Date, or POSIXct from read_csv.
  df <- df %>%
    mutate(
      Week_Start = as.Date(Week_Start, tz = "UTC"),
      Week_End   = as.Date(Week_End,   tz = "UTC")
    )
  
  # ---- Cutoff: drop any overlapping week starting AFTER 9 Dec 2019 ----
  # Keeps the last overlapping window whose Week_Start is still <= cutoff,
  # i.e. the final window that touches 9 Dec, and excludes everything
  # after it.
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
    # Safety net: never use a non-positive or non-finite variance/CI value,
    # even if one somehow made it into the CSV. A CI half-width must be > 0.
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

## ---- Plot builder ----------------------------------------------------------
## COLOUR CHANGE (UPDATE 4): aes(color = ...) is now Week_Index instead of
## H_norm, and the diverging red-yellow-blue scale is replaced with a
## continuous viridis scale, appropriate for a sequential/ordered quantity
## like week number rather than a "high vs low" contrast.
make_hc_plot <- function(df_long, bounds_cropped, label_points, var_label, with_ci) {
  
  p <- ggplot(df_long, aes(x = H, y = C, color = Week_Index)) +
    geom_point(size = 1.6, alpha = 0.8) +
    scale_color_viridis_c(
      name   = "Week",
      option = "D"
    ) +
    facet_wrap(vars(Measure), ncol = 2, scales = "free") +
    labs(
      title = bquote(.(var_label) ~ italic(H) %*% italic(C) ~
                       "plane"),
      x = expression(italic(H)),
      y = expression(italic(C))
    ) +
    hc_theme
  
  if (with_ci) {
    # Horizontal CI (H direction) -- available for all four measures.
    # H can never be negative, but H - Semi_H can dip below 0 when the
    # point estimate sits close to the boundary -- clamp the lower end.
    p <- p +
      geom_segment(
        data = df_long %>% filter(is.finite(Semi_H)),
        aes(x = pmax(0, H - Semi_H), xend = H + Semi_H, y = C, yend = C),
        linewidth = 0.3, alpha = 0.6
      )
    # Vertical CI (C direction) -- only available for Shannon.
    # Same reasoning: C - Semi_C can go negative even with a valid Semi_C.
    p <- p +
      geom_segment(
        data = df_long %>% filter(Measure == "Shannon", is.finite(Semi_C)),
        aes(x = H, xend = H, y = pmax(0, C - Semi_C), yend = C + Semi_C),
        linewidth = 0.3, alpha = 0.6
      )
  }
  
  p <- add_shannon_boundary(p, bounds_cropped)
  
  # ---- Ring the exact point being labeled -----------------------------
  # Makes the true target of each label unambiguous regardless of where
  # ggrepel nudges the text box -- without this, a crowded/compressed
  # facet (e.g. Fisher, which can have a skewed H_norm distribution) can
  # make it look like a label is pointing at the wrong neighbouring dot.
  p <- p +
    geom_point(
      data = label_points,
      aes(x = H, y = C),
      inherit.aes = FALSE,
      shape = 21, size = 3, stroke = 1,
      color = "black", fill = NA
    )
  
  # ---- Label the extreme points per facet -----------------------------
  # (Shannon/Tsallis/Renyi: low-H/high-C corner. Fisher: high-H/high-C
  #  corner -- see get_label_points().)
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
      min.segment.length = 0,   # always draw the leader line -- never
      # silently omit it, which is what made
      # labels look like they pointed at the
      # wrong dot
      box.padding = 0.4,
      point.padding = 0.3,
      max.overlaps = Inf,
      seed = 42
    )
  
  p
}

## ---- Run for each variable, save with and without CI ----------------------
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
  
  label_points <- get_label_points(df_long, n = N_LABEL)
  
  var_label <- if (var_name == "rsam_avg") "RSAM" else "Displacement"
  
  # -- without CI --
  p_noci <- make_hc_plot(df_long, bounds_cropped, label_points, var_label, with_ci = FALSE)
  print(p_noci)
  ggsave(
    filename = file.path(out_dir, sprintf("HC_Plane_%s_D%d_noCI_upto9Dec_byWeek.pdf", var_name, D)),
    plot = p_noci, width = 20, height = 16, units = "cm", device = "pdf"
  )
  
  # -- with CI --
  p_ci <- make_hc_plot(df_long, bounds_cropped, label_points, var_label, with_ci = TRUE)
  print(p_ci)
  ggsave(
    filename = file.path(out_dir, sprintf("HC_Plane_%s_D%d_withCI_upto9Dec_byWeek.pdf", var_name, D)),
    plot = p_ci, width = 20, height = 16, units = "cm", device = "pdf"
  )
  
  message("  -> saved noCI and withCI figures to ", out_dir)
}

message("Done.")

#####################################################################

## ============================================================================
## Weekly H-C (Entropy-Complexity) Plane Plots
## One figure PER MEASURE per variable (rsam_avg, displacement_avg_m):
## Shannon, Tsallis, Renyi, Fisher are now separate stand-alone plots rather
## than facets of one combined figure. Produced with AND without confidence
## intervals.
##
## [UPDATE 1-4 as before -- cutoff filter, extreme-point labelling, ggrepel
##  ring-marker fix, colour-by-week-number -- see prior version for details]
##
## UPDATE 5: Switched from one faceted figure (4 panels via facet_wrap) to
## FOUR SEPARATE single-panel figures per variable -- one per measure. This
## means:
##   - Each measure now gets its own PDF file (8 files per variable instead
##     of 2: 4 measures x {noCI, withCI}).
##   - The Shannon theoretical boundary curve is only drawn on the Shannon
##     plot (it was already Shannon-only in the faceted version, just less
##     obvious there since it only appeared in one of the four panels).
##   - Each plot now has its own title identifying both the variable AND
##     the measure (e.g. "RSAM Shannon H x C plane"), since there's no
##     shared facet strip to do that job anymore.
##   - Axis ranges are naturally independent per plot now (no more shared
##     facet_wrap(scales = "free") behaviour to think about -- each plot
##     is just its own ggplot object with its own default scale).
## ============================================================================
##
## ASSUMPTIONS (edit CONFIG below if wrong):
##  1. D (embedding dimension) here MUST match the D you used when you
##     generated the weekly CSVs, since it's only used to pick the correct
##     theoretical Shannon H-C boundary curve from StatOrdPattHxC::LinfLsup.
##  2. Point colour: Week_Index, on a continuous viridis scale -- shows
##     temporal evolution rather than duplicating the H position on the
##     x-axis.
##  3. CI layers: your weekly CSVs only contain a variance/CI for H in all
##     four measures (Semi_H_*), but a variance/CI for C only for Shannon
##     (Semi_C_Shannon). So the "with CI" plots show horizontal error bars
##     (H direction) on all four measures, and vertical error bars (C
##     direction) on the Shannon plot only.
##  4. One figure per variable PER MEASURE (see UPDATE 5) -- not split by
##     month, since your weekly script doesn't produce a month field.
##  5. Week_Start / Week_End in the CSV are already correct UTC dates/times
##     (per your confirmation) -- so cutoff filtering below simply compares
##     against them directly, forcing tz = "UTC" whenever they are coerced
##     to Date/POSIXct, so no local-timezone shift is ever introduced here.
##  6. "Min entropy / max complexity" is judged the same way as in your
##     extreme-weeks Excel script: C is rescaled to [0,1] per measure too
##     (C_norm), and score = H_norm - C_norm. The N_LABEL lowest-score
##     points per measure are labelled -- low-H/high-C corner points.
##  7. EXCEPTION -- Fisher: its high-H corner coincides with high C, not low
##     C, so its label corner uses score_fisher = H_norm + C_norm (largest
##     values) instead of score (smallest values).
## ============================================================================

library(readr)
library(tidyverse)
library(lubridate)
library(here)
library(StatOrdPattHxC)
library(scales)
library(ggrepel)   # install.packages("ggrepel") if not already installed

## ---- CONFIG ---------------------------------------------------------------
D <- 5   # <-- MUST match the embedding dimension used to generate the CSVs

# Only keep weekly (overlapping) windows starting on/before this date.
cutoff_date <- as.Date("2019-12-09", tz = "UTC")

# Number of extreme points to label per measure
N_LABEL <- 4

results_dir <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results"
out_dir     <- file.path(results_dir, "HC_Plane_Plots_single_2019")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  rsam_avg           = file.path(results_dir, "weekly_entropy_complexity_rsam_avg_2019.csv"),
  displacement_avg_m = file.path(results_dir, "weekly_entropy_complexity_displacement_avg_m_2019.csv")
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
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

## ---- Reshape one variable's weekly CSV into long H-C-CI format -----------
## (unchanged from before -- still builds all four measures into one long
##  data frame; the splitting into separate plots happens later, at the
##  plotting step, not here)
build_long_df <- function(path) {
  
  df <- read_csv(path, show_col_types = FALSE)
  
  df <- df %>%
    mutate(
      Week_Start = as.Date(Week_Start, tz = "UTC"),
      Week_End   = as.Date(Week_End,   tz = "UTC")
    )
  
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
      score        = H_norm - C_norm,
      score_fisher = H_norm + C_norm
    ) %>%
    ungroup()
  
  df_long
}

## ---- Pick the extreme points for ONE measure to label ---------------------
## (same logic as before, just applied to a single-measure slice instead of
##  grouped across all four at once)
get_label_points_one <- function(df_measure, measure_name, n = N_LABEL) {
  if (measure_name == "Fisher") {
    df_measure %>% slice_max(score_fisher, n = n, with_ties = FALSE)
  } else {
    df_measure %>% slice_min(score, n = n, with_ties = FALSE)
  }
}

## ---- Plot builder: ONE measure, ONE variable -------------------------------
## This replaces the old facet-based make_hc_plot(). It now takes a single
## measure's data (already filtered to one Measure) and returns one
## stand-alone ggplot object -- no facet_wrap involved.
make_hc_plot_single <- function(df_measure, bounds_cropped, label_points,
                                var_label, measure_name, with_ci) {
  
  p <- ggplot(df_measure, aes(x = H, y = C, color = Week_Index)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_viridis_c(
      name   = "Week",
      option = "D"
    ) +
    labs(
      title = bquote(.(var_label) ~ .(measure_name) ~ italic(H) %*% italic(C) ~
                       "plane"),
      x = expression(italic(H)),
      y = expression(italic(C))
    ) +
    hc_theme
  
  if (with_ci) {
    # Horizontal CI (H direction) -- available for every measure.
    p <- p +
      geom_segment(
        data = df_measure %>% filter(is.finite(Semi_H)),
        aes(x = pmax(0, H - Semi_H), xend = H + Semi_H, y = C, yend = C),
        linewidth = 0.3, alpha = 0.6
      )
    # Vertical CI (C direction) -- Shannon only (see ASSUMPTION 3 above).
    if (measure_name == "Shannon") {
      p <- p +
        geom_segment(
          data = df_measure %>% filter(is.finite(Semi_C)),
          aes(x = H, xend = H, y = pmax(0, C - Semi_C), yend = C + Semi_C),
          linewidth = 0.3, alpha = 0.6
        )
    }
  }
  
  # Theoretical boundary curve -- Shannon only, since LinfLsup is a Shannon
  # H-C boundary and doesn't apply to Tsallis/Renyi/Fisher's H-C geometry.
  if (measure_name == "Shannon") {
    p <- add_shannon_boundary(p, bounds_cropped)
  }
  
  # ---- Ring the exact point being labeled -----------------------------
  p <- p +
    geom_point(
      data = label_points,
      aes(x = H, y = C),
      inherit.aes = FALSE,
      shape = 21, size = 3.5, stroke = 1,
      color = "black", fill = NA
    )
  
  # ---- Label the extreme points -----------------------------
  p <- p +
    geom_text_repel(
      data = label_points,
      aes(x = H, y = C, label = paste0("week ", Week_Index)),
      inherit.aes = FALSE,
      color = "black",
      size = 3.2,
      fontface = "bold",
      segment.size = 0.3,
      segment.color = "grey30",
      min.segment.length = 0,
      box.padding = 0.4,
      point.padding = 0.3,
      max.overlaps = Inf,
      seed = 42
    )
  
  p
}

## ---- Run for each variable x each measure, save with and without CI ------
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
  
  var_label <- if (var_name == "rsam_avg") "RSAM" else "Displacement"
  
  # ---- Loop over each of the four measures, building a separate plot ----
  for (measure_name in measure_levels) {
    
    df_measure <- df_long %>% filter(Measure == measure_name)
    
    if (nrow(df_measure) == 0) {
      warning("No data for measure ", measure_name, " (", var_name, ") -- skipping.")
      next
    }
    
    label_points <- get_label_points_one(df_measure, measure_name, n = N_LABEL)
    
    # -- without CI --
    p_noci <- make_hc_plot_single(
      df_measure, bounds_cropped, label_points, var_label, measure_name, with_ci = FALSE
    )
    print(p_noci)
    ggsave(
      filename = file.path(
        out_dir,
        sprintf("HC_Plane_%s_%s_D%d_noCI_upto9Dec_byWeek.pdf", var_name, measure_name, D)
      ),
      plot = p_noci, width = 14, height = 12, units = "cm", device = "pdf"
    )
    
    # -- with CI --
    p_ci <- make_hc_plot_single(
      df_measure, bounds_cropped, label_points, var_label, measure_name, with_ci = TRUE
    )
    print(p_ci)
    ggsave(
      filename = file.path(
        out_dir,
        sprintf("HC_Plane_%s_%s_D%d_withCI_upto9Dec_byWeek.pdf", var_name, measure_name, D)
      ),
      plot = p_ci, width = 14, height = 12, units = "cm", device = "pdf"
    )
    
    message("  -> saved noCI and withCI figures for ", measure_name)
  }
}

message("Done.")