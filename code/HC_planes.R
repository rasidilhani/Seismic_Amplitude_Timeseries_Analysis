# HC plots for all — monthly, by variable, with confidence intervals
library(readxl)
library(tidyverse)
library(here)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────
D <- 5

# ── Read the combined long-format results ──────────────────────────────────
# NOTE: this assumes your R session's working directory / project root is
# C:\Users\UserA1\Documents\GitHub\Seismic_Amplitude_Timeseries_Analysis,
# so that here("results", ...) resolves to the "results" folder under it.
# If here() is pointing somewhere else, run here::here() on its own line to
# check, or call here::i_am("Seismic_Amplitude_Timeseries_Analysis.Rproj")
# (or any file that lives at the repo root) once at the top of the script.
input_path <- here("results", "WIZ_OrdinalPatterns_D5_Monthly_Quarterly.xlsx")

df_all <- read_xlsx(input_path, sheet = 1)

# Keep only the monthly rows — the quarterly rows (Q1..Q4) live in the same
# table via Period_Type and would otherwise show up as stray, miscolored
# points with no real calendar position.
df_all <- df_all %>%
  filter(Period_Type == "Monthly") %>%
  mutate(
    Period = factor(Period, levels = month.abb)   # "Jan".."Dec", in calendar order
  )

# ── Theoretical HC boundary for Shannon only ──────────────────────────────
data("LinfLsup")

bounds <- LinfLsup %>%
  filter(as.integer(as.character(Dimension)) == D)

bound_group_col <- "Side"

# ── Article-friendly colors, keyed to the actual Period values (Jan..Dec) ──
month_colors <- c(
  Jan = "#1B9E77",
  Feb = "#D95F02",
  Mar = "#7570B3",
  Apr = "#E7298A",
  May = "#66A61E",
  Jun = "#E6AB02",
  Jul = "#A6761D",
  Aug = "#666666",
  Sep = "#377EB8",
  Oct = "#E41A1C",
  Nov = "#4DAF4A",
  Dec = "#984EA3"
)

# ── Common plotting function ──────────────────────────────────────────────
# H_col/C_col: the entropy/complexity columns to plot
# Semi_H_col/Semi_C_col: their matching semi-length columns, used to draw a
#   95% CI cross (horizontal errorbar on H, vertical errorbar on C) around
#   each point. Pass NULL for either to skip that direction (Var_C_Shannon
#   is the only complexity-variance column in your sheet, so CI on C is only
#   available for the Shannon plot).
make_hc_plot <- function(df, H_col, C_col, Semi_H_col = NULL, Semi_C_col = NULL,
                         title_text, use_bounds = FALSE) {
  
  df <- df %>% filter(is.finite(.data[[H_col]]), is.finite(.data[[C_col]]))
  
  # Expand the axis ranges to leave room for the error bars, not just the points
  h_semi_max <- if (!is.null(Semi_H_col)) max(df[[Semi_H_col]], na.rm = TRUE) else 0
  c_semi_max <- if (!is.null(Semi_C_col)) max(df[[Semi_C_col]], na.rm = TRUE) else 0
  if (!is.finite(h_semi_max)) h_semi_max <- 0
  if (!is.finite(c_semi_max)) c_semi_max <- 0
  
  x_min <- min(df[[H_col]], na.rm = TRUE) - h_semi_max - 0.01
  x_max <- max(df[[H_col]], na.rm = TRUE) + h_semi_max + 0.01
  y_min <- max(0, min(df[[C_col]], na.rm = TRUE) - c_semi_max - 0.01)
  y_max <- max(df[[C_col]], na.rm = TRUE) + c_semi_max + 0.01
  
  p <- ggplot(df, aes(x = .data[[H_col]], y = .data[[C_col]], color = Period))
  
  # Vertical CI bar (uncertainty in C)
  if (!is.null(Semi_C_col)) {
    p <- p + geom_errorbar(
      aes(ymin = .data[[C_col]] - .data[[Semi_C_col]],
          ymax = .data[[C_col]] + .data[[Semi_C_col]]),
      width = 0, alpha = 0.5, linewidth = 0.4
    )
  }
  
  # Horizontal CI bar (uncertainty in H)
  if (!is.null(Semi_H_col)) {
    p <- p + geom_errorbarh(
      aes(xmin = .data[[H_col]] - .data[[Semi_H_col]],
          xmax = .data[[H_col]] + .data[[Semi_H_col]]),
      height = 0, alpha = 0.5, linewidth = 0.4
    )
  }
  
  p <- p +
    geom_point(size = 1, alpha = 0.85) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    scale_color_manual(values = month_colors, drop = FALSE) +
    labs(
      x     = expression(italic(H)),
      y     = expression(italic(C)),
      color = NULL,
      title = title_text
    ) +
    theme_bw(base_size = 11, base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.border     = element_rect(color = "black", linewidth = 0.5),
      legend.position  = "bottom",
      legend.text      = element_text(size = 7),
      axis.title       = element_text(size = 9),
      axis.text        = element_text(size = 8, color = "black"),
      plot.title       = element_text(size = 9, face = "bold", hjust = 0.5)
    ) +
    guides(
      color = guide_legend(
        nrow = 3,
        override.aes = list(size = 1, alpha = 0.8)
      )
    )
  
  if (use_bounds) {
    p <- p +
      geom_line(
        data        = bounds,
        aes(x = H, y = C, group = .data[[bound_group_col]]),
        color       = "grey35",
        linewidth   = 0.6,
        inherit.aes = FALSE
      )
  }
  
  return(p)
}

# ── Build and save one full set of 4 plots for a given variable ────────────
build_and_save_variable <- function(var_label, df_all, D) {
  
  df <- df_all %>% filter(Variable == var_label)
  
  p_shannon <- make_hc_plot(
    df,
    H_col = "H_Shannon", C_col = "C_Shannon",
    Semi_H_col = "Semi_H_Shannon", Semi_C_col = "Semi_C_Shannon",
    title_text = paste0("Shannon, D = ", D, " (", var_label, ")"),
    use_bounds = TRUE
  )
  
  p_renyi <- make_hc_plot(
    df,
    H_col = "H_Renyi", C_col = "C_Renyi",
    Semi_H_col = "Semi_H_Renyi", Semi_C_col = NULL,   # no Var_C_Renyi column available
    title_text = paste0("Renyi, D = ", D, " (", var_label, ")"),
    use_bounds = FALSE
  )
  
  p_tsallis <- make_hc_plot(
    df,
    H_col = "H_Tsallis", C_col = "C_Tsallis",
    Semi_H_col = "Semi_H_Tsallis", Semi_C_col = NULL, # no Var_C_Tsallis column available
    title_text = paste0("Tsallis, D = ", D, " (", var_label, ")"),
    use_bounds = FALSE
  )
  
  p_fisher <- make_hc_plot(
    df,
    H_col = "H_Fisher", C_col = "C_Fisher",
    Semi_H_col = "Semi_H_Fisher", Semi_C_col = NULL,  # no Var_C_Fisher column available
    title_text = paste0("Fisher, D = ", D, " (", var_label, ")"),
    use_bounds = FALSE
  )
  
  print(p_shannon)
  print(p_renyi)
  print(p_tsallis)
  print(p_fisher)
  
  out_dir <- here("results", "HC_Results_New")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(file.path(out_dir, paste0("HC_Plot_Shannon_", var_label, "_D5.pdf")),
         plot = p_shannon, width = 8, height = 7, units = "cm", device = "pdf")
  ggsave(file.path(out_dir, paste0("HC_Plot_Renyi_", var_label, "_D5.pdf")),
         plot = p_renyi, width = 8, height = 7, units = "cm", device = "pdf")
  ggsave(file.path(out_dir, paste0("HC_Plot_Tsallis_", var_label, "_D5.pdf")),
         plot = p_tsallis, width = 8, height = 7, units = "cm", device = "pdf")
  ggsave(file.path(out_dir, paste0("HC_Plot_Fisher_", var_label, "_D5.pdf")),
         plot = p_fisher, width = 8, height = 7, units = "cm", device = "pdf")
  
  cat("Plots saved for", var_label, "to", out_dir, "\n")
}

# ── Run for both variables ─────────────────────────────────────────────────
for (v in unique(df_all$Variable)) {
  build_and_save_variable(v, df_all, D)
}

# End of the code

###################################################################################
# Without confidence interval
# HC plots for all — monthly, by variable, no confidence intervals
library(readxl)
library(tidyverse)
library(here)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────
D <- 5

# ── Read the combined long-format results ──────────────────────────────────
# NOTE: this assumes your R session's working directory / project root is
# C:\Users\UserA1\Documents\GitHub\Seismic_Amplitude_Timeseries_Analysis,
# so that here("results", ...) resolves to the "results" folder under it.
input_path <- here("results", "WIZ_OrdinalPatterns_D5_Monthly_Quarterly.xlsx")

df_all <- read_xlsx(input_path, sheet = 1)

# Keep only the monthly rows — the quarterly rows (Q1..Q4) live in the same
# table via Period_Type and would otherwise show up as stray, miscolored
# points with no real calendar position.
df_all <- df_all %>%
  filter(Period_Type == "Monthly") %>%
  mutate(
    Period = factor(Period, levels = month.abb)   # "Jan".."Dec", in calendar order
  )

# ── Theoretical HC boundary for Shannon only ──────────────────────────────
data("LinfLsup")

bounds <- LinfLsup %>%
  filter(as.integer(as.character(Dimension)) == D)

bound_group_col <- "Side"

# ── Article-friendly colors, keyed to the actual Period values (Jan..Dec) ──
month_colors <- c(
  Jan = "#1B9E77",
  Feb = "#D95F02",
  Mar = "#7570B3",
  Apr = "#E7298A",
  May = "#66A61E",
  Jun = "#E6AB02",
  Jul = "#A6761D",
  Aug = "#666666",
  Sep = "#377EB8",
  Oct = "#E41A1C",
  Nov = "#4DAF4A",
  Dec = "#984EA3"
)

# ── Common plotting function ──────────────────────────────────────────────
make_hc_plot <- function(df, H_col, C_col, title_text, use_bounds = FALSE) {
  
  df <- df %>% filter(is.finite(.data[[H_col]]), is.finite(.data[[C_col]]))
  
  x_min <- min(df[[H_col]], na.rm = TRUE) - 0.01
  x_max <- max(df[[H_col]], na.rm = TRUE) + 0.01
  y_min <- max(0, min(df[[C_col]], na.rm = TRUE) - 0.01)
  y_max <- max(df[[C_col]], na.rm = TRUE) + 0.01
  
  p <- ggplot(df, aes(x = .data[[H_col]], y = .data[[C_col]], color = Period)) +
    geom_point(size = 1, alpha = 0.65) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    scale_color_manual(values = month_colors, drop = FALSE) +
    labs(
      x     = expression(italic(H)),
      y     = expression(italic(C)),
      color = NULL,
      title = title_text
    ) +
    theme_bw(base_size = 11, base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.border     = element_rect(color = "black", linewidth = 0.5),
      legend.position  = "bottom",
      legend.text      = element_text(size = 7),
      axis.title       = element_text(size = 9),
      axis.text        = element_text(size = 8, color = "black"),
      plot.title       = element_text(size = 9, face = "bold", hjust = 0.5)
    ) +
    guides(
      color = guide_legend(
        nrow = 2,
        override.aes = list(size = 1, alpha = 0.6)
      )
    )
  
  if (use_bounds) {
    p <- p +
      geom_line(
        data        = bounds,
        aes(x = H, y = C, group = .data[[bound_group_col]]),
        color       = "grey35",
        linewidth   = 0.6,
        inherit.aes = FALSE
      )
  }
  
  return(p)
}

# ── Build and save one full set of 4 plots for a given variable ────────────
build_and_save_variable <- function(var_label, df_all, D) {
  
  df <- df_all %>% filter(Variable == var_label)
  
  p_shannon <- make_hc_plot(
    df,
    H_col = "H_Shannon", C_col = "C_Shannon",
    title_text = paste0("Shannon, D = ", D, " (", var_label, ")"),
    use_bounds = TRUE
  )
  
  p_renyi <- make_hc_plot(
    df,
    H_col = "H_Renyi", C_col = "C_Renyi",
    title_text = paste0("Renyi, D = ", D, " (", var_label, ")"),
    use_bounds = FALSE
  )
  
  p_tsallis <- make_hc_plot(
    df,
    H_col = "H_Tsallis", C_col = "C_Tsallis",
    title_text = paste0("Tsallis, D = ", D, " (", var_label, ")"),
    use_bounds = FALSE
  )
  
  p_fisher <- make_hc_plot(
    df,
    H_col = "H_Fisher", C_col = "C_Fisher",
    title_text = paste0("Fisher, D = ", D, " (", var_label, ")"),
    use_bounds = FALSE
  )
  
  print(p_shannon)
  print(p_renyi)
  print(p_tsallis)
  print(p_fisher)
  
  out_dir <- here("results", "HC_Results_New")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(file.path(out_dir, paste0("HC_Plot_Shannon_", var_label, "_D5.pdf")),
         plot = p_shannon, width = 8, height = 7, units = "cm", device = "pdf")
  ggsave(file.path(out_dir, paste0("HC_Plot_Renyi_", var_label, "_D5.pdf")),
         plot = p_renyi, width = 8, height = 7, units = "cm", device = "pdf")
  ggsave(file.path(out_dir, paste0("HC_Plot_Tsallis_", var_label, "_D5.pdf")),
         plot = p_tsallis, width = 8, height = 7, units = "cm", device = "pdf")
  ggsave(file.path(out_dir, paste0("HC_Plot_Fisher_", var_label, "_D5.pdf")),
         plot = p_fisher, width = 8, height = 7, units = "cm", device = "pdf")
  
  cat("Plots saved for", var_label, "to", out_dir, "\n")
}

# ── Run for both variables ─────────────────────────────────────────────────
for (v in unique(df_all$Variable)) {
  build_and_save_variable(v, df_all, D)
}

# End of the code

###############################################################################
#With CI but separate plots for each month
# HC plots — one plot per month, per variable, per metric, with CIs
# Each point in a given plot is one year (2011-2022); colors distinguish years.
library(readxl)
library(tidyverse)
library(here)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────
D <- 5

# ── Read the combined long-format results ──────────────────────────────────
# Assumes the R session's working directory / project root is
# C:\Users\UserA1\Documents\GitHub\Seismic_Amplitude_Timeseries_Analysis,
# so here("results", ...) resolves to the "results" folder under it.
input_path <- here("results", "WIZ_OrdinalPatterns_D5_Monthly_Quarterly.xlsx")

df_all <- read_xlsx(input_path, sheet = 1)

df_all <- df_all %>%
  filter(Period_Type == "Monthly") %>%
  mutate(
    Period = factor(Period, levels = month.abb),  # "Jan".."Dec"
    Year   = factor(Year, levels = sort(unique(Year)))
  )

# ── Theoretical HC boundary for Shannon only ──────────────────────────────
data("LinfLsup")

bounds <- LinfLsup %>%
  filter(as.integer(as.character(Dimension)) == D)

bound_group_col <- "Side"

# ── Article-friendly colors, one per year (2011-2022) ──────────────────────
year_colors <- c(
  "2011" = "#1B9E77",
  "2012" = "#D95F02",
  "2013" = "#7570B3",
  "2014" = "#E7298A",
  "2015" = "#66A61E",
  "2016" = "#E6AB02",
  "2017" = "#A6761D",
  "2018" = "#666666",
  "2019" = "#377EB8",
  "2020" = "#E41A1C",
  "2021" = "#4DAF4A",
  "2022" = "#984EA3"
)

# ── Common plotting function (one month, one metric, one variable) ────────
# Each point is one year for that month; CI cross drawn from Semi_H_*/Semi_C_*.
make_hc_plot <- function(df, H_col, C_col, Semi_H_col = NULL, Semi_C_col = NULL,
                         title_text, use_bounds = FALSE) {
  
  df <- df %>% filter(is.finite(.data[[H_col]]), is.finite(.data[[C_col]]))
  
  h_semi_max <- if (!is.null(Semi_H_col)) max(df[[Semi_H_col]], na.rm = TRUE) else 0
  c_semi_max <- if (!is.null(Semi_C_col)) max(df[[Semi_C_col]], na.rm = TRUE) else 0
  if (!is.finite(h_semi_max)) h_semi_max <- 0
  if (!is.finite(c_semi_max)) c_semi_max <- 0
  
  x_min <- min(df[[H_col]], na.rm = TRUE) - h_semi_max - 0.01
  x_max <- max(df[[H_col]], na.rm = TRUE) + h_semi_max + 0.01
  y_min <- max(0, min(df[[C_col]], na.rm = TRUE) - c_semi_max - 0.01)
  y_max <- max(df[[C_col]], na.rm = TRUE) + c_semi_max + 0.01
  
  p <- ggplot(df, aes(x = .data[[H_col]], y = .data[[C_col]], color = Year))
  
  if (!is.null(Semi_C_col)) {
    p <- p + geom_errorbar(
      aes(ymin = .data[[C_col]] - .data[[Semi_C_col]],
          ymax = .data[[C_col]] + .data[[Semi_C_col]]),
      width = 0, alpha = 0.5, linewidth = 0.4
    )
  }
  
  if (!is.null(Semi_H_col)) {
    p <- p + geom_errorbarh(
      aes(xmin = .data[[H_col]] - .data[[Semi_H_col]],
          xmax = .data[[H_col]] + .data[[Semi_H_col]]),
      height = 0, alpha = 0.5, linewidth = 0.4
    )
  }
  
  p <- p +
    geom_point(size = 1.2, alpha = 0.85) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    scale_color_manual(values = year_colors, drop = FALSE) +
    labs(
      x     = expression(italic(H)),
      y     = expression(italic(C)),
      color = NULL,
      title = title_text
    ) +
    theme_bw(base_size = 11, base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.border     = element_rect(color = "black", linewidth = 0.5),
      legend.position  = "bottom",
      legend.text      = element_text(size = 7),
      axis.title       = element_text(size = 9),
      axis.text        = element_text(size = 8, color = "black"),
      plot.title       = element_text(size = 9, face = "bold", hjust = 0.5)
    ) +
    guides(
      color = guide_legend(
        nrow = 2,
        override.aes = list(size = 1, alpha = 0.8)
      )
    )
  
  if (use_bounds) {
    p <- p +
      geom_line(
        data        = bounds,
        aes(x = H, y = C, group = .data[[bound_group_col]]),
        color       = "grey35",
        linewidth   = 0.6,
        inherit.aes = FALSE
      )
  }
  
  return(p)
}

# ── Build and save the 4 metric plots for one (Variable, Month) pair ──────
build_and_save_month <- function(var_label, month_label, df_all, D) {
  
  df <- df_all %>% filter(Variable == var_label, Period == month_label)
  
  if (nrow(df) == 0) {
    cat("Skipping", var_label, month_label, "- no rows.\n")
    return(invisible(NULL))
  }
  
  p_shannon <- make_hc_plot(
    df,
    H_col = "H_Shannon", C_col = "C_Shannon",
    Semi_H_col = "Semi_H_Shannon", Semi_C_col = "Semi_C_Shannon",
    title_text = paste0("Shannon, D = ", D, " (", var_label, ", ", month_label, ")"),
    use_bounds = TRUE
  )
  
  p_renyi <- make_hc_plot(
    df,
    H_col = "H_Renyi", C_col = "C_Renyi",
    Semi_H_col = "Semi_H_Renyi", Semi_C_col = NULL,   # no Var_C_Renyi column available
    title_text = paste0("Renyi, D = ", D, " (", var_label, ", ", month_label, ")"),
    use_bounds = FALSE
  )
  
  p_tsallis <- make_hc_plot(
    df,
    H_col = "H_Tsallis", C_col = "C_Tsallis",
    Semi_H_col = "Semi_H_Tsallis", Semi_C_col = NULL, # no Var_C_Tsallis column available
    title_text = paste0("Tsallis, D = ", D, " (", var_label, ", ", month_label, ")"),
    use_bounds = FALSE
  )
  
  p_fisher <- make_hc_plot(
    df,
    H_col = "H_Fisher", C_col = "C_Fisher",
    Semi_H_col = "Semi_H_Fisher", Semi_C_col = NULL,  # no Var_C_Fisher column available
    title_text = paste0("Fisher, D = ", D, " (", var_label, ", ", month_label, ")"),
    use_bounds = FALSE
  )
  
  out_dir <- here("results", "HC_Results_New", "By_Month")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(file.path(out_dir, paste0("HC_Plot_Shannon_", month_label, "_", var_label, "_D5.pdf")),
         plot = p_shannon, width = 8, height = 7, units = "cm", device = "pdf")
  ggsave(file.path(out_dir, paste0("HC_Plot_Renyi_", month_label, "_", var_label, "_D5.pdf")),
         plot = p_renyi, width = 8, height = 7, units = "cm", device = "pdf")
  ggsave(file.path(out_dir, paste0("HC_Plot_Tsallis_", month_label, "_", var_label, "_D5.pdf")),
         plot = p_tsallis, width = 8, height = 7, units = "cm", device = "pdf")
  ggsave(file.path(out_dir, paste0("HC_Plot_Fisher_", month_label, "_", var_label, "_D5.pdf")),
         plot = p_fisher, width = 8, height = 7, units = "cm", device = "pdf")
  
  cat("Plots saved for", var_label, month_label, "to", out_dir, "\n")
}

# ── Run for every variable x every month ───────────────────────────────────
for (v in unique(df_all$Variable)) {
  for (m in month.abb) {
    build_and_save_month(v, m, df_all, D)
  }
}

# End of the code
################################################################
library(readxl)
library(tidyverse)
library(here)
library(StatOrdPattHxC)

D <- 5

input_path <- here("results", "WIZ_OrdinalPatterns_D5_Monthly_Quarterly.xlsx")

df_all <- read_xlsx(input_path, sheet = 1) %>%
  filter(Period_Type == "Monthly") %>%
  mutate(
    Period = factor(Period, levels = month.abb),
    Year = factor(Year, levels = 2011:2022)
  )

data("LinfLsup")

bounds_half <- bounds %>%
  group_by(H) %>%
  slice_max(C, n = 1, with_ties = FALSE) %>%
  ungroup()

year_colors <- c(
  "2011" = "#1B9E77",
  "2012" = "#D95F02",
  "2013" = "#7570B3",
  "2014" = "#E7298A",
  "2015" = "#66A61E",
  "2016" = "#E6AB02",
  "2017" = "#A6761D",
  "2018" = "#666666",
  "2019" = "#377EB8",
  "2020" = "#E41A1C",
  "2021" = "#4DAF4A",
  "2022" = "#984EA3"
)

make_hc_month_panel <- function(df, H_col, C_col,
                                Semi_H_col = NULL,
                                Semi_C_col = NULL,
                                title_text,
                                use_bounds = FALSE) {
  
  df <- df %>%
    filter(is.finite(.data[[H_col]]),
           is.finite(.data[[C_col]]))
  
  p <- ggplot(df, aes(x = .data[[H_col]],
                      y = .data[[C_col]],
                      colour = Year))
  
  if (!is.null(Semi_C_col)) {
    p <- p +
      geom_errorbar(
        aes(ymin = .data[[C_col]] - .data[[Semi_C_col]],
            ymax = .data[[C_col]] + .data[[Semi_C_col]]),
        width = 0,
        alpha = 0.45,
        linewidth = 0.25
      )
  }
  
  if (!is.null(Semi_H_col)) {
    p <- p +
      geom_errorbarh(
        aes(xmin = .data[[H_col]] - .data[[Semi_H_col]],
            xmax = .data[[H_col]] + .data[[Semi_H_col]]),
        height = 0,
        alpha = 0.45,
        linewidth = 0.25
      )
  }
  
  p <- p +
    geom_point(size = 1.3, alpha = 0.9) +
    facet_wrap(~ Period, ncol = 4) +
    scale_colour_manual(values = year_colors, drop = FALSE) +
    labs(
      title = title_text,
      x = expression(italic(H)),
      y = expression(italic(C)),
      colour = "Year"
    ) +
    theme_bw(base_family = "serif", base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 7, colour = "black"),
      axis.title = element_text(size = 10),
      panel.grid.minor = element_blank()
    ) +
    guides(colour = guide_legend(nrow = 2))
  
  if (use_bounds) {
    p <- p +
      geom_line(
        data = bounds_half,
        aes(x = H, y = C, group = Side),
        colour = "grey35",
        linewidth = 0.4,
        inherit.aes = FALSE
      )
  }
  
  return(p)
}

save_month_panel_plots <- function(var_label) {
  
  df_var <- df_all %>%
    filter(Variable == var_label)
  
  out_dir <- here("results", "HC_Results_New", "Monthly_Panel")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  p_shannon <- make_hc_month_panel(
    df_var,
    H_col = "H_Shannon",
    C_col = "C_Shannon",
    Semi_H_col = "Semi_H_Shannon",
    Semi_C_col = "Semi_C_Shannon",
    title_text = paste0("Monthly Shannon H-C Plane, D = ", D, " (", var_label, ")"),
    use_bounds = TRUE
  )
  
  p_renyi <- make_hc_month_panel(
    df_var,
    H_col = "H_Renyi",
    C_col = "C_Renyi",
    Semi_H_col = "Semi_H_Renyi",
    Semi_C_col = NULL,
    title_text = paste0("Monthly Rényi H-C Plane, D = ", D, " (", var_label, ")")
  )
  
  p_tsallis <- make_hc_month_panel(
    df_var,
    H_col = "H_Tsallis",
    C_col = "C_Tsallis",
    Semi_H_col = "Semi_H_Tsallis",
    Semi_C_col = NULL,
    title_text = paste0("Monthly Tsallis H-C Plane, D = ", D, " (", var_label, ")")
  )
  
  p_fisher <- make_hc_month_panel(
    df_var,
    H_col = "H_Fisher",
    C_col = "C_Fisher",
    Semi_H_col = "Semi_H_Fisher",
    Semi_C_col = NULL,
    title_text = paste0("Monthly Fisher H-C Plane, D = ", D, " (", var_label, ")")
  )
  
  ggsave(file.path(out_dir, paste0("Monthly_Panel_Shannon_", var_label, "_D5.pdf")),
         p_shannon, width = 18, height = 14, units = "cm")
  
  ggsave(file.path(out_dir, paste0("Monthly_Panel_Renyi_", var_label, "_D5.pdf")),
         p_renyi, width = 18, height = 14, units = "cm")
  
  ggsave(file.path(out_dir, paste0("Monthly_Panel_Tsallis_", var_label, "_D5.pdf")),
         p_tsallis, width = 18, height = 14, units = "cm")
  
  ggsave(file.path(out_dir, paste0("Monthly_Panel_Fisher_", var_label, "_D5.pdf")),
         p_fisher, width = 18, height = 14, units = "cm")
}

for (v in unique(df_all$Variable)) {
  save_month_panel_plots(v)
}


# ── Yearly HC plots: one plot per year, all months in one plot ─────────────

library(readxl)
library(tidyverse)
library(here)
library(StatOrdPattHxC)

D <- 5

input_path <- here("results", "WIZ_OrdinalPatterns_D5_Monthly_Quarterly.xlsx")

df_all <- read_xlsx(input_path, sheet = 1) %>%
  filter(Period_Type == "Monthly") %>%
  mutate(
    Period = factor(Period, levels = month.abb),
    Year   = factor(Year, levels = 2011:2022)
  )

data("LinfLsup")

bounds <- LinfLsup %>%
  filter(Dimension == as.character(D))

month_colors <- c(
  "Jan" = "#E41A1C",
  "Feb" = "#FF7F00",
  "Mar" = "#B3A214",
  "Apr" = "#66A61E",
  "May" = "#4DAF4A",
  "Jun" = "#1B9E77",
  "Jul" = "#00A6A6",
  "Aug" = "#377EB8",
  "Sep" = "#4C78A8",
  "Oct" = "#984EA3",
  "Nov" = "#F781BF",
  "Dec" = "#A65628"
)

make_hc_year_plot <- function(df, H_col, C_col,
                              Semi_H_col = NULL,
                              Semi_C_col = NULL,
                              title_text,
                              use_bounds = FALSE) {
  
  df <- df %>%
    filter(
      is.finite(.data[[H_col]]),
      is.finite(.data[[C_col]])
    )
  
  x_min <- min(df[[H_col]], na.rm = TRUE) - 0.01
  x_max <- max(df[[H_col]], na.rm = TRUE) + 0.01
  y_min <- max(0, min(df[[C_col]], na.rm = TRUE) - 0.01)
  y_max <- min(1, max(df[[C_col]], na.rm = TRUE) + 0.01)
  
  p <- ggplot(df, aes(x = .data[[H_col]],
                      y = .data[[C_col]],
                      colour = Period))
  
  if (use_bounds) {
    p <- p +
      geom_line(
        data = bounds,
        aes(x = H, y = C, group = Side),
        colour = "grey60",
        linewidth = 0.5,
        inherit.aes = FALSE
      )
  }
  
  if (!is.null(Semi_H_col)) {
    p <- p +
      geom_errorbarh(
        aes(xmin = .data[[H_col]] - .data[[Semi_H_col]],
            xmax = .data[[H_col]] + .data[[Semi_H_col]]),
        height = 0.003,
        linewidth = 0.4,
        alpha = 0.4
      )
  }
  
  if (!is.null(Semi_C_col)) {
    p <- p +
      geom_errorbar(
        aes(ymin = .data[[C_col]] - .data[[Semi_C_col]],
            ymax = .data[[C_col]] + .data[[Semi_C_col]]),
        width = 0.003,
        linewidth = 0.4,
        alpha = 0.7
      )
  }
  
  p +
    geom_point(
      shape = 21,
      aes(fill = Period),
      colour = "black",
      size = 3.5,
      stroke = 0.4
    ) +
    scale_colour_manual(values = month_colors, drop = FALSE) +
    scale_fill_manual(values = month_colors, drop = FALSE) +
    coord_cartesian(
      xlim = c(x_min, x_max),
      ylim = c(y_min, y_max)
    ) +
    labs(
      title = title_text,
      x = expression(italic(H)),
      y = expression(italic(C)),
      fill = "Month"
    ) +
    theme_bw(base_size = 11, base_family = "serif") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8, colour = "black"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(6, 6, 12, 6)
    ) +
    guides(
      colour = "none",
      fill = guide_legend(nrow = 3, byrow = TRUE)
    )
}
save_yearly_hc_plots <- function(var_label, year_label) {
  
  df_year <- df_all %>%
    filter(
      Variable == var_label,
      Year == year_label
    )
  
  out_dir <- here("results", "HC_Results_New", "Yearly_Plots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  p_shannon <- make_hc_year_plot(
    df_year,
    H_col = "H_Shannon",
    C_col = "C_Shannon",
    Semi_H_col = "Semi_H_Shannon",
    Semi_C_col = "Semi_C_Shannon",
    title_text = paste0("Shannon H-C Plane, D = ", D,
                        " (", var_label, ", ", year_label, ")"),
    use_bounds = TRUE
  )
  
  p_renyi <- make_hc_year_plot(
    df_year,
    H_col = "H_Renyi",
    C_col = "C_Renyi",
    Semi_H_col = "Semi_H_Renyi",
    Semi_C_col = NULL,
    title_text = paste0("Rényi H-C Plane, D = ", D,
                        " (", var_label, ", ", year_label, ")"),
    use_bounds = FALSE
  )
  
  p_tsallis <- make_hc_year_plot(
    df_year,
    H_col = "H_Tsallis",
    C_col = "C_Tsallis",
    Semi_H_col = "Semi_H_Tsallis",
    Semi_C_col = NULL,
    title_text = paste0("Tsallis H-C Plane, D = ", D,
                        " (", var_label, ", ", year_label, ")"),
    use_bounds = FALSE
  )
  
  p_fisher <- make_hc_year_plot(
    df_year,
    H_col = "H_Fisher",
    C_col = "C_Fisher",
    Semi_H_col = "Semi_H_Fisher",
    Semi_C_col = NULL,
    title_text = paste0("Fisher H-C Plane, D = ", D,
                        " (", var_label, ", ", year_label, ")"),
    use_bounds = FALSE
  )
  
  safe_var <- str_replace_all(var_label, "[^A-Za-z0-9_]", "_")
  
  ggsave(file.path(out_dir, paste0("Yearly_Shannon_", year_label, "_", safe_var, "_D5.pdf")),
         p_shannon, width = 15, height = 12, units = "cm")
  
  ggsave(file.path(out_dir, paste0("Yearly_Renyi_", year_label, "_", safe_var, "_D5.pdf")),
         p_renyi, width = 15, height = 12, units = "cm")
  
  ggsave(file.path(out_dir, paste0("Yearly_Tsallis_", year_label, "_", safe_var, "_D5.pdf")),
         p_tsallis, width = 15, height = 12, units = "cm")
  
  ggsave(file.path(out_dir, paste0("Yearly_Fisher_", year_label, "_", safe_var, "_D5.pdf")),
         p_fisher, width = 15, height = 12, units = "cm")
}

for (v in unique(df_all$Variable)) {
  for (y in levels(df_all$Year)) {
    save_yearly_hc_plots(v, y)
  }
}