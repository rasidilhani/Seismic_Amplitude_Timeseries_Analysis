##########################################################
# HC plots by variable, combining Shannon, Tsallis, Renyi, and Fisher
# Legend labelled by Year
# Generates both No-CI and With-CI versions in a single pass
# Shannon panel zoomed to H = 0.85-1.00 and C = 0.00-0.20
#
# Works for BOTH Monthly and Quarterly data.
# To switch between them, change ONLY the line below (period_mode).
# Everything else (sheet, filter, factor levels, loop, labels,
# output folders, filenames) is derived automatically.
##########################################################

library(readxl)
library(tidyverse)
library(here)
library(StatOrdPattHxC)
library(patchwork)

# ---- Master switch ----

period_mode <- "Quarterly"   # change to "Monthly" to switch

# ---- Parameters ----

D <- 5

input_path <- here("results", "WIZ_OrdinalPatterns_D5_Monthly_Quarterly.xlsx")

# ---- Settings derived from period_mode ----
# (Nothing below this block needs manual editing when switching modes)

if (period_mode == "Monthly") {
  
  data_sheet     <- 1
  period_levels  <- month.abb
  period_values  <- month.abb
  period_label   <- "Month"
  results_suffix <- "Monthly_By_Year"
  
} else if (period_mode == "Quarterly") {
  
  data_sheet     <- 2
  period_levels  <- c("Q1", "Q2", "Q3", "Q4")
  period_values  <- c("Q1", "Q2", "Q3", "Q4")
  period_label   <- "Quarter"
  results_suffix <- "Quarterly_By_Year"
  
} else {
  
  stop('period_mode must be either "Monthly" or "Quarterly"')
  
}

# ---- Read data ----

df_all <- read_xlsx(input_path, sheet = data_sheet)

df_all <- df_all %>%
  filter(
    Period_Type == period_mode,
    Year >= 2007,
    Year <= 2022
  ) %>%
  mutate(
    Period = factor(Period, levels = period_levels),
    Year = factor(Year)
  )

# ---- Theoretical HC boundary for Shannon only ----

data("LinfLsup")

bounds <- LinfLsup %>%
  filter(as.integer(as.character(Dimension)) == D)

bounds_shannon_zoom <- bounds %>%
  filter(
    H >= 0.85,
    H <= 1.00,
    C >= 0.00,
    C <= 0.20
  )

bound_group_col <- "Side"

# ---- Year colours ----

year_colors <- c(
  "2007" = "#66C2A5",
  "2008" = "#FC8D62",
  "2009" = "#8DA0CB",
  "2010" = "#E78AC3",
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
  "2022" = "#F781BF"
)

# ---- Common theme ----

hc_theme <- theme_bw(base_size = 11, base_family = "serif") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5)
  )

# ---- Function for one HC plane ----

make_single_hc_plot <- function(df,
                                H_col,
                                C_col,
                                Semi_H_col = NULL,
                                Semi_C_col = NULL,
                                measure_label,
                                add_boundary = FALSE,
                                show_ci = FALSE,
                                show_year_legend = FALSE) {
  
  plot_df <- df %>%
    filter(
      is.finite(.data[[H_col]]),
      is.finite(.data[[C_col]])
    )
  
  p <- ggplot(
    plot_df,
    aes(
      x = .data[[H_col]],
      y = .data[[C_col]],
      color = Year
    )
  )
  
  if (show_ci) {
    
    if (!is.null(Semi_C_col) && Semi_C_col %in% names(plot_df)) {
      p <- p +
        geom_errorbar(
          aes(
            ymin = .data[[C_col]] - .data[[Semi_C_col]],
            ymax = .data[[C_col]] + .data[[Semi_C_col]]
          ),
          width = 0,
          alpha = 0.5,
          linewidth = 0.4,
          show.legend = FALSE
        )
    }
    
    if (!is.null(Semi_H_col) && Semi_H_col %in% names(plot_df)) {
      p <- p +
        geom_errorbarh(
          aes(
            xmin = .data[[H_col]] - .data[[Semi_H_col]],
            xmax = .data[[H_col]] + .data[[Semi_H_col]]
          ),
          height = 0,
          alpha = 0.5,
          linewidth = 0.4,
          show.legend = FALSE
        )
    }
  }
  
  p <- p +
    geom_point(size = 2, alpha = 0.85, show.legend = show_year_legend) +
    scale_color_manual(values = year_colors, drop = FALSE) +
    labs(
      title = measure_label,
      x = expression(italic(H)),
      y = expression(italic(C)),
      color = "Year"
    ) +
    hc_theme +
    guides(
      color = guide_legend(
        nrow = 2,
        override.aes = list(size = 2, alpha = 1, linewidth = 0)
      )
    )
  
  if (!show_year_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  if (add_boundary) {
    p <- p +
      geom_line(
        data = bounds_shannon_zoom,
        aes(x = H, y = C, group = .data[[bound_group_col]]),
        color = "grey35",
        linewidth = 0.6,
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      coord_cartesian(
        xlim = c(0.85, 1.00),
        ylim = c(0.00, 0.20)
      ) +
      scale_x_continuous(
        breaks = c(0.85, 0.90, 0.95, 1.00)
      ) +
      scale_y_continuous(
        breaks = c(0.00, 0.05, 0.10, 0.15, 0.20)
      )
  }
  
  return(p)
}

# ---- Make one combined plot for a given period value (a month or a quarter) ----

make_period_plot <- function(df, period_value, var_label, D, show_ci = FALSE) {
  
  period_df <- df %>%
    filter(
      Period == period_value,
      Variable == var_label
    )
  
  p_shannon <- make_single_hc_plot(
    df = period_df,
    H_col = "H_Shannon",
    C_col = "C_Shannon",
    Semi_H_col = "Semi_H_Shannon",
    Semi_C_col = "Semi_C_Shannon",
    measure_label = "Shannon",
    add_boundary = TRUE,
    show_ci = show_ci,
    show_year_legend = TRUE
  )
  
  p_tsallis <- make_single_hc_plot(
    df = period_df,
    H_col = "H_Tsallis",
    C_col = "C_Tsallis",
    Semi_H_col = "Semi_H_Tsallis",
    Semi_C_col = NULL,
    measure_label = "Tsallis",
    add_boundary = FALSE,
    show_ci = show_ci,
    show_year_legend = FALSE
  )
  
  p_renyi <- make_single_hc_plot(
    df = period_df,
    H_col = "H_Renyi",
    C_col = "C_Renyi",
    Semi_H_col = "Semi_H_Renyi",
    Semi_C_col = NULL,
    measure_label = "Renyi",
    add_boundary = FALSE,
    show_ci = show_ci,
    show_year_legend = FALSE
  )
  
  p_fisher <- make_single_hc_plot(
    df = period_df,
    H_col = "H_Fisher",
    C_col = "C_Fisher",
    Semi_H_col = "Semi_H_Fisher",
    Semi_C_col = NULL,
    measure_label = "Fisher",
    add_boundary = FALSE,
    show_ci = show_ci,
    show_year_legend = FALSE
  )
  
  combined_plot <-
    (p_shannon + p_tsallis) /
    (p_renyi + p_fisher) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = bquote(
        .(period_value) ~ italic(H) %*% italic(C) ~ "planes," ~ italic(D) == .(D) ~
          "(" * .(var_label) * ")"
      ),
      theme = theme(
        plot.title = element_text(
          size = 13,
          face = "bold",
          hjust = 0.5,
          family = "serif"
        ),
        legend.position = "bottom"
      )
    ) &
    theme(legend.position = "bottom")
  
  return(combined_plot)
}

# ---- Output folders ----

out_dir_no_ci   <- here("results", paste0("HC_Results_", results_suffix, "_No_CI"))
out_dir_with_ci <- here("results", paste0("HC_Results_", results_suffix, "_With_CI"))

dir.create(out_dir_no_ci, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_with_ci, recursive = TRUE, showWarnings = FALSE)

# ---- Save plots for each variable, both No-CI and With-CI ----

for (v in unique(df_all$Variable)) {
  
  safe_var <- str_replace_all(v, "[^A-Za-z0-9]+", "_")
  
  for (pv in period_values) {
    
    p_no_ci <- make_period_plot(
      df = df_all,
      period_value = pv,
      var_label = v,
      D = D,
      show_ci = FALSE
    )
    
    ggsave(
      filename = file.path(
        out_dir_no_ci,
        paste0("HC_Combined_", pv, "_", safe_var, "_D", D, "_Years_2007_2022_No_CI.pdf")
      ),
      plot = p_no_ci,
      width = 16,
      height = 12,
      units = "cm",
      device = "pdf"
    )
    
    p_with_ci <- make_period_plot(
      df = df_all,
      period_value = pv,
      var_label = v,
      D = D,
      show_ci = TRUE
    )
    
    ggsave(
      filename = file.path(
        out_dir_with_ci,
        paste0("HC_Combined_", pv, "_", safe_var, "_D", D, "_Years_2007_2022_With_CI.pdf")
      ),
      plot = p_with_ci,
      width = 16,
      height = 12,
      units = "cm",
      device = "pdf"
    )
  }
  
  cat("Saved", length(period_values), period_label, "plots (No-CI and With-CI) for", v, "\n")
}

# End of code

