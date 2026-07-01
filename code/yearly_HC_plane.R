library(readxl)
library(tidyverse)
library(here)
library(StatOrdPattHxC)
library(patchwork)

# ---- Parameters ----

D <- 5

input_path <- here("results", "WIZ_OrdinalPatterns_D5_Monthly_Quarterly.xlsx")

# ---- Read data ----

df_all <- read_xlsx(input_path, sheet = 1)

df_all <- df_all %>%
  filter(
    Period_Type == "Monthly",
    Year >= 2011,
    Year <= 2021
  ) %>%
  mutate(
    Period = factor(Period, levels = month.abb),
    Year = factor(Year)
  )

# ---- Theoretical HC boundary for Shannon only ----

data("LinfLsup")

bounds <- LinfLsup %>%
  filter(as.integer(as.character(Dimension)) == D)

bounds_shannon_zoom <- bounds %>%
  filter(
    H >= 0.90,
    H <= 1.00,
    C >= 0.00,
    C <= 0.20
  )

bound_group_col <- "Side"

# ---- Month colours ----

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
                                show_month_legend = FALSE) {
  
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
      color = Period
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
    geom_point(size = 2, alpha = 0.85, show.legend = show_month_legend) +
    scale_color_manual(values = month_colors, drop = FALSE) +
    labs(
      title = measure_label,
      x = expression(italic(H)),
      y = expression(italic(C)),
      color = "Month"
    ) +
    hc_theme +
    guides(
      color = guide_legend(
        nrow = 2,
        override.aes = list(size = 2, alpha = 1, linewidth = 0)
      )
    )
  
  if (!show_month_legend) {
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
        xlim = c(0.90, 1.00),
        ylim = c(0.00, 0.20)
      ) +
      scale_x_continuous(
        breaks = c(0.90, 0.95, 1.00)
      ) +
      scale_y_continuous(
        breaks = c(0.00, 0.05, 0.10, 0.15, 0.20)
      )
  }
  
  return(p)
}

# ---- Make one combined yearly plot ----

make_year_plot <- function(df, year_value, var_label, D, show_ci = FALSE) {
  
  year_df <- df %>%
    filter(
      Year == year_value,
      Variable == var_label
    )
  
  p_shannon <- make_single_hc_plot(
    df = year_df,
    H_col = "H_Shannon",
    C_col = "C_Shannon",
    Semi_H_col = "Semi_H_Shannon",
    Semi_C_col = "Semi_C_Shannon",
    measure_label = "Shannon",
    add_boundary = TRUE,
    show_ci = show_ci,
    show_month_legend = TRUE
  )
  
  p_tsallis <- make_single_hc_plot(
    df = year_df,
    H_col = "H_Tsallis",
    C_col = "C_Tsallis",
    Semi_H_col = "Semi_H_Tsallis",
    Semi_C_col = NULL,
    measure_label = "Tsallis",
    add_boundary = FALSE,
    show_ci = show_ci,
    show_month_legend = FALSE
  )
  
  p_renyi <- make_single_hc_plot(
    df = year_df,
    H_col = "H_Renyi",
    C_col = "C_Renyi",
    Semi_H_col = "Semi_H_Renyi",
    Semi_C_col = NULL,
    measure_label = "Renyi",
    add_boundary = FALSE,
    show_ci = show_ci,
    show_month_legend = FALSE
  )
  
  p_fisher <- make_single_hc_plot(
    df = year_df,
    H_col = "H_Fisher",
    C_col = "C_Fisher",
    Semi_H_col = "Semi_H_Fisher",
    Semi_C_col = NULL,
    measure_label = "Fisher",
    add_boundary = FALSE,
    show_ci = show_ci,
    show_month_legend = FALSE
  )
  
  combined_plot <-
    (p_shannon + p_tsallis) /
    (p_renyi + p_fisher) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = bquote(
        .(as.character(year_value)) ~ italic(H) %*% italic(C) ~ "planes," ~
          italic(D) == .(D) ~ "(" * .(var_label) * ")"
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

out_dir_no_ci <- here("results", "HC_Results_Yearly_By_Month_No_CI")
out_dir_with_ci <- here("results", "HC_Results_Yearly_By_Month_With_CI")

dir.create(out_dir_no_ci, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_with_ci, recursive = TRUE, showWarnings = FALSE)

# ---- Save yearly plots for each variable ----

for (v in unique(df_all$Variable)) {
  
  safe_var <- str_replace_all(v, "[^A-Za-z0-9]+", "_")
  
  for (yr in levels(df_all$Year)) {
    
    p_no_ci <- make_year_plot(
      df = df_all,
      year_value = yr,
      var_label = v,
      D = D,
      show_ci = FALSE
    )
    
    ggsave(
      filename = file.path(
        out_dir_no_ci,
        paste0("HC_Combined_", yr, "_", safe_var, "_D", D, "_Months_No_CI.pdf")
      ),
      plot = p_no_ci,
      width = 16,
      height = 12,
      units = "cm",
      device = "pdf"
    )
    
    p_with_ci <- make_year_plot(
      df = df_all,
      year_value = yr,
      var_label = v,
      D = D,
      show_ci = TRUE
    )
    
    ggsave(
      filename = file.path(
        out_dir_with_ci,
        paste0("HC_Combined_", yr, "_", safe_var, "_D", D, "_Months_With_CI.pdf")
      ),
      plot = p_with_ci,
      width = 16,
      height = 12,
      units = "cm",
      device = "pdf"
    )
  }
  
  cat("Saved yearly plots with and without CI for", v, "\n")
}

# End of code