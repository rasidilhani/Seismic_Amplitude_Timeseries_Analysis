# HC plots for each month, by variable
# Each plot combines Shannon, Tsallis, Renyi, and Fisher as facets
# Legend is labelled by Year
# No confidence intervals

library(readxl)
library(tidyverse)
library(here)
library(StatOrdPattHxC)

# ---- Parameters ----

D <- 5

input_path <- here("results", "WIZ_OrdinalPatterns_D5_Monthly_Quarterly.xlsx")

# ---- Read data ----

df_all <- read_xlsx(input_path, sheet = 1)

# Keep only monthly data from 2011 to 2021
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
  filter(as.integer(as.character(Dimension)) == D) %>%
  filter(
    H >= 0.90,
    H <= 1.00,
    C >= 0.00,
    C <= 0.20
  )

bound_group_col <- "Side"

# ---- Year colours ----

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
  "2021" = "#4DAF4A"
)

# ---- Convert Shannon, Tsallis, Renyi, Fisher into one long table ----

make_long_measure_data <- function(df) {
  
  df %>%
    transmute(
      Year,
      Period,
      Variable,
      
      H_Shannon,
      C_Shannon,
      H_Tsallis,
      C_Tsallis,
      H_Renyi,
      C_Renyi,
      H_Fisher,
      C_Fisher
    ) %>%
    pivot_longer(
      cols = c(
        H_Shannon, C_Shannon,
        H_Tsallis, C_Tsallis,
        H_Renyi, C_Renyi,
        H_Fisher, C_Fisher
      ),
      names_to = c(".value", "Measure"),
      names_pattern = "([HC])_(.*)"
    ) %>%
    mutate(
      Measure = factor(
        Measure,
        levels = c("Shannon", "Tsallis", "Renyi", "Fisher")
      )
    ) %>%
    filter(
      is.finite(H),
      is.finite(C)
    )
}

# ---- Make one faceted monthly plot ----

make_month_plot <- function(df, month_name, var_label, D) {
  
  plot_df <- df %>%
    filter(
      Period == month_name,
      Variable == var_label
    ) %>%
    make_long_measure_data()
  
  shannon_bounds <- bounds %>%
    mutate(
      Measure = factor(
        "Shannon",
        levels = c("Shannon", "Tsallis", "Renyi", "Fisher")
      )
    )
  
  p <- ggplot(plot_df, aes(x = H, y = C, color = Year)) +
    geom_point(size = 2, alpha = 0.85) +
    geom_line(
      data = shannon_bounds,
      aes(x = H, y = C, group = .data[[bound_group_col]]),
      color = "grey35",
      linewidth = 0.6,
      inherit.aes = FALSE
    ) +
    facet_wrap(~ Measure, scales = "free", ncol = 2) +
    scale_color_manual(values = year_colors, drop = FALSE) +
    labs(
      title = bquote(
        .(month_name) ~ italic(H) %*% italic(C) ~ "planes," ~ italic(D) == .(D) ~
          "(" * .(var_label) * ")"
      ),
      x = expression(italic(H)),
      y = expression(italic(C)),
      color = "Year"
    ) +
    theme_bw(base_size = 11, base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.border = element_rect(color = "black", linewidth = 0.5),
      strip.background = element_rect(fill = "grey92", color = "black"),
      strip.text = element_text(size = 9, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8, color = "black"),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    ) +
    guides(
      color = guide_legend(
        nrow = 2,
        override.aes = list(size = 2, alpha = 1)
      )
    )
  
  return(p)
}

# ---- Save monthly plots for each variable ----

out_dir <- here("results", "HC_Results_Monthly_By_Year_No_CI")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (v in unique(df_all$Variable)) {
  
  for (m in month.abb) {
    
    p <- make_month_plot(
      df = df_all,
      month_name = m,
      var_label = v,
      D = D
    )
    
    print(p)
    
    safe_var <- str_replace_all(v, "[^A-Za-z0-9]+", "_")
    
    ggsave(
      filename = file.path(
        out_dir,
        paste0("HC_Combined_", m, "_", safe_var, "_D", D, "_Years_2011_2021.pdf")
      ),
      plot = p,
      width = 16,
      height = 12,
      units = "cm",
      device = "pdf"
    )
  }
  
  cat("Saved 12 monthly combined plots for", v, "to", out_dir, "\n")
}

# End of code
#######################################################################
# HC plots for each month, by variable
# One combined plot per month:
# Shannon, Tsallis, Renyi, and Fisher
# Legend labelled by Year
# WITH confidence intervals
# Shannon panel zoomed to H = 0.85-1.00 and C = 0.00-0.20

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
    H >= 0.85,
    H <= 1.00,
    C >= 0.00,
    C <= 0.20
  )

bound_group_col <- "Side"

# ---- Year colours ----

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
  "2021" = "#4DAF4A"
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

# ---- Function for one HC plane with confidence intervals ----

make_single_hc_plot <- function(df,
                                H_col,
                                C_col,
                                Semi_H_col = NULL,
                                Semi_C_col = NULL,
                                measure_label,
                                add_boundary = FALSE,
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
  
  # Vertical confidence interval for C
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
  
  # Horizontal confidence interval for H
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

# ---- Make one combined monthly plot ----

make_month_plot <- function(df, month_name, var_label, D) {
  
  month_df <- df %>%
    filter(
      Period == month_name,
      Variable == var_label
    )
  
  p_shannon <- make_single_hc_plot(
    df = month_df,
    H_col = "H_Shannon",
    C_col = "C_Shannon",
    Semi_H_col = "Semi_H_Shannon",
    Semi_C_col = "Semi_C_Shannon",
    measure_label = "Shannon",
    add_boundary = TRUE,
    show_year_legend = TRUE
  )
  
  p_tsallis <- make_single_hc_plot(
    df = month_df,
    H_col = "H_Tsallis",
    C_col = "C_Tsallis",
    Semi_H_col = "Semi_H_Tsallis",
    Semi_C_col = NULL,
    measure_label = "Tsallis",
    add_boundary = FALSE,
    show_year_legend = FALSE
  )
  
  p_renyi <- make_single_hc_plot(
    df = month_df,
    H_col = "H_Renyi",
    C_col = "C_Renyi",
    Semi_H_col = "Semi_H_Renyi",
    Semi_C_col = NULL,
    measure_label = "Renyi",
    add_boundary = FALSE,
    show_year_legend = FALSE
  )
  
  p_fisher <- make_single_hc_plot(
    df = month_df,
    H_col = "H_Fisher",
    C_col = "C_Fisher",
    Semi_H_col = "Semi_H_Fisher",
    Semi_C_col = NULL,
    measure_label = "Fisher",
    add_boundary = FALSE,
    show_year_legend = FALSE
  )
  
  combined_plot <-
    (p_shannon + p_tsallis) /
    (p_renyi + p_fisher) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = bquote(
        .(month_name) ~ italic(H) %*% italic(C) ~ "planes," ~ italic(D) == .(D) ~
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

# ---- Save monthly plots for each variable ----

out_dir <- here("results", "HC_Results_Monthly_By_Year_With_CI")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (v in unique(df_all$Variable)) {
  
  for (m in month.abb) {
    
    p <- make_month_plot(
      df = df_all,
      month_name = m,
      var_label = v,
      D = D
    )
    
    print(p)
    
    safe_var <- str_replace_all(v, "[^A-Za-z0-9]+", "_")
    
    ggsave(
      filename = file.path(
        out_dir,
        paste0("HC_Combined_", m, "_", safe_var, "_D", D, "_Years_2011_2021_With_CI.pdf")
      ),
      plot = p,
      width = 16,
      height = 12,
      units = "cm",
      device = "pdf"
    )
  }
  
  cat("Saved 12 monthly combined plots with CI for", v, "to", out_dir, "\n")
}

# End of code