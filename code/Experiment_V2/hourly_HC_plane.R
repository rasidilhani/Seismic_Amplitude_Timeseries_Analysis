##########################################################
# HC plots for HOURLY WIZ ordinal-pattern data (Jul-Dec 2019)
# Facet dimension = MEASURE (Shannon, Tsallis, Renyi, Fisher)
# No-CI only.
#
# Two outputs:
#   1) ALL-DATA plot: every hourly point (all months, all
#      hours), 4 facets (one per measure), colored by Month.
#   2) MONTHLY plots: one figure per month (e.g. July has
#      ~240 points = 10 sampled days x 24 hours), same 4
#      facets, colored by Hour-of-day to show diurnal
#      structure within that month.
##########################################################

library(readr)
library(tidyverse)
library(lubridate)
library(here)
library(StatOrdPattHxC)

# ---- Parameters ----

D <- 5

input_path <- here("results", "WIZ_OP_Hourly_Every3Days_July_Dec_2019.csv")
out_dir    <- here("results", "HC_Results_hourly")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Read data ----
# (Fisher_Info and C_Fisher are already merged into this file.)

df_all <- read_csv(input_path, show_col_types = FALSE)

month_levels <- c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
measure_levels <- c("Shannon", "Tsallis", "Renyi", "Fisher")

df_all <- df_all %>%
  mutate(
    Date       = mdy(Day),
    Month      = factor(month(Date, label = TRUE, abbr = TRUE),
                        levels = month_levels),
    # FIX: Hour_Start arrives as a full ISO-8601 UTC timestamp
    # (e.g. "2019-07-03T00:00:00Z"), which read_csv() parses as
    # POSIXct. The old `as.integer(Hour_Start)` converted the whole
    # datetime to Unix epoch seconds (a huge number, e.g. 1562112000)
    # instead of the hour-of-day -- every value then fell outside the
    # cut() breaks below and became NA, so Hour_Block was NA for every
    # row and the monthly plots had no color. hour() correctly
    # extracts the 0-23 hour component instead.
    Hour_Start = hour(ymd_hms(Hour_Start, tz = "UTC"))
  ) %>%
  filter(!is.na(Month), !is.na(Hour_Start))

# ---- Reshape to long format: one row per (record, measure) ----
# Shannon_H/Shannon_C, Tsallis_H/Tsallis_C, Renyi_H/Renyi_C,
# Fisher_H/Fisher_C -> pivoted into Measure + H + C columns.

df_long <- df_all %>%
  transmute(
    Month, Hour_Start, Day, Date,
    Shannon_H = H_Shannon,   Shannon_C = C_Shannon,
    Tsallis_H = H_Tsallis,   Tsallis_C = C_Tsallis,
    Renyi_H   = H_Renyi,     Renyi_C   = C_Renyi,
    Fisher_H  = Fisher_Info, Fisher_C  = C_Fisher
  ) %>%
  pivot_longer(
    cols = -c(Month, Hour_Start, Day, Date),
    names_to = c("Measure", ".value"),
    names_sep = "_"
  ) %>%
  mutate(Measure = factor(Measure, levels = measure_levels)) %>%
  filter(is.finite(H), is.finite(C))

# ---- Theoretical HC boundary for Shannon only ----
# Tagged with Measure = "Shannon" so facet_wrap only draws it on that facet.

data("LinfLsup")

bounds <- LinfLsup %>%
  filter(as.integer(as.character(Dimension)) == D)

bounds_shannon_zoom <- bounds %>%
  filter(
    H >= 0.85,
    H <= 1.00,
    C >= 0.00,
    C <= 0.20
  ) %>%
  mutate(Measure = factor("Shannon", levels = measure_levels))

bound_group_col <- "Side"

# ---- Month colours (categorical, used in the all-data plot) ----

month_colors <- c(
  "Jul" = "#440154",
  "Aug" = "darkgreen",
  "Sep" = "green",
  "Oct" = "blue",
  "Nov" = "black",
  "Dec" = "brown"
)

# ---- Hour-block colours (categorical, used in the monthly plots) ----
# Hour_Start binned into 4-hour blocks with hand-picked, high-contrast
# colors (rather than a continuous gradient) so points at opposite HC
# corners stay distinguishable.

hour_block_levels <- c("00-03", "04-07", "08-11", "12-15", "16-19", "20-23")

hour_block_colors <- c(
  "00-03" = "#08306B",
  "04-07" = "#4292C6",
  "08-11" = "#41AB5D",
  "12-15" = "#FDB462",
  "16-19" = "#E31A1C",
  "20-23" = "#6A3D9A"
)

df_long <- df_long %>%
  mutate(
    Hour_Block = cut(
      Hour_Start,
      breaks = c(-1, 3, 7, 11, 15, 19, 23),
      labels = hour_block_levels
    )
  )

# ---- Common theme ----

hc_theme <- theme_classic(base_family = "serif", base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

# ---- Shared boundary layer helper ----

add_shannon_boundary <- function(p) {
  p +
    geom_line(
      data = bounds_shannon_zoom,
      aes(x = H, y = C, group = .data[[bound_group_col]]),
      color = "grey35",
      linewidth = 0.5,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
}

# ---- 1) ALL-DATA plot: 4 facets by Measure, colored by Month ----

p_all <- ggplot(df_long, aes(x = H, y = C, color = Month)) +
  geom_point(size = 1.4, alpha = 0.75) +
  scale_color_manual(values = month_colors, drop = FALSE, name = "Month") +
  facet_wrap(vars(Measure), ncol = 2, scales = "free") +
  labs(
    #title = bquote(italic(H) %*% italic(C) ~ "planes by measure," ~
    #                 italic(D) == .(D) ~ "(WIZ, hourly, Jul-Dec 2019)"),
    x = expression(italic(H)),
    y = expression(italic(C))
  ) +
  hc_theme +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 2, alpha = 1)))

p_all <- add_shannon_boundary(p_all)

print(p_all)   # display in plot window before saving

ggsave(
  filename = file.path(out_dir,
                       sprintf("HC_ByMeasure_AllData_ColorByMonth_D%d_Jul_Dec_2019_No_CI.pdf", D)),
  plot = p_all,
  width = 20, height = 16, units = "cm", device = "pdf"
)

cat("Saved all-data by-measure figure to", out_dir, "\n")

# ---- 2) MONTHLY plots: one figure per month, 4 facets by Measure, colored by Hour ----

make_month_plot <- function(df, month_value, D) {
  
  month_df <- df %>% filter(Month == month_value)
  
  p <- ggplot(month_df, aes(x = H, y = C, color = Hour_Block)) +
    geom_point(size = 1.4, alpha = 0.75) +
    scale_color_manual(values = hour_block_colors, drop = FALSE, name = "Hour") +
    facet_wrap(vars(Measure), ncol = 2, scales = "free") +
    labs(
      title = bquote(.(as.character(month_value)) ~ italic(H) %*% italic(C) ~
                       "planes by measure," ~ italic(D) == .(D) ~
                       "(WIZ, hourly, 2019)"),
      x = expression(italic(H)),
      y = expression(italic(C))
    ) +
    hc_theme
  
  p <- add_shannon_boundary(p)
  
  return(p)
}

for (mv in month_levels) {
  
  month_data_available <- df_long %>% filter(Month == mv) %>% nrow()
  if (month_data_available == 0) next
  
  p_month <- make_month_plot(df_long, month_value = mv, D = D)
  
  print(p_month)   # display in plot window before saving
  
  ggsave(
    filename = file.path(out_dir,
                         sprintf("HC_ByMeasure_%s_ColorByHour_D%d_2019_No_CI.pdf", mv, D)),
    plot = p_month,
    width = 20, height = 16, units = "cm", device = "pdf"
  )
}

cat("Saved monthly by-measure figures to", out_dir, "\n")

# End of code

######################################################
#Check the extreme points

library(dplyr)

# Filter to December only (~240 hourly points)
df_dec <- df_long %>% filter(Month == "Dec")

# For each measure, find the rows with min/max H and min/max C
extreme_points <- df_dec %>%
  filter(Measure %in% c("Tsallis", "Renyi", "Fisher")) %>%
  group_by(Measure) %>%
  summarise(
    min_H_row = list(slice_min(pick(everything()), H, n = 1)),
    max_H_row = list(slice_max(pick(everything()), H, n = 1)),
    min_C_row = list(slice_min(pick(everything()), C, n = 1)),
    max_C_row = list(slice_max(pick(everything()), C, n = 1)),
    .groups = "drop"
  )

# Unpack into a single readable table
extreme_table <- bind_rows(
  df_dec %>% filter(Measure %in% c("Tsallis","Renyi","Fisher")) %>%
    group_by(Measure) %>% slice_min(H, n = 1) %>% mutate(Extreme = "Min H"),
  df_dec %>% filter(Measure %in% c("Tsallis","Renyi","Fisher")) %>%
    group_by(Measure) %>% slice_max(H, n = 1) %>% mutate(Extreme = "Max H"),
  df_dec %>% filter(Measure %in% c("Tsallis","Renyi","Fisher")) %>%
    group_by(Measure) %>% slice_min(C, n = 1) %>% mutate(Extreme = "Min C"),
  df_dec %>% filter(Measure %in% c("Tsallis","Renyi","Fisher")) %>%
    group_by(Measure) %>% slice_max(C, n = 1) %>% mutate(Extreme = "Max C")
) %>%
  select(Measure, Extreme, Day, Hour_Start, H, C) %>%
  arrange(Measure, Extreme)

print(extreme_table, n = Inf)