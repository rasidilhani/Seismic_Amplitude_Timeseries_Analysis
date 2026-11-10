# ============================================================
# DAILY H-C PLANE PLOTS - RSAM 2019
# Embedding Dimension D = 3
# ============================================================


# ============================================================
# 1. LOAD PACKAGES
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(here)
library(scales)
library(writexl)
library(StatOrdPattHxC)


# ============================================================
# 2. SETTINGS
# ============================================================

# Embedding dimension
D <- 3

# Number of extreme points to label
N_LABEL <- 5


# ============================================================
# 3. FILE LOCATIONS
# ============================================================

# Folder containing the daily results
results_dir <- here(
  "results",
  "Daily_OP_Analysis_2019"
)


# CSV input file
input_file <- here(
  "results",
  "Daily_OP_Analysis_2019",
  "daily_entropy_complexity_rsam_avg_2019.csv"
)


# Output folders
out_noCI <- here(
  "results",
  "Daily_OP_Analysis_2019",
  "HC_Plane_Plots_noCI"
)

out_withCI <- here(
  "results",
  "Daily_OP_Analysis_2019",
  "HC_Plane_Plots_withCI"
)


# Create output folders
dir.create(
  out_noCI,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  out_withCI,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. READ DAILY CSV FILE
# ============================================================

df <- read_csv(
  input_file,
  show_col_types = FALSE
)


# Add a simple daily index
df <- df %>%
  mutate(
    Day_Index = row_number()
  )


# ============================================================
# 5. CONVERT DATA TO LONG FORMAT
# ============================================================

df_long <- bind_rows(
  
  # Shannon
  df %>%
    transmute(
      Day_Index,
      Day_Start,
      Day_End,
      Measure = "Shannon",
      H = H_Shannon,
      C = C_Shannon,
      Semi_H = Semi_H_Shannon,
      Semi_C = Semi_C_Shannon
    ),
  
  # Tsallis
  df %>%
    transmute(
      Day_Index,
      Day_Start,
      Day_End,
      Measure = "Tsallis",
      H = H_Tsallis,
      C = C_Tsallis,
      Semi_H = Semi_H_Tsallis,
      Semi_C = NA_real_
    ),
  
  # Renyi
  df %>%
    transmute(
      Day_Index,
      Day_Start,
      Day_End,
      Measure = "Renyi",
      H = H_Renyi,
      C = C_Renyi,
      Semi_H = Semi_H_Renyi,
      Semi_C = NA_real_
    ),
  
  # Fisher
  df %>%
    transmute(
      Day_Index,
      Day_Start,
      Day_End,
      Measure = "Fisher",
      H = Fisher_Info,
      C = C_Fisher,
      Semi_H = Semi_H_Fisher,
      Semi_C = NA_real_
    )
)

# Set order of facets
df_long$Measure <- factor(
  df_long$Measure,
  levels = c(
    "Shannon",
    "Tsallis",
    "Renyi",
    "Fisher"
  )
)


# Remove missing or non-finite values
df_long <- df_long %>%
  filter(
    is.finite(H),
    is.finite(C)
  )


# ============================================================
# 6. NORMALISE H AND C
#    Used only to identify extreme points
# ============================================================

df_long <- df_long %>%
  
  group_by(Measure) %>%
  
  mutate(
    
    H_norm = rescale(
      H,
      to = c(0, 1)
    ),
    
    C_norm = rescale(
      C,
      to = c(0, 1)
    ),
    
    # Shannon, Tsallis and Renyi:
    # low score = low H / high C
    score = H_norm - C_norm,
    
    # Fisher:
    # high score = high H / high C
    score_fisher = H_norm + C_norm
  ) %>%
  
  ungroup()


# ============================================================
# 7. FIND EXTREME DAILY POINTS
# ============================================================

# Shannon, Tsallis and Renyi:
# select low-H / high-C points

label_non_fisher <- df_long %>%
  
  filter(
    Measure != "Fisher"
  ) %>%
  
  group_by(Measure) %>%
  
  slice_min(
    score,
    n = N_LABEL,
    with_ties = FALSE
  ) %>%
  
  ungroup()


# Fisher:
# select high-H / high-C points

label_fisher <- df_long %>%
  
  filter(
    Measure == "Fisher"
  ) %>%
  
  slice_max(
    score_fisher,
    n = N_LABEL,
    with_ties = FALSE
  )


# Combine labelled points
label_points <- bind_rows(
  label_non_fisher,
  label_fisher
)


# ============================================================
# 8. SHANNON THEORETICAL H-C BOUNDARY FOR D = 3
# ============================================================

data("LinfLsup")


shannon_bounds <- LinfLsup %>%
  
  filter(
    as.integer(
      as.character(Dimension)
    ) == D
  ) %>%
  
  mutate(
    Measure = factor(
      "Shannon",
      levels = levels(df_long$Measure)
    )
  )


# ============================================================
# 9. CROP SHANNON BOUNDARY TO RELEVANT DATA REGION
#    Same approach used in the previous H-C plane code
# ============================================================

crop_bounds_to_data <- function(
    bounds,
    shannon_df,
    pad_frac = 0.05
) {
  
  # Range of observed Shannon entropy and complexity
  H_range <- range(
    shannon_df$H,
    na.rm = TRUE
  )
  
  C_range <- range(
    shannon_df$C,
    na.rm = TRUE
  )
  
  
  # Small margin around the observed data
  H_pad <- diff(H_range) * pad_frac
  
  C_pad <- diff(C_range) * pad_frac
  
  
  # Safety check for very small ranges
  if (!is.finite(H_pad) || H_pad == 0) {
    H_pad <- 0.01
  }
  
  if (!is.finite(C_pad) || C_pad == 0) {
    C_pad <- 0.01
  }
  
  
  # Keep only the part of the theoretical boundary
  # relevant to the observed Shannon data
  bounds %>%
    
    filter(
      
      H >= H_range[1] - H_pad,
      
      H <= H_range[2] + H_pad,
      
      C >= max(
        0,
        C_range[1] - C_pad
      ),
      
      C <= C_range[2] + C_pad
    )
}


# Extract Shannon observations
shannon_data <- df_long %>%
  
  filter(
    Measure == "Shannon"
  )


# Crop theoretical boundary
shannon_bounds_cropped <- crop_bounds_to_data(
  shannon_bounds,
  shannon_data
)


# ============================================================
# 10. COMMON H-C PLOT
# ============================================================

p_base <- ggplot(
  df_long,
  aes(
    x = H,
    y = C,
    color = Day_Index
  )
) +
  
  geom_point(
    size = 1.8,
    alpha = 0.8
  ) +
  
  scale_color_viridis_c(
    name = "Day"
  ) +
  
  facet_wrap(
    ~ Measure,
    ncol = 2,
    scales = "free"
  ) +
  
  
  # ----------------------------------------------------------
# Cropped Shannon theoretical boundary
# ----------------------------------------------------------

geom_line(
  
  data = shannon_bounds_cropped,
  
  aes(
    x = H,
    y = C,
    group = Side
  ),
  
  inherit.aes = FALSE,
  
  color = "grey40",
  
  linewidth = 0.5
) +
  
  
  # ----------------------------------------------------------
# Circle extreme points
# ----------------------------------------------------------

geom_point(
  
  data = label_points,
  
  aes(
    x = H,
    y = C
  ),
  
  inherit.aes = FALSE,
  
  shape = 21,
  
  size = 3,
  
  stroke = 1,
  
  color = "black",
  
  fill = NA
) +
  
  
  # ----------------------------------------------------------
# Label extreme points
# ----------------------------------------------------------

geom_text_repel(
  
  data = label_points,
  
  aes(
    x = H,
    y = C,
    label = paste0(
      "day ",
      Day_Index
    )
  ),
  
  inherit.aes = FALSE,
  
  color = "black",
  
  size = 3,
  
  fontface = "bold",
  
  min.segment.length = 0,
  
  max.overlaps = Inf,
  
  seed = 42
) +
  
  
  labs(
    title = "Daily RSAM H-C Plane - 2019",
    x = expression(italic(H)),
    y = expression(italic(C))
  ) +
  
  
  theme_classic(
    base_family = "serif",
    base_size = 12
  ) +
  
  
  theme(
    
    strip.background =
      element_rect(
        fill = "grey95",
        color = "black"
      ),
    
    strip.text =
      element_text(
        face = "bold"
      ),
    
    legend.position = "bottom",
    
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      )
  )


# ============================================================
# 11. H-C PLOT WITHOUT CONFIDENCE INTERVALS
# ============================================================

p_noCI <- p_base


print(
  p_noCI
)


ggsave(
  filename = file.path(
    out_noCI,
    "HC_Plane_RSAM_D3_noCI.pdf"
  ),
  plot = p_noCI,
  width = 20,
  height = 16,
  units = "cm"
)


# ============================================================
# 12. H-C PLOT WITH CONFIDENCE INTERVALS
# ============================================================

p_withCI <- p_base +
  
  
  # ----------------------------------------------------------
# Horizontal CI for entropy
# ----------------------------------------------------------

geom_segment(
  
  data = df_long %>%
    filter(
      is.finite(Semi_H),
      Semi_H > 0
    ),
  
  aes(
    
    x = pmax(
      0,
      H - Semi_H
    ),
    
    xend = H + Semi_H,
    
    y = C,
    
    yend = C
  ),
  
  linewidth = 0.3,
  
  alpha = 0.6
) +
  
  
  # ----------------------------------------------------------
# Vertical CI for Shannon complexity
# ----------------------------------------------------------

geom_segment(
  
  data = df_long %>%
    filter(
      Measure == "Shannon",
      is.finite(Semi_C),
      Semi_C > 0
    ),
  
  aes(
    
    x = H,
    
    xend = H,
    
    y = pmax(
      0,
      C - Semi_C
    ),
    
    yend = C + Semi_C
  ),
  
  linewidth = 0.3,
  
  alpha = 0.6
)


print(
  p_withCI
)


ggsave(
  filename = file.path(
    out_withCI,
    "HC_Plane_RSAM_D3_withCI.pdf"
  ),
  plot = p_withCI,
  width = 20,
  height = 16,
  units = "cm"
)


# ============================================================
# 13. CHECK AND SAVE EXTREME DAILY POINTS
# ============================================================

extreme_points <- label_points %>%
  select(
    Day_Index,
    Day_Start,
    Day_End,
    Measure,
    H,
    C,
    score,
    score_fisher
  ) %>%
  arrange(
    Measure,
    Day_Index
  )


print(extreme_points)


write_xlsx(
  extreme_points,
  file.path(
    results_dir,
    "Extreme_Daily_Points_D3_2019.xlsx"
  )
)

# ============================================================
# 14. FINISHED
# ============================================================

message(
  "Daily D = 3 H-C plots completed."
)
