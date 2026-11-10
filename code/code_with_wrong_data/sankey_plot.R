library(tidyverse)
library(ggalluvial)
library(here)

# ── Read data ─────────────────────────────────────────────────────────────
DATA_DIR <- here("Data")

val_df <- read_csv(
  file.path(DATA_DIR, "VAL_WhakaariWhiteIsland.csv"),
  show_col_types = FALSE
)

# ── Colours for VAL levels ────────────────────────────────────────────────
val_colors <- c(
  "VAL 1" = "#4DAF4A",
  "VAL 2" = "#FFB000",
  "VAL 3" = "#E66101",
  "VAL 4" = "#B2182B"
)

# ── Sankey 1: VAL → system version → description ─────────────────────────
sankey_df <- val_df %>%
  transmute(
    VAL = paste("VAL", VAL),
    System = paste("System version", `VAL System Version`),
    Description = Description
  ) %>%
  count(
    VAL,
    System,
    Description,
    name = "Number_of_periods"
  )

p1 <- ggplot(
  sankey_df,
  aes(
    axis1 = VAL,
    axis2 = System,
    axis3 = Description,
    y = Number_of_periods
  )
) +
  geom_alluvium(
    aes(fill = VAL),
    width = 0.15,
    alpha = 0.75
  ) +
  geom_stratum(
    width = 0.15,
    fill = "grey95",
    color = "grey35",
    linewidth = 0.3
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 2.6,
    family = "serif"
  ) +
  scale_x_discrete(
    limits = c(
      "Volcanic Alert Level",
      "VAL system",
      "Activity description"
    ),
    expand = c(0.08, 0.08)
  ) +
  scale_fill_manual(values = val_colors) +
  labs(
    title = "VAL system and activity descriptions",
    x = NULL,
    y = "Number of recorded VAL periods",
    fill = NULL
  ) +
  theme_minimal(
    base_size = 9,
    base_family = "serif"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(color = "black"),
    legend.position = "bottom"
  )

print(p1)

# ── Create broader activity classes ───────────────────────────────────────
val_df <- val_df %>%
  mutate(
    Activity_Class = case_when(
      Description == "Moderate volcanic eruption" ~
        "Moderate eruption",
      
      Description %in% c(
        "Minor volcanic eruption",
        "Minor eruptive activity"
      ) ~
        "Minor eruption/activity",
      
      Description %in% c(
        "Signs of volcano unrest",
        "Minor volcanic unrest",
        "Moderate to heightened volcanic unrest"
      ) ~
        "Unrest",
      
      TRUE ~ "Other"
    ),
    Activity_Class = factor(
      Activity_Class,
      levels = c(
        "Unrest",
        "Minor eruption/activity",
        "Moderate eruption",
        "Other"
      )
    )
  )

# ── Sankey 2: VAL → activity class ────────────────────────────────────────
activity_sankey <- val_df %>%
  transmute(
    VAL = paste("VAL", VAL),
    Activity_Class
  ) %>%
  count(
    VAL,
    Activity_Class,
    name = "Number_of_periods"
  )

p2 <- ggplot(
  activity_sankey,
  aes(
    axis1 = VAL,
    axis2 = Activity_Class,
    y = Number_of_periods
  )
) +
  geom_alluvium(
    aes(fill = VAL),
    width = 0.2,
    alpha = 0.8
  ) +
  geom_stratum(
    width = 0.2,
    fill = "grey95",
    color = "grey30"
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3,
    family = "serif"
  ) +
  scale_x_discrete(
    limits = c(
      "Volcanic Alert Level",
      "Activity class"
    ),
    expand = c(0.15, 0.15)
  ) +
  scale_fill_manual(values = val_colors) +
  labs(
    title = "VAL and volcanic activity classification",
    x = NULL,
    y = "Number of recorded VAL periods",
    fill = NULL
  ) +
  theme_minimal(base_family = "serif") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p2)

ggsave(
  "VAL_System_Description_Sankey.pdf",
  plot = p1,
  width = 7,
  height = 4.5
)

ggsave(
  "VAL_Activity_Class_Sankey.pdf",
  plot = p2,
  width = 6,
  height = 4
)