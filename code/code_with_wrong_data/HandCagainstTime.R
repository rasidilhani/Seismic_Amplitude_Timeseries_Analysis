library(readr)
library(ggplot2)
library(dplyr)

# ---- Load data ----
HC_2007to2019 <- read_csv("results/Weekly_OP_Analysis_2007to2019/weekly_entropy_complexity_rsam_avg_2007to2019.csv")
df <- data.frame(HC_2007to2019)

# ---- Make sure Year is numeric/integer ----
df$Year <- as.integer(df$Year)

# ---- Build a continuous week index (1, 2, 3, ... like the Excel row axis) ----
df$week_index <- seq_len(nrow(df))

# ---- Year labels placed at the FIRST week_index of each year (like Excel's 1,29,57...) ----
year_starts <- df %>%
  group_by(Year) %>%
  summarise(first_index = min(week_index)) %>%
  arrange(Year)

year_breaks <- year_starts$first_index
year_labels <- year_starts$Year

# ---- Plot 1: Entropy across time, x-axis labeled by year ----
p_entropy <- ggplot(df, aes(x = week_index, y = H_Shannon)) +
  geom_line(color = "darkblue", linewidth = 0.5) +
  scale_x_continuous(breaks = year_breaks, labels = as.character(year_labels)) +
  labs(
    x = "Year",
    y = "Shannon Entropy"
    #title = "Shannon Entropy"
  ) +
  theme_bw(base_size = 13, base_family = "serif") +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# ---- Plot 2: Complexity across time, x-axis labeled by year ----
p_complexity <- ggplot(df, aes(x = week_index, y = C_Shannon)) +
  geom_line(color = "black", linewidth = 0.5) +
  scale_x_continuous(breaks = year_breaks, labels = as.character(year_labels)) +
  labs(
    x = "Year",
    y = "Complexity"
    #title = "Complexity"
  ) +
  theme_bw(base_size = 13, base_family = "serif") +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# ---- Print to check in RStudio viewer ----
print(p_entropy)
print(p_complexity)

# ---- Save as separate PDFs ----
ggsave("results/H_vs_year.pdf", p_entropy, width = 8, height = 5)
ggsave("results/C_vs_year.pdf", p_complexity, width = 8, height = 5)
