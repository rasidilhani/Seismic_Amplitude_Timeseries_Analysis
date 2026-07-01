library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)

# 1. Read 2011 data
df <- read_csv("data/WIZ_NZ_2011.csv")

# 2. Create date, year, and month columns
df <- df %>%
  mutate(
    datetime_nz = as.POSIXct(
      unix_timestamp,
      origin = "1970-01-01",
      tz = "Pacific/Auckland"
    ),
    year = year(datetime_nz),
    month = month(datetime_nz, label = TRUE, abbr = TRUE)
  )

# 3. Check number of observations in each month
df %>%
  count(month)

ggplot(df, aes(x = datetime_nz, y = rsam_avg)) +
  geom_line() +
  scale_y_log10() +
  facet_wrap(~ month, scales = "free_x", ncol = 3) +
  labs(
    title = "Monthly RSAM time series: 2011",
    x = "Date",
    y = "RSAM average"
  ) +
  theme_bw()

# Monthly summary Statistics
monthly_summary <- df %>%
  group_by(month) %>%
  summarise(
    n = n(),
    mean = mean(rsam_avg, na.rm = TRUE),
    median = median(rsam_avg, na.rm = TRUE),
    sd = sd(rsam_avg, na.rm = TRUE),
    min = min(rsam_avg, na.rm = TRUE),
    max = max(rsam_avg, na.rm = TRUE),
    cv = sd / mean,
    .groups = "drop"
  )

monthly_summary

# Boxplot comparison
ggplot(df, aes(x = month, y = rsam_avg)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(
    title = "Monthly RSAM boxplots: 2011",
    x = "Month",
    y = "RSAM average"
  ) +
  theme_bw()

# Density plot comparison
ggplot(df, aes(x = rsam_avg, colour = month)) +
  geom_density() +
  scale_x_log10() +
  labs(
    title = "Monthly RSAM density plots: 2011",
    x = "RSAM average",
    y = "Density"
  ) +
  theme_bw()

# Kruksal Wallis test
kruskal.test(rsam_avg ~ month, data = df)

# ACF for one month at a time
apr <- df %>% filter(month == "Apr")

acf(apr$rsam_avg, main = "ACF: April 2011")

nov <- df %>% filter(month == "Nov")

acf(nov$rsam_avg, main = "ACF: November 2011")

#ACF for all months automatically
months_list <- levels(df$month)

for (m in months_list) {
  x <- df %>%
    filter(month == m) %>%
    pull(rsam_avg)
  
  acf(x, main = paste("ACF:", m, "2011"))
}

#Spectral plot for all months
for (m in months_list) {
  x <- df %>%
    filter(month == m) %>%
    pull(rsam_avg)
  
  spec.pgram(x, main = paste("Spectrum:", m, "2011"))
}


######### Yearly comparison of summary statistics
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(purrr)
library(stringr)
library(here)

# 1. Read all yearly CSV files
files <- list.files(
  path = "data",
  pattern = "WIZ_NZ_[0-9]{4}\\.csv$",
  full.names = TRUE
)

seismic_all <- files %>%
  map_dfr(function(file) {
    year_value <- str_extract(file, "[0-9]{4}")
    
    read_csv(file, show_col_types = FALSE) %>%
      mutate(
        year = as.integer(year_value),
        datetime_nz = as.POSIXct(
          unix_timestamp,
          origin = "1970-01-01",
          tz = "Pacific/Auckland"
        ),
        month = month(datetime_nz, label = TRUE, abbr = TRUE)
      )
  })

# 2. Check data
glimpse(seismic_all)

# 3. Yearly summary table for Displacement
displacement_yearly_summary <- seismic_all %>%
  filter(year >= 2011, year <= 2022) %>%
  group_by(year) %>%
  summarise(
    Mean = mean(displacement_avg_m, na.rm = TRUE),
    Median = median(displacement_avg_m, na.rm = TRUE),
    S.d = sd(displacement_avg_m, na.rm = TRUE),
    Min = min(displacement_avg_m, na.rm = TRUE),
    Max = max(displacement_avg_m, na.rm = TRUE),
    Coefficient_of_variation =
      S.d / Mean,
    .groups = "drop"
  )

displacement_yearly_summary


# 4. Yearly summary table for RSAM
rsam_yearly_summary <- seismic_all %>%
  filter(year >= 2011, year <= 2022) %>%
  group_by(year) %>%
  summarise(
    Mean = mean(rsam_avg, na.rm = TRUE),
    Median = median(rsam_avg, na.rm = TRUE),
    S.d = sd(rsam_avg, na.rm = TRUE),
    Min = min(rsam_avg, na.rm = TRUE),
    Max = max(rsam_avg, na.rm = TRUE),
    Coefficient_of_variation =
      S.d / Mean,
    .groups = "drop"
  )

rsam_yearly_summary


# 5. Save yearly summary tables
write_csv(
  displacement_yearly_summary,
  "results/yearly_summary_statistics_displacement.csv"
)

write_csv(
  rsam_yearly_summary,
  "results/yearly_summary_statistics_rsam.csv"
)

# Monthly boxplot for Displacement
p_displacement_box <- ggplot(
  seismic_all,
  aes(x = month, y = displacement_avg_m)
) +
#  geom_boxplot(outlier.alpha = 0.3) +
  geom_boxplot(
    outlier.size = 0.5,
    outlier.alpha = 0.05
  )+
  labs(
    title = "Monthly distribution of displacement, 2011-2022",
    x = "Month",
    y = "Average displacement"
  ) +
  theme_bw(base_size = 12)

p_displacement_box

ggsave(
  "figures/monthly_boxplot_displacement_all_years.pdf",
  p_displacement_box,
  width = 8,
  height = 5
)

# Monthly boxplot for RSAM using all years
p_rsam_box <- ggplot(
  seismic_all,
  aes(x = month, y = rsam_avg)
) +
#  geom_boxplot(outlier.shape = NA) +
  geom_boxplot(
    outlier.size = 0.5,
    outlier.alpha = 0.05
  )+
  scale_y_log10() +
  labs(
    title = "Monthly distribution of RSAM, 2011-2022",
    x = "Month",
    y = "RSAM"
  ) +
  theme_bw(base_size = 12)

p_rsam_box

ggsave(
  "figures/monthly_boxplot_RSAM_all_years.pdf",
  p_rsam_box,
  width = 8,
  height = 5
)

# Boxplot + violin plot
# 1. Displacement 
violin_plot_displacement <- ggplot(seismic_all,
                           aes(month, displacement_avg_m)) +
  geom_violin(
    fill = "grey85",
    colour = "black",
    alpha = 0.7
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white"
  ) +
  scale_y_log10() +
  labs(
    title = "Monthly distribution of Displacement Average, 2011-2022",
    x = "Month",
    y = "Displacement_avg"
  ) +
  theme_bw(base_size = 12, base_family = "serif")

violin_plot_displacement

ggsave(
  "figures/monthly_violin_boxplot_displacement_all_years.pdf",
  violin_plot_displacement,
  width = 8,
  height = 5
)

#2> RSAM violin boxplot
violin_plot_rsam <- ggplot(seismic_all,
       aes(month, rsam_avg)) +
  geom_violin(
    fill = "grey85",
    colour = "black",
    alpha = 0.7
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white"
  ) +
  scale_y_log10() +
  labs(
    title = "Monthly distribution of RSAM Average, 2011-2022",
    x = "Month",
    y = "RSAM_avg"
  ) +
  theme_bw(base_size = 12, base_family = "serif")

violin_plot_rsam

ggsave(
  "figures/monthly_violin_boxplot_RSAM_all_years.pdf",
  violin_plot_rsam,
  width = 8,
  height = 5
)

#####################################################
# Density plots
p_displacement_density <- ggplot(
  seismic_all,
  aes(x = displacement_avg_m)
) +
  geom_density() +
  labs(
    title = "Density plot of displacement, 2011-2022",
    x = "Average displacement",
    y = "Density"
  ) +
  theme_bw(base_size = 12)

p_displacement_density

ggsave(
  "figures/density_displacement_all_years.pdf",
  p_displacement_density,
  width = 8,
  height = 5
)

# RSAM
p_rsam_density <- ggplot(
  seismic_all,
  aes(x = rsam_avg)
) +
  geom_density() +
  scale_x_log10() +
  labs(
    title = "Density plot of RSAM, 2011--2022",
    x = "RSAM",
    y = "Density"
  ) +
  theme_bw(base_size = 12)

p_rsam_density

ggsave(
  "figures/density_RSAM_all_years.pdf",
  p_rsam_density,
  width = 8,
  height = 5
)