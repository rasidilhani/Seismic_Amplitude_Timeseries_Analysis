#############################################
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)

WIZ_NZ_2011 <- read_csv("data/WIZ_NZ_2011.csv") %>%
  mutate(
    datetime_nz = as.POSIXct(
      unix_timestamp,
      origin = "1970-01-01",
      tz = "Pacific/Auckland"
    ),
    month = month(datetime_nz, label = TRUE, abbr = TRUE)
  )

a <- WIZ_NZ_2011 %>% filter(month == "Apr")
b <- WIZ_NZ_2011 %>% filter(month == "Nov")

p1 <- ggplot(a, aes(x = datetime_nz, y = rsam_avg)) +
  geom_line() +
#  scale_y_log10() +
  ggtitle("April 2011") +
  theme_bw()

p2 <- ggplot(b, aes(x = datetime_nz, y = rsam_avg)) +
  geom_line() +
#  scale_y_log10() +
  ggtitle("November 2011") +
  theme_bw()

p1 + p2



