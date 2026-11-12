# =============================================================================
# PROJECT 12 - Epidemiological Surveillance: Control Chart (Endemic Channel)
# =============================================================================
#
# DESCRIPTION
# Building a control chart to detect epidemics from the weekly incidence
# series of a disease (e.g., arboviruses), using non-epidemic years as the
# reference for the expected pattern (mean ± 2 standard deviations per
# epidemiological week). Also shows how to build the series from individual
# notifications (e.g., SINAN data) and how to query a public surveillance API.
#
# WHAT WAS DONE
#  - Building the query URL for the public InfoDengue API (optional)
#  - Weekly time series (xts) and boxplots by year to identify epidemic years
#  - Endemic channel per epidemiological week
#  - Overlay of the current-year incidence and flagging of weeks above the
#    upper limit
#  - Aggregation of individual notifications into epidemiological weeks and
#    incidence per 100,000 inhabitants
#
# TECHNIQUES / PACKAGES: tidyverse, lubridate, xts, ggplot2
# NOTE: simulated data; variable and object names anonymized.
# =============================================================================

library(tidyverse)
library(lubridate)

set.seed(2357)

# -----------------------------------------------------------------------------
# 1. (Optional) query the public InfoDengue API
# -----------------------------------------------------------------------------
build_url <- function(geocode, disease = "dengue", yr_start = 2012, yr_end = 2022) {
  paste0("https://info.dengue.mat.br/api/alertcity?",
         "geocode=", geocode, "&disease=", disease, "&format=csv",
         "&ew_start=1&ew_end=53&ey_start=", yr_start, "&ey_end=", yr_end)
}
# api_series <- readr::read_csv(build_url(geocode = 0000000))  # IBGE city code

# -----------------------------------------------------------------------------
# 2. Simulated weekly series (epidemic years: 2015, 2016, 2019, 2020, 2022)
# -----------------------------------------------------------------------------
epi_years <- c(2015, 2016, 2019, 2020, 2022)

wk_series <- tibble(wk_start = seq(as.Date("2012-01-01"), as.Date("2022-12-31"), by = "week")) |>
  mutate(
    yr       = year(wk_start),
    ew       = epiweek(wk_start),
    base_s   = 8 + 25 * exp(-((ew - 14)^2) / 60),                 # seasonality
    factor_e = if_else(yr %in% epi_years, runif(n(), 4, 9), 1),
    inc_100k = rgamma(n(), shape = 5, rate = 5 / (base_s * factor_e))
  )

if (requireNamespace("xts", quietly = TRUE)) {
  ts_inc <- xts::xts(wk_series$inc_100k, order.by = wk_series$wk_start)
  plot(ts_inc, main = "Weekly incidence", ylab = "Incidence per 100,000")
}

# Boxplot by year to identify atypical years
boxplot(inc_100k ~ yr, data = wk_series,
        xlab = "Year", ylab = "Incidence per 100,000",
        main = "Distribution of weekly incidence by year")
abline(h = c(300, 100, 50), lty = 2, lwd = 2, col = c("red", "orange", "blue"))

# -----------------------------------------------------------------------------
# 3. Endemic channel using non-epidemic years
# -----------------------------------------------------------------------------
ref_years <- setdiff(2012:2021, epi_years)

endemic_ch <- wk_series |>
  filter(yr %in% ref_years) |>
  group_by(ew) |>
  summarise(
    n_years = n(),
    mean_i  = mean(inc_100k, na.rm = TRUE),
    sd_i    = sd(inc_100k, na.rm = TRUE),
    upper   = mean_i + 2 * sd_i,
    lower   = pmax(mean_i - 2 * sd_i, 0),
    .groups = "drop"
  )

target_yr <- wk_series |>
  filter(yr == 2022) |>
  left_join(endemic_ch, by = "ew") |>
  mutate(alert = inc_100k > upper)

# weeks in an epidemic situation
target_yr |> filter(alert) |> select(ew, inc_100k, upper)

ggplot(endemic_ch, aes(x = ew)) +
  geom_col(data = target_yr, aes(y = inc_100k, fill = "Incidence in the year"), alpha = 0.45) +
  geom_line(aes(y = mean_i, color = "Mean incidence"), linewidth = 1.2) +
  geom_line(aes(y = upper,  color = "Upper limit"),    linewidth = 1.2) +
  geom_point(data = filter(target_yr, alert), aes(y = inc_100k),
             shape = 8, color = "red", size = 2) +
  scale_x_continuous(breaks = seq(1, 53, 2), expand = c(0, 0)) +
  scale_color_manual(name = "", values = c("Mean incidence" = "darkblue",
                                           "Upper limit" = "red")) +
  scale_fill_manual(name = "", values = c("Incidence in the year" = "deepskyblue")) +
  labs(x = "Epidemiological week", y = "Incidence per 100,000",
       title = "Control chart - year 2022") +
  theme_classic() +
  theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# 4. Building the series from individual notifications
# -----------------------------------------------------------------------------
notif_ind <- tibble(dt_notif = as.Date("2016-01-01") + sample(0:2900, 6000, TRUE))
pop_year  <- tibble(iso_yr = 2016:2023,
                    pop_ref = round(seq(1.9e6, 2.2e6, length.out = 8)))

ind_series <- notif_ind |>
  mutate(iso_yr = isoyear(dt_notif), ew = epiweek(dt_notif)) |>
  group_by(iso_yr, ew) |>
  summarise(wk_start = min(dt_notif), cases = n(), .groups = "drop") |>
  left_join(pop_year, by = "iso_yr") |>
  mutate(inc_100k = cases / pop_ref * 1e5) |>
  arrange(iso_yr, ew)

head(ind_series)
