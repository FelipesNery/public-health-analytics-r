# =============================================================================
# PROJECT 13 - Epidemiological Indicators: Incidence, Mortality and Years of
#              Life Lost (YLL)
# =============================================================================
#
# DESCRIPTION
# Calculating disease burden indicators from death records (Brazilian
# Mortality Information System - SIM/DATASUS format) for diseases of the
# circulatory system (ICD-10 chapter IX) in states of the Legal Amazon over
# two decades: mortality rates with 95% CI, temporal trends by sex,
# distribution by cause group, and YLL using state- and year-specific life
# expectancy.
#
# WHAT WAS DONE
#  - Grouping underlying causes of death (ICD-10) into cause groups
#  - Calculation of age at death and year of death
#  - Mortality rate per 100,000 with 95% CI (Poisson approximation)
#  - Smoothed temporal trend (LOESS) by sex
#  - YLL = life expectancy (state, year) - age at death
#  - Municipality ranking and percentage distribution by cause group
#  - Rate by cause group and year as a heat map
#
# TECHNIQUES / PACKAGES: tidyverse, lubridate, ggplot2, scales
# NOTE: simulated data; variable and object names anonymized.
# =============================================================================

library(tidyverse)
library(lubridate)

set.seed(3571)

# -----------------------------------------------------------------------------
# 1. Simulated data
# -----------------------------------------------------------------------------
states <- c("ST_A", "ST_B", "ST_C", "ST_D", "ST_E")
icd_pool <- c("I10", "I11", "I21", "I25", "I48", "I50", "I63", "I64", "I69", "I70", "I80")

n_d <- 15000
deaths_sm <- tibble(
  state_res = sample(states, n_d, TRUE, prob = c(.35, .25, .2, .12, .08)),
  mun_res   = paste0(state_res, "_M", sprintf("%02d", sample(1:30, n_d, TRUE))),
  dt_birth  = as.Date("1920-01-01") + sample(0:30000, n_d, TRUE),
  dt_death  = as.Date("2003-01-01") + sample(0:7660, n_d, TRUE),
  sex_d     = sample(c(1, 2), n_d, TRUE, prob = c(0.55, 0.45)),
  icd_cause = sample(icd_pool, n_d, TRUE)
) |>
  filter(dt_birth < dt_death)

pop_st <- expand_grid(state_res = states, yr_death = 2003:2023, sex_d = c(1, 2)) |>
  mutate(pop_ref = round(runif(n(), 4e5, 2.2e6)))

life_exp <- expand_grid(state_res = states, yr_death = 2003:2023) |>
  mutate(le_years = round(70 + 0.2 * (yr_death - 2003) + rnorm(n(), 0, 0.8), 1))

# -----------------------------------------------------------------------------
# 2. Derived variables and cause groups
# -----------------------------------------------------------------------------
deaths_sm <- deaths_sm |>
  mutate(
    yr_death  = year(dt_death),
    age_death = as.numeric(floor(time_length(interval(dt_birth, dt_death), "years"))),
    icd3      = substr(icd_cause, 1, 3),
    cause_grp = case_when(
      icd3 %in% c("I10", "I11", "I12", "I13", "I15")          ~ "Hypertensive diseases",
      icd3 %in% c("I20", "I21", "I22", "I23", "I24", "I25")   ~ "Ischaemic heart diseases",
      icd3 %in% sprintf("I%d", 30:52)                         ~ "Other forms of heart disease",
      icd3 %in% sprintf("I%d", 60:69)                         ~ "Cerebrovascular diseases",
      icd3 %in% sprintf("I%d", 70:79)                         ~ "Diseases of arteries",
      icd3 %in% sprintf("I%d", 80:89)                         ~ "Diseases of veins and lymphatics",
      TRUE                                                    ~ "Other circulatory disorders"
    )
  )

# -----------------------------------------------------------------------------
# 3. Mortality rate by sex and year with 95% CI
# -----------------------------------------------------------------------------
mort_rate <- deaths_sm |>
  count(state_res, yr_death, sex_d, name = "n_deaths") |>
  inner_join(pop_st, by = c("state_res", "yr_death", "sex_d")) |>
  group_by(yr_death, sex_d) |>
  summarise(n_deaths = sum(n_deaths), pop_ref = sum(pop_ref), .groups = "drop") |>
  mutate(
    rate_100k = n_deaths / pop_ref * 1e5,
    se_rate   = sqrt(n_deaths) / pop_ref * 1e5,
    ci_low    = rate_100k - 1.96 * se_rate,
    ci_high   = rate_100k + 1.96 * se_rate,
    sex_d     = factor(sex_d, levels = c(1, 2), labels = c("Male", "Female"))
  )

ggplot(mort_rate, aes(x = yr_death, y = rate_100k, color = sex_d, fill = sex_d)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_point(size = 1.8) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 1.2) +
  scale_color_manual(values = c("#4682b4", "#e9967a"), name = "Sex") +
  scale_fill_manual(values = c("#add8e6", "#ffdfd4"), name = "Sex") +
  scale_x_continuous(breaks = seq(2003, 2023, 2)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  labs(title = "Mortality from diseases of the circulatory system, by sex",
       x = "Year", y = "Mortality (per 100,000)") +
  theme_classic()

# -----------------------------------------------------------------------------
# 4. Years of life lost (YLL)
# -----------------------------------------------------------------------------
yll_df <- deaths_sm |>
  inner_join(life_exp, by = c("state_res", "yr_death")) |>
  mutate(yll_ind = le_years - age_death) |>
  filter(yll_ind > 0)                         # premature deaths only

# YLL by state and cause group
yll_df |>
  group_by(state_res, cause_grp) |>
  summarise(n_deaths = n(), yll_total = sum(yll_ind), yll_mean = mean(yll_ind),
            .groups = "drop") |>
  arrange(desc(yll_total))

# YLL by municipality and year (for maps / time series)
yll_mun <- yll_df |>
  group_by(yr_death, mun_res, cause_grp) |>
  summarise(yll_total = sum(yll_ind), .groups = "drop")
# write_csv(yll_mun, "yll_municipalities.csv")

# -----------------------------------------------------------------------------
# 5. Percentage distribution by cause group
# -----------------------------------------------------------------------------
deaths_sm |>
  count(cause_grp) |>
  mutate(pct = n / sum(n) * 100, lbl = paste0(round(pct, 1), "%")) |>
  ggplot(aes(x = reorder(cause_grp, pct), y = pct)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_text(aes(label = lbl), hjust = -0.1) +
  coord_flip() +
  expand_limits(y = 40) +
  labs(x = "Cause group", y = "%", title = "Deaths by cause group (ICD-10 chapter IX)") +
  theme_classic()

# -----------------------------------------------------------------------------
# 6. Top 10 municipalities by state
# -----------------------------------------------------------------------------
top_mun <- function(df, state_sel, k = 10) {
  df |>
    filter(state_res == state_sel) |>
    count(mun_res, name = "total") |>
    slice_max(total, n = k) |>
    mutate(mun_res = fct_reorder(mun_res, total))
}

top_mun(deaths_sm, "ST_A") |>
  ggplot(aes(x = total, y = mun_res)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = total), hjust = -0.1, size = 3.5) +
  labs(title = "Municipalities with the most cumulative deaths - ST_A",
       x = "Number of deaths", y = "Municipality") +
  theme_classic()

# -----------------------------------------------------------------------------
# 7. Heat map: rate by cause group and year
# -----------------------------------------------------------------------------
deaths_sm |>
  count(yr_death, cause_grp, name = "n_deaths") |>
  left_join(pop_st |> group_by(yr_death) |> summarise(pop_ref = sum(pop_ref)),
            by = "yr_death") |>
  mutate(rate_100k = round(n_deaths / pop_ref * 1e5, 2)) |>
  ggplot(aes(x = yr_death, y = cause_grp, fill = rate_100k)) +
  geom_tile(color = "white") +
  geom_text(aes(label = rate_100k), size = 2.3) +
  scale_fill_gradient(low = "#fff5eb", high = "#d94801", name = "Rate") +
  labs(x = "Year", y = "", title = "Mortality rate by cause group") +
  theme_minimal()
