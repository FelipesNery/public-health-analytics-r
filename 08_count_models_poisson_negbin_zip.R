# =============================================================================
# PROJECT 08 - Count Data Models: Poisson, Negative Binomial and
#              Zero-Inflated Models
# =============================================================================
#
# DESCRIPTION
# Modeling counts of cases/notifications as a function of environmental
# (air pollutants, precipitation, temperature), regional and temporal
# variables, using a population offset to estimate Incidence Rate Ratios
# (IRR). Includes assessment of overdispersion and excess zeros.
#
# WHAT WAS DONE
#  - Monthly aggregation of cases by stratum and rescaling of exposures
#  - Poisson regression with log(population) offset -> IRR
#  - Overdispersion test (dispersion ratio and Cameron-Trivedi test)
#  - Negative binomial (MASS::glm.nb) and likelihood ratio test
#  - Exposure x region interaction term and reference category
#  - Zero-Inflated Poisson and Zero-Inflated Negative Binomial models (pscl)
#  - Comparison by log-likelihood / AIC and Vuong test
#  - Diagnostics with simulated residuals (DHARMa)
#
# TECHNIQUES / PACKAGES: stats, MASS, pscl, lmtest, AER, DHARMa, gtsummary
# NOTE: simulated data; variable and object names anonymized.
# =============================================================================

library(tidyverse)
library(MASS)
library(gtsummary)
library(lmtest)

set.seed(6174)

# -----------------------------------------------------------------------------
# 1. Simulated data: monthly cases by region
# -----------------------------------------------------------------------------
reg_nm <- c("R_North", "R_South", "R_East", "R_West", "R_Center")

ct_month <- expand_grid(
  reg_u    = reg_nm,
  year_ref = 2015:2023,
  month_ref = 1:12
) |>
  mutate(
    pop_ref = case_when(reg_u == "R_North" ~ 1.8e6, reg_u == "R_South" ~ 2.9e6,
                        reg_u == "R_East"  ~ 4.1e6, reg_u == "R_West"  ~ 1.2e6,
                        TRUE ~ 2.2e6) * (1 + 0.01 * (year_ref - 2015)),
    pm_exp   = rgamma(n(), shape = 4, rate = 4) * 1e-9,        # original scale
    rain_exp = rgamma(n(), shape = 2, rate = 20),
    temp_exp = 26 + 3 * sin(2 * pi * month_ref / 12) + rnorm(n(), 0, 1),
    season   = factor(if_else(month_ref %in% c(12, 1:5), "Rainy", "Dry")),
    reg_u    = factor(reg_u),
    year_ref = factor(year_ref)
  )

# Rescaling exposures so coefficients are interpretable
ct_month <- ct_month |>
  mutate(pm_sc   = pm_exp * 1e9,
         rain_sc = rain_exp * 1000)

mu_sim <- with(ct_month, exp(-9.5 + 0.18 * pm_sc + 0.002 * rain_sc +
                               0.03 * temp_exp + 0.3 * (season == "Rainy")) * pop_ref)
ct_month$n_cs <- rnbinom(nrow(ct_month), mu = mu_sim, size = 3)    # overdispersed
ct_month <- ct_month |> mutate(log_pop = log(pop_ref),
                               inc_100k = n_cs / pop_ref * 1e5)

# -----------------------------------------------------------------------------
# 2. Poisson regression with offset
# -----------------------------------------------------------------------------
fit_ps <- glm(n_cs ~ pm_sc + rain_sc + temp_exp + season + reg_u + offset(log_pop),
              data = ct_month, family = poisson(link = "log"))
summary(fit_ps)

tbl_regression(fit_ps, exponentiate = TRUE,
               label = list(pm_sc    ~ "Particulate matter (scaled)",
                            rain_sc  ~ "Precipitation (scaled)",
                            temp_exp ~ "Mean temperature (°C)",
                            season   ~ "Season",
                            reg_u    ~ "Region")) |>
  modify_header(estimate ~ "**IRR**")

# -----------------------------------------------------------------------------
# 3. Overdispersion
# -----------------------------------------------------------------------------
# Pearson dispersion ratio (>1 indicates overdispersion)
sum(residuals(fit_ps, type = "pearson")^2) / df.residual(fit_ps)

if (requireNamespace("AER", quietly = TRUE)) AER::dispersiontest(fit_ps)

# -----------------------------------------------------------------------------
# 4. Negative binomial
# -----------------------------------------------------------------------------
fit_nb <- glm.nb(n_cs ~ pm_sc + rain_sc + temp_exp + season + reg_u + offset(log_pop),
                 data = ct_month)
summary(fit_nb)
fit_nb$theta

lrtest(fit_ps, fit_nb)                 # p < 0.05 -> NB better than Poisson
AIC(fit_ps, fit_nb)

tab_nb <- tbl_regression(fit_nb, exponentiate = TRUE) |>
  modify_header(estimate ~ "**IRR**") |> bold_p()
tab_nb

# Exposure x region interaction, with a defined reference region
ct_month$reg_u <- relevel(ct_month$reg_u, ref = "R_North")
fit_nb_int <- glm.nb(n_cs ~ pm_sc * reg_u + year_ref + offset(log_pop), data = ct_month)
tbl_regression(fit_nb_int, exponentiate = TRUE)

# Diagnostics with simulated residuals
if (requireNamespace("DHARMa", quietly = TRUE)) {
  res_sim <- DHARMa::simulateResiduals(fit_nb)
  plot(res_sim)
  DHARMa::testDispersion(res_sim)
  DHARMa::testZeroInflation(res_sim)
}

# -----------------------------------------------------------------------------
# 5. Zero-inflated models
# -----------------------------------------------------------------------------
n_z <- 300
zr_dt <- tibble(
  nx_a = rpois(n_z, 1),
  nx_b = rbinom(n_z, 1, 0.4),
  nx_c = sample(1:4, n_z, TRUE)
) |>
  mutate(
    struct_zero = rbinom(n_z, 1, plogis(0.5 - 0.8 * nx_c + 0.6 * nx_a)),
    cnt_y = if_else(struct_zero == 1, 0L,
                    rnbinom(n_z, mu = exp(-0.3 + 0.5 * nx_c - 0.6 * nx_a + 0.4 * nx_b),
                            size = 1.5))
  )

table(zr_dt$cnt_y)
ggplot(zr_dt, aes(x = cnt_y)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Count", y = "Frequency", title = "Distribution of the count outcome") +
  theme_minimal()

zp_ps <- glm(cnt_y ~ nx_a + nx_b + nx_c, data = zr_dt, family = poisson)
zp_nb <- glm.nb(cnt_y ~ nx_a + nx_b + nx_c, data = zr_dt)

# Point prediction: expected count for a specific profile
predict(zp_ps, newdata = data.frame(nx_a = 1, nx_b = 1, nx_c = 3), type = "response")

if (requireNamespace("pscl", quietly = TRUE)) {
  # part 1 (before |): count component
  # part 2 (after |): logit component for zero inflation
  zp_zip  <- pscl::zeroinfl(cnt_y ~ nx_a + nx_b + nx_c | nx_a + nx_c,
                            data = zr_dt, dist = "poisson")
  zp_zinb <- pscl::zeroinfl(cnt_y ~ nx_a + nx_b + nx_c | nx_a + nx_c,
                            data = zr_dt, dist = "negbin")
  summary(zp_zip)
  summary(zp_zinb)

  pscl::vuong(zp_ps, zp_zip)        # Poisson vs ZIP
  pscl::vuong(zp_nb, zp_zinb)       # NB vs ZINB

  data.frame(
    model  = c("Poisson", "NegBin", "ZIP", "ZINB"),
    logLik = c(logLik(zp_ps), logLik(zp_nb), logLik(zp_zip), logLik(zp_zinb)),
    AIC    = c(AIC(zp_ps), AIC(zp_nb), AIC(zp_zip), AIC(zp_zinb))
  )
}
