# =============================================================================
# PROJECT 05 - Simple and Multiple Linear Regression, Box-Cox and Stepwise
# =============================================================================
#
# DESCRIPTION
# Modeling a continuous outcome (e.g., emergency service response time,
# newborn length, sales volume) from numeric and categorical predictors, with
# residual diagnostics, Box-Cox transformation, variable selection and
# predictive validation.
#
# WHAT WAS DONE
#  - Univariable linear regressions combined into a single gtsummary table
#  - Multiple linear regression with dummy variables
#  - Diagnostics: normality of residuals (Shapiro-Francia), VIF, plots
#  - Box-Cox transformation (car::powerTransform) and R² comparison
#  - Stepwise selection with k = qchisq(0.95, 1) (equivalent to p < 0.05)
#  - Interaction term and prediction on a test set (RMSE, R², MAE)
#
# TECHNIQUES / PACKAGES: stats, car, nortest, gtsummary, fastDummies, caret
# NOTE: simulated data; variable and object names anonymized.
# =============================================================================

library(tidyverse)
library(gtsummary)
library(car)

set.seed(1414)

# -----------------------------------------------------------------------------
# 1. Simulated data
# -----------------------------------------------------------------------------
n_l <- 350
db_lin <- tibble(
  site_cat = factor(sample(c("Site1", "Site2", "Site3"), n_l, TRUE)),
  age_grp  = factor(sample(c("0-19", "20-59", "60+"), n_l, TRUE)),
  sex_f    = factor(sample(c("F", "M"), n_l, TRUE)),
  season   = factor(sample(c("Normal", "Critical"), n_l, TRUE)),
  dist_km  = runif(n_l, 2, 60)
) |>
  mutate(
    resp_y = 12 + 0.55 * dist_km +
      if_else(season == "Critical", 8, 0) +
      if_else(site_cat == "Site3", 5, 0) +
      rnorm(n_l, 0, 6),
    resp_y = pmax(resp_y, 1)
  )

# -----------------------------------------------------------------------------
# 2. Univariable regressions (one per predictor) in a single table
# -----------------------------------------------------------------------------
tab_uni <- db_lin |>
  tbl_uvregression(
    method = lm,
    y = resp_y,
    include = c(site_cat, age_grp, sex_f, season, dist_km)
  ) |>
  bold_labels() |>
  bold_p()
tab_uni

# -----------------------------------------------------------------------------
# 3. Multiple linear regression
# -----------------------------------------------------------------------------
fit_mult <- lm(resp_y ~ site_cat + age_grp + sex_f + season + dist_km, data = db_lin)
summary(fit_mult)

tab_mult <- tbl_regression(fit_mult,
                           label = list(site_cat ~ "Site",
                                        age_grp ~ "Age group",
                                        sex_f ~ "Sex",
                                        season ~ "Period",
                                        dist_km ~ "Distance (km)")) |>
  bold_labels() |> bold_p()

tbl_merge(list(tab_uni, tab_mult),
          tab_spanner = c("**Univariable**", "**Multivariable**"))

# Diagnostics
vif(fit_mult)                         # multicollinearity (VIF > 5 = warning)
par(mfrow = c(2, 2)); plot(fit_mult); par(mfrow = c(1, 1))
if (requireNamespace("nortest", quietly = TRUE)) nortest::sf.test(residuals(fit_mult))

# -----------------------------------------------------------------------------
# 4. Box-Cox transformation
# -----------------------------------------------------------------------------
fit_simple <- lm(resp_y ~ dist_km, data = db_lin)
shapiro.test(residuals(fit_simple))

lmb_bc <- powerTransform(db_lin$resp_y)
lmb_bc$lambda

db_lin <- db_lin |>
  mutate(resp_bc = (resp_y^lmb_bc$lambda - 1) / lmb_bc$lambda)

fit_bc <- lm(resp_bc ~ dist_km, data = db_lin)
summary(fit_bc)
shapiro.test(residuals(fit_bc))

tibble(R2_OLS = summary(fit_simple)$r.squared,
       R2_BoxCox = summary(fit_bc)$r.squared) |>
  round(4)

# Histogram of residuals with theoretical normal curve
tibble(res_bc = residuals(fit_bc)) |>
  ggplot(aes(x = res_bc)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 color = "grey50", fill = "grey90") +
  stat_function(fun = dnorm,
                args = list(mean = mean(residuals(fit_bc)), sd = sd(residuals(fit_bc))),
                color = "#440154", linewidth = 1.2) +
  labs(x = "Residuals (Box-Cox)", y = "Density") +
  theme_minimal()

# Prediction back-transformed to the original scale
pred_bc <- predict(fit_bc, newdata = db_lin)
db_lin$pred_orig <- (pred_bc * lmb_bc$lambda + 1)^(1 / lmb_bc$lambda)

# -----------------------------------------------------------------------------
# 5. Dummies + Stepwise + validation on a test set
# -----------------------------------------------------------------------------
idx_tr <- sample(seq_len(n_l), size = round(0.7 * n_l))
db_train <- db_lin[idx_tr, ]
db_test  <- db_lin[-idx_tr, ]

if (requireNamespace("fastDummies", quietly = TRUE)) {
  db_dummy <- fastDummies::dummy_columns(
    db_train |> select(resp_y, dist_km, site_cat, age_grp, sex_f, season),
    select_columns = c("site_cat", "age_grp", "sex_f", "season"),
    remove_selected_columns = TRUE,
    remove_first_dummy = TRUE
  )
  head(db_dummy)
}

fit_full <- lm(resp_y ~ site_cat + age_grp + sex_f + season * dist_km, data = db_train)
fit_step <- step(fit_full, k = qchisq(p = 0.05, df = 1, lower.tail = FALSE), trace = 0)
summary(fit_step)

pred_test <- predict(fit_step, newdata = db_test)
df_val    <- data.frame(obs = db_test$resp_y, pred = pred_test)

# Performance metrics
tibble(
  RMSE = sqrt(mean((df_val$obs - df_val$pred)^2)),
  R2   = cor(df_val$obs, df_val$pred)^2,
  MAE  = mean(abs(df_val$obs - df_val$pred))
)
# equivalent: caret::defaultSummary(df_val)

ggplot(df_val, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
  labs(x = "Observed", y = "Predicted", title = "Observed vs. Predicted (test set)") +
  theme_bw()
