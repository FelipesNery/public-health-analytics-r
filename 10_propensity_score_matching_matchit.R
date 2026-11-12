# =============================================================================
# PROJECT 10 - Propensity Score Matching
# =============================================================================
#
# DESCRIPTION
# Comparing a continuous outcome (e.g., emergency service response time)
# between two periods/conditions (e.g., drought vs. normal river levels),
# controlling for confounding by individual characteristics through
# propensity score matching, with exact matching on sex.
#
# WHAT WAS DONE
#  - Propensity score estimation (logit) with MatchIt
#  - 1:1 nearest-neighbor matching with exact matching on one covariate and
#    a caliper
#  - Balance assessment (standardized mean difference before/after)
#  - Comparison table in the matched sample (gtsummary)
#  - Wilcoxon test and effect size in the matched sample
#  - Estimation of the average effect with regression weighted by the
#    matching weights and pair-clustered robust standard errors
#
# TECHNIQUES / PACKAGES: MatchIt, gtsummary, rstatix, sandwich, lmtest
# NOTE: simulated data; variable and object names anonymized.
# =============================================================================

library(tidyverse)
library(gtsummary)
library(rstatix)

set.seed(7351)

# -----------------------------------------------------------------------------
# 1. Simulated data (treatment depends on covariates -> confounding)
# -----------------------------------------------------------------------------
n_m <- 800
obs_pm <- tibble(
  sex_pm  = factor(sample(c("F", "M"), n_m, TRUE)),
  age_pm  = round(runif(n_m, 1, 90)),
  dist_pm = rgamma(n_m, shape = 2, rate = 0.08)
) |>
  mutate(
    trt_pm = rbinom(n_m, 1, plogis(-1 + 0.03 * dist_pm + 0.01 * (age_pm - 45))),
    trt_pm = factor(trt_pm, levels = c(0, 1), labels = c("Normal", "Critical")),
    out_pm = 15 + 0.6 * dist_pm + 4 * (trt_pm == "Critical") + rnorm(n_m, 0, 6),
    agegrp_pm = cut(age_pm, c(0, 17, 59, Inf), labels = c("0-17", "18-59", "60+"))
  )

# Naive comparison (unadjusted)
obs_pm |> group_by(trt_pm) |> summarise(mean = mean(out_pm), n = n())

# -----------------------------------------------------------------------------
# 2. Matching
# -----------------------------------------------------------------------------
if (requireNamespace("MatchIt", quietly = TRUE)) {
  library(MatchIt)

  pm_obj <- matchit(trt_pm ~ sex_pm + age_pm + dist_pm,
                    data = obs_pm,
                    method = "nearest",
                    distance = "glm",
                    ratio = 1,
                    caliper = 0.2,
                    exact = ~ sex_pm)

  pm_obj
  summary(pm_obj, un = TRUE)             # balance before x after
  plot(summary(pm_obj))                  # love plot (SMD)
  plot(pm_obj, type = "jitter", interactive = FALSE)

  pairs_pm <- match.data(pm_obj)
  nrow(pairs_pm)

  # ---------------------------------------------------------------------------
  # 3. Table in the matched sample
  # ---------------------------------------------------------------------------
  pairs_pm |>
    select(trt_pm, sex_pm, agegrp_pm, dist_pm, out_pm) |>
    tbl_summary(by = trt_pm, missing = "no",
                label = list(sex_pm ~ "Sex", agegrp_pm ~ "Age group",
                             dist_pm ~ "Distance", out_pm ~ "Outcome (min)")) |>
    add_overall(col_label = "**Total**  \nN = {style_number(N)}") |>
    add_p() |>
    modify_header(label = "**Variables**", p.value = "**p-value**") |>
    modify_spanning_header(all_stat_cols() ~ "**Condition**") |>
    bold_labels()

  # ---------------------------------------------------------------------------
  # 4. Test in the matched sample
  # ---------------------------------------------------------------------------
  pairs_pm |> wilcox_test(out_pm ~ trt_pm)
  pairs_pm |> wilcox_effsize(out_pm ~ trt_pm)

  # ---------------------------------------------------------------------------
  # 5. Average treatment effect on the treated (ATT) with pair-clustered SE
  # ---------------------------------------------------------------------------
  fit_att <- lm(out_pm ~ trt_pm + sex_pm + age_pm + dist_pm,
                data = pairs_pm, weights = weights)
  if (requireNamespace("sandwich", quietly = TRUE)) {
    lmtest::coeftest(fit_att, vcov. = sandwich::vcovCL, cluster = ~ subclass)
  } else {
    summary(fit_att)
  }
}
