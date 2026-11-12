# =============================================================================
# PROJECT 06 - Logistic Regression, ROC Curve and Models for Rare Events
# =============================================================================
#
# DESCRIPTION
# Evaluating the performance of an early warning clinical score to predict
# clinical deterioration, death and ICU admission in hospitalized patients,
# and identifying factors associated with binary outcomes. Includes
# strategies for complete separation / rare events (Firth and bias
# reduction), which are common in small hospital samples.
#
# WHAT WAS DONE
#  - ROC curve, AUC with 95% CI and optimal cut-off (Youden index)
#  - Sensitivity, specificity, PPV and NPV at the cut-off
#  - Univariable (tbl_uvregression) and multivariable logistic regression (OR)
#  - Multicollinearity check (VIF) and goodness of fit
#  - Firth's penalized logistic regression (logistf)
#  - Bias-reduced logistic regression (brglm2::brglmFit)
#  - Univariable + multivariable table side by side (tbl_merge)
#
# TECHNIQUES / PACKAGES: stats, pROC, logistf, brglm2, gtsummary, car, flextable
# NOTE: simulated data; variable and object names anonymized.
# =============================================================================

library(tidyverse)
library(gtsummary)
library(flextable)

set.seed(8128)

# -----------------------------------------------------------------------------
# 1. Simulated data
# -----------------------------------------------------------------------------
n_p <- 260
pts_lg <- tibble(
  warn_score = sample(0:9, n_p, replace = TRUE, prob = c(8, 12, 14, 14, 12, 10, 8, 4, 2, 1)),
  age_p      = round(rnorm(n_p, 55, 16)),
  sex_p      = factor(sample(c("Male", "Female"), n_p, TRUE)),
  ward       = factor(sample(c("General ward", "WardB", "WardC"), n_p, TRUE, prob = c(.6, .25, .15))),
  comorb_p   = factor(sample(c("No", "Yes"), n_p, TRUE, prob = c(.55, .45))),
  stay_d     = round(rexp(n_p, 1 / 9))
) |>
  mutate(
    lp_det    = -3.2 + 0.55 * warn_score + 0.015 * (age_p - 55) + 0.5 * (comorb_p == "Yes"),
    out_det   = rbinom(n_p, 1, plogis(lp_det)),
    out_death = rbinom(n_p, 1, plogis(-5 + 0.6 * warn_score)),     # rare event
    score_bin = factor(if_else(warn_score >= 4, ">= 4", "< 4"), levels = c("< 4", ">= 4"))
  )

table(pts_lg$out_det)
table(pts_lg$out_death)

# -----------------------------------------------------------------------------
# 2. ROC curve of the score
# -----------------------------------------------------------------------------
if (requireNamespace("pROC", quietly = TRUE)) {
  roc_sc <- pROC::roc(pts_lg$out_det, pts_lg$warn_score, ci = TRUE)
  auc_sc <- pROC::auc(roc_sc)
  print(roc_sc$ci)                                  # AUC with 95% CI

  plot(roc_sc, col = "black", lwd = 2, legacy.axes = TRUE,
       main = "ROC curve - score vs. deterioration",
       xlab = "1 - Specificity", ylab = "Sensitivity")
  abline(a = 0, b = 1, lty = 2, col = "gray")
  text(0.35, 0.25, paste("AUC =", round(auc_sc, 3)))

  # Optimal cut-off (Youden)
  cut_yd <- pROC::coords(roc_sc, x = "best", best.method = "youden",
                         transpose = FALSE,
                         ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"))

  df_cut <- data.frame(
    Cutoff      = round(cut_yd$threshold, 2),
    Sensitivity = round(cut_yd$sensitivity, 2),
    Specificity = round(cut_yd$specificity, 2),
    PPV         = round(cut_yd$ppv, 2),
    NPV         = round(cut_yd$npv, 2),
    AUC         = round(as.numeric(auc_sc), 3)
  )

  flextable(df_cut) |>
    set_caption("Performance of the score in predicting deterioration") |>
    autofit() |> align(align = "center", part = "all") |> bold(part = "header")
}

# -----------------------------------------------------------------------------
# 3. Univariable and multivariable logistic regression
# -----------------------------------------------------------------------------
data_lg <- pts_lg |> select(out_det, warn_score, age_p, sex_p, ward, comorb_p, stay_d)

tab_uv <- data_lg |>
  tbl_uvregression(
    method = glm,
    y = out_det,
    method.args = list(family = binomial),
    exponentiate = TRUE
  ) |>
  bold_labels() |> bold_p()

fit_lg <- glm(out_det ~ warn_score + age_p + sex_p + ward + comorb_p + stay_d,
              data = data_lg, family = binomial(link = "logit"))
summary(fit_lg)
exp(cbind(OR = coef(fit_lg), confint.default(fit_lg)))

car::vif(fit_lg)

tab_mv <- tbl_regression(fit_lg, exponentiate = TRUE) |>
  bold_labels() |> bold_p()

tab_final <- tbl_merge(list(tab_uv, tab_mv),
                       tab_spanner = c("**Crude OR**", "**Adjusted OR**"))
tab_final
# tab_final |> as_flex_table() |> save_as_docx(path = "logistic_table.docx")

# McFadden's pseudo-R²
1 - fit_lg$deviance / fit_lg$null.deviance

# -----------------------------------------------------------------------------
# 4. Rare events / separation: Firth and bias reduction
# -----------------------------------------------------------------------------
if (requireNamespace("logistf", quietly = TRUE)) {
  fit_firth <- logistf::logistf(out_death ~ warn_score + age_p + sex_p + comorb_p, data = pts_lg)
  summary(fit_firth)
  exp(cbind(OR = coef(fit_firth), confint(fit_firth)))

  fit_firth_bin <- logistf::logistf(out_det ~ score_bin, data = pts_lg)
  summary(fit_firth_bin)
}

if (requireNamespace("brglm2", quietly = TRUE)) {
  fit_br <- glm(out_death ~ warn_score + sex_p + comorb_p + ward,
                data = pts_lg, family = binomial(), method = brglm2::brglmFit)
  summary(fit_br)
  tbl_regression(fit_br, exponentiate = TRUE)
}
