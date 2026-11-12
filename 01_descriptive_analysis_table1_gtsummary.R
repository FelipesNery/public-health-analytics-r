# =============================================================================
# PROJECT 01 - Descriptive Analysis and Stratified "Table 1" with gtsummary
# =============================================================================
#
# DESCRIPTION
# Building publication-ready descriptive tables ("Table 1") that compare
# sociodemographic and clinical characteristics across outcome groups. Used in
# several epidemiological studies (surveillance cohorts, hospital-based
# studies and field surveys).
#
# WHAT WAS DONE
#  - Recoding numeric codes into labelled factors
#  - Overall descriptive table and table stratified by outcome
#  - Median (IQR) for continuous and n / N (%) for categorical variables
#  - Fisher's exact test (categorical) and Wilcoxon test (continuous) via add_p()
#  - Removal of empty rows, header formatting and export to Word
#
# TECHNIQUES / PACKAGES: tidyverse, gtsummary, flextable, officer
# NOTE: simulated data; variable and object names anonymized.
# =============================================================================

library(tidyverse)
library(gtsummary)
library(flextable)

set.seed(2718)

# -----------------------------------------------------------------------------
# 1. Simulated data
# -----------------------------------------------------------------------------
n_obs <- 420

cohort_kx <- tibble(
  out_grp    = sample(c(0, 1), n_obs, replace = TRUE, prob = c(0.72, 0.28)),
  age_yr     = round(rnorm(n_obs, mean = 41, sd = 14)),
  sex_cod    = sample(c(1, 2), n_obs, replace = TRUE),
  area_cod   = sample(c(1, 2, 9), n_obs, replace = TRUE, prob = c(0.6, 0.35, 0.05)),
  edu_cod    = sample(1:4, n_obs, replace = TRUE),
  comorb_n   = rpois(n_obs, lambda = 0.8),
  stay_days  = round(rexp(n_obs, rate = 1 / 7))
)

# -----------------------------------------------------------------------------
# 2. Recoding (numeric codes -> labelled factors)
# -----------------------------------------------------------------------------
cohort_kx <- cohort_kx |>
  mutate(
    out_grp  = factor(out_grp, levels = c(0, 1), labels = c("No", "Yes")),
    sex_cod  = factor(sex_cod, levels = c(1, 2), labels = c("Male", "Female")),
    area_cod = na_if(area_cod, 9),                       # 9 = unknown -> NA
    area_cod = factor(area_cod, levels = c(1, 2), labels = c("Urban", "Rural")),
    edu_cod  = factor(edu_cod, levels = 1:4,
                      labels = c("No schooling", "Elementary",
                                 "High school", "College")),
    age_grp  = cut(age_yr, breaks = c(-Inf, 19, 39, 59, Inf),
                   labels = c("<20", "20-39", "40-59", "60+")),
    comorb_bin = factor(if_else(comorb_n > 0, "Yes", "No"))
  )

# -----------------------------------------------------------------------------
# 3. Variable labels
# -----------------------------------------------------------------------------
var_labels <- list(
  age_yr     ~ "Age (years)",
  age_grp    ~ "Age group (years)",
  sex_cod    ~ "Sex",
  area_cod   ~ "Area of residence",
  edu_cod    ~ "Education",
  comorb_bin ~ "Comorbidity",
  stay_days  ~ "Length of stay (days)"
)

# -----------------------------------------------------------------------------
# 4. Overall descriptive table (no stratification)
# -----------------------------------------------------------------------------
tab_overall <- cohort_kx |>
  select(-out_grp, -comorb_n) |>
  tbl_summary(missing = "no", label = var_labels) |>
  # drop category rows with no observations
  modify_table_body(~ .x |> filter(!(stat_0 %in% c("0 (0%)")))) |>
  bold_labels()

tab_overall

# -----------------------------------------------------------------------------
# 5. Table 1 stratified by outcome with hypothesis tests
# -----------------------------------------------------------------------------
tab_strat <- cohort_kx |>
  select(out_grp, age_yr, age_grp, sex_cod, area_cod, edu_cod,
         comorb_bin, stay_days) |>
  tbl_summary(
    by = out_grp,
    missing_text = "Missing",
    label = var_labels,
    statistic = list(
      all_continuous()  ~ "{median} ({p25}; {p75})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    digits = all_continuous() ~ 1
  ) |>
  add_overall(col_label = "**Total**  \nN = {style_number(N)}") |>
  add_p(
    test = list(all_categorical() ~ "fisher.test",
                all_continuous()  ~ "wilcox.test"),
    test.args = all_categorical() ~ list(workspace = 2e7, simulate.p.value = TRUE),
    pvalue_fun = label_style_pvalue(digits = 2)
  ) |>
  modify_header(label ~ "**Variables**", p.value ~ "**p-value**") |>
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Outcome**") |>
  modify_caption("**Table 1. Characteristics of the study population by outcome**") |>
  bold_labels() |>
  bold_p(t = 0.05)

tab_strat

# -----------------------------------------------------------------------------
# 6. Export to Word
# -----------------------------------------------------------------------------
tab_strat |>
  as_flex_table() |>
  autofit() |>
  save_as_docx(path = "table1_stratified.docx")
