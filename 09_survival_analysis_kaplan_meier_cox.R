# =============================================================================
# PROJECT 09 - Survival Analysis: Kaplan-Meier, Log-Rank and Cox
# =============================================================================
#
# DESCRIPTION
# Evaluating time to recurrence of an infectious disease in individuals with
# and without a coinfection, using surveillance databases previously joined
# by record linkage. Follow-up time was analyzed in 90-, 180- and 365-day
# windows with administrative censoring.
# (Master's thesis project - epidemiology of infectious diseases.)
#
# WHAT WAS DONE
#  - Calculation of time between notification dates and creation of follow-up
#    windows censored at the end of each window
#  - Descriptive tables by outcome for each window (Fisher's test)
#  - Kaplan-Meier curves with 95% CI, p-value and risk table (survminer)
#  - Log-rank test (survdiff)
#  - Multivariable Cox regression -> Hazard Ratios with 95% CI
#  - Check of the proportional hazards assumption (Schoenfeld residuals -
#    cox.zph) and formatted table with gt
#  - Subgroup sensitivity analyses
#
# TECHNIQUES / PACKAGES: survival, survminer, gtsummary, broom, gt, gridExtra
# NOTE: simulated data; variable and object names anonymized.
# =============================================================================

library(tidyverse)
library(survival)
library(survminer)
library(gtsummary)
library(broom)
library(gt)

set.seed(4669)

# -----------------------------------------------------------------------------
# 1. Simulated data
# -----------------------------------------------------------------------------
n_s <- 900
cohort_sv <- tibble(
  id_pers  = sprintf("P%04d", 1:n_s),
  coinf    = factor(sample(c("No", "Yes"), n_s, TRUE, prob = c(0.8, 0.2))),
  sex_s    = factor(sample(c("Female", "Male"), n_s, TRUE)),
  age_s    = round(runif(n_s, 5, 80)),
  area_s   = factor(sample(c("Urban", "Rural"), n_s, TRUE, prob = c(0.55, 0.45))),
  dt_start = as.Date("2012-01-01") + sample(0:2500, n_s, TRUE)
) |>
  mutate(
    hz_ind  = 0.0025 * exp(0.6 * (coinf == "Yes") + 0.25 * (area_s == "Rural")),
    t_evt   = round(rexp(n_s, rate = hz_ind)),
    t_cens  = round(runif(n_s, 30, 900)),
    dt_end  = dt_start + pmin(t_evt, t_cens),
    evt_rec = as.integer(t_evt <= t_cens)
  )

# -----------------------------------------------------------------------------
# 2. Follow-up time and windows with administrative censoring
# -----------------------------------------------------------------------------
cohort_sv <- cohort_sv |>
  mutate(t_days = as.numeric(difftime(dt_end, dt_start, units = "days")))

make_window <- function(df, time, event, w) {
  df |>
    mutate(
      "t_{w}"   := pmin({{ time }}, w),
      "evt_{w}" := if_else({{ event }} == 1 & {{ time }} <= w, 1L, 0L)
    )
}

cohort_sv <- cohort_sv |>
  make_window(t_days, evt_rec, 90) |>
  make_window(t_days, evt_rec, 180) |>
  make_window(t_days, evt_rec, 365)

# -----------------------------------------------------------------------------
# 3. Descriptive table by outcome (365-day window)
# -----------------------------------------------------------------------------
lbl_sv <- list(age_s ~ "Age (years)", sex_s ~ "Sex",
               area_s ~ "Area", coinf ~ "Coinfection")

cohort_sv |>
  mutate(evt_365 = factor(evt_365, labels = c("No", "Yes"))) |>
  select(evt_365, age_s, sex_s, area_s, coinf) |>
  tbl_summary(by = evt_365, label = lbl_sv,
              statistic = all_categorical() ~ "{n} / {N} ({p}%)") |>
  add_overall(col_label = "**Total**") |>
  add_p(test = all_categorical() ~ "fisher.test",
        test.args = all_categorical() ~ list(workspace = 2e7)) |>
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Recurrence within 365 days**") |>
  bold_labels() |> bold_p()

# -----------------------------------------------------------------------------
# 4. Kaplan-Meier + Log-Rank
# -----------------------------------------------------------------------------
km_365 <- survfit(Surv(t_365, evt_365) ~ coinf, data = cohort_sv)
summary(km_365, times = c(90, 180, 365))

survdiff(Surv(t_365, evt_365) ~ coinf, data = cohort_sv)        # Log-Rank

g_km <- ggsurvplot(
  km_365, data = cohort_sv,
  conf.int = TRUE, pval = TRUE, pval.coord = c(10, 0.6),
  break.time.by = 30,
  xlab = "Time (days)", ylab = "Recurrence-free probability",
  legend.title = "Coinfection", legend.labs = c("No", "Yes"),
  risk.table = "abs_pct", risk.table.title = "At risk - n (%)",
  risk.table.fontsize = 3, risk.table.y.text = FALSE,
  ggtheme = theme_bw(), ylim = c(0.5, 1),
  title = "Time to recurrence within 365 days"
)
g_km

g_comb <- gridExtra::arrangeGrob(g_km$plot, g_km$table, ncol = 1, heights = c(2, 1))
ggsave("km_recurrence_365.png", plot = g_comb, dpi = 300,
       width = 21, height = 20, units = "cm")

# Curves for the three windows in sequence
for (w in c(90, 180, 365)) {
  frm_w <- as.formula(paste0("Surv(t_", w, ", evt_", w, ") ~ coinf"))
  km_w  <- surv_fit(frm_w, data = cohort_sv)
  print(ggsurvplot(km_w, data = cohort_sv, pval = TRUE, conf.int = TRUE,
                   title = paste(w, "-day window")))
}

# -----------------------------------------------------------------------------
# 5. Cox regression
# -----------------------------------------------------------------------------
cx_365 <- coxph(Surv(t_365, evt_365) ~ coinf + sex_s + age_s + area_s, data = cohort_sv)
summary(cx_365)

tidy(cx_365, exponentiate = TRUE, conf.int = TRUE) |>
  transmute(Variable = term,
            HR = round(estimate, 2),
            `95% CI` = paste0(round(conf.low, 2), " - ", round(conf.high, 2)),
            `p-value` = signif(p.value, 3)) |>
  gt() |>
  tab_header(title = "Cox regression", subtitle = "Recurrence within 365 days")

tbl_regression(cx_365, exponentiate = TRUE) |> bold_p()

forest_cx <- ggforest(cx_365, data = as.data.frame(cohort_sv))
forest_cx

# -----------------------------------------------------------------------------
# 6. Proportional hazards assumption (Schoenfeld)
# -----------------------------------------------------------------------------
zph_365 <- cox.zph(cx_365)
zph_365

as.data.frame(zph_365$table) |>
  rownames_to_column("Variable") |>
  rename(Chi_square = chisq, p_value = p) |>
  mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value, 4)))) |>
  gt() |>
  tab_header(title = "Proportional hazards test (Schoenfeld residuals)")

ggcoxzph(zph_365)

# -----------------------------------------------------------------------------
# 7. Subgroup sensitivity analysis (e.g., rural area only)
# -----------------------------------------------------------------------------
cx_sub <- coxph(Surv(t_365, evt_365) ~ coinf + sex_s + age_s,
                data = filter(cohort_sv, area_s == "Rural"))
tbl_merge(
  list(tbl_regression(cx_365, exponentiate = TRUE),
       tbl_regression(cx_sub, exponentiate = TRUE)),
  tab_spanner = c("**Overall**", "**Rural area**")
)
