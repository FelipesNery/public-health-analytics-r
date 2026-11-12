# =============================================================================
# PROJECT 03 - Group Comparisons: t-test, Wilcoxon, ANOVA and Tukey
# =============================================================================
#
# DESCRIPTION
# Comparing continuous variables between two or more groups, choosing the test
# according to its assumptions (normality and homogeneity of variances),
# computing effect sizes and drawing boxplots annotated with p-values.
# Applied, for example, to compare emergency service response times across
# seasonal periods and anthropometric measures of newborns across clinical
# groups.
#
# WHAT WAS DONE
#  - Descriptive statistics by group (n, mean, median, SD, quartiles)
#  - Shapiro-Wilk (normality) and Levene (homoscedasticity) tests
#  - Student's / Welch's t-test for 2 groups
#  - Mann-Whitney (Wilcoxon) test + effect size r (rstatix/coin)
#  - One-way ANOVA, critical F and Tukey's multiple comparisons (HSD)
#  - Kruskal-Wallis + Dunn's test as the non-parametric alternative
#  - Boxplots with p-values (ggpubr)
#
# TECHNIQUES / PACKAGES: stats, rstatix, coin, car, ggpubr, tidyverse
# NOTE: simulated data; variable and object names anonymized.
# =============================================================================

library(tidyverse)
library(rstatix)
library(ggpubr)
library(car)

set.seed(3141)

# -----------------------------------------------------------------------------
# 1. Simulated data
# -----------------------------------------------------------------------------
sample_zt <- tibble(
  period_cat = factor(sample(c("Normal", "Critical"), 300, TRUE, prob = c(0.6, 0.4))),
  clin_grp   = factor(sample(c("Control", "GroupA", "GroupB", "GroupA+B"),
                             300, TRUE)),
  time_min   = rgamma(300, shape = 3, rate = 0.08),
  measure_kg = rnorm(300, 3.1, 0.45)
) |>
  mutate(
    time_min   = time_min + if_else(period_cat == "Critical", 6, 0),
    measure_kg = measure_kg - case_when(clin_grp == "GroupA"   ~ 0.40,
                                        clin_grp == "GroupA+B" ~ 0.55,
                                        TRUE ~ 0)
  )

# -----------------------------------------------------------------------------
# 2. Descriptive statistics by group
# -----------------------------------------------------------------------------
desc_grp <- sample_zt |>
  group_by(clin_grp) |>
  summarise(
    n_obs  = n(),
    mean   = mean(measure_kg),
    median = median(measure_kg),
    sd     = sd(measure_kg),
    q1     = quantile(measure_kg, 0.25, type = 5),
    q3     = quantile(measure_kg, 0.75, type = 5),
    min    = min(measure_kg),
    max    = max(measure_kg)
  )
desc_grp

# -----------------------------------------------------------------------------
# 3. Assumptions
# -----------------------------------------------------------------------------
sample_zt |> group_by(period_cat) |> shapiro_test(time_min)
sample_zt |> group_by(clin_grp)   |> shapiro_test(measure_kg)

leveneTest(time_min ~ period_cat, data = sample_zt)
leveneTest(measure_kg ~ clin_grp, data = sample_zt)

# -----------------------------------------------------------------------------
# 4. Two groups - t-test and Mann-Whitney
# -----------------------------------------------------------------------------
t.test(time_min ~ period_cat, data = sample_zt, var.equal = FALSE)  # Welch

mw_test <- sample_zt |>
  wilcox_test(time_min ~ period_cat) |>
  add_significance()
mw_test

# effect size r (small ~0.1, medium ~0.3, large ~0.5)
mw_eff <- sample_zt |> wilcox_effsize(time_min ~ period_cat)
mw_eff

# Annotated boxplot
mw_test <- mw_test |> add_xy_position(x = "period_cat")

ggboxplot(sample_zt, x = "period_cat", y = "time_min", fill = "period_cat") +
  stat_pvalue_manual(mw_test, label = "p = {p}") +
  annotate("text", x = 1.5, y = max(sample_zt$time_min) * 0.9,
           label = paste0("r = ", round(mw_eff$effsize, 3))) +
  labs(x = "Period", y = "Time (min)",
       title = "Distribution of time by period") +
  theme(legend.position = "none")

# -----------------------------------------------------------------------------
# 5. Three or more groups - ANOVA + Tukey
# -----------------------------------------------------------------------------
# H0: all means are equal | H1: at least one mean differs
fit_aov <- aov(measure_kg ~ clin_grp, data = sample_zt)
summary(fit_aov)

# Critical F (df1 = k - 1 ; df2 = N - k)
k_grp <- nlevels(sample_zt$clin_grp)
qf(0.95, df1 = k_grp - 1, df2 = nrow(sample_zt) - k_grp)

# Normality of residuals
shapiro.test(residuals(fit_aov))

# Multiple comparisons
tk_res <- TukeyHSD(fit_aov, conf.level = 0.95)
tk_res

par(mar = c(5, 9, 4, 2))
plot(tk_res, las = 1, cex.axis = 0.7)

# -----------------------------------------------------------------------------
# 6. Non-parametric alternative - Kruskal-Wallis + Dunn
# -----------------------------------------------------------------------------
sample_zt |> kruskal_test(measure_kg ~ clin_grp)
sample_zt |> kruskal_effsize(measure_kg ~ clin_grp)
post_dunn <- sample_zt |> dunn_test(measure_kg ~ clin_grp, p.adjust.method = "bonferroni")
post_dunn

post_dunn <- post_dunn |> add_xy_position(x = "clin_grp")
ggboxplot(sample_zt, x = "clin_grp", y = "measure_kg", fill = "clin_grp") +
  stat_pvalue_manual(post_dunn, hide.ns = TRUE) +
  labs(x = "Group", y = "Measure (kg)") +
  theme(legend.position = "none")
