# =============================================================================
# PROJECT 02 - Contingency Tables, Association Tests and Risk Measures
# =============================================================================
#
# DESCRIPTION
# Assessing the association between exposure and outcome in 2x2 and RxC
# tables, including paired data, and computing the classic epidemiological
# effect measures (Risk Ratio and Odds Ratio) with 95% CI.
#
# WHAT WAS DONE
#  - Building 2x2 tables from aggregated counts and from long-format data
#  - Chi-square test of independence (checking expected counts < 5)
#  - Fisher's exact test as an alternative for small samples
#  - McNemar's test for paired data (manual calculation + built-in function)
#  - Risk Ratio (RR) and Odds Ratio (OR) with 95% CI (epitools)
#  - Mosaic plot and stacked bar chart to visualize the association
#
# TECHNIQUES / PACKAGES: stats, epitools, tidyverse
# NOTE: simulated data; variable and object names anonymized.
# =============================================================================

library(tidyverse)
library(epitools)

set.seed(1618)

# -----------------------------------------------------------------------------
# 1. 2x2 table from aggregated counts
# -----------------------------------------------------------------------------
mtx_ab <- matrix(c(21, 139, 204, 451), nrow = 2, byrow = TRUE)
rownames(mtx_ab) <- c("Exposed", "Unexposed")
colnames(mtx_ab) <- c("Outcome +", "Outcome -")
mtx_ab

# Chi-square test of independence
chi_res <- chisq.test(mtx_ab)
chi_res
chi_res$expected            # check assumption (expected counts >= 5)

# Fisher's exact test (alternative when expected counts < 5)
fisher.test(mtx_ab)

# -----------------------------------------------------------------------------
# 2. Table from individual-level data (long format)
# -----------------------------------------------------------------------------
rec_long <- tibble(
  exp_fct = sample(c("exposed", "unexposed"), 600, TRUE, prob = c(0.3, 0.7))
) |>
  mutate(out_fct = if_else(
    runif(n()) < if_else(exp_fct == "exposed", 0.25, 0.12),
    "positive", "negative")
  ) |>
  mutate(
    exp_num = case_when(exp_fct == "exposed" ~ 1, TRUE ~ 0),
    out_num = case_when(out_fct == "positive" ~ 1, TRUE ~ 0)
  )

tab_cross <- table(rec_long$exp_fct, rec_long$out_fct)
tab_cross
prop.table(tab_cross, margin = 1) * 100      # row percentages
chisq.test(tab_cross, correct = TRUE)        # with Yates' correction

# -----------------------------------------------------------------------------
# 3. Paired data - McNemar's test
# -----------------------------------------------------------------------------
# Discordant pairs: b and c
pair_b <- 41
pair_c <- 18

# Manual statistic (with continuity correction)
mcn_manual <- (abs(pair_b - pair_c) - 1)^2 / (pair_b + pair_c)
mcn_manual
pchisq(mcn_manual, df = 1, lower.tail = FALSE)

mtx_pair <- matrix(c(11, pair_b, pair_c, 77), nrow = 2, byrow = TRUE,
                   dimnames = list(Case = c("Exp+", "Exp-"),
                                   Control = c("Exp+", "Exp-")))
mcnemar.test(mtx_pair)

# -----------------------------------------------------------------------------
# 4. Effect measures: Risk Ratio and Odds Ratio
# -----------------------------------------------------------------------------
# epitools expects: rows = exposure (reference in 1st row), columns = outcome
# (reference in 1st column). rev = "both" flips the table into that layout.
rr_res <- riskratio(mtx_ab, method = "wald", rev = "both", conf.level = 0.95)
rr_res$measure

or_res <- oddsratio(mtx_ab, method = "wald", rev = "both")
or_res$measure

# From individual-level vectors
riskratio(x = rec_long$exp_num, y = rec_long$out_num, method = "wald")$measure
oddsratio(x = rec_long$exp_num, y = rec_long$out_num, method = "wald")$measure

# -----------------------------------------------------------------------------
# 5. Visualization
# -----------------------------------------------------------------------------
mosaicplot(t(mtx_ab), color = c("#2b8cbe", "#a6bddb"),
           main = "Association between exposure and outcome")

rec_long |>
  count(exp_fct, out_fct) |>
  group_by(exp_fct) |>
  mutate(pct = n / sum(n) * 100) |>
  ggplot(aes(x = exp_fct, y = pct, fill = out_fct)) +
  geom_col(position = "stack") +
  geom_text(aes(label = paste0(round(pct, 1), "%")),
            position = position_stack(vjust = 0.5)) +
  labs(x = "Exposure", y = "%", fill = "Outcome") +
  theme_classic()
