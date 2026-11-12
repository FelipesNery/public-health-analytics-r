# =============================================================================
# PROJECT 04 - Correlation Analysis between Continuous Variables
# =============================================================================
#
# DESCRIPTION
# Investigating linear and monotonic relationships between continuous
# variables (e.g., anthropometric measures, air pollutants and climate
# variables), with significance tests and several ways of visualizing the
# correlation matrix.
#
# WHAT WAS DONE
#  - Pearson and Spearman correlation matrices
#  - Significance test (cor.test) for pairs of variables
#  - Matrix of p-values with rstatix
#  - Correlation maps (corrplot) and scatterplot matrix (GGally)
#  - Scatterplot with regression line, R² and p-value (ggpubr)
#
# TECHNIQUES / PACKAGES: stats, rstatix, corrplot, GGally, ggpubr, tidyverse
# NOTE: simulated data; variable and object names anonymized.
# =============================================================================

library(tidyverse)
library(rstatix)
library(corrplot)
library(ggpubr)

set.seed(5772)

# -----------------------------------------------------------------------------
# 1. Simulated data
# -----------------------------------------------------------------------------
n_rec <- 150
data_cr <- tibble(
  var_wk  = round(runif(n_rec, 30, 42), 1),
  grp_bin = factor(if_else(var_wk < 37, "Yes", "No"))
) |>
  mutate(
    var_m = 0.18 * var_wk - 3.6 + rnorm(n_rec, 0, 0.35),
    var_l = 1.3  * var_wk + 0.5 + rnorm(n_rec, 0, 2),
    var_h = 0.85 * var_wk + 3   + rnorm(n_rec, 0, 1.5)
  )

num_vars <- data_cr |> select(var_wk, var_m, var_l, var_h)

# -----------------------------------------------------------------------------
# 2. Correlation matrices
# -----------------------------------------------------------------------------
cm_pearson  <- cor(num_vars, method = "pearson")
cm_spearman <- cor(num_vars, method = "spearman")
round(cm_pearson, 2)
round(cm_spearman, 2)

# Test for a single pair
cor.test(data_cr$var_wk, data_cr$var_m, method = "pearson")
cor.test(data_cr$var_wk, data_cr$var_m, method = "spearman", exact = FALSE)

# All pairs with p-values (tidy format)
num_vars |> cor_test(method = "pearson")
num_vars |> cor_mat() |> cor_get_pval()

# -----------------------------------------------------------------------------
# 3. Visualizations
# -----------------------------------------------------------------------------
corrplot(cm_pearson)
corrplot.mixed(cm_pearson, lower = "circle", upper = "number", tl.col = "black")

if (requireNamespace("GGally", quietly = TRUE)) {
  GGally::ggpairs(data_cr, columns = 1:4,
                  lower = list(continuous = "smooth"))
  GGally::ggpairs(data_cr, columns = 1:4, aes(colour = grp_bin))
}

# Scatterplot with regression equation, R² and p-value
ggplot(data_cr, aes(x = var_wk, y = var_m)) +
  geom_point(size = 2.5, pch = 21, col = "brown4", fill = "coral") +
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE) +
  stat_regline_equation(label.y.npc = 0.95) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~")),
           label.y.npc = 0.85) +
  labs(x = "Variable X", y = "Variable Y") +
  theme_light()
