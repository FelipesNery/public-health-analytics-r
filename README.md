# Statistical and Epidemiological Analysis in R

![R](https://img.shields.io/badge/R-4.6-276DC3?logo=r&logoColor=white)
![Field](https://img.shields.io/badge/field-epidemiology%20%7C%20public%20health-2e8b57)
![Data](https://img.shields.io/badge/data-simulated-lightgrey)

A collection of code adapted from real **public health and epidemiology** research projects: peer-reviewed articles, my master's thesis, statistical consulting, and disease surveillance projects.

> **Data privacy:** all original data have been replaced with **simulated data** (`set.seed`), and the names of objects, variables, places, and institutions have been changed. Each script runs on its own and gives reproducible results.

---

## Projects

| # | Project | Main techniques |
|---|---------|-----------------|
| 01 | [Descriptive analysis and Table 1]() | `gtsummary`, Fisher/Wilcoxon tests, export to Word |
| 02 | [Contingency tables and risk measures](02_contingency_tables_association_risk.R) | Chi-square, Fisher's exact, McNemar, Risk Ratio, Odds Ratio |
| 03 | [Group comparisons](03_group_comparison_ttest_wilcoxon_anova_tukey.R) | Shapiro-Wilk, Levene, t-test, Mann-Whitney + effect size, ANOVA + Tukey, Kruskal-Wallis + Dunn |
| 04 | [Correlation analysis](04_correlation_analysis.R) | Pearson, Spearman, `corrplot`, `GGally` |
| 05 | [Linear regression, Box-Cox and stepwise](05_linear_regression_boxcox_stepwise.R) | `lm`, VIF, Box-Cox transformation, stepwise selection, train/test validation (RMSE, R², MAE) |
| 06 | [Logistic regression and ROC curve](06_logistic_regression_roc_firth.R) | Crude/adjusted OR, AUC, Youden index, Firth regression, `brglm2` |
| 07 | [Classification / machine learning](07_classification_machine_learning.R) | Logit, multinomial, neural network (`nnet`), random forest, confusion matrix |
| 08 | [Count models](08_count_models_poisson_negbin_zip.R) | Poisson with offset (IRR), negative binomial, ZIP/ZINB, Vuong test, `DHARMa` |
| 09 | [Survival analysis](09_survival_analysis_kaplan_meier_cox.R) | Kaplan-Meier, log-rank, Cox regression, Schoenfeld residuals |
| 10 | [Propensity score matching](10_propensity_score_matching_matchit.R) | `MatchIt`, balance diagnostics (SMD), ATT with robust SE |
| 11 | [Probabilistic record linkage and deduplication](11_probabilistic_record_linkage_deduplication.R) | Name standardization, blocking, Soundex, Jaro-Winkler, EpiLink weights |
| 12 | [Epidemic control chart (endemic channel)](12_surveillance_epidemic_control_chart.R) | Epidemiological weeks, mean ± 2 SD, epidemic alerts |
| 13 | [Mortality indicators and YLL](13_mortality_indicators_incidence_yll.R) | Rates with 95% CI, LOESS trends, years of life lost |
| 14 | [ETL for large health databases](14_etl_datasus_duckdb_parquet.R) | DATASUS FTP, `.dbc` files, `DuckDB` (SQL), `Parquet`/`arrow` |

---

## Tech stack

| Area | Packages |
|------|----------|
| Data wrangling | `tidyverse`, `lubridate`, `stringr` |
| Tables | `gtsummary`, `flextable`, `gt` |
| Inference | `rstatix`, `coin`, `car`, `epitools` |
| Modeling | `MASS`, `pscl`, `logistf`, `brglm2`, `nnet`, `survival`, `MatchIt` |
| Machine learning | `caret`, `randomForest`, `NeuralNetTools`, `pROC` |
| Diagnostics | `DHARMa`, `nortest`, `lmtest`, `sandwich`, `AER` |
| Visualization | `ggplot2`, `ggpubr`, `survminer`, `corrplot`, `GGally` |
| Big data | `duckdb`, `DBI`, `arrow`, `read.dbc` |
| Record linkage | `RecordLinkage` |

---

## How to run

```r
packages <- c("tidyverse", "gtsummary", "broom.helpers", "flextable", "gt",
              "epitools", "rstatix", "coin", "ggpubr", "car", "corrplot",
              "GGally", "nortest", "fastDummies", "caret", "pROC", "logistf",
              "brglm2", "nnet", "NeuralNetTools", "randomForest", "MASS",
              "pscl", "lmtest", "AER", "DHARMa", "survival", "survminer",
              "broom", "gridExtra", "MatchIt", "sandwich", "RecordLinkage",
              "xts", "arrow", "duckdb", "DBI", "RCurl", "curl", "read.dbc")

install.packages(setdiff(packages, rownames(installed.packages())))

source("09_survival_analysis_kaplan_meier_cox.R")  # or any other script
```

All scripts were tested with R 4.6.1.
