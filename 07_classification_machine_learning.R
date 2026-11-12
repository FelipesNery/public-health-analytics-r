# =============================================================================
# PROJECT 07 - Supervised Classification: Logit, Multinomial, Neural Network
#              and Random Forest
# =============================================================================
#
# DESCRIPTION
# Building and comparing classification models for diagnosis (binary and
# multiclass) from clinical/laboratory markers, with a train/test split,
# variable selection and evaluation through confusion matrices.
#
# WHAT WAS DONE
#  - Random train (70%) / test (30%) split
#  - Binary logistic regression + stepwise, with accuracy comparison
#  - Variable importance (caret::varImp) in a lollipop chart
#  - Neural network with one hidden layer (nnet) and architecture plot
#  - Random Forest with imputation (na.roughfix) and importance (Gini)
#  - Multinomial logistic regression (3-class diagnosis)
#  - Confusion matrices: accuracy, sensitivity, specificity, Kappa
#
# TECHNIQUES / PACKAGES: stats, nnet, randomForest, caret, NeuralNetTools
# NOTE: simulated data; variable and object names anonymized.
# =============================================================================

library(tidyverse)
library(nnet)

set.seed(9973)

# -----------------------------------------------------------------------------
# 1. Simulated data (binary diagnosis from 8 markers)
# -----------------------------------------------------------------------------
n_c <- 500
lab_ml <- as_tibble(matrix(rnorm(n_c * 8), ncol = 8,
                           dimnames = list(NULL, paste0("mk_", 1:8))))
lab_ml <- lab_ml |>
  mutate(
    cls_bin = rbinom(n_c, 1, plogis(-0.4 + 1.3 * mk_1 - 0.9 * mk_3 + 0.7 * mk_5)),
    cls_bin = factor(cls_bin, levels = c(0, 1))
  )
lab_ml$mk_7[sample(n_c, 25)] <- NA            # some missing values

# -----------------------------------------------------------------------------
# 2. Train / test
# -----------------------------------------------------------------------------
idx_ml   <- sample(seq_len(n_c), round(0.7 * n_c))
train_ml <- lab_ml[idx_ml, ]  |> drop_na()
test_ml  <- lab_ml[-idx_ml, ] |> drop_na()

conf_mtx <- function(pred, obs) {
  tb <- table(Predicted = pred, Observed = obs)
  tp <- tb["1", "1"]; tn <- tb["0", "0"]; fp <- tb["1", "0"]; fn <- tb["0", "1"]
  list(table = tb,
       accuracy = (tp + tn) / sum(tb),
       sensitivity = tp / (tp + fn),
       specificity = tn / (tn + fp))
}

# -----------------------------------------------------------------------------
# 3. Logistic regression + stepwise
# -----------------------------------------------------------------------------
fit_lgt <- glm(cls_bin ~ ., data = train_ml, family = binomial)
summary(fit_lgt)
car::vif(fit_lgt)

test_ml <- test_ml |>
  mutate(p_lgt   = predict(fit_lgt, newdata = test_ml, type = "response"),
         cls_lgt = factor(if_else(p_lgt > 0.5, 1, 0), levels = c(0, 1)))
res_lgt <- conf_mtx(test_ml$cls_lgt, test_ml$cls_bin)
res_lgt

fit_stp <- step(fit_lgt, k = qchisq(0.05, df = 1, lower.tail = FALSE), trace = 0)
summary(fit_stp)
cat("Variables removed:", length(coef(fit_lgt)) - length(coef(fit_stp)), "\n")

test_ml <- test_ml |>
  mutate(p_stp   = predict(fit_stp, newdata = test_ml, type = "response"),
         cls_stp = factor(if_else(p_stp > 0.5, 1, 0), levels = c(0, 1)))
res_stp <- conf_mtx(test_ml$cls_stp, test_ml$cls_bin)

tibble(Full = res_lgt$accuracy, Stepwise = res_stp$accuracy)

# Variable importance (|z statistic|)
if (requireNamespace("caret", quietly = TRUE)) {
  imp_v <- caret::varImp(fit_stp)
  caret::confusionMatrix(test_ml$cls_stp, test_ml$cls_bin, positive = "1")
} else {
  imp_v <- data.frame(Overall = abs(summary(fit_stp)$coefficients[-1, "z value"]))
}

ggplot(imp_v, aes(x = reorder(rownames(imp_v), Overall), y = Overall)) +
  geom_segment(aes(xend = rownames(imp_v), y = 0, yend = Overall), color = "skyblue") +
  geom_point(color = "blue", size = 4, alpha = 0.6) +
  coord_flip() +
  labs(x = "Variable", y = "Importance") +
  theme_light()

# -----------------------------------------------------------------------------
# 4. Neural network (1 hidden layer)
# -----------------------------------------------------------------------------
sel_vars <- names(coef(fit_stp))[-1]
frm_nn   <- as.formula(paste("cls_bin ~", paste(sel_vars, collapse = " + ")))

fit_nn <- nnet(frm_nn, data = train_ml, size = 6, decay = 0.01,
               maxit = 500, trace = FALSE)
if (requireNamespace("NeuralNetTools", quietly = TRUE)) NeuralNetTools::plotnet(fit_nn)

pred_nn <- factor(predict(fit_nn, newdata = test_ml, type = "class"), levels = c(0, 1))
res_nn  <- conf_mtx(pred_nn, test_ml$cls_bin)
res_nn

# -----------------------------------------------------------------------------
# 5. Random Forest (with NA imputation)
# -----------------------------------------------------------------------------
if (requireNamespace("randomForest", quietly = TRUE)) {
  fit_rf <- randomForest::randomForest(cls_bin ~ ., data = lab_ml[idx_ml, ],
                                       na.action = randomForest::na.roughfix,
                                       ntree = 500, importance = TRUE)
  print(fit_rf)                                # OOB error
  randomForest::importance(fit_rf)
  randomForest::varImpPlot(fit_rf, main = "Variable importance - RF")

  pred_rf <- predict(fit_rf, newdata = test_ml)
  res_rf  <- conf_mtx(pred_rf, test_ml$cls_bin)

  tibble(Logit = res_lgt$accuracy, Stepwise = res_stp$accuracy,
         NeuralNet = res_nn$accuracy, RandomForest = res_rf$accuracy)
}

# -----------------------------------------------------------------------------
# 6. Multinomial logistic regression (3 classes)
# -----------------------------------------------------------------------------
n_m <- 400
diag_mn <- tibble(
  age_m  = rnorm(n_m, 70, 8),
  sex_m  = factor(sample(c("F", "M"), n_m, TRUE)),
  cog_sc = rnorm(n_m, 25, 3),
  bio_a  = rnorm(n_m, 600, 150),
  bio_b  = rnorm(n_m, 80, 25)
) |>
  mutate(
    eta_b = -0.5 - 0.25 * (cog_sc - 25) + 0.004 * (bio_b - 80),
    eta_c = -1.5 - 0.55 * (cog_sc - 25) - 0.004 * (bio_a - 600),
    pr_a  = 1 / (1 + exp(eta_b) + exp(eta_c)),
    pr_b  = exp(eta_b) * pr_a,
    u_m   = runif(n_m),
    cls_mn = case_when(u_m < pr_a ~ "Control",
                       u_m < pr_a + pr_b ~ "Mild",
                       TRUE ~ "Severe"),
    cls_mn = relevel(factor(cls_mn), ref = "Control")
  ) |>
  select(cls_mn, age_m, sex_m, cog_sc, bio_a, bio_b)

fit_mn <- multinom(cls_mn ~ ., data = diag_mn, trace = FALSE)
summary(fit_mn)

# Relative risk ratios (RRR) and p-values (Wald test)
exp(coef(fit_mn))
z_mn <- summary(fit_mn)$coefficients / summary(fit_mn)$standard.errors
round(2 * (1 - pnorm(abs(z_mn))), 4)

diag_mn$pred_mn <- predict(fit_mn, newdata = diag_mn, type = "class")
table(Predicted = diag_mn$pred_mn, Observed = diag_mn$cls_mn)
mean(diag_mn$pred_mn == diag_mn$cls_mn)        # accuracy
