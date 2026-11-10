# 10-centre-only.R
# Sensitivity analysis S8: the centre as the only predictor.
# Each woman's predicted risk is the outcome proportion of her centre,
# estimated in the training folds (20 repeats of 5-fold cross-validation).
# Run from the repository root:  Rscript R/10-centre-only.R
# Input : data/processed/analysis.csv (from 01-clean.R)
# Output: results/tables/tableS-centre-only.csv

source("R/utils.R")

d <- load_analysis()
y <- d$y_primary
centre <- d$centre

# AUC as the Mann-Whitney statistic (ties count one half)
auc <- function(y, p) {
  r <- rank(p)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

set.seed(2026)
aucs <- unlist(lapply(1:20, function(k) {
  fold <- sample(rep(1:5, length.out = length(y)))
  sapply(1:5, function(j) {
    train <- fold != j
    rate <- tapply(y[train], centre[train], mean)
    auc(y[!train], rate[centre[!train]])
  })
}))

res <- data.frame(
  analysis = "S8 centre only",
  auc_mean = mean(aucs),
  auc_p2.5 = unname(quantile(aucs, 0.025)),
  auc_p97.5 = unname(quantile(aucs, 0.975)),
  auc_apparent = auc(y, ave(y, centre))
)
print(round(res[-1], 3))
write.csv(res, "results/tables/tableS-centre-only.csv", row.names = FALSE)
