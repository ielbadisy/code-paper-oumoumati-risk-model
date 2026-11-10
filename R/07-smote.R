# 07-smote.R
# Sensitivity analysis S4 (protocol 7.6): SMOTE applied to each completed
# training set inside the outer folds, compared with the same learners and
# folds without rebalancing (repeats 1-5 of 02-model.R).
# Run from the repository root:  OMP_NUM_THREADS=1 Rscript R/07-smote.R
# Input : data/processed/analysis.csv, results/cv-predictions.rds
# Output: results/smote-predictions.rds, results/tables/tableS-smote.csv

source("R/utils.R")
suppressPackageStartupMessages(library(parallel))

n_repeats <- as.integer(Sys.getenv("N_SENS_REPEATS", 5))
n_imp     <- as.integer(Sys.getenv("N_IMP", 5))
n_cores   <- min(detectCores() - 2, 14)
seed0     <- 20261007   # same folds and imputation seeds as 02-model.R
learners  <- c("glmnet", "ranger", "bart")

d <- load_analysis()

# SMOTE to a 1:1 class ratio (k = 5 neighbours) on one completed training set;
# synthetic binary, count and education values are rounded back to observed codes
smote_train <- function(z, seed) {
  set.seed(seed)
  lev <- levels(z$education)
  z$education <- as.integer(z$education)
  s <- as.data.frame(themis::smote(z, "y", k = 5, over_ratio = 1))
  for (v in c(binary_vars, "gravidity", "parity", "education")) s[[v]] <- round(s[[v]])
  s$education <- factor(lev[s$education], levels = lev)
  s[c(predictors, "y")]
}

jobs <- expand.grid(fold = 1:5, rep = seq_len(n_repeats))
folds <- lapply(seq_len(n_repeats), function(r) stratified_folds(d$y, 5, seed0 + r))
res <- mclapply(seq_len(nrow(jobs)), function(j) {
  r <- jobs$rep[j]; k <- jobs$fold[j]
  test_i <- which(folds[[r]] == k)
  seed <- seed0 + 1000 * r + k
  imp <- impute_fold(d[-test_i, ], d[test_i, ], predictors, n_imp, seed)
  trains <- lapply(seq_along(imp$train), function(i) smote_train(imp$train[[i]], seed + i))
  do.call(rbind, lapply(learners, function(l) {
    fits <- fit_mi(trains, learner_set[[l]], predictors, seed)
    data.frame(rep = r, fold = k, learner = l, row = test_i, y = d$y[test_i],
               p = predict_mi(fits, imp$test))
  }))
}, mc.cores = n_cores, mc.preschedule = FALSE)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop("failed folds: ", paste(which(bad), collapse = ", "))
smote <- do.call(rbind, res)
smote$analysis <- "S4_smote"
saveRDS(smote, "results/smote-predictions.rds")

ref <- readRDS("results/cv-predictions.rds")
ref <- ref[ref$rep <= n_repeats & ref$learner %in% learners, ]
ref$analysis <- "no_rebalancing"
both <- rbind(ref[names(smote)], smote)

per_rep <- do.call(rbind, lapply(split(both, list(both$analysis, both$learner, both$rep),
                                       drop = TRUE),
  function(s) cbind(analysis = s$analysis[1], learner = s$learner[1], rep = s$rep[1],
                    performance(s$y, s$p), mean_risk = mean(s$p))))
tab <- summarise_repeats(per_rep, by = c("analysis", "learner"))
fmt <- function(k) sprintf("%.3f (%.3f to %.3f)", tab[[k]], tab[[paste0(k, "_lo")]],
                           tab[[paste0(k, "_hi")]])
out <- data.frame(analysis = tab$analysis, learner = tab$learner)
for (k in c("auc", "citl", "slope", "oe", "scaled_brier", "ece", "mean_risk")) out[[k]] <- fmt(k)
out <- out[order(out$learner, out$analysis != "no_rebalancing"), ]
write.csv(out, "results/tables/tableS-smote.csv", row.names = FALSE)
print(out, row.names = FALSE)
