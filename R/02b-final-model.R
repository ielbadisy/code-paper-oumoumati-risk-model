# 02b-final-model.R
# Learner selection and final model (protocol sections 7.4-7.5).
# Run from the repository root:  OMP_NUM_THREADS=1 Rscript R/02b-final-model.R
# Input : data/processed/analysis.csv, results/cv-predictions.rds (02-model.R)
# Output: results/final-model.rds   final learner refitted on all women

source("R/utils.R")

n_imp <- as.integer(Sys.getenv("N_IMP", 5))
seed0 <- 20261007

d <- load_analysis()
cv_pred <- readRDS("results/cv-predictions.rds")

per_rep <- do.call(rbind, lapply(split(cv_pred, list(cv_pred$learner, cv_pred$rep), drop = TRUE),
  function(s) cbind(learner = s$learner[1], rep = s$rep[1], performance(s$y, s$p))))
cv_summary <- summarise_repeats(per_rep)
rule_selected <- select_learner(cv_summary)

# v0.4 overrode the rule (random forest) with BART; from v0.5 the rule itself
# selects BART (random forest slope above 1.2), so no override is applied
final_learner <- rule_selected
message("rule-selected learner: ", rule_selected, "; final learner: ", final_learner)

imp_all <- impute_fold(d, NULL, predictors, n_imp, seed0)$train
final_fits <- lapply(fit_mi(imp_all, learner_set[[final_learner]], predictors, seed0),
                     serialisable)
p_full <- predict_mi(final_fits, imp_all)
stopifnot(cor(p_full, as.integer(d$y == "yes")) > 0)  # sign check (protocol 7.5)
saveRDS(list(learner = final_learner, rule_selected = rule_selected, fits = final_fits,
             data = imp_all, cv_summary = cv_summary), "results/final-model.rds")
fm <- readRDS("results/final-model.rds")
stopifnot(all.equal(predict_mi(fm$fits, fm$data), p_full))  # fits survive saveRDS
