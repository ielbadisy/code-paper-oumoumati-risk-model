# 05-audit-2025.R
# Audit of the 2025 analysis (protocol section 9), for the supplement.
# Run from the repository root:  OMP_NUM_THREADS=1 Rscript R/05-audit-2025.R
# Input : data/processed/analysis.csv
#         data/audit-2025/maternalrisk_data_rabat.csv  (the 2025 analysis file)
# Output: results/tables/audit-composite.csv     old vs corrected composite
#         results/tables/audit-leakage.csv       single-variable AUC, old vs new outcome
#         results/tables/audit-label-baseline.csv
#         results/tables/audit-upsampling.csv    effect of up-sampling on predicted risk

source("R/utils.R")

d <- load_analysis()
seed0 <- 20261007
fmt_auc <- function(y01, x) {
  ok <- !is.na(x)
  auc(factor(ifelse(y01[ok] == 1, "yes", "no"), levels = c("no", "yes")), x[ok])
}

# ---- 1. composite outcome: old vs corrected ----------------------------------

comp <- c(c_anaemia = "Anaemia", c_scarred_uterus = "Scarred uterus",
          c_gdm = "Gestational diabetes", c_gest_htn = "Gestational hypertension",
          c_preeclampsia = "Pre-eclampsia", c_eclampsia = "Eclampsia",
          c_haemorrhage = "Haemorrhage", c_infection = "Infection",
          c_preterm = "Preterm birth", c_cardiopathy = "Cardiopathy")
audit_comp <- data.frame(
  component = unname(comp),
  n = sapply(names(comp), function(v) sum(d[[v]] == 1, na.rm = TRUE)),
  in_2025_composite = "yes",
  in_corrected_composite = ifelse(names(comp) %in% c("c_anaemia", "c_scarred_uterus",
                                                     "c_cardiopathy"), "no", "yes"),
  only_event_in_2025 = sapply(names(comp), function(v)
    sum(d[[v]] == 1 & d$y_primary == 0, na.rm = TRUE))
)
audit_comp <- rbind(audit_comp, data.frame(
  component = c("Any (2025 composite)", "Any (corrected composite)"),
  n = c(sum(d$y_old == 1), sum(d$y_primary == 1)),
  in_2025_composite = "", in_corrected_composite = "", only_event_in_2025 = NA))
audit_comp$percent <- round(100 * audit_comp$n / nrow(d), 1)
write.csv(audit_comp, "results/tables/audit-composite.csv", row.names = FALSE)

# ---- 2. leakage: predictors that define part of the old outcome --------------

leak <- data.frame(
  predictor = c("Haemoglobin (lower = higher risk)", "Gestational diabetes flag",
                "Prior cesarean"),
  auc_2025_outcome = c(fmt_auc(d$y_old, -d$hb_gdl), fmt_auc(d$y_old, d$c_gdm),
                       fmt_auc(d$y_old, d$prior_caesarean)),
  auc_corrected_outcome = c(fmt_auc(d$y_primary, -d$hb_gdl), NA,
                            fmt_auc(d$y_primary, d$prior_caesarean))
)
leak[-1] <- lapply(leak[-1], round, 3)
write.csv(leak, "results/tables/audit-leakage.csv", row.names = FALSE)

# ---- 3. label coding: baseline of the 2025 explainer -------------------------

old <- read.csv("data/audit-2025/maternalrisk_data_rabat.csv")
old <- old[!is.na(old$outcome_complication), ]
base <- data.frame(
  quantity = c("women in 2025 file", "outcome = 1 (complication)",
               "outcome = 0 (no complication)", "proportion outcome = 1",
               "proportion outcome = 0", "baseline reported in 2025 explainer"),
  value = c(nrow(old), sum(old$outcome_complication == 1),
            sum(old$outcome_complication == 0),
            round(mean(old$outcome_complication == 1), 3),
            round(mean(old$outcome_complication == 0), 3), 0.544)
)
write.csv(base, "results/tables/audit-label-baseline.csv", row.names = FALSE)

# ---- 4. up-sampling: predicted risk vs observed proportion -------------------

# This audit was run with the first version of the analysis (2026-10-07): the
# 22 prespecified predictors and single kNN imputation (missknn). Both are
# defined here so that the audit reproduces the reported values.
audit_binary <- c("rural", "consanguinity", "wanted_pregnancy", "prior_abortion",
                  "prior_adverse_outcome", "prior_caesarean", "family_diabetes",
                  "family_hypertension", "personal_medical_history",
                  "rhesus_negative")
predictors <- c("age", "gravidity", "parity", "bmi", "booking_ga_weeks",
                "sbp", "dbp", "hb_gdl", "platelets_k", audit_binary, "education")
round_binaries <- function(x, vars) {
  for (v in intersect(vars, names(x))) x[[v]] <- as.numeric(x[[v]] >= 0.5)
  x
}
# impute the training set, then the test set with the imputed training rows as
# donors; the outcome never enters the imputation
impute_split <- function(train, test, vars, seed) {
  knn <- function(x) as.data.frame(missknn::missknn(x, k = 5, seed = seed,
                                                     parallel_cores = 1))
  tr <- knn(train[vars])
  stacked <- rbind(tr, test[vars])
  te <- knn(stacked)[(nrow(tr) + 1):nrow(stacked), , drop = FALSE]
  list(train = round_binaries(tr, audit_binary), test = round_binaries(te, audit_binary))
}

# Logistic regression on the corrected outcome, 5-fold CV, in-fold imputation;
# with and without random up-sampling of events in the training folds.
folds <- stratified_folds(d$y, 5, seed0)
f <- as.formula(paste("y ~", paste(predictors, collapse = " + ")))
res <- do.call(rbind, lapply(1:5, function(k) {
  te <- which(folds == k)
  imp <- impute_split(d[-te, ], d[te, ], predictors, seed0 + k)
  tr <- cbind(imp$train, y = d$y[-te])
  set.seed(seed0 + k)
  ev <- which(tr$y == "yes")
  up <- rbind(tr, tr[sample(ev, sum(tr$y == "no") - length(ev), replace = TRUE), ])
  data.frame(row = te, y = d$y[te],
             p_plain = prob_yes(fit(f, tr, model = "glm"), imp$test),
             p_upsampled = prob_yes(fit(f, up, model = "glm"), imp$test))
}))
ups <- rbind(cbind(training = "no up-sampling", performance(res$y, res$p_plain),
                   mean_predicted = mean(res$p_plain)),
             cbind(training = "random up-sampling", performance(res$y, res$p_upsampled),
                   mean_predicted = mean(res$p_upsampled)))
ups$observed <- mean(res$y == "yes")
ups[-1] <- lapply(ups[-1], round, 3)
write.csv(ups, "results/tables/audit-upsampling.csv", row.names = FALSE)

print(audit_comp, row.names = FALSE); print(leak, row.names = FALSE)
print(base, row.names = FALSE); print(ups, row.names = FALSE)
