# 02-model.R
# Development and internal validation of the complication risk model (protocol section 7).
# Run from the repository root:  OMP_NUM_THREADS=1 Rscript R/02-model.R
# Set N_REPEATS=1 N_SENS_REPEATS=1 N_IMP=2 for a quick smoke run.
# Input : data/processed/analysis.csv (from 01-clean.R)
# Output: results/cv-predictions.rds      repeated 5-fold CV, out-of-fold risks
#         results/cv-specs.rds            tuned hyperparameters per fold
#         results/sens-predictions.rds    sensitivity analyses S1-S3, S5-S7 (S4: 07-smote.R)
# Selection and final model: 02b-final-model.R; interpretation: 03-interpret.R

source("R/utils.R")
suppressPackageStartupMessages(library(parallel))

n_repeats      <- as.integer(Sys.getenv("N_REPEATS", 20))
n_sens_repeats <- as.integer(Sys.getenv("N_SENS_REPEATS", 5))
n_imp          <- as.integer(Sys.getenv("N_IMP", 5))
n_cores        <- min(detectCores() - 2, 14)
seed0          <- 20261007
dir.create("results", showWarnings = FALSE)

d <- load_analysis()
stopifnot(sum(d$y == "yes") == sum(d$y_primary == 1))

# ---- repeated cross-validation ----------------------------------------------
# one job per outer fold: impute once (m completed sets), then fit every learner

run_cv <- function(data, vars, repeats, v = 5, learners = names(learner_set)) {
  jobs <- expand.grid(fold = seq_len(v), rep = seq_len(repeats))
  folds <- lapply(seq_len(repeats), function(r) stratified_folds(data$y, v, seed0 + r))
  res <- mclapply(seq_len(nrow(jobs)), function(j) {
    r <- jobs$rep[j]; k <- jobs$fold[j]
    test_i <- which(folds[[r]] == k)
    seed <- seed0 + 1000 * r + k
    imp <- impute_fold(data[-test_i, ], data[test_i, ], vars, n_imp, seed)
    out <- lapply(learners, function(l) {
      t0 <- Sys.time()
      fits <- tryCatch(fit_mi(imp$train, learner_set[[l]], vars, seed),
                       error = function(e) {
                         message(l, " rep ", r, " fold ", k, ": ", conditionMessage(e))
                         NULL })
      p <- if (is.null(fits)) NA_real_ else predict_mi(fits, imp$test)
      list(pred = data.frame(rep = r, fold = k, learner = l, row = test_i,
                             y = data$y[test_i], p = p),
           spec = data.frame(rep = r, fold = k, learner = l,
                             spec = if (is.null(fits)) NA_character_ else
                               paste(deparse(fits[[1]]$spec[names(learner_set[[l]]$grid)]),
                                     collapse = ""),
                             secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    })
    list(pred = do.call(rbind, lapply(out, `[[`, "pred")),
         spec = do.call(rbind, lapply(out, `[[`, "spec")))
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) stop("failed folds: ", paste(which(bad), collapse = ", "))
  list(pred = do.call(rbind, lapply(res, `[[`, "pred")),
       spec = do.call(rbind, lapply(res, `[[`, "spec")))
}

message("repeated CV: ", n_repeats, " x 5 folds x ", length(learner_set),
        " learners, m = ", n_imp, " imputations")
cv <- run_cv(d, predictors, n_repeats)
saveRDS(cv$pred, "results/cv-predictions.rds")
saveRDS(cv$spec, "results/cv-specs.rds")
print(aggregate(secs ~ learner, cv$spec, sum))

# ---- sensitivity analyses (protocol 7.6) -------------------------------------

d_s2 <- d
d_s2$hb_gdl[!is.na(d_s2$hb_gdl_month) & d_s2$hb_gdl_month > 4] <- NA
d_s2$platelets_k[!is.na(d_s2$platelets_k_month) & d_s2$platelets_k_month > 4] <- NA

# S6: narrower outcome, gestational diabetes or hypertensive disorder
d_s6 <- d
y6 <- with(d, (c_gdm %in% 1) | (c_gest_htn %in% 1) | (c_preeclampsia %in% 1) |
              (c_eclampsia %in% 1))
d_s6$y <- factor(ifelse(y6, "yes", "no"), levels = c("no", "yes"))

sens <- list(
  S1_glycaemia   = list(data = d, vars = c(predictors, "glycaemia_gl")),
  S2_early_labs  = list(data = d_s2, vars = predictors),
  S3_booked_le20 = list(data = d[!is.na(d$booking_ga_weeks) & d$booking_ga_weeks <= 20, ],
                        vars = predictors),
  S5_process     = list(data = d, vars = c(predictors, process_vars)),
  S6_gdm_hdp     = list(data = d_s6, vars = predictors),
  S7_no_rural    = list(data = d, vars = setdiff(predictors, "rural"))
)
sens_pred <- do.call(rbind, lapply(names(sens), function(s) {
  message("sensitivity ", s, ": n = ", nrow(sens[[s]]$data),
          ", events = ", sum(sens[[s]]$data$y == "yes"))
  out <- run_cv(sens[[s]]$data, sens[[s]]$vars, n_sens_repeats)$pred
  out$analysis <- s
  out
}))
saveRDS(sens_pred, "results/sens-predictions.rds")

message("done; final model: R/02b-final-model.R")
