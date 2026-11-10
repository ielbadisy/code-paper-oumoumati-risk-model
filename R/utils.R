# utils.R
# Shared definitions for 02-model.R, 03-interpret.R and 04-results.R (protocol sections 4 and 7).

# one thread per process: parallelism is across folds (02-model.R)
Sys.setenv(OMP_NUM_THREADS = 1, R_RANGER_NUM_THREADS = 1)
options(ranger.num.threads = 1)

# mimar is not attached: its fit(), complete() and evaluate() would mask funcml
suppressPackageStartupMessages({
  library(funcml)
  library(ranger)   # predict() on fits read back from rds
  library(xgboost)
})

# ---- predictors (protocol section 4) ----------------------------------------

binary_vars <- c("rural", "consanguinity", "wanted_pregnancy", "prior_abortion",
                 "prior_adverse_outcome", "prior_caesarean", "family_diabetes",
                 "family_hypertension", "personal_medical_history",
                 "rhesus_negative", "surgical_history", "prior_htag",
                 "prior_haemorrhage", "pallor", "vaginal_discharge")
numeric_vars <- c("age", "gravidity", "parity", "bmi", "booking_ga_weeks",
                  "sbp", "dbp", "hb_gdl", "platelets_k")
factor_vars <- "education"
predictors <- c(numeric_vars, binary_vars, factor_vars)
# care-process markers, sensitivity analysis S5 only (protocol 7.6)
process_vars <- c("serology_done", "rubella_immune", "hb_missing", "platelets_missing")

labels <- c(
  age = "Age (years)", gravidity = "Gravidity", parity = "Parity",
  bmi = "BMI (kg/m2)", booking_ga_weeks = "Gestational age at booking (weeks)",
  sbp = "Systolic BP (mmHg)", dbp = "Diastolic BP (mmHg)",
  hb_gdl = "Hemoglobin (g/dL)", platelets_k = "Platelets (10^3/mm3)",
  glycaemia_gl = "Glycemia (g/L)",
  rural = "Rural residence", consanguinity = "Consanguinity",
  wanted_pregnancy = "Wanted pregnancy", prior_abortion = "Prior abortion",
  prior_adverse_outcome = "Prior adverse obstetric outcome",
  prior_caesarean = "Prior cesarean", family_diabetes = "Family history of diabetes",
  family_hypertension = "Family history of hypertension",
  personal_medical_history = "Personal medical history",
  rhesus_negative = "Rhesus negative", education = "Education",
  surgical_history = "Surgical history", prior_htag = "Prior hypertensive disorder",
  prior_haemorrhage = "Prior delivery hemorrhage", pallor = "Conjunctival pallor",
  vaginal_discharge = "Abnormal vaginal discharge", serology_done = "Serologies done",
  rubella_immune = "Rubella immune", hb_missing = "Hemoglobin not measured",
  platelets_missing = "Platelets not measured"
)

load_analysis <- function(path = "data/processed/analysis.csv") {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$education <- factor(d$education, levels = c("none", "primary", "middle", "secondary"))
  d$y <- factor(ifelse(d$y_primary == 1, "yes", "no"), levels = c("no", "yes"))
  d$hb_missing <- as.integer(is.na(d$hb_gdl))
  d$platelets_missing <- as.integer(is.na(d$platelets_k))
  d
}

# ---- learners (protocol 7.1) ------------------------------------------------
# grid: tuned by inner 5-fold CV on log loss; spec: fixed settings passed to fit()

learner_set <- list(
  glm      = list(model = "glm"),
  glmnet   = list(model = "glmnet",
                  grid = expand.grid(alpha = c(0, 0.5, 1),
                                     lambda = 10^seq(-4, -0.5, length.out = 8))),
  earth    = list(model = "earth",
                  grid = expand.grid(degree = 1:2, nprune = c(3, 5, 8, 12))),
  kknn     = list(model = "kknn", grid = expand.grid(k = c(25, 50, 100))),
  svm      = list(model = "e1071_svm",
                  grid = expand.grid(cost = c(0.25, 1, 4),
                                     gamma = c(0.005, 0.02, 0.05))),
  ranger   = list(model = "ranger",
                  grid = expand.grid(num.trees = 500, mtry = c(2, 4, 8),
                                     min.node.size = c(10, 25, 50))),
  xgboost  = list(model = "xgboost",
                  grid = expand.grid(max_depth = c(1, 2, 3), eta = c(0.03, 0.1),
                                     nrounds = c(100, 300))),
  lightgbm = list(model = "lightgbm",
                  grid = expand.grid(num_leaves = c(4, 8, 16),
                                     learning_rate = c(0.03, 0.1),
                                     nrounds = c(100, 300))),
  bart     = list(model = "bart"),
  densemlp = list(model = "densemlp",
                  grid = expand.grid(lr = c(1e-3, 3e-3, 1e-2), dropout = c(0.1, 0.3)),
                  spec = list(hidden_units = c(32L, 16L), batch_norm = FALSE,
                              ensemble = 5L, epochs = 200L, patience = 30L)),
  superlearner = list(model = "superlearner",
                      spec = list(learners = c("glmnet", "earth", "e1071_svm",
                                               "ranger", "bart"),
                                  meta_model = "glmnet"))
)
simplicity_order <- names(learner_set)

# ---- folds, imputation, fitting ---------------------------------------------

stratified_folds <- function(y, v, seed) {
  set.seed(seed)
  f <- integer(length(y))
  for (lev in unique(y)) {
    i <- which(y == lev)
    f[i] <- sample(rep_len(seq_len(v), length(i)))
  }
  f
}

# Multiple imputation inside a training set (protocol 7.2): chained equations
# with predictive mean matching (mimar), outcome included as an imputation
# predictor. Test rows are imputed by stacking them under the observed training
# rows with their outcome set to missing, so the test outcome is never used;
# only the test rows of that second imputation are kept.
# Education (ordered, 4 levels) is imputed as its integer code: pmm draws only
# observed codes, and mimar 1.4.0 polyreg fails on a factor with exactly one
# missing value.
impute_mi <- function(x, m, seed) {
  lev <- levels(x$education)
  x$education <- as.integer(x$education)
  imp <- mimar::impute(x, m = m, imputer = "pmm", maxit = 5, seed = seed,
                       progress = FALSE)
  lapply(seq_len(m), function(i) {
    z <- as.data.frame(mimar::complete(imp, i))
    z$education <- factor(lev[z$education], levels = lev)
    z
  })
}

impute_fold <- function(train, test, vars, m, seed) {
  tr <- train[c(vars, "y")]
  out <- list(train = impute_mi(tr, m, seed))
  if (is.null(test)) return(out)
  te <- test[c(vars, "y")]
  te$y <- factor(NA, levels = levels(train$y))
  stacked <- rbind(tr, te)
  keep <- (nrow(tr) + 1):nrow(stacked)
  out$test <- lapply(impute_mi(stacked, m, seed + 1), function(z) z[keep, vars, drop = FALSE])
  out
}

model_formula <- function(vars) as.formula(paste("y ~", paste(vars, collapse = " + ")))

# tuned (or fixed-spec) fit on one completed training set
fit_learner <- function(train, learner, vars, seed, inner_v = 5) {
  f <- model_formula(vars)
  if (is.null(learner$grid)) {
    return(fit(f, train, model = learner$model, spec = learner$spec, seed = seed))
  }
  args <- c(list(data = train, formula = f, model = learner$model, grid = learner$grid,
                 resampling = cv(inner_v, seed = seed), metric = "logloss", seed = seed),
            learner$spec)
  do.call(tune, args)$fit_best
}

# hyperparameters tuned on the first completed training set, then the learner
# is refitted with those settings on each of the m completed sets
fit_mi <- function(trains, learner, vars, seed) {
  first <- fit_learner(trains[[1]], learner, vars, seed)
  c(list(first), lapply(trains[-1], function(z)
    fit(model_formula(vars), z, model = learner$model, spec = first$spec, seed = seed)))
}

prob_yes <- function(fit, newdata) {
  p <- predict(fit, newdata, type = "prob")
  stopifnot(identical(colnames(p), c("no", "yes")))
  p[, "yes"]
}

# dbarts keeps the BART sampler state outside R; accessing it stores it in the
# fit so that a fit saved with saveRDS() still predicts after readRDS()
serialisable <- function(fit) {
  b <- fit$state$state
  if (inherits(b, "bart") && !is.null(b$fit)) invisible(b$fit$state)
  fit
}

# predicted risk averaged over the m imputations
predict_mi <- function(fits, tests) {
  rowMeans(sapply(seq_along(fits), function(i) prob_yes(fits[[i]], tests[[i]])))
}

# ---- performance (protocol 7.3) ---------------------------------------------

clip <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)

performance <- function(y, p) {
  stopifnot(is.factor(y), identical(levels(y), c("no", "yes")))
  yy <- as.integer(y == "yes")
  lp <- qlogis(clip(p))
  brier <- mean((yy - p)^2)
  prev <- mean(yy)
  data.frame(
    auc = auc(y, p),
    citl = unname(coef(glm(yy ~ 1 + offset(lp), family = binomial))[1]),
    slope = unname(coef(glm(yy ~ lp, family = binomial))[2]),
    oe = prev / mean(p),
    brier = brier,
    scaled_brier = 1 - brier / (prev * (1 - prev)),
    ece = ece(y, p, positive = "yes"),
    logloss = -mean(yy * log(clip(p)) + (1 - yy) * log(1 - clip(p)))
  )
}

# pooled-per-repeat metrics -> mean and 2.5-97.5 percentiles across repeats
summarise_repeats <- function(m, by = "learner") {
  metrics <- setdiff(names(m), c(by, "rep"))
  out <- lapply(split(m, m[by], drop = TRUE), function(s) {
    r <- data.frame(s[1, by, drop = FALSE], row.names = NULL)
    for (k in metrics) {
      r[[k]] <- mean(s[[k]])
      r[[paste0(k, "_lo")]] <- unname(quantile(s[[k]], 0.025))
      r[[paste0(k, "_hi")]] <- unname(quantile(s[[k]], 0.975))
    }
    r
  })
  do.call(rbind, out)
}

# selection rule (protocol 7.4)
select_learner <- function(summary) {
  ok <- summary[summary$slope >= 0.8 & summary$slope <= 1.2, ]
  if (nrow(ok) == 0) return("glmnet")
  best <- max(summary$auc)
  ok <- ok[ok$auc >= best - 0.02, ]
  if (nrow(ok) == 0) return("glmnet")
  simplicity_order[simplicity_order %in% ok$learner][1]
}
