# 03-interpret.R
# Interpretation of the final model (protocol section 7.5): funcml interpret()
# methods (partial dependence and ALE dropped 2026-10-09). estimate() is in 03b-estimate.R.
# Run from the repository root:  OMP_NUM_THREADS=1 Rscript R/03-interpret.R
# Lower N_CORES if memory is short.
# Input : results/final-model.rds (from 02b-final-model.R)
# Output: results/interpret.rds   global and local explanations

source("R/utils.R")
suppressPackageStartupMessages(library(parallel))

seed0   <- 20261007
n_cores <- as.integer(Sys.getenv("N_CORES", min(detectCores() - 2, 14)))  # each worker holds a BART copy (~1.5 GB)
fm_path <- "results/final-model.rds"
fm <- readRDS(fm_path)

# explanations use the fit on the first completed dataset (of m imputations)
fit1 <- fm$fits[[1]]
d_full <- fm$data[[1]]
p_full <- prob_yes(fit1, d_full)
stopifnot(cor(p_full, as.integer(d_full$y == "yes")) > 0)  # sign check

f <- model_formula(predictors)
cont <- c("age", "bmi", "sbp", "dbp", "hb_gdl", "platelets_k", "booking_ga_weeks", "parity")
ip <- function(method, ...) interpret(fit1, d_full, f, method = method, class_level = "yes",
                                      seed = seed0, ...)

# BART fits are not safe to reuse across forks: each child reads its own copy
in_child <- function(jobs, fun) {
  res <- mclapply(jobs, function(j) {
    fm_child <- readRDS(fm_path)
    fun(j, fm_child$fits[[1]], fm_child$data[[1]])
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) stop(res[[which(bad)[1]]])
  res
}

# SHAP for many women, computed in chunks of rows and reassembled
shap_rows_of <- function(rows) {
  chunks <- split(rows, cut(seq_along(rows), min(n_cores, length(rows)), labels = FALSE))
  parts <- in_child(chunks, function(r, fit, z)
    interpret(fit, z, f, method = "shap", class_level = "yes", newdata = z[r, ],
              nsamples = 100, seed = seed0))
  out <- parts[[1]]
  offset <- 0L
  out$result <- do.call(rbind, lapply(parts, function(x) {
    x$result$observation <- x$result$observation + offset
    offset <<- offset + x$diagnostics$observations
    x$result
  }))
  out$diagnostics$observations <- length(rows)
  out$diagnostics$baseline <- out$diagnostics$prediction <- NA_real_
  out$newdata_ref <- d_full[rows, , drop = FALSE]
  out
}

# ---- global explanations -----------------------------------------------------

message("global explanations")
set.seed(seed0)
shap_rows <- c(which(d_full$y == "yes"), sample(which(d_full$y == "no"), 200))
vip <- ip("vip", metric = "auc", nsim = 20)
top6 <- head(vip$result$scores$feature[order(-vip$result$scores$importance)], 6)

interp <- list(
  learner = fm$learner,
  # agreement between the explained fit and the m-imputation average
  cor_mi = cor(p_full, predict_mi(fm$fits, fm$data)),
  vip = vip,
  permute = ip("permute", metric = "brier", compare = "ratio", nsim = 20),
  shap = shap_rows_of(shap_rows),
  ice = ip("ice", features = cont),
  profile = ip("profile", features = cont),
  interaction = ip("interaction", features = top6, nsamples = 100, grid_size = 10),
  surrogate = ip("surrogate"),
  # apparent (in-sample); the reported curves are out of fold (04-results.R)
  calibration = ip("calibration"),
  dca = ip("dca")
)

# ---- local explanations: women at the 10th, 50th and 90th risk percentiles ----

message("local explanations")
examples <- vapply(c(p10 = 0.10, p50 = 0.50, p90 = 0.90),
                   function(q) which.min(abs(p_full - quantile(p_full, q))), integer(1))
interp$examples <- data.frame(example = names(examples), row = unname(examples),
                              risk = p_full[examples], y = d_full$y[examples])
# funcml 0.9.5: "local", "lime" and "local_model" call the same weighted local
# linear surrogate, so it is computed once, keeping every coefficient (the
# default keeps 3 and drops zero contributions, often leaving none);
# "ceteris_paribus" ignores newdata, so the woman's row is passed as the data
# with the grid taken from the whole cohort
cp_grid <- lapply(setNames(cont, cont), function(v)
  seq(min(d_full[[v]]), max(d_full[[v]]), length.out = 51))
interp$local <- lapply(examples, function(i) {
  w <- d_full[i, , drop = FALSE]
  list(lime = ip("lime", newdata = w, k = length(predictors)),
       shap = ip("shap", newdata = w, nsamples = 100),
       ceteris_paribus = interpret(fit1, w, f, method = "ceteris_paribus", features = cont,
                                   grid = cp_grid, class_level = "yes", seed = seed0))
})
saveRDS(interp, "results/interpret.rds")
message("done: ", fm$learner)
