# 03b-estimate.R
# Standardised risk differences for each binary predictor (protocol section 7.5):
# plug-in g-computation on the final model with funcml estimate(); the
# bootstrap refits BART n_boot times. Model-based, not causal effects.
# Run from the repository root:  OMP_NUM_THREADS=1 Rscript R/03b-estimate.R
# Set N_BOOT=20 for a quick run; lower N_CORES if memory is short.
# Input : results/final-model.rds (from 02b-final-model.R)
# Output: results/estimate.rds; one file per job in results/estimate-jobs/ (resumable)

source("R/utils.R")
suppressPackageStartupMessages(library(parallel))

seed0   <- 20261007
n_boot  <- as.integer(Sys.getenv("N_BOOT", 200))
n_cores <- as.integer(Sys.getenv("N_CORES", 10))  # each worker holds a BART copy (~1.5 GB)
min_n   <- 10  # minimum women with and without the predictor in the target population
fm_path <- "results/final-model.rds"
job_dir <- "results/estimate-jobs"
dir.create(job_dir, showWarnings = FALSE)

d_full <- readRDS(fm_path)$data[[1]]
f <- model_formula(predictors)
subgroups <- list(all = rep(TRUE, nrow(d_full)), nulliparous = d_full$parity == 0,
                  age35 = d_full$age >= 35)

jobs <- expand.grid(var = binary_vars, estimand = c("ATE", "ATT", "CATE_nulliparous",
                                                    "CATE_age35", "IATE"),
                    stringsAsFactors = FALSE)
jobs$population <- ifelse(grepl("^CATE_", jobs$estimand), sub("CATE_", "", jobs$estimand), "all")
jobs$n_with <- mapply(function(v, g) sum(d_full[[v]][subgroups[[g]]] == 1), jobs$var, jobs$population)
jobs$n_without <- mapply(function(v, g) sum(d_full[[v]][subgroups[[g]]] == 0), jobs$var, jobs$population)
jobs$run <- jobs$n_with >= min_n & jobs$n_without >= min_n
jobs$file <- file.path(job_dir, paste0(jobs$var, "-", jobs$estimand, ".rds"))
message("estimate(): ", sum(jobs$run), " of ", nrow(jobs), " jobs eligible, n_boot = ", n_boot)

todo <- which(jobs$run & !file.exists(jobs$file))
res <- mclapply(todo, function(i) {
  j <- jobs[i, ]
  fit <- readRDS(fm_path)$fits[[1]]  # BART fits are not safe to reuse across forks
  estimand <- sub("_.*", "", j$estimand)
  nd <- if (estimand == "CATE") d_full[subgroups[[j$population]], ] else NULL
  e <- try(estimate(d_full, f, fit = fit, treatment = j$var, estimand = estimand, newdata = nd,
                    interval = if (estimand == "IATE") "normal" else "bootstrap",
                    n_boot = n_boot, seed = seed0), silent = TRUE)
  if (!inherits(e, "try-error")) saveRDS(e, j$file)
  if (inherits(e, "try-error")) paste(j$var, j$estimand, ":", conditionMessage(attr(e, "condition"))) else NA
}, mc.cores = n_cores, mc.preschedule = FALSE)
err <- na.omit(unlist(res))
if (length(err)) message("failed jobs:\n", paste(err, collapse = "\n"))

jobs$done <- file.exists(jobs$file)
est <- lapply(jobs$file[jobs$done], readRDS)
names(est) <- paste(jobs$var, jobs$estimand, sep = ":")[jobs$done]
saveRDS(list(jobs = jobs, estimates = est, subgroups = subgroups, min_n = min_n),
        "results/estimate.rds")
message("done: ", sum(jobs$done), " of ", sum(jobs$run), " eligible jobs")
