# 08-logistic.R
# Unpenalised multivariable logistic regression on the 27 parameters (the
# reference learner of 02-model.R), fitted on each of the m completed data sets
# of the whole cohort and pooled by Rubin's rules; crude odds ratios likewise.
# Run from the repository root:  Rscript R/08-logistic.R
# Input : results/final-model.rds (completed data sets from 02b-final-model.R)
# Output: results/tables/table4-logistic.csv

source("R/utils.R")

imps <- readRDS("results/final-model.rds")$data

# continuous predictors on a clinically readable scale
scale_by <- c(age = 5, bmi = 5, booking_ga_weeks = 4, sbp = 10, dbp = 10,
              platelets_k = 50)
unit_label <- c(age = "per 5 years", bmi = "per 5 kg/m2", booking_ga_weeks = "per 4 weeks",
                sbp = "per 10 mmHg", dbp = "per 10 mmHg", platelets_k = "per 50 x 10^3/mm3",
                hb_gdl = "per 1 g/dL", gravidity = "per pregnancy", parity = "per birth")
rescale <- function(z) {
  for (v in names(scale_by)) z[[v]] <- z[[v]] / scale_by[[v]]
  z$y <- as.integer(z$y == "yes")
  z
}
imps <- lapply(imps, rescale)

# Rubin's rules on the log-odds scale
pool <- function(fits) {
  b <- sapply(fits, coef)
  v <- sapply(fits, function(m) diag(vcov(m)))
  m <- ncol(b)
  qbar <- rowMeans(b)
  ubar <- rowMeans(v)
  bvar <- apply(b, 1, var)
  tvar <- ubar + (1 + 1 / m) * bvar
  r <- (1 + 1 / m) * bvar / ubar
  df <- (m - 1) * (1 + 1 / r)^2
  q <- qt(0.975, df)
  data.frame(term = names(qbar), or = exp(qbar), lo = exp(qbar - q * sqrt(tvar)),
             hi = exp(qbar + q * sqrt(tvar)), p = 2 * pt(-abs(qbar / sqrt(tvar)), df))
}

adj <- pool(lapply(imps, function(z) glm(model_formula(predictors), binomial, z)))
crude <- do.call(rbind, lapply(predictors, function(v)
  pool(lapply(imps, function(z) glm(model_formula(v), binomial, z)))[-1, ]))

fmt <- function(x) sprintf("%.2f (%.2f to %.2f)", x$or, x$lo, x$hi)
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
terms <- adj$term[-1]
var_of <- ifelse(terms %in% predictors, terms, "education")
level <- ifelse(var_of == "education", sub("^education", "", terms), "")

d1 <- readRDS("results/final-model.rds")$data[[1]]
n_with <- vapply(var_of, function(v)
  if (v %in% binary_vars) sum(d1[[v]] == 1) else NA_real_, numeric(1))
ev_with <- vapply(var_of, function(v)
  if (v %in% binary_vars) sum(d1[[v]] == 1 & d1$y == "yes") else NA_real_, numeric(1))
# complete separation (every woman with the characteristic had the outcome):
# the odds ratio is not estimable
separated <- !is.na(n_with) & ev_with == n_with

out <- data.frame(
  predictor = ifelse(level == "", labels[var_of],
                     paste0(labels["education"], ": ", level, " (vs none)")),
  unit = ifelse(var_of %in% names(unit_label), unit_label[var_of], ""),
  women_with = n_with,
  events_with = ev_with,
  crude_or = ifelse(separated, "not estimable", fmt(crude[match(terms, crude$term), ])),
  adjusted_or = ifelse(separated, "not estimable", fmt(adj[-1, ])),
  adjusted_p = ifelse(separated, "", fmt_p(adj$p[-1]))
)
# intercept on the log-odds scale (reference: all binary predictors 0, education
# none, continuous predictors 0 on the scales above)
out <- rbind(data.frame(predictor = "Intercept (log odds)", unit = "", women_with = NA,
                        events_with = NA, crude_or = "",
                        adjusted_or = sprintf("%.2f (%.2f to %.2f)", log(adj$or[1]),
                                              log(adj$lo[1]), log(adj$hi[1])),
                        adjusted_p = ""), out)
write.csv(out, "results/tables/table4-logistic.csv", row.names = FALSE)
message("done: ", nrow(out), " terms, m = ", length(imps))
