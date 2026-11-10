# 00-sample-size.R
# Minimum sample size for a binary-outcome prediction model (Riley et al. criteria).
# Run from the repository root:  Rscript R/00-sample-size.R
# Output: results/tables/sample-size.csv
# Inputs fixed in protocol.md, section 6.

suppressPackageStartupMessages(library(testflow))

n_params <- 27      # candidate predictor parameters (protocol section 4)
phi      <- 0.133   # outcome proportion (103 / 776)
c_stat   <- 0.70    # anticipated C statistic
n_avail  <- 776
events   <- 103

# Cox-Snell R2 implied by the C statistic and prevalence: simulate a
# binormal linear predictor with that AUC, fit a logistic model, take its LR.
set.seed(20261007)
n_sim <- 1e6
y  <- rbinom(n_sim, 1, phi)
x  <- rnorm(n_sim, mean = sqrt(2) * qnorm(c_stat) * y)
lr <- with(glm(y ~ x, family = binomial), null.deviance - deviance)
r2_cs <- 1 - exp(-lr / n_sim)

# criterion 1: expected uniform shrinkage of at least 0.9
s1 <- 0.9
n1 <- n_params / ((s1 - 1) * log(1 - r2_cs / s1))

# criterion 2: optimism in Nagelkerke R2 of at most 0.05
ll_null <- phi * log(phi) + (1 - phi) * log(1 - phi)
max_r2  <- 1 - exp(2 * ll_null)
s2 <- r2_cs / (r2_cs + 0.05 * max_r2)
n2 <- n_params / ((s2 - 1) * log(1 - r2_cs / s2))

# criterion 3: overall risk estimated within +/- 0.05 (95% Wilson interval)
n3 <- suppressWarnings(sample_size_precision(
  endpoint = "binary", design = "one_sample", width = 0.05, p = phi,
  method = "wilson"))$n

# expected shrinkage achieved with the available sample (inverting criterion 1)
f <- function(s) n_params / ((s - 1) * log(1 - r2_cs / s)) - n_avail
s_avail <- uniroot(f, c(r2_cs + 1e-6, 0.999))$root

out <- data.frame(
  quantity = c("anticipated C statistic", "anticipated Cox-Snell R2",
               "max Cox-Snell R2", "criterion 1 n (shrinkage 0.9)",
               "criterion 2 n (optimism 0.05)", "criterion 3 n (risk +/- 0.05)",
               "required n (max)", "required events", "available n",
               "available events", "events per parameter",
               "expected shrinkage with available n"),
  value = c(c_stat, r2_cs, max_r2, ceiling(n1), ceiling(n2), n3,
            max(ceiling(c(n1, n2, n3))), ceiling(max(c(n1, n2, n3)) * phi),
            n_avail, events, events / n_params, s_avail)
)
out$value <- round(out$value, 3)

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "results/tables/sample-size.csv", row.names = FALSE)
print(out, row.names = FALSE)
