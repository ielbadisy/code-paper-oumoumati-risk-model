# 04-results.R
# Tables and figures for the manuscript (protocol section 8). Figures are the
# funcml plot methods (compare, roc, calibration, dca, permute, shap)
# applied to out-of-fold predictions and to the interpret() objects.
# Run from the repository root:  Rscript R/04-results.R
# Input : data/processed/analysis.csv, results/*.rds (02-model.R, 03-interpret.R)
# Output: results/tables/table1-characteristics.csv
#         results/tables/table2-cv-performance.csv
#         results/tables/table3-classification.csv
#         results/tables/tableS-sensitivity.csv
#         results/tables/tableS-importance.csv
#         results/tables/tableS-tuned-specs.csv
#         results/figures/fig1-compare.pdf, fig2-roc.pdf, fig3-calibration.pdf,
#         fig4-dca.pdf, fig5-importance.pdf,
#         figS-roc-learners.pdf, figS-calibration-learners.pdf, figS-dca-learners.pdf,
#         figS-ice-profile.pdf, figS-interaction-surrogate.pdf, figS-local.pdf,
#         figS-apparent.pdf, fig6-estimate.pdf; tables/tableS-estimate.csv,
#         tables/tableS-dca-final.csv

source("R/utils.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

d <- load_analysis()
cv_pred <- readRDS("results/cv-predictions.rds")
cv_spec <- readRDS("results/cv-specs.rds")
sens_pred <- readRDS("results/sens-predictions.rds")
fm <- readRDS("results/final-model.rds")
interp <- readRDS("results/interpret.rds")
est <- readRDS("results/estimate.rds")
selected <- fm$learner
learner_names <- c(glm = "Logistic regression", glmnet = "Elastic net",
                   earth = "MARS", kknn = "k-nearest neighbors",
                   svm = "Support vector machine", ranger = "Random forest",
                   xgboost = "XGBoost", lightgbm = "LightGBM", bart = "BART",
                   densemlp = "Dense MLP", superlearner = "Super learner")
stopifnot(setequal(names(learner_names), names(learner_set)))

fmt <- function(x, k = 2) formatC(x, format = "f", digits = k)
write_table <- function(x, name) write.csv(x, file.path("results/tables", name), row.names = FALSE)
save_fig <- function(g, name, w, h) ggsave(file.path("results/figures", name), g, width = w, height = h)

failed <- aggregate(is.na(p) ~ learner, cv_pred, sum)
if (any(failed[[2]] > 0)) { message("CV predictions missing (failed fits):"); print(failed) }

# ---- Table 1: characteristics by primary outcome -----------------------------

describe <- function(x, g) {
  out <- function(v) {
    if (all(v[!is.na(v)] %in% c(0, 1))) {
      sprintf("%d (%s)", sum(v == 1, na.rm = TRUE), fmt(100 * mean(v == 1, na.rm = TRUE), 1))
    } else {
      q <- quantile(v, c(0.5, 0.25, 0.75), na.rm = TRUE)
      sprintf("%s (%s-%s)", fmt(q[1], 1), fmt(q[2], 1), fmt(q[3], 1))
    }
  }
  c(all = out(x), no = out(x[g == "no"]), yes = out(x[g == "yes"]),
    missing = sprintf("%d (%s)", sum(is.na(x)), fmt(100 * mean(is.na(x)), 1)))
}
rows <- list()
for (v in c(numeric_vars, binary_vars, "glycaemia_gl", process_vars)) {
  rows[[v]] <- c(variable = unname(labels[v]), describe(d[[v]], d$y))
}
for (lev in levels(d$education)) {
  rows[[paste0("education_", lev)]] <- c(
    variable = paste("Education:", lev),
    describe(as.numeric(d$education == lev), d$y))
}
table1 <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(table1) <- c("variable",
                   sprintf("all (n = %d)", nrow(d)),
                   sprintf("no complication (n = %d)", sum(d$y == "no")),
                   sprintf("complication (n = %d)", sum(d$y == "yes")),
                   "missing, n (%)")
write_table(table1, "table1-characteristics.csv")

# ---- Table 2: repeated cross-validation performance --------------------------
# metrics per repeat on pooled out-of-fold risks; mean and 2.5-97.5 percentiles

per_rep <- function(pred, by) {
  s <- split(pred, pred[by], drop = TRUE)
  do.call(rbind, lapply(s, function(x) {
    x <- x[!is.na(x$p), ]
    cbind(x[1, by, drop = FALSE], performance(x$y, x$p))
  }))
}
cv_rep <- per_rep(cv_pred, c("learner", "rep"))
cv_sum <- summarise_repeats(cv_rep)
cv_sum <- cv_sum[match(simplicity_order, cv_sum$learner), ]
metrics <- c("auc", "citl", "slope", "oe", "brier", "scaled_brier", "ece", "logloss")
show_ci <- function(s, k = 3) {
  out <- data.frame(learner = learner_names[s$learner])
  for (m in metrics) out[[m]] <- sprintf("%s (%s to %s)", fmt(s[[m]], k),
                                         fmt(s[[paste0(m, "_lo")]], k),
                                         fmt(s[[paste0(m, "_hi")]], k))
  out
}

# risk per woman averaged over the 20 repeats: curves and DeLong interval
oof <- aggregate(p ~ learner + row, data = cv_pred, FUN = mean)
oof$y <- d$y[oof$row]
oof_by <- split(oof, oof$learner)[simplicity_order]
delong <- do.call(rbind, lapply(names(oof_by), function(l) {
  a <- auc_ci(oof_by[[l]]$y, oof_by[[l]]$p, positive = "yes")
  data.frame(learner = l, auc_delong = sprintf("%s (%s to %s)", fmt(a$auc, 3),
                                               fmt(a$conf_low, 3), fmt(a$conf_high, 3)))
}))
table2 <- show_ci(cv_sum)
table2$auc_delong_mean_risk <- delong$auc_delong
table2$selected <- ifelse(simplicity_order == selected, "final",
                          ifelse(simplicity_order == fm$rule_selected, "rule", ""))
write_table(table2, "table2-cv-performance.csv")

specs <- aggregate(spec ~ learner, cv_spec, function(s) {
  t <- sort(table(s), decreasing = TRUE)
  sprintf("%s (%d of %d folds)", names(t)[1], t[1], length(s))
})
specs$learner <- learner_names[specs$learner]
write_table(specs, "tableS-tuned-specs.csv")

# ---- Figure 1: funcml learner comparison --------------------------------------

cmp_metrics <- c(auc = "AUC", slope = "Calibration slope", scaled_brier = "Scaled Brier",
                 ece = "ECE")
cmp <- do.call(rbind, lapply(names(cmp_metrics), function(m) data.frame(
  model = learner_names[cv_sum$learner], metric = cmp_metrics[[m]],
  mean = cv_sum[[m]], conf_low = cv_sum[[paste0(m, "_lo")]],
  conf_high = cv_sum[[paste0(m, "_hi")]])))
cmp$metric <- factor(cmp$metric, levels = cmp_metrics)
# same layout as plot.funcml_compare, learners on the y axis for legibility
cmp$model <- factor(cmp$model, levels = rev(learner_names[simplicity_order]))
save_fig(ggplot(cmp, aes(mean, model)) +
           geom_linerange(aes(xmin = conf_low, xmax = conf_high), colour = "#0072B2",
                          linewidth = 0.8) +
           geom_point(size = 1.8) +
           facet_wrap(~metric, nrow = 1, scales = "free_x") +
           labs(x = "Mean over 20 repeats (2.5-97.5 percentiles)", y = NULL) +
           theme_funcml() + theme(axis.text.x = element_text(size = 7)),
         "fig1-compare.pdf", 11, 3.8)

# ---- Figures 2-4: ROC, calibration, decision curve ----------------------------
# funcml objects built from out-of-fold risks (not in-sample interpret() calls)

as_calibration <- function(y, p, bins = 10) {
  structure(list(method = "calibration",
                 result = list(curve = calibration_curve(y, p, bins = bins, positive = "yes"),
                               prob = p, truth = y, positive = "yes",
                               ece = ece(y, p, bins = bins, positive = "yes"),
                               mce = mce(y, p, bins = bins, positive = "yes")),
                 diagnostics = list(bins = bins, strategy = "quantile", positive = "yes")),
            class = "funcml_calibration")
}
thr <- seq(0.01, 0.50, by = 0.01)
as_dca <- function(y, p) {
  structure(list(method = "dca",
                 result = list(curve = dca(y, p, positive = "yes", thresholds = thr),
                               positive = "yes", prevalence = mean(y == "yes"))),
            class = "funcml_dca")
}
rocs <- lapply(oof_by, function(x) roc_curve(x$y, x$p, positive = "yes"))
cals <- lapply(oof_by, function(x) as_calibration(x$y, x$p))
dcas <- lapply(oof_by, function(x) as_dca(x$y, x$p))
write_table(dcas[[selected]]$result$curve, "tableS-dca-final.csv")
panel_title <- function(l) labs(title = learner_names[[l]])

roc_all <- do.call(rbind, lapply(names(rocs), function(l) data.frame(
  learner = sprintf("%s (%.3f)", learner_names[[l]], rocs[[l]]$auc),
  fpr = 1 - rocs[[l]]$curve$specificity, tpr = rocs[[l]]$curve$sensitivity)))
roc_all <- roc_all[order(roc_all$learner, roc_all$fpr, roc_all$tpr), ]
g_roc_all <- ggplot(roc_all, aes(fpr, tpr, colour = learner)) +
  geom_abline(linetype = 2, colour = "grey70") + geom_step(linewidth = 0.5) +
  coord_equal() + labs(x = "1 - Specificity", y = "Sensitivity", colour = "Learner (AUC)",
                       title = "All learners") + theme_funcml()
save_fig(plot(rocs[[selected]]) + labs(title = sprintf("%s (AUC = %.3f)",
           learner_names[[selected]], rocs[[selected]]$auc)) | g_roc_all,
         "fig2-roc.pdf", 11, 4.8)
save_fig(wrap_plots(lapply(names(rocs), function(l) plot(rocs[[l]]) +
           labs(title = sprintf("%s (%.3f)", learner_names[[l]], rocs[[l]]$auc))), ncol = 4),
         "figS-roc-learners.pdf", 12, 9)

save_fig(plot(cals[[selected]]) + panel_title(selected) |
           plot(cals[[selected]], style = "histogram") + panel_title(selected),
         "fig3-calibration.pdf", 9, 4.5)
save_fig(wrap_plots(lapply(names(cals), function(l) plot(cals[[l]]) + panel_title(l)), ncol = 4),
         "figS-calibration-learners.pdf", 12, 10)

save_fig(plot(dcas[[selected]]) + panel_title(selected) + coord_cartesian(ylim = c(-0.02, 0.14)),
         "fig4-dca.pdf", 5.5, 4.2)
save_fig(wrap_plots(lapply(names(dcas), function(l) plot(dcas[[l]]) + panel_title(l) +
           coord_cartesian(ylim = c(-0.02, 0.14)) + theme(legend.position = "none")),
           ncol = 4) + plot_layout(guides = "collect"),
         "figS-dca-learners.pdf", 12, 9)

# ---- Table 3: classification at fixed thresholds (selected learner) -----------

sel <- cv_pred[cv_pred$learner == selected & !is.na(cv_pred$p), ]
class_rows <- list()
for (t in c(0.10, 0.15, 0.20)) {
  r <- do.call(rbind, lapply(split(sel, sel$rep), function(x) {
    pos <- x$p >= t; yy <- x$y == "yes"
    data.frame(flagged = mean(pos), sensitivity = mean(pos[yy]),
               specificity = mean(!pos[!yy]), ppv = mean(yy[pos]),
               npv = mean(!yy[!pos]),
               net_benefit = mean(pos & yy) - mean(pos & !yy) * t / (1 - t))
  }))
  class_rows[[as.character(t)]] <- data.frame(
    threshold = t, t(sapply(r, function(v) sprintf("%s (%s to %s)",
      fmt(mean(v, na.rm = TRUE), 3), fmt(quantile(v, 0.025, na.rm = TRUE), 3),
      fmt(quantile(v, 0.975, na.rm = TRUE), 3)))))
}
write_table(do.call(rbind, class_rows), "table3-classification.csv")

# ---- Supplementary: sensitivity analyses -------------------------------------

sens_rep <- per_rep(sens_pred, c("analysis", "learner", "rep"))
sens_sum <- do.call(rbind, lapply(split(sens_rep, sens_rep$analysis), function(s) {
  r <- summarise_repeats(s[setdiff(names(s), "analysis")])
  r <- r[match(intersect(simplicity_order, r$learner), r$learner), ]
  cbind(analysis = s$analysis[1], show_ci(r))
}))
sens_n <- tapply(sens_pred$row, sens_pred$analysis, function(r) length(unique(r)))
sens_e <- sapply(split(sens_pred, sens_pred$analysis),
                 function(x) sum(x$y[!duplicated(x$row)] == "yes"))
sens_sum$n <- sens_n[sens_sum$analysis]
sens_sum$events <- sens_e[sens_sum$analysis]
# k-nearest neighbors can predict risks of exactly 0 or 1, so the logit offset of
# calibration-in-the-large is infinite (S6): reported as not estimable
blown <- abs(as.numeric(sub(" .*", "", sens_sum$citl))) > 100
sens_sum$citl[blown] <- "not estimable"
write_table(sens_sum, "tableS-sensitivity.csv")

# ---- Figures 5-6: funcml interpret() plots (final model) ----------------------

relabel <- function(x) {
  lab <- function(f) ifelse(f %in% names(labels), labels[f], f)
  # local SHAP and LIME print "variable = value" strings
  lab_df <- function(z) {
    for (cc in intersect(c("feature_label", "feature.value"), names(z)))
      z[[cc]] <- paste0(lab(z$feature), " = ", sub("^.*? = ", "", z[[cc]]))
    z$feature <- lab(z$feature)
    z
  }
  r <- x$result
  if (is.data.frame(r) && "feature" %in% names(r)) r <- lab_df(r)
  if (is.list(r) && !is.data.frame(r)) for (k in names(r)) {
    if (is.data.frame(r[[k]]) && "feature" %in% names(r[[k]])) r[[k]] <- lab_df(r[[k]])
    if (is.data.frame(r[[k]])) for (cc in intersect(c("feature_x", "feature_y"), names(r[[k]])))
      r[[k]][[cc]] <- lab(r[[k]][[cc]])
    if (is.matrix(r[[k]]) && !is.null(colnames(r[[k]]))) colnames(r[[k]]) <- lab(colnames(r[[k]]))
  }
  x$result <- r
  x
}
vip <- relabel(interp$vip); shap <- relabel(interp$shap)
permute <- relabel(interp$permute)

imp <- Reduce(function(a, b) merge(a, b, by = "feature"), list(
  setNames(interp$vip$result$scores[c("feature", "importance")],
           c("feature", "permutation_auc_loss")),
  setNames(interp$permute$result$scores[c("feature", "importance")],
           c("feature", "permutation_brier_ratio")),
  setNames(aggregate(abs(shap) ~ feature, interp$shap$result, mean),
           c("feature", "mean_abs_shap"))))
imp <- merge(imp, setNames(interp$interaction$result$results[c("feature", "interaction")],
                           c("feature", "friedman_h")), by = "feature", all.x = TRUE)
imp <- imp[order(-imp$permutation_auc_loss), ]
write_table(imp, "tableS-importance.csv")

save_fig(plot(vip) + labs(title = "Permutation importance (loss in AUC)") |
           plot(shap) + labs(title = "SHAP values (risk scale)"),
         "fig5-importance.pdf", 13, 6)

# ---- Supplementary interpretation figures (remaining interpret() methods) ------

save_fig(plot(relabel(interp$ice)) + labs(title = "Individual conditional expectation") |
           plot(relabel(interp$profile)) + labs(title = "Profiles"),
         "figS-ice-profile.pdf", 13, 6)
save_fig(plot(relabel(interp$interaction)) + labs(title = "Interaction strength (Friedman H)") |
           plot(interp$surrogate) + labs(title = sprintf("Global surrogate tree (fidelity %.2f)",
                                                         interp$surrogate$result$fidelity)),
         "figS-interaction-surrogate.pdf", 14, 6)
ex_title <- function(k) with(interp$examples[interp$examples$example == k, ],
  sprintf("%s risk woman (predicted %.2f, outcome %s)", c(p10 = "Low", p50 = "Median",
          p90 = "High")[[k]], risk, y))
save_fig(wrap_plots(unlist(lapply(names(interp$local), function(k) {
  l <- interp$local[[k]]
  list(plot(relabel(l$shap)) + labs(title = ex_title(k), subtitle = "SHAP"),
       plot(relabel(l$lime)) + labs(title = " ", subtitle = sprintf(
         "Local surrogate (LIME), fidelity %.2f", l$lime$result$fidelity)),
       plot(relabel(l$ceteris_paribus)) + labs(title = " ", subtitle = "Ceteris paribus"))
}), recursive = FALSE), ncol = 3, widths = c(1, 1, 1.4)), "figS-local.pdf", 17, 15)
save_fig(plot(interp$calibration) + labs(title = "Apparent calibration (final fit, all women)") |
           plot(interp$dca) + labs(title = "Apparent decision curve") +
           coord_cartesian(ylim = c(-0.02, 0.14)),
         "figS-apparent.pdf", 11, 4.5)

# ---- estimate(): standardised risk differences (binary predictors) -------------

jobs <- est$jobs[est$jobs$done, ]
est_tab <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(i) {
  e <- est$estimates[[paste(jobs$var[i], jobs$estimand[i], sep = ":")]]
  # IATE: median and quartiles of the individual differences
  ci <- if (jobs$estimand[i] == "IATE") quantile(e$effects$effect, c(0.5, 0.25, 0.75))
        else c(e$estimate, e$conf_int[["lower"]], e$conf_int[["upper"]])
  data.frame(var = jobs$var[i], estimand = jobs$estimand[i], estimate = ci[[1]],
             lower = ci[[2]], upper = ci[[3]], n_with = jobs$n_with[i])
}))
est_tab$predictor <- labels[est_tab$var]
est_out <- reshape(transform(est_tab, value = sprintf("%s (%s to %s)", fmt(estimate, 3),
                                                      fmt(lower, 3), fmt(upper, 3)))[
                     c("predictor", "estimand", "value")],
                   idvar = "predictor", timevar = "estimand", direction = "wide")
names(est_out) <- sub("^value[.]", "", names(est_out))
est_out <- merge(est_out, setNames(est_tab[est_tab$estimand == "ATE", c("predictor", "n_with")],
                                   c("predictor", "women_with")), by = "predictor", all.x = TRUE)
write_table(est_out, "tableS-estimate.csv")
ate <- est_tab[est_tab$estimand != "IATE", ]
ate$estimand <- factor(ate$estimand, levels = c("ATE", "ATT", "CATE_nulliparous", "CATE_age35"),
                       labels = c("ATE", "ATT", "CATE: nulliparous", "CATE: age 35+"))
ord <- ate[ate$estimand == "ATE", ]
ate$predictor <- factor(ate$predictor, levels = ord$predictor[order(ord$estimate)])
save_fig(ggplot(ate, aes(estimate, predictor)) +
           geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
           geom_pointrange(aes(xmin = lower, xmax = upper), colour = "#0072B2", size = 0.3) + facet_wrap(~estimand, nrow = 1) +
           labs(x = "Difference in predicted risk, predictor set to 1 vs 0 (bootstrap 95% CI)",
                y = NULL) + theme_funcml(),
         "fig6-estimate.pdf", 13, 5.5)

message("selected learner: ", selected, "; explained fit vs MI average r = ",
        round(interp$cor_mi, 3))
print(table2[c("learner", "auc", "slope", "citl", "scaled_brier", "ece")], row.names = FALSE)
