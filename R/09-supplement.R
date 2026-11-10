# 09-supplement.R
# LaTeX tables of Additional file 1 from the result tables (no new analysis).
# Run from the repository root:  Rscript R/09-supplement.R
# Input : results/tables/*.csv
# Output: results/supplement-tables/*.tex

suppressPackageStartupMessages(library(knitr))
source("R/utils.R")

out <- "results/supplement-tables"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) read.csv(file.path("results/tables", f), check.names = FALSE,
                           stringsAsFactors = FALSE)
lab <- function(f) ifelse(f %in% names(labels), labels[f], f)
us <- function(x) {
  x <- gsub("aemorrhage", "emorrhage", x); x <- gsub("aemia", "emia", x)
  x <- gsub("aesarean", "esarean", x); x <- gsub("ycaemia", "ycemia", x)
  gsub("centre", "center", x)
}
tex <- function(d, file, align = NULL, digits = 3, font = "\\footnotesize", stack = FALSE) {
  d[] <- lapply(d, function(v) if (is.character(v)) us(v) else v)
  # stack = TRUE puts each interval under its estimate (wide interval tables)
  if (stack) d[] <- lapply(d, function(v)
    if (is.character(v)) sub("^(\\S+) \\((.*)\\)$", "\\\\makecell{\\1\\\\\\\\(\\2)}", v) else v)
  k <- kable(d, format = "latex", booktabs = TRUE, align = align, digits = digits,
             linesep = "", escape = !stack, row.names = FALSE)
  k <- gsub("10\\\\textasciicircum\\{\\}3", "$10^3$", k)
  k <- gsub("kg/m2", "kg/m$^2$", k, fixed = TRUE)
  writeLines(c(font, k), file.path(out, file))
}
ci_to <- function(x) gsub(" to ", ", ", x)

# ---- flow, outcome components, missingness -----------------------------------

comp_lab <- c(c_gdm = "Gestational diabetes", c_gest_htn = "Gestational hypertension",
              c_preeclampsia = "Pre-eclampsia", c_eclampsia = "Eclampsia",
              c_haemorrhage = "Hemorrhage", c_infection = "Infection",
              c_cardiopathy = "Cardiopathy", c_anaemia = "Anemia (not counted)",
              c_scarred_uterus = "Scarred uterus (not counted)", c_preterm = "Preterm birth",
              y_primary = "Primary composite outcome",
              y_hdp = "Hypertensive disorders (secondary outcome)",
              y_old = "Composite of the earlier analysis")
oc <- rd("outcome-components.csv")
oc <- data.frame(Item = comp_lab[oc$component], `Women, n` = oc$events,
                 `%` = sprintf("%.1f", oc$pct), check.names = FALSE)
tex(oc, "outcome-components.tex", c("l", "r", "r"))

fl <- rd("flow.csv")
tex(data.frame(Step = fl$step, n = fl$n), "flow.tex", c("l", "r"))

mi <- rd("missingness.csv")
mi <- mi[mi$variable %in% c(predictors, "height_cm", "weight_kg", "glycaemia_gl") &
           mi$n_missing > 0, ]
mi$variable <- lab(mi$variable)
mi$variable[mi$variable == "height_cm"] <- "Height (cm)"
mi$variable[mi$variable == "weight_kg"] <- "Weight (kg)"
mi <- mi[order(-mi$n_missing), ]
tex(data.frame(Variable = mi$variable, `Missing, n` = mi$n_missing,
               `%` = sprintf("%.1f", mi$pct), check.names = FALSE),
    "missingness.tex", c("l", "r", "r"))

# ---- center profile and sample size --------------------------------------------

cp <- rd("tableS-centre-profile.csv")
names(cp) <- c("Item", "MOUAZZA", "SAFSAF", "Sidi Abderrazak", "Sidi Yahya El Gharb", "All")
tex(cp, "center-profile.tex", c("l", rep("r", 5)), digits = 1)

ss <- rd("sample-size.csv")
ss$value <- ifelse(ss$value >= 100, format(ss$value, big.mark = ","), as.character(ss$value))
tex(setNames(ss, c("Quantity", "Value")), "sample-size.tex", c("l", "r"))

# ---- tuned hyperparameters -----------------------------------------------------

sp <- rd("tableS-tuned-specs.csv")
sp$spec <- sub("^structure\\(list\\(\\), names = character\\(0\\)\\)", "no tuning (defaults)",
               sp$spec)
sp$spec <- gsub("^list\\((.*)\\) \\(", "\\1 (", sp$spec)
sp$spec <- gsub("([0-9]+)L\\b", "\\1", sp$spec)
sp$spec <- gsub("0.0316227766016838", "0.032", sp$spec, fixed = TRUE)
sp$spec <- sub("\\((\\d+) of 100 folds\\)", "[\\1 of 100 folds]", sp$spec)
tex(setNames(sp, c("Learner", "Most frequent tuned setting [folds selecting it]")),
    "tuned-specs.tex", c("l", "l"))

# ---- sensitivity analyses S1-S7 -------------------------------------------------

an_lab <- c(S1_glycaemia = "S1 glycemia added", S2_early_labs = "S2 labs months 1-4",
            S3_booked_le20 = "S3 booked by 20 weeks", S5_process = "S5 care-process markers",
            S6_gdm_hdp = "S6 GDM or HDP outcome", S7_no_rural = "S7 no rural residence")
se <- rd("tableS-sensitivity.csv")
for (a in names(an_lab)) {
  z <- se[se$analysis == a, ]
  tex(data.frame(Learner = z$learner, AUC = ci_to(z$auc), Slope = ci_to(z$slope),
                 CITL = ci_to(z$citl), `O/E` = ci_to(z$oe),
                 `Scaled Brier` = ci_to(z$scaled_brier), ECE = ci_to(z$ece),
                 check.names = FALSE),
      paste0("sens-", tolower(sub("_.*", "", a)), ".tex"), c("l", rep("c", 6)),
      font = "\\footnotesize\\setlength{\\tabcolsep}{2.5pt}", stack = TRUE)
}
sm <- rd("tableS-smote.csv")
sm$learner <- c(bart = "BART", glmnet = "Elastic net", ranger = "Random forest")[sm$learner]
sm$analysis <- ifelse(sm$analysis == "S4_smote", "SMOTE", "None")
tex(data.frame(Learner = sm$learner, Rebalancing = sm$analysis, AUC = ci_to(sm$auc),
               Slope = ci_to(sm$slope), CITL = ci_to(sm$citl), `O/E` = ci_to(sm$oe),
               `Scaled Brier` = ci_to(sm$scaled_brier), `Mean risk` = ci_to(sm$mean_risk),
               check.names = FALSE),
    "smote.tex", c("l", "l", rep("c", 6)), font = "\\footnotesize\\setlength{\\tabcolsep}{2pt}", stack = TRUE)

# ---- combined importance, apparent decision curve ---------------------------------

im <- rd("tableS-importance.csv")
tex(data.frame(Predictor = lab(im$feature), `Loss in AUC` = sprintf("%.4f", im$permutation_auc_loss),
               `Brier ratio` = sprintf("%.3f", im$permutation_brier_ratio),
               `Mean |SHAP|` = sprintf("%.4f", im$mean_abs_shap),
               `Friedman H` = ifelse(is.na(im$friedman_h), "--", sprintf("%.2f", im$friedman_h)),
               check.names = FALSE),
    "importance.tex", c("l", rep("r", 4)))

dc <- rd("tableS-dca-final.csv")
dc <- dc[round(dc$threshold * 100) %% 5 == 0 & dc$threshold <= 0.5, ]
tex(data.frame(Threshold = sprintf("%.2f", dc$threshold), Model = sprintf("%.3f", dc$model),
               `Treat all` = sprintf("%.3f", dc$treat_all), `Treat none` = "0",
               check.names = FALSE),
    "dca-final.tex", rep("r", 4))

# ---- audit of the earlier analysis --------------------------------------------------

ac <- rd("audit-composite.csv")
tex(data.frame(Component = ac$component, n = ac$n, `In earlier composite` = ac$in_2025_composite,
               `In corrected composite` = ac$in_corrected_composite,
               `Only event (earlier)` = ifelse(is.na(ac$only_event_in_2025), "",
                                               ac$only_event_in_2025),
               `%` = ac$percent, check.names = FALSE),
    "audit-composite.tex", c("l", "r", "c", "c", "r", "r"), digits = 1)
al <- rd("audit-leakage.csv")
tex(data.frame(Predictor = al$predictor,
               `AUC, earlier outcome` = sprintf("%.3f", al$auc_2025_outcome),
               `AUC, corrected outcome` = ifelse(is.na(al$auc_corrected_outcome), "--",
                                                 sprintf("%.3f", al$auc_corrected_outcome)),
               check.names = FALSE),
    "audit-leakage.tex", c("l", "r", "r"))
au <- rd("audit-upsampling.csv")
tex(data.frame(Training = au$training, AUC = au$auc, CITL = au$citl, Slope = au$slope,
               `O/E` = au$oe, `Scaled Brier` = au$scaled_brier,
               `Mean predicted` = au$mean_predicted, Observed = au$observed,
               check.names = FALSE),
    "audit-upsampling.tex", c("l", rep("r", 7)))
ab <- rd("audit-label-baseline.csv")
tex(setNames(ab, c("Quantity", "Value")), "audit-baseline.tex", c("l", "r"))
hb <- rd("hb-by-anaemia.csv")
tex(data.frame(Anemia = ifelse(hb$c_anaemia == 1, "yes", "no"), n = hb$n,
               `Median (IQR) hemoglobin, g/dL` = sprintf("%.1f (%.1f-%.1f)", hb$median, hb$q1, hb$q3),
               Range = sprintf("%.1f-%.1f", hb$min, hb$max), check.names = FALSE),
    "hb-anemia.tex", c("l", "r", "r", "r"))
message("supplement tables written to ", out)
