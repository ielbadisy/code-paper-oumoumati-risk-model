# 06-centre-profile.R
# Outcome components and predictor recording by centre (protocol section 5.3):
# evidence for differential ascertainment across centres.
# Run from the repository root:  Rscript R/06-centre-profile.R
# Input : data/processed/analysis.csv (from 01-clean.R)
# Output: results/tables/tableS-centre-profile.csv

source("R/utils.R")

d <- load_analysis()
n_sum <- function(x) sum(x, na.rm = TRUE)

rows <- list(
  "Women"                          = function(x) nrow(x),
  "Primary composite"              = function(x) n_sum(x$y_primary),
  "Gestational diabetes"           = function(x) n_sum(x$c_gdm),
  "Gestational hypertension"       = function(x) n_sum(x$c_gest_htn),
  "Pre-eclampsia or eclampsia"     = function(x) n_sum(x$c_preeclampsia | x$c_eclampsia),
  "Haemorrhage"                    = function(x) n_sum(x$c_haemorrhage),
  "Infection"                      = function(x) n_sum(x$c_infection),
  "Preterm birth"                  = function(x) n_sum(x$c_preterm),
  "Rural residence"                = function(x) n_sum(x$rural),
  "Wanted pregnancy"               = function(x) n_sum(x$wanted_pregnancy),
  "Family history of diabetes"     = function(x) n_sum(x$family_diabetes),
  "Personal medical history"       = function(x) n_sum(x$personal_medical_history),
  "Booking GA, median (weeks)"     = function(x) round(median(x$booking_ga_weeks, na.rm = TRUE), 1),
  "Haemoglobin not recorded"       = function(x) n_sum(is.na(x$hb_gdl)),
  "Platelets not recorded"         = function(x) n_sum(is.na(x$platelets_k))
)

by_centre <- split(d, d$centre)
tab <- data.frame(item = names(rows), stringsAsFactors = FALSE)
for (cc in names(by_centre)) tab[[cc]] <- vapply(rows, function(f) f(by_centre[[cc]]), numeric(1))
tab$all <- vapply(rows, function(f) f(d), numeric(1))

write.csv(tab, "results/tables/tableS-centre-profile.csv", row.names = FALSE)
print(tab, row.names = FALSE)
