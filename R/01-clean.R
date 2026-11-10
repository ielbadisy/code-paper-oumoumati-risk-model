# 01-clean.R
# Build the analysis dataset from the four per-centre GAR workbooks.
# Run from the repository root:  Rscript R/01-clean.R
# Input : data/raw/BDD GAR AI *.xlsm  (one workbook per health centre)
# Output: data/processed/analysis.csv   one row per woman
#         results/tables/flow.csv          participant flow
#         results/tables/outcome-components.csv
#         results/tables/missingness.csv
#         results/tables/hb-by-anaemia.csv  (audit of the 2025 model)
# Rules are the ones written in protocol.md, section 4-6.

suppressPackageStartupMessages({
  library(readxl)
  library(basetable)
})

raw_dir  <- "data/raw"
out_data <- "data/processed"
out_tab  <- "results/tables"
dir.create(out_data, showWarnings = FALSE, recursive = TRUE)
dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)

delivery_year <- 2024  # all deliveries in the registry took place in 2024

# ---- helpers ---------------------------------------------------------------

files <- list.files(raw_dir, pattern = "^BDD GAR AI .*\\.xlsm$", full.names = TRUE)
centre_of <- function(f) {
  x <- sub("^BDD GAR AI\\s+", "", tools::file_path_sans_ext(basename(f)))
  trimws(sub("\\(RECTIF\\)", "", x))
}

# one sheet from every workbook, stacked, all columns as text
read_sheet <- function(sheet) {
  parts <- lapply(files, function(f) {
    x <- suppressMessages(suppressWarnings(
      read_excel(f, sheet = sheet, col_types = "text")))
    x <- as.data.frame(x)
    x <- x[, !grepl("^\\.\\.\\.|^autres$", names(x)), drop = FALSE]
    names(x)[1] <- "id"
    x$centre <- centre_of(f)
    x[!is.na(x$id) & grepl("^[0-9]+$", trimws(x$id)), , drop = FALSE]
  })
  out <- do.call(rbindfill, parts)
  out$id <- as.integer(trimws(out$id))
  out
}

# first decimal number in a string ("11,3 g/dl" -> 11.3, "36SA+1J" -> 36)
firstnum <- function(x) {
  x <- sub(",", ".", as.character(x), fixed = TRUE)
  m <- regmatches(x, regexpr("[0-9]+(\\.[0-9]+)?", x))
  out <- rep(NA_real_, length(x))
  out[grepl("[0-9]", x)] <- as.numeric(m)
  out
}

# numbers with thousands separators ("200 000/mm3" -> 200000)
spacednum <- function(x) firstnum(gsub("(?<=[0-9])[ .](?=[0-9]{3}\\b)", "", x, perl = TRUE))

within_range <- function(x, lo, hi) ifelse(!is.na(x) & x >= lo & x <= hi, x, NA)

has <- function(x, pattern) {
  x <- lower(removeaccents(as.character(x)))
  !is.na(x) & grepl(pattern, x)
}

yesno <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ifelse(x == "OUI", 1L, ifelse(x == "NON", 0L, NA_integer_))
}

# ---- outcome sheet ---------------------------------------------------------

pred <- read_sheet("VARIABLES PREDECTIVES")
comp <- c(gdm = "pred_diabeteGestationnel", gest_htn = "pred_hag",
          preeclampsia = "pred_preeclampsie", eclampsia = "pred_eclampsie",
          haemorrhage = "pred_hemorragie", infection = "pred_infection",
          cardiopathy = "pred_cardiopathies", anaemia = "pred_anemie",
          scarred_uterus = "pred_uterusCicatriciel", preterm = "pred_premature")
y <- data.frame(centre = pred$centre, id = pred$id)
for (k in names(comp)) y[[paste0("c_", k)]] <- yesno(pred[[comp[[k]]]])
y <- y[!duplicated(y[c("centre", "id")]), ]

# primary composite: incident complications only (protocol 5.1)
primary <- c("c_gdm", "c_gest_htn", "c_preeclampsia", "c_eclampsia",
             "c_haemorrhage", "c_infection", "c_preterm")
y$y_primary <- as.integer(rowSums(y[primary], na.rm = TRUE) > 0)
y$y_primary[rowSums(!is.na(y[primary])) == 0] <- NA
y$y_hdp <- as.integer(rowSums(y[c("c_gest_htn", "c_preeclampsia", "c_eclampsia")],
                              na.rm = TRUE) > 0)
# 2025 composite (all ten flags), kept for the audit only
y$y_old <- as.integer(rowSums(y[paste0("c_", names(comp))], na.rm = TRUE) > 0)

# ---- identification --------------------------------------------------------

idt <- read_sheet("IDENTIFICATION")
age_raw <- firstnum(idt$ident_age)
age <- ifelse(!is.na(age_raw) & age_raw > 1900, delivery_year - age_raw, age_raw)
edu <- lower(removeaccents(idt$ident_niveauInstruction))
X <- data.frame(
  centre = idt$centre, id = idt$id,
  age = within_range(age, 13, 55),
  education = casewhen(list(
    none      = grepl("sans", edu),
    primary   = grepl("primaire", edu),
    middle    = grepl("college", edu),
    secondary = grepl("lycee|bac|univ|diplome", edu))),
  rural = ifelse(idt$ident_adresse %in% c("Rural", "Urbain"),
                 as.integer(idt$ident_adresse == "Rural"), NA_integer_),
  consanguinity = yesno(idt$ident_consanguinite),
  wanted_pregnancy = yesno(idt$ident_grossesseDesiree)
)
X <- X[!duplicated(X[c("centre", "id")]), ]

# ---- obstetric and medical history ----------------------------------------

acc <- read_sheet("INFOS ACCH ANTER")
acc <- acc[!duplicated(acc[c("centre", "id")]), ]
hist <- data.frame(
  centre = acc$centre, id = acc$id,
  gravidity = within_range(firstnum(acc$infaccant_gestation), 1, 15),
  parity    = within_range(firstnum(acc$infaccant_parite), 0, 12)
)
hist$parity[!is.na(hist$parity) & grepl("^-", acc$infaccant_parite)] <- NA
hist$gravidity[!is.na(hist$gravidity) & grepl("^-", acc$infaccant_gestation)] <- NA

flag_by_woman <- function(d, value, pattern) {
  hit <- has(d[[value]], pattern)
  agg <- unique(d[c("centre", "id")])
  key <- paste(d$centre, d$id)
  agg$flag <- as.integer(paste(agg$centre, agg$id) %in% key[hit])
  agg
}
add_flag <- function(base, d, value, pattern, name) {
  f <- flag_by_woman(d, value, pattern)
  names(f)[3] <- name
  out <- merge(base, f, by = c("centre", "id"), all.x = TRUE)
  out[[name]][is.na(out[[name]])] <- 0L  # absent from the sheet: none recorded
  out
}

obs <- read_sheet("ANTECEDENTS OBSTETRICAUX")
hist <- add_flag(hist, obs, "antobs_nature", "avortement", "prior_abortion")
hist <- add_flag(hist, obs, "antobs_nature", "mort|premature|hemorragie",
                 "prior_adverse_outcome")

dac <- read_sheet("DEROUL ACCH ANTER")
hist <- add_flag(hist, dac, "deraccant_modaliteExtraction", "cesarienne",
                 "prior_caesarean_hist")
hist <- add_flag(hist, dac, "deraccant_typeComplicationAccouchement",
                 "htag|eclam|toxemie", "prior_htag_dac")
hist <- add_flag(hist, dac, "deraccant_typeComplicationAccouchement",
                 "hemorr|hgie", "prior_haemorrhage")

fam <- read_sheet("ANTECEDENTS HEREDITAIRES")
hist <- add_flag(hist, fam, "anther_type", "diab", "family_diabetes")
hist <- add_flag(hist, fam, "anther_type", "hta", "family_hypertension")

per <- read_sheet("ANTECEDENTS FEMME")
hist <- add_flag(hist, per, "antf_type", "medicaux|htag|anemie",
                 "personal_medical_history")
hist <- add_flag(hist, per, "antf_type", "chirurg", "surgical_history")
hist <- add_flag(hist, per, "antf_type", "htag", "prior_htag_antf")
hist$prior_htag <- pmax(hist$prior_htag_dac, hist$prior_htag_antf)
hist$prior_htag_dac <- hist$prior_htag_antf <- NULL

# ---- current pregnancy (fixed at booking) ----------------------------------

gro <- read_sheet("INFOS GROSS ACTUELLE")
gro <- gro[!duplicated(gro[c("centre", "id")]), ]
h <- firstnum(gro$infgrsact_taille)
h <- ifelse(!is.na(h) & h < 3, h * 100, h)
cur <- data.frame(
  centre = gro$centre, id = gro$id,
  height_cm = within_range(h, 130, 190),
  rhesus_negative = ifelse(gro$infgrsact_rhesus %in% c("RH+", "RH-"),
                           as.integer(gro$infgrsact_rhesus == "RH-"), NA_integer_)
)

# ---- baseline clinical exam: earliest valid visit with each measure --------

ex <- read_sheet("EXAMENS CLINIQUES")
ex$month <- within_range(suppressWarnings(as.numeric(ex$excl_mois)), 1, 10)
ga <- firstnum(ex$excl_ageGrossesse)
gd <- ifelse(grepl("\\+", ex$excl_ageGrossesse),
             firstnum(sub(".*\\+", "", ex$excl_ageGrossesse)), 0)
ex$ga_weeks <- within_range(ga + ifelse(is.na(gd), 0, gd) / 7, 4, 43)
ex$weight <- within_range(firstnum(ex$excl_poids), 35, 150)
bp <- regmatches(ex$excl_ta, regexec("^\\s*([0-9]+)\\s*/\\s*([0-9]+)", ex$excl_ta))
sbp <- sapply(bp, function(m) if (length(m) == 3) as.numeric(m[2]) else NA)
dbp <- sapply(bp, function(m) if (length(m) == 3) as.numeric(m[3]) else NA)
ex$sbp <- within_range(ifelse(sbp < 30, sbp * 10, sbp), 70, 220)  # cmHg -> mmHg
ex$dbp <- within_range(ifelse(dbp < 30, dbp * 10, dbp), 40, 140)
ex$sbp[is.na(ex$dbp)] <- NA
ex$dbp[is.na(ex$sbp)] <- NA
ex$order <- ifelse(is.na(ex$ga_weeks), ex$month * 4.3, ex$ga_weeks)

first_measure <- function(d, value, name) {
  d <- d[!is.na(d[[value]]) & !is.na(d$order), c("centre", "id", "order", value)]
  d <- firstby(d, by = c("centre", "id"), order = "order")
  d <- as.data.frame(d)
  names(d)[names(d) == value] <- name
  names(d)[names(d) == "order"] <- paste0(name, "_ga")
  d
}
base <- Reduce(function(a, b) merge(a, b, by = c("centre", "id"), all = TRUE), list(
  first_measure(ex, "weight", "weight_kg"),
  first_measure(transform(ex, sbp = sbp), "sbp", "sbp"),
  first_measure(transform(ex, dbp = dbp), "dbp", "dbp")
))
base$dbp_ga <- NULL
names(base)[names(base) == "weight_kg_ga"] <- "booking_ga_weeks"

# findings of the first recorded clinical exam (protocol 4.2)
ex1 <- as.data.frame(firstby(ex[!is.na(ex$order), ], by = c("centre", "id"), order = "order"))
conj <- lower(removeaccents(ex1$excl_etatConjonctives))
spec <- lower(removeaccents(ex1$excl_speculum))
exam1 <- data.frame(
  centre = ex1$centre, id = ex1$id,
  pallor = ifelse(grepl("decolor", conj), 1L, ifelse(grepl("normal", conj), 0L, NA_integer_)),
  vaginal_discharge = as.integer(!is.na(spec) & grepl("perte|leucor|fetid|candid", spec))
)

# ---- baseline biology -------------------------------------------------------

bio <- read_sheet("EXAMENS BIOLOGIQUES")
bio$month <- within_range(suppressWarnings(as.numeric(bio$exbio_mois)), 1, 10)
bio$order <- bio$month
bio$hb <- within_range(firstnum(bio$exbio_hemoglobine), 5, 18)
g <- firstnum(bio$exbio_bilanGlycemique)
g <- ifelse(!is.na(g) & g > 20, g / 100, g)  # mg/dL -> g/L
bio$glycaemia <- within_range(g, 0.4, 3)
p <- spacednum(bio$exbio_plaquettes)
p <- ifelse(!is.na(p) & p > 2000, p / 1000, p)  # /mm3 -> 10^3/mm3
bio$platelets <- within_range(p, 20, 800)
biob <- Reduce(function(a, b) merge(a, b, by = c("centre", "id"), all = TRUE), list(
  first_measure(bio, "hb", "hb_gdl"),
  first_measure(bio, "glycaemia", "glycaemia_gl"),
  first_measure(bio, "platelets", "platelets_k")
))
names(biob) <- sub("_ga$", "_month", names(biob))

# care-process markers from the first biology record (protocol 4.3, S5 only);
# records without a valid month are taken last
bio$order1 <- ifelse(is.na(bio$month), 99, bio$month)
bio1 <- as.data.frame(firstby(bio, by = c("centre", "id"), order = "order1"))
syph <- lower(removeaccents(bio1$exbio_syphilis))
rub <- lower(removeaccents(bio1$exbio_rubeole))
proc <- data.frame(
  centre = bio1$centre, id = bio1$id,
  serology_done = as.integer(!is.na(syph) & !grepl("^nr", syph)),
  rubella_immune = ifelse(grepl("^\\+|posit|imm", rub), 1L,
                          ifelse(grepl("^-|neg|absence", rub), 0L, NA_integer_))
)

# ---- assemble ---------------------------------------------------------------

keys <- c("centre", "id")
d <- Reduce(function(a, b) merge(a, b, by = keys, all.x = TRUE),
            list(X, hist, cur, base, exam1, biob, proc))
d <- merge(d, y, by = keys, all.x = TRUE)
d$bmi <- round(d$weight_kg / (d$height_cm / 100)^2, 1)
d$bmi <- within_range(d$bmi, 13, 60)
d$nulliparous <- as.integer(d$parity == 0)
d$prior_caesarean <- as.integer(d$prior_caesarean_hist == 1 |
                                (!is.na(d$c_scarred_uterus) & d$c_scarred_uterus == 1))
d$prior_caesarean_hist <- NULL
d$serology_done[is.na(d$serology_done)] <- 0L  # no biology record: none done

# ---- participant flow and exclusions ---------------------------------------

flow <- data.frame(step = "registered women", n = nrow(d))
d <- d[!is.na(d$y_primary), ]
flow <- rbind(flow, data.frame(step = "with outcome record (analysis set)", n = nrow(d)))
flow <- rbind(flow, data.frame(step = "primary composite events", n = sum(d$y_primary)))
d <- orderrows(d, by = keys)

# ---- summaries ---------------------------------------------------------------

cc <- grep("^c_|^y_", names(d), value = TRUE)
comp_tab <- data.frame(component = cc,
                       events = sapply(d[cc], function(v) sum(v == 1, na.rm = TRUE)),
                       n = sapply(d[cc], function(v) sum(!is.na(v))))
comp_tab$pct <- round(100 * comp_tab$events / comp_tab$n, 1)

miss <- data.frame(variable = names(d),
                   n_missing = sapply(d, function(v) sum(is.na(v))))
miss$pct <- round(100 * miss$n_missing / nrow(d), 1)

# audit: haemoglobin against the anaemia flag (explains the 2025 Hb signal)
hb <- d[!is.na(d$hb_gdl) & !is.na(d$c_anaemia), ]
hb_audit <- do.call(rbind, lapply(split(hb$hb_gdl, hb$c_anaemia), function(v)
  data.frame(n = length(v), median = median(v), q1 = quantile(v, 0.25),
             q3 = quantile(v, 0.75), min = min(v), max = max(v))))
hb_audit <- cbind(c_anaemia = rownames(hb_audit), hb_audit)

btwrite(d, file.path(out_data, "analysis.csv"))
btwrite(flow, file.path(out_tab, "flow.csv"))
btwrite(comp_tab, file.path(out_tab, "outcome-components.csv"))
btwrite(miss, file.path(out_tab, "missingness.csv"))
btwrite(hb_audit, file.path(out_tab, "hb-by-anaemia.csv"))

print(flow)
print(comp_tab, row.names = FALSE)
print(miss[miss$n_missing > 0, ], row.names = FALSE)
print(hb_audit)
