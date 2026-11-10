# Analysis code: booking-visit risk of pregnancy complications (OUMOUMATI)

R code for the paper **"Predicting pregnancy complications from booking-visit
data in Moroccan primary care: development and internal validation of a risk
model comparing eleven learners"** (submitted to _BMC Pregnancy and
Childbirth_).

The registry data are **not** included (see [`data/README.md`](data/README.md)).
The code reproduces every table and figure of the paper and its Additional
file 1 when the data are placed in `data/`.

## Scripts

Run from the repository root, in this order:

| Script | Purpose | Output |
|---|---|---|
| `R/00-sample-size.R` | Minimum sample size (Riley et al.) | `results/tables/sample-size.csv` |
| `R/01-clean.R` | Cleaning, booking values, outcome definition | `data/processed/analysis.csv`, flow, components, missingness |
| `R/02-model.R` | Eleven learners, 20 x 5-fold CV, nested imputation and tuning; S1-S3, S5-S7 | `results/cv-predictions.rds`, `results/sens-predictions.rds`, `results/cv-specs.rds` |
| `R/02b-final-model.R` | Final BART fit on all women (m = 5) | `results/final-model.rds` |
| `R/03-interpret.R` | Global and local explanations | `results/interpret.rds` |
| `R/03b-estimate.R` | Standardized risk differences (bootstrap) | `results/estimate.rds` |
| `R/04-results.R` | Tables and figures of the paper | `results/tables/`, `results/figures/` |
| `R/05-audit-2025.R` | Audit of the earlier analysis | `results/tables/audit-*.csv` |
| `R/06-centre-profile.R` | Recording by center | `results/tables/tableS-centre-profile.csv` |
| `R/07-smote.R` | Sensitivity analysis S4 (SMOTE) | `results/smote-predictions.rds` |
| `R/08-logistic.R` | Crude and adjusted odds ratios (Rubin's rules) | `results/tables/table4-logistic.csv` |
| `R/09-supplement.R` | LaTeX tables of Additional file 1 | `results/supplement-tables/` |
| `R/10-centre-only.R` | Sensitivity analysis S8: the centre as the only predictor | `results/tables/` |

```sh
Rscript R/00-sample-size.R
Rscript R/01-clean.R
OMP_NUM_THREADS=1 Rscript R/02-model.R
OMP_NUM_THREADS=1 Rscript R/02b-final-model.R
OMP_NUM_THREADS=1 Rscript R/03-interpret.R
OMP_NUM_THREADS=1 Rscript R/03b-estimate.R
Rscript R/04-results.R
OMP_NUM_THREADS=1 Rscript R/05-audit-2025.R
Rscript R/06-centre-profile.R
OMP_NUM_THREADS=1 Rscript R/07-smote.R
Rscript R/08-logistic.R
Rscript R/09-supplement.R
Rscript R/10-centre-only.R
```

`R/utils.R` holds the predictor list, labels, learner grids, and helpers, and
is sourced by every script. Model scripts parallelize across folds and should
be run with `OMP_NUM_THREADS=1`.

## Software

Results in the paper were obtained with R 4.5.1 and:

| Package | Version | Role |
|---|---|---|
| funcml | 0.9.5 | fitting, tuning, prediction, metrics, interpretation, figures |
| mimar | 1.4.0 | multiple imputation |
| basetable | 1.4.2 | data management |
| testflow | 1.0.0 | sample size |
| themis | 1.0.3 | SMOTE |
| glmnet | 5.0 | elastic net |
| earth | 5.3.5 | MARS |
| kknn | 1.4.1 | k-nearest neighbors |
| e1071 | 1.7-17 | support vector machine |
| ranger | 0.18.0 | random forest |
| xgboost | 3.2.1.1 | gradient boosting |
| lightgbm | 4.6.0 | gradient boosting |
| dbarts | 0.9-33 | BART |
| densemlp | 0.7.1 | multilayer perceptron |
| readxl, ggplot2, patchwork, knitr | | reading workbooks, figures, tables |

All are on CRAN. Older versions can be installed with
`remotes::install_version()`.

## License

Code: MIT (see `LICENSE.md`).
