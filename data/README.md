# Data (not included)

The OUMOUMATI registry data contain individual clinical records and are not
shared in this repository. They are available from the OUMOUMATI investigators
on reasonable request and with the approval of the ethics committee.

To rerun the analysis, place the files here:

| Path | Content |
|---|---|
| `data/raw/BDD GAR AI <center>.xlsm` | One registry workbook per health center (four files) |
| `data/audit-2025/maternalrisk_data_rabat.csv` | Analysis file of the earlier (2025) analysis, used only by `R/05-audit-2025.R` |

`R/01-clean.R` writes the analysis data set to `data/processed/analysis.csv`.
Nothing under `data/` except this file is tracked by git.
