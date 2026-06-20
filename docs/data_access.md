# Data access and privacy

## MIMIC-IV

The MIMIC-IV source data are publicly accessible through PhysioNet only to
credentialed users who complete the required training and data use agreement.
This repository provides SQL and analysis code but does not redistribute
patient-level MIMIC-IV data or derived patient-level CSV files.

Users with approved MIMIC-IV access can run `sql/mimic_landmark72_extraction.sql`
inside their authorized MIMIC-IV environment and then continue the local
analysis workflow.

## Single-center exploratory cohort

The Linyi exploratory cohort contains institutional clinical data. Patient-level
records, raw Excel files, imputed files, and MICE objects are not released here
because of privacy and ethics restrictions.

Only aggregate results are included in `results/tables/` and `results/figures/`.

## Files intentionally excluded from GitHub

- Patient identifiers and local hospital records.
- Raw or imputed patient-level CSV/XLSX files.
- MICE `mids` objects and other RDS/RData files.
- Journal submission files, cover letters, proof responses, and similarity
  reports.
- Large journal-formatted TIFF files.

