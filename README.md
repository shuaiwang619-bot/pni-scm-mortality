# PNI-SCM mortality analysis

This repository contains public reproducibility materials for a 72-hour landmark
analysis of whether septic cardiomyopathy (SCM) modifies the association between
early prognostic nutritional index (PNI) and mortality in adults with sepsis.

## What is included

- `code/`: R scripts used for cohort summaries, Cox interaction models,
  diagnostics, sensitivity analyses, component analyses, and figures.
- `code/imputation/`: MICE imputation scripts for authorized local use.
- `sql/`: MIMIC-IV landmark cohort extraction SQL.
- `results/tables/`: aggregate, non-patient-level result tables.
- `results/figures/`: compressed PNG versions of the manuscript figures.
- `docs/`: data access notes, variable definitions, analysis workflow, and
  upload manifest.

## What is not included

This public release intentionally does not include any patient-level datasets,
raw clinical files, imputed patient-level files, RDS/MICE objects, manuscript
submission documents, cover letters, proof responses, similarity reports, or
journal-formatted TIFF source files.

## Data access

The MIMIC-IV data are available to credentialed users through PhysioNet under
the MIMIC data use agreement. Derived patient-level MIMIC datasets cannot be
redistributed here.

The single-center Linyi exploratory cohort contains institutional clinical data
and is not publicly released because of privacy and ethics restrictions.

See `docs/data_access.md` for details.

## Reproducibility notes

All stochastic analysis code uses seed `2025`. The public tables and figures in
`results/` are aggregate outputs only.

Authorized users who wish to rerun patient-level analyses should place their
locally approved datasets in `data_private/` using the filenames documented in
`data_private/README.md`. That directory is ignored by Git and should never be
uploaded.

## Main model

The primary Cox models used:

```r
Surv(time, event) ~ pni_low5_c * scm + covariates
```

where `pni_low5 = -PNI / 5` and `pni_low5_c` is mean-centered. Therefore, the
PNI effect estimates are interpreted per 5-point lower PNI. Interaction
estimates below 1 indicate attenuation of the lower-PNI mortality association
among patients with SCM.

## License

No reuse license has been selected yet. Add a license before making a public
repository if redistribution or reuse terms should be explicit.
