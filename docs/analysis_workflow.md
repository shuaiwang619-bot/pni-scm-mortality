# Analysis workflow

## 1. Cohort extraction

Run `sql/mimic_landmark72_extraction.sql` in an authorized MIMIC-IV database
environment. The SQL defines the 72-hour landmark cohort, SCM ascertainment,
PNI calculation, outcomes, and baseline covariates.

## 2. Missing-data handling

Imputation scripts are provided in `code/imputation/` for local authorized use.
Exposure components, outcomes, follow-up times, and SCM status should not be
imputed for the public analysis workflow.

## 3. Main analyses

The main MIMIC-IV Cox interaction model is implemented in:

- `code/00_通用设置与函数_seed2025.R`
- `code/03_主模型分析_seed2025.R`

The interaction model estimates PNI-associated mortality gradients separately
for No SCM and SCM strata, and tests the interaction using likelihood-ratio
tests.

## 4. Diagnostics and sensitivity analyses

The repository includes scripts for:

- covariate adjustment sensitivity analyses,
- proportional-hazards diagnostics,
- VIF and collinearity diagnostics,
- natural cubic spline checks,
- PNI winsorization sensitivity,
- absolute-risk translation for PNI 30 versus PNI 40,
- PNI component comparisons, and
- exploratory external consistency analyses.

## 5. Public outputs

Aggregate tables are in `results/tables/`. Compressed display figures are in
`results/figures/`.

Do not commit `data_private/`, local SQL exports, or generated patient-level
output folders.

