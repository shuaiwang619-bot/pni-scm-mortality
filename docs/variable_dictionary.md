# Variable dictionary

## Core exposures

- `pni`: Prognostic nutritional index, calculated as albumin in g/L plus
  5 times lymphocyte count in 10^9/L.
- `pni_low5`: Lower-PNI coding, calculated as `-pni / 5`.
- `pni_low5_c`: Mean-centered `pni_low5`; used in interaction models.
- `scm`: Septic cardiomyopathy indicator.

## SCM definition

In MIMIC-IV, SCM was defined from pre-landmark TTE-derived LVEF as either:

- any eligible LVEF less than 50%, or
- an LVEF decrease greater than 10 percentage points with post-decline LVEF
  less than 55%.

In the single-center exploratory cohort, SCM was defined using available
LVEF less than 50% because serial eligible LVEF measurements were not
consistently available.

## Time origin and landmark

- Index time: ICU admission.
- Landmark time: ICU admission plus 72 hours.
- MIMIC-IV PNI window: ICU admission -24 hours to +48 hours.
- MIMIC-IV SCM ascertainment window: ICU admission -24 hours to +72 hours.
- Follow-up for the primary time-to-event analyses began at the 72-hour
  landmark.

## Outcomes

- `event_28d`, `time_28d`: 28-day mortality event and follow-up time after the
  72-hour landmark.
- `event_90d`, `time_90d`: 90-day mortality event and follow-up time after the
  72-hour landmark.
- Exploratory external endpoint: in-hospital death endpoint as available in the
  single-center cohort.

## Main covariates

The fully adjusted MIMIC-IV model included age, sex, BMI, SOFA score,
log-transformed lactate, log-transformed creatinine, coronary artery disease,
diabetes mellitus, hypertension, chronic kidney disease, COPD, mechanical
ventilation, and vasopressor use.

