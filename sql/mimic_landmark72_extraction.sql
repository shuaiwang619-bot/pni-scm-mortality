\set ON_ERROR_STOP on
\timing on
\pset pager off

/*
PNI-SCM mortality project: 72-hour landmark extraction.

Time origin:
  ICU admission (intime).

Primary SCM ascertainment:
  TTE-LVEF measured from ICU admission -24h to ICU admission +72h,
  bounded by the hospital admission/discharge interval.

Primary PNI:
  albumin(g/L) + 5 * lymphocyte(k/uL), using ICU admission -24h to +48h.
  No albumin-lymphocyte simultaneity rule is applied.

This archive version also extracts Table 1 baseline variables and modeling
transforms needed downstream:
  BMI, coronary artery disease, diabetes mellitus, hypertension, chronic kidney
  disease, COPD, WBC, platelets, BUN, pni_low5, albumin_low5, lymph_low05,
  log_lactate, and log_creatinine.
*/

CREATE INDEX IF NOT EXISTS lvef_measurements_subject_datetime_tte_idx
ON mimiciv_echo.lvef_measurements_20260531 (subject_id, measurement_datetime)
WHERE test_type = 'tte';

DROP TABLE IF EXISTS public.pni_scm_landmark72_20260605;

DROP TABLE IF EXISTS pg_temp.tmp_first_icu;
CREATE TEMP TABLE tmp_first_icu AS
SELECT
  ie.*,
  ROW_NUMBER() OVER (PARTITION BY ie.subject_id ORDER BY ie.intime, ie.stay_id) AS rn
FROM mimiciv_icu.icustays ie;

CREATE INDEX tmp_first_icu_stay_idx ON tmp_first_icu (stay_id);
CREATE INDEX tmp_first_icu_hadm_idx ON tmp_first_icu (hadm_id);
ANALYZE tmp_first_icu;

DROP TABLE IF EXISTS pg_temp.tmp_base_pre_exclusion;
CREATE TEMP TABLE tmp_base_pre_exclusion AS
SELECT
  fi.subject_id,
  fi.hadm_id,
  fi.stay_id,
  adm.admittime,
  adm.dischtime,
  adm.deathtime,
  adm.hospital_expire_flag,
  adm.admission_type,
  adm.admission_location,
  adm.discharge_location,
  adm.insurance,
  adm.race,
  fi.intime AS icu_intime,
  fi.outtime AS icu_outtime,
  EXTRACT(EPOCH FROM (fi.outtime - fi.intime)) / 3600.0 AS icu_los_hours,
  p.anchor_age AS age,
  p.gender,
  p.dod,
  s.suspected_infection_time,
  s.sofa_time,
  s.sofa_score,
  fi.intime + INTERVAL '72 hours' AS landmark_time,
  fi.intime - INTERVAL '24 hours' AS echo_window_start,
  fi.intime + INTERVAL '72 hours' AS echo_window_end,
  fi.intime - INTERVAL '24 hours' AS pni_window_start,
  fi.intime + INTERVAL '48 hours' AS pni_window_end
FROM tmp_first_icu fi
JOIN mimiciv_hosp.admissions adm
  ON adm.hadm_id = fi.hadm_id
JOIN mimiciv_hosp.patients p
  ON p.subject_id = fi.subject_id
JOIN mimiciv_derived.sepsis3 s
  ON s.stay_id = fi.stay_id
 AND s.sepsis3 = true
WHERE fi.rn = 1
  AND p.anchor_age >= 18
  AND EXTRACT(EPOCH FROM (fi.outtime - fi.intime)) / 3600.0 > 24;

CREATE INDEX tmp_base_pre_exclusion_stay_idx ON tmp_base_pre_exclusion (stay_id);
CREATE INDEX tmp_base_pre_exclusion_hadm_idx ON tmp_base_pre_exclusion (hadm_id);
CREATE INDEX tmp_base_pre_exclusion_subject_idx ON tmp_base_pre_exclusion (subject_id);
ANALYZE tmp_base_pre_exclusion;

DROP TABLE IF EXISTS pg_temp.tmp_exclusions;
CREATE TEMP TABLE tmp_exclusions AS
SELECT
  hadm_id,
  BOOL_OR(
    (icd_version = 10 AND (
      icd_code LIKE 'I50%'
      OR icd_code IN ('I110', 'I130', 'I132', 'I255')
    ))
    OR
    (icd_version = 9 AND (
      icd_code LIKE '428%'
      OR icd_code IN ('39891', '40201', '40211', '40291', '40401', '40403', '40411', '40413', '40491', '40493')
    ))
  ) AS has_heart_failure,
  BOOL_OR(
    (icd_version = 10 AND (
      icd_code LIKE 'I42%'
      OR icd_code LIKE 'I43%'
    ))
    OR
    (icd_version = 9 AND icd_code LIKE '425%')
  ) AS has_cardiomyopathy,
  BOOL_OR(
    (icd_version = 10 AND (
      icd_code LIKE 'I40%'
      OR icd_code LIKE 'I41%'
      OR icd_code = 'I514'
    ))
    OR
    (icd_version = 9 AND (
      icd_code LIKE '422%'
      OR icd_code = '4290'
    ))
  ) AS has_myocarditis
FROM mimiciv_hosp.diagnoses_icd
GROUP BY hadm_id;

CREATE INDEX tmp_exclusions_hadm_idx ON tmp_exclusions (hadm_id);
ANALYZE tmp_exclusions;

DROP TABLE IF EXISTS pg_temp.tmp_base;
CREATE TEMP TABLE tmp_base AS
SELECT
  b.*,
  COALESCE(x.has_heart_failure, false) AS excluded_heart_failure,
  COALESCE(x.has_cardiomyopathy, false) AS excluded_cardiomyopathy,
  COALESCE(x.has_myocarditis, false) AS excluded_myocarditis,
  CASE
    WHEN b.deathtime IS NOT NULL THEN b.deathtime
    WHEN b.dod IS NOT NULL THEN b.dod::timestamp
  END AS death_datetime_est,
  (
    (b.deathtime IS NOT NULL AND b.deathtime <= b.landmark_time)
    OR
    (b.deathtime IS NULL AND b.dod IS NOT NULL AND b.dod < b.landmark_time::date)
  ) AS died_before_landmark,
  (
    NOT (
      (b.deathtime IS NOT NULL AND b.deathtime <= b.landmark_time)
      OR
      (b.deathtime IS NULL AND b.dod IS NOT NULL AND b.dod < b.landmark_time::date)
    )
  ) AS alive_at_landmark,
  (b.dischtime > b.landmark_time) AS in_hospital_at_landmark
FROM tmp_base_pre_exclusion b
LEFT JOIN tmp_exclusions x
  ON x.hadm_id = b.hadm_id
WHERE COALESCE(x.has_heart_failure, false) = false
  AND COALESCE(x.has_cardiomyopathy, false) = false
  AND COALESCE(x.has_myocarditis, false) = false;

CREATE INDEX tmp_base_stay_idx ON tmp_base (stay_id);
CREATE INDEX tmp_base_hadm_idx ON tmp_base (hadm_id);
CREATE INDEX tmp_base_subject_idx ON tmp_base (subject_id);
ANALYZE tmp_base;

DROP TABLE IF EXISTS pg_temp.tmp_ef_raw;
CREATE TEMP TABLE tmp_ef_raw AS
SELECT
  b.subject_id,
  b.hadm_id,
  b.stay_id,
  e.measurement_id,
  e.measurement_datetime,
  e.lvef_lower,
  e.lvef_upper,
  e.biplane_lvef,
  e.lvef_3d,
  e.lvef_value,
  e.lvef_source
FROM tmp_base b
JOIN mimiciv_echo.lvef_measurements_20260531 e
  ON e.subject_id = b.subject_id
 AND e.test_type = 'tte'
 AND e.lvef_value IS NOT NULL
 AND e.measurement_datetime BETWEEN b.admittime AND b.dischtime
 AND e.measurement_datetime >= b.echo_window_start
 AND e.measurement_datetime <= b.echo_window_end;

CREATE INDEX tmp_ef_raw_stay_time_idx ON tmp_ef_raw (stay_id, measurement_datetime, measurement_id);
ANALYZE tmp_ef_raw;

DROP TABLE IF EXISTS pg_temp.tmp_ef_per_time;
CREATE TEMP TABLE tmp_ef_per_time AS
SELECT DISTINCT ON (stay_id, measurement_datetime)
  stay_id,
  measurement_datetime,
  measurement_id,
  lvef_lower,
  lvef_upper,
  biplane_lvef,
  lvef_3d,
  lvef_value,
  lvef_source
FROM tmp_ef_raw
ORDER BY
  stay_id,
  measurement_datetime,
  CASE lvef_source
    WHEN 'biplane_lvef' THEN 1
    WHEN 'lvef_3d' THEN 2
    WHEN 'lvef_midpoint' THEN 3
    WHEN 'lvef_lower' THEN 4
    WHEN 'lvef_upper' THEN 5
    ELSE 9
  END,
  measurement_id;

CREATE INDEX tmp_ef_per_time_stay_time_idx ON tmp_ef_per_time (stay_id, measurement_datetime);
ANALYZE tmp_ef_per_time;

DROP TABLE IF EXISTS pg_temp.tmp_ef_seq;
CREATE TEMP TABLE tmp_ef_seq AS
SELECT
  stay_id,
  measurement_datetime,
  measurement_id,
  lvef_value,
  lvef_source,
  MAX(lvef_value) OVER (
    PARTITION BY stay_id
    ORDER BY measurement_datetime
    ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING
  ) AS prior_max_lvef
FROM tmp_ef_per_time;

CREATE INDEX tmp_ef_seq_stay_idx ON tmp_ef_seq (stay_id);
ANALYZE tmp_ef_seq;

DROP TABLE IF EXISTS pg_temp.tmp_ef_summary;
CREATE TEMP TABLE tmp_ef_summary AS
WITH ef_flags AS (
  SELECT
    stay_id,
    COUNT(*) AS ef_timepoint_count,
    MIN(lvef_value) AS ef_min,
    MAX(lvef_value) AS ef_max,
    BOOL_OR(lvef_value < 50) AS ef_lt50
  FROM tmp_ef_per_time
  GROUP BY stay_id
),
drop_flags AS (
  SELECT
    stay_id,
    MAX(prior_max_lvef - lvef_value) AS ef_drop_max,
    BOOL_OR(prior_max_lvef - lvef_value > 10) AS ef_drop_gt10_raw,
    BOOL_OR(prior_max_lvef - lvef_value > 10 AND lvef_value < 55) AS ef_drop_gt10_to_lt55,
    BOOL_OR(prior_max_lvef - lvef_value > 10 AND lvef_value < 50) AS ef_drop_gt10_to_lt50
  FROM tmp_ef_seq
  GROUP BY stay_id
),
lowest_ef AS (
  SELECT DISTINCT ON (stay_id)
    stay_id,
    measurement_id AS ef_measurement_id,
    measurement_datetime AS ef_measurement_datetime,
    lvef_lower,
    lvef_upper,
    biplane_lvef,
    lvef_3d,
    lvef_value AS lvef_selected,
    lvef_source AS lvef_source_selected
  FROM tmp_ef_per_time
  ORDER BY stay_id, lvef_value, measurement_datetime, measurement_id
)
SELECT
  f.stay_id,
  f.ef_timepoint_count,
  f.ef_min,
  f.ef_max,
  f.ef_lt50,
  COALESCE(d.ef_drop_max, NULL) AS ef_drop_max,
  COALESCE(d.ef_drop_gt10_raw, false) AS ef_drop_gt10_raw,
  COALESCE(d.ef_drop_gt10_to_lt55, false) AS ef_drop_gt10_to_lt55,
  COALESCE(d.ef_drop_gt10_to_lt50, false) AS ef_drop_gt10_to_lt50,
  (f.ef_lt50 OR COALESCE(d.ef_drop_gt10_to_lt55, false)) AS scm_landmark72,
  (f.ef_lt50 OR COALESCE(d.ef_drop_gt10_to_lt50, false)) AS scm_strict_landmark72,
  le.ef_measurement_id,
  le.ef_measurement_datetime,
  le.lvef_lower,
  le.lvef_upper,
  le.biplane_lvef,
  le.lvef_3d,
  le.lvef_selected,
  le.lvef_source_selected
FROM ef_flags f
LEFT JOIN drop_flags d
  ON d.stay_id = f.stay_id
LEFT JOIN lowest_ef le
  ON le.stay_id = f.stay_id;

CREATE INDEX tmp_ef_summary_stay_idx ON tmp_ef_summary (stay_id);
ANALYZE tmp_ef_summary;

DROP TABLE IF EXISTS pg_temp.tmp_pni_base;
CREATE TEMP TABLE tmp_pni_base AS
SELECT
  b.subject_id,
  b.hadm_id,
  b.stay_id,
  b.pni_window_start,
  b.pni_window_end
FROM tmp_base b;

CREATE INDEX tmp_pni_base_hadm_idx ON tmp_pni_base (subject_id, hadm_id);
CREATE INDEX tmp_pni_base_stay_idx ON tmp_pni_base (stay_id);
ANALYZE tmp_pni_base;

DROP TABLE IF EXISTS pg_temp.tmp_hosp_pni_labs;
CREATE TEMP TABLE tmp_hosp_pni_labs AS
SELECT
  b.stay_id,
  le.specimen_id,
  le.charttime,
  le.itemid,
  di.label,
  le.valueuom,
  le.valuenum
FROM tmp_pni_base b
JOIN mimiciv_hosp.labevents le
  ON le.subject_id = b.subject_id
 AND le.hadm_id = b.hadm_id
 AND le.charttime >= b.pni_window_start
 AND le.charttime < b.pni_window_end
JOIN mimiciv_hosp.d_labitems di
  ON di.itemid = le.itemid
WHERE le.valuenum IS NOT NULL
  AND le.itemid IN (
    50862, 53085, 53138, 52022,
    51133, 52769, 53157,
    51244, 51245, 53188,
    51300, 51301, 51755, 51756
  );

CREATE INDEX tmp_hosp_pni_labs_idx ON tmp_hosp_pni_labs (stay_id, itemid, charttime);
ANALYZE tmp_hosp_pni_labs;

DROP TABLE IF EXISTS pg_temp.tmp_icu_pni_labs;
CREATE TEMP TABLE tmp_icu_pni_labs AS
SELECT
  b.stay_id,
  ce.charttime,
  ce.itemid,
  di.label,
  ce.valueuom,
  ce.valuenum
FROM tmp_pni_base b
JOIN mimiciv_icu.chartevents ce
  ON ce.stay_id = b.stay_id
 AND ce.charttime >= b.pni_window_start
 AND ce.charttime < b.pni_window_end
JOIN mimiciv_icu.d_items di
  ON di.itemid = ce.itemid
WHERE ce.valuenum IS NOT NULL
  AND ce.itemid IN (227456, 220574, 229358, 225641, 220546);

CREATE INDEX tmp_icu_pni_labs_idx ON tmp_icu_pni_labs (stay_id, itemid, charttime);
ANALYZE tmp_icu_pni_labs;

DROP TABLE IF EXISTS pg_temp.tmp_pni_unit_diagnostics;
CREATE TEMP TABLE tmp_pni_unit_diagnostics AS
WITH src AS (
  SELECT
    'hosp_labevents' AS source_table,
    itemid,
    label,
    valueuom,
    valuenum,
    stay_id
  FROM tmp_hosp_pni_labs
  UNION ALL
  SELECT
    'icu_chartevents' AS source_table,
    itemid,
    label,
    valueuom,
    valuenum,
    stay_id
  FROM tmp_icu_pni_labs
)
SELECT
  source_table,
  itemid,
  label,
  valueuom,
  COUNT(*) AS rows_n,
  COUNT(DISTINCT stay_id) AS stays_n,
  ROUND(MIN(valuenum)::numeric, 4) AS min_value,
  ROUND(PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY valuenum)::numeric, 4) AS p25_value,
  ROUND(PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY valuenum)::numeric, 4) AS median_value,
  ROUND(PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY valuenum)::numeric, 4) AS p75_value,
  ROUND(MAX(valuenum)::numeric, 4) AS max_value
FROM src
GROUP BY source_table, itemid, label, valueuom
ORDER BY source_table, label, itemid, valueuom;

DROP TABLE IF EXISTS pg_temp.tmp_albumin_candidates;
CREATE TEMP TABLE tmp_albumin_candidates AS
SELECT
  stay_id,
  charttime,
  CASE
    WHEN itemid IN (50862, 53085) THEN 'hosp_albumin_primary_50862_53085'
    WHEN itemid = 53138 THEN 'hosp_albumin_rescue_53138'
    WHEN itemid = 52022 THEN 'hosp_albumin_blood_gas_52022'
  END AS albumin_source,
  valuenum * 10.0 AS albumin_g_l,
  CASE
    WHEN itemid IN (50862, 53085) THEN 1
    WHEN itemid = 53138 THEN 2
    WHEN itemid = 52022 THEN 3
  END AS source_rank
FROM tmp_hosp_pni_labs
WHERE itemid IN (50862, 53085, 53138, 52022)
  AND valuenum * 10.0 BETWEEN 5 AND 80

UNION ALL

SELECT
  stay_id,
  charttime,
  CASE
    WHEN itemid = 227456 THEN 'icu_albumin_227456'
    WHEN itemid = 220574 THEN 'icu_zalbumin_220574'
  END AS albumin_source,
  CASE
    WHEN valuenum BETWEEN 0.5 AND 8 THEN valuenum * 10.0
    WHEN valuenum BETWEEN 5 AND 80 THEN valuenum
  END AS albumin_g_l,
  CASE
    WHEN itemid = 227456 THEN 4
    WHEN itemid = 220574 THEN 5
  END AS source_rank
FROM tmp_icu_pni_labs
WHERE itemid IN (227456, 220574)
  AND (
    valuenum BETWEEN 0.5 AND 8
    OR valuenum BETWEEN 5 AND 80
  );

DELETE FROM tmp_albumin_candidates
WHERE albumin_source IS NULL OR albumin_g_l IS NULL;

DROP TABLE IF EXISTS pg_temp.tmp_first_albumin;
CREATE TEMP TABLE tmp_first_albumin AS
SELECT DISTINCT ON (stay_id)
  stay_id,
  charttime AS albumin_charttime,
  ROUND(albumin_g_l::numeric, 2) AS albumin_g_l,
  albumin_source
FROM tmp_albumin_candidates
WHERE albumin_g_l BETWEEN 5 AND 80
ORDER BY stay_id, source_rank, charttime;

CREATE INDEX tmp_first_albumin_stay_idx ON tmp_first_albumin (stay_id);
ANALYZE tmp_first_albumin;

DROP TABLE IF EXISTS pg_temp.tmp_lymph_candidates;
CREATE TEMP TABLE tmp_lymph_candidates AS
WITH hosp_direct AS (
  SELECT DISTINCT ON (stay_id)
    stay_id,
    charttime,
    'hosp_direct_absolute' AS lymphocyte_source,
    CASE WHEN itemid = 52769 THEN valuenum / 1000.0 ELSE valuenum END AS lymphocyte_k_ul,
    1 AS source_rank
  FROM tmp_hosp_pni_labs
  WHERE itemid IN (51133, 52769, 53157)
    AND CASE WHEN itemid = 52769 THEN valuenum / 1000.0 ELSE valuenum END BETWEEN 0 AND 20
  ORDER BY stay_id, charttime, itemid
),
hosp_same_specimen AS (
  SELECT DISTINCT ON (stay_id)
    stay_id,
    charttime,
    'hosp_same_specimen_wbc_x_pct' AS lymphocyte_source,
    wbc_k_ul * lymphocyte_pct / 100.0 AS lymphocyte_k_ul,
    2 AS source_rank
  FROM (
    SELECT
      stay_id,
      specimen_id,
      MIN(charttime) AS charttime,
      MAX(valuenum) FILTER (
        WHERE itemid IN (51300, 51301, 51755, 51756)
          AND valuenum BETWEEN 0.1 AND 200
      ) AS wbc_k_ul,
      MAX(valuenum) FILTER (
        WHERE itemid IN (51244, 51245, 53188)
          AND valuenum BETWEEN 0 AND 100
      ) AS lymphocyte_pct
    FROM tmp_hosp_pni_labs
    WHERE specimen_id IS NOT NULL
    GROUP BY stay_id, specimen_id
  ) d
  WHERE wbc_k_ul IS NOT NULL
    AND lymphocyte_pct IS NOT NULL
    AND wbc_k_ul * lymphocyte_pct / 100.0 BETWEEN 0 AND 20
  ORDER BY stay_id, charttime
),
hosp_same_charttime AS (
  SELECT DISTINCT ON (stay_id)
    stay_id,
    charttime,
    'hosp_same_charttime_wbc_x_pct' AS lymphocyte_source,
    wbc_k_ul * lymphocyte_pct / 100.0 AS lymphocyte_k_ul,
    3 AS source_rank
  FROM (
    SELECT
      stay_id,
      charttime,
      MAX(valuenum) FILTER (
        WHERE itemid IN (51300, 51301, 51755, 51756)
          AND valuenum BETWEEN 0.1 AND 200
      ) AS wbc_k_ul,
      MAX(valuenum) FILTER (
        WHERE itemid IN (51244, 51245, 53188)
          AND valuenum BETWEEN 0 AND 100
      ) AS lymphocyte_pct
    FROM tmp_hosp_pni_labs
    GROUP BY stay_id, charttime
  ) d
  WHERE wbc_k_ul IS NOT NULL
    AND lymphocyte_pct IS NOT NULL
    AND wbc_k_ul * lymphocyte_pct / 100.0 BETWEEN 0 AND 20
  ORDER BY stay_id, charttime
),
icu_direct AS (
  SELECT DISTINCT ON (stay_id)
    stay_id,
    charttime,
    'icu_direct_absolute_229358' AS lymphocyte_source,
    CASE
      WHEN valuenum BETWEEN 0 AND 20 THEN valuenum
      WHEN valuenum BETWEEN 20 AND 20000 THEN valuenum / 1000.0
    END AS lymphocyte_k_ul,
    4 AS source_rank
  FROM tmp_icu_pni_labs
  WHERE itemid = 229358
    AND (
      valuenum BETWEEN 0 AND 20
      OR valuenum BETWEEN 20 AND 20000
    )
  ORDER BY stay_id, charttime
),
icu_same_charttime AS (
  SELECT DISTINCT ON (stay_id)
    stay_id,
    charttime,
    'icu_same_charttime_wbc_x_pct' AS lymphocyte_source,
    wbc_k_ul * lymphocyte_pct / 100.0 AS lymphocyte_k_ul,
    5 AS source_rank
  FROM (
    SELECT
      stay_id,
      charttime,
      MAX(valuenum) FILTER (WHERE itemid = 220546 AND valuenum BETWEEN 0.1 AND 200) AS wbc_k_ul,
      MAX(valuenum) FILTER (WHERE itemid = 225641 AND valuenum BETWEEN 0 AND 100) AS lymphocyte_pct
    FROM tmp_icu_pni_labs
    GROUP BY stay_id, charttime
  ) d
  WHERE wbc_k_ul IS NOT NULL
    AND lymphocyte_pct IS NOT NULL
    AND wbc_k_ul * lymphocyte_pct / 100.0 BETWEEN 0 AND 20
  ORDER BY stay_id, charttime
)
SELECT * FROM hosp_direct
UNION ALL SELECT * FROM hosp_same_specimen
UNION ALL SELECT * FROM hosp_same_charttime
UNION ALL SELECT * FROM icu_direct
UNION ALL SELECT * FROM icu_same_charttime;

DELETE FROM tmp_lymph_candidates
WHERE lymphocyte_source IS NULL OR lymphocyte_k_ul IS NULL;

DROP TABLE IF EXISTS pg_temp.tmp_first_lymphocyte;
CREATE TEMP TABLE tmp_first_lymphocyte AS
SELECT DISTINCT ON (stay_id)
  stay_id,
  charttime AS lymphocyte_charttime,
  ROUND(lymphocyte_k_ul::numeric, 3) AS lymphocyte_k_ul,
  lymphocyte_source
FROM tmp_lymph_candidates
WHERE lymphocyte_k_ul BETWEEN 0 AND 20
ORDER BY stay_id, source_rank, charttime;

CREATE INDEX tmp_first_lymphocyte_stay_idx ON tmp_first_lymphocyte (stay_id);
ANALYZE tmp_first_lymphocyte;

DROP TABLE IF EXISTS pg_temp.tmp_lactate_raw;
CREATE TEMP TABLE tmp_lactate_raw AS
WITH lab_lactate AS (
  SELECT
    b.stay_id,
    MAX(le.valuenum) FILTER (WHERE le.valuenum BETWEEN 0.1 AND 30) AS lactate_lab_max
  FROM tmp_base b
  JOIN mimiciv_hosp.labevents le
    ON le.hadm_id = b.hadm_id
   AND le.charttime >= b.icu_intime - INTERVAL '6 hours'
   AND le.charttime < b.icu_intime + INTERVAL '24 hours'
   AND le.itemid IN (50813, 52442, 53154)
   AND le.valuenum IS NOT NULL
  GROUP BY b.stay_id
),
bg_lactate AS (
  SELECT
    b.stay_id,
    MAX(bg.lactate) FILTER (WHERE bg.lactate BETWEEN 0.1 AND 30) AS lactate_bg_max
  FROM tmp_base b
  JOIN mimiciv_derived.bg bg
    ON bg.hadm_id = b.hadm_id
   AND bg.charttime >= b.icu_intime - INTERVAL '6 hours'
   AND bg.charttime < b.icu_intime + INTERVAL '24 hours'
  GROUP BY b.stay_id
)
SELECT
  b.stay_id,
  GREATEST(
    COALESCE(fd.lactate_max, ll.lactate_lab_max, bl.lactate_bg_max),
    COALESCE(ll.lactate_lab_max, fd.lactate_max, bl.lactate_bg_max),
    COALESCE(bl.lactate_bg_max, fd.lactate_max, ll.lactate_lab_max)
  ) AS lactate_max
FROM tmp_base b
LEFT JOIN mimiciv_derived.first_day_bg fd
  ON fd.stay_id = b.stay_id
LEFT JOIN lab_lactate ll
  ON ll.stay_id = b.stay_id
LEFT JOIN bg_lactate bl
  ON bl.stay_id = b.stay_id;

CREATE INDEX tmp_lactate_raw_stay_idx ON tmp_lactate_raw (stay_id);
ANALYZE tmp_lactate_raw;

DROP TABLE IF EXISTS pg_temp.tmp_bmi_rescue;
CREATE TEMP TABLE tmp_bmi_rescue AS
WITH base_subjects AS (
  SELECT DISTINCT subject_id
  FROM tmp_base
),
same_stay_height AS (
  SELECT
    b.stay_id,
    percentile_cont(0.5) WITHIN GROUP (ORDER BY h.height::numeric) AS height_cm
  FROM tmp_base b
  JOIN mimiciv_derived.height h
    ON h.stay_id = b.stay_id
   AND h.height BETWEEN 120 AND 220
  GROUP BY b.stay_id
),
same_hadm_height AS (
  SELECT
    b.stay_id,
    percentile_cont(0.5) WITHIN GROUP (
      ORDER BY CASE
        WHEN ce.itemid = 226707 THEN ce.valuenum * 2.54
        WHEN ce.itemid = 226730 THEN ce.valuenum
      END
    ) AS height_cm
  FROM tmp_base b
  JOIN mimiciv_icu.chartevents ce
    ON ce.hadm_id = b.hadm_id
   AND ce.valuenum IS NOT NULL
   AND (
     (ce.itemid = 226707 AND ce.valuenum BETWEEN 48 AND 87)
     OR (ce.itemid = 226730 AND ce.valuenum BETWEEN 120 AND 220)
   )
  GROUP BY b.stay_id
),
all_subject_height AS (
  SELECT
    h.subject_id,
    h.height::numeric AS height_cm
  FROM mimiciv_derived.height h
  JOIN base_subjects bs
    ON bs.subject_id = h.subject_id
  WHERE h.height BETWEEN 120 AND 220

  UNION ALL

  SELECT
    ce.subject_id,
    CASE
      WHEN ce.itemid = 226707 THEN ce.valuenum * 2.54
      WHEN ce.itemid = 226730 THEN ce.valuenum
    END AS height_cm
  FROM mimiciv_icu.chartevents ce
  JOIN base_subjects bs
    ON bs.subject_id = ce.subject_id
  WHERE ce.valuenum IS NOT NULL
    AND (
      (ce.itemid = 226707 AND ce.valuenum BETWEEN 48 AND 87)
      OR (ce.itemid = 226730 AND ce.valuenum BETWEEN 120 AND 220)
    )

  UNION ALL

  SELECT
    omr.subject_id,
    omr.result_value::numeric * 2.54 AS height_cm
  FROM mimiciv_hosp.omr omr
  JOIN base_subjects bs
    ON bs.subject_id = omr.subject_id
  WHERE omr.result_name IN ('Height', 'Height (Inches)')
    AND omr.result_value ~ '^ *[0-9]+([.][0-9]+)? *$'
    AND omr.result_value::numeric BETWEEN 48 AND 87
),
subject_height AS (
  SELECT
    subject_id,
    percentile_cont(0.5) WITHIN GROUP (ORDER BY height_cm) AS height_cm
  FROM all_subject_height
  WHERE height_cm BETWEEN 120 AND 220
  GROUP BY subject_id
),
weight AS (
  SELECT
    b.stay_id,
    COALESCE(
      CASE WHEN fdw.weight_admit BETWEEN 20 AND 300 THEN fdw.weight_admit END,
      CASE WHEN fdw.weight BETWEEN 20 AND 300 THEN fdw.weight END,
      CASE WHEN wd.weight BETWEEN 20 AND 300 THEN wd.weight END
    ) AS weight_kg
  FROM tmp_base b
  LEFT JOIN mimiciv_derived.first_day_weight fdw
    ON fdw.stay_id = b.stay_id
  LEFT JOIN LATERAL (
    SELECT w.weight
    FROM mimiciv_derived.weight_durations w
    WHERE w.stay_id = b.stay_id
      AND w.weight BETWEEN 20 AND 300
    ORDER BY
      CASE WHEN w.weight_type = 'admit' THEN 0 ELSE 1 END,
      ABS(EXTRACT(EPOCH FROM (w.starttime - b.icu_intime)))
    LIMIT 1
  ) wd
    ON true
)
SELECT
  b.stay_id,
  CASE
    WHEN w.weight_kg BETWEEN 20 AND 300
     AND COALESCE(ssh.height_cm, shh.height_cm, subj.height_cm) BETWEEN 120 AND 220
    THEN w.weight_kg / POWER(COALESCE(ssh.height_cm, shh.height_cm, subj.height_cm) / 100.0, 2)
  END AS bmi
FROM tmp_base b
LEFT JOIN same_stay_height ssh
  ON ssh.stay_id = b.stay_id
LEFT JOIN same_hadm_height shh
  ON shh.stay_id = b.stay_id
LEFT JOIN subject_height subj
  ON subj.subject_id = b.subject_id
LEFT JOIN weight w
  ON w.stay_id = b.stay_id;

CREATE INDEX tmp_bmi_rescue_stay_idx ON tmp_bmi_rescue (stay_id);
ANALYZE tmp_bmi_rescue;

DROP TABLE IF EXISTS pg_temp.tmp_comorbidity_flags;
CREATE TEMP TABLE tmp_comorbidity_flags AS
SELECT
  d.hadm_id,
  BOOL_OR(
    ch.myocardial_infarct = 1
    OR (d.icd_version = 10 AND d.icd_code SIMILAR TO 'I2[0-5]%')
    OR (d.icd_version = 9 AND d.icd_code SIMILAR TO '41[0-4]%')
  ) AS coronary_artery_disease,
  BOOL_OR(ch.diabetes_without_cc = 1 OR ch.diabetes_with_cc = 1) AS diabetes_mellitus,
  BOOL_OR(
    (d.icd_version = 10 AND d.icd_code SIMILAR TO 'I1[0-6]%')
    OR (d.icd_version = 9 AND d.icd_code SIMILAR TO '40[1-5]%')
  ) AS hypertension,
  BOOL_OR(ch.renal_disease = 1) AS chronic_kidney_disease,
  BOOL_OR(ch.chronic_pulmonary_disease = 1) AS copd
FROM mimiciv_hosp.diagnoses_icd d
LEFT JOIN mimiciv_derived.charlson ch
  ON ch.hadm_id = d.hadm_id
GROUP BY d.hadm_id;

CREATE INDEX tmp_comorbidity_flags_hadm_idx ON tmp_comorbidity_flags (hadm_id);
ANALYZE tmp_comorbidity_flags;

DROP TABLE IF EXISTS pg_temp.tmp_support;
CREATE TEMP TABLE tmp_support AS
WITH vent AS (
  SELECT DISTINCT b.stay_id
  FROM tmp_base b
  JOIN mimiciv_derived.ventilation v
    ON v.stay_id = b.stay_id
   AND v.ventilation_status = 'InvasiveVent'
   AND v.starttime < b.icu_intime + INTERVAL '24 hours'
   AND COALESCE(v.endtime, v.starttime) >= b.icu_intime
),
vaso AS (
  SELECT DISTINCT b.stay_id
  FROM tmp_base b
  JOIN mimiciv_derived.vasoactive_agent va
    ON va.stay_id = b.stay_id
   AND va.starttime < b.icu_intime + INTERVAL '24 hours'
   AND COALESCE(va.endtime, va.starttime) >= b.icu_intime
)
SELECT
  b.stay_id,
  CASE WHEN vent.stay_id IS NOT NULL THEN 1 ELSE 0 END AS mechanical_vent,
  CASE WHEN vaso.stay_id IS NOT NULL THEN 1 ELSE 0 END AS vasopressor
FROM tmp_base b
LEFT JOIN vent
  ON vent.stay_id = b.stay_id
LEFT JOIN vaso
  ON vaso.stay_id = b.stay_id;

CREATE INDEX tmp_support_stay_idx ON tmp_support (stay_id);
ANALYZE tmp_support;

CREATE TABLE public.pni_scm_landmark72_20260605 AS
SELECT
  b.subject_id,
  b.hadm_id,
  b.stay_id,
  b.age,
  b.gender,
  CASE WHEN b.gender = 'F' THEN 1 ELSE 0 END AS female,
  ROUND(bmi.bmi::numeric, 2) AS bmi,
  b.admittime,
  b.dischtime,
  b.deathtime,
  b.dod,
  b.icu_intime,
  b.icu_outtime,
  b.icu_los_hours,
  b.suspected_infection_time,
  b.sofa_time,
  ROUND(b.sofa_score::numeric, 2) AS sofa,
  aps.apsiii,
  b.landmark_time,
  b.echo_window_start,
  b.echo_window_end,
  b.pni_window_start,
  b.pni_window_end,
  b.hospital_expire_flag AS hospital_death_raw,
  b.died_before_landmark,
  CASE WHEN b.alive_at_landmark THEN 1 ELSE 0 END AS alive_at_landmark,
  CASE WHEN b.in_hospital_at_landmark THEN 1 ELSE 0 END AS in_hospital_at_landmark,
  CASE WHEN e.stay_id IS NOT NULL THEN 1 ELSE 0 END AS echo_pre_landmark_available,
  e.ef_timepoint_count,
  ROUND(e.ef_min::numeric, 2) AS ef_min,
  ROUND(e.ef_max::numeric, 2) AS ef_max,
  ROUND(e.ef_drop_max::numeric, 2) AS ef_drop_max,
  e.ef_lt50,
  e.ef_drop_gt10_raw,
  e.ef_drop_gt10_to_lt55,
  e.ef_drop_gt10_to_lt50,
  CASE WHEN e.stay_id IS NULL THEN NULL ELSE CASE WHEN e.scm_landmark72 THEN 1 ELSE 0 END END AS scm,
  CASE WHEN e.stay_id IS NULL THEN NULL ELSE CASE WHEN e.scm_strict_landmark72 THEN 1 ELSE 0 END END AS scm_strict_lvef_lt50,
  e.ef_measurement_id,
  e.ef_measurement_datetime,
  ROUND(e.lvef_selected::numeric, 2) AS lvef_selected,
  e.lvef_source_selected,
  a.albumin_g_l,
  a.albumin_charttime,
  a.albumin_source,
  l.lymphocyte_k_ul,
  l.lymphocyte_charttime,
  l.lymphocyte_source,
  CASE
    WHEN a.albumin_g_l IS NOT NULL AND l.lymphocyte_k_ul IS NOT NULL
    THEN ROUND((a.albumin_g_l + 5.0 * l.lymphocyte_k_ul)::numeric, 2)
  END AS pni,
  CASE WHEN a.albumin_g_l IS NOT NULL AND l.lymphocyte_k_ul IS NOT NULL THEN 1 ELSE 0 END AS pni_complete,
  ROUND(lr.lactate_max::numeric, 2) AS lactate,
  ROUND(fd_lab.creatinine_max::numeric, 2) AS creatinine,
  ROUND(fd_lab.wbc_max::numeric, 2) AS wbc,
  ROUND(fd_lab.platelets_min::numeric, 2) AS platelets,
  ROUND((fd_lab.bun_max * 0.357)::numeric, 2) AS bun_mmol_l,
  sp.mechanical_vent,
  sp.vasopressor,
  b.admission_type,
  b.admission_location,
  b.discharge_location,
  b.insurance,
  b.race,
  CASE WHEN COALESCE(cf.coronary_artery_disease, false) THEN 1 ELSE 0 END AS coronary_artery_disease,
  CASE WHEN COALESCE(cf.diabetes_mellitus, false) THEN 1 ELSE 0 END AS diabetes_mellitus,
  CASE WHEN COALESCE(cf.hypertension, false) THEN 1 ELSE 0 END AS hypertension,
  CASE WHEN COALESCE(cf.chronic_kidney_disease, false) THEN 1 ELSE 0 END AS chronic_kidney_disease,
  CASE WHEN COALESCE(cf.copd, false) THEN 1 ELSE 0 END AS copd,
  b.excluded_heart_failure,
  b.excluded_cardiomyopathy,
  b.excluded_myocarditis,
  CASE
    WHEN b.deathtime IS NOT NULL AND b.deathtime > b.landmark_time
    THEN EXTRACT(EPOCH FROM (b.deathtime - b.landmark_time)) / 86400.0
    WHEN b.deathtime IS NULL AND b.dod IS NOT NULL AND b.dod >= b.landmark_time::date
    THEN (b.dod - b.landmark_time::date)::numeric + 0.5
  END AS death_days_from_landmark,
  CASE
    WHEN b.deathtime IS NOT NULL AND b.deathtime > b.landmark_time THEN 1
    WHEN b.deathtime IS NULL AND b.dod IS NOT NULL AND b.dod >= b.landmark_time::date THEN 1
    ELSE 0
  END AS death_after_landmark,
  CASE
    WHEN b.in_hospital_at_landmark
     AND b.hospital_expire_flag = 1
     AND COALESCE(b.deathtime, b.dischtime) > b.landmark_time
    THEN 1 ELSE 0
  END AS hosp_event,
  CASE
    WHEN b.in_hospital_at_landmark
     AND b.hospital_expire_flag = 1
     AND COALESCE(b.deathtime, b.dischtime) > b.landmark_time
    THEN GREATEST(EXTRACT(EPOCH FROM (COALESCE(b.deathtime, b.dischtime) - b.landmark_time)) / 86400.0, 0.01)
    WHEN b.in_hospital_at_landmark
    THEN GREATEST(EXTRACT(EPOCH FROM (b.dischtime - b.landmark_time)) / 86400.0, 0.01)
  END AS hosp_time_days,
  CASE
    WHEN (
      CASE
        WHEN b.deathtime IS NOT NULL AND b.deathtime > b.landmark_time
        THEN EXTRACT(EPOCH FROM (b.deathtime - b.landmark_time)) / 86400.0
        WHEN b.deathtime IS NULL AND b.dod IS NOT NULL AND b.dod >= b.landmark_time::date
        THEN (b.dod - b.landmark_time::date)::numeric + 0.5
      END
    ) <= 28 THEN 1 ELSE 0
  END AS event_28d,
  CASE
    WHEN (
      CASE
        WHEN b.deathtime IS NOT NULL AND b.deathtime > b.landmark_time
        THEN EXTRACT(EPOCH FROM (b.deathtime - b.landmark_time)) / 86400.0
        WHEN b.deathtime IS NULL AND b.dod IS NOT NULL AND b.dod >= b.landmark_time::date
        THEN (b.dod - b.landmark_time::date)::numeric + 0.5
      END
    ) <= 28
    THEN GREATEST(
      CASE
        WHEN b.deathtime IS NOT NULL AND b.deathtime > b.landmark_time
        THEN EXTRACT(EPOCH FROM (b.deathtime - b.landmark_time)) / 86400.0
        WHEN b.deathtime IS NULL AND b.dod IS NOT NULL AND b.dod >= b.landmark_time::date
        THEN (b.dod - b.landmark_time::date)::numeric + 0.5
      END,
      0.01
    )
    ELSE 28.0
  END AS time_28d,
  CASE
    WHEN (
      CASE
        WHEN b.deathtime IS NOT NULL AND b.deathtime > b.landmark_time
        THEN EXTRACT(EPOCH FROM (b.deathtime - b.landmark_time)) / 86400.0
        WHEN b.deathtime IS NULL AND b.dod IS NOT NULL AND b.dod >= b.landmark_time::date
        THEN (b.dod - b.landmark_time::date)::numeric + 0.5
      END
    ) <= 90 THEN 1 ELSE 0
  END AS event_90d,
  CASE
    WHEN (
      CASE
        WHEN b.deathtime IS NOT NULL AND b.deathtime > b.landmark_time
        THEN EXTRACT(EPOCH FROM (b.deathtime - b.landmark_time)) / 86400.0
        WHEN b.deathtime IS NULL AND b.dod IS NOT NULL AND b.dod >= b.landmark_time::date
        THEN (b.dod - b.landmark_time::date)::numeric + 0.5
      END
    ) <= 90
    THEN GREATEST(
      CASE
        WHEN b.deathtime IS NOT NULL AND b.deathtime > b.landmark_time
        THEN EXTRACT(EPOCH FROM (b.deathtime - b.landmark_time)) / 86400.0
        WHEN b.deathtime IS NULL AND b.dod IS NOT NULL AND b.dod >= b.landmark_time::date
        THEN (b.dod - b.landmark_time::date)::numeric + 0.5
      END,
      0.01
    )
    ELSE 90.0
  END AS time_90d,
  CASE
    WHEN a.albumin_g_l IS NOT NULL AND l.lymphocyte_k_ul IS NOT NULL
    THEN ROUND((-(a.albumin_g_l + 5.0 * l.lymphocyte_k_ul) / 5.0)::numeric, 6)
  END AS pni_low5,
  CASE
    WHEN a.albumin_g_l IS NOT NULL THEN ROUND((-a.albumin_g_l / 5.0)::numeric, 6)
  END AS albumin_low5,
  CASE
    WHEN l.lymphocyte_k_ul IS NOT NULL THEN ROUND((-l.lymphocyte_k_ul / 0.5)::numeric, 6)
  END AS lymph_low05,
  CASE
    WHEN lr.lactate_max > 0 THEN LN(lr.lactate_max)
  END AS log_lactate,
  CASE
    WHEN fd_lab.creatinine_max > 0 THEN LN(fd_lab.creatinine_max)
  END AS log_creatinine
FROM tmp_base b
LEFT JOIN tmp_ef_summary e
  ON e.stay_id = b.stay_id
LEFT JOIN tmp_bmi_rescue bmi
  ON bmi.stay_id = b.stay_id
LEFT JOIN tmp_first_albumin a
  ON a.stay_id = b.stay_id
LEFT JOIN tmp_first_lymphocyte l
  ON l.stay_id = b.stay_id
LEFT JOIN mimiciv_derived.apsiii aps
  ON aps.stay_id = b.stay_id
LEFT JOIN mimiciv_derived.first_day_lab fd_lab
  ON fd_lab.stay_id = b.stay_id
LEFT JOIN tmp_lactate_raw lr
  ON lr.stay_id = b.stay_id
LEFT JOIN tmp_support sp
  ON sp.stay_id = b.stay_id
LEFT JOIN tmp_comorbidity_flags cf
  ON cf.hadm_id = b.hadm_id;

CREATE INDEX pni_scm_landmark72_20260605_stay_idx ON public.pni_scm_landmark72_20260605 (stay_id);
CREATE INDEX pni_scm_landmark72_20260605_subject_idx ON public.pni_scm_landmark72_20260605 (subject_id);
ANALYZE public.pni_scm_landmark72_20260605;

DROP TABLE IF EXISTS pg_temp.tmp_landmark_flow_counts;
CREATE TEMP TABLE tmp_landmark_flow_counts AS
SELECT '01_first_icu_adult_sepsis3_los_gt24' AS step, COUNT(*) AS n
FROM tmp_base_pre_exclusion
UNION ALL
SELECT '02_after_hf_cardiomyopathy_myocarditis_exclusions', COUNT(*)
FROM tmp_base
UNION ALL
SELECT '03_alive_at_72h_landmark', COUNT(*)
FROM public.pni_scm_landmark72_20260605
WHERE alive_at_landmark = 1
UNION ALL
SELECT '04_in_hospital_at_72h_landmark', COUNT(*)
FROM public.pni_scm_landmark72_20260605
WHERE alive_at_landmark = 1 AND in_hospital_at_landmark = 1
UNION ALL
SELECT '05_alive_with_pre_landmark_tte_lvef', COUNT(*)
FROM public.pni_scm_landmark72_20260605
WHERE alive_at_landmark = 1 AND echo_pre_landmark_available = 1
UNION ALL
SELECT '06_alive_with_tte_and_pni_complete_28d_90d_cohort', COUNT(*)
FROM public.pni_scm_landmark72_20260605
WHERE alive_at_landmark = 1 AND echo_pre_landmark_available = 1 AND pni_complete = 1
UNION ALL
SELECT '07_alive_in_hospital_with_tte_and_pni_complete_hospital_cohort', COUNT(*)
FROM public.pni_scm_landmark72_20260605
WHERE alive_at_landmark = 1 AND in_hospital_at_landmark = 1 AND echo_pre_landmark_available = 1 AND pni_complete = 1;

DROP TABLE IF EXISTS pg_temp.tmp_landmark_selection_summary;
CREATE TEMP TABLE tmp_landmark_selection_summary AS
SELECT
  COUNT(*) AS cohort_after_exclusion_n,
  COUNT(*) FILTER (WHERE died_before_landmark) AS died_before_landmark_n,
  COUNT(*) FILTER (WHERE alive_at_landmark = 1) AS alive_at_landmark_n,
  COUNT(*) FILTER (WHERE alive_at_landmark = 1 AND in_hospital_at_landmark = 1) AS in_hospital_at_landmark_n,
  COUNT(*) FILTER (WHERE alive_at_landmark = 1 AND echo_pre_landmark_available = 1) AS pre_landmark_tte_n,
  COUNT(*) FILTER (WHERE alive_at_landmark = 1 AND echo_pre_landmark_available = 1 AND scm = 1) AS scm_n,
  COUNT(*) FILTER (WHERE alive_at_landmark = 1 AND echo_pre_landmark_available = 1 AND pni_complete = 1) AS pni_complete_n,
  COUNT(*) FILTER (WHERE alive_at_landmark = 1 AND echo_pre_landmark_available = 1 AND pni_complete = 1 AND scm = 1) AS scm_pni_complete_n,
  COUNT(*) FILTER (WHERE alive_at_landmark = 1 AND echo_pre_landmark_available = 1 AND pni_complete = 1 AND event_28d = 1) AS death_28d_n,
  COUNT(*) FILTER (WHERE alive_at_landmark = 1 AND echo_pre_landmark_available = 1 AND pni_complete = 1 AND event_90d = 1) AS death_90d_n,
  COUNT(*) FILTER (WHERE alive_at_landmark = 1 AND in_hospital_at_landmark = 1 AND echo_pre_landmark_available = 1 AND pni_complete = 1 AND hosp_event = 1) AS hospital_death_n
FROM public.pni_scm_landmark72_20260605;

DROP TABLE IF EXISTS pg_temp.tmp_landmark_missingness;
CREATE TEMP TABLE tmp_landmark_missingness AS
WITH d AS (
  SELECT *
  FROM public.pni_scm_landmark72_20260605
  WHERE alive_at_landmark = 1
    AND echo_pre_landmark_available = 1
)
SELECT 'albumin_g_l' AS variable, COUNT(*) AS n, COUNT(albumin_g_l) AS nonmissing_n, COUNT(*) - COUNT(albumin_g_l) AS missing_n, ROUND(100.0 * (COUNT(*) - COUNT(albumin_g_l)) / COUNT(*), 2) AS missing_pct FROM d
UNION ALL SELECT 'lymphocyte_k_ul', COUNT(*), COUNT(lymphocyte_k_ul), COUNT(*) - COUNT(lymphocyte_k_ul), ROUND(100.0 * (COUNT(*) - COUNT(lymphocyte_k_ul)) / COUNT(*), 2) FROM d
UNION ALL SELECT 'pni', COUNT(*), COUNT(pni), COUNT(*) - COUNT(pni), ROUND(100.0 * (COUNT(*) - COUNT(pni)) / COUNT(*), 2) FROM d
UNION ALL SELECT 'sofa', COUNT(*), COUNT(sofa), COUNT(*) - COUNT(sofa), ROUND(100.0 * (COUNT(*) - COUNT(sofa)) / COUNT(*), 2) FROM d
UNION ALL SELECT 'apsiii', COUNT(*), COUNT(apsiii), COUNT(*) - COUNT(apsiii), ROUND(100.0 * (COUNT(*) - COUNT(apsiii)) / COUNT(*), 2) FROM d
UNION ALL SELECT 'lactate', COUNT(*), COUNT(lactate), COUNT(*) - COUNT(lactate), ROUND(100.0 * (COUNT(*) - COUNT(lactate)) / COUNT(*), 2) FROM d
UNION ALL SELECT 'creatinine', COUNT(*), COUNT(creatinine), COUNT(*) - COUNT(creatinine), ROUND(100.0 * (COUNT(*) - COUNT(creatinine)) / COUNT(*), 2) FROM d
UNION ALL SELECT 'bmi', COUNT(*), COUNT(bmi), COUNT(*) - COUNT(bmi), ROUND(100.0 * (COUNT(*) - COUNT(bmi)) / COUNT(*), 2) FROM d
UNION ALL SELECT 'wbc', COUNT(*), COUNT(wbc), COUNT(*) - COUNT(wbc), ROUND(100.0 * (COUNT(*) - COUNT(wbc)) / COUNT(*), 2) FROM d
UNION ALL SELECT 'platelets', COUNT(*), COUNT(platelets), COUNT(*) - COUNT(platelets), ROUND(100.0 * (COUNT(*) - COUNT(platelets)) / COUNT(*), 2) FROM d
UNION ALL SELECT 'bun_mmol_l', COUNT(*), COUNT(bun_mmol_l), COUNT(*) - COUNT(bun_mmol_l), ROUND(100.0 * (COUNT(*) - COUNT(bun_mmol_l)) / COUNT(*), 2) FROM d
UNION ALL SELECT 'mechanical_vent', COUNT(*), COUNT(mechanical_vent), COUNT(*) - COUNT(mechanical_vent), ROUND(100.0 * (COUNT(*) - COUNT(mechanical_vent)) / COUNT(*), 2) FROM d
UNION ALL SELECT 'vasopressor', COUNT(*), COUNT(vasopressor), COUNT(*) - COUNT(vasopressor), ROUND(100.0 * (COUNT(*) - COUNT(vasopressor)) / COUNT(*), 2) FROM d;

DROP TABLE IF EXISTS pg_temp.tmp_landmark_echo_selection;
CREATE TEMP TABLE tmp_landmark_echo_selection AS
WITH d AS (
  SELECT
    CASE WHEN echo_pre_landmark_available = 1 THEN 'pre_landmark_tte_available' ELSE 'no_pre_landmark_tte' END AS group_name,
    *
  FROM public.pni_scm_landmark72_20260605
  WHERE alive_at_landmark = 1
)
SELECT
  group_name,
  COUNT(*) AS n,
  ROUND(AVG(age)::numeric, 2) AS age_mean,
  ROUND(AVG(CASE WHEN gender = 'F' THEN 1 ELSE 0 END)::numeric, 3) AS female_pct_fraction,
  ROUND(AVG(sofa)::numeric, 2) AS sofa_mean,
  ROUND(AVG(CASE WHEN event_28d = 1 THEN 1 ELSE 0 END)::numeric, 3) AS death_28d_fraction,
  ROUND(AVG(CASE WHEN event_90d = 1 THEN 1 ELSE 0 END)::numeric, 3) AS death_90d_fraction,
  ROUND(AVG(CASE WHEN hospital_death_raw = 1 THEN 1 ELSE 0 END)::numeric, 3) AS hospital_death_raw_fraction,
  COUNT(*) FILTER (WHERE pni_complete = 1) AS pni_complete_n
FROM d
GROUP BY group_name
ORDER BY group_name;

DROP TABLE IF EXISTS pg_temp.tmp_landmark_pni_selection;
CREATE TEMP TABLE tmp_landmark_pni_selection AS
WITH d AS (
  SELECT
    CASE WHEN pni_complete = 1 THEN 'pni_complete' ELSE 'pni_missing' END AS group_name,
    *
  FROM public.pni_scm_landmark72_20260605
  WHERE alive_at_landmark = 1
    AND echo_pre_landmark_available = 1
)
SELECT
  group_name,
  COUNT(*) AS n,
  ROUND(AVG(age)::numeric, 2) AS age_mean,
  ROUND(AVG(CASE WHEN gender = 'F' THEN 1 ELSE 0 END)::numeric, 3) AS female_pct_fraction,
  ROUND(AVG(sofa)::numeric, 2) AS sofa_mean,
  ROUND(AVG(CASE WHEN scm = 1 THEN 1 ELSE 0 END)::numeric, 3) AS scm_fraction,
  ROUND(AVG(CASE WHEN event_28d = 1 THEN 1 ELSE 0 END)::numeric, 3) AS death_28d_fraction,
  ROUND(AVG(CASE WHEN event_90d = 1 THEN 1 ELSE 0 END)::numeric, 3) AS death_90d_fraction
FROM d
GROUP BY group_name
ORDER BY group_name;

SELECT * FROM tmp_landmark_flow_counts ORDER BY step;
SELECT * FROM tmp_landmark_selection_summary;
SELECT * FROM tmp_landmark_missingness ORDER BY missing_pct DESC, variable;
SELECT * FROM tmp_landmark_echo_selection;
SELECT * FROM tmp_landmark_pni_selection;

\copy (SELECT * FROM public.pni_scm_landmark72_20260605 ORDER BY subject_id, hadm_id, stay_id) TO 'outputs/pni_scm_landmark72_20260605/landmark72_all_after_exclusions_20260605.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8')
\copy (SELECT * FROM public.pni_scm_landmark72_20260605 WHERE alive_at_landmark = 1 AND echo_pre_landmark_available = 1 AND pni_complete = 1 ORDER BY subject_id, hadm_id, stay_id) TO 'outputs/pni_scm_landmark72_20260605/landmark72_analysis_dataset_20260605.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8')
\copy (SELECT * FROM tmp_landmark_flow_counts ORDER BY step) TO 'outputs/pni_scm_landmark72_20260605/landmark72_flow_counts_20260605.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8')
\copy (SELECT * FROM tmp_landmark_selection_summary) TO 'outputs/pni_scm_landmark72_20260605/landmark72_selection_summary_20260605.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8')
\copy (SELECT * FROM tmp_landmark_missingness ORDER BY missing_pct DESC, variable) TO 'outputs/pni_scm_landmark72_20260605/landmark72_missingness_20260605.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8')
\copy (SELECT * FROM tmp_landmark_echo_selection) TO 'outputs/pni_scm_landmark72_20260605/landmark72_echo_selection_20260605.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8')
\copy (SELECT * FROM tmp_landmark_pni_selection) TO 'outputs/pni_scm_landmark72_20260605/landmark72_pni_selection_20260605.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8')
\copy (SELECT * FROM tmp_pni_unit_diagnostics) TO 'outputs/pni_scm_landmark72_20260605/landmark72_pni_unit_diagnostics_20260605.csv' WITH (FORMAT csv, HEADER true, ENCODING 'UTF8')
