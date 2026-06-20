library(readxl)

source("00_通用设置与函数_seed2025.R", encoding = "UTF-8")

m <- read.csv(mimic_path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
d <- read_excel(domestic_path)

mimic <- data.frame(
  cohort = "MIMIC-IV",
  scm = to_num(m$scm),
  age = to_num(m$age),
  female = to_num(m$female),
  bmi = to_num(m$bmi),
  sofa = to_num(m$sofa),
  lvef = to_num(m$lvef_selected),
  albumin = to_num(m$albumin_g_l),
  lymphocyte = to_num(m$lymphocyte_k_ul),
  pni = to_num(m$pni),
  lactate = to_num(m$lactate),
  creatinine = to_num(m$creatinine) * 88.4,
  wbc = to_num(m$wbc),
  platelet = to_num(m$platelets),
  bun = to_num(m$bun_mmol_l),
  cad = to_num(m$coronary_artery_disease),
  diabetes = to_num(m$diabetes_mellitus),
  hypertension = to_num(m$hypertension),
  ckd = to_num(m$chronic_kidney_disease),
  copd = to_num(m$copd),
  ventilation = to_num(m$mechanical_vent),
  vasopressor = to_num(m$vasopressor),
  event28 = to_num(m$event_28d),
  event90 = to_num(m$event_90d),
  hospital_or_terminal = NA_real_
)

sex_d <- trimws(as.character(d$性别))
domestic <- data.frame(
  cohort = "Domestic",
  scm = ifelse(to_num(d$射血分数) < 50, 1, 0),
  age = to_num(d$年龄),
  female = ifelse(sex_d == "女", 1, ifelse(sex_d == "男", 0, NA)),
  bmi = to_num(d$bmi),
  sofa = to_num(d$sofa),
  lvef = to_num(d$射血分数),
  albumin = to_num(d$白蛋白),
  lymphocyte = to_num(d$淋巴细胞计数),
  pni = to_num(d$白蛋白) + 5 * to_num(d$淋巴细胞计数),
  lactate = to_num(d$lac),
  creatinine = to_num(d$cr),
  wbc = to_num(d$wbc),
  platelet = to_num(d$plt),
  bun = to_num(d$bun),
  cad = to_binary(d$冠心病),
  diabetes = to_binary(d$糖尿病),
  hypertension = to_binary(d$高血压),
  ckd = to_binary(d$慢性肾病),
  copd = to_binary(d$慢阻肺),
  ventilation = to_binary(d$呼吸机),
  vasopressor = to_binary(d$升压药),
  event28 = NA_real_,
  event90 = NA_real_,
  hospital_or_terminal = to_binary(d$结局)
)

variables <- data.frame(
  Section = c(
    "Demographics and severity", "Demographics and severity", "Demographics and severity", "Demographics and severity",
    "Cardiac and nutritional markers", "Cardiac and nutritional markers", "Cardiac and nutritional markers", "Cardiac and nutritional markers",
    "Laboratory markers", "Laboratory markers", "Laboratory markers", "Laboratory markers", "Laboratory markers",
    "Comorbidities", "Comorbidities", "Comorbidities", "Comorbidities", "Comorbidities",
    "Organ support", "Organ support",
    "Outcomes; descriptive only", "Outcomes; descriptive only", "Outcomes; descriptive only"
  ),
  Variable = c(
    "Age, years", "Female sex, n (%)", "BMI, kg/m2", "SOFA score",
    "LVEF, %", "Albumin, g/L", "Lymphocyte count, x10^9/L", "PNI",
    "Lactate, mmol/L", "Creatinine, umol/L", "WBC, x10^9/L", "Platelet count, x10^9/L", "BUN, mmol/L",
    "Coronary artery disease, n (%)", "Diabetes mellitus, n (%)", "Hypertension, n (%)", "Chronic kidney disease, n (%)", "COPD, n (%)",
    "Mechanical ventilation, n (%)", "Vasopressor use, n (%)",
    "28-day mortality, n (%)", "90-day mortality, n (%)", "In-hospital death or terminal discharge, n (%)"
  ),
  Name = c(
    "age", "female", "bmi", "sofa",
    "lvef", "albumin", "lymphocyte", "pni",
    "lactate", "creatinine", "wbc", "platelet", "bun",
    "cad", "diabetes", "hypertension", "ckd", "copd",
    "ventilation", "vasopressor",
    "event28", "event90", "hospital_or_terminal"
  ),
  Type = c(
    "cont", "cat", "cont", "cont",
    "cont", "cont", "cont", "cont",
    "cont", "cont", "cont", "cont", "cont",
    "cat", "cat", "cat", "cat", "cat",
    "cat", "cat",
    "cat", "cat", "cat"
  ),
  Digits = c(1, NA, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  stringsAsFactors = FALSE
)

rows <- list()
for (i in seq_len(nrow(variables))) {
  v <- variables[i, ]
  blank_p <- grepl("Outcomes", v$Section)
  sm <- summarize_by_scm(mimic, v$Name, v$Type, digits = ifelse(is.na(v$Digits), 1, v$Digits), blank_p = blank_p)
  sd <- summarize_by_scm(domestic, v$Name, v$Type, digits = ifelse(is.na(v$Digits), 1, v$Digits), blank_p = blank_p)
  rows[[length(rows) + 1]] <- data.frame(
    Section = v$Section,
    Variable = v$Variable,
    MIMIC_Overall = sm$Overall,
    MIMIC_No_SCM = sm$No_SCM,
    MIMIC_SCM = sm$SCM,
    MIMIC_P_value = sm$P_value,
    Domestic_Overall = sd$Overall,
    Domestic_No_SCM = sd$No_SCM,
    Domestic_SCM = sd$SCM,
    Domestic_P_value = sd$P_value,
    stringsAsFactors = FALSE
  )
}

table1 <- do.call(rbind, rows)
cat("Table 1 check: MIMIC n=", nrow(mimic), ", SCM=", sum(mimic$scm == 1), "; Domestic n=", nrow(domestic), ", SCM=", sum(domestic$scm == 1), "\n", sep = "")
print_table(table1)


