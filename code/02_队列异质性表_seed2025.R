library(readxl)

source("00_通用设置与函数_seed2025.R", encoding = "UTF-8")

m <- read.csv(mimic_path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
d <- read_excel(domestic_path)

mimic <- data.frame(
  age = to_num(m$age),
  female = to_num(m$female),
  bmi = to_num(m$bmi),
  sofa = to_num(m$sofa),
  scm = to_num(m$scm),
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
  vasopressor = to_num(m$vasopressor)
)

sex_d <- trimws(as.character(d$性别))
domestic <- data.frame(
  age = to_num(d$年龄),
  female = ifelse(sex_d == "女", 1, ifelse(sex_d == "男", 0, NA)),
  bmi = to_num(d$bmi),
  sofa = to_num(d$sofa),
  scm = ifelse(to_num(d$射血分数) < 50, 1, 0),
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
  vasopressor = to_binary(d$升压药)
)

vars <- data.frame(
  Variable = c(
    "Age, years", "Female sex, n (%)", "BMI, kg/m2", "SOFA score", "SCM prevalence, n (%)",
    "LVEF, %", "Albumin, g/L", "Lymphocyte count, x10^9/L", "PNI",
    "Lactate, mmol/L", "Creatinine, umol/L", "WBC, x10^9/L", "Platelet count, x10^9/L", "BUN, mmol/L",
    "Coronary artery disease, n (%)", "Diabetes mellitus, n (%)", "Hypertension, n (%)", "Chronic kidney disease, n (%)", "COPD, n (%)",
    "Mechanical ventilation, n (%)", "Vasopressor use, n (%)"
  ),
  Name = c(
    "age", "female", "bmi", "sofa", "scm",
    "lvef", "albumin", "lymphocyte", "pni",
    "lactate", "creatinine", "wbc", "platelet", "bun",
    "cad", "diabetes", "hypertension", "ckd", "copd",
    "ventilation", "vasopressor"
  ),
  Type = c("cont", "cat", "cont", "cont", "cat", "cont", "cont", "cont", "cont", "cont", "cont", "cont", "cont", "cont", "cat", "cat", "cat", "cat", "cat", "cat", "cat"),
  Digits = c(1, NA, 1, 1, NA, 1, 1, 1, 1, 1, 0, 1, 1, 1, NA, NA, NA, NA, NA, NA, NA),
  stringsAsFactors = FALSE
)

rows <- list()
for (i in seq_len(nrow(vars))) {
  v <- vars[i, ]
  a <- mimic[[v$Name]]
  b <- domestic[[v$Name]]
  if (v$Type == "cont") {
    rows[[i]] <- data.frame(
      Variable = v$Variable,
      MIMIC = fmt_cont(a, v$Digits),
      Domestic = fmt_cont(b, v$Digits),
      P_value = fmt_p(p_cont(a, b)),
      Absolute_SMD = sprintf("%.3f", smd_cont(a, b)),
      stringsAsFactors = FALSE
    )
  } else {
    rows[[i]] <- data.frame(
      Variable = v$Variable,
      MIMIC = fmt_cat(a),
      Domestic = fmt_cat(b),
      P_value = fmt_p(p_cat(a, b)),
      Absolute_SMD = sprintf("%.3f", smd_cat(a, b)),
      stringsAsFactors = FALSE
    )
  }
}

heterogeneity_table <- do.call(rbind, rows)
cat("Supplementary Table 1 check: MIMIC n=", nrow(mimic), "; Domestic n=", nrow(domestic), "\n", sep = "")
print_table(heterogeneity_table)


