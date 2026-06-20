set.seed(2025)

library(readxl)
library(survival)
library(officer)
library(flextable)

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
analysis_dir <- dirname(script_dir)

domestic_candidates <- list.files(analysis_dir, pattern = "\\.xlsx$", full.names = TRUE)
if (length(domestic_candidates) == 0) stop("Domestic cohort XLSX not found in analysis directory.")
domestic_path <- domestic_candidates[1]

subdirs <- list.dirs(analysis_dir, recursive = FALSE, full.names = TRUE)
docx_counts <- vapply(subdirs, function(x) length(list.files(x, pattern = "\\.docx$")), integer(1))
if (length(docx_counts) == 0 || max(docx_counts) == 0) {
  model_result_dir <- analysis_dir
} else {
  model_result_dir <- subdirs[which.max(docx_counts)]
}

to_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9eE+\\.-]", "", trimws(as.character(x)))))

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt_hr <- function(hr, lo, hi) {
  if (any(is.na(c(hr, lo, hi)))) return("NA")
  sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)
}

lincomb <- function(fit, terms) {
  b <- coef(fit)
  v <- vcov(fit)
  l <- rep(0, length(b))
  names(l) <- names(b)
  missing_terms <- setdiff(terms, names(b))
  if (length(missing_terms) > 0) stop(paste("Missing term(s):", paste(missing_terms, collapse = ", ")))
  l[terms] <- 1
  beta <- sum(l * b)
  se <- sqrt(as.numeric(t(l) %*% v %*% l))
  data.frame(
    hr = exp(beta),
    conf.low = exp(beta - 1.96 * se),
    conf.high = exp(beta + 1.96 * se),
    p.value = 2 * pnorm(abs(beta / se), lower.tail = FALSE)
  )
}

fit_one_cox <- function(data, time_col, event_col, covars) {
  rhs_terms <- c("pni_low5_c * scm", covars)
  f <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", paste(rhs_terms, collapse = " + ")))
  vars <- unique(all.vars(f))
  d <- data[complete.cases(data[, vars]), ]
  d <- d[d[[time_col]] > 0 & !is.na(d[[event_col]]), ]

  fit <- coxph(f, data = d, ties = "efron", x = TRUE)
  reduced_f <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ",
    paste(c("pni_low5_c", "scm", covars), collapse = " + ")
  ))
  reduced <- coxph(reduced_f, data = d, ties = "efron", x = TRUE)
  lrt <- anova(reduced, fit, test = "LRT")

  int_name <- intersect(c("pni_low5_c:scm", "scm:pni_low5_c"), names(coef(fit)))[1]
  no_scm <- lincomb(fit, "pni_low5_c")
  scm <- lincomb(fit, c("pni_low5_c", int_name))
  interaction <- lincomb(fit, int_name)

  data.frame(
    N = nrow(d),
    Deaths = sum(d[[event_col]] == 1, na.rm = TRUE),
    No_SCM_HR_95CI = fmt_hr(no_scm$hr, no_scm$conf.low, no_scm$conf.high),
    No_SCM_P = fmt_p(no_scm$p.value),
    SCM_HR_95CI = fmt_hr(scm$hr, scm$conf.low, scm$conf.high),
    SCM_P = fmt_p(scm$p.value),
    Interaction_HR_95CI = fmt_hr(interaction$hr, interaction$conf.low, interaction$conf.high),
    Interaction_LRT_P = fmt_p(lrt$`Pr(>|Chi|)`[2]),
    stringsAsFactors = FALSE
  )
}

raw <- read_excel(domestic_path)
if (ncol(raw) < 50) stop("Domestic cohort file has fewer columns than expected.")

sex_text <- trimws(as.character(raw[[5]]))
status_text <- trimws(as.character(raw[[48]]))

dat <- data.frame(
  age = to_num(raw[[6]]),
  female = ifelse(grepl("\u5973", sex_text), 1, ifelse(grepl("\u7537", sex_text), 0, NA_real_)),
  bmi = to_num(raw[[7]]),
  sofa = to_num(raw[[19]]),
  lactate = to_num(raw[[24]]),
  creatinine = to_num(raw[[22]]),
  lvef = to_num(raw[[25]]),
  lymphocyte = to_num(raw[[28]]),
  albumin = to_num(raw[[33]]),
  hospital_days = to_num(raw[[49]]),
  stringsAsFactors = FALSE
)

dat$pni <- dat$albumin + 5 * dat$lymphocyte
dat$scm <- ifelse(dat$lvef < 50, 1, 0)
dat$death_event <- ifelse(grepl("\u597d\u8f6c", status_text), 0, 1)
dat$pni_low5 <- -dat$pni / 5
dat$pni_low5_c <- dat$pni_low5 - mean(dat$pni_low5, na.rm = TRUE)
dat$log_lactate <- ifelse(dat$lactate > 0, log(dat$lactate), NA_real_)
dat$log_creatinine <- ifelse(dat$creatinine > 0, log(dat$creatinine), NA_real_)

dat$time_inhospital <- ifelse(dat$hospital_days > 0, dat$hospital_days, NA_real_)
dat$event_inhospital <- dat$death_event
dat$time_28_proxy <- ifelse(!is.na(dat$time_inhospital), pmin(dat$time_inhospital, 28), NA_real_)
dat$event_28_proxy <- ifelse(dat$death_event == 1 & dat$time_inhospital <= 28, 1, 0)

models <- list(
  "Model 1: unadjusted" = character(0),
  "Model 2: age, sex, BMI" = c("age", "female", "bmi"),
  "Model 3: limited clinical adjustment" = c("age", "female", "bmi", "sofa", "log_lactate", "log_creatinine")
)

scenarios <- list(
  "In-hospital Cox using hospital days" = c(time = "time_inhospital", event = "event_inhospital"),
  "28-day proxy Cox censored at 28 hospital days" = c(time = "time_28_proxy", event = "event_28_proxy")
)

rows <- list()
for (scenario in names(scenarios)) {
  time_col <- scenarios[[scenario]]["time"]
  event_col <- scenarios[[scenario]]["event"]
  for (model_name in names(models)) {
    result <- fit_one_cox(dat, time_col, event_col, models[[model_name]])
    result$Scenario <- scenario
    result$Model <- model_name
    rows[[length(rows) + 1]] <- result
  }
}

result_table <- do.call(rbind, rows)
result_table <- result_table[, c(
  "Scenario", "Model", "N", "Deaths",
  "No_SCM_HR_95CI", "No_SCM_P",
  "SCM_HR_95CI", "SCM_P",
  "Interaction_HR_95CI", "Interaction_LRT_P"
)]

section <- prop_section(
  page_size = page_size(orient = "landscape"),
  page_margins = page_mar(top = 0.55, bottom = 0.55, left = 0.45, right = 0.45)
)

doc <- read_docx()
doc <- body_set_default_section(doc, section)
doc <- body_add_fpar(
  doc,
  fpar(ftext(
    "Supplementary Table S5. Exploratory Cox HR analysis in the domestic cohort",
    prop = fp_text(font.size = 11, bold = TRUE)
  ))
)
doc <- body_add_par(
  doc,
  sprintf(
    "Domestic cohort: n=%d, deaths=%d, SCM=%d, No SCM=%d.",
    nrow(dat), sum(dat$death_event == 1, na.rm = TRUE), sum(dat$scm == 1, na.rm = TRUE), sum(dat$scm == 0, na.rm = TRUE)
  ),
  style = "Normal"
)

ft <- flextable(result_table)
ft <- set_header_labels(
  ft,
  Scenario = "Scenario",
  Model = "Model",
  N = "N",
  Deaths = "Deaths",
  No_SCM_HR_95CI = "No SCM HR (95% CI)",
  No_SCM_P = "No SCM P",
  SCM_HR_95CI = "SCM HR (95% CI)",
  SCM_P = "SCM P",
  Interaction_HR_95CI = "Interaction HR (95% CI)",
  Interaction_LRT_P = "LRT P for interaction"
)
ft <- theme_booktabs(ft)
ft <- fontsize(ft, size = 7, part = "all")
ft <- bold(ft, part = "header")
ft <- align(ft, align = "center", part = "all")
ft <- align(ft, j = c("Scenario", "Model"), align = "left", part = "body")
ft <- padding(ft, padding = 2, part = "all")
ft <- width(ft, j = "Scenario", width = 1.9)
ft <- width(ft, j = "Model", width = 1.6)
ft <- width(ft, j = c("N", "Deaths"), width = 0.35)
ft <- width(ft, j = c("No_SCM_HR_95CI", "SCM_HR_95CI", "Interaction_HR_95CI"), width = 1.05)
ft <- width(ft, j = c("No_SCM_P", "SCM_P", "Interaction_LRT_P"), width = 0.55)
ft <- merge_v(ft, j = "Scenario")
ft <- valign(ft, valign = "center", part = "body")
ft <- autofit(ft)
ft <- fit_to_width(ft, max_width = 9.4)

doc <- body_add_flextable(doc, ft)
doc <- body_add_par(
  doc,
  "Note. This exploratory Cox analysis used hospital days as the time scale and should not be interpreted as formal validation of the MIMIC-IV 28-day or 90-day time-to-event models. HRs are per 5-point decrease in PNI. The interaction HR represents the ratio of the PNI-associated HR in the SCM group to that in the No SCM group; values below 1 indicate attenuation in SCM. The 28-day proxy analysis censored follow-up at 28 hospital days. P values for interaction were obtained by likelihood-ratio tests comparing models with and without the 5-point lower PNI x SCM interaction. Model 1 was unadjusted. Model 2 adjusted for age, sex, and BMI. Model 3 adjusted for age, sex, BMI, SOFA score, log-transformed lactate, and log-transformed creatinine. BMI, body mass index; CI, confidence interval; HR, hazard ratio; PNI, prognostic nutritional index; SCM, septic cardiomyopathy; SOFA, Sequential Organ Failure Assessment.",
  style = "Normal"
)

out_path <- file.path(model_result_dir, "\u56fd\u5185\u961f\u5217\u63a2\u7d22\u6027Cox_HR\u8868.docx")
print(doc, target = out_path)

cat("Domestic exploratory Cox HR analysis completed.\n")
cat("Source file:", domestic_path, "\n")
cat("Saved DOCX:", out_path, "\n")
cat("Cohort counts: N =", nrow(dat), "; deaths =", sum(dat$death_event == 1, na.rm = TRUE),
    "; SCM =", sum(dat$scm == 1, na.rm = TRUE), "; No SCM =", sum(dat$scm == 0, na.rm = TRUE), "\n")
print(result_table, row.names = FALSE)


