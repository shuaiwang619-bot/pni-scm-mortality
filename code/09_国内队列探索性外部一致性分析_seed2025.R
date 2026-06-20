set.seed(2025)

library(readxl)
library(logistf)
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

fmt_or <- function(or, lo, hi) {
  if (any(is.na(c(or, lo, hi)))) return("NA")
  sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
}

extract_direct_or <- function(fit, term) {
  if (!term %in% names(coef(fit))) stop(paste("Missing term:", term))
  data.frame(
    OR = exp(coef(fit)[term]),
    low = exp(fit$ci.lower[term]),
    high = exp(fit$ci.upper[term]),
    p = fit$prob[term],
    row.names = NULL
  )
}

fit_firth_model <- function(data, covars) {
  rhs <- c("pni_low5_c * scm_fac", covars)
  f <- as.formula(paste("death_event ~", paste(rhs, collapse = " + ")))
  vars <- unique(all.vars(f))
  d <- data[complete.cases(data[, vars]), ]

  d$scm_fac <- factor(ifelse(d$scm == 1, "SCM", "No SCM"), levels = c("No SCM", "SCM"))
  fit_no_scm_ref <- logistf(f, data = d, pl = TRUE)

  d$scm_fac <- factor(ifelse(d$scm == 1, "SCM", "No SCM"), levels = c("SCM", "No SCM"))
  fit_scm_ref <- logistf(f, data = d, pl = TRUE)

  int_name <- intersect(c("pni_low5_c:scm_facSCM", "scm_facSCM:pni_low5_c"), names(coef(fit_no_scm_ref)))[1]
  if (is.na(int_name)) stop("Interaction term not found.")

  int_idx <- which(names(coef(fit_no_scm_ref)) == int_name)
  int_test <- logistftest(fit_no_scm_ref, test = int_idx)

  no_scm <- extract_direct_or(fit_no_scm_ref, "pni_low5_c")
  scm <- extract_direct_or(fit_scm_ref, "pni_low5_c")
  interaction <- extract_direct_or(fit_no_scm_ref, int_name)

  data.frame(
    N = nrow(d),
    Events = sum(d$death_event == 1, na.rm = TRUE),
    No_SCM_OR_95CI = fmt_or(no_scm$OR, no_scm$low, no_scm$high),
    No_SCM_P = fmt_p(no_scm$p),
    SCM_OR_95CI = fmt_or(scm$OR, scm$low, scm$high),
    SCM_P = fmt_p(scm$p),
    Interaction_OR_95CI = fmt_or(interaction$OR, interaction$low, interaction$high),
    Interaction_LRT_P = fmt_p(as.numeric(int_test$prob)),
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
  stringsAsFactors = FALSE
)

dat$pni <- dat$albumin + 5 * dat$lymphocyte
dat$scm <- ifelse(dat$lvef < 50, 1, 0)
dat$death_event <- ifelse(grepl("\u597d\u8f6c", status_text), 0, 1)
dat$pni_low5 <- -dat$pni / 5
dat$pni_low5_c <- dat$pni_low5 - mean(dat$pni_low5, na.rm = TRUE)
dat$log_lactate <- ifelse(dat$lactate > 0, log(dat$lactate), NA_real_)
dat$log_creatinine <- ifelse(dat$creatinine > 0, log(dat$creatinine), NA_real_)
dat$scm_fac <- factor(ifelse(dat$scm == 1, "SCM", "No SCM"), levels = c("No SCM", "SCM"))

models <- list(
  "Model 1: unadjusted" = character(0),
  "Model 2: age, sex, BMI" = c("age", "female", "bmi"),
  "Model 3: limited clinical adjustment" = c("age", "female", "bmi", "sofa", "log_lactate", "log_creatinine")
)

rows <- list()
for (m in names(models)) {
  result <- fit_firth_model(dat, models[[m]])
  result$Model <- m
  rows[[length(rows) + 1]] <- result
}

result_table <- do.call(rbind, rows)
result_table <- result_table[, c(
  "Model", "N", "Events",
  "No_SCM_OR_95CI", "No_SCM_P",
  "SCM_OR_95CI", "SCM_P",
  "Interaction_OR_95CI", "Interaction_LRT_P"
)]

doc <- read_docx()
doc <- body_add_fpar(
  doc,
  fpar(ftext(
    "Supplementary Table S4. Exploratory external consistency analysis in the domestic cohort",
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
  Model = "Model",
  N = "N",
  Events = "Deaths",
  No_SCM_OR_95CI = "No SCM OR (95% CI)",
  No_SCM_P = "No SCM P",
  SCM_OR_95CI = "SCM OR (95% CI)",
  SCM_P = "SCM P",
  Interaction_OR_95CI = "Interaction OR (95% CI)",
  Interaction_LRT_P = "LRT P for interaction"
)
ft <- theme_booktabs(ft)
ft <- fontsize(ft, size = 8, part = "all")
ft <- bold(ft, part = "header")
ft <- align(ft, align = "center", part = "all")
ft <- align(ft, j = "Model", align = "left", part = "body")
ft <- padding(ft, padding = 3, part = "all")
ft <- width(ft, j = "Model", width = 1.55)
ft <- width(ft, j = c("N", "Events"), width = 0.42)
ft <- width(ft, j = c("No_SCM_OR_95CI", "SCM_OR_95CI", "Interaction_OR_95CI"), width = 1.05)
ft <- width(ft, j = c("No_SCM_P", "SCM_P", "Interaction_LRT_P"), width = 0.65)
ft <- autofit(ft)
ft <- fit_to_width(ft, max_width = 7.2)

doc <- body_add_flextable(doc, ft)
doc <- body_add_par(
  doc,
  "Note. Firth penalized logistic regression was used for this exploratory analysis because of the limited domestic sample size. ORs are per 5-point decrease in PNI. The interaction OR represents the ratio of the PNI-associated OR in the SCM group to that in the No SCM group; values below 1 indicate attenuation in SCM. P values for interaction were obtained using penalized likelihood-ratio tests. Model 1 was unadjusted. Model 2 adjusted for age, sex, and BMI. Model 3 adjusted for age, sex, BMI, SOFA score, log-transformed lactate, and log-transformed creatinine. BMI, body mass index; CI, confidence interval; OR, odds ratio; PNI, prognostic nutritional index; SCM, septic cardiomyopathy; SOFA, Sequential Organ Failure Assessment.",
  style = "Normal"
)

out_path <- file.path(model_result_dir, "\u56fd\u5185\u961f\u5217\u63a2\u7d22\u6027\u5916\u90e8\u4e00\u81f4\u6027\u5206\u6790.docx")
print(doc, target = out_path)

cat("Domestic exploratory external consistency analysis completed.\n")
cat("Source file:", domestic_path, "\n")
cat("Saved DOCX:", out_path, "\n")
cat("Cohort counts: N =", nrow(dat), "; deaths =", sum(dat$death_event == 1, na.rm = TRUE),
    "; SCM =", sum(dat$scm == 1, na.rm = TRUE), "; No SCM =", sum(dat$scm == 0, na.rm = TRUE), "\n")
print(result_table, row.names = FALSE)


