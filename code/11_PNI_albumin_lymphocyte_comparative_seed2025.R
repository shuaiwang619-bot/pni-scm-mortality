set.seed(2025)

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

mimic_candidates <- list.files(analysis_dir, pattern = "^mimic.*\\.csv$", full.names = TRUE)
if (length(mimic_candidates) == 0) stop("MIMIC CSV not found in analysis directory.")
mimic_path <- mimic_candidates[1]

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
    beta = beta,
    se = se,
    hr = exp(beta),
    conf.low = exp(beta - 1.96 * se),
    conf.high = exp(beta + 1.96 * se),
    p.value = 2 * pnorm(abs(beta / se), lower.tail = FALSE)
  )
}

interaction_name <- function(marker, fit) {
  out <- intersect(c(paste0(marker, ":scm"), paste0("scm:", marker)), names(coef(fit)))[1]
  if (is.na(out)) stop(paste("Interaction term not found for", marker))
  out
}

summarize_marker_from_fit <- function(fit, d, event_col, marker, model_label, marker_label, lrt_p) {
  int_name <- interaction_name(marker, fit)
  no_scm <- lincomb(fit, marker)
  scm <- lincomb(fit, c(marker, int_name))
  inter <- lincomb(fit, int_name)
  data.frame(
    Model = model_label,
    Marker = marker_label,
    N = nrow(d),
    Events = sum(d[[event_col]] == 1, na.rm = TRUE),
    No_SCM_HR_95CI = fmt_hr(no_scm$hr, no_scm$conf.low, no_scm$conf.high),
    No_SCM_P = fmt_p(no_scm$p.value),
    SCM_HR_95CI = fmt_hr(scm$hr, scm$conf.low, scm$conf.high),
    SCM_P = fmt_p(scm$p.value),
    Interaction_HR_95CI = fmt_hr(inter$hr, inter$conf.low, inter$conf.high),
    Interaction_LRT_P = fmt_p(lrt_p),
    AIC = sprintf("%.1f", AIC(fit)),
    C_index = sprintf("%.3f", as.numeric(summary(fit)$concordance[1])),
    stringsAsFactors = FALSE
  )
}

fit_single_marker <- function(data, time_col, event_col, marker, model_label, marker_label, covars) {
  f <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ",
    paste(c(paste0(marker, " * scm"), covars), collapse = " + ")
  ))
  vars <- unique(all.vars(f))
  d <- data[complete.cases(data[, vars]), ]
  d <- d[d[[time_col]] > 0 & !is.na(d[[event_col]]), ]
  fit <- coxph(f, data = d, ties = "efron", x = TRUE)

  reduced_f <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ",
    paste(c(marker, "scm", covars), collapse = " + ")
  ))
  reduced <- coxph(reduced_f, data = d, ties = "efron", x = TRUE)
  lrt_p <- anova(reduced, fit, test = "LRT")$`Pr(>|Chi|)`[2]

  summarize_marker_from_fit(fit, d, event_col, marker, model_label, marker_label, lrt_p)
}

fit_combined_markers <- function(data, time_col, event_col, covars) {
  f <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ",
    paste(c("albumin_low5_c * scm", "lymph_low05_c * scm", covars), collapse = " + ")
  ))
  vars <- unique(all.vars(f))
  d <- data[complete.cases(data[, vars]), ]
  d <- d[d[[time_col]] > 0 & !is.na(d[[event_col]]), ]
  fit <- coxph(f, data = d, ties = "efron", x = TRUE)

  reduced_albumin_f <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ",
    paste(c("albumin_low5_c", "scm", "lymph_low05_c * scm", covars), collapse = " + ")
  ))
  reduced_lymph_f <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ",
    paste(c("lymph_low05_c", "scm", "albumin_low5_c * scm", covars), collapse = " + ")
  ))
  reduced_albumin <- coxph(reduced_albumin_f, data = d, ties = "efron", x = TRUE)
  reduced_lymph <- coxph(reduced_lymph_f, data = d, ties = "efron", x = TRUE)
  lrt_albumin <- anova(reduced_albumin, fit, test = "LRT")$`Pr(>|Chi|)`[2]
  lrt_lymph <- anova(reduced_lymph, fit, test = "LRT")$`Pr(>|Chi|)`[2]

  rbind(
    summarize_marker_from_fit(fit, d, event_col, "albumin_low5_c", "Albumin + lymphocyte", "Albumin, per 5 g/L lower", lrt_albumin),
    summarize_marker_from_fit(fit, d, event_col, "lymph_low05_c", "Albumin + lymphocyte", "Lymphocyte, per 0.5 x10^9/L lower", lrt_lymph)
  )
}

dat <- read.csv(mimic_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))

num_cols <- c(
  "age", "female", "bmi", "sofa", "scm", "pni", "albumin_g_l", "lymphocyte_k_ul",
  "lactate", "creatinine", "coronary_artery_disease", "diabetes_mellitus",
  "hypertension", "chronic_kidney_disease", "copd", "mechanical_vent", "vasopressor",
  "event_28d", "time_28d", "event_90d", "time_90d"
)
for (cc in intersect(num_cols, names(dat))) dat[[cc]] <- to_num(dat[[cc]])

dat$pni_low5 <- -dat$pni / 5
dat$albumin_low5 <- -dat$albumin_g_l / 5
dat$lymph_low05 <- -dat$lymphocyte_k_ul / 0.5
dat$pni_low5_c <- dat$pni_low5 - mean(dat$pni_low5, na.rm = TRUE)
dat$albumin_low5_c <- dat$albumin_low5 - mean(dat$albumin_low5, na.rm = TRUE)
dat$lymph_low05_c <- dat$lymph_low05 - mean(dat$lymph_low05, na.rm = TRUE)
dat$log_lactate <- ifelse(dat$lactate > 0, log(dat$lactate), NA_real_)
dat$log_creatinine <- ifelse(dat$creatinine > 0, log(dat$creatinine), NA_real_)

covars <- c(
  "age", "female", "bmi", "sofa", "log_lactate", "log_creatinine",
  "coronary_artery_disease", "diabetes_mellitus", "hypertension",
  "chronic_kidney_disease", "copd", "mechanical_vent", "vasopressor"
)

outcomes <- list(
  "28-day mortality" = c(time = "time_28d", event = "event_28d"),
  "90-day mortality" = c(time = "time_90d", event = "event_90d")
)

rows <- list()
for (outcome_label in names(outcomes)) {
  time_col <- outcomes[[outcome_label]]["time"]
  event_col <- outcomes[[outcome_label]]["event"]
  outcome_rows <- rbind(
    fit_single_marker(dat, time_col, event_col, "pni_low5_c", "PNI", "PNI, per 5-point lower", covars),
    fit_single_marker(dat, time_col, event_col, "albumin_low5_c", "Albumin", "Albumin, per 5 g/L lower", covars),
    fit_single_marker(dat, time_col, event_col, "lymph_low05_c", "Lymphocyte", "Lymphocyte, per 0.5 x10^9/L lower", covars),
    fit_combined_markers(dat, time_col, event_col, covars)
  )
  outcome_rows$Outcome <- outcome_label
  rows[[length(rows) + 1]] <- outcome_rows
}

result_table <- do.call(rbind, rows)
result_table <- result_table[, c(
  "Outcome", "Model", "Marker", "N", "Events",
  "No_SCM_HR_95CI", "No_SCM_P",
  "SCM_HR_95CI", "SCM_P",
  "Interaction_HR_95CI", "Interaction_LRT_P",
  "AIC", "C_index"
)]

section <- prop_section(
  page_size = page_size(orient = "landscape"),
  page_margins = page_mar(top = 0.55, bottom = 0.55, left = 0.45, right = 0.45)
)

make_ft <- function(tbl) {
  ft <- flextable(tbl[, c(
    "Model", "Marker", "N", "Events",
    "No_SCM_HR_95CI", "No_SCM_P",
    "SCM_HR_95CI", "SCM_P",
    "Interaction_HR_95CI", "Interaction_LRT_P",
    "AIC", "C_index"
  )])
  ft <- set_header_labels(
    ft,
    Model = "Model",
    Marker = "Marker",
    N = "N",
    Events = "Events",
    No_SCM_HR_95CI = "No SCM HR (95% CI)",
    No_SCM_P = "No SCM P",
    SCM_HR_95CI = "SCM HR (95% CI)",
    SCM_P = "SCM P",
    Interaction_HR_95CI = "Interaction HR (95% CI)",
    Interaction_LRT_P = "LRT P for interaction",
    AIC = "AIC",
    C_index = "C-index"
  )
  ft <- theme_booktabs(ft)
  ft <- fontsize(ft, size = 6.8, part = "all")
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "center", part = "all")
  ft <- align(ft, j = c("Model", "Marker"), align = "left", part = "body")
  ft <- padding(ft, padding = 2, part = "all")
  ft <- width(ft, j = "Model", width = 1.05)
  ft <- width(ft, j = "Marker", width = 1.55)
  ft <- width(ft, j = c("N", "Events"), width = 0.35)
  ft <- width(ft, j = c("No_SCM_HR_95CI", "SCM_HR_95CI", "Interaction_HR_95CI"), width = 0.95)
  ft <- width(ft, j = c("No_SCM_P", "SCM_P", "Interaction_LRT_P"), width = 0.52)
  ft <- width(ft, j = c("AIC", "C_index"), width = 0.48)
  ft <- merge_v(ft, j = "Model")
  ft <- valign(ft, valign = "center", part = "body")
  ft <- autofit(ft)
  fit_to_width(ft, max_width = 9.4)
}

doc <- read_docx()
doc <- body_set_default_section(doc, section)
doc <- body_add_fpar(
  doc,
  fpar(ftext(
    "Supplementary Table S6. Comparative prognostic analyses of PNI and its individual components in the MIMIC-IV cohort",
    prop = fp_text(font.size = 11, bold = TRUE)
  ))
)
doc <- body_add_par(
  doc,
  "Fully adjusted Cox models were fitted separately for PNI, albumin, lymphocyte count, and albumin plus lymphocyte count.",
  style = "Normal"
)

for (outcome_label in names(outcomes)) {
  doc <- body_add_par(doc, outcome_label, style = "Normal")
  doc <- body_add_flextable(doc, make_ft(result_table[result_table$Outcome == outcome_label, ]))
  doc <- body_add_par(doc, "", style = "Normal")
}

doc <- body_add_par(
  doc,
  "Note. HRs are scaled as 5-point lower PNI, 5 g/L lower albumin, or 0.5 x10^9/L lower lymphocyte count. All markers were centered before interaction modeling. Models were adjusted for age, sex, BMI, SOFA score, log-transformed lactate, log-transformed creatinine, coronary artery disease, diabetes mellitus, hypertension, chronic kidney disease, COPD, mechanical ventilation, and vasopressor use. PNI, albumin, and lymphocyte count were not included together except in the dedicated albumin plus lymphocyte model. In the combined model, LRT P values test each component-specific interaction while retaining the other component and its interaction. Lower AIC and higher C-index indicate better model fit/discrimination. AIC, Akaike information criterion; BMI, body mass index; CI, confidence interval; COPD, chronic obstructive pulmonary disease; HR, hazard ratio; PNI, prognostic nutritional index; SCM, septic cardiomyopathy; SOFA, Sequential Organ Failure Assessment.",
  style = "Normal"
)

out_path <- file.path(model_result_dir, "\u0050\u004e\u0049\u4e0e\u5355\u72ec\u6210\u5206\u6bd4\u8f83\u5206\u6790.docx")
print(doc, target = out_path)

cat("PNI vs albumin vs lymphocyte comparative analysis completed.\n")
cat("Source file:", mimic_path, "\n")
cat("Saved DOCX:", out_path, "\n")
print(result_table, row.names = FALSE)


