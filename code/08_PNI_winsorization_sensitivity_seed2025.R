set.seed(2025)

library(survival)
library(splines)
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
    hr = exp(beta),
    conf.low = exp(beta - 1.96 * se),
    conf.high = exp(beta + 1.96 * se),
    p.value = 2 * pnorm(abs(beta / se), lower.tail = FALSE)
  )
}

dat <- read.csv(mimic_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))

num_cols <- c(
  "age", "female", "bmi", "sofa", "scm", "pni", "lactate", "creatinine",
  "coronary_artery_disease", "diabetes_mellitus", "hypertension",
  "chronic_kidney_disease", "copd", "mechanical_vent", "vasopressor",
  "event_28d", "time_28d", "event_90d", "time_90d"
)
for (cc in intersect(num_cols, names(dat))) dat[[cc]] <- to_num(dat[[cc]])

dat$log_lactate <- ifelse(dat$lactate > 0, log(dat$lactate), NA_real_)
dat$log_creatinine <- ifelse(dat$creatinine > 0, log(dat$creatinine), NA_real_)

pni_limits <- as.numeric(quantile(dat$pni, probs = c(0.01, 0.99), na.rm = TRUE, type = 2))
dat$pni_winsor <- pmin(pmax(dat$pni, pni_limits[1]), pni_limits[2])
dat$pni_low5 <- -dat$pni / 5
dat$pni_low5_c <- dat$pni_low5 - mean(dat$pni_low5, na.rm = TRUE)
dat$pni_winsor_low5 <- -dat$pni_winsor / 5
dat$pni_winsor_low5_c <- dat$pni_winsor_low5 - mean(dat$pni_winsor_low5, na.rm = TRUE)

n_low <- sum(dat$pni < pni_limits[1], na.rm = TRUE)
n_high <- sum(dat$pni > pni_limits[2], na.rm = TRUE)

covars <- c(
  "age", "female", "bmi", "sofa", "log_lactate", "log_creatinine",
  "coronary_artery_disease", "diabetes_mellitus", "hypertension",
  "chronic_kidney_disease", "copd", "mechanical_vent", "vasopressor"
)

outcomes <- list(
  "28-day mortality" = c(time = "time_28d", event = "event_28d"),
  "90-day mortality" = c(time = "time_90d", event = "event_90d")
)

fit_linear_interaction <- function(data, time_col, event_col, pni_low5_col) {
  rhs <- c(paste0(pni_low5_col, " * scm"), covars)
  f <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", paste(rhs, collapse = " + ")))
  vars <- unique(all.vars(f))
  d <- data[complete.cases(data[, vars]), ]
  d <- d[d[[time_col]] > 0 & !is.na(d[[event_col]]), ]

  fit <- coxph(f, data = d, ties = "efron", x = TRUE)
  reduced_f <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ",
    paste(c(pni_low5_col, "scm", covars), collapse = " + ")
  ))
  reduced <- coxph(reduced_f, data = d, ties = "efron", x = TRUE)
  lrt <- anova(reduced, fit, test = "LRT")

  int_name <- intersect(c(paste0(pni_low5_col, ":scm"), paste0("scm:", pni_low5_col)), names(coef(fit)))[1]
  no_scm <- lincomb(fit, c(pni_low5_col))
  scm <- lincomb(fit, c(pni_low5_col, int_name))
  interaction <- lincomb(fit, c(int_name))

  list(
    N = nrow(d),
    Events = sum(d[[event_col]] == 1, na.rm = TRUE),
    No_SCM = no_scm,
    SCM = scm,
    Interaction = interaction,
    Linear_interaction_P = lrt$`Pr(>|Chi|)`[2]
  )
}

fit_spline_tests <- function(data, time_col, event_col, pni_col) {
  vars <- c(time_col, event_col, pni_col, "scm", covars)
  d <- data[complete.cases(data[, vars]), ]
  d <- d[d[[time_col]] > 0 & !is.na(d[[event_col]]), ]

  linear_f <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ", pni_col, " * scm + ",
    paste(covars, collapse = " + ")
  ))
  spline_main_f <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ns(", pni_col, ", df = 3) + scm + ",
    paste(covars, collapse = " + ")
  ))
  spline_int_f <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ns(", pni_col, ", df = 3) * scm + ",
    paste(covars, collapse = " + ")
  ))

  linear_fit <- coxph(linear_f, data = d, ties = "efron", x = TRUE)
  spline_main_fit <- coxph(spline_main_f, data = d, ties = "efron", x = TRUE)
  spline_int_fit <- coxph(spline_int_f, data = d, ties = "efron", x = TRUE)

  nonlin_lrt <- anova(linear_fit, spline_int_fit, test = "LRT")
  interaction_lrt <- anova(spline_main_fit, spline_int_fit, test = "LRT")

  list(
    Nonlinearity_P = nonlin_lrt$`Pr(>|Chi|)`[2],
    Spline_interaction_P = interaction_lrt$`Pr(>|Chi|)`[2]
  )
}

build_row <- function(outcome_label, pni_handling, pni_col, pni_low5_col) {
  time_col <- outcomes[[outcome_label]]["time"]
  event_col <- outcomes[[outcome_label]]["event"]
  lin <- fit_linear_interaction(dat, time_col, event_col, pni_low5_col)
  spl <- fit_spline_tests(dat, time_col, event_col, pni_col)

  data.frame(
    Outcome = outcome_label,
    PNI_handling = pni_handling,
    N = lin$N,
    Events = lin$Events,
    No_SCM_HR_95CI = fmt_hr(lin$No_SCM$hr, lin$No_SCM$conf.low, lin$No_SCM$conf.high),
    SCM_HR_95CI = fmt_hr(lin$SCM$hr, lin$SCM$conf.low, lin$SCM$conf.high),
    Interaction_HR_95CI = fmt_hr(lin$Interaction$hr, lin$Interaction$conf.low, lin$Interaction$conf.high),
    Linear_interaction_P = fmt_p(lin$Linear_interaction_P),
    Nonlinearity_P = fmt_p(spl$Nonlinearity_P),
    Spline_interaction_P = fmt_p(spl$Spline_interaction_P),
    stringsAsFactors = FALSE
  )
}

rows <- list()
for (label in names(outcomes)) {
  rows[[length(rows) + 1]] <- build_row(label, "Original PNI", "pni", "pni_low5_c")
  rows[[length(rows) + 1]] <- build_row(label, "Winsorized PNI", "pni_winsor", "pni_winsor_low5_c")
}
result_table <- do.call(rbind, rows)

section <- prop_section(
  page_size = page_size(orient = "landscape"),
  page_margins = page_mar(top = 0.55, bottom = 0.55, left = 0.45, right = 0.45)
)

doc <- read_docx()
doc <- body_set_default_section(doc, section)
doc <- body_add_fpar(
  doc,
  fpar(ftext(
    "Supplementary Table S3. Sensitivity analysis using winsorized PNI",
    prop = fp_text(font.size = 16, bold = TRUE)
  ))
)
doc <- body_add_par(
  doc,
  sprintf(
    "PNI was winsorized at the 1st and 99th percentiles (%.2f and %.2f). No patients were removed; %d low-end and %d high-end observations were capped.",
    pni_limits[1], pni_limits[2], n_low, n_high
  ),
  style = "Normal"
)

ft <- flextable(result_table)
ft <- set_header_labels(
  ft,
  Outcome = "Outcome",
  PNI_handling = "PNI handling",
  N = "N",
  Events = "Events",
  No_SCM_HR_95CI = "No SCM HR (95% CI)",
  SCM_HR_95CI = "SCM HR (95% CI)",
  Interaction_HR_95CI = "Interaction HR (95% CI)",
  Linear_interaction_P = "Linear interaction P",
  Nonlinearity_P = "Nonlinearity P",
  Spline_interaction_P = "Spline interaction P"
)
ft <- theme_booktabs(ft)
ft <- fontsize(ft, size = 8, part = "all")
ft <- bold(ft, part = "header")
ft <- align(ft, align = "center", part = "all")
ft <- align(ft, j = c("Outcome", "PNI_handling"), align = "left", part = "body")
ft <- padding(ft, padding = 3, part = "all")
ft <- width(ft, j = "Outcome", width = 1.05)
ft <- width(ft, j = "PNI_handling", width = 1.05)
ft <- width(ft, j = c("N", "Events"), width = 0.45)
ft <- width(ft, j = c("No_SCM_HR_95CI", "SCM_HR_95CI", "Interaction_HR_95CI"), width = 1.05)
ft <- width(ft, j = c("Linear_interaction_P", "Nonlinearity_P", "Spline_interaction_P"), width = 0.85)
ft <- autofit(ft)

doc <- body_add_flextable(doc, ft)
doc <- body_add_par(
  doc,
  "Note. Winsorization capped values outside the 1st and 99th percentiles without excluding patients. HRs are per 5-point decrease in PNI. Linear interaction P values were obtained by likelihood-ratio tests comparing fully adjusted Cox models with and without the 5-point lower PNI x SCM interaction. Nonlinearity P values compare the natural cubic spline model with the linear PNI model, both including the SCM interaction. Spline interaction P values test whether the PNI spline differed by SCM status. Models were adjusted as in Model 3: age, sex, BMI, SOFA score, log-transformed lactate, log-transformed creatinine, coronary artery disease, diabetes mellitus, hypertension, chronic kidney disease, COPD, mechanical ventilation, and vasopressor use.",
  style = "Normal"
)

out_path <- file.path(model_result_dir, "PNI_winsorization_sensitivity_seed2025.docx")
print(doc, target = out_path)

cat("PNI winsorization sensitivity completed.\n")
cat("Winsorization limits:", sprintf("%.4f", pni_limits[1]), sprintf("%.4f", pni_limits[2]), "\n")
print(result_table, row.names = FALSE)
cat("Saved DOCX:", out_path, "\n")


