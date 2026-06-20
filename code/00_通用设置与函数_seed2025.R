set.seed(2025)

project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE)
analysis_dir <- project_dir
model_result_dir <- file.path(project_dir, "results")
mimic_path <- file.path(project_dir, "data_private", "mimic队列插补后.csv")
domestic_path <- file.path(project_dir, "data_private", "临沂队列插补后.xlsx")

to_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9eE+\\.-]", "", trimws(as.character(x)))))

to_binary <- function(x) {
  z <- trimws(as.character(x))
  out <- rep(NA_real_, length(z))
  out[z %in% c("1", "1.0", "是", "有", "死亡", "自动出院")] <- 1
  out[z %in% c("0", "0.0", "否", "无", "存活", "好转")] <- 0
  idx <- is.na(out) & !is.na(z) & z != ""
  out[idx] <- suppressWarnings(as.numeric(z[idx]))
  out
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt_hr <- function(hr, lo, hi) {
  if (any(is.na(c(hr, lo, hi)))) return("NA")
  sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)
}

fmt_cont <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  qs <- as.numeric(quantile(x, probs = c(0.5, 0.25, 0.75), na.rm = TRUE, type = 2))
  sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"), qs[1], qs[2], qs[3])
}

fmt_cat <- function(x) {
  x <- x[!is.na(x)]
  n <- sum(x == 1, na.rm = TRUE)
  N <- length(x)
  sprintf("%d (%.1f)", n, 100 * n / N)
}

p_cont <- function(a, b) {
  suppressWarnings(tryCatch(wilcox.test(a, b, exact = FALSE)$p.value, error = function(e) NA_real_))
}

p_cat <- function(a, b) {
  tab <- rbind(
    c(sum(a == 1, na.rm = TRUE), sum(a == 0, na.rm = TRUE)),
    c(sum(b == 1, na.rm = TRUE), sum(b == 0, na.rm = TRUE))
  )
  if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) return(NA_real_)
  expected <- suppressWarnings(chisq.test(tab)$expected)
  if (any(expected < 5)) {
    suppressWarnings(fisher.test(tab)$p.value)
  } else {
    suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
  }
}

smd_cont <- function(a, b) {
  a <- a[!is.na(a)]
  b <- b[!is.na(b)]
  sp <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  abs((mean(a) - mean(b)) / sp)
}

smd_cat <- function(a, b) {
  p1 <- mean(a == 1, na.rm = TRUE)
  p2 <- mean(b == 1, na.rm = TRUE)
  denom <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
  if (denom == 0) return(NA_real_)
  abs((p1 - p2) / denom)
}

summarize_by_scm <- function(dat, name, type, digits = 1, blank_p = FALSE) {
  x <- dat[[name]]
  g <- dat$scm
  if (type == "cont") {
    data.frame(
      Overall = fmt_cont(x, digits),
      No_SCM = fmt_cont(x[g == 0], digits),
      SCM = fmt_cont(x[g == 1], digits),
      P_value = if (blank_p) "-" else fmt_p(p_cont(x[g == 0], x[g == 1])),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Overall = fmt_cat(x),
      No_SCM = fmt_cat(x[g == 0]),
      SCM = fmt_cat(x[g == 1]),
      P_value = if (blank_p) "-" else fmt_p(p_cat(x[g == 0], x[g == 1])),
      stringsAsFactors = FALSE
    )
  }
}

load_mimic_model_data <- function() {
  dat <- read.csv(mimic_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
  num_cols <- c(
    "age", "female", "bmi", "sofa", "scm", "pni", "lactate", "creatinine",
    "coronary_artery_disease", "diabetes_mellitus", "hypertension",
    "chronic_kidney_disease", "copd", "mechanical_vent", "vasopressor",
    "alive_at_landmark", "event_28d", "time_28d", "event_90d", "time_90d"
  )
  for (cc in intersect(num_cols, names(dat))) dat[[cc]] <- to_num(dat[[cc]])
  dat$pni_low5 <- -dat$pni / 5
  dat$pni_low5_c <- dat$pni_low5 - mean(dat$pni_low5, na.rm = TRUE)
  dat$log_lactate <- ifelse(dat$lactate > 0, log(dat$lactate), NA_real_)
  dat$log_creatinine <- ifelse(dat$creatinine > 0, log(dat$creatinine), NA_real_)
  dat$event_28_90 <- ifelse(dat$time_90d > 28 & dat$event_90d == 1, 1, 0)
  dat$time_28_90 <- ifelse(dat$time_90d > 28, dat$time_90d - 28, NA_real_)
  dat
}

main_covars <- list(
  "Model 1: unadjusted" = character(0),
  "Model 2: age, sex, BMI" = c("age", "female", "bmi"),
  "Model 3: fully adjusted" = c(
    "age", "female", "bmi", "sofa", "log_lactate", "log_creatinine",
    "coronary_artery_disease", "diabetes_mellitus", "hypertension",
    "chronic_kidney_disease", "copd", "mechanical_vent", "vasopressor"
  )
)

sensitivity_covars <- list(
  "Sensitivity 1: full model without organ support" = c(
    "age", "female", "bmi", "sofa", "log_lactate", "log_creatinine",
    "coronary_artery_disease", "diabetes_mellitus", "hypertension",
    "chronic_kidney_disease", "copd"
  ),
  "Sensitivity 2: demographics and comorbidities" = c(
    "age", "female", "bmi",
    "coronary_artery_disease", "diabetes_mellitus", "hypertension",
    "chronic_kidney_disease", "copd"
  )
)

outcomes <- list(
  "28-day mortality" = c(time = "time_28d", event = "event_28d"),
  "90-day mortality" = c(time = "time_90d", event = "event_90d")
)

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

fit_one_cox <- function(data, time_col, event_col, covars) {
  rhs_terms <- c("pni_low5_c * scm", covars)
  f <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", paste(rhs_terms, collapse = " + ")))
  vars <- unique(all.vars(f))
  d <- data[complete.cases(data[, vars]), ]
  d <- d[d[[time_col]] > 0 & !is.na(d[[event_col]]), ]
  fit <- survival::coxph(f, data = d, ties = "efron", x = TRUE)
  reduced_f <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ",
    paste(c("pni_low5_c", "scm", covars), collapse = " + ")
  ))
  reduced <- survival::coxph(reduced_f, data = d, ties = "efron", x = TRUE)
  lrt <- anova(reduced, fit, test = "LRT")
  list(fit = fit, data = d, lrt_p = lrt$`Pr(>|Chi|)`[2])
}

summarize_cox_fit <- function(fit_object, outcome_label, model_label, event_col) {
  fit <- fit_object$fit
  d <- fit_object$data
  int_name <- intersect(c("pni_low5_c:scm", "scm:pni_low5_c"), names(coef(fit)))[1]
  no_scm <- lincomb(fit, c("pni_low5_c"))
  scm <- lincomb(fit, c("pni_low5_c", int_name))
  interaction <- lincomb(fit, c(int_name))
  data.frame(
    Outcome = outcome_label,
    Model = model_label,
    N = nrow(d),
    Events = sum(d[[event_col]] == 1, na.rm = TRUE),
    No_SCM_HR_95CI = fmt_hr(no_scm$hr, no_scm$conf.low, no_scm$conf.high),
    No_SCM_P = fmt_p(no_scm$p.value),
    SCM_HR_95CI = fmt_hr(scm$hr, scm$conf.low, scm$conf.high),
    SCM_P = fmt_p(scm$p.value),
    Interaction_HR_95CI = fmt_hr(interaction$hr, interaction$conf.low, interaction$conf.high),
    Interaction_LRT_P = fmt_p(fit_object$lrt_p),
    row.names = NULL
  )
}

print_table <- function(x) {
  print(x, row.names = FALSE)
}



