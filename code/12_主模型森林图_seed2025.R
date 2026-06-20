set.seed(2025)

library(survival)
library(readxl)
library(logistf)
library(ggplot2)

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

domestic_preferred <- file.path(analysis_dir, "\u4e34\u6c82\u961f\u5217\u63d2\u8865\u540e.xlsx")
if (file.exists(domestic_preferred)) {
  domestic_path <- domestic_preferred
} else {
  domestic_candidates <- list.files(analysis_dir, pattern = "\\.xlsx$", full.names = TRUE)
  if (length(domestic_candidates) == 0) stop("Domestic cohort XLSX not found in analysis directory.")
  domestic_path <- domestic_candidates[1]
}

subdirs <- list.dirs(analysis_dir, recursive = FALSE, full.names = TRUE)
docx_counts <- vapply(subdirs, function(x) length(list.files(x, pattern = "\\.docx$")), integer(1))
if (length(docx_counts) == 0 || max(docx_counts) == 0) {
  model_result_dir <- analysis_dir
} else {
  model_result_dir <- subdirs[which.max(docx_counts)]
}

to_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9eE+\\.-]", "", trimws(as.character(x)))))

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

fmt_hr <- function(hr, lo, hi) sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)

fit_model3 <- function(data, time_col, event_col, covars) {
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
  lrt_p <- anova(reduced, fit, test = "LRT")$`Pr(>|Chi|)`[2]
  list(fit = fit, data = d, lrt_p = lrt_p)
}

build_cox_rows <- function(fit_object, outcome_label, event_col) {
  fit <- fit_object$fit
  int_name <- intersect(c("pni_low5_c:scm", "scm:pni_low5_c"), names(coef(fit)))[1]
  no_scm <- lincomb(fit, "pni_low5_c")
  scm <- lincomb(fit, c("pni_low5_c", int_name))
  interaction <- lincomb(fit, int_name)

  data.frame(
    outcome = outcome_label,
    estimand = c("No SCM", "SCM", "Interaction"),
    measure = "HR",
    hr = c(no_scm$hr, scm$hr, interaction$hr),
    lo = c(no_scm$conf.low, scm$conf.low, interaction$conf.low),
    hi = c(no_scm$conf.high, scm$conf.high, interaction$conf.high),
    label = c(
      fmt_hr(no_scm$hr, no_scm$conf.low, no_scm$conf.high),
      fmt_hr(scm$hr, scm$conf.low, scm$conf.high),
      fmt_hr(interaction$hr, interaction$conf.low, interaction$conf.high)
    ),
    events = sum(fit_object$data[[event_col]] == 1, na.rm = TRUE),
    n = nrow(fit_object$data),
    interaction_p = fit_object$lrt_p,
    stringsAsFactors = FALSE
  )
}

extract_direct_or <- function(fit, term) {
  if (!term %in% names(coef(fit))) stop(paste("Missing term:", term))
  data.frame(
    estimate = exp(coef(fit)[term]),
    conf.low = exp(fit$ci.lower[term]),
    conf.high = exp(fit$ci.upper[term]),
    p.value = fit$prob[term],
    row.names = NULL
  )
}

build_domestic_firth_rows <- function(data, outcome_label) {
  covars <- c("age", "female", "bmi", "sofa", "log_lactate", "log_creatinine")
  data$scm_fac <- factor(ifelse(data$scm == 1, "SCM", "No SCM"), levels = c("No SCM", "SCM"))
  rhs <- c("pni_low5_c * scm_fac", covars)
  f <- as.formula(paste("death_event ~", paste(rhs, collapse = " + ")))
  vars <- unique(all.vars(f))
  d <- data[complete.cases(data[, vars]), ]

  fit_no_scm_ref <- logistf(f, data = d, pl = TRUE)

  d$scm_fac <- factor(ifelse(d$scm == 1, "SCM", "No SCM"), levels = c("SCM", "No SCM"))
  fit_scm_ref <- logistf(f, data = d, pl = TRUE)

  int_name <- intersect(c("pni_low5_c:scm_facSCM", "scm_facSCM:pni_low5_c"), names(coef(fit_no_scm_ref)))[1]
  if (is.na(int_name)) stop("Domestic interaction term not found.")
  int_idx <- which(names(coef(fit_no_scm_ref)) == int_name)
  int_test <- logistftest(fit_no_scm_ref, test = int_idx)

  no_scm <- extract_direct_or(fit_no_scm_ref, "pni_low5_c")
  scm <- extract_direct_or(fit_scm_ref, "pni_low5_c")
  interaction <- extract_direct_or(fit_no_scm_ref, int_name)

  data.frame(
    outcome = outcome_label,
    estimand = c("No SCM", "SCM", "Interaction"),
    measure = "OR",
    hr = c(no_scm$estimate, scm$estimate, interaction$estimate),
    lo = c(no_scm$conf.low, scm$conf.low, interaction$conf.low),
    hi = c(no_scm$conf.high, scm$conf.high, interaction$conf.high),
    label = c(
      fmt_hr(no_scm$estimate, no_scm$conf.low, no_scm$conf.high),
      fmt_hr(scm$estimate, scm$conf.low, scm$conf.high),
      fmt_hr(interaction$estimate, interaction$conf.low, interaction$conf.high)
    ),
    events = sum(d$death_event == 1, na.rm = TRUE),
    n = nrow(d),
    interaction_p = as.numeric(int_test$prob),
    stringsAsFactors = FALSE
  )
}

load_domestic_model_data <- function(path) {
  raw <- read_excel(path)
  if (ncol(raw) < 50) stop("Domestic cohort file has fewer columns than expected.")
  sex_text <- trimws(as.character(raw[[5]]))
  status_text <- trimws(as.character(raw[[48]]))
  d <- data.frame(
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
  d$pni <- d$albumin + 5 * d$lymphocyte
  d$scm <- ifelse(d$lvef < 50, 1, 0)
  d$death_event <- ifelse(grepl("\u597d\u8f6c", status_text), 0, 1)
  d$pni_low5 <- -d$pni / 5
  d$pni_low5_c <- d$pni_low5 - mean(d$pni_low5, na.rm = TRUE)
  d$log_lactate <- ifelse(d$lactate > 0, log(d$lactate), NA_real_)
  d$log_creatinine <- ifelse(d$creatinine > 0, log(d$creatinine), NA_real_)
  d
}

dat <- read.csv(mimic_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
num_cols <- c(
  "age", "female", "bmi", "sofa", "scm", "pni", "lactate", "creatinine",
  "coronary_artery_disease", "diabetes_mellitus", "hypertension",
  "chronic_kidney_disease", "copd", "mechanical_vent", "vasopressor",
  "event_28d", "time_28d", "event_90d", "time_90d"
)
for (cc in intersect(num_cols, names(dat))) dat[[cc]] <- to_num(dat[[cc]])

dat$pni_low5 <- -dat$pni / 5
dat$pni_low5_c <- dat$pni_low5 - mean(dat$pni_low5, na.rm = TRUE)
dat$log_lactate <- ifelse(dat$lactate > 0, log(dat$lactate), NA_real_)
dat$log_creatinine <- ifelse(dat$creatinine > 0, log(dat$creatinine), NA_real_)

covars_mimic <- c(
  "age", "female", "bmi", "sofa", "log_lactate", "log_creatinine",
  "coronary_artery_disease", "diabetes_mellitus", "hypertension",
  "chronic_kidney_disease", "copd", "mechanical_vent", "vasopressor"
)
fit28 <- fit_model3(dat, "time_28d", "event_28d", covars_mimic)
fit90 <- fit_model3(dat, "time_90d", "event_90d", covars_mimic)
domestic_dat <- load_domestic_model_data(domestic_path)

plot_dat <- rbind(
  build_cox_rows(fit28, "MIMIC-IV: 28-day mortality", "event_28d"),
  build_cox_rows(fit90, "MIMIC-IV: 90-day mortality", "event_90d"),
  build_domestic_firth_rows(domestic_dat, "Domestic: exploratory in-hospital analysis")
)

plot_dat$row <- c(9.2, 8.2, 7.2, 5.5, 4.5, 3.5, 1.8, 0.8, -0.2)
plot_dat$group <- factor(plot_dat$estimand, levels = c("No SCM", "SCM", "Interaction"))
plot_dat$outcome_label <- paste0(plot_dat$outcome, " (", format(plot_dat$events, big.mark = ","), " events)")

fmt_interaction_p <- function(p) {
  if (is.na(p)) return("interaction P=NA")
  if (p < 0.001) return("interaction P<0.001")
  paste0("interaction P=", sprintf("%.3f", p))
}

header_dat <- data.frame(
  x = 0.36,
  y = c(9.9, 6.2, 2.5),
  label = c(
    paste0("MIMIC-IV: 72-h landmark 28-day mortality, n=", format(plot_dat$n[1], big.mark = ","), ", events=", plot_dat$events[1],
           ", ", fmt_interaction_p(plot_dat$interaction_p[1])),
    paste0("MIMIC-IV: 72-h landmark 90-day mortality, n=", format(plot_dat$n[4], big.mark = ","), ", events=", plot_dat$events[4],
           ", ", fmt_interaction_p(plot_dat$interaction_p[4])),
    paste0("Domestic exploratory: death endpoint, n=", format(plot_dat$n[7], big.mark = ","), ", events=", plot_dat$events[7],
           ", ", fmt_interaction_p(plot_dat$interaction_p[7]))
  )
)

axis_breaks <- c(0.5, 0.75, 1.0, 1.25, 1.5)

p <- ggplot(plot_dat, aes(x = hr, y = row, xmin = lo, xmax = hi, color = group, shape = group)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, color = "grey45") +
  geom_errorbar(orientation = "y", width = 0, linewidth = 0.85) +
  geom_point(size = 2.8) +
  geom_text(aes(x = 0.36, label = estimand), hjust = 0, color = "black", size = 3.35, family = "serif") +
  geom_text(aes(x = 1.55, label = label), hjust = 0, color = "black", size = 3.25, family = "serif") +
  geom_text(data = header_dat, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            hjust = 0, fontface = "bold", size = 3.45, family = "serif") +
  annotate("text", x = 0.36, y = 10.55, label = "Effect estimate", hjust = 0,
           fontface = "bold", size = 3.45, family = "serif") +
  annotate("text", x = 1.55, y = 10.55, label = "Estimate (95% CI)", hjust = 0,
           fontface = "bold", size = 3.45, family = "serif") +
  scale_x_log10(
    limits = c(0.35, 1.88),
    breaks = axis_breaks,
    labels = c("0.50", "0.75", "1.00", "1.25", "1.50")
  ) +
  scale_y_continuous(limits = c(-0.85, 10.85), breaks = NULL) +
  scale_color_manual(values = c("No SCM" = "#1F4E79", "SCM" = "#6E6E6E", "Interaction" = "#B23A2E")) +
  scale_shape_manual(values = c("No SCM" = 16, "SCM" = 15, "Interaction" = 18)) +
  labs(
    x = "Effect estimate per 5-point lower PNI (log scale)",
    y = NULL
  ) +
  theme_classic(base_family = "serif") +
  theme(
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.45),
    axis.ticks.x = element_line(color = "black", linewidth = 0.45),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, margin = margin(t = 7)),
    plot.margin = margin(8, 10, 8, 10)
  )

out_path <- file.path(model_result_dir, "\u53cc\u961f\u5217\u63a2\u7d22\u6027\u68ee\u6797\u56fe.tiff")
tiff(out_path, width = 7.6, height = 5.3, units = "in", res = 600, compression = "lzw")
print(p)
dev.off()

cat("Dual-cohort forest plot completed.\n")
cat("Source file:", mimic_path, "\n")
cat("Domestic file:", domestic_path, "\n")
cat("Saved TIFF:", out_path, "\n")
print(plot_dat[, c("outcome", "estimand", "hr", "lo", "hi", "label")], row.names = FALSE)


