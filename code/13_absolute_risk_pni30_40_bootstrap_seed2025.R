set.seed(2025)

code_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
common_file <- list.files(code_dir, pattern = "^00_.*seed2025[.]R$", full.names = TRUE)
if (length(common_file) != 1) {
  stop("Expected exactly one 00_*_seed2025.R file in the code directory.")
}

source(common_file, encoding = "UTF-8")

analysis_dir <- normalizePath(file.path(code_dir, ".."), winslash = "/", mustWork = TRUE)
model_result_dir <- file.path(analysis_dir, basename(model_result_dir))
if (!dir.exists(model_result_dir)) dir.create(model_result_dir, recursive = TRUE)

mimic_candidates <- list.files(analysis_dir, pattern = "^mimic.*[.]csv$", full.names = TRUE)
if (length(mimic_candidates) < 1) {
  stop("Could not find the MIMIC analysis CSV in the analysis directory.")
}
mimic_path <- mimic_candidates[1]

boot_B <- as.integer(Sys.getenv("PNI_ABS_BOOT_B", "1000"))
if (is.na(boot_B) || boot_B < 100) stop("PNI_ABS_BOOT_B must be an integer >= 100.")

pni_contrast <- c(low_risk = 40, high_risk = 30)
full_covars <- main_covars[["Model 3: fully adjusted"]]

dat <- load_mimic_model_data()
dat$pni_low5 <- -dat$pni / 5
pni_low5_center <- mean(dat$pni_low5, na.rm = TRUE)
dat$pni_low5_c <- dat$pni_low5 - pni_low5_center

fmt_num <- function(x, digits = 1) {
  ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
}

fmt_est_ci <- function(est, lo, hi, multiplier = 100, digits = 1) {
  if (any(is.na(c(est, lo, hi)))) return("NA")
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
    multiplier * est, multiplier * lo, multiplier * hi
  )
}

prepare_model_data <- function(data, event_col, covars) {
  vars <- unique(c(event_col, "pni_low5_c", "scm", covars))
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop(paste("Missing variable(s):", paste(missing_vars, collapse = ", ")))
  }
  d <- data[complete.cases(data[, vars]), vars, drop = FALSE]
  d <- d[!is.na(d[[event_col]]) & d[[event_col]] %in% c(0, 1), , drop = FALSE]
  d[[event_col]] <- as.integer(d[[event_col]] == 1)
  d
}

build_formula <- function(event_col, covars) {
  rhs <- paste(c("pni_low5_c * scm", covars), collapse = " + ")
  as.formula(paste(event_col, "~", rhs))
}

standardized_risk <- function(fit, data, scm_value, pni_value) {
  newdata <- data
  newdata$scm <- scm_value
  newdata$pni_low5_c <- (-pni_value / 5) - pni_low5_center
  mean(stats::predict(fit, newdata = newdata, type = "response"), na.rm = TRUE)
}

estimate_absolute_effects <- function(data, event_col, covars) {
  f <- build_formula(event_col, covars)
  fit <- stats::glm(
    f,
    data = data,
    family = stats::binomial(),
    control = stats::glm.control(maxit = 50)
  )

  risk_no_scm_pni40 <- standardized_risk(fit, data, 0, pni_contrast["low_risk"])
  risk_no_scm_pni30 <- standardized_risk(fit, data, 0, pni_contrast["high_risk"])
  risk_scm_pni40 <- standardized_risk(fit, data, 1, pni_contrast["low_risk"])
  risk_scm_pni30 <- standardized_risk(fit, data, 1, pni_contrast["high_risk"])

  rd_no_scm <- risk_no_scm_pni30 - risk_no_scm_pni40
  rd_scm <- risk_scm_pni30 - risk_scm_pni40

  c(
    risk_no_scm_pni40 = risk_no_scm_pni40,
    risk_no_scm_pni30 = risk_no_scm_pni30,
    rd_no_scm = rd_no_scm,
    risk_scm_pni40 = risk_scm_pni40,
    risk_scm_pni30 = risk_scm_pni30,
    rd_scm = rd_scm,
    rd_difference_scm_minus_no_scm = rd_scm - rd_no_scm
  )
}

bootstrap_effects <- function(data, event_col, covars, B) {
  point <- estimate_absolute_effects(data, event_col, covars)
  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(point))
  colnames(boot_mat) <- names(point)

  for (b in seq_len(B)) {
    idx <- sample.int(nrow(data), size = nrow(data), replace = TRUE)
    boot_data <- data[idx, , drop = FALSE]
    boot_mat[b, ] <- tryCatch(
      estimate_absolute_effects(boot_data, event_col, covars),
      error = function(e) rep(NA_real_, length(point))
    )
    if (b %% 100 == 0) {
      cat("Bootstrap", b, "of", B, "for", event_col, "\n")
    }
  }

  ci <- t(apply(
    boot_mat,
    2,
    function(x) stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  ))
  colnames(ci) <- c("ci_low", "ci_high")

  list(
    point = point,
    ci = ci,
    successful_bootstrap = sum(stats::complete.cases(boot_mat)),
    boot_mat = boot_mat
  )
}

make_outcome_rows <- function(outcome_label, event_col) {
  model_data <- prepare_model_data(dat, event_col, full_covars)
  effects <- bootstrap_effects(model_data, event_col, full_covars, boot_B)
  p <- effects$point
  ci <- effects$ci

  numeric_rows <- data.frame(
    Outcome = outcome_label,
    SCM = c("No SCM", "SCM"),
    N_model = nrow(model_data),
    Events = sum(model_data[[event_col]] == 1, na.rm = TRUE),
    Bootstrap_resamples = boot_B,
    Successful_bootstrap = effects$successful_bootstrap,
    Risk_PNI40 = c(p["risk_no_scm_pni40"], p["risk_scm_pni40"]),
    Risk_PNI40_CI_low = c(ci["risk_no_scm_pni40", "ci_low"], ci["risk_scm_pni40", "ci_low"]),
    Risk_PNI40_CI_high = c(ci["risk_no_scm_pni40", "ci_high"], ci["risk_scm_pni40", "ci_high"]),
    Risk_PNI30 = c(p["risk_no_scm_pni30"], p["risk_scm_pni30"]),
    Risk_PNI30_CI_low = c(ci["risk_no_scm_pni30", "ci_low"], ci["risk_scm_pni30", "ci_low"]),
    Risk_PNI30_CI_high = c(ci["risk_no_scm_pni30", "ci_high"], ci["risk_scm_pni30", "ci_high"]),
    Risk_difference_PNI30_minus_PNI40 = c(p["rd_no_scm"], p["rd_scm"]),
    Risk_difference_CI_low = c(ci["rd_no_scm", "ci_low"], ci["rd_scm", "ci_low"]),
    Risk_difference_CI_high = c(ci["rd_no_scm", "ci_high"], ci["rd_scm", "ci_high"]),
    row.names = NULL,
    check.names = FALSE
  )

  formatted_rows <- data.frame(
    Outcome = outcome_label,
    SCM = c("No SCM", "SCM"),
    N = nrow(model_data),
    Events = sum(model_data[[event_col]] == 1, na.rm = TRUE),
    `Adjusted risk at PNI 40, % (95% CI)` = c(
      fmt_est_ci(p["risk_no_scm_pni40"], ci["risk_no_scm_pni40", "ci_low"], ci["risk_no_scm_pni40", "ci_high"]),
      fmt_est_ci(p["risk_scm_pni40"], ci["risk_scm_pni40", "ci_low"], ci["risk_scm_pni40", "ci_high"])
    ),
    `Adjusted risk at PNI 30, % (95% CI)` = c(
      fmt_est_ci(p["risk_no_scm_pni30"], ci["risk_no_scm_pni30", "ci_low"], ci["risk_no_scm_pni30", "ci_high"]),
      fmt_est_ci(p["risk_scm_pni30"], ci["risk_scm_pni30", "ci_low"], ci["risk_scm_pni30", "ci_high"])
    ),
    `Risk difference, percentage points (95% CI)` = c(
      fmt_est_ci(p["rd_no_scm"], ci["rd_no_scm", "ci_low"], ci["rd_no_scm", "ci_high"]),
      fmt_est_ci(p["rd_scm"], ci["rd_scm", "ci_low"], ci["rd_scm", "ci_high"])
    ),
    check.names = FALSE
  )

  interaction_row <- data.frame(
    Outcome = outcome_label,
    Contrast = "Risk difference in SCM minus risk difference in No SCM",
    Estimate = p["rd_difference_scm_minus_no_scm"],
    CI_low = ci["rd_difference_scm_minus_no_scm", "ci_low"],
    CI_high = ci["rd_difference_scm_minus_no_scm", "ci_high"],
    Formatted = fmt_est_ci(
      p["rd_difference_scm_minus_no_scm"],
      ci["rd_difference_scm_minus_no_scm", "ci_low"],
      ci["rd_difference_scm_minus_no_scm", "ci_high"]
    ),
    row.names = NULL,
    check.names = FALSE
  )

  list(
    numeric_rows = numeric_rows,
    formatted_rows = formatted_rows,
    interaction_row = interaction_row
  )
}

pni_distribution <- do.call(rbind, lapply(c(0, 1), function(g) {
  x <- dat$pni[dat$scm == g]
  qs <- stats::quantile(x, probs = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99), na.rm = TRUE)
  data.frame(
    SCM = ifelse(g == 1, "SCM", "No SCM"),
    N = sum(!is.na(x)),
    PNI_min = min(x, na.rm = TRUE),
    PNI_p01 = qs[1],
    PNI_p05 = qs[2],
    PNI_p25 = qs[3],
    PNI_median = qs[4],
    PNI_p75 = qs[5],
    PNI_p95 = qs[6],
    PNI_p99 = qs[7],
    PNI_max = max(x, na.rm = TRUE),
    PNI30_within_p05_p95 = 30 >= qs[2] & 30 <= qs[6],
    PNI40_within_p05_p95 = 40 >= qs[2] & 40 <= qs[6],
    row.names = NULL,
    check.names = FALSE
  )
}))

outcome_results <- lapply(names(outcomes), function(label) {
  make_outcome_rows(label, outcomes[[label]]["event"])
})

numeric_table <- do.call(rbind, lapply(outcome_results, `[[`, "numeric_rows"))
formatted_table <- do.call(rbind, lapply(outcome_results, `[[`, "formatted_rows"))
interaction_table <- do.call(rbind, lapply(outcome_results, `[[`, "interaction_row"))

numeric_path <- file.path(model_result_dir, "absolute_risk_PNI30_vs_40_bootstrap_seed2025_numeric.csv")
formatted_path <- file.path(model_result_dir, "absolute_risk_PNI30_vs_40_bootstrap_seed2025_formatted.csv")
interaction_path <- file.path(model_result_dir, "absolute_risk_PNI30_vs_40_bootstrap_seed2025_interaction.csv")
pni_dist_path <- file.path(model_result_dir, "absolute_risk_PNI_distribution_seed2025.csv")

utils::write.csv(numeric_table, numeric_path, row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(formatted_table, formatted_path, row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(interaction_table, interaction_path, row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(pni_distribution, pni_dist_path, row.names = FALSE, fileEncoding = "UTF-8")

if (requireNamespace("officer", quietly = TRUE) && requireNamespace("flextable", quietly = TRUE)) {
  docx_path <- file.path(model_result_dir, "absolute_risk_PNI30_vs_40_bootstrap_seed2025.docx")
  doc <- officer::read_docx()
  doc <- officer::body_add_par(
    doc,
    "Secondary absolute-risk translation: PNI 30 vs PNI 40",
    style = "heading 1"
  )
  doc <- officer::body_add_par(
    doc,
    paste0(
      "Fixed-time logistic models were fitted for 28-day and 90-day mortality. ",
      "Risks were standardized over the analytic covariate distribution after setting SCM status ",
      "and PNI values. Bootstrap resamples: ", boot_B, "."
    ),
    style = "Normal"
  )
  ft <- flextable::flextable(formatted_table)
  ft <- flextable::autofit(ft)
  doc <- flextable::body_add_flextable(doc, ft)
  doc <- officer::body_add_par(doc, "Absolute-scale interaction contrast", style = "heading 2")
  ft2 <- flextable::flextable(interaction_table)
  ft2 <- flextable::autofit(ft2)
  doc <- flextable::body_add_flextable(doc, ft2)
  print(doc, target = docx_path)
  cat("DOCX written:", docx_path, "\n")
} else {
  cat("officer/flextable not available; skipped DOCX output.\n")
}

cat("Numeric CSV written:", numeric_path, "\n")
cat("Formatted CSV written:", formatted_path, "\n")
cat("Interaction CSV written:", interaction_path, "\n")
cat("PNI distribution CSV written:", pni_dist_path, "\n")
cat("PNI lower-5 center:", sprintf("%.8f", pni_low5_center), "\n")
cat("Analysis complete.\n")


