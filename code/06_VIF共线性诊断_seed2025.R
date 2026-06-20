library(survival)

source("00_通用设置与函数_seed2025.R", encoding = "UTF-8")

dat <- load_mimic_model_data()

vif_from_fit <- function(fit) {
  mm <- model.matrix(fit)
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
  }
  vif <- sapply(seq_len(ncol(mm)), function(j) {
    y <- mm[, j]
    x <- mm[, -j, drop = FALSE]
    if (ncol(x) == 0 || stats::var(y) == 0) return(NA_real_)
    r2 <- tryCatch(summary(lm(y ~ x))$r.squared, error = function(e) NA_real_)
    if (is.na(r2) || r2 >= 1) return(Inf)
    1 / (1 - r2)
  })
  names(vif) <- colnames(mm)
  vif
}

label_vif_term <- function(x) {
  labels <- c(
    pni_low5_c = "PNI, per 5-point decrease, centered",
    scm = "SCM",
    age = "Age",
    female = "Female sex",
    bmi = "BMI",
    sofa = "SOFA score",
    log_lactate = "Log-transformed lactate",
    log_creatinine = "Log-transformed creatinine",
    coronary_artery_disease = "Coronary artery disease",
    diabetes_mellitus = "Diabetes mellitus",
    hypertension = "Hypertension",
    chronic_kidney_disease = "Chronic kidney disease",
    copd = "COPD",
    mechanical_vent = "Mechanical ventilation",
    vasopressor = "Vasopressor use",
    `pni_low5_c:scm` = "PNI, per 5-point decrease, centered x SCM",
    `scm:pni_low5_c` = "PNI, per 5-point decrease, centered x SCM"
  )
  out <- labels[x]
  out[is.na(out)] <- x[is.na(out)]
  unname(out)
}

rows <- list()
for (outcome_label in names(outcomes)) {
  fit_obj <- fit_one_cox(
    dat,
    outcomes[[outcome_label]]["time"],
    outcomes[[outcome_label]]["event"],
    main_covars[["Model 3: fully adjusted"]]
  )
  vif <- vif_from_fit(fit_obj$fit)
  rows[[length(rows) + 1]] <- data.frame(
    Outcome = outcome_label,
    Variable = label_vif_term(names(vif)),
    VIF = sprintf("%.2f", vif),
    row.names = NULL
  )
}

vif_table <- do.call(rbind, rows)
cat("VIF diagnosis: fully adjusted model after centering pni_low5.\n")
print_table(vif_table)


