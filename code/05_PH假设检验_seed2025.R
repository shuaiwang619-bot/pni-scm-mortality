library(survival)

source("00_通用设置与函数_seed2025.R", encoding = "UTF-8")

dat <- load_mimic_model_data()

rows <- list()
for (outcome_label in names(outcomes)) {
  fit_obj <- fit_one_cox(
    dat,
    outcomes[[outcome_label]]["time"],
    outcomes[[outcome_label]]["event"],
    main_covars[["Model 3: fully adjusted"]]
  )
  z <- cox.zph(fit_obj$fit, transform = "km")
  tab <- as.data.frame(z$table)
  rows[[length(rows) + 1]] <- data.frame(
    Outcome = outcome_label,
    Term = rownames(tab),
    Chisq = sprintf("%.3f", tab$chisq),
    DF = tab$df,
    P_value = vapply(tab$p, fmt_p, character(1)),
    row.names = NULL
  )
}

ph_table <- do.call(rbind, rows)
cat("PH assumption test: Schoenfeld residual test using cox.zph.\n")
print_table(ph_table)


