library(survival)

source("00_通用设置与函数_seed2025.R", encoding = "UTF-8")

dat <- load_mimic_model_data()

covariate_rows <- list()
for (outcome_label in names(outcomes)) {
  for (model_label in names(sensitivity_covars)) {
    fit_obj <- fit_one_cox(
      dat,
      outcomes[[outcome_label]]["time"],
      outcomes[[outcome_label]]["event"],
      sensitivity_covars[[model_label]]
    )
    covariate_rows[[length(covariate_rows) + 1]] <- summarize_cox_fit(
      fit_obj,
      outcome_label,
      model_label,
      outcomes[[outcome_label]]["event"]
    )
  }
}

covariate_sensitivity_table <- do.call(rbind, covariate_rows)

time_rows <- list()
time_rows[[1]] <- summarize_cox_fit(
  fit_one_cox(dat, "time_28d", "event_28d", main_covars[["Model 3: fully adjusted"]]),
  "0-28 days after landmark",
  "Model 3: fully adjusted",
  "event_28d"
)

dat_28day_survivors <- dat[dat$time_90d > 28, ]
time_rows[[2]] <- summarize_cox_fit(
  fit_one_cox(dat_28day_survivors, "time_28_90", "event_28_90", main_covars[["Model 3: fully adjusted"]]),
  "28-90 days after landmark",
  "Model 3: among 28-day survivors",
  "event_28_90"
)

time_sensitivity_table <- do.call(rbind, time_rows)

cat("A. Covariate-adjustment sensitivity analysis.\n")
print_table(covariate_sensitivity_table)
cat("\nB. Time-window sensitivity analysis: early and late post-landmark intervals.\n")
print_table(time_sensitivity_table)


