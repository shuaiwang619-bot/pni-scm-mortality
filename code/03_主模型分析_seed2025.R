library(survival)

source("00_通用设置与函数_seed2025.R", encoding = "UTF-8")

dat <- load_mimic_model_data()

rows <- list()
for (outcome_label in names(outcomes)) {
  for (model_label in names(main_covars)) {
    fit_obj <- fit_one_cox(
      dat,
      outcomes[[outcome_label]]["time"],
      outcomes[[outcome_label]]["event"],
      main_covars[[model_label]]
    )
    rows[[length(rows) + 1]] <- summarize_cox_fit(
      fit_obj,
      outcome_label,
      model_label,
      outcomes[[outcome_label]]["event"]
    )
  }
}

main_hr_table <- do.call(rbind, rows)
cat("Main Cox model: HRs are per 5-point decrease in PNI; PNI was centered for interaction modeling.\n")
print_table(main_hr_table)


