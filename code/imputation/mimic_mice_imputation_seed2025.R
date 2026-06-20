# PNI-SCM mortality project
# MIMIC-IV 72-hour landmark cohort: R mice imputation archive script
# Random seed is fixed at 2025 by project convention.

set.seed(2025)

required_pkgs <- c("mice")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "))
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
work_dir <- script_dir

csv_files <- list.files(work_dir, pattern = "\\.csv$", full.names = TRUE)
raw_mimic_name <- "\u754c\u6807\u5206\u6790mimic\u961f\u5217.csv"
input_file <- file.path(work_dir, raw_mimic_name)
if (!file.exists(input_file)) {
  input_file <- csv_files[
    grepl("mimic", tolower(basename(csv_files))) &
      !grepl("MICE|missingness|flow|seed2025|mids|log|R_mice", basename(csv_files), ignore.case = TRUE)
  ]
}
if (length(input_file) != 1) {
  stop("Could not uniquely locate the raw MIMIC cohort CSV.")
}

output_file <- file.path(work_dir, "mimic_landmark72_R_mice_seed2025.csv")
missingness_file <- file.path(work_dir, "mimic_landmark72_R_mice_seed2025_missingness.csv")
mids_file <- file.path(work_dir, "mimic_landmark72_R_mice_seed2025_mids.rds")
logged_events_file <- file.path(work_dir, "mimic_landmark72_R_mice_seed2025_logged_events.csv")
log_file <- file.path(work_dir, "mimic_landmark72_R_mice_seed2025_log.txt")

df <- read.csv(input_file, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE, check.names = FALSE)

target_vars <- c("bmi", "lactate", "wbc")
required_cols <- c(
  "subject_id", "hadm_id", "stay_id",
  target_vars,
  "age", "female", "sofa", "apsiii", "scm",
  "lvef_selected", "albumin_g_l", "lymphocyte_k_ul", "pni",
  "creatinine", "platelets", "bun_mmol_l",
  "mechanical_vent", "vasopressor",
  "coronary_artery_disease", "diabetes_mellitus", "hypertension",
  "chronic_kidney_disease", "copd",
  "event_28d", "time_28d", "event_90d", "time_90d", "hosp_event"
)
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Input file is missing required columns: ", paste(missing_cols, collapse = ", "))
}

range_rules <- data.frame(
  variable = c("bmi", "lactate", "wbc"),
  lower = c(10, 0.1, 0.1),
  upper = c(80, 30, 200),
  stringsAsFactors = FALSE
)
extreme_cleaning <- data.frame(
  variable = character(),
  lower = numeric(),
  upper = numeric(),
  extreme_set_missing = integer(),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(range_rules))) {
  v <- range_rules$variable[i]
  lower <- range_rules$lower[i]
  upper <- range_rules$upper[i]
  values <- suppressWarnings(as.numeric(df[[v]]))
  extreme <- !is.na(values) & (values < lower | values > upper)
  df[[v]] <- values
  df[[v]][extreme] <- NA
  extreme_cleaning <- rbind(
    extreme_cleaning,
    data.frame(
      variable = v,
      lower = lower,
      upper = upper,
      extreme_set_missing = sum(extreme),
      stringsAsFactors = FALSE
    )
  )
}

before_missing <- data.frame(
  variable = target_vars,
  n = nrow(df),
  missing_before = vapply(target_vars, function(x) sum(is.na(df[[x]])), integer(1)),
  stringsAsFactors = FALSE
)
before_missing$missing_before_pct <- round(100 * before_missing$missing_before / before_missing$n, 2)

imp_vars <- required_cols
imp_dat <- df[, imp_vars]

# Coerce variables used by mice to numeric where clinically binary or continuous.
for (v in imp_vars) {
  imp_dat[[v]] <- suppressWarnings(as.numeric(imp_dat[[v]]))
}

method <- mice::make.method(imp_dat)
method[] <- ""
method["bmi"] <- "pmm"
method["lactate"] <- "pmm"
method["wbc"] <- "pmm"

predictor_matrix <- mice::make.predictorMatrix(imp_dat)
predictor_matrix[,] <- 0
predictors <- setdiff(imp_vars, c("subject_id", "hadm_id", "stay_id", "pni"))
for (target in target_vars) {
  predictor_matrix[target, setdiff(predictors, target)] <- 1
}

imp <- mice::mice(
  imp_dat,
  m = 10,
  maxit = 5,
  method = method,
  predictorMatrix = predictor_matrix,
  seed = 2025,
  printFlag = FALSE
)

completed <- mice::complete(imp, action = 1)
for (target in target_vars) {
  df[[target]] <- completed[[target]]
}

# Recompute transforms after imputation. Exposure components are not imputed here.
df$log_lactate <- ifelse(!is.na(df$lactate) & df$lactate > 0, log(df$lactate), NA)
df$log_creatinine <- ifelse(!is.na(df$creatinine) & df$creatinine > 0, log(df$creatinine), NA)
df$pni_low5 <- ifelse(!is.na(df$pni), -df$pni / 5, NA)
df$albumin_low5 <- ifelse(!is.na(df$albumin_g_l), -df$albumin_g_l / 5, NA)
df$lymph_low05 <- ifelse(!is.na(df$lymphocyte_k_ul), -df$lymphocyte_k_ul / 0.5, NA)

after_missing <- data.frame(
  variable = target_vars,
  missing_after = vapply(target_vars, function(x) sum(is.na(df[[x]])), integer(1)),
  stringsAsFactors = FALSE
)
missingness <- merge(before_missing, after_missing, by = "variable", all = TRUE)
missingness$missing_after_pct <- round(100 * missingness$missing_after / missingness$n, 2)

write.csv(df, output_file, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(missingness, missingness_file, row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(imp, mids_file)

logged_events <- imp$loggedEvents
if (is.null(logged_events)) {
  logged_events <- data.frame(
    it = integer(),
    im = integer(),
    dep = character(),
    meth = character(),
    out = character()
  )
}
write.csv(logged_events, logged_events_file, row.names = FALSE, fileEncoding = "UTF-8")

log_lines <- c(
  "MIMIC-IV 72-hour landmark R mice imputation archive",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Working directory: ", normalizePath(work_dir, winslash = "/", mustWork = FALSE)),
  paste0("Input file: ", basename(input_file)),
  paste0("Output file: ", basename(output_file)),
  paste0("mids RDS: ", basename(mids_file)),
  paste0("mice logged events file: ", basename(logged_events_file)),
  paste0("mice logged events count: ", nrow(logged_events)),
  "Seed: 2025",
  "Package: mice",
  "Imputation model: m=10, maxit=5, predictive mean matching for bmi, lactate, and wbc",
  "Predictor note: PNI was not used as an imputation predictor because it is derived from albumin and lymphocyte.",
  "Extreme-value cleaning before imputation: values outside prespecified physiologic ranges were set to missing.",
  "Not imputed: IDs, outcomes, follow-up times, SCM, PNI, albumin, lymphocyte, and complete covariates",
  "",
  "Extreme-value cleaning counts:",
  paste(capture.output(print(extreme_cleaning, row.names = FALSE)), collapse = "\n"),
  "",
  "Missingness before and after:",
  paste(capture.output(print(missingness, row.names = FALSE)), collapse = "\n")
)
con <- file(log_file, open = "w", encoding = "UTF-8")
writeLines(enc2utf8(log_lines), con = con)
close(con)

cat("Completed R mice imputation.\n")
cat("Input:", input_file, "\n")
cat("Output:", output_file, "\n")
cat("Missingness:", missingness_file, "\n")
cat("mids:", mids_file, "\n")
cat("Logged events:", logged_events_file, "\n")
cat("Log:", log_file, "\n")

