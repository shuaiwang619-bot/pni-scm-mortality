# PNI-SCM mortality project
# Linyi domestic cohort: R mice imputation archive script
# Random seed is fixed at 2025 by project convention.

set.seed(2025)

required_pkgs <- c("mice", "readxl", "rpart")
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

input_file <- file.path(work_dir, "\u533b\u9662\u6570\u636e\u8113\u6bd2\u75c7\u961f\u5217\u539f\u59cb\u6587\u4ef6.xlsx")
sheet_name <- "Sheet1"

completed_csv <- file.path(work_dir, "linyi_domestic_mice_seed2025_completed.csv")
imputed_cells_csv <- file.path(work_dir, "linyi_domestic_mice_seed2025_imputed_cells.csv")
missingness_csv <- file.path(work_dir, "linyi_domestic_mice_seed2025_missingness.csv")
mids_file <- file.path(work_dir, "linyi_domestic_mice_seed2025_mids.rds")
logged_events_csv <- file.path(work_dir, "linyi_domestic_mice_seed2025_logged_events.csv")
log_file <- file.path(work_dir, "linyi_domestic_mice_seed2025_log.txt")

raw_tbl <- readxl::read_excel(input_file, sheet = sheet_name, col_types = "text", .name_repair = "unique")
raw <- as.data.frame(raw_tbl, stringsAsFactors = FALSE, check.names = FALSE)
names(raw) <- trimws(names(raw))
raw[] <- lapply(raw, function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA
  x
})

no_impute_vars <- c("\u5e8f\u53f7", "\u4f4f\u9662\u53f7", "\u59d3\u540d", "icu\u5929\u6570")

numeric_vars <- c(
  "\u5e74\u9f84", "bmi", "\u5fc3\u7387", "\u547c\u5438", "\u5e73\u5747\u52a8\u8109\u538b", "\u4f53\u6e29",
  "apache", "sofa", "wbc", "plt", "cr", "bun", "lac", "\u5c04\u8840\u5206\u6570",
  "\u808c\u9499\u86cb\u767dT", "ntprobnp", "\u6dcb\u5df4\u7ec6\u80de\u8ba1\u6570", "c\u53cd\u5e94\u86cb\u767d",
  "d\u4e8c\u805a\u4f53", "\u808c\u9178\u6fc0\u9176\u540c\u5de5\u9176", "\u964d\u9499\u7d20\u539f", "\u767d\u86cb\u767d",
  "\u603b\u80c6\u7ea2\u7d20", "\u7ed3\u5408\u80c6\u7ea2\u7d20", "\u6e38\u79bb\u80c6\u7ea2\u7d20",
  "\u8c37\u8349\u8f6c\u6c28\u9176", "\u8c37\u4e19\u8f6c\u6c28\u9176", "il_6", "\u6c27\u5408\u6307\u6570",
  "\u4f4f\u9662\u5929\u6570"
)
numeric_vars <- intersect(numeric_vars, names(raw))

model_vars <- setdiff(names(raw), no_impute_vars)
dat <- raw[, model_vars, drop = FALSE]

to_numeric <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", x)))
}

for (v in model_vars) {
  if (v %in% numeric_vars) {
    dat[[v]] <- to_numeric(dat[[v]])
  } else {
    dat[[v]] <- factor(dat[[v]])
  }
}

missing_before_all <- vapply(raw, function(x) sum(is.na(x)), integer(1))
target_vars <- names(missing_before_all)[missing_before_all > 0]
target_vars <- setdiff(target_vars, no_impute_vars)

if (length(target_vars) == 0) {
  stop("No eligible missing variables to impute.")
}

constant_vars <- names(dat)[vapply(dat, function(x) length(unique(x[!is.na(x)])) <= 1, logical(1))]

method <- mice::make.method(dat)
method[] <- ""
for (v in target_vars) {
  if (v %in% names(dat)) {
    if (v %in% numeric_vars) {
      method[v] <- "pmm"
    } else {
      method[v] <- "cart"
    }
  }
}
method[constant_vars] <- ""

predictor_matrix <- mice::make.predictorMatrix(dat)
predictor_matrix[,] <- 1
diag(predictor_matrix) <- 0
if (length(constant_vars) > 0) {
  predictor_matrix[, constant_vars] <- 0
  predictor_matrix[constant_vars, ] <- 0
}
if ("\u79d1\u5ba4" %in% colnames(predictor_matrix)) {
  predictor_matrix[, "\u79d1\u5ba4"] <- 0
}

imp <- mice::mice(
  dat,
  m = 10,
  maxit = 5,
  method = method,
  predictorMatrix = predictor_matrix,
  seed = 2025,
  printFlag = FALSE
)

completed_dat <- mice::complete(imp, action = 1)
completed_raw <- raw
imputed_cells <- data.frame(
  excel_row = integer(),
  variable = character(),
  imputed_value = character(),
  stringsAsFactors = FALSE
)

for (v in target_vars) {
  miss_idx <- which(is.na(raw[[v]]))
  if (length(miss_idx) == 0 || !v %in% names(completed_dat)) {
    next
  }
  values <- completed_dat[[v]][miss_idx]
  if (is.factor(values)) {
    values <- as.character(values)
  }
  if (v %in% numeric_vars) {
    values_chr <- as.character(signif(as.numeric(values), 10))
  } else {
    values_chr <- as.character(values)
  }
  completed_raw[[v]][miss_idx] <- values_chr
  imputed_cells <- rbind(
    imputed_cells,
    data.frame(
      excel_row = miss_idx + 1L,
      variable = v,
      imputed_value = values_chr,
      stringsAsFactors = FALSE
    )
  )
}

missing_after_all <- vapply(completed_raw, function(x) sum(is.na(x)), integer(1))
missingness <- data.frame(
  variable = names(raw),
  missing_before = as.integer(missing_before_all),
  missing_after = as.integer(missing_after_all),
  imputed_by_mice = names(raw) %in% target_vars,
  excluded_from_imputation = names(raw) %in% no_impute_vars,
  stringsAsFactors = FALSE
)
missingness$n <- nrow(raw)
missingness$missing_before_pct <- round(100 * missingness$missing_before / missingness$n, 2)
missingness$missing_after_pct <- round(100 * missingness$missing_after / missingness$n, 2)

write.csv(completed_raw, completed_csv, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(imputed_cells, imputed_cells_csv, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(missingness, missingness_csv, row.names = FALSE, fileEncoding = "UTF-8")
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
write.csv(logged_events, logged_events_csv, row.names = FALSE, fileEncoding = "UTF-8")

imputed_summary <- aggregate(
  excel_row ~ variable,
  data = imputed_cells,
  FUN = length
)
names(imputed_summary)[2] <- "imputed_cells"

log_lines <- c(
  "Linyi domestic cohort R mice imputation archive",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Working directory: ", normalizePath(work_dir, winslash = "/", mustWork = FALSE)),
  paste0("Input file: ", basename(input_file)),
  paste0("Sheet: ", sheet_name),
  "Seed: 2025",
  "Package: mice",
  "Imputation model: m=10, maxit=5; numeric variables use predictive mean matching; categorical variables use CART.",
  "Excluded from imputation: 序号, 住院号, 姓名, icu天数.",
  paste0("mice logged events count: ", nrow(logged_events)),
  "",
  "Imputed-cell summary:",
  paste(capture.output(print(imputed_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "Missingness before and after:",
  paste(capture.output(print(missingness[missingness$missing_before > 0 | missingness$missing_after > 0, ], row.names = FALSE)), collapse = "\n")
)
con <- file(log_file, open = "w", encoding = "UTF-8")
writeLines(enc2utf8(log_lines), con = con)
close(con)

cat("Completed Linyi R mice imputation.\n")
cat("Input:", input_file, "\n")
cat("Completed CSV:", completed_csv, "\n")
cat("Imputed cells:", imputed_cells_csv, "\n")
cat("Missingness:", missingness_csv, "\n")
cat("mids:", mids_file, "\n")
cat("Logged events:", logged_events_csv, "\n")
cat("Log:", log_file, "\n")

