set.seed(2025)

library(survival)
library(splines)
library(ggplot2)
library(patchwork)

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) return(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()
project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
analysis_dir <- project_dir
model_result_dir <- file.path(project_dir, "results")
mimic_path <- file.path(project_dir, "data_private", "mimic队列插补后.csv")

dir.create(model_result_dir, showWarnings = FALSE, recursive = TRUE)

to_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9eE+\\.-]", "", trimws(as.character(x)))))

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
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

covars <- c(
  "age", "female", "bmi", "sofa", "log_lactate", "log_creatinine",
  "coronary_artery_disease", "diabetes_mellitus", "hypertension",
  "chronic_kidney_disease", "copd", "mechanical_vent", "vasopressor"
)

outcomes <- list(
  "28-day mortality" = c(time = "time_28d", event = "event_28d"),
  "90-day mortality" = c(time = "time_90d", event = "event_90d")
)

typical_covariates <- function(d) {
  out <- data.frame(row.names = 1)
  for (cc in covars) {
    x <- d[[cc]]
    if (all(x %in% c(0, 1), na.rm = TRUE)) {
      out[[cc]] <- as.numeric(mean(x, na.rm = TRUE) >= 0.5)
    } else {
      out[[cc]] <- median(x, na.rm = TRUE)
    }
  }
  out
}

fit_rcs_outcome <- function(data, time_col, event_col, outcome_label) {
  vars <- c(time_col, event_col, "pni", "scm", covars)
  d <- data[complete.cases(data[, vars]), ]
  d <- d[d[[time_col]] > 0 & !is.na(d[[event_col]]), ]

  linear_formula <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ pni * scm + ",
    paste(covars, collapse = " + ")
  ))
  spline_formula <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ns(pni, df = 3) * scm + ",
    paste(covars, collapse = " + ")
  ))
  spline_main_formula <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ns(pni, df = 3) + scm + ",
    paste(covars, collapse = " + ")
  ))

  fit_linear <- coxph(linear_formula, data = d, ties = "efron", x = TRUE)
  fit_spline_main <- coxph(spline_main_formula, data = d, ties = "efron", x = TRUE)
  fit_spline <- coxph(spline_formula, data = d, ties = "efron", x = TRUE)
  lrt <- anova(fit_linear, fit_spline, test = "LRT")
  p_nonlinear <- lrt$`Pr(>|Chi|)`[2]
  interaction_lrt <- anova(fit_spline_main, fit_spline, test = "LRT")
  p_interaction <- interaction_lrt$`Pr(>|Chi|)`[2]

  xlim <- as.numeric(quantile(d$pni, probs = c(0.01, 0.99), na.rm = TRUE, type = 2))
  grid <- seq(xlim[1], xlim[2], length.out = 160)
  cov_base <- typical_covariates(d)

  pred_rows <- list()
  for (scm_value in c(0, 1)) {
    nd <- data.frame(pni = grid, scm = scm_value)
    for (cc in covars) nd[[cc]] <- cov_base[[cc]][1]

    ref <- data.frame(pni = 35, scm = scm_value)
    for (cc in covars) ref[[cc]] <- cov_base[[cc]][1]

    mm <- model.matrix(delete.response(terms(fit_spline)), data = nd)
    mm_ref <- model.matrix(delete.response(terms(fit_spline)), data = ref)
    keep <- names(coef(fit_spline))
    mm <- mm[, keep, drop = FALSE]
    mm_ref <- mm_ref[, keep, drop = FALSE]

    contrast <- sweep(mm, 2, mm_ref[1, ], "-")
    beta <- coef(fit_spline)
    vc <- vcov(fit_spline)
    lp <- as.numeric(contrast %*% beta)
    se <- sqrt(rowSums((contrast %*% vc) * contrast))

    pred_rows[[length(pred_rows) + 1]] <- data.frame(
      Outcome = outcome_label,
      pni = grid,
      scm = ifelse(scm_value == 1, "SCM", "No SCM"),
      HR = exp(lp),
      lower = exp(lp - 1.96 * se),
      upper = exp(lp + 1.96 * se),
      p_nonlinear = p_nonlinear,
      p_interaction = p_interaction,
      N = nrow(d),
      Events = sum(d[[event_col]] == 1, na.rm = TRUE)
    )
  }

  list(
    prediction = do.call(rbind, pred_rows),
    rug = data.frame(
      Outcome = outcome_label,
      pni = d$pni,
      scm = ifelse(d$scm == 1, "SCM", "No SCM"),
      stringsAsFactors = FALSE
    ),
    diagnostics = data.frame(
      Outcome = outcome_label,
      N = nrow(d),
      Events = sum(d[[event_col]] == 1, na.rm = TRUE),
      Reference_PNI = 35,
      Plot_PNI_min_1pct = xlim[1],
      Plot_PNI_max_99pct = xlim[2],
      LRT_P_for_nonlinearity = p_nonlinear,
      LRT_P_for_spline_interaction = p_interaction,
      stringsAsFactors = FALSE
    )
  )
}

results <- list()
for (label in names(outcomes)) {
  results[[label]] <- fit_rcs_outcome(
    dat,
    outcomes[[label]]["time"],
    outcomes[[label]]["event"],
    label
  )
}

pred <- do.call(rbind, lapply(results, `[[`, "prediction"))
rug_data <- do.call(rbind, lapply(results, `[[`, "rug"))
diag_table <- do.call(rbind, lapply(results, `[[`, "diagnostics"))

make_panel <- function(data, rug_data, title_label) {
  sub <- data[data$Outcome == title_label, ]
  sub_rug <- rug_data[
    rug_data$Outcome == title_label &
      rug_data$pni >= min(sub$pni, na.rm = TRUE) &
      rug_data$pni <= max(sub$pni, na.rm = TRUE),
  ]
  ptxt <- paste0(
    "N=", unique(sub$N),
    ", events=", unique(sub$Events),
    "; nonlinearity P=", fmt_p(unique(sub$p_nonlinear)),
    "; spline interaction P=", fmt_p(unique(sub$p_interaction))
  )
  ggplot(sub, aes(x = pni, y = HR, color = scm, fill = scm)) +
    geom_hline(yintercept = 1, linewidth = 0.4, linetype = "dashed", color = "gray45") +
    geom_vline(xintercept = 35, linewidth = 0.4, linetype = "dotted", color = "gray45") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.14, linewidth = 0, color = NA) +
    geom_line(linewidth = 1.15) +
    geom_rug(
      data = sub_rug[sub_rug$scm == "No SCM", ],
      aes(x = pni, color = scm),
      inherit.aes = FALSE,
      sides = "b",
      alpha = 0.12,
      linewidth = 0.25,
      length = grid::unit(0.025, "npc")
    ) +
    geom_rug(
      data = sub_rug[sub_rug$scm == "SCM", ],
      aes(x = pni, color = scm),
      inherit.aes = FALSE,
      sides = "b",
      alpha = 0.50,
      linewidth = 0.35,
      length = grid::unit(0.045, "npc")
    ) +
    scale_color_manual(values = c("No SCM" = "#0B5D66", "SCM" = "#B2472D")) +
    scale_fill_manual(values = c("No SCM" = "#0B5D66", "SCM" = "#B2472D")) +
    scale_y_log10() +
    labs(
      title = title_label,
      subtitle = ptxt,
      x = "PNI",
      y = "Adjusted HR within strata (reference PNI = 35)",
      color = NULL,
      fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

p28 <- make_panel(pred, rug_data, "28-day mortality")
p90 <- make_panel(pred, rug_data, "90-day mortality")

fig <- (p28 + p90) +
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = paste(
      "Rug marks show the PNI distribution by SCM status. Shaded areas indicate 95% confidence intervals.",
      "HRs are estimated within each SCM stratum using PNI = 35 as the reference.",
      "Estimates at PNI extremes, particularly high PNI among SCM patients, should be interpreted cautiously.",
      sep = "\n"
    )
  ) &
  theme(
    legend.position = "bottom",
    plot.caption = element_text(size = 8, hjust = 0)
  )

tiff_path <- file.path(model_result_dir, "Figure_PNI_RCS_28_90day_seed2025.tiff")

ggsave(tiff_path, fig, width = 10, height = 5.8, dpi = 600, compression = "lzw")

cat("PNI RCS analysis completed.\n")
print(diag_table)
cat("Saved TIFF:", tiff_path, "\n")



