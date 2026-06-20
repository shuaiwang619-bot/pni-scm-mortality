set.seed(2025)

library(survival)
library(ggplot2)
library(grid)

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
analysis_dir <- dirname(script_dir)
model_root <- file.path(analysis_dir, "\u6a21\u578b\u7ed3\u679c")
result_dir <- file.path(model_root, "\u65b0\u5efa\u6587\u4ef6\u5939")
if (!dir.exists(result_dir)) result_dir <- model_root
out_tiff <- file.path(result_dir, "Figure 4-\u6a21\u578b\u5171\u7ebf\u6027\u8bca\u65ad.tiff")

mimic_candidates <- list.files(analysis_dir, pattern = "^mimic.*\\.csv$", full.names = TRUE)
if (length(mimic_candidates) == 0) stop("MIMIC CSV not found in analysis directory.")
mimic_path <- mimic_candidates[1]

to_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9eE+\\.-]", "", trimws(as.character(x)))))

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

covars <- c(
  "age", "female", "bmi", "sofa", "log_lactate", "log_creatinine",
  "coronary_artery_disease", "diabetes_mellitus", "hypertension",
  "chronic_kidney_disease", "copd", "mechanical_vent", "vasopressor"
)

fit_model3 <- function(data, time_col, event_col) {
  f <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ ",
    paste(c("pni_low5_c * scm", covars), collapse = " + ")
  ))
  vars <- unique(all.vars(f))
  d <- data[complete.cases(data[, vars]), ]
  d <- d[d[[time_col]] > 0 & !is.na(d[[event_col]]), ]
  fit <- coxph(f, data = d, ties = "efron", x = TRUE)
  list(fit = fit, data = d)
}

vif_from_fit <- function(fit) {
  mm <- model.matrix(fit)
  mm <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
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

label_map <- c(
  pni_low5_c = "PNI-low5, centered",
  scm = "SCM",
  age = "Age",
  female = "Sex",
  bmi = "BMI",
  sofa = "SOFA score",
  log_lactate = "Log lactate",
  log_creatinine = "Log creatinine",
  coronary_artery_disease = "Coronary artery disease",
  diabetes_mellitus = "Diabetes mellitus",
  hypertension = "Hypertension",
  chronic_kidney_disease = "Chronic kidney disease",
  copd = "COPD",
  mechanical_vent = "Mechanical ventilation",
  vasopressor = "Vasopressor use",
  `pni_low5_c:scm` = "PNI-low5 x SCM",
  `scm:pni_low5_c` = "PNI-low5 x SCM"
)

short_label_map <- c(
  pni_low5_c = "PNI-low5",
  scm = "SCM",
  age = "Age",
  female = "Sex",
  bmi = "BMI",
  sofa = "SOFA",
  log_lactate = "log Lac",
  log_creatinine = "log Cr",
  coronary_artery_disease = "CAD",
  diabetes_mellitus = "DM",
  hypertension = "HTN",
  chronic_kidney_disease = "CKD",
  copd = "COPD",
  mechanical_vent = "MV",
  vasopressor = "Vaso",
  `pni_low5_c:scm` = "PNI x SCM",
  `scm:pni_low5_c` = "PNI x SCM"
)

label_terms <- function(x, short = FALSE) {
  mp <- if (short) short_label_map else label_map
  out <- mp[x]
  out[is.na(out)] <- x[is.na(out)]
  unname(out)
}

fit28 <- fit_model3(dat, "time_28d", "event_28d")
vif_values <- vif_from_fit(fit28$fit)

vif_df <- data.frame(
  term = names(vif_values),
  variable = label_terms(names(vif_values)),
  vif = as.numeric(vif_values),
  stringsAsFactors = FALSE
)
vif_df <- vif_df[order(vif_df$vif, decreasing = TRUE), ]
vif_df$variable <- factor(vif_df$variable, levels = rev(vif_df$variable))

mm <- model.matrix(fit28$fit)
mm <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
colnames(mm) <- label_terms(colnames(mm), short = TRUE)
cor_mat <- suppressWarnings(cor(mm, method = "spearman", use = "pairwise.complete.obs"))
cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
names(cor_df) <- c("var1", "var2", "rho")
term_levels <- colnames(mm)
cor_df$var1 <- factor(cor_df$var1, levels = term_levels)
cor_df$var2 <- factor(cor_df$var2, levels = rev(term_levels))

base_theme <- theme_minimal(base_family = "sans") +
  theme(
    text = element_text(color = "#1F2630"),
    plot.title = element_text(face = "bold", size = 11, margin = margin(b = 8)),
    axis.title = element_text(face = "bold", size = 8.5),
    axis.text = element_text(size = 7.2, color = "#1F2630"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

vif_plot <- ggplot(vif_df, aes(x = vif, y = variable)) +
  geom_vline(xintercept = 5, linetype = "dashed", linewidth = 0.35, color = "#9C2F2F") +
  geom_col(width = 0.66, fill = "#4169E1") +
  geom_text(aes(label = sprintf("%.2f", vif)), hjust = -0.12, size = 2.35, color = "#1F2630") +
  annotate("text", x = 5, y = 0.65, label = "VIF = 5 reference", angle = 90,
           hjust = 0, vjust = -0.35, size = 2.25, color = "#9C2F2F") +
  scale_x_continuous(limits = c(0, 5.45), breaks = 0:5, expand = c(0, 0)) +
  labs(title = "A. Variance inflation factors", x = "VIF", y = NULL) +
  base_theme +
  theme(
    panel.grid.major.x = element_line(color = "#E0E6EF", linewidth = 0.35),
    plot.margin = margin(5, 10, 5, 5)
  )

heat_plot <- ggplot(cor_df, aes(x = var1, y = var2, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(abs(rho) >= 0.40 & var1 != var2, sprintf("%.2f", rho), "")),
            size = 1.85, color = "#1F2630") +
  scale_fill_gradient2(
    low = "#2F6DB3", mid = "white", high = "#B23A48",
    midpoint = 0, limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Spearman r"
  ) +
  coord_fixed() +
  labs(title = "B. Pairwise correlations among model terms", x = NULL, y = NULL) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.3),
    axis.text.y = element_text(size = 6.3),
    legend.position = "right",
    legend.title = element_text(size = 7.0, face = "bold"),
    legend.text = element_text(size = 6.4),
    plot.margin = margin(5, 5, 5, 8)
  )

tiff(out_tiff, width = 11.2, height = 6.6, units = "in", res = 600, compression = "lzw")
grid.newpage()
grid.rect(gp = gpar(fill = "#EAF0F7", col = NA))

grid.roundrect(
  x = 0.035, y = 0.085, width = 0.930, height = 0.875,
  just = c("left", "bottom"), r = unit(0.012, "snpc"),
  gp = gpar(fill = "white", col = "#D7DEE8", lwd = 0.7)
)

pushViewport(viewport(x = 0.055, y = 0.130, width = 0.385, height = 0.805, just = c("left", "bottom")))
grid.draw(ggplotGrob(vif_plot))
popViewport()

pushViewport(viewport(x = 0.455, y = 0.130, width = 0.495, height = 0.805, just = c("left", "bottom")))
grid.draw(ggplotGrob(heat_plot))
popViewport()

grid.text(
  "Main model terms only; LVEF, albumin and lymphocyte count were not included in the mortality model.",
  x = 0.5, y = 0.058, just = c(0.5, 0.5),
  gp = gpar(fontsize = 7.5, fontfamily = "sans", col = "#5D6876")
)
grid.text(
  "VIFs were calculated from the fully adjusted Cox model matrix after centering PNI-low5; the same covariate structure was used for 28-day and 90-day mortality.",
  x = 0.5, y = 0.033, just = c(0.5, 0.5),
  gp = gpar(fontsize = 7.0, fontfamily = "sans", col = "#5D6876")
)
dev.off()

message("Saved: ", out_tiff)


