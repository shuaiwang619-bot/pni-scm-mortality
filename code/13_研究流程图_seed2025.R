set.seed(2025)

library(grid)

project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE)
preview_mode <- identical(Sys.getenv("FIG1_PREVIEW"), "1")
result_dir <- file.path(project_dir, "results")
if (!dir.exists(result_dir)) {
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
}

out_tiff <- if (preview_mode) {
  preview_out <- Sys.getenv("FIG1_PREVIEW_OUT")
  if (nzchar(preview_out)) preview_out else file.path(tempdir(), "figure1_preview.tiff")
} else {
  file.path(result_dir, "研究流程图_含MIMIC纳排.tiff")
}
submission_tiff <- file.path(project_dir, "results", "Figure 1-流程图.tiff")

wrap_label <- function(x, width = 32) {
  parts <- unlist(strsplit(x, "\n", fixed = TRUE))
  paste(unlist(lapply(parts, function(z) strwrap(z, width = width))), collapse = "\n")
}

draw_box <- function(x, y, w, h, label, fill = "white", border = "black",
                     fontsize = 7.4, fontface = "plain", r = unit(0.10, "snpc")) {
  grid.roundrect(
    x = x, y = y, width = w, height = h, r = r,
    gp = gpar(fill = fill, col = border, lwd = 1.0)
  )
  grid.text(
    wrap_label(label),
    x = x, y = y,
    gp = gpar(fontsize = fontsize, fontface = fontface, fontfamily = "serif", lineheight = 0.92)
  )
}

draw_note <- function(x, y, label, fontsize = 6.4, col = "#444444", hjust = 0) {
  grid.text(
    wrap_label(label, width = 25),
    x = x, y = y, just = c(hjust, 0.5),
    gp = gpar(fontsize = fontsize, col = col, fontfamily = "serif", lineheight = 0.92)
  )
}

draw_arrow <- function(x0, y0, x1, y1, col = "black", lwd = 0.9) {
  grid.lines(
    x = unit(c(x0, x1), "npc"),
    y = unit(c(y0, y1), "npc"),
    gp = gpar(col = col, lwd = lwd),
    arrow = arrow(length = unit(0.08, "inches"), type = "closed")
  )
}

draw_elbow <- function(x0, y0, x1, y1, via_x = NULL, col = "black", lwd = 0.9) {
  if (is.null(via_x)) via_x <- (x0 + x1) / 2
  grid.lines(
    x = unit(c(x0, via_x, via_x, x1), "npc"),
    y = unit(c(y0, y0, y1, y1), "npc"),
    gp = gpar(col = col, lwd = lwd),
    arrow = arrow(length = unit(0.08, "inches"), type = "closed")
  )
}

draw_panel_title <- function(x, y, label) {
  grid.text(label, x = x, y = y, just = c(0, 0.5),
            gp = gpar(fontsize = 11, fontface = "bold", fontfamily = "serif"))
}

draw_backplate <- function(x, y, w, h, fill, border = NA) {
  grid.roundrect(
    x = x, y = y, width = w, height = h, r = unit(0.012, "snpc"),
    gp = gpar(fill = fill, col = border, lwd = 0.6)
  )
}

tiff(
  out_tiff,
  width = 11.0, height = 7.0, units = "in",
  res = if (preview_mode) 180 else 600,
  compression = "lzw"
)
grid.newpage()
grid.rect(gp = gpar(fill = "white", col = NA))

draw_backplate(0.245, 0.500, 0.430, 0.850, "#FBFBFA")
draw_backplate(0.720, 0.500, 0.445, 0.850, "#FBFBFA")
draw_backplate(0.675, 0.510, 0.300, 0.655, "#F6F9FB")
draw_backplate(0.925, 0.690, 0.145, 0.330, "#FFF9E6")
draw_backplate(0.925, 0.330, 0.145, 0.120, "#FFF9E6")

draw_panel_title(0.04, 0.965, "A. MIMIC-IV cohort construction")
draw_panel_title(0.52, 0.965, "B. Analysis workflow")

grid.lines(x = unit(c(0.49, 0.49), "npc"), y = unit(c(0.06, 0.94), "npc"),
           gp = gpar(col = "#BFBFBF", lwd = 0.8))

left_x <- 0.210
left_w <- 0.265
left_h <- 0.058
left_y <- c(0.895, 0.805, 0.715, 0.625, 0.535, 0.445, 0.355, 0.265, 0.155)
left_heights <- c(rep(left_h, 6), 0.074, 0.074, 0.095)
left_fontsizes <- c(rep(7.5, 6), 6.9, 6.9, 7.0)
left_labels <- c(
  "All ICU stays in MIMIC-IV\nn = 94,458",
  "First ICU stay per patient\nn = 65,366",
  "ICU length of stay >24 h\nn = 51,837",
  "Sepsis-3 positive ICU stays\nn = 25,105",
  "No prior heart failure,\ncardiomyopathy, or myocarditis\nn = 17,813",
  "Alive at 72-h landmark\nn = 17,076",
  "Clinically obtained TTE-LVEF\nbefore or at 72-h landmark\n(ICU -24 h to +72 h)\nn = 4,669",
  "Albumin and lymphocyte available\nfor early PNI (ICU -24 h to +48 h)\nn = 1,976",
  "Final MIMIC-IV analytic cohort\nn = 1,976\nNo SCM = 1,736; SCM = 240\n28-day deaths = 506; 90-day deaths = 633"
)
removed_labels <- c(
  "Excluded duplicate ICU stays\nn = 29,092",
  "Age <18 years\nn = 0",
  "ICU stay <=24 h\nn = 13,529",
  "No Sepsis-3 diagnosis\nn = 26,732",
  "Prior cardiac disease excluded\nn = 7,292",
  "Died before landmark\nn = 737",
  "No eligible TTE-LVEF\nn = 12,407",
  "PNI components unavailable\nn = 2,693"
)

for (i in seq_along(left_y)) {
  fill <- if (i == length(left_y)) "#EAF2F8" else "white"
  draw_box(left_x, left_y[i], left_w, left_heights[i],
           left_labels[i], fill = fill, fontsize = left_fontsizes[i],
           fontface = if (i == length(left_y)) "bold" else "plain")
  if (i < length(left_y)) {
    y0 <- left_y[i] - left_heights[i] / 2
    y1 <- left_y[i + 1] + left_heights[i + 1] / 2
    draw_arrow(left_x, y0, left_x, y1)
    draw_note(0.365, (y0 + y1) / 2, removed_labels[i], fontsize = 6.2)
  }
}

draw_box(0.610, 0.890, 0.160, 0.056, "MIMIC-IV analytic cohort\nn = 1,976", fill = "#EAF2F8", fontsize = 7.1, fontface = "bold")
draw_box(0.925, 0.920, 0.135, 0.052, "Single-center external\nsepsis cohort\nn = 235", fill = "#FFF2CC", fontsize = 6.0, fontface = "bold")
draw_box(0.925, 0.835, 0.135, 0.066, "SCM status and PNI derived\nPNI components complete-case\ncovariates harmonized", fill = "#FFF2CC", fontsize = 5.3)
draw_box(0.925, 0.750, 0.135, 0.072, "External exploratory cohort\nNo SCM = 177; SCM = 58\nin-hospital deaths = 77", fill = "#FFF2CC", fontsize = 5.4, fontface = "bold")

draw_arrow(0.925, 0.894, 0.925, 0.868)
draw_arrow(0.925, 0.802, 0.925, 0.786)

draw_box(0.685, 0.775, 0.180, 0.052, "Baseline characteristics\nTable 1", fill = "white", fontsize = 7.1, fontface = "bold")
draw_box(0.685, 0.660, 0.185, 0.058, "Primary Cox models\n28-day and 90-day mortality", fill = "white", fontsize = 7.0, fontface = "bold")

draw_arrow(0.610, 0.862, 0.685, 0.801)
draw_arrow(0.858, 0.750, 0.792, 0.750)
draw_arrow(0.685, 0.749, 0.685, 0.689)

draw_box(0.685, 0.505, 0.180, 0.116,
         "Main interaction analysis\n72-h landmark follow-up\ncentered lower PNI x SCM\nfully adjusted Model 3",
         fill = "#F7F7F7", fontsize = 6.9, fontface = "bold")
draw_arrow(0.685, 0.631, 0.685, 0.563)

draw_box(0.545, 0.505, 0.105, 0.078, "Diagnostics\nSchoenfeld PH\nVIF and Spearman\nSpline checks", fill = "white", fontsize = 5.9)
draw_arrow(0.595, 0.505, 0.630, 0.505)

draw_box(0.815, 0.505, 0.115, 0.080, "Sensitivity analyses\nAdjustment sets\nTime windows\nWinsorization", fill = "white", fontsize = 5.9)
draw_arrow(0.775, 0.505, 0.758, 0.505)

draw_box(0.575, 0.330, 0.125, 0.064, "Nonlinear analysis\nNatural cubic spline\nPNI reference = 35", fill = "white", fontsize = 6.0)
draw_elbow(0.655, 0.447, 0.575, 0.363, via_x = 0.655)

draw_box(0.685, 0.330, 0.130, 0.064, "Absolute-risk translation\nPNI 30 vs 40\nbootstrap risk differences", fill = "white", fontsize = 5.9)
draw_arrow(0.685, 0.447, 0.685, 0.363)

draw_box(0.925, 0.330, 0.135, 0.070, "External exploratory\nconsistency\nFirth logistic regression\npenalized LRT", fill = "#FFF2CC", fontsize = 5.3)
draw_arrow(0.925, 0.714, 0.925, 0.365)

draw_box(0.685, 0.180, 0.220, 0.066, "Component comparison\nPNI vs albumin vs lymphocyte\nAIC and C-index", fill = "white", fontsize = 6.5)
draw_arrow(0.685, 0.298, 0.685, 0.213)

draw_elbow(0.575, 0.298, 0.620, 0.213, via_x = 0.575)
draw_elbow(0.740, 0.298, 0.750, 0.213, via_x = 0.740)

grid.text(
  "Index time, ICU admission; landmark, ICU admission +72 h. PNI = albumin (g/L) + 5 x lymphocyte count (10^9/L).\nPNI components were complete-case; eligible non-exposure covariates were imputed by MICE.\nMIMIC-IV SCM: LVEF <50% or >10-point LVEF decline to <55%; external cohort SCM: available LVEF <50%.",
  x = 0.04, y = 0.033, just = c(0, 0.5),
  gp = gpar(fontsize = 5.4, col = "#555555", fontfamily = "serif", lineheight = 0.9)
)

dev.off()

if (!preview_mode && dir.exists(dirname(submission_tiff))) {
  file.copy(out_tiff, submission_tiff, overwrite = TRUE)
}

cat("Study flow chart saved to:", out_tiff, "\n")
cat("Submission copy saved to:", submission_tiff, "\n")



