set.seed(2025)

library(grid)

project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE)
result_dir <- file.path(project_dir, "results")
if (!dir.exists(result_dir)) {
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
}

stats_file <- file.path(project_dir, "results", "tables", "missingness_summary.csv")
out_tiff <- file.path(result_dir, "双队列缺失值统计图.tiff")

lines <- readLines(stats_file, encoding = "UTF-8")
mimic_start <- grep("^mimic", lines, ignore.case = TRUE)[1]
domestic_start <- grep("^LinYincohort", lines, ignore.case = TRUE)[1]
if (is.na(mimic_start) || is.na(domestic_start)) {
  stop("Could not identify cohort sections in missingness CSV.")
}

mimic_text <- paste(lines[(mimic_start + 1):(domestic_start - 1)], collapse = "\n")
domestic_text <- paste(lines[(domestic_start + 1):length(lines)], collapse = "\n")

mimic <- read.csv(text = mimic_text, stringsAsFactors = FALSE, check.names = FALSE)
domestic <- read.csv(text = domestic_text, stringsAsFactors = FALSE, check.names = FALSE)

mimic_plot <- data.frame(
  variable = mimic$table1_variable,
  missing_pct = as.numeric(mimic$missing_pct),
  stringsAsFactors = FALSE
)
mimic_plot <- mimic_plot[mimic_plot$missing_pct > 0, ]
mimic_plot <- mimic_plot[order(mimic_plot$missing_pct, decreasing = TRUE), ]

domestic_label_map <- c(
  "il_6" = "IL-6",
  "c反应蛋白" = "C-reactive protein",
  "sofa" = "SOFA score",
  "apache" = "APACHE score",
  "lac" = "Lactate",
  "氧合指数" = "Oxygenation index",
  "肌酸激酶同工酶" = "CK-MB",
  "ntprobnp" = "NT-proBNP",
  "d二聚体" = "D-dimer",
  "ards" = "ARDS",
  "bmi" = "BMI"
)

domestic_raw_variable <- domestic$variable
domestic_labels <- domestic_label_map[domestic_raw_variable]
domestic_labels[is.na(domestic_labels)] <- domestic_raw_variable[is.na(domestic_labels)]

domestic_plot <- data.frame(
  variable = domestic_labels,
  missing_pct = as.numeric(domestic$missing_pct),
  stringsAsFactors = FALSE
)
domestic_plot <- domestic_plot[domestic_plot$missing_pct > 0, ]
domestic_plot <- domestic_plot[order(domestic_plot$missing_pct, decreasing = TRUE), ]

fmt_pct <- function(x) sprintf("%.2f%%", x)

draw_card <- function(x, y, w, h, title, data, bar_color, x_max = 32) {
  grid.roundrect(
    x = x, y = y, width = w, height = h, r = unit(0.012, "snpc"),
    gp = gpar(fill = "white", col = "#D7DEE8", lwd = 0.7)
  )

  left_pad <- 0.025 * w
  label_w <- 0.185 * w
  right_pad <- 0.025 * w
  top_pad <- 0.105 * h
  bottom_pad <- 0.125 * h

  plot_x0 <- x - w / 2 + left_pad + label_w
  plot_x1 <- x + w / 2 - right_pad
  plot_y0 <- y - h / 2 + bottom_pad
  plot_y1 <- y + h / 2 - top_pad
  plot_w <- plot_x1 - plot_x0
  plot_h <- plot_y1 - plot_y0

  grid.text(
    title,
    x = x - w / 2 + 0.016 * w,
    y = y + h / 2 - 0.040 * h,
    just = c(0, 1),
    gp = gpar(fontsize = 9.2, fontface = "bold", fontfamily = "sans", col = "#1F2630")
  )

  breaks <- seq(0, 30, by = 5)
  for (b in breaks) {
    xb <- plot_x0 + plot_w * b / x_max
    grid.lines(
      x = unit(c(xb, xb), "npc"),
      y = unit(c(plot_y0, plot_y1), "npc"),
      gp = gpar(col = ifelse(b == 0, "#9BA7B5", "#E0E6EF"), lwd = ifelse(b == 0, 0.75, 0.55))
    )
  }

  x_ref <- plot_x0 + plot_w * 30 / x_max
  grid.lines(
    x = unit(c(x_ref, x_ref), "npc"),
    y = unit(c(plot_y0, plot_y1), "npc"),
    gp = gpar(col = "#9C2F2F", lwd = 0.8, lty = "dashed")
  )
  grid.text(
    "30% reference",
    x = x_ref + 0.004,
    y = plot_y1 + 0.030 * h,
    just = c(0, 0),
    gp = gpar(fontsize = 6.5, fontfamily = "sans", col = "#9C2F2F")
  )

  n <- nrow(data)
  if (n == 0) return(invisible(NULL))
  row_h <- plot_h / n
  for (i in seq_len(n)) {
    yc <- plot_y1 - (i - 0.5) * row_h
    if (i %% 2 == 1) {
      grid.rect(
        x = x,
        y = yc,
        width = w - 2 * left_pad,
        height = row_h * 0.94,
        gp = gpar(fill = "#F3F7FB", col = NA)
      )
    }
    grid.text(
      data$variable[i],
      x = x - w / 2 + left_pad + 0.006 * w,
      y = yc,
      just = c(0, 0.5),
      gp = gpar(fontsize = 7.6, fontfamily = "sans", col = "#1F2630")
    )
    bar_x1 <- plot_x0 + plot_w * data$missing_pct[i] / x_max
    bar_w <- max(bar_x1 - plot_x0, 0.004)
    grid.roundrect(
      x = plot_x0 + bar_w / 2,
      y = yc,
      width = bar_w,
      height = min(row_h * 0.34, 0.020),
      r = unit(0.004, "snpc"),
      gp = gpar(fill = bar_color, col = NA)
    )
    grid.text(
      fmt_pct(data$missing_pct[i]),
      x = bar_x1 + 0.006,
      y = yc,
      just = c(0, 0.5),
      gp = gpar(fontsize = 7.0, fontfamily = "sans", col = "#1F2630")
    )
  }

  grid.lines(
    x = unit(c(plot_x0, plot_x1), "npc"),
    y = unit(c(plot_y0, plot_y0), "npc"),
    gp = gpar(col = "#8794A4", lwd = 0.75)
  )
  for (b in breaks) {
    xb <- plot_x0 + plot_w * b / x_max
    grid.lines(
      x = unit(c(xb, xb), "npc"),
      y = unit(c(plot_y0, plot_y0 - 0.008), "npc"),
      gp = gpar(col = "#8794A4", lwd = 0.55)
    )
    grid.text(
      as.character(b),
      x = xb,
      y = plot_y0 - 0.018,
      just = c(0.5, 1),
      gp = gpar(fontsize = 6.8, fontfamily = "sans", col = "#344054")
    )
  }
}

tiff(out_tiff, width = 8.8, height = 6.35, units = "in", res = 600, compression = "lzw")
grid.newpage()
grid.rect(gp = gpar(fill = "#EAF0F7", col = NA))

grid.text(
  "Variable-wise missingness rates before imputation in the MIMIC-IV and domestic exploratory cohorts",
  x = 0.5, y = 0.965, just = c(0.5, 0.5),
  gp = gpar(fontsize = 12.0, fontface = "bold", fontfamily = "sans", col = "#1F2630")
)
grid.text(
  "Only variables with nonzero missingness are shown; dashed line indicates a 30% reference threshold.",
  x = 0.5, y = 0.932, just = c(0.5, 0.5),
  gp = gpar(fontsize = 7.6, fontfamily = "sans", col = "#5D6876")
)

draw_card(
  x = 0.5, y = 0.715, w = 0.88, h = 0.285,
  title = "MIMIC-IV cohort (n = 1,976)",
  data = mimic_plot,
  bar_color = "#4169E1"
)

draw_card(
  x = 0.5, y = 0.390, w = 0.88, h = 0.440,
  title = "Domestic exploratory cohort (n = 235)",
  data = domestic_plot,
  bar_color = "#D2812B"
)

grid.text(
  "Missingness (%)",
  x = 0.5, y = 0.070,
  gp = gpar(fontsize = 8.6, fontface = "bold", fontfamily = "sans", col = "#1F2630")
)

grid.roundrect(
  x = 0.5, y = 0.030, width = 0.34, height = 0.030,
  r = unit(0.006, "snpc"),
  gp = gpar(fill = "#F5F8FC", col = "#DDE4EE", lwd = 0.4)
)
grid.roundrect(x = 0.39, y = 0.030, width = 0.014, height = 0.010, r = unit(0.003, "snpc"),
               gp = gpar(fill = "#4169E1", col = NA))
grid.text("MIMIC-IV", x = 0.405, y = 0.030, just = c(0, 0.5),
          gp = gpar(fontsize = 7.2, fontfamily = "sans", col = "#1F2630"))
grid.roundrect(x = 0.515, y = 0.030, width = 0.014, height = 0.010, r = unit(0.003, "snpc"),
               gp = gpar(fill = "#D2812B", col = NA))
grid.text("Domestic", x = 0.530, y = 0.030, just = c(0, 0.5),
          gp = gpar(fontsize = 7.2, fontfamily = "sans", col = "#1F2630"))

dev.off()

cat("Missingness figure saved to:", out_tiff, "\n")



