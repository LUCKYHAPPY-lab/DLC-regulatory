# ============================================================
# Fig. 3a
# Left axis: -log10(BBFDR) for the lead cis-eQTL.
# Right axis: kernel density of distances to the TSS.
# Input: precomputed scatterplot data.
# Output: use a separate directory to preserve the original figure.
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript scripts/plot_fig3a_lead_eqtl_tss_distance.R SCATTER_TSV OUTDIR"
if (length(args) == 1L && args[1] == "--help") {
  cat(usage, "\n")
  quit(save = "no")
}
if (length(args) != 2L) stop(usage, call. = FALSE)
input_file <- args[1]
out_dir <- args[2]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(out_dir, "Fig3a_distance_density_right_axis.pdf")
out_png <- file.path(out_dir, "Fig3a_distance_density_right_axis.png")
out_density <- file.path(out_dir, "Fig3a_density_curve.tsv")

# -------------------------
# 1. Read the precomputed scatterplot data
# -------------------------

if (!file.exists(input_file)) {
  stop("Scatterplot input file not found: ", input_file)
}

dat <- read.delim(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c("Gene.ID", "pos_to_TSS", "BBFDR")
missing_cols <- setdiff(required_cols, names(dat))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

dat$pos_to_TSS <- as.numeric(dat$pos_to_TSS)
dat$BBFDR <- as.numeric(dat$BBFDR)

dat <- dat[
  is.finite(dat$pos_to_TSS) &
    is.finite(dat$BBFDR) &
    dat$BBFDR >= 0 &
    dat$BBFDR <= 1,
]

dat$distance_Mb <- dat$pos_to_TSS / 1e6
dat <- dat[abs(dat$distance_Mb) <= 1, ]

if (nrow(dat) < 2) {
  stop("At least two valid records are required.")
}

if (anyDuplicated(dat$Gene.ID)) {
  stop("Duplicate Gene.ID values found. Provide one lead eQTL per gene.")
}

dat$neglog10_BBFDR <- -log10(pmax(dat$BBFDR, 1e-300))

# -------------------------
# 2. Estimate kernel density
# Use the original default settings of density().
# -------------------------

den <- density(
  dat$distance_Mb,
  from = -1,
  to = 1,
  n = 512,
  na.rm = TRUE
)

write.table(
  data.frame(
    distance_to_TSS_Mb = den$x,
    density_per_Mb = den$y
  ),
  file = out_density,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------
# 3. Plot settings
# -------------------------

point_color <- "#D97A2B"
curve_color <- "#2C9FC2"

point_alpha <- 0.55
point_cex <- 0.34
curve_lwd <- 2.2

fig_width <- 4.2
fig_height <- 3.6

left_max <- max(dat$neglog10_BBFDR) * 1.12

# Leave room for the legend above the plot.
right_max <- max(den$y) / 0.82

# -------------------------
# 4. Plotting function
# -------------------------

draw_plot <- function() {

  par(
    mar = c(4.2, 4.4, 0.8, 4.2),
    mgp = c(2.5, 0.65, 0),
    tcl = -0.25,
    las = 1,
    family = "Helvetica"
  )

  # Draw the scatterplot against the left axis.
  plot(
    dat$distance_Mb,
    dat$neglog10_BBFDR,
    pch = 16,
    cex = point_cex,
    col = adjustcolor(point_color, alpha.f = point_alpha),
    xlim = c(-1, 1),
    ylim = c(0, left_max),
    xaxs = "i",
    yaxs = "i",
    xlab = "Distance to TSS (Mb)",
    ylab = expression(-log[10](BBFDR)),
    axes = FALSE,
    bty = "n"
  )

  axis(
    1,
    at = seq(-1, 1, by = 0.5),
    labels = c("-1.0", "-0.5", "0", "0.5", "1.0"),
    cex.axis = 0.9
  )

  axis(2, cex.axis = 0.9)

  # Draw the unscaled density values against the right axis.
  par(new = TRUE)

  plot(
    den$x,
    den$y,
    type = "l",
    col = curve_color,
    lwd = curve_lwd,
    xlim = c(-1, 1),
    ylim = c(0, right_max),
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    xlab = "",
    ylab = "",
    bty = "n"
  )

  density_ticks <- pretty(c(0, right_max), n = 5)
  density_ticks <- density_ticks[
    density_ticks >= 0 & density_ticks <= right_max
  ]

  axis(
    4,
    at = density_ticks,
    labels = format(
      density_ticks,
      trim = TRUE,
      scientific = FALSE
    ),
    col = curve_color,
    col.axis = curve_color,
    col.ticks = curve_color,
    cex.axis = 0.9
  )

  mtext(
    "Density",
    side = 4,
    line = 2.5,
    col = curve_color,
    las = 0
  )

  box()

  legend(
    "topright",
    legend = c("Lead cis-eQTL", "Density"),
    col = c(point_color, curve_color),
    pch = c(16, NA),
    lty = c(NA, 1),
    lwd = c(NA, curve_lwd),
    pt.cex = 0.8,
    cex = 0.76,
    bty = "n",
    inset = 0.02,
    y.intersp = 1.1
  )
}

# -------------------------
# 5. Save a transparent vector PDF
# -------------------------

pdf(
  out_pdf,
  width = fig_width,
  height = fig_height,
  family = "Helvetica",
  bg = "transparent",
  useDingbats = FALSE
)

draw_plot()
dev.off()

# -------------------------
# 6. Save a transparent PNG
# -------------------------

png(
  out_png,
  width = fig_width,
  height = fig_height,
  units = "in",
  res = 600,
  bg = "transparent"
)

draw_plot()
dev.off()

cat(
  "\nCompleted!\n",
  "Sites plotted: ", nrow(dat), "\n",
  "PDF: ", out_pdf, "\n",
  "PNG: ", out_png, "\n",
  "Kernel density data: ", out_density, "\n",
  sep = ""
)
