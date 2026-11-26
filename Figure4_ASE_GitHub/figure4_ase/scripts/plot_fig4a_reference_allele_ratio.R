args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript scripts/plot_fig4a_reference_allele_ratio.R SAMPLE_MAP COUNT_DIR OUTDIR"
if (length(args) == 1L && args[1] == "--help") {
  cat(usage, "\n")
  quit(save = "no")
}
if (length(args) != 3L) stop(usage, call. = FALSE)

sample_map <- args[1]
count_dir <- args[2]
outdir <- args[3]
if (!dir.exists(count_dir)) stop("Count directory does not exist: ", count_dir)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

library(ggplot2)

# Read the sample map linking ear tag IDs to sample IDs.
sample_info <- read.table(
  sample_map,
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Extract ear tag IDs.
if (!"Ear_ID" %in% names(sample_info)) stop("Sample map must contain Ear_ID.")
ear_id <- sample_info$Ear_ID

# Initialize storage for per-sample records.
all_data_list <- list()
data_index <- 1

# Read the records for each sample.
for (i in ear_id) {

  file_pattern <- paste0(
    paste0(count_dir, "/"),
    i,
    "*ASEReadCounter_binom_test.txt"
  )

  # Find files matching the ear tag ID.
  matched_files <- Sys.glob(file_pattern)

  # Skip samples with no matching file.
  if (length(matched_files) == 0) {
    warning("No matching file for sample: ", i)
    next
  }

  # Read all matching files when a sample has multiple matches.
  for (input_file in matched_files) {

    data <- read.table(
      input_file,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    # Check that the ref_ratio column is present.
    if (!"ref_ratio" %in% colnames(data)) {
      warning(
        "Skipping file without a ref_ratio column: ",
        input_file
      )
      next
    }

    # Convert the reference-allele ratios to numeric values.
    data$ref_ratio <- suppressWarnings(
      as.numeric(data$ref_ratio)
    )

    # Exclude missing values and ratios outside [0, 1].
    data <- data[
      !is.na(data$ref_ratio) &
        data$ref_ratio >= 0 &
        data$ref_ratio <= 1,
      ,
      drop = FALSE
    ]

    if (nrow(data) == 0) {
      warning(
        "Skipping file with no valid ref_ratio values: ",
        input_file
      )
      next
    }

    # Attach sample information.
    data$Ear_ID <- i
    data$source_file <- basename(input_file)

    # Store the records for pooling.
    all_data_list[[data_index]] <- data
    data_index <- data_index + 1
  }
}

# Check that at least one file contributed valid records.
if (length(all_data_list) == 0) {
  stop("No valid records were loaded. Check the sample IDs and input paths.")
}

# Pool the records across samples.
all_data <- do.call(
  rbind,
  all_data_list
)

# Summarize ratios to investigate the peak at 0.5.
n_total <- nrow(all_data)
n_equal_half <- sum(all_data$ref_ratio == 0.5)
percentage_equal_half <- n_equal_half / n_total * 100

message("Total valid records: ", n_total)
message(
  "Records with ref_ratio exactly equal to 0.5: ",
  n_equal_half,
  "（",
  round(percentage_equal_half, 2),
  "%）"
)

# Plot the pooled reference-allele ratios.
sp <- ggplot(
  all_data,
  aes(x = ref_ratio)
) +

  # Histogram
  geom_histogram(
    aes(
      y = after_stat(density),

      # Color each bin according to its center.
      fill = after_stat(
        ifelse(
          x < 0.495,
          "Alternative allele higher",
          ifelse(
            x <= 0.505,
            "Balanced",
            "Reference allele higher"
          )
        )
      )
    ),

    # Use a bin width of 0.01.
    binwidth = 0.01,

    # Center the middle bin at 0.5.
    center = 0.5,

    colour = "black",
    linewidth = 0.25
  ) +

  # Set the colors for the three bin groups.
  scale_fill_manual(
    name = NULL,
    values = c(
      "Alternative allele higher" = "#0057FF",
      "Balanced" = "grey70",
      "Reference allele higher" = "#FF1A1A"
    ),
    breaks = c(
      "Alternative allele higher",
      "Balanced",
      "Reference allele higher"
    ),
    labels = c(
      "Reference allele ratio < 0.5",
      "Balanced bin around 0.5",
      "Reference allele ratio > 0.5"
    )
  ) +

  # Draw the density curve as a black line without a fill.
  geom_density(
    colour = "black",
    fill = NA,
    linewidth = 0.7,
    adjust = 1,
    trim = TRUE
  ) +

  # Mark the reference ratio of 0.5.
  geom_vline(
    xintercept = 0.5,
    colour = "grey30",
    linetype = "dashed",
    linewidth = 0.5
  ) +

  # Configure the x-axis.
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = sprintf("%.2f", seq(0, 1, by = 0.25)),
    limits = c(0, 1),
    expand = expansion(mult = c(0.01, 0.01))
  ) +

  # Set the axis titles.
  labs(
    x = "Reference allele ratio",
    y = "Density"
  ) +

  # Set the plot theme.
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(
      colour = "black",
      size = 12
    ),
    axis.text = element_text(
      colour = "black",
      size = 10
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    legend.position = "top",
    legend.text = element_text(size = 9),
    legend.key.width = grid::unit(0.8, "cm"),
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )

# Use a distinct output filename to preserve the original plot.
output_file <- file.path(
  outdir,
  "histogram_density_of_all_sample_balanced_bin.pdf"
)

ggsave(
  filename = output_file,
  plot = sp,
  width = 8,
  height = 6,
  units = "in",
  device = cairo_pdf
)

message("Figure saved to: ", output_file)