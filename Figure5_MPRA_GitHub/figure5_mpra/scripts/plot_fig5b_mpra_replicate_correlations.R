#!/usr/bin/env Rscript

# Version: correlation_v4 (2/3-decimal labels + improved correlation gradient)

# ============================================================================
# Four MPRA datasets: barcode counts -> oligo counts -> DESeq2 + summit shift
# -> two log2 transformations -> Pearson correlation heatmaps
#
# Correlation transformations generated for every heatmap:
#   1) log2(normalized oligo count + 1)
#   2) log2(normalized oligo count), with Oligos containing a zero in any
#      sample of that heatmap excluded before transformation
#
# Additional library-level heatmaps:
#   Plasmid (3 replicates) -> PK15 RNA (3) -> PSC RNA (3)
#   one PDF for the 3'UTR library and one PDF for the CRE library
#
# Additional joint-normalization RNA heatmaps:
#   for each library, jointly normalize the 3 Plasmid + 3 PK15 RNA + 3 PSC RNA
#   samples, extract the six RNA columns, and draw a clustered 6 x 6 heatmap.
#
# The official MPRASuite functions used here are:
#   1) oligoIsolate(): sum barcode counts for the same Oligo
#   2) tagNorm(method = "ss"): DESeq2 normalization + summit-shift adjustment
#
# Heatmap colour rule:
#   continuous light-blue -> near-white -> light-red -> vivid-red gradient
#   r >= 0.60 starts from light red; values above 0.60 become progressively redder
# ============================================================================

options(stringsAsFactors = FALSE, scipen = 999, digits = 15)

# ----------------------------------------------------------------------------
# 1. Paths and parameters
# ----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage: Rscript scripts/plot_fig5b_mpra_replicate_correlations.R",
  "MPRA_MODEL CRE_PK15 CRE_PSC UTR_PK15 UTR_PSC OUTDIR"
)
if (length(args) == 1L && args[1] == "--help") {
  cat(usage, "\n")
  quit(save = "no")
}
if (length(args) != 6L) stop(usage, call. = FALSE)
MPRA_MODEL_FILE <- args[1]
OUTPUT_DIR <- args[6]

# This is a colour anchor, not a binary threshold.
CORRELATION_COLOR_ANCHOR <- 0.60

# Heatmap labels: keep both publication-friendly and detailed versions.
# Source matrices always retain full precision.
CORRELATION_LABEL_FORMAT_3 <- "%.3f"
CORRELATION_LABEL_FORMAT_2 <- "%.2f"

# Save the barcode-summed raw oligo count matrices as intermediate files.
WRITE_RAW_OLIGO_COUNTS <- TRUE

# Usually unnecessary because the normalized matrices are already saved.
# Change to TRUE if the log2-transformed matrices are also needed.
WRITE_LOG2_MATRICES <- FALSE

# File tags are deliberately explicit so the +1 and no-pseudocount outputs
# cannot be confused downstream.
CORRELATION_TRANSFORMS <- list(
  log2norm_plus1 = list(
    key = "log2norm_plus1",
    label = "log2(normalized oligo count + 1)",
    add_one = TRUE
  ),
  log2norm_no_plus1 = list(
    key = "log2norm_no_plus1",
    label = "log2(normalized oligo count)",
    add_one = FALSE
  )
)

DATASETS <- list(
  `3UTR_PK15` = list(
    input = args[4],
    cell_type = "PK15",
    title = "3'UTR library - PK15"
  ),
  `3UTR_PSC` = list(
    input = args[5],
    cell_type = "PSC",
    title = "3'UTR library - PSC"
  ),
  `CRE_PK15` = list(
    input = args[2],
    cell_type = "PK15",
    title = "CRE library - PK15"
  ),
  `CRE_PSC` = list(
    input = args[3],
    cell_type = "PSC",
    title = "CRE library - PSC"
  )
)


# ----------------------------------------------------------------------------
# 2. Initial checks and official MPRASuite code
# ----------------------------------------------------------------------------

if (!file.exists(MPRA_MODEL_FILE)) {
  stop("Official MPRAmodel.R was not found: ", MPRA_MODEL_FILE)
}

missing_inputs <- vapply(
  DATASETS,
  function(x) !file.exists(x$input),
  logical(1)
)

if (any(missing_inputs)) {
  stop(
    "The following input files were not found:\n",
    paste(
      vapply(DATASETS[missing_inputs], `[[`, character(1), "input"),
      collapse = "\n"
    )
  )
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

write.table(
  data.frame(
    Transform_file_tag = vapply(
      CORRELATION_TRANSFORMS,
      `[[`,
      character(1),
      "key"
    ),
    Transformation = vapply(
      CORRELATION_TRANSFORMS,
      `[[`,
      character(1),
      "label"
    ),
    Zero_count_handling = c(
      "All Oligos retained; zero becomes log2(0 + 1) = 0",
      paste0(
        "Exclude Oligos with a zero in any sample displayed in the ",
        "heatmap; no pseudocount is added"
      )
    ),
    Primary_or_supplementary = c(
      "Primary; matches the article public QC-code transformation",
      "Supplementary sensitivity comparison"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ),
  file = file.path(OUTPUT_DIR, "Correlation_transform_definitions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

message("Loading the official MPRASuite code: ", MPRA_MODEL_FILE)
source(MPRA_MODEL_FILE)

required_functions <- c("oligoIsolate", "tagNorm")
missing_functions <- required_functions[
  !vapply(required_functions, exists, logical(1), mode = "function")
]

if (length(missing_functions) > 0L) {
  stop(
    "The following required functions were not found after sourcing MPRAmodel.R: ",
    paste(missing_functions, collapse = ", ")
  )
}

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("The R package DESeq2 is required but is not installed.")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The R package ggplot2 is required but is not installed.")
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop(
    "The R package pheatmap is required for the combined RNA heatmaps. ",
    "Install it with: install.packages('pheatmap')"
  )
}


# ----------------------------------------------------------------------------
# 3. Helper functions
# ----------------------------------------------------------------------------

read_barcode_count_file <- function(path, cell_type) {
  message("Reading barcode-level count file: ", path)

  # data.table is used only for faster file reading when available. It does not
  # alter the MPRASuite aggregation or normalization method.
  if (requireNamespace("data.table", quietly = TRUE)) {
    dat <- data.table::fread(
      file = path,
      sep = "\t",
      header = TRUE,
      data.table = FALSE,
      showProgress = TRUE
    )
  } else {
    dat <- read.delim(
      file = path,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      quote = "",
      comment.char = ""
    )
  }

  names(dat) <- trimws(names(dat))

  if (anyDuplicated(names(dat))) {
    duplicated_names <- unique(names(dat)[duplicated(names(dat))])
    stop(
      "Duplicated column names in ", path, ": ",
      paste(duplicated_names, collapse = ", ")
    )
  }

  required_columns <- c("Barcode", "Oligo")
  missing_columns <- setdiff(required_columns, names(dat))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required column(s) in ", path, ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  sample_columns <- setdiff(names(dat), required_columns)
  if (length(sample_columns) < 2L) {
    stop("Fewer than two sample count columns were found in: ", path)
  }

  sample_names_lower <- tolower(sample_columns)
  dna_columns <- sample_columns[grepl("^plasmid_", sample_names_lower)]
  rna_columns <- sample_columns[
    grepl(paste0("^", tolower(cell_type), "_"), sample_names_lower)
  ]
  unknown_columns <- setdiff(sample_columns, c(dna_columns, rna_columns))

  if (length(dna_columns) < 2L) {
    stop(
      "Fewer than two plasmid columns were detected in ", path,
      ". Detected: ", paste(dna_columns, collapse = ", ")
    )
  }

  if (length(rna_columns) < 2L) {
    stop(
      "Fewer than two ", cell_type, " RNA columns were detected in ", path,
      ". Detected: ", paste(rna_columns, collapse = ", ")
    )
  }

  if (length(unknown_columns) > 0L) {
    stop(
      "Unexpected sample column(s) in ", path, ": ",
      paste(unknown_columns, collapse = ", "),
      "\nExpected only plasmid_* and ", cell_type, "_* columns."
    )
  }

  dat$Barcode <- as.character(dat$Barcode)
  dat$Oligo <- trimws(as.character(dat$Oligo))

  if (any(is.na(dat$Oligo) | dat$Oligo == "")) {
    stop("Missing or empty Oligo IDs were found in: ", path)
  }

  # Validate and coerce every count column to non-negative integer counts.
  for (column_name in sample_columns) {
    original_values <- dat[[column_name]]
    numeric_values <- suppressWarnings(
      as.numeric(as.character(original_values))
    )

    original_text <- trimws(as.character(original_values))
    conversion_failed <- is.na(numeric_values) &
      !is.na(original_values) & original_text != ""

    if (any(conversion_failed)) {
      stop(
        "Non-numeric count values were found in column ", column_name,
        " of file: ", path
      )
    }

    if (anyNA(numeric_values)) {
      stop(
        "NA count values were found in column ", column_name,
        " of file: ", path
      )
    }

    if (any(numeric_values < 0)) {
      stop(
        "Negative count values were found in column ", column_name,
        " of file: ", path
      )
    }

    if (any(abs(numeric_values - round(numeric_values)) > 1e-8)) {
      stop(
        "Non-integer count values were found in column ", column_name,
        " of file: ", path
      )
    }

    numeric_values <- round(numeric_values)
    if (max(numeric_values, na.rm = TRUE) <= .Machine$integer.max) {
      numeric_values <- as.integer(numeric_values)
    }
    dat[[column_name]] <- numeric_values
  }

  # oligoIsolate() assumes Barcode is column 1 and Oligo is column 2.
  dat <- dat[, c("Barcode", "Oligo", sample_columns), drop = FALSE]

  list(
    data = dat,
    sample_columns = sample_columns,
    dna_columns = dna_columns,
    rna_columns = rna_columns
  )
}


run_official_oligo_isolate <- function(barcode_data, dataset_name) {
  # The official function writes an intermediate file to ./results/. A private
  # temporary directory is used here; the returned official oligo matrix is
  # saved below with a stable, descriptive filename.
  temporary_directory <- tempfile(
    pattern = paste0("MPRASuite_oligoIsolate_", dataset_name, "_")
  )
  dir.create(
    file.path(temporary_directory, "results"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  original_working_directory <- getwd()
  on.exit({
    setwd(original_working_directory)
    unlink(temporary_directory, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  setwd(temporary_directory)
  oligo_counts <- oligoIsolate(barcode_data, dataset_name)

  as.data.frame(oligo_counts, check.names = FALSE)
}


write_matrix_with_oligo_column <- function(matrix_or_df, output_file) {
  output_data <- data.frame(
    Oligo = rownames(matrix_or_df),
    as.data.frame(matrix_or_df, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  connection <- gzfile(output_file, open = "wt")
  on.exit(close(connection), add = TRUE)

  write.table(
    output_data,
    file = connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}


write_square_matrix <- function(square_matrix, output_file) {
  output_data <- data.frame(
    Sample = rownames(square_matrix),
    as.data.frame(square_matrix, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write.table(
    output_data,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}


calculate_log2_correlation <- function(normalized_counts,
                                       transform_specification,
                                       context_label) {
  if (!is.matrix(normalized_counts)) {
    normalized_counts <- as.matrix(normalized_counts)
  }
  if (nrow(normalized_counts) < 2L || ncol(normalized_counts) < 2L) {
    stop(
      context_label,
      " requires at least two Oligos and two samples for correlation."
    )
  }
  if (anyNA(normalized_counts) ||
      any(!is.finite(normalized_counts)) ||
      any(normalized_counts < 0)) {
    stop("Invalid normalized counts supplied for ", context_label, ".")
  }

  total_oligos <- nrow(normalized_counts)

  if (isTRUE(transform_specification$add_one)) {
    keep_rows <- rep(TRUE, total_oligos)
    transformed_counts <- log2(normalized_counts + 1)
  } else {
    # A genuine log2(count) transformation cannot represent zero. Exclude an
    # Oligo if any sample displayed in this heatmap has a zero. This produces
    # one common positive-count Oligo set for all pairwise correlations and
    # avoids both -Inf values and pair-specific feature sets.
    keep_rows <- rowSums(normalized_counts <= 0) == 0L
    if (sum(keep_rows) < 2L) {
      stop(
        "Fewer than two Oligos have positive counts in every sample for ",
        context_label, " using ", transform_specification$label, "."
      )
    }
    transformed_counts <- log2(
      normalized_counts[keep_rows, , drop = FALSE]
    )
  }

  if (anyNA(transformed_counts) || any(!is.finite(transformed_counts))) {
    stop(
      "NA/Inf values remain after ", transform_specification$label,
      " for ", context_label, "."
    )
  }

  sample_standard_deviations <- apply(
    transformed_counts,
    2,
    stats::sd
  )
  if (any(!is.finite(sample_standard_deviations)) ||
      any(sample_standard_deviations == 0)) {
    failed_samples <- names(sample_standard_deviations)[
      !is.finite(sample_standard_deviations) |
        sample_standard_deviations == 0
    ]
    stop(
      "Zero or invalid variance after ", transform_specification$label,
      " for ", context_label, ": ",
      paste(failed_samples, collapse = ", ")
    )
  }

  correlation_matrix <- stats::cor(
    transformed_counts,
    method = "pearson",
    use = "everything"
  )
  if (anyNA(correlation_matrix) || any(!is.finite(correlation_matrix))) {
    stop(
      "The Pearson correlation matrix contains NA/Inf values for ",
      context_label, " using ", transform_specification$label, "."
    )
  }

  list(
    transformed_counts = transformed_counts,
    correlation_matrix = correlation_matrix,
    kept_oligo_ids = rownames(normalized_counts)[keep_rows],
    total_oligos = total_oligos,
    oligos_used = sum(keep_rows),
    oligos_excluded_due_to_zero = total_oligos - sum(keep_rows),
    transform_key = transform_specification$key,
    transform_label = transform_specification$label
  )
}


get_correlation_gradient_spec <- function(colour_anchor = 0.60,
                                     label_format = CORRELATION_LABEL_FORMAT_3) {
  if (colour_anchor <= 0 || colour_anchor >= 1) {
    stop("colour_anchor must be between 0 and 1.")
  }

  # Bright, publication-style continuous gradient. The 0.80 value is a colour
  # anchor, not a binary cut-off: values below it move continuously from light
  # blue through near-white and pale pink; values at and above it progress
  # through light coral to vivid red without using a dark maroon endpoint.
  list(
    colours = c(
      "#DCEFFD",
      "#F7FBFF",
      "#FFFFFF",
      "#FDE0E0",
      "#F7A6A6",
      "#F06164"
    ),
    positions = c(
      0,
      0.35,
      0.55,
      0.60,
      0.80,
      1
    )
  )
}


make_correlation_heatmap <- function(correlation_matrix,
                                     plot_title,
                                     transform_label,
                                     colour_anchor = 0.80,
                                     label_format = "%.3f") {
  if (colour_anchor <= 0 || colour_anchor >= 1) {
    stop("colour_anchor must be between 0 and 1.")
  }

  correlation_long <- as.data.frame(
    as.table(correlation_matrix),
    stringsAsFactors = FALSE
  )
  names(correlation_long) <- c("Sample_y", "Sample_x", "Correlation")

  sample_order <- colnames(correlation_matrix)
  correlation_long$Sample_x <- factor(
    correlation_long$Sample_x,
    levels = sample_order
  )
  correlation_long$Sample_y <- factor(
    correlation_long$Sample_y,
    levels = rev(sample_order)
  )

  # Values close to the 0.80 anchor remain light red. Progressively higher
  # correlations become brighter red, without using a dark maroon endpoint.
  dark_text_cutoff <- (colour_anchor + 1) / 2
  correlation_long$Text_group <- ifelse(
    !is.na(correlation_long$Correlation) &
      correlation_long$Correlation >= dark_text_cutoff,
    "white",
    "black"
  )

  correlation_long$Label <- ifelse(
    is.na(correlation_long$Correlation),
    "NA",
    sprintf(label_format, correlation_long$Correlation)
  )

  text_values <- c("black" = "black", "white" = "white")

  gradient_specification <- get_correlation_gradient_spec(colour_anchor)

  ggplot2::ggplot(
    correlation_long,
    ggplot2::aes(x = Sample_x, y = Sample_y, fill = Correlation)
  ) +
    # "size" is used instead of "linewidth" for compatibility with older
    # ggplot2 versions commonly installed in MPRASuite environments.
    ggplot2::geom_tile(color = "white", size = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = Label, color = Text_group),
      size = 3.5,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradientn(
      colours = gradient_specification$colours,
      values = gradient_specification$positions,
      limits = c(0, 1),
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
      oob = scales::squish,
      na.value = "#BDBDBD",
      name = "Pearson correlation coefficient (r)"
    ) +
    ggplot2::scale_color_manual(values = text_values) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = plot_title,
      subtitle = paste0(
        "Pearson r of ", transform_label,
        "; DESeq2 + summit-shift normalization"
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 0,
        vjust = 0,
        color = "black",
        size = 10
      ),
      axis.text.y = ggplot2::element_text(color = "black", size = 10),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 9),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(12, 12, 12, 12)
    )
}


save_combined_heatmaps <- function(plot_list, output_pdf, output_png) {
  draw_plots <- function() {
    grid::grid.newpage()
    layout <- grid::grid.layout(nrow = 2, ncol = 2)
    grid::pushViewport(grid::viewport(layout = layout))

    for (i in seq_along(plot_list)) {
      row_number <- ((i - 1L) %/% 2L) + 1L
      column_number <- ((i - 1L) %% 2L) + 1L
      print(
        plot_list[[i]],
        vp = grid::viewport(
          layout.pos.row = row_number,
          layout.pos.col = column_number
        )
      )
    }
  }

  grDevices::pdf(
    file = output_pdf,
    width = 15,
    height = 13,
    useDingbats = FALSE
  )
  draw_plots()
  grDevices::dev.off()

  grDevices::png(
    filename = output_png,
    width = 4500,
    height = 3900,
    res = 300,
    bg = "white"
  )
  draw_plots()
  grDevices::dev.off()
}


make_pheatmap_gradient <- function(number_of_colours = 200,
                                   colour_anchor = 0.80) {
  if (number_of_colours < 10L) {
    stop("number_of_colours must be at least 10.")
  }
  if (colour_anchor <= 0 || colour_anchor >= 1) {
    stop("colour_anchor must be between 0 and 1.")
  }

  gradient_specification <- get_correlation_gradient_spec(colour_anchor)

  # Use the same colour anchors and positions as the ggplot2 heatmaps, so all
  # heatmap families share one continuous colour mapping.
  palette_function <- scales::gradient_n_pal(
    colours = gradient_specification$colours,
    values = gradient_specification$positions,
    space = "Lab"
  )
  palette_function(seq(0, 1, length.out = number_of_colours))
}


save_rna_clustered_heatmap <- function(correlation_matrix,
                                       library_label,
                                       normalization_scope_title,
                                       transform_label,
                                       output_pdf,
                                       colour_anchor = 0.80) {
  sample_names <- colnames(correlation_matrix)
  cell_line <- ifelse(
    grepl("^PK15_", sample_names, ignore.case = TRUE),
    "PK15",
    ifelse(
      grepl("^PSC_", sample_names, ignore.case = TRUE),
      "PSC",
      NA_character_
    )
  )
  if (anyNA(cell_line)) {
    stop(
      "Unable to assign CellLine annotation to sample(s): ",
      paste(sample_names[is.na(cell_line)], collapse = ", ")
    )
  }

  annotation_data <- data.frame(
    CellLine = factor(cell_line, levels = c("PK15", "PSC")),
    row.names = sample_names,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  annotation_colours <- list(
    CellLine = c(PK15 = "#2C7FB8", PSC = "#F16913")
  )

  correlation_colours <- make_pheatmap_gradient(
    number_of_colours = 200,
    colour_anchor = colour_anchor
  )
  correlation_breaks <- seq(
    from = 0,
    to = 1,
    length.out = length(correlation_colours) + 1L
  )

  # Correlations below zero, if present, receive the lowest light-blue colour;
  # the printed values and clustering still use the original correlations.
  colour_matrix <- pmin(pmax(correlation_matrix, 0), 1)
  dimnames(colour_matrix) <- dimnames(correlation_matrix)
  displayed_numbers <- matrix(
    sprintf(CORRELATION_LABEL_FORMAT_3, correlation_matrix),
    nrow = nrow(correlation_matrix),
    ncol = ncol(correlation_matrix),
    dimnames = dimnames(correlation_matrix)
  )

  clustering_distance_matrix <- 1 - correlation_matrix
  clustering_distance_matrix[clustering_distance_matrix < 0] <- 0
  diag(clustering_distance_matrix) <- 0
  clustering_distance <- stats::as.dist(clustering_distance_matrix)

  pheatmap::pheatmap(
    mat = colour_matrix,
    color = correlation_colours,
    breaks = correlation_breaks,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = clustering_distance,
    clustering_distance_cols = clustering_distance,
    clustering_method = "complete",
    annotation_row = annotation_data,
    annotation_col = annotation_data,
    annotation_colors = annotation_colours,
    display_numbers = displayed_numbers,
    number_color = "#333333",
    fontsize = 10,
    fontsize_number = 7.5,
    border_color = "white",
    angle_col = 45,
    treeheight_row = 45,
    treeheight_col = 45,
    legend = TRUE,
    legend_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
    legend_labels = format(
      c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
      trim = TRUE,
      nsmall = 1
    ),
    main = paste0(
      library_label, " library: RNA replicate correlation\n",
      normalization_scope_title, "; ", transform_label
    ),
    filename = output_pdf,
    width = 8.8,
    height = 7.8,
    silent = TRUE
  )
}


make_rna_only_clustered_heatmap <- function(library_label,
                                            pk15_normalized_counts,
                                            psc_normalized_counts,
                                            pk15_rna_columns,
                                            psc_rna_columns,
                                            output_directory,
                                            colour_anchor = 0.80) {
  missing_pk15_columns <- setdiff(
    pk15_rna_columns,
    colnames(pk15_normalized_counts)
  )
  missing_psc_columns <- setdiff(
    psc_rna_columns,
    colnames(psc_normalized_counts)
  )
  if (length(missing_pk15_columns) > 0L) {
    stop(
      "Missing PK15 RNA columns for ", library_label, ": ",
      paste(missing_pk15_columns, collapse = ", ")
    )
  }
  if (length(missing_psc_columns) > 0L) {
    stop(
      "Missing PSC RNA columns for ", library_label, ": ",
      paste(missing_psc_columns, collapse = ", ")
    )
  }

  # These RNA columns originate from two independent DESeq2/SS fits. Match by
  # Oligo ID before concatenating them so every sample uses a common feature
  # universe. This output is retained as a sensitivity analysis.
  common_oligos <- intersect(
    rownames(pk15_normalized_counts),
    rownames(psc_normalized_counts)
  )
  if (length(common_oligos) < 2L) {
    stop(
      "Fewer than two common Oligos were found between PK15 and PSC for ",
      library_label
    )
  }

  combined_rna_counts <- cbind(
    pk15_normalized_counts[
      common_oligos,
      pk15_rna_columns,
      drop = FALSE
    ],
    psc_normalized_counts[
      common_oligos,
      psc_rna_columns,
      drop = FALSE
    ]
  )
  if (anyDuplicated(colnames(combined_rna_counts))) {
    stop("Duplicated RNA sample names in the combined ", library_label, " data.")
  }

  summary_rows <- list()
  for (transform_name in names(CORRELATION_TRANSFORMS)) {
    transform_specification <- CORRELATION_TRANSFORMS[[transform_name]]
    result <- calculate_log2_correlation(
      normalized_counts = combined_rna_counts,
      transform_specification = transform_specification,
      context_label = paste0(
        library_label,
        " RNA from separate PK15/PSC normalization"
      )
    )

    matrix_output <- file.path(
      output_directory,
      paste0(
        library_label,
        "_RNA_from_separate_celltype_normalization_Pearson_correlation_",
        result$transform_key,
        ".tsv"
      )
    )
    write_square_matrix(result$correlation_matrix, matrix_output)

    if (WRITE_LOG2_MATRICES) {
      transformed_output <- file.path(
        output_directory,
        paste0(
          library_label,
          "_RNA_from_separate_celltype_normalization_",
          result$transform_key,
          "_counts.tsv.gz"
        )
      )
      write_matrix_with_oligo_column(
        result$transformed_counts,
        transformed_output
      )
    }

    pdf_output <- file.path(
      output_directory,
      paste0(
        library_label,
        "_RNA_from_separate_celltype_normalization_correlation_",
        "clustered_heatmap_",
        result$transform_key,
        "_gradient.pdf"
      )
    )
    save_rna_clustered_heatmap(
      correlation_matrix = result$correlation_matrix,
      library_label = library_label,
      normalization_scope_title = "separate PK15 and PSC normalization",
      transform_label = result$transform_label,
      output_pdf = pdf_output,
      colour_anchor = colour_anchor
    )

    message(
      library_label, " separate-normalization RNA heatmap (",
      result$transform_key, ") used ",
      format(result$oligos_used, big.mark = ","), " of ",
      format(result$total_oligos, big.mark = ","), " common Oligos."
    )
    message("RNA-only clustered heatmap: ", pdf_output)

    summary_rows[[transform_name]] <- data.frame(
      Library = library_label,
      Normalization_scope = "Separate_PK15_and_PSC",
      Transformation = result$transform_label,
      Transform_file_tag = result$transform_key,
      PK15_oligos = nrow(pk15_normalized_counts),
      PSC_oligos = nrow(psc_normalized_counts),
      Common_oligos_available = length(common_oligos),
      Oligos_used_for_correlation = result$oligos_used,
      Oligos_excluded_due_to_zero = result$oligos_excluded_due_to_zero,
      PK15_only_oligos = nrow(pk15_normalized_counts) - length(common_oligos),
      PSC_only_oligos = nrow(psc_normalized_counts) - length(common_oligos),
      Number_of_PK15_RNA_samples = length(pk15_rna_columns),
      Number_of_PSC_RNA_samples = length(psc_rna_columns),
      Correlation_matrix_file = matrix_output,
      Heatmap_PDF = pdf_output,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  summary_table <- do.call(rbind, summary_rows)
  rownames(summary_table) <- NULL
  summary_table
}


make_library_all_sample_heatmap <- function(library_label,
                                             pk15_raw_oligo_counts,
                                             psc_raw_oligo_counts,
                                             pk15_dna_columns,
                                             psc_dna_columns,
                                             pk15_rna_columns,
                                             psc_rna_columns,
                                             output_directory,
                                             colour_anchor = 0.80) {
  if (length(pk15_dna_columns) != 3L || length(psc_dna_columns) != 3L) {
    stop(
      library_label,
      " requires exactly three plasmid replicates in each source file. ",
      "Detected PK15-file plasmids: ", length(pk15_dna_columns),
      "; PSC-file plasmids: ", length(psc_dna_columns), "."
    )
  }
  if (length(pk15_rna_columns) != 3L || length(psc_rna_columns) != 3L) {
    stop(
      library_label,
      " requires exactly three RNA replicates per cell line. ",
      "Detected PK15 RNA: ", length(pk15_rna_columns),
      "; PSC RNA: ", length(psc_rna_columns), "."
    )
  }

  if (!setequal(pk15_dna_columns, psc_dna_columns)) {
    stop(
      "The plasmid sample names differ between the PK15 and PSC files for ",
      library_label, ".\nPK15 file: ",
      paste(pk15_dna_columns, collapse = ", "),
      "\nPSC file: ", paste(psc_dna_columns, collapse = ", ")
    )
  }

  missing_pk15_columns <- setdiff(
    c(pk15_dna_columns, pk15_rna_columns),
    colnames(pk15_raw_oligo_counts)
  )
  missing_psc_columns <- setdiff(
    c(psc_dna_columns, psc_rna_columns),
    colnames(psc_raw_oligo_counts)
  )
  if (length(missing_pk15_columns) > 0L) {
    stop(
      "Missing columns in the ", library_label, " PK15 matrix: ",
      paste(missing_pk15_columns, collapse = ", ")
    )
  }
  if (length(missing_psc_columns) > 0L) {
    stop(
      "Missing columns in the ", library_label, " PSC matrix: ",
      paste(missing_psc_columns, collapse = ", ")
    )
  }

  # PK15 and PSC must be compared using the same Oligos. Keep the PK15 row
  # order but retain only constructs also present in the PSC count matrix.
  common_oligos <- intersect(
    rownames(pk15_raw_oligo_counts),
    rownames(psc_raw_oligo_counts)
  )
  if (length(common_oligos) < 2L) {
    stop(
      "Fewer than two common Oligos were found between PK15 and PSC for ",
      library_label, "."
    )
  }

  # Both cell-type count files contain the same physical plasmid libraries.
  # Confirm their oligo-level counts are identical before retaining one copy;
  # silently choosing between discordant plasmid tables would be arbitrary.
  plasmid_columns <- pk15_dna_columns
  pk15_plasmid_counts <- as.matrix(
    pk15_raw_oligo_counts[
      common_oligos,
      plasmid_columns,
      drop = FALSE
    ]
  )
  psc_plasmid_counts <- as.matrix(
    psc_raw_oligo_counts[
      common_oligos,
      plasmid_columns,
      drop = FALSE
    ]
  )

  plasmid_difference <- abs(pk15_plasmid_counts - psc_plasmid_counts)
  if (any(plasmid_difference != 0)) {
    mismatch_rows <- rowSums(plasmid_difference != 0) > 0
    stop(
      "Plasmid oligo counts are not identical between the PK15 and PSC ",
      "source files for ", library_label, ". Found ", sum(mismatch_rows),
      " common Oligos with at least one mismatch; maximum absolute count ",
      "difference = ", max(plasmid_difference),
      ". Please confirm both files were generated from the same plasmid ",
      "sequencing data before making the combined heatmap."
    )
  }

  combined_raw_counts <- cbind(
    pk15_plasmid_counts,
    as.matrix(
      pk15_raw_oligo_counts[
        common_oligos,
        pk15_rna_columns,
        drop = FALSE
      ]
    ),
    as.matrix(
      psc_raw_oligo_counts[
        common_oligos,
        psc_rna_columns,
        drop = FALSE
      ]
    )
  )
  if (max(combined_raw_counts) > .Machine$integer.max) {
    stop(
      "A combined raw count exceeds R's integer limit for ",
      library_label, "."
    )
  }
  storage.mode(combined_raw_counts) <- "integer"

  expected_sample_order <- c(
    plasmid_columns,
    pk15_rna_columns,
    psc_rna_columns
  )
  if (!identical(colnames(combined_raw_counts), expected_sample_order)) {
    stop("Unexpected sample order in the combined ", library_label, " matrix.")
  }
  if (anyDuplicated(colnames(combined_raw_counts))) {
    stop("Duplicated sample names in the combined ", library_label, " matrix.")
  }
  if (anyNA(combined_raw_counts) || any(combined_raw_counts < 0)) {
    stop("Invalid raw counts in the combined ", library_label, " matrix.")
  }

  # Normalize the nine samples together, rather than mixing normalized values
  # from two independent DESeq2 fits. tagNorm() fits one DESeq2 model and then
  # applies a separate summit-shift correction to PK15 and PSC relative to the
  # same three plasmid samples.
  combined_condition <- c(
    rep("DNA", length(plasmid_columns)),
    rep("PK15", length(pk15_rna_columns)),
    rep("PSC", length(psc_rna_columns))
  )
  combined_condition_data <- data.frame(
    condition = combined_condition,
    row.names = expected_sample_order,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  combined_attributes <- data.frame(
    ID = common_oligos,
    ctrl_exp = "MPRA_PROJ",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(combined_attributes) <- combined_attributes$ID

  message(
    "Running joint Plasmid + PK15 + PSC DESeq2/SS normalization for ",
    library_label, " ..."
  )
  combined_normalization <- tagNorm(
    countsData = combined_raw_counts,
    conditionData = combined_condition_data,
    attributesData = combined_attributes,
    exclList = c(),
    method = "ss",
    negCtrlName = "__NO_NEGATIVE_CONTROLS__",
    upDisp = TRUE,
    prior = FALSE
  )
  if (!is.list(combined_normalization) ||
      length(combined_normalization) < 2L) {
    stop("Unexpected tagNorm() return value for combined ", library_label, ".")
  }

  combined_dds_ss <- combined_normalization[[1]]
  combined_dds_before_ss <- combined_normalization[[2]]
  combined_normalized_counts <- as.matrix(
    DESeq2::counts(combined_dds_ss, normalized = TRUE)
  )
  combined_normalized_counts <- combined_normalized_counts[
    common_oligos,
    expected_sample_order,
    drop = FALSE
  ]

  if (any(!is.finite(combined_normalized_counts)) ||
      any(combined_normalized_counts < 0)) {
    stop("Invalid joint normalized counts for ", library_label, ".")
  }

  normalized_output <- file.path(
    output_directory,
    paste0(
      library_label,
      "_Plasmid_PK15_PSC_oligo_DESeq2_ss_normalized_counts.tsv.gz"
    )
  )
  write_matrix_with_oligo_column(
    combined_normalized_counts,
    normalized_output
  )

  joint_rna_columns <- c(pk15_rna_columns, psc_rna_columns)
  joint_rna_normalized_counts <- combined_normalized_counts[
    ,
    joint_rna_columns,
    drop = FALSE
  ]
  joint_rna_normalized_output <- file.path(
    output_directory,
    paste0(
      library_label,
      "_RNA_extracted_from_joint9_DESeq2_ss_normalized_counts.tsv.gz"
    )
  )
  write_matrix_with_oligo_column(
    joint_rna_normalized_counts,
    joint_rna_normalized_output
  )

  size_factors_before <- stats::setNames(
    as.numeric(DESeq2::sizeFactors(combined_dds_before_ss)),
    colnames(combined_dds_before_ss)
  )
  size_factors_after <- stats::setNames(
    as.numeric(DESeq2::sizeFactors(combined_dds_ss)),
    colnames(combined_dds_ss)
  )
  size_factor_output <- file.path(
    output_directory,
    paste0(
      library_label,
      "_Plasmid_PK15_PSC_DESeq2_ss_size_factors.tsv"
    )
  )
  write.table(
    data.frame(
      Library = library_label,
      Sample = expected_sample_order,
      Group = c(
        rep("Plasmid", length(plasmid_columns)),
        rep("PK15", length(pk15_rna_columns)),
        rep("PSC", length(psc_rna_columns))
      ),
      DESeq2_size_factor_before_SS = as.numeric(
        size_factors_before[expected_sample_order]
      ),
      Final_size_factor_after_SS = as.numeric(
        size_factors_after[expected_sample_order]
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    file = size_factor_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  group_names <- c(
    rep("Plasmid", length(plasmid_columns)),
    rep("PK15", length(pk15_rna_columns)),
    rep("PSC", length(psc_rna_columns))
  )
  annotation_data <- data.frame(
    Sample_group = factor(
      group_names,
      levels = c("Plasmid", "PK15", "PSC")
    ),
    row.names = expected_sample_order,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  annotation_colours <- list(
    Sample_group = c(
      Plasmid = "#4D4D4D",
      PK15 = "#2C7FB8",
      PSC = "#F16913"
    )
  )

  correlation_colours <- make_pheatmap_gradient(
    number_of_colours = 200,
    colour_anchor = colour_anchor
  )
  correlation_breaks <- seq(
    from = 0,
    to = 1,
    length.out = length(correlation_colours) + 1L
  )
  group_gaps <- c(
    length(plasmid_columns),
    length(plasmid_columns) + length(pk15_rna_columns)
  )

  all_sample_summary_rows <- list()
  joint_rna_summary_rows <- list()

  for (transform_name in names(CORRELATION_TRANSFORMS)) {
    transform_specification <- CORRELATION_TRANSFORMS[[transform_name]]

    # Calculate the 9 x 9 matrix first. For the no-+1 version, positivity is
    # required across all nine samples. The same transformed Oligo matrix is
    # then used for the six-RNA submatrix, guaranteeing that the RNA-only
    # correlations exactly match the RNA block of the corresponding 9 x 9
    # correlation matrix (apart from printed rounding).
    all_sample_result <- calculate_log2_correlation(
      normalized_counts = combined_normalized_counts,
      transform_specification = transform_specification,
      context_label = paste0(
        library_label,
        " jointly normalized Plasmid/PK15/PSC samples"
      )
    )
    correlation_matrix <- all_sample_result$correlation_matrix

    correlation_output <- file.path(
      output_directory,
      paste0(
        library_label,
        "_Plasmid_PK15_PSC_joint9_Pearson_correlation_",
        all_sample_result$transform_key,
        ".tsv"
      )
    )
    write_square_matrix(correlation_matrix, correlation_output)

    if (WRITE_LOG2_MATRICES) {
      transformed_output <- file.path(
        output_directory,
        paste0(
          library_label,
          "_Plasmid_PK15_PSC_joint9_",
          all_sample_result$transform_key,
          "_counts.tsv.gz"
        )
      )
      write_matrix_with_oligo_column(
        all_sample_result$transformed_counts,
        transformed_output
      )
    }

    colour_matrix <- pmin(pmax(correlation_matrix, 0), 1)
    dimnames(colour_matrix) <- dimnames(correlation_matrix)
    displayed_numbers <- matrix(
      sprintf(CORRELATION_LABEL_FORMAT_3, correlation_matrix),
      nrow = nrow(correlation_matrix),
      ncol = ncol(correlation_matrix),
      dimnames = dimnames(correlation_matrix)
    )
    heatmap_output <- file.path(
      output_directory,
      paste0(
        library_label,
        "_Plasmid_PK15_PSC_joint9_correlation_heatmap_",
        all_sample_result$transform_key,
        "_gradient.pdf"
      )
    )

    # Fixed article-like order: Plasmid -> PK15 -> PSC. White gaps divide the
    # three 3 x 3 within-group blocks; off-diagonal cells show cross-group r.
    pheatmap::pheatmap(
      mat = colour_matrix,
      color = correlation_colours,
      breaks = correlation_breaks,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      gaps_row = group_gaps,
      gaps_col = group_gaps,
      annotation_row = annotation_data,
      annotation_col = annotation_data,
      annotation_colors = annotation_colours,
      display_numbers = displayed_numbers,
      number_color = "#333333",
      fontsize = 10,
      fontsize_number = 7.5,
      border_color = "white",
      angle_col = 45,
      legend = TRUE,
      legend_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
      legend_labels = format(
        c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1),
        trim = TRUE,
        nsmall = 1
      ),
      main = paste0(
        library_label,
        " library: jointly normalized Plasmid, PK15 RNA and PSC RNA\n",
        all_sample_result$transform_label
      ),
      filename = heatmap_output,
      width = 9.8,
      height = 9.2,
      silent = TRUE
    )

    # Extract the six RNA columns from the already jointly normalized and
    # transformed nine-sample matrix. No second normalization is performed.
    joint_rna_transformed_counts <- all_sample_result$transformed_counts[
      ,
      joint_rna_columns,
      drop = FALSE
    ]
    joint_rna_correlation <- stats::cor(
      joint_rna_transformed_counts,
      method = "pearson",
      use = "everything"
    )
    expected_rna_submatrix <- correlation_matrix[
      joint_rna_columns,
      joint_rna_columns,
      drop = FALSE
    ]
    if (!isTRUE(all.equal(
      joint_rna_correlation,
      expected_rna_submatrix,
      tolerance = 1e-12,
      check.attributes = TRUE
    ))) {
      stop(
        "Internal consistency check failed: the ", library_label,
        " joint-normalization RNA correlation matrix is not identical to ",
        "the RNA submatrix of the 9-sample correlation matrix."
      )
    }

    joint_rna_correlation_output <- file.path(
      output_directory,
      paste0(
        library_label,
        "_RNA_from_joint9_normalization_Pearson_correlation_",
        all_sample_result$transform_key,
        ".tsv"
      )
    )
    write_square_matrix(
      joint_rna_correlation,
      joint_rna_correlation_output
    )

    joint_rna_heatmap_output <- file.path(
      output_directory,
      paste0(
        library_label,
        "_RNA_from_joint9_normalization_correlation_clustered_heatmap_",
        all_sample_result$transform_key,
        "_gradient.pdf"
      )
    )
    save_rna_clustered_heatmap(
      correlation_matrix = joint_rna_correlation,
      library_label = library_label,
      normalization_scope_title = "RNA extracted from joint 9-sample normalization",
      transform_label = all_sample_result$transform_label,
      output_pdf = joint_rna_heatmap_output,
      colour_anchor = colour_anchor
    )

    message(
      library_label, " joint 9-sample heatmap (",
      all_sample_result$transform_key, ") used ",
      format(all_sample_result$oligos_used, big.mark = ","), " of ",
      format(all_sample_result$total_oligos, big.mark = ","),
      " common Oligos."
    )
    message("Library-level all-sample heatmap: ", heatmap_output)
    message("RNA extracted from joint-9 heatmap: ", joint_rna_heatmap_output)

    all_sample_summary_rows[[transform_name]] <- data.frame(
      Library = library_label,
      Normalization_scope = "Joint_9_samples",
      Heatmap_scope = "Plasmid_PK15_PSC_9x9",
      Transformation = all_sample_result$transform_label,
      Transform_file_tag = all_sample_result$transform_key,
      PK15_oligos = nrow(pk15_raw_oligo_counts),
      PSC_oligos = nrow(psc_raw_oligo_counts),
      Common_oligos_available = length(common_oligos),
      Oligos_used_for_correlation = all_sample_result$oligos_used,
      Oligos_excluded_due_to_zero =
        all_sample_result$oligos_excluded_due_to_zero,
      Zero_filter_scope = ifelse(
        isTRUE(transform_specification$add_one),
        "Not_applicable",
        "All_9_joint_samples"
      ),
      PK15_only_oligos =
        nrow(pk15_raw_oligo_counts) - length(common_oligos),
      PSC_only_oligos =
        nrow(psc_raw_oligo_counts) - length(common_oligos),
      Number_of_plasmid_samples = length(plasmid_columns),
      Number_of_PK15_RNA_samples = length(pk15_rna_columns),
      Number_of_PSC_RNA_samples = length(psc_rna_columns),
      Joint_normalized_counts_file = normalized_output,
      Correlation_matrix_file = correlation_output,
      Size_factor_file = size_factor_output,
      Heatmap_PDF = heatmap_output,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    joint_rna_summary_rows[[transform_name]] <- data.frame(
      Library = library_label,
      Normalization_scope = "Joint_9_samples_then_extract_RNA",
      Heatmap_scope = "PK15_PSC_RNA_6x6",
      Transformation = all_sample_result$transform_label,
      Transform_file_tag = all_sample_result$transform_key,
      Common_oligos_available = length(common_oligos),
      Oligos_used_for_correlation = all_sample_result$oligos_used,
      Oligos_excluded_due_to_zero =
        all_sample_result$oligos_excluded_due_to_zero,
      Zero_filter_scope = ifelse(
        isTRUE(transform_specification$add_one),
        "Not_applicable",
        "All_9_joint_samples"
      ),
      Number_of_PK15_RNA_samples = length(pk15_rna_columns),
      Number_of_PSC_RNA_samples = length(psc_rna_columns),
      Joint_9sample_normalized_counts_file = normalized_output,
      Extracted_6RNA_normalized_counts_file = joint_rna_normalized_output,
      Correlation_matrix_file = joint_rna_correlation_output,
      Heatmap_PDF = joint_rna_heatmap_output,
      Exact_RNA_submatrix_of_9x9 = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  all_sample_summary <- do.call(rbind, all_sample_summary_rows)
  rownames(all_sample_summary) <- NULL
  joint_rna_summary <- do.call(rbind, joint_rna_summary_rows)
  rownames(joint_rna_summary) <- NULL

  list(
    all_sample_summary = all_sample_summary,
    joint_rna_summary = joint_rna_summary
  )
}


off_diagonal_values <- function(correlation_matrix) {
  correlation_matrix[upper.tri(correlation_matrix, diag = FALSE)]
}


safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}


safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}


# ----------------------------------------------------------------------------
# 4. Run all four datasets
# ----------------------------------------------------------------------------

all_heatmaps <- stats::setNames(
  lapply(names(CORRELATION_TRANSFORMS), function(x) list()),
  names(CORRELATION_TRANSFORMS)
)
all_size_factors <- list()
all_summaries <- list()
all_normalized_counts <- list()
all_rna_columns <- list()
all_raw_oligo_counts <- list()
all_dna_columns <- list()

for (dataset_name in names(DATASETS)) {
  specification <- DATASETS[[dataset_name]]

  message("\n", paste(rep("=", 78), collapse = ""))
  message("Processing dataset: ", dataset_name)
  message(paste(rep("=", 78), collapse = ""))

  imported <- read_barcode_count_file(
    path = specification$input,
    cell_type = specification$cell_type
  )

  barcode_data <- imported$data
  sample_columns <- imported$sample_columns
  dna_columns <- imported$dna_columns
  rna_columns <- imported$rna_columns

  message("Barcode-level rows: ", format(nrow(barcode_data), big.mark = ","))
  message("Detected plasmid samples: ", paste(dna_columns, collapse = ", "))
  message(
    "Detected ", specification$cell_type, " RNA samples: ",
    paste(rna_columns, collapse = ", ")
  )

  # Official MPRASuite barcode-to-oligo aggregation.
  message("Aggregating barcode counts with official oligoIsolate() ...")
  oligo_counts <- run_official_oligo_isolate(
    barcode_data = barcode_data,
    dataset_name = dataset_name
  )

  # Restore the exact original sample-column order.
  missing_after_aggregation <- setdiff(sample_columns, colnames(oligo_counts))
  if (length(missing_after_aggregation) > 0L) {
    stop(
      "Sample columns missing after oligoIsolate() for ", dataset_name, ": ",
      paste(missing_after_aggregation, collapse = ", ")
    )
  }
  oligo_counts <- oligo_counts[, sample_columns, drop = FALSE]

  if (anyDuplicated(rownames(oligo_counts))) {
    stop("Duplicated Oligo row names remain after aggregation: ", dataset_name)
  }

  if (anyNA(oligo_counts)) {
    stop("NA values were generated after oligo aggregation: ", dataset_name)
  }

  sample_totals <- colSums(as.matrix(oligo_counts))
  if (any(sample_totals <= 0)) {
    stop(
      "At least one sample has a total oligo count of zero in ", dataset_name,
      ": ", paste(names(sample_totals)[sample_totals <= 0], collapse = ", ")
    )
  }

  message(
    "Oligos after barcode summation: ",
    format(nrow(oligo_counts), big.mark = ",")
  )

  if (WRITE_RAW_OLIGO_COUNTS) {
    raw_oligo_output <- file.path(
      OUTPUT_DIR,
      paste0(dataset_name, "_oligo_raw_counts.tsv.gz")
    )
    write_matrix_with_oligo_column(oligo_counts, raw_oligo_output)
  }

  # Sample metadata for the official MPRASuite functions. Plasmids are DNA;
  # PK15/PSC columns are the corresponding RNA condition.
  condition_vector <- ifelse(
    sample_columns %in% dna_columns,
    "DNA",
    specification$cell_type
  )
  condition_data <- data.frame(
    condition = condition_vector,
    row.names = sample_columns,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # tagNorm(method="ss") needs ID and ctrl_exp fields. The standard ss method
  # uses all oligos and does not require negative-control annotations, so a
  # minimal attributes table is sufficient for this normalization-only task.
  normalization_attributes <- data.frame(
    ID = rownames(oligo_counts),
    ctrl_exp = "MPRA_PROJ",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(normalization_attributes) <- normalization_attributes$ID

  message("Running official DESeq2 + summit-shift normalization ...")
  normalization_result <- tagNorm(
    countsData = oligo_counts,
    conditionData = condition_data,
    attributesData = normalization_attributes,
    exclList = c(),
    method = "ss",
    negCtrlName = "__NO_NEGATIVE_CONTROLS__",
    upDisp = TRUE,
    prior = FALSE
  )

  if (!is.list(normalization_result) || length(normalization_result) < 2L) {
    stop(
      "Unexpected return value from tagNorm() for ", dataset_name,
      ". The local MPRAmodel.R may not match the official interface."
    )
  }

  dds_summit_shift <- normalization_result[[1]]
  dds_before_summit_shift <- normalization_result[[2]]

  normalized_counts <- as.matrix(
    DESeq2::counts(dds_summit_shift, normalized = TRUE)
  )
  normalized_counts <- normalized_counts[, sample_columns, drop = FALSE]

  if (!identical(rownames(normalized_counts), rownames(oligo_counts))) {
    stop("Oligo order changed during normalization: ", dataset_name)
  }

  if (any(!is.finite(normalized_counts)) || any(normalized_counts < 0)) {
    stop("Invalid normalized counts were generated for: ", dataset_name)
  }

  normalized_output <- file.path(
    OUTPUT_DIR,
    paste0(dataset_name, "_oligo_DESeq2_ss_normalized_counts.tsv.gz")
  )
  write_matrix_with_oligo_column(normalized_counts, normalized_output)

  # No activity/count-QC/emVAR filtering is applied before correlation. The
  # +1 branch reproduces the article's public QC-code transformation. The
  # no-+1 branch is an explicitly labelled supplementary comparison.
  dataset_summary_rows <- list()
  for (transform_name in names(CORRELATION_TRANSFORMS)) {
    transform_specification <- CORRELATION_TRANSFORMS[[transform_name]]
    correlation_result <- calculate_log2_correlation(
      normalized_counts = normalized_counts,
      transform_specification = transform_specification,
      context_label = dataset_name
    )

    if (WRITE_LOG2_MATRICES) {
      log2_output <- file.path(
        OUTPUT_DIR,
        paste0(
          dataset_name,
          "_oligo_",
          correlation_result$transform_key,
          "_counts.tsv.gz"
        )
      )
      write_matrix_with_oligo_column(
        correlation_result$transformed_counts,
        log2_output
      )
    }

    correlation_matrix <- correlation_result$correlation_matrix
    correlation_output <- file.path(
      OUTPUT_DIR,
      paste0(
        dataset_name,
        "_Pearson_correlation_",
        correlation_result$transform_key,
        ".tsv"
      )
    )
    write_square_matrix(correlation_matrix, correlation_output)

    heatmap_plot <- make_correlation_heatmap(
      correlation_matrix = correlation_matrix,
      plot_title = specification$title,
      transform_label = correlation_result$transform_label,
      colour_anchor = CORRELATION_COLOR_ANCHOR,
      label_format = "%.3f"
    )

    colour_anchor_text <- sprintf("%.2f", CORRELATION_COLOR_ANCHOR)
    heatmap_pdf <- file.path(
      OUTPUT_DIR,
      paste0(
        dataset_name,
        "_correlation_heatmap_",
        correlation_result$transform_key,
        "_r",
        colour_anchor_text,
        ".pdf"
      )
    )
    heatmap_png <- file.path(
      OUTPUT_DIR,
      paste0(
        dataset_name,
        "_correlation_heatmap_",
        correlation_result$transform_key,
        "_r",
        colour_anchor_text,
        ".png"
      )
    )

    ggplot2::ggsave(
      filename = heatmap_pdf,
      plot = heatmap_plot,
      width = 7.4,
      height = 6.6,
      units = "in",
      device = "pdf",
      useDingbats = FALSE,
      bg = "white"
    )
    ggplot2::ggsave(
      filename = heatmap_png,
      plot = heatmap_plot,
      width = 7.4,
      height = 6.6,
      units = "in",
      dpi = 400,
      bg = "white"
    )

    plasmid_correlation <- correlation_matrix[
      dna_columns,
      dna_columns,
      drop = FALSE
    ]
    rna_correlation <- correlation_matrix[
      rna_columns,
      rna_columns,
      drop = FALSE
    ]
    cross_correlation <- correlation_matrix[
      dna_columns,
      rna_columns,
      drop = FALSE
    ]

    dataset_summary_rows[[transform_name]] <- data.frame(
      Dataset = dataset_name,
      Transformation = correlation_result$transform_label,
      Transform_file_tag = correlation_result$transform_key,
      Input_file = specification$input,
      Barcode_rows = nrow(barcode_data),
      Unique_barcodes = length(unique(barcode_data$Barcode)),
      Unique_oligos_before_aggregation = length(unique(barcode_data$Oligo)),
      Oligo_rows_after_aggregation = nrow(oligo_counts),
      All_zero_oligo_rows = sum(rowSums(oligo_counts) == 0),
      Oligos_used_for_correlation = correlation_result$oligos_used,
      Oligos_excluded_due_to_zero =
        correlation_result$oligos_excluded_due_to_zero,
      Number_of_samples = length(sample_columns),
      Number_of_plasmid_samples = length(dna_columns),
      Number_of_RNA_samples = length(rna_columns),
      Minimum_all_off_diagonal_r = safe_min(
        off_diagonal_values(correlation_matrix)
      ),
      Maximum_all_off_diagonal_r = safe_max(
        off_diagonal_values(correlation_matrix)
      ),
      Minimum_within_plasmid_r = safe_min(
        off_diagonal_values(plasmid_correlation)
      ),
      Minimum_within_RNA_r = safe_min(
        off_diagonal_values(rna_correlation)
      ),
      Minimum_plasmid_RNA_r = safe_min(as.vector(cross_correlation)),
      Maximum_plasmid_RNA_r = safe_max(as.vector(cross_correlation)),
      Normalized_counts_file = normalized_output,
      Correlation_matrix_file = correlation_output,
      Heatmap_PDF = heatmap_pdf,
      Heatmap_PNG = heatmap_png,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    all_heatmaps[[transform_name]][[dataset_name]] <- heatmap_plot
    message(
      "Heatmap (", correlation_result$transform_key, "): ", heatmap_pdf,
      "; Oligos used = ",
      format(correlation_result$oligos_used, big.mark = ",")
    )
  }

  size_factors_before <- stats::setNames(
    as.numeric(DESeq2::sizeFactors(dds_before_summit_shift)),
    colnames(dds_before_summit_shift)
  )
  size_factors_after <- stats::setNames(
    as.numeric(DESeq2::sizeFactors(dds_summit_shift)),
    colnames(dds_summit_shift)
  )

  size_factor_table <- data.frame(
    Dataset = dataset_name,
    Sample = sample_columns,
    Condition = condition_vector,
    DESeq2_size_factor_before_SS = as.numeric(
      size_factors_before[sample_columns]
    ),
    Final_size_factor_after_SS = as.numeric(
      size_factors_after[sample_columns]
    ),
    SS_multiplier = as.numeric(
      size_factors_after[sample_columns] /
        size_factors_before[sample_columns]
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  all_size_factors[[dataset_name]] <- size_factor_table
  all_summaries[[dataset_name]] <- do.call(rbind, dataset_summary_rows)
  rownames(all_summaries[[dataset_name]]) <- NULL
  all_normalized_counts[[dataset_name]] <- normalized_counts
  all_rna_columns[[dataset_name]] <- rna_columns
  all_raw_oligo_counts[[dataset_name]] <- oligo_counts
  all_dna_columns[[dataset_name]] <- dna_columns

  message("Finished: ", dataset_name)
  message("Normalized counts: ", normalized_output)

  rm(
    imported,
    barcode_data,
    oligo_counts,
    normalization_attributes,
    normalization_result,
    dds_summit_shift,
    dds_before_summit_shift,
    normalized_counts,
    correlation_matrix,
    dataset_summary_rows
  )
  invisible(gc())
}


# ----------------------------------------------------------------------------
# 5. Combined reports
# ----------------------------------------------------------------------------

# Retain the RNA-only sensitivity analysis based on independently normalized
# PK15 and PSC matrices. Each library now receives both log2 transformations.
rna_heatmap_summaries <- list(
  `3UTR` = make_rna_only_clustered_heatmap(
    library_label = "3UTR",
    pk15_normalized_counts = all_normalized_counts[["3UTR_PK15"]],
    psc_normalized_counts = all_normalized_counts[["3UTR_PSC"]],
    pk15_rna_columns = all_rna_columns[["3UTR_PK15"]],
    psc_rna_columns = all_rna_columns[["3UTR_PSC"]],
    output_directory = OUTPUT_DIR,
    colour_anchor = CORRELATION_COLOR_ANCHOR
  ),
  `CRE` = make_rna_only_clustered_heatmap(
    library_label = "CRE",
    pk15_normalized_counts = all_normalized_counts[["CRE_PK15"]],
    psc_normalized_counts = all_normalized_counts[["CRE_PSC"]],
    pk15_rna_columns = all_rna_columns[["CRE_PK15"]],
    psc_rna_columns = all_rna_columns[["CRE_PSC"]],
    output_directory = OUTPUT_DIR,
    colour_anchor = CORRELATION_COLOR_ANCHOR
  )
)

rna_heatmap_summary_report <- do.call(rbind, rna_heatmap_summaries)
rownames(rna_heatmap_summary_report) <- NULL

write.table(
  rna_heatmap_summary_report,
  file = file.path(
    OUTPUT_DIR,
    "RNA_from_separate_celltype_normalization_heatmap_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# For each library, perform one joint 9-sample normalization. The function
# returns both the fixed-order 9 x 9 heatmaps and the newly requested 6 x 6
# RNA-only heatmaps extracted from that same joint-normalized matrix.
library_joint_normalization_results <- list(
  `3UTR` = make_library_all_sample_heatmap(
    library_label = "3UTR",
    pk15_raw_oligo_counts = all_raw_oligo_counts[["3UTR_PK15"]],
    psc_raw_oligo_counts = all_raw_oligo_counts[["3UTR_PSC"]],
    pk15_dna_columns = all_dna_columns[["3UTR_PK15"]],
    psc_dna_columns = all_dna_columns[["3UTR_PSC"]],
    pk15_rna_columns = all_rna_columns[["3UTR_PK15"]],
    psc_rna_columns = all_rna_columns[["3UTR_PSC"]],
    output_directory = OUTPUT_DIR,
    colour_anchor = CORRELATION_COLOR_ANCHOR
  ),
  `CRE` = make_library_all_sample_heatmap(
    library_label = "CRE",
    pk15_raw_oligo_counts = all_raw_oligo_counts[["CRE_PK15"]],
    psc_raw_oligo_counts = all_raw_oligo_counts[["CRE_PSC"]],
    pk15_dna_columns = all_dna_columns[["CRE_PK15"]],
    psc_dna_columns = all_dna_columns[["CRE_PSC"]],
    pk15_rna_columns = all_rna_columns[["CRE_PK15"]],
    psc_rna_columns = all_rna_columns[["CRE_PSC"]],
    output_directory = OUTPUT_DIR,
    colour_anchor = CORRELATION_COLOR_ANCHOR
  )
)

library_all_sample_heatmap_summary_report <- do.call(
  rbind,
  lapply(
    library_joint_normalization_results,
    `[[`,
    "all_sample_summary"
  )
)
rownames(library_all_sample_heatmap_summary_report) <- NULL

write.table(
  library_all_sample_heatmap_summary_report,
  file = file.path(
    OUTPUT_DIR,
    "Library_joint9_Plasmid_PK15_PSC_heatmap_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

joint_rna_heatmap_summary_report <- do.call(
  rbind,
  lapply(
    library_joint_normalization_results,
    `[[`,
    "joint_rna_summary"
  )
)
rownames(joint_rna_heatmap_summary_report) <- NULL

write.table(
  joint_rna_heatmap_summary_report,
  file = file.path(
    OUTPUT_DIR,
    "RNA_from_joint9_normalization_heatmap_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

size_factor_report <- do.call(rbind, all_size_factors)
rownames(size_factor_report) <- NULL

write.table(
  size_factor_report,
  file = file.path(OUTPUT_DIR, "All_4_datasets_DESeq2_ss_size_factors.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

summary_report <- do.call(rbind, all_summaries)
rownames(summary_report) <- NULL

write.table(
  summary_report,
  file = file.path(OUTPUT_DIR, "All_4_datasets_correlation_run_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

for (transform_name in names(CORRELATION_TRANSFORMS)) {
  transform_tag <- CORRELATION_TRANSFORMS[[transform_name]]$key
  save_combined_heatmaps(
    plot_list = all_heatmaps[[transform_name]],
    output_pdf = file.path(
      OUTPUT_DIR,
      paste0(
        "All_4_datasets_correlation_heatmaps_",
        transform_tag,
        "_r",
        sprintf("%.2f", CORRELATION_COLOR_ANCHOR),
        ".pdf"
      )
    ),
    output_png = file.path(
      OUTPUT_DIR,
      paste0(
        "All_4_datasets_correlation_heatmaps_",
        transform_tag,
        "_r",
        sprintf("%.2f", CORRELATION_COLOR_ANCHOR),
        ".png"
      )
    )
  )
}

capture.output(
  sessionInfo(),
  file = file.path(OUTPUT_DIR, "R_sessionInfo.txt")
)

message("\n", paste(rep("=", 78), collapse = ""))
message(
  "All four datasets, separate-normalization RNA heatmaps, joint 9-sample ",
  "heatmaps, and RNA-from-joint9 heatmaps completed successfully for both ",
  "log2 transformations."
)
message("Output directory: ", OUTPUT_DIR)
message("Heatmap labels: 3 decimals retained; v4 also generates 2-decimal display-ready versions.")
message(paste(rep("=", 78), collapse = ""))
