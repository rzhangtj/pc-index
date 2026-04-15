#!/usr/bin/env Rscript

# Compute PC-index, mean CP10K, and percent-positive for each gene
# from a cell-by-gene CP10K matrix.
#
# Definition used here:
# PC-index = max_k min(sorted_expression_k, 100 * k / n_cells)
#
# This is the continuous version of the prevalence-CP10K intersection.
# For example, a PC-index of 21 means that at least 21% of cells
# express the gene at >= 21 CP10K.
#
# Also supports direct input from Seurat objects via
# pc_index_seurat().
#
# Input:
# - matrix/data.frame: rows = cells, columns = genes
# - values = CP10K
#
# Output:
# - data.frame with one row per gene and columns:
#     gene, n_cells, mean_cp10k, pct_positive, pc_index

compute_pc_index <- function(values) {
  x <- as.numeric(values)
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  x <- sort(x, decreasing = TRUE)
  prevalence_pct <- 100 * seq_along(x) / length(x)
  max(pmin(x, prevalence_pct))
}

summarize_matrix <- function(mat) {
  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }

  if (!is.matrix(mat)) {
    stop("Input must be a matrix or data.frame with rows = cells and columns = genes.")
  }

  n_cells <- nrow(mat)
  genes <- colnames(mat)

  if (is.null(genes)) {
    genes <- paste0("gene_", seq_len(ncol(mat)))
  }

  mean_cp10k <- apply(mat, 2, function(x) mean(as.numeric(x), na.rm = TRUE))
  pct_positive <- apply(mat, 2, function(x) mean(as.numeric(x) > 0, na.rm = TRUE) * 100)
  pc_index <- apply(mat, 2, compute_pc_index)

  out <- data.frame(
    gene = genes,
    n_cells = n_cells,
    mean_cp10k = as.numeric(mean_cp10k),
    pct_positive = as.numeric(pct_positive),
    pc_index = as.numeric(pc_index),
    stringsAsFactors = FALSE
  )

  out <- out[order(out$pc_index, decreasing = TRUE), ]
  rownames(out) <- NULL
  out
}

pc_index_seurat <- function(seurat_obj, assay = NULL, slot = "data") {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("SeuratObject package is required for pc_index_seurat().")
  }

  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(seurat_obj)
  }

  expr <- SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = slot)

  if (inherits(expr, "dgCMatrix")) {
    expr <- as.matrix(expr)
  } else if (!is.matrix(expr)) {
    expr <- as.matrix(expr)
  }

  # Seurat matrices are typically gene-by-cell; transpose to cell-by-gene
  expr <- t(expr)
  summarize_matrix(expr)
}

read_table <- function(path, sep = NULL, row_names = TRUE) {
  if (is.null(sep)) {
    if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) {
      sep <- "\t"
    } else {
      sep <- ","
    }
  }

  df <- read.table(
    path,
    sep = sep,
    header = TRUE,
    row.names = if (row_names) 1 else NULL,
    check.names = FALSE
  )
  as.data.frame(df)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 2) {
    cat("Usage: Rscript pc_index.R input.csv output.csv [sep] [transpose]\n")
    quit(status = 1)
  }

  input <- args[1]
  output <- args[2]
  sep <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NULL
  transpose <- if (length(args) >= 4) tolower(args[4]) %in% c("true", "t", "1", "yes") else FALSE

  df <- read_table(input, sep = sep, row_names = TRUE)

  if (transpose) {
    df <- as.data.frame(t(df), check.names = FALSE)
  }

  summary_df <- summarize_matrix(df)
  write.csv(summary_df, output, row.names = FALSE)

  cat(sprintf("Done. Wrote %d genes to: %s\n", nrow(summary_df), output))
}

if (sys.nframe() == 0) {
  main()
}
