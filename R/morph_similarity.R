#' Pairwise similarity among morphological traits
#'
#' @description
#' Computes pairwise correlations among all morphological trait pairs within
#' one or more trait blocks of a morphometric matrix produced by
#' [morph_matrix_setting()], and exports correlation heatmaps, Excel result
#' tables, and, optionally, annotated scatter plots for user-selected trait
#' pairs.
#'
#' @details
#' `analysis_data` is expected to be the morphometric matrix returned by
#' [morph_matrix_setting()], i.e. a numeric data frame of specimens (rows) by
#' morphological traits (columns), carrying a `"base_cols"` attribute — a
#' named list assigning trait columns to the `"veg"`, `"flo"`, and `"fru"`
#' structural blocks. This attribute is attached automatically by
#' [morph_matrix_setting()]; if `analysis_data` does not carry it, the function
#' stops with an informative error.
#'
#' When `sim_blocks = "all"`, the function runs the analysis on each
#' individual block (`"veg"`, `"flo"`, `"fru"`) as well as on every pairwise
#' and three-way combination of blocks (`"vegflo"`, `"vegfru"`, `"flofru"`,
#' `"vegflofru"`). Alternatively, `sim_blocks` can be a character vector
#' naming a subset of the blocks present in `base_cols` (e.g. `c("veg",
#' "flo")`), in which case only those individual blocks are analyzed.
#'
#' For each block, all pairwise combinations of traits are correlated using
#' [stats::cor.test()] with the method given in `sim_method`. Pairs with
#' fewer than `sim_min_n` complete (non-`NA`) observations are skipped.
#' Resulting p-values are adjusted across all pairs within a block using
#' [stats::p.adjust()] with the method given in `sim_p_adjust`, and pairs
#' with an adjusted p-value below `sim_alpha` are flagged as significant.
#'
#' For each block, the function produces:
#' \itemize{
#'   \item a correlation heatmap (PDF and/or JPEG, per `file_formats`);
#'   \item an Excel spreadsheet of pairwise correlation statistics, sorted by
#'   absolute correlation strength.
#' }
#'
#' A combined Excel spreadsheet across all analyzed blocks is also written.
#'
#' If `sim_pairs` is supplied, an annotated scatter plot (with a linear fit,
#' correlation coefficient, significance, and sample size) is additionally
#' produced for each specified trait pair, regardless of block.
#'
#' @param analysis_data A numeric data frame of morphological measurements,
#' typically the matrix returned by [morph_matrix_setting()], carrying a
#' `"base_cols"` attribute (see Details).
#'
#' @param sim_blocks Character. Either `"all"` (default), to analyze every
#' individual block plus all block combinations, or a character vector naming
#' a subset of blocks in `base_cols` (e.g. `"veg"`, `c("veg", "flo")`).
#'
#' @param sim_method Character string giving the correlation method passed to
#' [stats::cor.test()]. One of `"pearson"` (default), `"spearman"`, or
#' `"kendall"`.
#'
#' @param sim_p_adjust Character string giving the p-value adjustment method
#' passed to [stats::p.adjust()] (e.g. `"BH"`, `"bonferroni"`, `"holm"`,
#' `"none"`). Default is `"BH"`.
#'
#' @param sim_alpha Numeric. Adjusted p-value threshold below which a
#' correlation is flagged as significant. Default is `0.05`.
#'
#' @param sim_min_n Integer. Minimum number of complete pairwise observations
#' required to test a trait pair. Pairs with fewer observations are skipped.
#' Default is `5`.
#'
#' @param sim_pairs Optional list of character vectors of length two, each
#' giving the names of two trait columns in `analysis_data` for which an
#' annotated scatter plot should be produced (e.g.
#' `list(c("petiole_length/PETIlng", "rachis_length/RACHlng"))`). If `NULL`
#' (default), only heatmaps are produced.
#'
#' @param output_dir Character string giving the directory where outputs are
#' saved. If `NULL` (default), a folder named `"Figs_Similarity/<DDMonYYYY>"`
#' is created in the working directory, using the current date.
#'
#' @param file_formats Character vector specifying which image formats to
#' save plots in. One or more of `"pdf"`, `"jpeg"`. Default is
#' `c("pdf", "jpeg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' per-block cluster summaries, skipped blocks, and the output directory to
#' the console.
#'
#' @return
#' Invisibly, a named list with three elements:
#'
#' \describe{
#'   \item{by_block}{A named list of data frames, one per analyzed block, each
#'   containing the pairwise correlation statistics for that block.}
#'   \item{combined}{A single data frame combining the results of all
#'   analyzed blocks, sorted by absolute correlation strength, or `NULL` if no
#'   block produced results.}
#'   \item{output_dir}{The directory where all outputs were saved.}
#' }
#'
#' As a side effect, correlation heatmaps, per-block and combined Excel
#' result tables, and (if requested) scatter plots are written to
#' `output_dir`.
#'
#' @seealso [morph_matrix_setting()]
#'
#' @examples
#' \dontrun{
#' morph <- morph_matrix_setting(
#'   xlsx_path = "morphological_data.xlsx",
#'   taxon_col = "species",
#'   trait_name_type = "code",
#'   veg_first = "petiole_length/PETIlng",
#'   flo_first = "inflorescence_length/INFLlng",
#'   fru_first = "fruit_stipe_length/FRSTlng"
#' )
#'
#' # Run pairwise similarity across every block and block combination
#' sim <- morph_similarity(morph, sim_method = "spearman")
#'
#' # Restrict to a single block and add a scatter plot for one trait pair
#' sim_veg <- morph_similarity(
#'   morph,
#'   sim_blocks = "veg",
#'   sim_pairs  = list(c("petiole_length/PETIlng", "leaf_length/LEAFlng"))
#' )
#'
#' # Inspect combined results
#' sim$combined
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_point geom_smooth annotate
#'   labs theme_bw theme element_text element_blank scale_fill_gradient2
#' @importFrom cowplot save_plot
#' @importFrom openxlsx write.xlsx
#' @importFrom stats cor.test p.adjust setNames na.omit
#'
#' @export

morph_similarity <- function(analysis_data,
                                sim_blocks = "all",
                                sim_method = c("pearson", "spearman", "kendall"),
                                sim_p_adjust = "BH",
                                sim_alpha = 0.05,
                                sim_min_n = 5,
                                sim_pairs = NULL,
                                output_dir = NULL,
                                file_formats = c("pdf", "jpeg"),
                                verbose = TRUE) {


# ---- Validate input ---------------------------------------------------------

if (!is.data.frame(analysis_data)) {
  stop("`analysis_data` must be a data frame.")
}

if (!"taxon" %in% names(analysis_data)) {
  stop("`analysis_data` must contain a `taxon` column")
}

if (is.null(rownames(analysis_data))) {
  stop("`analysis_data` must have specimen IDs as row names.")
}

base_cols <- attr(analysis_data, "base_cols")

if (!is.list(base_cols) || is.null(names(base_cols)) || any(names(base_cols) == "")) {
  stop(
    "`analysis_data` has no valid \"base_cols\" attribute. ",
    "Build `analysis_data` with `morph_matrix_setting()`, which attaches ",
    "this attribute automatically."
  )
}

sim_method <- match.arg(sim_method)

if (!is.logical(verbose) || length(verbose) != 1L) {
  stop("`verbose` must be TRUE or FALSE.")
}

if (verbose) {
  message("Running similarity analysis (", sim_method, " distance)")
}

# ---- Output folder -----------------------------------------------------------

if (is.null(output_dir)) {
  output_dir <- file.path("Figs_Similarity", format(Sys.time(), "%d%b%Y"))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Resolve blocks to run ----------------------------------------------------

if (identical(sim_blocks, "all")) {
  combo_defs <- list(
    vegflo    = c("veg", "flo"),
    vegfru    = c("veg", "fru"),
    flofru    = c("flo", "fru"),
    vegflofru = c("veg", "flo", "fru")
  )

  combos <- lapply(combo_defs, function(parts) {
    if (!all(parts %in% names(base_cols))) return(NULL)
    unlist(base_cols[parts], use.names = FALSE)
  })

  run_blocks <- c(base_cols, combos)
  run_blocks <- run_blocks[!vapply(run_blocks, is.null, logical(1))]
} else {
  requested <- intersect(sim_blocks, names(base_cols))
  if (length(requested) == 0) {
    stop(
      "No valid blocks in 'sim_blocks'. Use \"all\" or any of: ",
      paste(names(base_cols), collapse = ", ")
    )
  }
  run_blocks <- base_cols[requested]
}

message("Running Pairwise Similarity Analysis...")
message("  Blocks to run: ", paste(names(run_blocks), collapse = ", "))

# ---- Run all blocks ------------------------------------------------------------

all_results <- Filter(
  Negate(is.null),
  mapply(
    .run_similarity_block,
    names(run_blocks), run_blocks,
    MoreArgs = list(
      analysis_data = analysis_data,
      sim_method    = sim_method,
      sim_min_n     = sim_min_n,
      sim_p_adjust  = sim_p_adjust,
      sim_alpha     = sim_alpha,
      output_dir    = output_dir,
      file_formats  = file_formats
    ),
    SIMPLIFY = FALSE
  )
)

# ---- Save combined Excel --------------------------------------------------------

combined <- NULL

if (length(all_results) > 0) {
  combined        <- do.call(rbind, all_results)
  combined$trait1 <- .clean_label(combined$trait1)
  combined$trait2 <- .clean_label(combined$trait2)
  combined        <- combined[order(abs(combined$cor), decreasing = TRUE), ]

  openxlsx::write.xlsx(
    combined,
    file.path(output_dir, "Similarity_all_results.xlsx")
  )
}

# ---- Optional scatter plots for specific pairs -----------------------------------

if (!is.null(sim_pairs)) {
  message("\n  Generating scatter plots for specified pairs...")

  for (pair in sim_pairs) {
    .scatter_pair(
      pair          = pair,
      analysis_data = analysis_data,
      sim_method    = sim_method,
      sim_min_n     = sim_min_n,
      output_dir    = output_dir,
      file_formats  = file_formats
    )
  }
}

message("\nDone. Outputs saved in: ", output_dir)

invisible(list(
  by_block   = all_results,
  combined   = combined,
  output_dir = output_dir
))
}

#' Clean morphological trait labels for plotting
#'
#' @description
#' Internal helper that strips a `"/TRAITcode"` suffix (if present) from a
#' \pkg{moRphoTaxa}-style trait name and replaces underscores with spaces, for
#' use as a human-readable axis or heatmap label.
#'
#' @param x Character vector of trait names.
#'
#' @return A character vector of cleaned trait labels.
#'
#' @keywords internal
#'
#' @noRd
.clean_label <- function(x) gsub("_", " ", gsub("/.*", "", x))

#' Run a pairwise similarity analysis for one trait block
#'
#' @description
#' Internal helper used by [morph_similarity()] to compute all pairwise
#' correlations within a single trait block, save the resulting heatmap and
#' Excel table, and return the block's result table.
#'
#' @param block_name Character. Name of the trait block (used in messages and
#' output file names).
#' @param block_columns Character vector of trait column names belonging to
#' this block.
#' @param analysis_data The morphometric data frame passed to
#' [morph_similarity()].
#' @param sim_method,sim_min_n,sim_p_adjust,sim_alpha,output_dir,file_formats
#' Passed through from [morph_similarity()].
#'
#' @return A data frame of pairwise correlation statistics for this block, or
#' `NULL` (invisibly) if the block has fewer than two traits or no pair meets
#' `sim_min_n`.
#'
#' @keywords internal
#'
#' @noRd
.run_similarity_block <- function(block_name, block_columns,
                                analysis_data, sim_method, sim_min_n,
                                sim_p_adjust, sim_alpha,
                                output_dir, file_formats) {

message("\n  --- Block: ", block_name, " ---")

valid_cols <- intersect(block_columns, names(analysis_data))
block_data <- analysis_data[, valid_cols, drop = FALSE]

if (ncol(block_data) < 2) {
  warning("Skipping '", block_name, "': fewer than 2 traits.")
  return(invisible(NULL))
}

trait_labels <- stats::setNames(.clean_label(valid_cols), valid_cols)

# ---- All pairwise correlations ---------------------------------------------

pairs <- utils::combn(valid_cols, 2, simplify = FALSE)

results <- lapply(pairs, function(pair) {
  x <- block_data[[pair[1]]]
  y <- block_data[[pair[2]]]
  valid <- !is.na(x) & !is.na(y)
  if (sum(valid) < sim_min_n) return(NULL)

  tryCatch({
    res <- stats::cor.test(x[valid], y[valid], method = sim_method)
    data.frame(
      trait1 = pair[1],
      trait2 = pair[2],
      n = sum(valid),
      cor = res$estimate,
      t_stat = if (!is.null(res$statistic)) res$statistic else NA,
      df = if (!is.null(res$parameter)) res$parameter else NA,
      p_value = res$p.value
    )
  }, error = function(e) NULL)
})

results <- do.call(rbind, Filter(Negate(is.null), results))
if (is.null(results) || nrow(results) == 0) return(invisible(NULL))

results$p_adjusted <- stats::p.adjust(results$p_value, method = sim_p_adjust)
results$significant <- results$p_adjusted < sim_alpha
results$block <- block_name

n_sig <- sum(results$significant, na.rm = TRUE)
cat("  Pairs tested:", nrow(results),
    "| Significant (adj. p <", sim_alpha, "):", n_sig, "\n")

# ---- Correlation matrix for heatmap -----------------------------------------

cor_matrix <- matrix(NA, nrow = length(valid_cols), ncol = length(valid_cols),
                     dimnames = list(valid_cols, valid_cols))
diag(cor_matrix) <- 1

for (i in seq_len(nrow(results))) {
  t1 <- results$trait1[i]
  t2 <- results$trait2[i]
  cor_matrix[t1, t2] <- results$cor[i]
  cor_matrix[t2, t1] <- results$cor[i]
}

# ---- Heatmap -----------------------------------------------------------------

heat_df <- as.data.frame(as.table(cor_matrix))
names(heat_df) <- c("trait1", "trait2", "cor")
heat_df$trait1 <- factor(trait_labels[as.character(heat_df$trait1)],
                         levels = trait_labels[valid_cols])
heat_df$trait2 <- factor(trait_labels[as.character(heat_df$trait2)],
                         levels = trait_labels[valid_cols])

p_heat <- ggplot2::ggplot(heat_df, ggplot2::aes(x = .data[["trait1"]], y = .data[["trait2"]], fill = .data[["cor"]])) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::scale_fill_gradient2(low      = "#d73027",
                                mid      = "white",
                                high     = "#1a9850",
                                midpoint = 0,
                                limits   = c(-1, 1),
                                name     = "r") +
  ggplot2::labs(x = NULL, y = NULL,
                title = paste0("Correlation Matrix \u2014 ", block_name)) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
                 axis.text.y  = ggplot2::element_text(size = 7),
                 panel.border = ggplot2::element_blank(),
                 panel.grid   = ggplot2::element_blank())

plot_size <- max(5, length(valid_cols) * 0.35)

if ("pdf" %in% file_formats) {
  cowplot::save_plot(
    file.path(output_dir, paste0("Similarity_heatmap_", block_name, ".pdf")),
    p_heat, ncol = 1, nrow = 1,
    base_height = plot_size, base_aspect_ratio = 1
  )
}
if ("jpeg" %in% file_formats) {
  cowplot::save_plot(
    file.path(output_dir, paste0("Similarity_heatmap_", block_name, ".jpeg")),
    p_heat, ncol = 1, nrow = 1,
    base_height = plot_size, base_aspect_ratio = 1
  )
}

# ---- Save Excel ---------------------------------------------------------------

out <- results[, c("trait1", "trait2", "n", "cor",
                          "p_value", "p_adjusted", "significant")]
out$trait1 <- trait_labels[out$trait1]
out$trait2 <- trait_labels[out$trait2]
out <- out[order(abs(out$cor), decreasing = TRUE), ]

openxlsx::write.xlsx(
  out,
  file.path(output_dir, paste0("Similarity_", block_name, ".xlsx"))
)

results
}

#' Plot a single trait-pair scatter plot with correlation annotation
#'
#' @description
#' Internal helper used by [morph_similarity()] to produce an annotated
#' scatter plot (linear fit, correlation coefficient, significance, and
#' sample size) for one user-specified pair of trait columns.
#'
#' @param pair Character vector of length two giving the two trait column
#' names to plot.
#' @param analysis_data The morphometric data frame passed to
#' [morph_similarity()].
#' @param sim_method,sim_min_n,output_dir,file_formats Passed through from
#' [morph_similarity()].
#'
#' @return `NULL`, invisibly. Called for the side effect of saving a scatter
#' plot to `output_dir`.
#'
#' @keywords internal
#'
#' @noRd
.scatter_pair <- function(pair, analysis_data, sim_method, sim_min_n,
                          output_dir, file_formats) {

  col1 <- pair[1]
  col2 <- pair[2]
  missing <- setdiff(c(col1, col2), names(analysis_data))

  if (length(missing) > 0) {
    warning("Skipping pair \u2014 column(s) not found: ",
            paste(missing, collapse = ", "))
    return(invisible(NULL))
  }

  df <- stats::na.omit(analysis_data[, c(col1, col2)])
  colnames(df) <- c("x", "y")

  if (nrow(df) < sim_min_n) {
    warning("Skipping pair (", col1, " ~ ", col2,
            "): fewer than ", sim_min_n, " complete observations.")
    return(invisible(NULL))
  }

  res <- stats::cor.test(df$x, df$y, method = sim_method)
  sig_label <- if (res$p.value < 0.001) "p < 0.001" else
    paste0("p = ", round(res$p.value, 3))

  p_scatter <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]])) +
    ggplot2::geom_point(alpha = 0.6, size = 2, color = "steelblue") +
    ggplot2::geom_smooth(method  = "lm", formula = y ~ x,
                         se        = TRUE,
                         color     = "red",
                         linetype  = "dashed",
                         linewidth = 0.8) +
    ggplot2::annotate("text",
                      x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
                      label = paste0("r = ", round(res$estimate, 3),
                                     "\n", sig_label,
                                     "\nn = ", nrow(df)),
                      size = 3.5) +
    ggplot2::labs(x = .clean_label(col1),
                  y = .clean_label(col2)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

    pair_name <- paste0(gsub("[^a-zA-Z0-9]", "_", col1),
                      "_vs_",
                      gsub("[^a-zA-Z0-9]", "_", col2))

  if ("pdf" %in% file_formats) {
    cowplot::save_plot(
      file.path(output_dir, paste0("Scatter_", pair_name, ".pdf")),
      p_scatter, ncol = 1, nrow = 1,
      base_height = 6, base_aspect_ratio = 1.3
    )
  }
  if ("jpeg" %in% file_formats) {
    cowplot::save_plot(
      file.path(output_dir, paste0("Scatter_", pair_name, ".jpeg")),
      p_scatter, ncol = 1, nrow = 1,
      base_height = 6, base_aspect_ratio = 1.3
    )
  }

  cat("  Scatter saved:", .clean_label(col1), "~", .clean_label(col2),
      "| r =", round(res$estimate, 3), "|", sig_label, "\n")

  invisible(NULL)
}
