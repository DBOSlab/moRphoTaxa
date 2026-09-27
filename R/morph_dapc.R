#' Run discriminant analysis of principal components (DAPC) on morphological
#' trait blocks
#'
#' @description
#' Performs a discriminant analysis of principal components (DAPC), via
#' [stats::prcomp()] followed by [adegenet::dapc()], on morphological trait
#' blocks derived from a morphometric matrix, and generates a standard set
#' of diagnostic figures and spreadsheets for each trait block analyzed. The
#' function is designed to work directly on the objects produced by
#' [morph_matrix_setting()].
#'
#' For each requested trait block, the function fits a PCA on the scaled
#' morphological traits, retains a subset of principal components, runs
#' [adegenet::dapc()] on the retained components using taxon as the
#' grouping variable, saves a discriminant-scores scatterplot (LD1 x LD2)
#' with confidence ellipses per taxon, saves a membership-probability
#' barplot showing each specimen's assignment probability to every taxon,
#' and saves an Excel workbook with the posterior membership probabilities
#' and discriminant scores.
#'
#' @details
#' `morph_dapc()` expects `analysis_data` as produced by
#' [morph_matrix_setting()]: a numeric data frame with `id` and `taxon`
#' columns plus one column per morphological trait, with specimen IDs as row
#' names. The trait blocks (`"veg"`, `"flo"`, `"fru"`, or any subset thereof)
#' are not passed in separately — they are read directly from the
#' `"base_cols"` attribute that [morph_matrix_setting()] attaches to
#' `analysis_data`. If `analysis_data` has no `"base_cols"` attribute (e.g.
#' it was not built with [morph_matrix_setting()], or the attribute was
#' dropped by an intervening subsetting operation), the function stops with
#' an informative error.
#'
#' The `dapc_blocks` argument controls which trait sets are analyzed:
#'
#' \itemize{
#'   \item `"all"`: runs DAPC on each individual block in `base_cols`
#'   (typically `"veg"`, `"flo"`, `"fru"`), as well as on every pairwise and
#'   full combination of these blocks (`"vegflo"`, `"vegfru"`, `"flofru"`,
#'   `"vegflofru"`);
#'   \item a character vector naming one or more blocks in `base_cols`: runs
#'   DAPC only on the requested block(s), individually.
#' }
#'
#' For each trait block, the grouping variable is the `taxon` column of
#' `analysis_data`. Rows with a missing taxon are dropped, and any taxon
#' with fewer than `dapc_min_n` specimens (after that step) is removed from
#' the block entirely, with a message identifying the removed taxa. Blocks
#' left with fewer than 2 taxa are skipped with a warning.
#'
#' Within each retained block, missing trait values are imputed with the
#' column mean, and any trait with zero variance is dropped. A block is
#' skipped, with a warning, if fewer than 2 traits remain, or if the number
#' of specimens does not exceed the number of taxa. Remaining traits are
#' centred and scaled ([base::scale()]) before a principal component
#' analysis is fit with [stats::prcomp()] (`scale. = FALSE`, since the
#' traits are already scaled).
#'
#' The number of principal components retained for the discriminant step is
#' controlled by `dapc_n_pca`: if `NULL` (default), the smallest number of
#' components explaining at least 80% of the cumulative variance is used
#' (with a minimum of 2 and a maximum of `nrow - 1`); otherwise the supplied
#' value is used, capped at the number of components available. The number
#' of discriminant functions is controlled by `dapc_n_da`: if `NULL`
#' (default), `n_groups - 1` is used (the maximum possible); otherwise the
#' supplied value is used, capped at `n_groups - 1`. Discriminant functions
#' are then fit with [adegenet::dapc()] on the retained principal
#' components.
#'
#' Taxon colours are taken from the `"taxon_colors"` attribute attached to
#' `analysis_data` by [morph_matrix_setting()], so that each taxon keeps the
#' same colour used in other plots (e.g. [morph_pca()], [morph_cva()],
#' [morph_boxplots()]). If this attribute is missing, or lacks an entry for
#' a taxon retained in a given block, colours are generated automatically
#' from the `viridis` palette for the affected taxa, with a warning.
#'
#' Output files are written to a date-specific directory inside `Figs_DAPC`.
#' For each analyzed block, the following files are generated:
#'
#' \itemize{
#'   \item a discriminant-scores scatterplot (LD1 x LD2, or LD1 alone if
#'   only one discriminant function is available), with a normal confidence
#'   ellipse per taxon at level `dapc_ellipse` (`"DAPC_scatter_"` prefix);
#'   \item a membership-probability barplot, with one panel per taxon and
#'   one bar per specimen, showing the posterior probability of assignment
#'   to each taxon (`"DAPC_membership_"` prefix);
#'   \item an Excel workbook with the posterior membership probabilities and
#'   discriminant scores (`"DAPC_"` prefix, `.xlsx`).
#' }
#'
#' Plots are saved in the formats specified by `file_formats`. Once every
#' requested block has been processed, a `"DAPC_summary.xlsx"` file
#' collecting one row per successfully analyzed block is written to the same
#' output directory.
#'
#' @param analysis_data A numeric data frame containing morphological traits
#' as columns and specimens as rows, with `id` and `taxon` columns and
#' specimen identifiers as row names. Typically the `analysis_data` element
#' returned by [morph_matrix_setting()]. If this object carries a
#' `"taxon_colors"` attribute, those colours are reused; otherwise a
#' `viridis` palette is generated automatically.
#'
#' @param dapc_blocks Character. Either `"all"`, to run DAPC on every
#' individual block in `base_cols` plus all pairwise and full combinations
#' of these blocks, or a character vector naming one or more blocks in
#' `base_cols` to analyze individually. Default is `"all"`.
#'
#' @param dapc_n_pca Integer or `NULL`. Number of principal components to
#' retain before the discriminant step. If `NULL` (default), the number of
#' components needed to explain at least 80% of the cumulative variance is
#' used automatically for each block.
#'
#' @param dapc_n_da Integer or `NULL`. Number of discriminant functions to
#' retain. If `NULL` (default), `n_groups - 1` is used for each block (the
#' maximum possible number of discriminant axes).
#'
#' @param dapc_min_n Integer. Minimum number of specimens a taxon must have,
#' within a given trait block, to be retained in the DAPC for that block.
#' Taxa with fewer specimens are dropped block-by-block. Default is `3`.
#'
#' @param dapc_ellipse Numeric. Confidence level (between 0 and 1) used for
#' the per-taxon normal confidence ellipses drawn on the discriminant-scores
#' scatterplot. Default is `0.95`.
#'
#' @param file_formats Character vector specifying which image formats to
#' save. One or more of `"pdf"` and `"jpeg"`. Default is
#' `c("pdf", "jpeg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' including removed taxa, skipped blocks, retained components, and the
#' output directory, to the console.
#'
#' @return
#' Invisibly returns a data frame summarizing the blocks successfully
#' analyzed, with one row per block and columns `block`, `n`, `n_groups`,
#' `n_pca`, and `n_da`; or invisibly returns `NULL` if no block could be
#' analyzed.
#'
#' The function is primarily called for its side effects: discriminant-score
#' scatterplots, membership-probability barplots, and per-block and summary
#' Excel workbooks are written to disk for each analyzed trait block.
#'
#' @seealso
#' [morph_matrix_setting()],
#' [morph_pca()],
#' [morph_cva()],
#' [morph_boxplots()]
#'
#' @examples
#' \dontrun{
#' # Build the morphometric matrix
#' res <- morph_matrix_setting(
#'   xlsx_path = "output_data/all_data.xlsx",
#'   taxon_col = "taxon",
#'   trait_name_type = "code",
#'   veg_first = "petiole_length/PETIlng",
#'   flo_first = "inflorescence_length/INFLlng",
#'   fru_first = "fruit_stipe_length/FRSTlng"
#' )
#'
#' # Run DAPC on every individual block and every block combination
#' morph_dapc(
#'   analysis_data = res,
#'   dapc_blocks = "all"
#' )
#'
#' # Run DAPC on the vegetative block only, retaining 5 PCs
#' morph_dapc(
#'   analysis_data = res,
#'   dapc_blocks = "veg",
#'   dapc_n_pca = 5
#' )
#'
#' # Run DAPC on selected blocks, requiring at least 5 specimens per taxon,
#' # saving PDF output only
#' morph_dapc(
#'   analysis_data = res,
#'   dapc_blocks = c("veg", "flo"),
#'   dapc_min_n = 5,
#'   file_formats = "pdf"
#' )
#' }
#'
#' @importFrom stats prcomp var setNames
#' @importFrom adegenet dapc
#' @importFrom ggplot2 aes element_blank element_text facet_grid geom_col
#'   geom_hline geom_point geom_vline ggplot labs scale_color_manual
#'   scale_fill_manual scale_y_continuous stat_ellipse theme theme_bw vars
#' @importFrom rlang .data
#' @importFrom cowplot save_plot
#' @importFrom viridis viridis
#' @importFrom openxlsx write.xlsx
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr arrange desc
#' @importFrom tidyselect all_of
#' @importFrom grid unit
#'
#' @export

morph_dapc <- function(analysis_data,
                       dapc_blocks = "all",
                       dapc_n_pca = NULL,
                       dapc_n_da = NULL,
                       dapc_min_n = 3,
                       dapc_ellipse = 0.95,
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

if (!is.character(dapc_blocks) || length(dapc_blocks) < 1) {
  stop("`dapc_blocks` must be \"all\" or a character vector of block names.")
}

if (!is.null(dapc_n_pca) &&
    (!is.numeric(dapc_n_pca) || length(dapc_n_pca) != 1L || dapc_n_pca < 1)) {
  stop("`dapc_n_pca` must be NULL or a positive integer.")
}

if (!is.null(dapc_n_da) &&
    (!is.numeric(dapc_n_da) || length(dapc_n_da) != 1L || dapc_n_da < 1)) {
  stop("`dapc_n_da` must be NULL or a positive integer.")
}

if (!is.numeric(dapc_min_n) ||
    length(dapc_min_n) != 1L ||
    dapc_min_n < 1) {
  stop("`dapc_min_n` must be a positive integer.")
}

if (!is.numeric(dapc_ellipse) ||
    length(dapc_ellipse) != 1L ||
    dapc_ellipse <= 0 || dapc_ellipse >= 1) {
  stop("`dapc_ellipse` must be a single numeric value between 0 and 1.")
}

file_formats <- match.arg(
  file_formats,
  choices = c("pdf", "jpeg"),
  several.ok = TRUE
)

if (!is.logical(verbose) || length(verbose) != 1L) {
  stop("`verbose` must be TRUE or FALSE.")
}

# ---- Species vector ----------------------------------------------------------

species_all <- stats::setNames(
  as.character(analysis_data$taxon),
  rownames(analysis_data)
)

# ---- Output folder ------------------------------------------------------------

output_dir <- file.path("Figs_DAPC", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (verbose) {
  message("Running DAPC analysis")
}

# ---- Resolve blocks ----------------------------------------------------------

if (identical(dapc_blocks, "all")) {
  combo_defs <- list(
    vegflo = c("veg", "flo"),
    vegfru = c("veg", "fru"),
    flofru = c("flo", "fru"),
    vegflofru = c("veg", "flo", "fru")
  )

  combos <- lapply(combo_defs, function(parts) {
    if (!all(parts %in% names(base_cols))) return(NULL)
    unlist(base_cols[parts], use.names = FALSE)
  })

  run_blocks <- c(base_cols, combos)

# Drop combinations whose required source block(s) are missing
run_blocks <- run_blocks[!vapply(run_blocks, is.null, logical(1))]
} else {
requested <- intersect(dapc_blocks, names(base_cols))
if (length(requested) == 0) {
  stop(
    "No valid blocks in `dapc_blocks`. Use \"all\" or any of: ",
    paste(names(base_cols), collapse = ", ")
  )
}
run_blocks <- base_cols[requested]
}

if (verbose) {
  message("  Blocks to run: ", paste(names(run_blocks), collapse = ", "))
}

# ---- Helper: save a plot in the requested formats ----------------------------

save_outputs <- function(filename_stub, plot_obj, base_height, base_aspect_ratio, nrow = 1) {
  if ("pdf" %in% file_formats) {
    cowplot::save_plot(
      file.path(output_dir, paste0(filename_stub, ".pdf")),
      plot_obj,
      ncol = 1,
      nrow = nrow,
      base_height = base_height,
      base_aspect_ratio = base_aspect_ratio,
      base_width = NULL
    )
  }

  if ("jpeg" %in% file_formats) {
    cowplot::save_plot(
      file.path(output_dir, paste0(filename_stub, ".jpeg")),
      plot_obj,
      ncol = 1,
      nrow = nrow,
      base_height = base_height,
      base_aspect_ratio = base_aspect_ratio,
      base_width = NULL
    )
  }
}

# ---- Helper: run DAPC for one block -------------------------------------------

run_dapc_block <- function(block_name, block_columns) {

if (verbose) {
  message("\n  --- Block: ", block_name, " ---")
}

valid_cols <- intersect(block_columns, names(analysis_data))
block_data <- analysis_data[, valid_cols, drop = FALSE]

# Attach grouping variable
group_vals <- analysis_data[rownames(block_data), "taxon"]
block_data$group <- as.factor(group_vals)
block_data <- block_data[!is.na(block_data$group), ]
block_data$group <- droplevels(block_data$group)

# Remove small groups
group_counts <- table(block_data$group)
small_groups <- names(group_counts)[group_counts < dapc_min_n]
if (length(small_groups) > 0) {
  if (verbose) {
    message("  Removing small groups: ", paste(small_groups, collapse = ", "))
  }
  block_data <- block_data[!block_data$group %in% small_groups, ]
  block_data$group <- droplevels(block_data$group)
}

groups <- block_data$group
n_groups <- length(levels(groups))

if (n_groups < 2) {
  warning("Skipping '", block_name, "': fewer than 2 groups.")
  return(NULL)
}

trait_data <- block_data[, valid_cols, drop = FALSE]

# Impute column means for NAs
trait_data[] <- lapply(trait_data, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

# Remove zero-variance columns
keep_cols <- vapply(trait_data, function(x) {
  v <- stats::var(x, na.rm = TRUE)
  !is.na(v) && v > 0
}, logical(1))
trait_data <- trait_data[, keep_cols, drop = FALSE]
trait_mat <- scale(trait_data)

if (nrow(trait_mat) < n_groups + 1 || ncol(trait_mat) < 2) {
  warning("Skipping '", block_name, "': not enough data.")
  return(NULL)
}

if (verbose) {
  message(
    "  Specimens: ", nrow(trait_mat),
    " | Traits: ", ncol(trait_mat),
    " | Groups: ", n_groups
  )
}

# ---- PCA step --------------------------------------------------------------

pca_res <- stats::prcomp(trait_mat, scale. = FALSE)
var_prop <- summary(pca_res)$importance[2, ]
cum_var <- cumsum(var_prop)

if (is.null(dapc_n_pca)) {
  n_pca <- max(2, min(which(cum_var >= 0.80)[1], nrow(trait_mat) - 1))
} else {
  n_pca <- min(dapc_n_pca, ncol(pca_res$x))
}

if (verbose) {
  message(
    "  PCs retained: ", n_pca,
    " (", round(cum_var[n_pca] * 100, 1), "% variance)"
  )
}

pca_scores <- pca_res$x[, seq_len(n_pca), drop = FALSE]

# ---- DA step ---------------------------------------------------------------

n_da <- if (is.null(dapc_n_da)) n_groups - 1 else min(dapc_n_da, n_groups - 1)

dapc_res <- tryCatch(
  adegenet::dapc(
    pca_scores, grp = groups,
    n.pca = n_pca, n.da = n_da
  ),
  error = function(e) {
    warning("DAPC failed for '", block_name, "': ", conditionMessage(e))
    NULL
  }
)
if (is.null(dapc_res)) return(NULL)

# ---- Print results ---------------------------------------------------------

if (verbose) {
  cat("\n  Block      :", block_name, "\n")
  cat("  Groups     :", paste(levels(groups), collapse = ", "), "\n")
  cat("  PCs used   :", n_pca, "\n")
  cat("  DAs used   :", n_da, "\n")
}

# ---- Taxon colours -------------------------------------------------------

color_vector <- attr(analysis_data, "taxon_colors")

if (is.null(color_vector)) {

  color_vector <- stats::setNames(
    viridis::viridis(n_groups),
    levels(groups)
  )

} else {

  missing_sp <- setdiff(levels(groups), names(color_vector))

  if (length(missing_sp) > 0) {
    warning(
      "No stored colour found for: ", paste(missing_sp, collapse = ", "),
      ". Assigning automatic colours for these taxa."
    )
    color_vector <- c(
      color_vector,
      stats::setNames(viridis::viridis(length(missing_sp)), missing_sp)
    )
  }

  color_vector <- color_vector[levels(groups)]
}

# ---- Scatter plot (LD1 vs LD2) ---------------------------------------------

da_scores <- as.data.frame(dapc_res$ind.coord)
da_scores$group <- groups
da_scores$taxon <- as.character(species_all[rownames(da_scores)])

x_ax <- "LD1"
y_ax <- if (n_da >= 2) "LD2" else "LD1"

p_scatter <- ggplot2::ggplot(
  da_scores,
  ggplot2::aes(
    x = .data[[x_ax]],
    y = if (n_da >= 2) .data[[y_ax]] else 0
  )
) +
  ggplot2::stat_ellipse(
    ggplot2::aes(color = .data[["group"]]),
    type = "norm", level = dapc_ellipse,
    linetype = "dashed", linewidth = 0.8,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = .data[["group"]]), shape = 21, size = 3,
    colour = "white", alpha = 0.8
  ) +
  ggplot2::scale_fill_manual(values = color_vector) +
  ggplot2::scale_color_manual(values = color_vector) +
  ggplot2::geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::labs(
    x = x_ax, y = if (n_da >= 2) y_ax else "",
    fill = "taxon",
    title = paste0("DAPC \u2014 ", block_name)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(face = "italic")
  )

save_outputs(
  paste0("DAPC_scatter_", block_name), p_scatter,
  base_height = 7, base_aspect_ratio = 1.3
)

# ---- Membership probability plot -------------------------------------------

proba_df <- as.data.frame(dapc_res$posterior)
proba_df$id <- rownames(proba_df)
proba_df$group <- as.character(groups)

proba_long <- tidyr::pivot_longer(
  proba_df,
  cols = -tidyselect::all_of(c("id", "group")),
  names_to = "assigned_group",
  values_to = "probability"
)
proba_long <- dplyr::arrange(
  proba_long,
  .data[["group"]], dplyr::desc(.data[["probability"]])
)

proba_long$id <- factor(
  proba_long$id,
  levels = unique(proba_long$id[order(proba_long$group)])
)

p_membership <- ggplot2::ggplot(
  proba_long,
  ggplot2::aes(
    x = .data[["id"]], y = .data[["probability"]],
    fill = .data[["assigned_group"]]
  )
) +
  ggplot2::geom_col(width = 1) +
  ggplot2::scale_fill_manual(values = color_vector) +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::facet_grid(
    cols = ggplot2::vars(.data[["group"]]),
    scales = "free_x", space = "free_x"
  ) +
  ggplot2::labs(
    x = "Specimens", y = "Membership Probability",
    fill = "taxon",
    title = paste0("DAPC Membership Probabilities \u2014 ", block_name)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    panel.spacing = grid::unit(0.5, "lines"),
    strip.text = ggplot2::element_text(face = "italic", size = 8),
    legend.text = ggplot2::element_text(face = "italic"),
    panel.grid = ggplot2::element_blank()
  )

save_outputs(
  paste0("DAPC_membership_", block_name), p_membership,
  base_height = 5, base_aspect_ratio = 2
)

# ---- Save Excel ------------------------------------------------------------

openxlsx::write.xlsx(
  list(
    membership = proba_df,
    da_scores = da_scores
  ),
  file.path(output_dir, paste0("DAPC_", block_name, ".xlsx"))
)

if (verbose) {
  message("  Done: ", block_name)
}

data.frame(
  block = block_name,
  n = nrow(trait_mat),
  n_groups = n_groups,
  n_pca = n_pca,
  n_da = n_da
)
}

# ---- Run all blocks ----------------------------------------------------------

all_results <- Filter(
  Negate(is.null),
  mapply(
    run_dapc_block,
    names(run_blocks), run_blocks,
    SIMPLIFY = FALSE
  )
)

# ---- Save summary ------------------------------------------------------------

if (length(all_results) > 0) {
  summary_df <- do.call(rbind, all_results)
  openxlsx::write.xlsx(summary_df, file.path(output_dir, "DAPC_summary.xlsx"))
} else {
  summary_df <- NULL
}

if (verbose) {
  message("\nDAPC analysis done. Outputs saved in: ", output_dir)
}

invisible(summary_df)
}
