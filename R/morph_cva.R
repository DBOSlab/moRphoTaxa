#' Run canonical variate analysis (CVA) on morphological trait blocks
#'
#' @description
#' Performs a canonical variate analysis (CVA), via a one-way MANOVA,
#' on morphological trait blocks derived from a morphometric matrix, and
#' generates a standard set of diagnostic figures and spreadsheets for each
#' trait block analyzed. The function is designed to work directly on the
#' objects produced by [morph_matrix_setting()].
#'
#' For each requested trait block, the function fits a one-way MANOVA on the
#' (scaled) morphological traits by taxon, derives canonical variates from
#' the MANOVA sums-of-squares-and-cross-products matrices, saves a canonical-scores scatterplot (Can1 x Can2)
#' with confidence ellipses per taxon, saves a barplot of standardized
#' canonical coefficients for the first canonical axis, saves a barplot of
#' structure coefficients (trait-axis correlations) for the first canonical
#' axis, and saves an Excel workbook with the block summary, canonical
#' scores, standardized coefficients, and structure coefficients.
#'
#' @details
#' `morph_cva()` expects `analysis_data` as produced by
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
#' The `cva_blocks` argument controls which trait sets are analyzed:
#'
#' \itemize{
#'   \item `"all"`: runs CVA on each individual block in `base_cols`
#'   (typically `"veg"`, `"flo"`, `"fru"`), as well as on every pairwise and
#'   full combination of these blocks (`"vegflo"`, `"vegfru"`, `"flofru"`,
#'   `"vegflofru"`);
#'   \item a character vector naming one or more blocks in `base_cols`: runs
#'   CVA only on the requested block(s), individually.
#' }
#'
#' For each trait block, the grouping variable is the `taxon` column of
#' `analysis_data`. Rows with a missing taxon are dropped, and any taxon
#' with fewer than `cva_min_n` specimens (after that step) is removed from
#' the block entirely, with a message identifying the removed taxa. Blocks
#' left with fewer than 2 taxa are skipped with a warning.
#'
#' Within each retained block, missing trait values are imputed with the
#' column mean, and any trait with zero variance is dropped. A block is
#' skipped, with a warning, if fewer than 2 traits remain, or if the number
#' of specimens does not exceed the number of taxa. Remaining traits are
#' centred and scaled ([base::scale()]) before the one-way
#' [stats::manova()] is fit and canonical variates are extracted from its
#' within- and between-group SSCP matrices. The computation is the one used
#' by `candisc::candisc()` for a one-way design (same axes, scores,
#' standardized and structure coefficients), implemented internally so that
#' \pkg{moRphoTaxa} does not depend on \pkg{candisc} and its 3D-graphics
#' dependency \pkg{rgl}.
#'
#' Taxon colours are taken from the `"taxon_colors"` attribute attached to
#' `analysis_data` by [morph_matrix_setting()], so that each taxon keeps the
#' same colour used in other plots (e.g. [morph_pca()],
#' [morph_boxplots()]). If this attribute is missing, or lacks an entry for
#' a taxon retained in a given block, colours are generated automatically
#' from the `viridis` palette for the affected taxa, with a warning.
#'
#' Output files are written to a date-specific directory inside `Figs_CVA`.
#' For each analyzed block, the following files are generated:
#'
#' \itemize{
#'   \item a canonical-scores scatterplot (Can1 x Can2, or Can1 alone if only
#'   one canonical axis is available), with a normal confidence ellipse per
#'   taxon at level `cva_ellipse` (`"CVA_scatter_"` prefix);
#'   \item a barplot of standardized canonical coefficients for Can1
#'   (`"CVA_coefficients_"` prefix);
#'   \item a barplot of structure coefficients for Can1, with reference lines
#'   at |r| = 0.3 (`"CVA_structure_"` prefix);
#'   \item an Excel workbook with the block summary, canonical scores,
#'   standardized coefficients, and structure coefficients
#'   (`"CVA_"` prefix, `.xlsx`).
#' }
#'
#' Plots are saved in the formats specified by `file_formats`. Once every
#' requested block has been processed, a `"CVA_summary.xlsx"` file
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
#' @param cva_blocks Character. Either `"all"`, to run CVA on every
#' individual block in `base_cols` plus all pairwise and full combinations
#' of these blocks, or a character vector naming one or more blocks in
#' `base_cols` to analyze individually. Default is `"all"`.
#'
#' @param cva_ellipse Numeric. Confidence level (between 0 and 1) used for
#' the per-taxon normal confidence ellipses drawn on the canonical-scores
#' scatterplot. Default is `0.95`.
#'
#' @param cva_min_n Integer. Minimum number of specimens a taxon must have,
#' within a given trait block, to be retained in the CVA for that block.
#' Taxa with fewer specimens are dropped block-by-block. Default is `3`.
#'
#' @param file_formats Character vector specifying which image formats to
#' save. One or more of `"pdf"` and `"jpeg"`. Default is
#' `c("pdf", "jpeg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' including removed taxa, skipped blocks, canonical-axis summaries, and the
#' output directory, to the console.
#'
#' @return
#' Invisibly returns a data frame summarizing the blocks successfully
#' analyzed, with one row per block and columns `block`, `n`, `n_groups`,
#' `can1_pct`, and `can1_r`; or invisibly returns `NULL` if no block could be
#' analyzed.
#'
#' The function is primarily called for its side effects: canonical-scores
#' scatterplots, coefficient and structure barplots, and per-block and
#' summary Excel workbooks are written to disk for each analyzed trait
#' block.
#'
#' @seealso
#' [morph_matrix_setting()],
#' [morph_pca()],
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
#' # Run CVA on every individual block and every block combination
#' morph_cva(
#'   analysis_data = res,
#'   cva_blocks = "all"
#' )
#'
#' # Run CVA on the vegetative block only, with a 90% confidence ellipse
#' morph_cva(
#'   analysis_data = res,
#'   cva_blocks = "veg",
#'   cva_ellipse = 0.90
#' )
#'
#' # Run CVA on selected blocks, requiring at least 5 specimens per taxon,
#' # saving PDF output only
#' morph_cva(
#'   analysis_data = res,
#'   cva_blocks = c("veg", "flo"),
#'   cva_min_n = 5,
#'   file_formats = "pdf"
#' )
#' }
#'
#' @importFrom stats manova var setNames
#' @importFrom ggplot2 aes coord_flip element_blank element_text geom_col
#'   geom_hline geom_point geom_vline ggplot ggtitle labs scale_color_manual
#'   scale_fill_manual stat_ellipse theme theme_bw
#' @importFrom rlang .data
#' @importFrom cowplot save_plot
#' @importFrom viridis viridis
#' @importFrom openxlsx write.xlsx
#'
#' @export

morph_cva <- function(analysis_data,
                      cva_blocks = "all",
                      cva_ellipse = 0.95,
                      cva_min_n = 3,
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

if (!is.character(cva_blocks) || length(cva_blocks) < 1) {
  stop("`cva_blocks` must be \"all\" or a character vector of block names.")
}

if (!is.numeric(cva_ellipse) ||
    length(cva_ellipse) != 1L ||
    cva_ellipse <= 0 || cva_ellipse >= 1) {
  stop("`cva_ellipse` must be a single numeric value between 0 and 1.")
}

if (!is.numeric(cva_min_n) ||
    length(cva_min_n) != 1L ||
    cva_min_n < 1) {
  stop("`cva_min_n` must be a positive integer.")
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

output_dir <- file.path("Figs_CVA", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (verbose) {
  message("Running CVA analysis")
}

# ---- Resolve blocks ----------------------------------------------------------

if (identical(cva_blocks, "all")) {
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
requested <- intersect(cva_blocks, names(base_cols))
if (length(requested) == 0) {
  stop(
    "No valid blocks in `cva_blocks`. Use \"all\" or any of: ",
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

# ---- Helper: run CVA for one block -------------------------------------------

run_cva_block <- function(block_name, block_columns) {

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
small_groups <- names(group_counts)[group_counts < cva_min_n]
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
kept_cols_names <- valid_cols[keep_cols]

if (ncol(trait_data) < 2 || nrow(trait_data) < n_groups + 1) {
  warning("Skipping '", block_name, "': not enough data.")
  return(NULL)
}

if (verbose) {
  message(
    "  Specimens: ", nrow(trait_data),
    " | Traits: ", ncol(trait_data),
    " | Groups: ", n_groups
  )
}

# Scale traits
trait_scaled <- as.data.frame(scale(trait_data))
trait_scaled$group <- groups

# ---- Fit MANOVA + CVA ------------------------------------------------------

manova_fit <- tryCatch(
  stats::manova(
    as.matrix(trait_scaled[, kept_cols_names, drop = FALSE]) ~ group,
    data = trait_scaled
  ),
  error = function(e) {
    warning("MANOVA failed for '", block_name, "': ", conditionMessage(e))
    NULL
  }
)
if (is.null(manova_fit)) return(NULL)

cva_res <- tryCatch(
  .candisc_fit(manova_fit),
  error = function(e) {
    warning("CVA failed for '", block_name, "': ", conditionMessage(e))
    NULL
  }
)
if (is.null(cva_res)) return(NULL)

# ---- Print results ---------------------------------------------------------

n_axes <- length(cva_res$canrsq)
pct_var <- round(cva_res$pct, 1)  # already a percentage
can_corr <- round(sqrt(cva_res$canrsq), 3)

if (verbose) {
  cat("\n  Block          :", block_name, "\n")
  cat("  Specimens      :", nrow(trait_data), "\n")
  cat("  Groups         :", paste(levels(groups), collapse = ", "), "\n")
  cat("  Canonical axes :", n_axes, "\n")
  for (i in seq_len(n_axes)) {
    cat(sprintf(
      "  Can%d: %.1f%% variance | r = %.3f\n",
      i, pct_var[i], can_corr[i]
    ))
  }
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

# ---- CVA scatter plot (Can1 vs Can2) ---------------------------------------

scores_df <- as.data.frame(cva_res$scores)
scores_df$group <- groups
scores_df$taxon <- as.character(species_all[rownames(scores_df)])

x_ax <- "Can1"
y_ax <- if (n_axes >= 2) "Can2" else "Can1"

x_lab <- paste0("Can1 (", pct_var[1], "% | r = ", can_corr[1], ")")
y_lab <- if (n_axes >= 2) {
  paste0("Can2 (", pct_var[2], "% | r = ", can_corr[2], ")")
} else {
  "Can2"
}

p_scatter <- ggplot2::ggplot(
  scores_df,
  ggplot2::aes(
    x = .data[[x_ax]],
    y = if (n_axes >= 2) .data[[y_ax]] else 0
  )
) +
  {
    if (n_axes >= 2) {
      ggplot2::stat_ellipse(
        ggplot2::aes(color = .data[["group"]]),
        level = cva_ellipse,
        type = "norm",
        linetype = "dashed",
        linewidth = 0.8,
        show.legend = FALSE
      )
    }
  } +
  ggplot2::geom_point(
    ggplot2::aes(fill = .data[["group"]]), shape = 21, size = 3,
    colour = "white", alpha = 0.8) +
  ggplot2::scale_fill_manual(values = color_vector) +
  ggplot2::scale_color_manual(values = color_vector) +
  ggplot2::geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::labs(
    x = x_lab, y = y_lab, fill = "taxon",
    title = paste0("CVA \u2014 ", block_name)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(face = "italic")
  )

save_outputs(
  paste0("CVA_scatter_", block_name), p_scatter,
  base_height = 7, base_aspect_ratio = 1.3
)

# ---- Standardized coefficients plot ----------------------------------------

std_coef <- as.data.frame(cva_res$coeffs.std)
std_coef$trait <- gsub("_", " ", gsub("/.*", "", rownames(std_coef)))
std_coef <- std_coef[order(abs(std_coef$Can1), decreasing = TRUE), ]
std_coef$trait <- factor(std_coef$trait, levels = rev(std_coef$trait))

p_coef <- ggplot2::ggplot(
  std_coef,
  ggplot2::aes(
    x = .data[["trait"]],
    y = .data[["Can1"]],
    fill = .data[["Can1"]] > 0
  )
) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(
    values = c("TRUE" = "#1a9850", "FALSE" = "#d73027"),
    guide = "none"
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL, y = "Standardized Canonical Coefficient (Can1)",
    title = paste0("CVA Coefficients \u2014 ", block_name)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 7)
  )

save_outputs(
  paste0("CVA_coefficients_", block_name), p_coef,
  base_height = max(4, nrow(std_coef) * 0.3),
  base_aspect_ratio = 1.5
)

# ---- Structure coefficients (trait-axis correlations) ----------------------

struct_coef <- as.data.frame(cva_res$structure)
struct_coef$trait <- gsub("_", " ", gsub("/.*", "", rownames(struct_coef)))
struct_coef <- struct_coef[order(abs(struct_coef$Can1), decreasing = TRUE), ]
struct_coef$trait <- factor(struct_coef$trait, levels = rev(struct_coef$trait))

p_struct <- ggplot2::ggplot(
  struct_coef,
  ggplot2::aes(
    x = .data[["trait"]],
    y = .data[["Can1"]],
    fill = .data[["Can1"]] > 0
  )
) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(
    values = c("TRUE" = "#1a9850", "FALSE" = "#d73027"),
    guide = "none"
  ) +
  ggplot2::geom_hline(
    yintercept = c(-0.3, 0.3), linetype = "dashed",
    color = "gray50", linewidth = 0.5
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL, y = "Structure Coefficient (Can1)",
    title = paste0("CVA Structure \u2014 ", block_name),
    subtitle = "Dashed lines = |r| = 0.3 (minimum meaningful loading)"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 7)
  )

save_outputs(
  paste0("CVA_structure_", block_name), p_struct,
  base_height = max(4, nrow(struct_coef) * 0.3),
  base_aspect_ratio = 1.5
)

# ---- Save Excel ------------------------------------------------------------

openxlsx::write.xlsx(
  list(
    summary = data.frame(
      block = block_name,
      n = nrow(trait_data),
      n_groups = n_groups,
      n_axes = n_axes,
      can1_pct = pct_var[1],
      can1_r = can_corr[1],
      can2_pct = if (n_axes >= 2) pct_var[2] else NA,
      can2_r = if (n_axes >= 2) can_corr[2] else NA
    ),
    scores = scores_df,
    std_coefficients = as.data.frame(cva_res$coeffs.std),
    structure_coefs = as.data.frame(cva_res$structure)
  ),
  file.path(output_dir, paste0("CVA_", block_name, ".xlsx"))
)

if (verbose) {
  message("  Done: ", block_name)
}

data.frame(
  block = block_name,
  n = nrow(trait_data),
  n_groups = n_groups,
  can1_pct = pct_var[1],
  can1_r = can_corr[1]
)
}

# ---- Run all blocks ----------------------------------------------------------

all_results <- Filter(
  Negate(is.null),
  mapply(
    run_cva_block,
    names(run_blocks), run_blocks,
    SIMPLIFY = FALSE
  )
)

# ---- Save summary ------------------------------------------------------------

if (length(all_results) > 0) {
  summary_df <- do.call(rbind, all_results)
  openxlsx::write.xlsx(summary_df, file.path(output_dir, "CVA_summary.xlsx"))
} else {
  summary_df <- NULL
}

if (verbose) {
  message("\nCVA analysis done. Outputs saved in: ", output_dir)
}

invisible(summary_df)
}


# Canonical discriminant analysis for a one-way MANOVA fit.
# Same algorithm as candisc::candisc() (type-II SSCP matrices for a single
# term): canonical axes are the eigenvectors of H relative to E, where E is
# the within-group and H the between-group SSCP matrix.
.candisc_fit <- function(fit) {
  mf <- stats::model.frame(fit)
  Y  <- stats::model.response(mf)
  group <- mf[["group"]]

  E <- crossprod(stats::residuals(fit))
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  H <- crossprod(Yc) - E
  dfe <- fit$df.residual
  dfh <- nlevels(droplevels(group)) - 1

  # E = t(Tm) %*% Tm
  eE <- eigen(E, symmetric = TRUE)
  Tm <- t(eE$vectors %*% diag(sqrt(eE$values), nrow = length(eE$values)))
  eInv <- solve(Tm)
  dc <- eigen(t(eInv) %*% H %*% eInv, symmetric = TRUE)

  ndim <- min(dfh, sum(dc$values > 0))
  if (ndim < 1) stop("no canonical dimensions with positive eigenvalues")
  cn <- paste0("Can", seq_len(ndim))

  coeffs.raw <- (eInv %*% dc$vectors * sqrt(dfe))[, seq_len(ndim), drop = FALSE]
  dimnames(coeffs.raw) <- list(colnames(Y), cn)
  coeffs.std <- diag(sqrt(diag(E / dfe)), nrow = ncol(Y)) %*% coeffs.raw
  dimnames(coeffs.std) <- list(colnames(Y), cn)

  scores <- Yc %*% coeffs.raw
  colnames(scores) <- cn

  list(
    eigenvalues = dc$values,
    canrsq      = dc$values[seq_len(ndim)] / (1 + dc$values[seq_len(ndim)]),
    pct         = 100 * dc$values / sum(dc$values),
    rank        = ndim,
    coeffs.raw  = coeffs.raw,
    coeffs.std  = coeffs.std,
    structure   = stats::cor(Yc, scores),
    scores      = data.frame(group = group, as.data.frame(scores),
                             row.names = rownames(Y))
  )
}
