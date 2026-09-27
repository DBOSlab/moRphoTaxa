#' Linear discriminant analysis of morphological trait blocks
#'
#' @description
#' Performs linear discriminant analysis (LDA) on morphological trait blocks
#' derived from a morphometric matrix, and generates a standard set of
#' diagnostic figures and tables for each trait block analyzed. The function
#' is designed to work directly on the objects produced by
#' [morph_matrix_setting()].
#'
#' For each requested trait block, the function imputes missing trait values,
#' standardizes the traits, fits an LDA with [MASS::lda()] against a chosen
#' grouping variable, saves a scatterplot of specimens on the requested
#' discriminant axes with confidence ellipses per group, saves a bar plot of
#' trait loadings on the first discriminant axis, and saves the specimen
#' scores and trait loadings to spreadsheets.
#'
#' @details
#' `morph_lda()` expects `analysis_data` as produced by
#' [morph_matrix_setting()]: a numeric data frame with `id` and `taxon`
#' columns plus one column per morphological trait, with specimen IDs as row
#' names. The trait blocks (`"veg"`, `"flo"`, `"fru"`, or any subset thereof)
#' are not passed in separately — they are read directly from the
#' `"base_cols"` attribute that [morph_matrix_setting()] attaches to
#' `analysis_data`. If `analysis_data` has no `"base_cols"` attribute (e.g.
#' it was not built with [morph_matrix_setting()], or the attribute was
#' dropped by an intervening subsetting operation), the function stops with an
#' informative error.
#'
#' The `lda_blocks` argument controls which trait sets are analyzed:
#'
#' \itemize{
#'   \item `"all"`: runs LDA on each individual block in `base_cols`
#'   (typically `"veg"`, `"flo"`, `"fru"`), as well as on every pairwise and
#'   full combination of these blocks (`"vegflo"`, `"vegfru"`, `"flofru"`,
#'   `"vegflofru"`). Combinations whose source blocks are absent from
#'   `base_cols` are silently skipped;
#'   \item a character vector naming one or more blocks in `base_cols`: runs
#'   LDA only on the requested block(s), individually.
#' }
#'
#' The grouping variable discriminated between is read from the column of
#' `analysis_data` named in `lda_group` — by default `"taxon"`, the grouping
#' column that [morph_matrix_setting()] always attaches. Passing any other
#' name only works if that column has been added to `analysis_data`
#' beforehand; `morph_lda()` does not read from any other object. Within
#' each block, specimens with a missing value for `lda_group` are removed,
#' as are specimens with no observation for any trait in the block; group
#' levels with fewer than `lda_min_n` remaining specimens are then dropped
#' (with a message naming them), since a group needs enough specimens to
#' contribute a stable within-group covariance estimate. The block is
#' skipped with a warning if fewer than 2 group levels remain after this
#' filtering, or if no trait varies among the retained specimens.
#'
#' LDA requires complete cases. Remaining missing values are imputed
#' trait-by-trait according to `lda_impute`: `"mean"` (default) fills a
#' trait's missing values with that trait's overall mean across all retained
#' specimens, matching the behaviour of the original analysis script;
#' `"group_mean"` fills them with the mean of the specimen's own group
#' instead, which is usually preferable for discriminant analysis since a
#' trait that is systematically higher or lower in one group is not pulled
#' toward the grand mean for the specimens missing it — but note this uses
#' group identity to fill in trait values before fitting a model that
#' predicts group identity from those same traits, so it can make group
#' separation look better than it would from complete data alone, and is
#' best treated with that caveat when a block has substantial missingness.
#' If a group has no non-missing values for a trait at all, that trait's
#' overall mean is used for that group's specimens instead, with a warning.
#' Traits that are constant among the retained specimens after imputation
#' are then dropped, since they cannot contribute to group separation; the
#' remaining traits are centered and scaled with [base::scale()] before
#' fitting, so that loadings are comparable across traits measured in
#' different units.
#'
#' Taxon colours are taken from the `"taxon_colors"` attribute attached to
#' `analysis_data` by [morph_matrix_setting()] when `lda_group = "taxon"`,
#' so that groups keep the same colour used in other plots (e.g.
#' [morph_pca()], [morph_kmeans()]). If this attribute is missing, lacks an
#' entry for a group present in a given block, or `lda_group` is not
#' `"taxon"`, colours are generated automatically from the `viridis`
#' palette for the groups analyzed, with a warning in the first two cases.
#'
#' `lda_lds` names the discriminant axes to plot (default `c("LD1", "LD2")`).
#' If a block yields fewer discriminant functions than requested — which
#' happens whenever the number of retained groups or traits is small, since
#' at most `min(n_groups - 1, n_traits)` axes exist — the function falls
#' back to the first two available axes (or the only one) and reports this.
#' When only one discriminant axis is available (the common case of exactly
#' two groups), a scatterplot with ellipses cannot be drawn in two
#' dimensions; the function instead saves a one-dimensional density plot of
#' that axis, with one curve per group and a rug of individual scores, so
#' the block still produces a usable figure rather than a degenerate one.
#'
#' Output files are written to a date-specific directory inside `Figs.LDA`.
#' For each analyzed block, a scores scatterplot (or density plot, in the
#' single-axis case), a loadings bar plot for the first discriminant axis,
#' and — when `save_xlsx = TRUE` — spreadsheets of specimen scores and
#' trait loadings, are saved with names prefixed `"LDA_plot_"`,
#' `"LDA_loadings_"` and `"LDA_scores_"`/`"LDA_loadings_"` respectively and
#' suffixed by the block name and `lda_group`. A single `"LDA_summary.xlsx"`
#' combining the summary rows of all analyzed blocks is written once every
#' block has been processed.
#'
#' @param analysis_data A numeric data frame containing morphological traits
#' as columns and specimens as rows, with `id` and `taxon` columns and
#' specimen identifiers as row names. Typically the object returned by
#' [morph_matrix_setting()]. If it carries a `"taxon_colors"` attribute and
#' `lda_group = "taxon"`, those colours are reused; otherwise a `viridis`
#' palette is generated automatically.
#'
#' @param lda_blocks Character. Either `"all"`, to run LDA on every
#' individual block in `base_cols` plus all pairwise and full combinations of
#' these blocks, or a character vector naming one or more blocks in
#' `base_cols` to analyze individually. Default is `"all"`.
#'
#' @param lda_group Character string naming the column of `analysis_data` to
#' discriminate between. Default is `"taxon"`.
#'
#' @param lda_lds Character vector of length 1 or 2 naming the discriminant
#' axes to plot, e.g. `c("LD1", "LD2")` (default). Falls back to the first
#' available axis or axes when a block does not have enough discriminant
#' functions for what is requested.
#'
#' @param lda_ellipse Numeric between 0 and 1. Confidence level of the
#' ellipses drawn on the two-dimensional scores plot. Default is `0.95`.
#' Unused when only one discriminant axis is available for a block.
#'
#' @param lda_min_n Integer. Minimum number of specimens a group must have,
#' within a given block, to be retained. Groups with fewer specimens are
#' dropped before fitting. Default is `3`.
#'
#' @param lda_impute Character string, `"mean"` (default) or `"group_mean"`,
#' controlling how missing trait values are imputed before fitting (see
#' Details).
#'
#' @param save_xlsx Logical. If `TRUE` (default), writes spreadsheets of
#' specimen scores, trait loadings, and the combined block summary to the
#' output directory.
#'
#' @param file_formats Character vector specifying which image formats to
#' save. One or more of `"pdf"` and `"jpeg"`. Default is
#' `c("pdf", "jpeg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' dropped groups, per-block summaries, skipped blocks, and the output
#' directory to the console.
#'
#' @return
#' Invisibly returns a summary data frame with one row per analyzed block, or
#' `NULL` if no block could be analyzed. The columns are:
#' \describe{
#'   \item{`block`}{Name of the trait block.}
#'   \item{`group`}{The grouping variable discriminated between (`lda_group`).}
#'   \item{`n`}{Number of specimens entering the LDA.}
#'   \item{`n_groups`}{Number of group levels retained.}
#'   \item{`n_lds`}{Number of discriminant functions available.}
#'   \item{`ld1_var`}{Percentage of between-group variance explained by the
#'   first discriminant axis.}
#' }
#'
#' The function is also called for its side effects: scores and loadings
#' plots, and (when `save_xlsx = TRUE`) spreadsheets, are written to disk for
#' each analyzed trait block.
#'
#' @seealso
#' [morph_matrix_setting()],
#' [morph_pca()],
#' [morph_kmeans()],
#' [morph_nmm()],
#' [morph_anova()],
#' [morph_kruskal()]
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
#' # LDA of taxon identity for every individual block and every combination
#' lda_summary <- morph_lda(
#'   analysis_data = res,
#'   lda_blocks = "all"
#' )
#'
#' lda_summary
#'
#' # Vegetative block only, with group-mean imputation
#' morph_lda(
#'   analysis_data = res,
#'   lda_blocks = "veg",
#'   lda_impute = "group_mean"
#' )
#'
#' # Discriminate a different grouping column added beforehand
#' res$locality <- specimen_localities[rownames(res)]
#' morph_lda(
#'   analysis_data = res,
#'   lda_group = "locality"
#' )
#' }
#'
#' @importFrom stats setNames var predict
#' @importFrom MASS lda
#' @importFrom ggplot2 aes coord_flip element_blank element_text geom_col
#'   geom_density geom_hline geom_point geom_rug geom_vline ggplot labs
#'   scale_color_manual scale_fill_manual stat_ellipse theme theme_bw
#' @importFrom cowplot save_plot
#' @importFrom openxlsx write.xlsx
#' @importFrom viridis viridis
#' @importFrom rlang .data
#'
#' @export

morph_lda <- function(analysis_data,
                      lda_blocks = "all",
                      lda_group = "taxon",
                      lda_lds = c("LD1", "LD2"),
                      lda_ellipse = 0.95,
                      lda_min_n = 3,
                      lda_impute = c("mean", "group_mean"),
                      save_xlsx = TRUE,
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

if (!is.character(lda_blocks) || length(lda_blocks) < 1) {
  stop("`lda_blocks` must be \"all\" or a character vector of block names.")
}

if (!is.character(lda_group) || length(lda_group) != 1L) {
  stop("`lda_group` must be a single column name in `analysis_data`.")
}

if (!lda_group %in% names(analysis_data)) {
  stop("Group column '", lda_group, "' not found in `analysis_data`.")
}

if (!is.character(lda_lds) || length(lda_lds) < 1 || length(lda_lds) > 2) {
  stop("`lda_lds` must be a character vector of length 1 or 2, e.g. c(\"LD1\", \"LD2\").")
}

if (!is.numeric(lda_ellipse) ||
    length(lda_ellipse) != 1L ||
    lda_ellipse <= 0 ||
    lda_ellipse >= 1) {
  stop("`lda_ellipse` must be a single number between 0 and 1.")
}

if (!is.numeric(lda_min_n) || length(lda_min_n) != 1L || lda_min_n < 2) {
  stop("`lda_min_n` must be a single integer >= 2.")
}

lda_impute <- match.arg(lda_impute)

if (!is.logical(save_xlsx) || length(save_xlsx) != 1L) {
  stop("`save_xlsx` must be TRUE or FALSE.")
}

file_formats <- match.arg(
  file_formats,
  choices = c("pdf", "jpeg"),
  several.ok = TRUE
)

if (!is.logical(verbose) || length(verbose) != 1L) {
  stop("`verbose` must be TRUE or FALSE.")
}

# ---- Species vector (always taxon, regardless of lda_group) ------------------

species_all <- stats::setNames(
  as.character(analysis_data$taxon),
  rownames(analysis_data)
)

# ---- Output folder ------------------------------------------------------------

output_dir <- file.path("Figs.LDA", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (verbose) {
  message("Running LDA")
}

# ---- Resolve blocks ----------------------------------------------------------

if (identical(lda_blocks, "all")) {
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
  requested <- intersect(lda_blocks, names(base_cols))
  if (length(requested) == 0) {
    stop(
      "No valid blocks in `lda_blocks`. Use \"all\" or any of: ",
      paste(names(base_cols), collapse = ", ")
    )
  }
  run_blocks <- base_cols[requested]
}

if (verbose) {
  message("  Blocks to run: ", paste(names(run_blocks), collapse = ", "))
}

# ---- Helper: trait display label ---------------------------------------------

trait_label <- function(x) {
  gsub("_", " ", gsub("/.*", "", x))
}

# ---- Helper: discriminant axis index and % variance label --------------------

ld_index <- function(ld_name) as.integer(gsub("^LD", "", ld_name))

ld_label <- function(ld_name, prop_var) {
  paste0(ld_name, " (", round(prop_var[ld_index(ld_name)], 1), "%)")
}

# ---- Helper: taxon/group colours ----------------------------------------------

resolve_colors <- function(levels_present) {

color_vector <- if (identical(lda_group, "taxon")) {
  attr(analysis_data, "taxon_colors")
} else {
  NULL
}

if (is.null(color_vector)) {
  return(stats::setNames(viridis::viridis(length(levels_present)), levels_present))
}

missing_lv <- setdiff(levels_present, names(color_vector))

if (length(missing_lv) > 0) {
  warning(
    "No stored colour found for: ", paste(missing_lv, collapse = ", "),
    ". Assigning automatic colours for these groups."
  )
  color_vector[missing_lv] <- viridis::viridis(length(missing_lv))
}

color_vector[levels_present]
}

# ---- Helper: impute missing trait values --------------------------------------

impute_traits <- function(trait_data, group) {

as.data.frame(lapply(trait_data, function(x) {

if (!anyNA(x)) return(x)

if (lda_impute == "mean") {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}

# "group_mean": fill within each group; fall back to the overall mean
# for a group with no non-missing values for this trait at all
overall_mean <- mean(x, na.rm = TRUE)

for (lv in levels(group)) {
  idx <- which(group == lv & is.na(x))
  if (length(idx) == 0) next

  lv_mean <- mean(x[group == lv], na.rm = TRUE)

  if (is.nan(lv_mean)) {
    warning(
      "Group '", lv, "' has no non-missing values for one trait; ",
      "using the overall mean for its specimens instead."
    )
    lv_mean <- overall_mean
  }

  x[idx] <- lv_mean
}

x
}))
}

# ---- Helper: drop zero-variance traits -----------------------------------------

drop_constant <- function(df) {
  keep <- vapply(df, function(x) {
    v <- stats::var(x, na.rm = TRUE)
    !is.na(v) && v > 0
  }, logical(1))

  df[, keep, drop = FALSE]
}

# ---- Helper: save a plot in the requested formats -------------------------------

save_outputs <- function(filename_stub, plot_obj, base_height, base_aspect_ratio) {
  if ("pdf" %in% file_formats) {
    cowplot::save_plot(
      file.path(output_dir, paste0(filename_stub, ".pdf")),
      plot_obj,
      ncol = 1,
      nrow = 1,
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
      nrow = 1,
      base_height = base_height,
      base_aspect_ratio = base_aspect_ratio,
      base_width = NULL
    )
  }
}

# ---- Helper: LDA for one trait block --------------------------------------------

run_lda_block <- function(block_name, block_columns) {

valid_cols <- intersect(block_columns, names(analysis_data))
block_data <- analysis_data[, valid_cols, drop = FALSE]

# Grouping variable, read directly from analysis_data so it stays aligned
# with block_data by specimen ID (row name) regardless of row order
block_data[[lda_group]] <- as.factor(analysis_data[rownames(block_data), lda_group])

# Drop specimens with a missing group or with no observation for any trait
block_data <- block_data[!is.na(block_data[[lda_group]]), , drop = FALSE]
block_data <- block_data[rowSums(!is.na(block_data[, valid_cols, drop = FALSE])) > 0, , drop = FALSE]

# Drop group levels with too few remaining specimens
group_counts <- table(block_data[[lda_group]])
small_groups <- names(group_counts)[group_counts < lda_min_n]

if (length(small_groups) > 0) {
  if (verbose) {
    message(
      "  Removing groups with fewer than ", lda_min_n, " specimens: ",
      paste(small_groups, collapse = ", ")
    )
  }
  block_data <- block_data[!block_data[[lda_group]] %in% small_groups, , drop = FALSE]
  block_data[[lda_group]] <- droplevels(block_data[[lda_group]])
}

n_groups <- length(unique(block_data[[lda_group]]))

if (n_groups < 2) {
  warning("Skipping '", block_name, "': fewer than 2 groups remaining.")
  return(NULL)
}

if (verbose) {
  message("Running LDA for: ", block_name)
}

# ---- Impute, then drop constant traits -------------------------------------

trait_data <- impute_traits(block_data[, valid_cols, drop = FALSE], block_data[[lda_group]])
trait_data <- drop_constant(trait_data)

if (ncol(trait_data) < 1) {
  warning("Skipping '", block_name, "': no variable traits after filtering.")
  return(NULL)
}

trait_scaled <- as.data.frame(scale(trait_data))

if (verbose) {
  message(
    "  Specimens: ", nrow(trait_scaled),
    " | Traits: ", ncol(trait_scaled),
    " | Groups: ", n_groups
  )
}

# ---- Fit LDA ------------------------------------------------------------------

lda_model <- tryCatch(
  MASS::lda(trait_scaled, grouping = block_data[[lda_group]]),
  error = function(e) {
    warning("LDA failed for '", block_name, "': ", conditionMessage(e))
    NULL
  }
)

if (is.null(lda_model)) return(NULL)

lda_result <- stats::predict(lda_model)

# ---- Plot data -----------------------------------------------------------------

lda_df <- as.data.frame(lda_result$x)
lda_df$group <- block_data[[lda_group]]
lda_df$taxon <- species_all[rownames(lda_df)]
lda_df$id <- rownames(lda_df)

prop_var <- lda_model$svd^2 / sum(lda_model$svd^2) * 100

available_lds <- colnames(lda_result$x)
missing_lds <- setdiff(lda_lds, available_lds)

if (length(missing_lds) > 0) {
  if (verbose) {
    message(
      "  Requested LD(s) not available: ", paste(missing_lds, collapse = ", "),
      ". Available: ", paste(available_lds, collapse = ", "),
      ". Using first available axis/axes."
    )
  }
  plot_lds <- available_lds[seq_len(min(2, length(available_lds)))]
} else {
  plot_lds <- lda_lds
}

if (verbose) {
  message("  Group sizes: ", paste(table(lda_df$group), collapse = " | "))
  message("  Discriminant functions available: ", length(available_lds))
  for (i in seq_along(prop_var)) {
    message(sprintf("  LD%d: %.1f%% variance explained", i, prop_var[i]))
  }
}

# ---- Colours -----------------------------------------------------------------

color_vector <- resolve_colors(sort(levels(lda_df$group)))

# ---- Scores plot -----------------------------------------------------------

fname_suffix <- paste0(block_name, "_", lda_group)

if (length(plot_lds) >= 2) {

  xld <- plot_lds[1]
  yld <- plot_lds[2]

  p_lda <- ggplot2::ggplot(
    lda_df,
    ggplot2::aes(x = .data[[xld]], y = .data[[yld]])
  ) +
    ggplot2::stat_ellipse(
      ggplot2::aes(color = .data[["group"]]),
      level = lda_ellipse,
      type = "norm",
      linetype = "dashed",
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data[["group"]]),
      shape = 21, size = 3, colour = "white", alpha = 0.8
    ) +
    ggplot2::scale_fill_manual(values = color_vector) +
    ggplot2::scale_color_manual(values = color_vector) +
    ggplot2::geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, color = "gray40", linetype = "dashed") +
    ggplot2::labs(
      x = ld_label(xld, prop_var),
      y = ld_label(yld, prop_var),
      fill = lda_group,
      title = paste0("LDA - ", block_name)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(face = "italic")
    )

  save_outputs(
    paste0("LDA_plot_", fname_suffix),
    p_lda,
    base_height = 7,
    base_aspect_ratio = 1.3
  )

} else {

  # Only one discriminant axis available (e.g. exactly 2 groups): a 2-D
  # ellipse plot is degenerate here, so show a 1-D density instead.
  xld <- plot_lds[1]

  p_lda <- ggplot2::ggplot(
    lda_df,
    ggplot2::aes(x = .data[[xld]])
  ) +
    ggplot2::geom_density(
      ggplot2::aes(fill = .data[["group"]], color = .data[["group"]]),
      alpha = 0.4
    ) +
    ggplot2::geom_rug(
      ggplot2::aes(color = .data[["group"]]),
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = color_vector) +
    ggplot2::scale_color_manual(values = color_vector) +
    ggplot2::geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
    ggplot2::labs(
      x = ld_label(xld, prop_var),
      y = "Density",
      fill = lda_group,
      title = paste0("LDA - ", block_name, " (single discriminant axis)")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  save_outputs(
    paste0("LDA_plot_", fname_suffix),
    p_lda,
    base_height = 5,
    base_aspect_ratio = 1.6
  )
}

# ---- Loadings plot: trait contributions to LD1 ------------------------------

loadings_df <- as.data.frame(lda_model$scaling)
loadings_df$trait_raw <- rownames(loadings_df)
loadings_df$trait <- trait_label(loadings_df$trait_raw)
loadings_df <- loadings_df[order(abs(loadings_df$LD1), decreasing = TRUE), ]
loadings_df$trait <- factor(loadings_df$trait, levels = rev(loadings_df$trait))

p_load <- ggplot2::ggplot(
  loadings_df,
  ggplot2::aes(
    x = .data[["trait"]],
    y = .data[["LD1"]],
    fill = .data[["LD1"]] > 0
  )
) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(
    values = c("TRUE" = "#1a9850", "FALSE" = "#d73027"),
    guide = "none"
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL, y = "LD1 loading",
    title = paste0("LDA loadings - ", block_name)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 7)
  )

save_outputs(
  paste0("LDA_loadings_", fname_suffix),
  p_load,
  base_height = max(4, nrow(loadings_df) * 0.3),
  base_aspect_ratio = 1.5
)

# ---- Save spreadsheets ------------------------------------------------------

if (save_xlsx) {

  scores_out <- lda_df[, c("id", "taxon", "group", available_lds), drop = FALSE]
  names(scores_out)[names(scores_out) == "group"] <- lda_group

  openxlsx::write.xlsx(
    scores_out,
    file.path(output_dir, paste0("LDA_scores_", fname_suffix, ".xlsx"))
  )

  openxlsx::write.xlsx(
    loadings_df[, c("trait_raw", "trait", setdiff(names(loadings_df), c("trait_raw", "trait")))],
    file.path(output_dir, paste0("LDA_loadings_", fname_suffix, ".xlsx"))
  )
}

if (verbose) {
  message("  Done: ", block_name)
}

# ---- Block summary -----------------------------------------------------------

data.frame(
  block = block_name,
  group = lda_group,
  n = nrow(lda_df),
  n_groups = n_groups,
  n_lds = length(available_lds),
  ld1_var = round(prop_var[1], 2),
  row.names = NULL,
  stringsAsFactors = FALSE
)
}

# ---- Run all blocks -------------------------------------------------------------

all_results <- Filter(
  Negate(is.null),
  mapply(
    run_lda_block,
    names(run_blocks),
    run_blocks,
    SIMPLIFY = FALSE
  )
)

if (length(all_results) == 0) {
  warning("No trait block could be analyzed.")
  return(invisible(NULL))
}

lda_summary <- do.call(rbind, all_results)
rownames(lda_summary) <- NULL

if (save_xlsx) {
  openxlsx::write.xlsx(lda_summary, file.path(output_dir, "LDA_summary.xlsx"))
}

if (verbose) {
  message("LDA analysis done. Outputs saved in: ", output_dir)
}

invisible(lda_summary)
}
