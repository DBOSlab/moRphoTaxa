#' One-way ANOVA of morphological traits across a grouping variable
#'
#' @description
#' Runs a one-way ANOVA for every morphological trait in each requested trait
#' block, testing whether trait means differ across the levels of one or more
#' grouping variables, and generates a summary bar plot and spreadsheet for
#' each block x group combination. The function is designed to work directly
#' on the objects produced by [morph_matrix_setting()].
#'
#' For each trait block and each grouping variable, the function fits
#' `trait ~ group` with [stats::aov()] for every trait in the block,
#' collects the F statistic and p-value, adjusts p-values for multiple
#' testing across the traits tested in that block, and saves a horizontal bar
#' plot of F statistics (coloured by significance) together with a
#' spreadsheet of the results.
#'
#' @details
#' `morph_anova()` expects `analysis_data` as produced by
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
#' Grouping variables are read from columns of `analysis_data` itself, named
#' in `anova_groups` — by default just `"taxon"`, the grouping column that
#' [morph_matrix_setting()] always attaches. Passing any other name only
#' works if that column has been added to `analysis_data` beforehand (e.g. a
#' locality or growth-form column merged in before calling this function);
#' `morph_anova()` does not read from any other object.
#'
#' The `anova_blocks` argument controls which trait sets are analyzed:
#'
#' \itemize{
#'   \item `"all"`: runs ANOVA on each individual block in `base_cols`
#'   (typically `"veg"`, `"flo"`, `"fru"`), as well as on every pairwise and
#'   full combination of these blocks (`"vegflo"`, `"vegfru"`, `"flofru"`,
#'   `"vegflofru"`). Combinations whose source blocks are absent from
#'   `base_cols` are silently skipped;
#'   \item a character vector naming one or more blocks in `base_cols`: runs
#'   ANOVA only on the requested block(s), individually.
#' }
#'
#' Within each block x group combination, specimens with a missing value for
#' the grouping variable are removed; the combination is skipped with a
#' warning if fewer than 3 specimens or fewer than 2 group levels remain.
#' Each trait is then tested individually: rows with a missing value for that
#' particular trait or the grouping variable are omitted by [stats::aov()]'s
#' default `na.action`, and a trait is skipped (silently, since this is
#' expected for traits scored only in part of a block) if it has fewer than 3
#' non-missing values. A trait whose model fails to fit — for instance
#' because it is constant within every group — is also skipped, since no F
#' statistic exists in that case. This is an omnibus per-trait test: it
#' indicates whether group means differ but does not identify which pairs of
#' groups drive that difference; TukeyHSD or another post-hoc procedure can be
#' run on the individual `aov` fits if that is needed.
#'
#' P-values are adjusted with [stats::p.adjust()] using `anova_p_adjust`,
#' separately within each block x group combination — that is, the
#' correction accounts for the number of traits tested in that block, not for
#' the total number of tests across all blocks and groups in the call. A
#' trait is flagged `significant` when its adjusted p-value is below
#' `anova_alpha`.
#'
#' Output files are written to a date-specific directory inside `Figs.ANOVA`.
#' For each block x group combination, a bar plot of F statistics (traits
#' sorted from highest to lowest, coloured by significance, with a reference
#' line at F = 1) is saved in the formats specified by `file_formats`, named
#' `"ANOVA_<block>_<group>"`. When `save_xlsx = TRUE`, a spreadsheet with
#' that combination's results is written alongside it, and a single
#' `"ANOVA_all_results.xlsx"` combining every block x group combination is
#' written once all combinations have been processed.
#'
#' @param analysis_data A numeric data frame containing morphological traits
#' as columns and specimens as rows, with `id` and `taxon` columns and
#' specimen identifiers as row names. Typically the object returned by
#' [morph_matrix_setting()].
#'
#' @param anova_blocks Character. Either `"all"`, to test every individual
#' block in `base_cols` plus all pairwise and full combinations of these
#' blocks, or a character vector naming one or more blocks in `base_cols` to
#' test individually. Default is `"all"`.
#'
#' @param anova_groups Character vector of column names in `analysis_data` to
#' use as grouping variables, one ANOVA run per name. Default is `"taxon"`.
#'
#' @param anova_p_adjust Character. Method passed to [stats::p.adjust()] for
#' correcting p-values across the traits tested within each block x group
#' combination. One of `"bonferroni"`, `"holm"`, `"hochberg"`, `"hommel"`,
#' `"BH"` (default), `"BY"`, or `"none"`.
#'
#' @param anova_alpha Numeric between 0 and 1. Adjusted-p-value threshold
#' below which a trait is flagged as significant. Default is `0.05`.
#'
#' @param save_xlsx Logical. If `TRUE` (default), writes a spreadsheet of
#' results for each block x group combination, plus a combined spreadsheet
#' of every combination, to the output directory.
#'
#' @param file_formats Character vector specifying which image formats to
#' save. One or more of `"pdf"` and `"jpeg"`. Default is
#' `c("pdf", "jpeg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' per-combination counts of traits tested and found significant, skipped
#' combinations, and the output directory to the console.
#'
#' @return
#' Invisibly returns a data frame combining the results of every analyzed
#' block x group combination, or `NULL` if none could be analyzed. The
#' columns are:
#' \describe{
#'   \item{`trait`}{Name of the morphological trait.}
#'   \item{`F_statistic`}{F statistic of the one-way ANOVA for that trait.}
#'   \item{`p_value`}{Unadjusted p-value.}
#'   \item{`p_adjusted`}{P-value adjusted with `anova_p_adjust`, within the
#'   block x group combination.}
#'   \item{`significant`}{`TRUE` when `p_adjusted < anova_alpha`.}
#'   \item{`block`}{Name of the trait block.}
#'   \item{`group`}{Name of the grouping variable.}
#' }
#'
#' The function is also called for its side effects: bar plots and (when
#' `save_xlsx = TRUE`) spreadsheets of results are written to disk for each
#' analyzed block x group combination.
#'
#' @seealso
#' [morph_matrix_setting()],
#' [morph_pca()],
#' [morph_kmeans()],
#' [morph_nmm()]
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
#' # Test every individual block and every block combination against taxon
#' anova_res <- morph_anova(
#'   analysis_data = res,
#'   anova_blocks = "all"
#' )
#'
#' # Traits that differ significantly among taxa
#' subset(anova_res, significant)
#'
#' # Vegetative block only, with a stricter correction
#' morph_anova(
#'   analysis_data = res,
#'   anova_blocks = "veg",
#'   anova_p_adjust = "bonferroni"
#' )
#'
#' # Against a custom grouping column added beforehand
#' res$locality <- specimen_localities[rownames(res)]
#' morph_anova(
#'   analysis_data = res,
#'   anova_groups = "locality"
#' )
#' }
#'
#' @importFrom stats aov p.adjust as.formula
#' @importFrom ggplot2 aes element_blank element_text geom_col geom_hline
#'   ggplot coord_flip labs scale_fill_manual theme theme_bw
#' @importFrom cowplot save_plot
#' @importFrom openxlsx write.xlsx
#' @importFrom rlang .data
#'
#' @export

morph_anova <- function(analysis_data,
                        anova_blocks = "all",
                        anova_groups = "taxon",
                        anova_p_adjust = "BH",
                        anova_alpha = 0.05,
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

if (!is.character(anova_blocks) || length(anova_blocks) < 1) {
  stop("`anova_blocks` must be \"all\" or a character vector of block names.")
}

if (!is.character(anova_groups) || length(anova_groups) < 1) {
  stop("`anova_groups` must be a character vector of column names in `analysis_data`.")
}

missing_groups <- setdiff(anova_groups, names(analysis_data))

if (length(missing_groups) > 0) {
  stop(
    "Group column(s) not found in `analysis_data`: ",
    paste(missing_groups, collapse = ", ")
  )
}

anova_p_adjust <- match.arg(anova_p_adjust, choices = stats::p.adjust.methods)

if (!is.numeric(anova_alpha) ||
    length(anova_alpha) != 1L ||
    anova_alpha <= 0 ||
    anova_alpha >= 1) {
  stop("`anova_alpha` must be a single number between 0 and 1.")
}

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

# ---- Output folder ------------------------------------------------------------

output_dir <- file.path("Figs.ANOVA", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (verbose) {
  message("Running ANOVA")
}

# ---- Resolve blocks ----------------------------------------------------------

if (identical(anova_blocks, "all")) {
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
  requested <- intersect(anova_blocks, names(base_cols))
  if (length(requested) == 0) {
    stop(
      "No valid blocks in `anova_blocks`. Use \"all\" or any of: ",
      paste(names(base_cols), collapse = ", ")
    )
  }
  run_blocks <- base_cols[requested]
}

if (verbose) {
  message("Blocks to run:", paste(names(run_blocks), collapse = ", "))
}

# ---- Helper: save a plot in the requested formats ----------------------------

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

# ---- Helper: ANOVA for one block x group combination -------------------------

run_anova_block <- function(block_name, block_columns, group_var) {

valid_cols <- intersect(block_columns, names(analysis_data))
block_data <- analysis_data[, valid_cols, drop = FALSE]

# Grouping variable, read directly from analysis_data so it stays aligned
# with block_data by specimen ID (row name) regardless of row order
block_data$group__ <- as.factor(analysis_data[rownames(block_data), group_var])

# Drop specimens with a missing group
block_data <- block_data[!is.na(block_data$group__), , drop = FALSE]

if (nrow(block_data) < 3 || length(unique(block_data$group__)) < 2) {
  warning(
    "Skipping '", block_name, "' x '", group_var,
    "': fewer than 3 specimens or fewer than 2 group levels after removing ",
    "missing groups."
  )
  return(NULL)
}

if (verbose) {
  message("Running ANOVA for: ", block_name, " x ", group_var)
}

trait_cols <- setdiff(names(block_data), "group__")

# ---- Run ANOVA per trait ---------------------------------------------------

results <- lapply(trait_cols, function(trait) {
  y <- block_data[[trait]]
  if (sum(!is.na(y)) < 3) return(NULL)

  tryCatch({
    fit <- stats::aov(
      stats::as.formula(paste0("`", trait, "` ~ group__")),
      data = block_data
    )
    sm <- summary(fit)[[1]]
    data.frame(
      trait = trait,
      F_statistic = sm[["F value"]][1],
      p_value = sm[["Pr(>F)"]][1]
    )
  }, error = function(e) NULL)
})

results <- do.call(rbind, Filter(Negate(is.null), results))

if (is.null(results) || nrow(results) == 0) {
  warning(
    "Skipping '", block_name, "' x '", group_var,
    "': no trait had enough data for a valid ANOVA."
  )
  return(NULL)
}

# ---- Adjust p-values -------------------------------------------------------

results$p_adjusted <- stats::p.adjust(results$p_value, method = anova_p_adjust)
results$significant <- results$p_adjusted < anova_alpha
results$block <- block_name
results$group <- group_var

n_sig <- sum(results$significant, na.rm = TRUE)

if (verbose) {
  message(
    "  Traits tested: ", nrow(results),
    " | Significant (adj. p < ", anova_alpha, "): ", n_sig
  )
}

# ---- Plot: F-statistic bar chart, coloured by significance -----------------

plot_df <- results[order(results$F_statistic, decreasing = TRUE), ]
plot_df$trait <- factor(plot_df$trait, levels = plot_df$trait)
plot_df$sig_label <- ifelse(plot_df$significant, "sig", "n.s.")

p_bar <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = .data[["trait"]],
    y = .data[["F_statistic"]],
    fill = .data[["sig_label"]]
  )
) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(
    values = c("n.s." = "gray70", "sig" = "#2ca25f"),
    labels = c(
      "n.s." = "n.s.",
      "sig" = paste0("adj. p < ", anova_alpha)
    )
  ) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL,
    y = "F statistic",
    fill = "",
    title = paste0("ANOVA - ", block_name, " | ", group_var)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 7)
  )

fname <- paste0("ANOVA_", block_name, "_", group_var)

save_outputs(
  fname,
  p_bar,
  base_height = max(4, nrow(plot_df) * 0.3),
  base_aspect_ratio = 1.5
)

# ---- Save spreadsheet -------------------------------------------------------

if (save_xlsx) {
  openxlsx::write.xlsx(
    results[, c("trait", "F_statistic", "p_value", "p_adjusted", "significant")],
    file.path(output_dir, paste0(fname, ".xlsx"))
  )
}

if (verbose) {
  message("  Done: ", block_name, " x ", group_var)
}

results
}

# ---- Run all block x group combinations ---------------------------------------

all_results <- list()

for (group_var in anova_groups) {
  for (block_name in names(run_blocks)) {
    key <- paste0(block_name, "_", group_var)
    all_results[[key]] <- run_anova_block(block_name, run_blocks[[block_name]], group_var)
  }
}

combined <- do.call(rbind, Filter(Negate(is.null), all_results))

if (is.null(combined) || nrow(combined) == 0) {
  warning("No block x group combination could be analyzed.")
  return(invisible(NULL))
}

rownames(combined) <- NULL

if (save_xlsx) {
  openxlsx::write.xlsx(combined, file.path(output_dir, "ANOVA_all_results.xlsx"))
}

if (verbose) {
  message("ANOVA analysis done. Outputs saved in: ", output_dir)
}

invisible(combined)
}
