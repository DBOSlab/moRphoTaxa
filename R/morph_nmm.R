#' Model-based clustering of morphological trait blocks with normal mixture models
#'
#' @description
#' Fits normal mixture models (NMM) to morphological trait blocks derived from
#' a morphometric matrix, and generates a standard set of diagnostic figures
#' and tables for each trait block analyzed. The function is designed to work
#' directly on the objects produced by [morph_matrix_setting()].
#'
#' For each requested trait block, the function fits a PCA to rank traits by
#' their contribution to the first principal component, retains the top
#' contributing traits, fits a Gaussian mixture model with
#' [mclust::Mclust()] — which selects both the covariance structure and the
#' number of components by BIC — and saves a BIC plot, a classification plot,
#' an uncertainty plot, and a spreadsheet of group assignments.
#'
#' @details
#' `morph_nmm()` expects `analysis_data` as produced by
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
#' The `nmm_blocks` argument controls which trait sets are analyzed:
#'
#' \itemize{
#'   \item `"all"`: runs NMM on each individual block in `base_cols`
#'   (typically `"veg"`, `"flo"`, `"fru"`), as well as on every pairwise and
#'   full combination of these blocks (`"vegflo"`, `"vegfru"`, `"flofru"`,
#'   `"vegflofru"`). Combinations whose source blocks are absent from
#'   `base_cols` are silently skipped;
#'   \item a character vector naming one or more blocks in `base_cols`: runs
#'   NMM only on the requested block(s), individually.
#' }
#'
#' Within each block, rows with missing values are removed via
#' [stats::na.omit()] and traits with zero variance across the retained
#' specimens are dropped, since they cannot contribute to a correlation-based
#' PCA nor to a mixture model. Blocks left with fewer than 2 traits, or with
#' fewer than `nmm_min_n` complete specimens, are skipped with a warning.
#'
#' Trait selection is based on PCA: a PCA is fit with [stats::prcomp()] using
#' `scale. = TRUE`, the squared loadings on the first principal component are
#' rescaled to percentages, and the `nmm_top_traits` highest-contributing
#' traits are retained for the mixture model. This step keeps the
#' dimensionality of the mixture model manageable relative to the number of
#' specimens; note that it is unsupervised and makes no use of the `taxon`
#' column, so the resulting groups remain independent of the current
#' taxonomic hypothesis.
#'
#' The mixture model is fit with [mclust::Mclust()]. By default the number of
#' components and the covariance parameterization are both selected by BIC
#' over the full set of `mclust` models; `nmm_G` and `nmm_model_names`
#' restrict that search when a particular number of groups or a particular
#' family of models is required. Model fitting is wrapped in [tryCatch()], so
#' a block whose model fails to converge produces a warning and is skipped
#' rather than aborting the run. `nmm_seed` is applied before each block for
#' reproducibility.
#'
#' For each analyzed block, the following figures are generated:
#'
#' \itemize{
#'   \item a BIC plot across covariance models and numbers of components
#'   ([factoextra::fviz_mclust_bic()]);
#'   \item a classification plot showing specimens projected onto the first
#'   two principal components of the retained traits, coloured by mixture
#'   component;
#'   \item an uncertainty plot, in which point size reflects the uncertainty
#'   of each specimen's assignment.
#' }
#'
#' Confidence ellipses are drawn on the classification and uncertainty plots
#' only when more than one component is retained and every component contains
#' at least 4 specimens; otherwise the ellipses are suppressed and a message
#' is emitted, since they cannot be estimated from fewer points. Both plots
#' are additionally wrapped in [tryCatch()], so a plotting failure does not
#' interrupt the block.
#'
#' Output files are written to a date-specific directory inside `Figs.NMM`.
#' Figures are saved in the formats specified by `file_formats`, with names
#' prefixed by figure type (`"NMM_BIC_"`, `"NMM_classification_"`,
#' `"NMM_uncertainty_"`) and suffixed by the trait block name. When
#' `save_xlsx = TRUE`, a per-block spreadsheet (`"NMM_groups_<block>.xlsx"`)
#' holding the retained trait values, taxon, assigned group and assignment
#' uncertainty for every specimen is written to the same directory, along
#' with a single `"NMM_summary.xlsx"` combining the summary rows of all
#' analyzed blocks.
#'
#' @param analysis_data A numeric data frame containing morphological traits
#' as columns and specimens as rows, with `id` and `taxon` columns and
#' specimen identifiers as row names. Typically the object returned by
#' [morph_matrix_setting()].
#'
#' @param nmm_blocks Character. Either `"all"`, to run NMM on every individual
#' block in `base_cols` plus all pairwise and full combinations of these
#' blocks, or a character vector naming one or more blocks in `base_cols` to
#' analyze individually. Default is `"all"`.
#'
#' @param nmm_top_traits Integer. Number of highest-contributing traits on
#' principal component 1 retained for the mixture model. Default is `10`. If a
#' block contains fewer traits, all of them are used.
#'
#' @param nmm_min_n Integer. Minimum number of complete specimens required to
#' analyze a block. Blocks with fewer specimens, before or after subsetting to
#' the top traits, are skipped with a warning. Default is `10`.
#'
#' @param nmm_G Integer vector giving the numbers of mixture components to
#' evaluate, passed to [mclust::Mclust()]. If `NULL` (default), the `mclust`
#' default range is used and the number of components is selected by BIC.
#'
#' @param nmm_model_names Character vector of `mclust` model names (covariance
#' parameterizations) to evaluate, e.g. `c("EII", "VII", "VVV")`, passed to
#' [mclust::Mclust()]. If `NULL` (default), all applicable models are
#' evaluated.
#'
#' @param ellipse_level Numeric between 0 and 1. Confidence level of the
#' ellipses drawn on the classification plot. Default is `0.4`.
#'
#' @param nmm_seed Integer. Random seed applied before fitting each block, for
#' reproducible results. Default is `123`.
#'
#' @param save_xlsx Logical. If `TRUE` (default), writes a per-block
#' spreadsheet of specimen group assignments and a combined summary
#' spreadsheet to the output directory.
#'
#' @param file_formats Character vector specifying which image formats to
#' save. One or more of `"pdf"` and `"jpeg"`. Default is
#' `c("pdf", "jpeg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' the traits retained per block, per-block model summaries, skipped blocks,
#' and the output directory to the console.
#'
#' @return
#' Invisibly returns a summary data frame with one row per analyzed block, or
#' `NULL` if no block could be analyzed. The columns are:
#' \describe{
#'   \item{`block`}{Name of the trait block.}
#'   \item{`n`}{Number of specimens entering the mixture model.}
#'   \item{`model`}{Covariance parameterization selected by `mclust`.}
#'   \item{`k`}{Number of mixture components retained.}
#'   \item{`bic`}{BIC of the selected model.}
#'   \item{`top_traits`}{Retained traits, separated by `"|"`.}
#' }
#'
#' The function is also called for its side effects: BIC, classification and
#' uncertainty plots, and (when `save_xlsx = TRUE`) spreadsheets of group
#' assignments, are written to disk for each analyzed trait block.
#'
#' @seealso
#' [morph_matrix_setting()],
#' [morph_pca()],
#' [morph_kmeans()]
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
#' # Fit NMMs to every individual block and every block combination
#' nmm <- morph_nmm(
#'   analysis_data = res,
#'   nmm_blocks = "all"
#' )
#'
#' # Inspect the model selected for each block
#' nmm
#'
#' # Vegetative block only, using the 15 top PC1-contributing traits
#' morph_nmm(
#'   analysis_data = res,
#'   nmm_blocks = "veg",
#'   nmm_top_traits = 15
#' )
#'
#' # Restrict the search to 2-4 components, saving PDF output only and
#' # skipping the spreadsheets
#' morph_nmm(
#'   analysis_data = res,
#'   nmm_blocks = c("veg", "flo"),
#'   nmm_G = 2:4,
#'   save_xlsx = FALSE,
#'   file_formats = "pdf"
#' )
#' }
#'
#' @importFrom stats prcomp na.omit setNames var
#' @importFrom mclust Mclust mclustBIC
#' @importFrom cowplot save_plot
#' @importFrom ggplot2 labs theme_bw
#' @importFrom openxlsx write.xlsx
#'
#' @export

morph_nmm <- function(analysis_data,
                      nmm_blocks = "all",
                      nmm_top_traits = 10,
                      nmm_min_n = 10,
                      nmm_G = NULL,
                      nmm_model_names = NULL,
                      ellipse_level = 0.4,
                      nmm_seed = 123,
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

if (!is.character(nmm_blocks) || length(nmm_blocks) < 1) {
stop("`nmm_blocks` must be \"all\" or a character vector of block names.")
}

if (!is.numeric(nmm_top_traits) || length(nmm_top_traits) != 1L || nmm_top_traits < 2) {
stop("`nmm_top_traits` must be a single integer >= 2.")
}

if (!is.numeric(nmm_min_n) || length(nmm_min_n) != 1L || nmm_min_n < 3) {
stop("`nmm_min_n` must be a single integer >= 3.")
}

if (!is.null(nmm_G) && (!is.numeric(nmm_G) || any(nmm_G < 1))) {
stop("`nmm_G` must be NULL or a vector of positive integers.")
}

if (!is.null(nmm_model_names) && !is.character(nmm_model_names)) {
stop("`nmm_model_names` must be NULL or a character vector of mclust model names.")
}

if (!is.numeric(ellipse_level) ||
  length(ellipse_level) != 1L ||
  ellipse_level <= 0 ||
  ellipse_level >= 1) {
stop("`ellipse_level` must be a single number between 0 and 1.")
}

if (!is.numeric(nmm_seed) || length(nmm_seed) != 1L) {
stop("`nmm_seed` must be a single number.")
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

# ---- mclust internals --------------------------------------------------------

# `mclust::Mclust()` rewrites its own call into a call to `mclustBIC()` and
# evaluates it in the *calling* frame, so `mclustBIC` must be reachable from
# inside this function. Binding it locally makes the call resolve both when
# the package is installed and when this file is merely sourced.
mclustBIC <- mclust::mclustBIC

# ---- Species vector ----------------------------------------------------------

species_all <- stats::setNames(
as.character(analysis_data$taxon),
rownames(analysis_data)
)

# ---- Output folder ------------------------------------------------------------

output_dir <- file.path("Figs.NMM", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (verbose) {
message("Running normal mixture models")
}

# ---- Resolve blocks ----------------------------------------------------------

if (identical(nmm_blocks, "all")) {
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
requested <- intersect(nmm_blocks, names(base_cols))
if (length(requested) == 0) {
  stop(
    "No valid blocks in `nmm_blocks`. Use \"all\" or any of: ",
    paste(names(base_cols), collapse = ", ")
  )
}
run_blocks <- base_cols[requested]
}

if (verbose) {
message("  Blocks to run: ", paste(names(run_blocks), collapse = ", "))
}

# ---- Helper: drop zero-variance traits ---------------------------------------

drop_constant <- function(df) {
keep <- vapply(df, function(x) {
  v <- stats::var(as.numeric(x), na.rm = TRUE)
  !is.na(v) && v > 0
}, logical(1))

df[, keep, drop = FALSE]
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

# ---- Helper: NMM for one trait block -----------------------------------------

run_nmm_block <- function(block_name, block_columns) {

valid_cols <- intersect(block_columns, names(analysis_data))
block_data <- analysis_data[, valid_cols, drop = FALSE]

n_all <- nrow(block_data)
p_all <- ncol(block_data)

# Complete cases first, then drop traits that are constant among them
block_data <- stats::na.omit(block_data)

if (nrow(block_data) < nmm_min_n) {
  warning(
    "Skipping '", block_name, "': only ", nrow(block_data), " of ", n_all,
    " specimens are scored for all ", p_all, " traits in this block ",
    "(required = ", nmm_min_n, "). Combined blocks often have no complete ",
    "specimens; consider running the blocks separately, or dropping the ",
    "traits responsible for most of the missing data."
  )
  return(NULL)
}

block_data <- drop_constant(block_data)

if (ncol(block_data) < 2) {
  warning(
    "Skipping '", block_name, "': fewer than 2 traits vary among the ",
    nrow(block_data), " complete specimens of this block."
  )
  return(NULL)
}

if (verbose) {
  message("Running NMM for: ", block_name)
  message("  Specimens: ", nrow(block_data), " | Traits: ", ncol(block_data))
}

# ---- PCA to rank traits by PC1 contribution --------------------------------

res_pca <- stats::prcomp(block_data, scale. = TRUE)

loadings1 <- res_pca$rotation[, "PC1"]
contrib1 <- (loadings1^2 / sum(loadings1^2)) * 100

n_top <- min(nmm_top_traits, length(contrib1))
top_traits <- names(sort(contrib1, decreasing = TRUE))[seq_len(n_top)]

if (verbose) {
  message(
    "  Top ", n_top, " traits by PC1 contribution: ",
    paste(top_traits, collapse = ", ")
  )
}

# ---- NMM data --------------------------------------------------------------

df_nmm <- drop_constant(block_data[, top_traits, drop = FALSE])

if (ncol(df_nmm) < 2) {
  warning(
    "Skipping '", block_name, "': fewer than 2 variable traits among the ",
    "top PC1 contributors."
  )
  return(NULL)
}

if (nrow(df_nmm) < nmm_min_n) {
  warning(
    "Skipping '", block_name, "': not enough specimens for the retained ",
    "traits (n = ", nrow(df_nmm), ", required = ", nmm_min_n, ")."
  )
  return(NULL)
}

top_traits <- names(df_nmm)
sp <- species_all[rownames(df_nmm)]

# ---- Fit the mixture model -------------------------------------------------

set.seed(nmm_seed)

mc_nmm <- tryCatch(
  mclust::Mclust(
    df_nmm,
    G = nmm_G,
    modelNames = nmm_model_names,
    verbose = FALSE
  ),
  error = function(e) {
    warning(
      "NMM failed for block '", block_name, "': ", conditionMessage(e)
    )
    NULL
  }
)

if (is.null(mc_nmm)) return(NULL)

n_groups <- mc_nmm$G

if (verbose) {
  message("  Model selected: ", mc_nmm$modelName)
  message("  Optimal number of groups: ", n_groups)
  message("  BIC: ", round(mc_nmm$bic, 3))
  message(
    "  Group sizes: ",
    paste(table(mc_nmm$classification), collapse = " | ")
  )
}

# Ellipses need more than one component and >= 4 specimens per component
group_sizes <- table(mc_nmm$classification)
can_ellipse <- (n_groups > 1) && all(group_sizes >= 4)
ellipse_type_use <- if (can_ellipse) "norm" else "none"

if (!can_ellipse && verbose) {
  message(
    "  Note: skipping confidence ellipses (groups = ", n_groups,
    ", smallest group = ", min(group_sizes), ")."
  )
}

# ---- BIC plot --------------------------------------------------------------

.quiet_load_namespaces("factoextra")
p_bic <- tryCatch(
  factoextra::fviz_mclust_bic(mc_nmm) +
    ggplot2::labs(title = paste0("NMM BIC - ", block_name)) +
    ggplot2::theme_bw(),
  error = function(e) {
    warning(
      "  BIC plot failed for block '", block_name, "': ",
      conditionMessage(e)
    )
    NULL
  }
)

if (!is.null(p_bic)) {
  save_outputs(
    paste0("NMM_BIC_", block_name),
    p_bic,
    base_height = 6,
    base_aspect_ratio = 1.4
  )
}

# ---- Classification plot ---------------------------------------------------

p_class <- tryCatch(
  factoextra::fviz_mclust(
    mc_nmm,
    what = "classification",
    ellipse.type = ellipse_type_use,
    ellipse.level = ellipse_level
  ) +
    ggplot2::labs(
      title = paste0("NMM classification (k = ", n_groups, ") - ", block_name)
    ) +
    ggplot2::theme_bw(),
  error = function(e) {
    warning(
      "  Classification plot failed for block '", block_name, "': ",
      conditionMessage(e)
    )
    NULL
  }
)

if (!is.null(p_class)) {
  save_outputs(
    paste0("NMM_classification_", block_name),
    p_class,
    base_height = 7,
    base_aspect_ratio = 1.3
  )
}

# ---- Uncertainty plot ------------------------------------------------------

p_unc <- tryCatch(
  factoextra::fviz_mclust(
    mc_nmm,
    what = "uncertainty",
    ellipse.type = ellipse_type_use
  ) +
    ggplot2::labs(title = paste0("NMM uncertainty - ", block_name)) +
    ggplot2::theme_bw(),
  error = function(e) {
    warning(
      "  Uncertainty plot failed for block '", block_name, "': ",
      conditionMessage(e)
    )
    NULL
  }
)

if (!is.null(p_unc)) {
  save_outputs(
    paste0("NMM_uncertainty_", block_name),
    p_unc,
    base_height = 7,
    base_aspect_ratio = 1.3
  )
}

# ---- Specimen assignments --------------------------------------------------

out_df <- as.data.frame(df_nmm)
out_df <- cbind(
  id = rownames(df_nmm),
  taxon = as.character(sp[rownames(df_nmm)]),
  out_df,
  stringsAsFactors = FALSE
)

out_df$nmm_group <- as.integer(mc_nmm$classification)
out_df$uncertainty <- mc_nmm$uncertainty
out_df$block <- block_name
rownames(out_df) <- NULL

if (save_xlsx) {
  openxlsx::write.xlsx(
    out_df,
    file.path(output_dir, paste0("NMM_groups_", block_name, ".xlsx"))
  )
}

if (verbose) {
  message("  Done: ", block_name)
}

# ---- Block summary ---------------------------------------------------------

data.frame(
  block = block_name,
  n = nrow(df_nmm),
  model = mc_nmm$modelName,
  k = n_groups,
  bic = round(mc_nmm$bic, 3),
  top_traits = paste(top_traits, collapse = "|"),
  row.names = NULL,
  stringsAsFactors = FALSE
)
}

# ---- Run all blocks ----------------------------------------------------------

all_results <- Filter(
Negate(is.null),
mapply(
  run_nmm_block,
  names(run_blocks),
  run_blocks,
  SIMPLIFY = FALSE
)
)

if (length(all_results) == 0) {
warning("No trait block could be analyzed.")
return(invisible(NULL))
}

nmm_summary <- do.call(rbind, all_results)
rownames(nmm_summary) <- NULL

if (save_xlsx) {
openxlsx::write.xlsx(
  nmm_summary,
  file.path(output_dir, "NMM_summary.xlsx")
)
}

if (verbose) {
message("NMM analysis done. Outputs saved in: ", output_dir)
}

invisible(nmm_summary)
}
