#' Non-metric multidimensional scaling of morphological trait blocks
#'
#' @description
#' Performs non-metric multidimensional scaling (NMDS) on morphological trait
#' blocks derived from a morphometric matrix, and generates an ordination
#' plot annotated with stress value and quality rating for each trait block
#' analyzed. The function is designed to work directly on the objects
#' produced by [morph_matrix_setting()].
#'
#' @details
#' `morph_nmds()` expects `analysis_data` as produced by
#' [morph_matrix_setting()]: a numeric data frame with `id` and `taxon`
#' columns plus one column per morphological trait, with specimen IDs as row
#' names. The trait blocks (`"veg"`, `"flo"`, `"fru"`, or any subset thereof)
#' are not passed in separately — they are read directly from the
#' `"base_cols"` attribute that [morph_matrix_setting()] attaches to
#' `analysis_data`. If `analysis_data` has no `"base_cols"` attribute, the
#' function stops with an informative error.
#'
#' The `nmds_blocks` argument controls which trait sets are analyzed:
#'
#' \itemize{
#'   \item `"all"`: runs NMDS on each individual block in `base_cols`
#'   (typically `"veg"`, `"flo"`, `"fru"`), as well as on every pairwise and
#'   full combination of these blocks (`"vegflo"`, `"vegfru"`, `"flofru"`,
#'   `"vegflofru"`) that can be fully built from the blocks present;
#'   \item a character vector naming one or more blocks in `base_cols`: runs
#'   NMDS only on the requested block(s), individually.
#' }
#'
#' For each trait block, specimens with no observations across the block's
#' traits are removed. When `nmds_dist = "gower"`, the Gower dissimilarity
#' ([cluster::daisy()]) is used, which tolerates missing values, so
#' specimens with partial trait coverage are retained. For any other
#' distance method, rows with any missing values are additionally removed
#' via [stats::na.omit()] before computing the distance matrix with
#' [vegan::vegdist()]. Blocks with fewer than 3 remaining specimens are
#' skipped with a warning.
#'
#' NMDS is fit with [vegan::metaMDS()] using the requested number of
#' dimensions (`nmds_k`) and maximum number of random starts (`nmds_trymax`).
#' A random seed (`nmds_seed`) is set immediately before each block's NMDS
#' run so that results are reproducible across calls. If `metaMDS()` fails
#' to converge or otherwise errors for a block, that block is skipped with a
#' warning.
#'
#' Stress values are classified for interpretation as `"excellent"` (< 0.05),
#' `"good"` (< 0.10), `"acceptable"` (< 0.20), or `"poor"` (>= 0.20).
#'
#' Taxon colours are taken from the `"taxon_colors"` attribute attached to
#' `analysis_data` by [morph_matrix_setting()], so that each taxon keeps the
#' same colour used in other plots (e.g. [morph_boxplots()], [morph_pca()]).
#' If this attribute is missing, or lacks an entry for a taxon present in a
#' given block, colours are generated automatically from the `viridis`
#' palette for the affected taxa, with a warning.
#'
#' For each analyzed block, an NMDS ordination scatter plot (NMDS1 vs NMDS2)
#' annotated with the stress value and its quality rating is generated.
#' Output files are written to a date-specific directory inside `Figs_NMDS`,
#' in the formats specified by `file_formats`.
#'
#' @param analysis_data A numeric data frame containing morphological traits
#' as columns and specimens as rows, with `id` and `taxon` columns and
#' specimen identifiers as row names. Typically the output of
#' [morph_matrix_setting()]. If this object carries a `"taxon_colors"`
#' attribute, those colours are reused; otherwise a `viridis` palette is
#' generated automatically.
#'
#' @param nmds_blocks Character. Either `"all"`, to run NMDS on every
#' individual block in `base_cols` plus all pairwise and full combinations of
#' these blocks, or a character vector naming one or more blocks in
#' `base_cols` to analyze individually. Default is `"all"`.
#'
#' @param nmds_dist Character string giving the distance/dissimilarity method.
#' `"gower"` (default) uses [cluster::daisy()] and tolerates missing values;
#' any other value (e.g. `"bray"`, `"jaccard"`, `"euclidean"`, `"manhattan"`)
#' is passed to [vegan::vegdist()], which requires complete cases.
#'
#' @param nmds_k Integer. Number of NMDS dimensions to fit. Default is `2`.
#'
#' @param nmds_trymax Integer. Maximum number of random starts used by
#' [vegan::metaMDS()] to find a stable solution. Default is `100`.
#'
#' @param nmds_seed Integer. Random seed set before each block's NMDS run, for
#' reproducibility. Default is `42`.
#'
#' @param file_formats Character vector specifying which image formats to
#' save. One or more of `"pdf"` and `"jpeg"`. Default is
#' `c("pdf", "jpeg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' including skipped blocks, per-block stress values, and the output
#' directory, to the console.
#'
#' @return
#' Invisibly returns a data frame summarizing, for each successfully analyzed
#' block, the block name, number of specimens retained, distance method,
#' stress value, and stress quality rating. Returns an empty data frame
#' (invisibly) if no block could be analyzed.
#'
#' As a side effect, one NMDS ordination plot per analyzed block is written
#' to disk.
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
#' # Run NMDS on every individual block and every block combination
#' morph_nmds(analysis_data = res, nmds_blocks = "all")
#'
#' # Run NMDS on the vegetative block only, using Bray-Curtis distance
#' morph_nmds(
#'   analysis_data = res,
#'   nmds_blocks = "veg",
#'   nmds_dist = "bray"
#' )
#'
#' # Run NMDS on selected blocks, saving PDF output only
#' summary_nmds <- morph_nmds(
#'   analysis_data = res,
#'   nmds_blocks = c("veg", "flo"),
#'   file_formats = "pdf"
#' )
#' summary_nmds
#' }
#'
#' @importFrom stats setNames na.omit
#' @importFrom ggplot2 aes annotate element_blank element_text geom_hline
#'   geom_point geom_vline ggplot labs scale_fill_manual theme theme_bw
#' @importFrom cluster daisy
#' @importFrom vegan vegdist metaMDS
#' @importFrom cowplot save_plot
#' @importFrom viridis viridis
#'
#' @export

morph_nmds <- function(analysis_data,
                       nmds_blocks = "all",
                       nmds_dist = "gower",
                       nmds_k = 2,
                       nmds_trymax = 100,
                       nmds_seed = 42,
                       file_formats = c("pdf", "jpeg"),
                       verbose = TRUE) {

# ---- Validate input ----------------------------------------------------------

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

if (!is.character(nmds_blocks) || length(nmds_blocks) < 1) {
stop("`nmds_blocks` must be \"all\" or a character vector of block names.")
}

if (!is.character(nmds_dist) || length(nmds_dist) != 1L) {
stop("`nmds_dist` must be a single character string.")
}

if (!is.numeric(nmds_k) || length(nmds_k) != 1L || nmds_k < 1) {
stop("`nmds_k` must be a positive integer.")
}

if (!is.numeric(nmds_trymax) || length(nmds_trymax) != 1L || nmds_trymax < 1) {
stop("`nmds_trymax` must be a positive integer.")
}

if (!is.numeric(nmds_seed) || length(nmds_seed) != 1L) {
stop("`nmds_seed` must be a single number.")
}

file_formats <- match.arg(
file_formats,
choices = c("pdf", "jpeg"),
several.ok = TRUE
)

if (!is.logical(verbose) || length(verbose) != 1L) {
stop("`verbose` must be TRUE or FALSE.")
}

nmds_min_n <- 3

# ---- Species vector -----------------------------------------------------------

species_all <- stats::setNames(
as.character(analysis_data$taxon),
rownames(analysis_data)
)

# ---- Output folder -------------------------------------------------------------

output_dir <- file.path("Figs_NMDS", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (verbose) {
message("Running NMDS...")
}

# ---- Resolve which blocks to run -----------------------------------------------

if (identical(nmds_blocks, "all")) {
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
requested <- intersect(nmds_blocks, names(base_cols))
if (length(requested) == 0) {
  stop(
    "No valid blocks in `nmds_blocks`. Use \"all\" or any of: ",
    paste(names(base_cols), collapse = ", ")
  )
}
run_blocks <- base_cols[requested]
}

if (verbose) {
message("  Blocks to run: ", paste(names(run_blocks), collapse = ", "))
}

col_names <- names(analysis_data)

# ---- Helper: run NMDS for one block --------------------------------------------

run_nmds_block <- function(block_name, block_columns) {

if (verbose) {
  message("\n  --- Block: ", block_name, " ---")
}

valid_cols <- intersect(block_columns, col_names)
block_data <- analysis_data[, valid_cols, drop = FALSE]
sp         <- species_all[rownames(block_data)]

# Remove all-NA rows
block_data <- block_data[rowSums(!is.na(block_data)) > 0, , drop = FALSE]
sp <- sp[rownames(block_data)]

if (nrow(block_data) < nmds_min_n) {
  warning(
    "Skipping '", block_name, "': not enough specimens (n = ",
    nrow(block_data), ")."
  )
  return(invisible(NULL))
}

if (verbose) {
  message("  Specimens: ", nrow(block_data), " | Traits: ", ncol(block_data))
}

# ---- Distance matrix --------------------------------------------------------

if (nmds_dist == "gower") {
  dist_mat <- stats::as.dist(cluster::daisy(block_data, metric = "gower"))
} else {
  block_data <- stats::na.omit(block_data)
  sp <- sp[rownames(block_data)]

  if (nrow(block_data) < nmds_min_n) {
    warning(
      "Skipping '", block_name, "': not enough complete specimens for '",
      nmds_dist, "'. Use `nmds_dist = \"gower\"` for incomplete data."
    )
    return(invisible(NULL))
  }

  dist_mat <- vegan::vegdist(block_data, method = nmds_dist)
}

# ---- NMDS ---------------------------------------------------------------

set.seed(nmds_seed)
nmds_result <- tryCatch(
  vegan::metaMDS(dist_mat, k = nmds_k, trymax = nmds_trymax, trace = FALSE),
  error = function(e) {
    warning("NMDS failed for block '", block_name, "': ", conditionMessage(e))
    return(NULL)
  }
)

if (is.null(nmds_result)) return(invisible(NULL))

# ---- Results --------------------------------------------------------------

stress <- round(nmds_result$stress, 4)
stress_label <- if (stress < 0.05) "excellent" else
  if (stress < 0.10) "good" else
    if (stress < 0.20) "acceptable" else "poor"

if (verbose) {
  cat("\n  Block     :", block_name, "\n")
  cat("  Specimens :", nrow(block_data), "\n")
  cat("  Distance  :", nmds_dist, "\n")
  cat("  Stress    :", stress, "(", stress_label, ")\n\n")
}

# ---- Plot data --------------------------------------------------------------

nmds_df <- as.data.frame(nmds_result$points)
colnames(nmds_df) <- paste0("NMDS", seq_len(ncol(nmds_df)))
nmds_df$taxon <- as.character(sp[rownames(nmds_df)])

# ---- Taxon colours -------------------------------------------------------

color_vector <- attr(analysis_data, "taxon_colors")

if (is.null(color_vector)) {

  color_vector <- stats::setNames(
    viridis::viridis(length(unique(nmds_df$taxon))),
    sort(unique(nmds_df$taxon))
  )

} else {

  missing_sp <- setdiff(unique(nmds_df$taxon), names(color_vector))

  if (length(missing_sp) > 0) {
    warning(
      "No stored colour found for: ", paste(missing_sp, collapse = ", "),
      ". Assigning automatic colours for these taxa."
    )
    color_vector[missing_sp] <- viridis::viridis(length(missing_sp))
  }

  color_vector <- color_vector[sort(unique(nmds_df$taxon))]
}

# ---- Plot ---------------------------------------------------------------

p_nmds <- ggplot2::ggplot(nmds_df, ggplot2::aes(x = .data[["NMDS1"]], y = .data[["NMDS2"]])) +
  ggplot2::geom_point(ggplot2::aes(fill = .data[["taxon"]]), shape = 21, size = 3, colour = "white") +
  ggplot2::scale_fill_manual(values = color_vector) +
  ggplot2::annotate(
    "text",
    x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
    label = paste0("Stress = ", stress, " (", stress_label, ")"),
    size = 3.5
  ) +
  ggplot2::geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::labs(
    x = "NMDS1", y = "NMDS2", fill = "",
    title = paste0("NMDS \u2014 ", block_name)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(face = "italic")
  )

if ("pdf" %in% file_formats) {
  cowplot::save_plot(
    file.path(output_dir, paste0("NMDS_plot_", block_name, ".pdf")),
    p_nmds, ncol = 1, nrow = 1, base_height = 7, base_aspect_ratio = 1.3
  )
}

if ("jpeg" %in% file_formats) {
  cowplot::save_plot(
    file.path(output_dir, paste0("NMDS_plot_", block_name, ".jpeg")),
    p_nmds, ncol = 1, nrow = 1, base_height = 7, base_aspect_ratio = 1.3
  )
}

# ---- Return summary -----------------------------------------------------

data.frame(block = block_name,
           n = nrow(block_data),
           distance = nmds_dist,
           stress   = stress,
           quality  = stress_label)
}

# ---- Run all blocks -------------------------------------------------------------

all_results <- Filter(
Negate(is.null),
mapply(
  run_nmds_block,
  names(run_blocks), run_blocks,
  SIMPLIFY = FALSE
)
)

summary_df <- if (length(all_results) > 0) {
do.call(rbind, all_results)
} else {
data.frame(
  block = character(0), n = integer(0), distance = character(0),
  stress = numeric(0), quality = character(0)
)
}

if (verbose) {
message("\nDone. Outputs saved in: ", output_dir)
}

invisible(summary_df)
}
