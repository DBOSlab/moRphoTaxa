#' Run PCA on morphological trait blocks
#'
#' @description
#' Performs principal component analysis (PCA) on morphological trait blocks
#' derived from a morphometric matrix, and generates a standard set of
#' diagnostic figures for each trait block analyzed. The function is designed
#' to work directly on the objects produced by [morph_matrix_setting()].
#'
#' For each requested trait block, the function fits a PCA, saves a biplot
#' with marginal density panels, saves variable-contribution plots for the
#' first two principal components, and — when at least three components are
#' available — saves a panel of pairwise PC scatterplots (PC1 x PC2, PC1 x
#' PC3, PC2 x PC3).
#'
#' @details
#' `morph_pca()` expects `analysis_data` as produced by
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
#' The `pca_blocks` argument controls which trait sets are analyzed:
#'
#' \itemize{
#'   \item `"all"`: runs PCA on each individual block in `base_cols`
#'   (typically `"veg"`, `"flo"`, `"fru"`), as well as on every pairwise and
#'   full combination of these blocks (`"vegflo"`, `"vegfru"`, `"flofru"`,
#'   `"vegflofru"`);
#'   \item a character vector naming one or more blocks in `base_cols`: runs
#'   PCA only on the requested block(s), individually.
#' }
#'
#' For each trait block, rows with missing values in the relevant columns are
#' removed via [stats::na.omit()] prior to fitting the PCA. Blocks with fewer
#' than 3 specimens or fewer than 2 traits after this filtering step are
#' skipped with a warning.
#'
#' PCA is fit with [stats::prcomp()] using `scale. = TRUE` (i.e., on the
#' correlation matrix). Taxon colours are taken from the `"taxon_colors"`
#' attribute attached to `analysis_data` by [morph_matrix_setting()], so that
#' each taxon keeps the same colour used in other plots (e.g.
#' [morph_boxplots()]). If this attribute is missing, or lacks an entry for a
#' taxon present in a given block, colours are generated automatically from
#' the `viridis` palette for the affected taxa, with a warning.
#'
#' For each analyzed block, the following figures are generated:
#'
#' \itemize{
#'   \item a PCA biplot (individuals and variables) with marginal density
#'   panels along each axis, one density per taxon;
#'   \item a variable-contribution barplot for principal component 1;
#'   \item a variable-contribution barplot for principal component 2;
#'   \item when at least 3 components are available, a three-panel figure
#'   with PC1 x PC2, PC1 x PC3, and PC2 x PC3 scatterplots.
#' }
#'
#' Output files are written to a date-specific directory inside `Figs.PCA`.
#' Each figure is saved in the formats specified by `file_formats`. File
#' names are prefixed by figure type (e.g. `"PCA_"`, `"Dim1_contrib_variables_"`,
#' `"Dim2_contrib_variables_"`, `"PCA_panels_"`) and suffixed by the trait
#' block name.
#'
#' @param analysis_data A numeric data frame containing morphological traits
#' as columns and specimens as rows, with `id` and `taxon` columns and
#' specimen identifiers as row names. Typically the `analysis_data` element
#' returned by [morph_matrix_setting()]. If this object carries a
#' `"taxon_colors"` attribute, those colours are reused; otherwise a
#' `viridis` palette is generated automatically.
#'
#' @param pca_blocks Character. Either `"all"`, to run PCA on every individual
#' block in `base_cols` plus all pairwise and full combinations of these
#' blocks, or a character vector naming one or more blocks in `base_cols` to
#' analyze individually. Default is `"all"`.
#'
#' @param top_contrib Integer. Number of top-contributing variables to display
#' in the variable-contribution plots for principal components 1 and 2.
#' Default is `10`.
#'
#' @param file_formats Character vector specifying which image formats to
#' save. One or more of `"pdf"` and `"jpeg"`. Default is
#' `c("pdf", "jpeg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' including skipped blocks and the output directory, to the console.
#'
#' @return
#' Invisibly returns `NULL`.
#'
#' The function is primarily called for its side effects: PCA biplots,
#' variable-contribution plots, and (when applicable) pairwise PC scatterplot
#' panels are written to disk for each analyzed trait block.
#'
#' @seealso
#' [morph_matrix_setting()],
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
#' # Run PCA on every individual block and every block combination
#' morph_pca(
#'   analysis_data = res,
#'   pca_blocks = "all"
#' )
#'
#' # Run PCA on the vegetative block only
#' morph_pca(
#'   analysis_data = res,
#'   pca_blocks = "veg"
#' )
#'
#' # Run PCA on selected blocks, showing the top 15 contributing variables,
#' # saving PDF output only
#' morph_pca(
#'   analysis_data = res,
#'   pca_blocks = c("veg", "flo"),
#'   top_contrib = 15,
#'   file_formats = "pdf"
#' )
#' }
#'
#' @importFrom stats prcomp na.omit setNames
#' @importFrom ggplot2 aes element_blank element_text geom_hline geom_point
#'   geom_vline ggplot ggtitle guide_colorbar guide_legend guides labs
#'   scale_fill_manual theme theme_bw
#' @importFrom cowplot save_plot plot_grid
#' @importFrom viridis viridis
#' @importFrom ggside geom_xsidedensity geom_ysidedensity theme_ggside_void
#'
#' @export

morph_pca <- function(analysis_data,
                      pca_blocks = "all",
                      top_contrib = 10,
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

if (!is.character(pca_blocks) || length(pca_blocks) < 1) {
  stop("`pca_blocks` must be \"all\" or a character vector of block names.")
}

if (!is.numeric(top_contrib) ||
    length(top_contrib) != 1L ||
    top_contrib < 1) {
  stop("`top_contrib` must be a positive integer.")
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

output_dir <- file.path("Figs.PCA", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (verbose) {
  message("Running PCA analysis")
}

# ---- Resolve blocks ----------------------------------------------------------

if (identical(pca_blocks, "all")) {
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
  requested <- intersect(pca_blocks, names(base_cols))
  if (length(requested) == 0) {
    stop(
      "No valid blocks in `pca_blocks`. Use \"all\" or any of: ",
      paste(names(base_cols), collapse = ", ")
    )
  }
  run_blocks <- base_cols[requested]
}

if (verbose) {
  message("Blocks to run:", paste(names(run_blocks), collapse = ", "))
}

# ---- Helper: PC scatter plot -------------------------------------------------

make_pc_plot <- function(df, x, y, color_vector) {
  ggplot2::ggplot(df,
    ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[["taxon"]])
  ) +
    ggplot2::geom_point(shape = 21, size = 3, colour = "white") +
    ggplot2::scale_fill_manual(values = color_vector) +
    ggplot2::labs(x = x, y = y, fill = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(face = "italic")
    ) +
    ggplot2::geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, color = "gray40", linetype = "dashed")
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

# ---- PCA + plots for each trait set ------------------------------------------

for (set_name in names(run_blocks)) {

valid_cols <- intersect(run_blocks[[set_name]], names(analysis_data))
d <- stats::na.omit(analysis_data[, valid_cols, drop = FALSE])

if (nrow(d) < 3 || ncol(d) < 2) {
  warning(
    "Skipping '", set_name, "': not enough data (n = ", nrow(d),
    ", p = ", ncol(d), ")."
  )
  next
}

sp <- species_all[rownames(d)]

if (verbose) {
  message("Running PCA for: ", set_name)
}

res.pca <- stats::prcomp(d, scale. = TRUE)

# ---- Taxon colours -------------------------------------------------------
color_vector <- attr(analysis_data, "taxon_colors")

if (is.null(color_vector)) {

  color_vector <- stats::setNames(
    viridis::viridis(length(unique(sp))),
    sort(unique(sp))
  )

} else {

  missing_sp <- setdiff(unique(sp), names(color_vector))

  if (length(missing_sp) > 0) {
    warning(
      "No stored colour found for: ", paste(missing_sp, collapse = ", "),
      ". Assigning automatic colours for these taxa."
    )
    color_vector[missing_sp] <- viridis::viridis(length(missing_sp))
  }

  color_vector <- color_vector[sort(unique(sp))]
}

# ---- Biplot with marginal densities -----------------------------------

.quiet_load_namespaces("factoextra")
p <- factoextra::fviz_pca_biplot(
  res.pca, geom = c("point", "text"),
  label = "var",
  alpha.var = "contrib", col.var = "contrib",
  fill.ind = sp, col.ind = "white",
  pointshape = 21, pointsize = 3, pointwidth = 0.5,
  gradient.cols = "RdBu"
) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(order = 0),
    fill = ggplot2::guide_legend(order = 1)
  ) +
  ggplot2::scale_fill_manual(values = color_vector) +
  ggplot2::labs(fill = "") +
  ggplot2::ggtitle("") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(face = "italic")
  ) +
  ggplot2::geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0, color = "gray40", linetype = "dashed") +
  ggside::geom_xsidedensity(ggplot2::aes(fill = sp), alpha = 0.5, lwd = 0.1, show.legend = FALSE) +
  ggside::geom_ysidedensity(ggplot2::aes(fill = sp), alpha = 0.5, lwd = 0.1, show.legend = FALSE) +
  ggside::theme_ggside_void()

save_outputs(paste0("PCA_", set_name), p, base_height = 8.5, base_aspect_ratio = 1.3)

# ---- Variable-contribution plots ---------------------------------------

p_contrib1 <- factoextra::fviz_contrib(res.pca, choice = "var", axes = 1, top = top_contrib)
p_contrib2 <- factoextra::fviz_contrib(res.pca, choice = "var", axes = 2, top = top_contrib)

save_outputs(paste0("Dim1_contrib_variables_", set_name), p_contrib1, base_height = 8.5, base_aspect_ratio = 1.3)
save_outputs(paste0("Dim2_contrib_variables_", set_name), p_contrib2, base_height = 8.5, base_aspect_ratio = 1.3)

# ---- PC1xPC2, PC1xPC3, PC2xPC3 scatter panels ---------------------------

if (ncol(res.pca$x) >= 3) {

  pca_df <- as.data.frame(res.pca$x)
  pca_df$taxon <- sp[rownames(pca_df)]

  p_all_pc <- cowplot::plot_grid(
    make_pc_plot(pca_df, "PC1", "PC2", color_vector),
    make_pc_plot(pca_df, "PC1", "PC3", color_vector),
    make_pc_plot(pca_df, "PC2", "PC3", color_vector),
    labels = c("A", "B", "C"), ncol = 1, nrow = 3, align = "hv"
  )

  save_outputs(paste0("PCA_panels_", set_name), p_all_pc, base_height = 7, base_aspect_ratio = 1.3, nrow = 3)
}

if (verbose) {
  message("  Done: ", set_name)
}
}

if (verbose) {
  message("PCA analysis done. Outputs saved in: ", output_dir)
}

invisible(NULL)
}
