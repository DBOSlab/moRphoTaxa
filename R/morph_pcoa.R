#' Run PCoA on morphological trait blocks
#'
#' @description
#' Performs principal coordinates analysis (PCoA) on morphological trait
#' blocks derived from a morphometric matrix, and generates a standard set of
#' diagnostic figures for each trait block analyzed. The function is designed
#' to work directly on the objects produced by [morph_matrix_setting()], and
#' mirrors the structure of [morph_pca()].
#'
#' For each requested trait block, the function computes a specimen-by-specimen
#' distance matrix, fits a PCoA, saves an ordination biplot with marginal
#' density panels and vectors for the traits most correlated with the first
#' two axes, saves a barplot of the correlation of the top traits with each of
#' these axes, saves a barplot of the variation explained by the first axes,
#' and, when at least three axes are available, saves a panel of pairwise
#' scatterplots (PCo1 x PCo2, PCo1 x PCo3, PCo2 x PCo3).
#'
#' @details
#' `morph_pcoa()` expects `analysis_data` as produced by
#' [morph_matrix_setting()]: a numeric data frame with `id` and `taxon`
#' columns plus one column per morphological trait, with specimen IDs as row
#' names. The trait blocks (`"veg"`, `"flo"`, `"fru"`, or any subset thereof)
#' are not passed in separately: they are read directly from the
#' `"base_cols"` attribute that [morph_matrix_setting()] attaches to
#' `analysis_data`. If `analysis_data` has no `"base_cols"` attribute (e.g.
#' it was not built with [morph_matrix_setting()], or the attribute was
#' dropped by an intervening subsetting operation), the function stops with
#' an informative error.
#'
#' The `pcoa_blocks` argument controls which trait sets are analyzed:
#'
#' \itemize{
#'   \item `"all"`: runs PCoA on each individual block in `base_cols`
#'   (typically `"veg"`, `"flo"`, `"fru"`), as well as on every pairwise and
#'   full combination of these blocks (`"vegflo"`, `"vegfru"`, `"flofru"`,
#'   `"vegflofru"`);
#'   \item a character vector naming one or more blocks in `base_cols`: runs
#'   PCoA only on the requested block(s), individually.
#' }
#'
#' For each trait block, rows with missing values in the relevant columns are
#' removed via [stats::na.omit()], and traits with zero variance are dropped.
#' Blocks with fewer than 3 specimens or fewer than 2 traits after this
#' filtering, or with fewer than 2 ordination axes, are skipped with a
#' warning.
#'
#' \strong{Distance.} The distance between specimens is set by `dist_method`:
#'
#' \itemize{
#'   \item `"gower"` (default): Gower distance for continuous traits, i.e. the
#'   mean, over traits, of the absolute difference between two specimens
#'   divided by the range of the trait. All traits have the same weight
#'   regardless of their units;
#'   \item `"euclidean"`: Euclidean distance on centred and scaled traits. The
#'   resulting PCoA is equivalent to the PCA computed by [morph_pca()].
#' }
#'
#' \strong{Ordination.} The PCoA is computed with [ape::pcoa()]. Distances
#' such as Gower's are not always Euclidean, in which case some eigenvalues are
#' negative. `correction` controls how this is handled: `"none"` (default)
#' uses only the axes with positive eigenvalues; `"lingoes"` and `"cailliez"`
#' apply the corresponding additive-constant corrections. The percentage of
#' variation shown on each axis is its eigenvalue divided by the sum of the
#' positive eigenvalues (corrected eigenvalues when a correction is used). When
#' negative eigenvalues are present and `correction = "none"`, a message
#' reports their magnitude.
#'
#' \strong{Traits.} PCoA is computed from distances, so it has no trait
#' loadings of its own. As in the usual PCoA biplot, each trait is instead
#' correlated ([stats::cor()]) with the ordination axes: the vectors in the
#' biplot are drawn for the `top_vars` traits with the highest squared
#' correlation with PCo1 and PCo2 combined, with lengths proportional to the
#' correlations; and the two barplots show the correlation of the `top_vars`
#' traits most correlated with PCo1 and with PCo2, respectively.
#'
#' Taxon colours are taken from the `"taxon_colors"` attribute attached to
#' `analysis_data` by [morph_matrix_setting()], so that each taxon keeps the
#' same colour used in other plots (e.g. [morph_pca()], [morph_dapc()]). If
#' this attribute is missing, or lacks an entry for a taxon present in a given
#' block, colours are generated automatically from the `viridis` palette for
#' the affected taxa, with a warning.
#'
#' For each analyzed block, the following figures are generated:
#'
#' \itemize{
#'   \item a PCoA biplot (specimens and trait vectors) with marginal density
#'   panels along each axis, one density per taxon;
#'   \item a barplot of the correlation of the top traits with PCo1;
#'   \item a barplot of the correlation of the top traits with PCo2;
#'   \item a barplot of the percentage of variation explained by the first
#'   (up to 10) axes;
#'   \item when at least 3 axes are available, a three-panel figure with
#'   PCo1 x PCo2, PCo1 x PCo3, and PCo2 x PCo3 scatterplots.
#' }
#'
#' Output files are written to a date-specific directory inside `Figs.PCoA`.
#' Each figure is saved in the formats specified by `file_formats`. File
#' names are prefixed by figure type (`"PCoA_"`, `"Dim1_corr_variables_"`,
#' `"Dim2_corr_variables_"`, `"PCoA_eigenvalues_"`, `"PCoA_panels_"`) and
#' suffixed by the trait block name.
#'
#' @param analysis_data A numeric data frame containing morphological traits
#' as columns and specimens as rows, with `id` and `taxon` columns and
#' specimen identifiers as row names. Typically the `analysis_data` element
#' returned by [morph_matrix_setting()]. If this object carries a
#' `"taxon_colors"` attribute, those colours are reused; otherwise a
#' `viridis` palette is generated automatically.
#'
#' @param pcoa_blocks Character. Either `"all"`, to run PCoA on every
#' individual block in `base_cols` plus all pairwise and full combinations of
#' these blocks, or a character vector naming one or more blocks in
#' `base_cols` to analyze individually. Default is `"all"`.
#'
#' @param dist_method Character. Distance between specimens: `"gower"`
#' (default) or `"euclidean"` (on scaled traits). See Details.
#'
#' @param correction Character. Correction for negative eigenvalues: `"none"`
#' (default), `"lingoes"` or `"cailliez"`. See Details.
#'
#' @param top_vars Integer. Number of traits shown as vectors in the biplot
#' and as bars in the trait-correlation plots for PCo1 and PCo2. Default is
#' `10`.
#'
#' @param file_formats Character vector specifying which image formats to
#' save. One or more of `"pdf"` and `"jpeg"`. Default is
#' `c("pdf", "jpeg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' including skipped blocks, negative eigenvalues and the output directory, to
#' the console.
#'
#' @return
#' Invisibly returns a named list with one element per analyzed block. Each
#' element is a list with:
#' \describe{
#'   \item{`scores`}{A data frame with the ordination coordinates of each
#'   specimen (`PCo1`, `PCo2`, ...) and its `taxon`, with specimen IDs as row
#'   names.}
#'   \item{`eigenvalues`}{A data frame with the `axis`, `eigenvalue` and
#'   `relative` (proportion of the sum of positive eigenvalues) of each axis.}
#'   \item{`n_specimens`, `n_traits`}{Number of specimens and traits used.}
#' }
#' The list is empty if no block could be analyzed.
#'
#' The function is primarily called for its side effects: PCoA biplots,
#' trait-correlation plots, eigenvalue plots and (when applicable) pairwise
#' axis scatterplot panels are written to disk for each analyzed trait block.
#'
#' @seealso
#' [morph_matrix_setting()],
#' [morph_pca()],
#' [morph_dapc()]
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
#' # Run PCoA (Gower distance) on every individual block and every combination
#' morph_pcoa(
#'   analysis_data = res,
#'   pcoa_blocks = "all"
#' )
#'
#' # Vegetative block only, with a Lingoes correction for negative eigenvalues
#' morph_pcoa(
#'   analysis_data = res,
#'   pcoa_blocks = "veg",
#'   correction = "lingoes"
#' )
#'
#' # Selected blocks, Euclidean distance on scaled traits, top 15 traits,
#' # saving PDF output only
#' pco <- morph_pcoa(
#'   analysis_data = res,
#'   pcoa_blocks = c("veg", "flo"),
#'   dist_method = "euclidean",
#'   top_vars = 15,
#'   file_formats = "pdf"
#' )
#'
#' # Ordination coordinates of the vegetative block
#' head(pco$veg$scores)
#' }
#'
#' @importFrom stats na.omit setNames var cor dist
#' @importFrom ape pcoa
#' @importFrom ggplot2 aes coord_flip element_blank element_text geom_col
#'   geom_hline geom_point geom_segment geom_text geom_vline ggplot labs
#'   scale_fill_manual theme theme_bw
#' @importFrom cowplot save_plot plot_grid
#' @importFrom viridis viridis
#' @importFrom ggside geom_xsidedensity geom_ysidedensity theme_ggside_void
#' @importFrom grid arrow unit
#' @importFrom rlang .data
#'
#' @export

morph_pcoa <- function(analysis_data,
                       pcoa_blocks = "all",
                       dist_method = c("gower", "euclidean"),
                       correction = c("none", "lingoes", "cailliez"),
                       top_vars = 10,
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

if (!is.character(pcoa_blocks) || length(pcoa_blocks) < 1) {
  stop("`pcoa_blocks` must be \"all\" or a character vector of block names.")
}

dist_method <- match.arg(dist_method)
correction  <- match.arg(correction)

if (!is.numeric(top_vars) ||
    length(top_vars) != 1L ||
    top_vars < 1) {
  stop("`top_vars` must be a positive integer.")
}

file_formats <- match.arg(
  file_formats,
  choices = c("pdf", "jpeg"),
  several.ok = TRUE
)

if (!is.logical(verbose) || length(verbose) != 1L) {
  stop("`verbose` must be TRUE or FALSE.")
}

if (verbose) {
  message("Running PCoA analysis (", dist_method, " distance)")
}

# ---- Species vector ----------------------------------------------------------

species_all <- stats::setNames(
  as.character(analysis_data$taxon),
  rownames(analysis_data)
)

# ---- Output folder ------------------------------------------------------------

output_dir <- file.path("Figs.PCoA", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Resolve blocks ----------------------------------------------------------

if (identical(pcoa_blocks, "all")) {
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
  requested <- intersect(pcoa_blocks, names(base_cols))
  if (length(requested) == 0) {
    stop(
      "No valid blocks in `pcoa_blocks`. Use \"all\" or any of: ",
      paste(names(base_cols), collapse = ", ")
    )
  }
  run_blocks <- base_cols[requested]
}

if (verbose) {
  message("Blocks to run: ", paste(names(run_blocks), collapse = ", "))
}

# ---- Helper: PCo scatter plot -------------------------------------------------

make_pco_plot <- function(df, x, y, color_vector, lab_x, lab_y) {
  ggplot2::ggplot(df,
                  ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[["taxon"]])
  ) +
    ggplot2::geom_point(shape = 21, size = 3, colour = "white") +
    ggplot2::scale_fill_manual(values = color_vector) +
    ggplot2::labs(x = lab_x, y = lab_y, fill = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(face = "italic")
    ) +
    ggplot2::geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, color = "gray40", linetype = "dashed")
}

# ---- Helper: trait-correlation barplot ----------------------------------------

make_corr_plot <- function(r, axis_label, top_n) {
ord <- order(abs(r), decreasing = TRUE)[seq_len(min(top_n, length(r)))]
df <- data.frame(variable = names(r)[ord], r = unname(r[ord]))
df$variable <- factor(df$variable, levels = rev(df$variable))
df$sign <- ifelse(df$r >= 0, "positive", "negative")

ggplot2::ggplot(df, ggplot2::aes(x = .data[["variable"]], y = .data[["r"]],
                                 fill = .data[["sign"]])) +
  ggplot2::geom_col(show.legend = FALSE) +
  ggplot2::scale_fill_manual(values = c(positive = "#2166ac", negative = "#b2182b")) +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = paste0("Correlation with ", axis_label)) +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
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

# ---- PCoA + plots for each trait set -------------------------------------------

results <- list()

for (set_name in names(run_blocks)) {

  valid_cols <- intersect(run_blocks[[set_name]], names(analysis_data))
  d <- stats::na.omit(analysis_data[, valid_cols, drop = FALSE])

  if (nrow(d) < 3) {
    warning(
      "Skipping '", set_name, "': fewer than 3 specimens with complete data ",
      "(n = ", nrow(d), ")."
    )
    next
  }

  # Drop traits with no variation (undefined for range/scale standardisation)
  keep_trait <- vapply(d, function(x) {
    v <- stats::var(x)
    !is.na(v) && v > 0
  }, logical(1))
  d <- d[, keep_trait, drop = FALSE]

  if (ncol(d) < 2) {  warning(
    "Skipping '", set_name, "': not enough data (n = ", nrow(d),
    ", p = ", ncol(d), ")."
  )
  next
}

sp <- species_all[rownames(d)]

if (verbose) {
  message("Running PCoA for: ", set_name)
}

# ---- Distance matrix -------------------------------------------------------
d_mat <- as.matrix(d)

if (dist_method == "gower") {
  rng   <- apply(d_mat, 2, function(x) diff(range(x)))
  z     <- sweep(sweep(d_mat, 2, apply(d_mat, 2, min), "-"), 2, rng, "/")
  dmat  <- stats::dist(z, method = "manhattan") / ncol(z)
} else {
  dmat  <- stats::dist(scale(d_mat))
}

# ---- PCoA -------------------------------------------------------------------
pc <- tryCatch(
  ape::pcoa(dmat, correction = correction),
  error = function(e) {
    warning("PCoA failed for '", set_name, "': ", conditionMessage(e))
    NULL
  }
)

if (is.null(pc)) next

if (correction == "none") {
  scores <- pc$vectors
  eig    <- pc$values$Eigenvalues
} else {
  scores <- pc$vectors.cor
  eig    <- pc$values$Corr_eig
}

if (is.null(scores) || ncol(scores) < 2) {
  warning("Skipping '", set_name, "': fewer than 2 ordination axes.")
  next
}

rel_eig <- eig / sum(eig[eig > 0])

if (verbose && correction == "none" && any(eig < -1e-8)) {
  message(sprintf(
    "  %d negative eigenvalue(s), summing to %.1f%% of the positive sum (consider `correction`).",
    sum(eig < -1e-8), 100 * abs(sum(eig[eig < -1e-8])) / sum(eig[eig > 0])
  ))
}

colnames(scores) <- paste0("PCo", seq_len(ncol(scores)))
rownames(scores) <- rownames(d)

pco_df <- as.data.frame(scores)
pco_df$taxon <- sp[rownames(pco_df)]

lab_x <- sprintf("PCo1 (%.1f%%)", 100 * rel_eig[1])
lab_y <- sprintf("PCo2 (%.1f%%)", 100 * rel_eig[2])

# ---- Taxon colours -----------------------------------------------------------
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

# ---- Trait correlations with PCo1 and PCo2 -----------------------------------
trait_r <- suppressWarnings(stats::cor(d_mat, scores[, 1:2, drop = FALSE]))

# ---- Biplot with marginal densities ------------------------------------------
top_n   <- min(top_vars, nrow(trait_r))
top_idx <- order(rowSums(trait_r^2), decreasing = TRUE)[seq_len(top_n)]
mult    <- 0.85 * max(abs(scores[, 1:2]))

arrows_df <- data.frame(
  variable = rownames(trait_r)[top_idx],
  xend = trait_r[top_idx, 1] * mult,
  yend = trait_r[top_idx, 2] * mult
)

p <- ggplot2::ggplot(
  pco_df,
  ggplot2::aes(x = .data[["PCo1"]], y = .data[["PCo2"]], fill = .data[["taxon"]])
) +
  ggplot2::geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::geom_segment(
    data = arrows_df,
    ggplot2::aes(x = 0, y = 0, xend = .data[["xend"]], yend = .data[["yend"]]),
    arrow = grid::arrow(length = grid::unit(0.2, "cm")),
    colour = "gray30", linewidth = 0.4, inherit.aes = FALSE
  ) +
  ggplot2::geom_text(
    data = arrows_df,
    ggplot2::aes(x = .data[["xend"]], y = .data[["yend"]], label = .data[["variable"]]),
    colour = "gray30", size = 3, vjust = -0.4, inherit.aes = FALSE
  ) +
  ggplot2::geom_point(shape = 21, size = 3, colour = "white") +
  ggplot2::scale_fill_manual(values = color_vector) +
  ggplot2::labs(x = lab_x, y = lab_y, fill = "") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(face = "italic")
  ) +
  ggside::geom_xsidedensity(ggplot2::aes(fill = .data[["taxon"]]),
                            alpha = 0.5, lwd = 0.1, show.legend = FALSE) +
  ggside::geom_ysidedensity(ggplot2::aes(fill = .data[["taxon"]]),
                            alpha = 0.5, lwd = 0.1, show.legend = FALSE) +
  ggside::theme_ggside_void()

save_outputs(paste0("PCoA_", set_name), p, base_height = 8.5, base_aspect_ratio = 1.3)

# ---- Trait-correlation plots ---------------------------------------------------
p_corr1 <- make_corr_plot(trait_r[, 1], "PCo1", top_vars)
p_corr2 <- make_corr_plot(trait_r[, 2], "PCo2", top_vars)

save_outputs(paste0("Dim1_corr_variables_", set_name), p_corr1, base_height = 8.5, base_aspect_ratio = 1.3)
save_outputs(paste0("Dim2_corr_variables_", set_name), p_corr2, base_height = 8.5, base_aspect_ratio = 1.3)

# ---- Variation explained by the first axes ---------------------------------------
n_show <- min(10, ncol(scores))

eig_df <- data.frame(
  axis = factor(paste0("PCo", seq_len(n_show)), levels = paste0("PCo", seq_len(n_show))),
  pct  = 100 * rel_eig[seq_len(n_show)]
)

p_eig <- ggplot2::ggplot(eig_df, ggplot2::aes(x = .data[["axis"]], y = .data[["pct"]])) +
  ggplot2::geom_col(fill = "gray50") +
  ggplot2::labs(x = NULL, y = "Variation explained (%)") +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())

save_outputs(paste0("PCoA_eigenvalues_", set_name), p_eig, base_height = 5, base_aspect_ratio = 1.5)

# ---- PCo1xPCo2, PCo1xPCo3, PCo2xPCo3 scatter panels --------------------------------
if (ncol(scores) >= 3) {

  lab_3 <- sprintf("PCo3 (%.1f%%)", 100 * rel_eig[3])

  p_all_pc <- cowplot::plot_grid(
    make_pco_plot(pco_df, "PCo1", "PCo2", color_vector, lab_x, lab_y),
    make_pco_plot(pco_df, "PCo1", "PCo3", color_vector, lab_x, lab_3),
    make_pco_plot(pco_df, "PCo2", "PCo3", color_vector, lab_y, lab_3),
    labels = c("A", "B", "C"), ncol = 1, nrow = 3, align = "hv"
  )

  save_outputs(paste0("PCoA_panels_", set_name), p_all_pc,
               base_height = 7, base_aspect_ratio = 1.3, nrow = 3)
}

results[[set_name]] <- list(
  scores = pco_df,
  eigenvalues = data.frame(
    axis = paste0("PCo", seq_along(eig)),
    eigenvalue = eig,
    relative = rel_eig
  ),
  n_specimens = nrow(d),
  n_traits = ncol(d)
)

if (verbose) {
  message("  Done: ", set_name)
}
}

if (verbose) {
  message("PCoA analysis done. Outputs saved in: ", output_dir)
}

invisible(results)
}
