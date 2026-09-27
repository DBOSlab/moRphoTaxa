#' Hierarchical clustering dendrograms of morphological trait blocks
#'
#' @description
#' Builds hierarchical clustering dendrograms for morphological trait blocks
#' derived from a morphometric matrix, with automated selection of the
#' linkage method and number of clusters, and generates diagnostic and
#' dendrogram figures for each trait block analyzed. The function is designed
#' to work directly on the objects produced by [morph_matrix_setting()].
#'
#' @details
#' `morph_dendrogram()` expects `analysis_data` as produced by
#' [morph_matrix_setting()]: a numeric data frame with `id` and `taxon`
#' columns plus one column per morphological trait, with specimen IDs as row
#' names. The trait blocks (`"veg"`, `"flo"`, `"fru"`, or any subset thereof)
#' are not passed in separately — they are read directly from the
#' `"base_cols"` attribute that [morph_matrix_setting()] attaches to
#' `analysis_data`. If `analysis_data` has no `"base_cols"` attribute, the
#' function stops with an informative error.
#'
#' The `dendro_blocks` argument controls which trait sets are analyzed:
#'
#' \itemize{
#'   \item `"all"`: runs the analysis on each individual block in `base_cols`
#'   (typically `"veg"`, `"flo"`, `"fru"`), as well as on every pairwise and
#'   full combination of these blocks (`"vegflo"`, `"vegfru"`, `"flofru"`,
#'   `"vegflofru"`) that can be fully built from the blocks present;
#'   \item a character vector naming one or more blocks in `base_cols`: runs
#'   the analysis only on the requested block(s), individually.
#' }
#'
#' For each trait block, rows with missing values in the block's columns are
#' removed via [stats::na.omit()]. Blocks with fewer than 3 remaining
#' specimens or fewer than 2 traits are skipped with a warning.
#'
#' When a block has more than 2 traits, the traits used for clustering are
#' first narrowed to those contributing more than the average expected
#' contribution to PC1 (`100 / p`, where `p` is the number of traits in the
#' block) — the same reference-line criterion used by
#' [factoextra::fviz_contrib()] to flag "significant" variables. If fewer
#' than 2 traits pass this criterion, the top 2 contributing traits are used
#' instead, so that at least 2 traits are always retained for clustering.
#'
#' A distance matrix is computed on the retained traits with
#' [vegan::vegdist()] using `dist_method`. Each of `linkage_candidates` is fit
#' with [stats::hclust()], and the linkage whose cophenetic distances
#' correlate best with the original distance matrix (Sokal & Rohlf's
#' cophenetic correlation) is selected automatically. The number of clusters
#' `k` is then chosen automatically, from 2 up to `min(k_max, n - 1)`, as the
#' value maximizing the average silhouette width ([cluster::silhouette()]).
#'
#' For each analyzed block, the function produces:
#'
#' \itemize{
#'   \item a bar plot of cophenetic correlation by candidate linkage method;
#'   \item a line plot of average silhouette width by candidate `k`, with the
#'   selected `k` marked;
#'   \item a vertical dendrogram, with branches colored by the statistical
#'   cluster assignment (`k` groups) and leaf labels colored by known taxon,
#'   so the two can be visually compared;
#'   \item a horizontal version of the same dendrogram.
#' }
#'
#' Diagnostic plots are saved in the formats specified by `file_formats`.
#' Dendrograms are base-R plots and are always saved as both PDF and JPEG,
#' since `file_formats` controls [ggplot2::ggplot()]-based diagnostic output
#' only.
#'
#' Leaf labels in each dendrogram are colored by known taxon. These taxon
#' colours are taken from the `"taxon_colors"` attribute attached to
#' `analysis_data` by [morph_matrix_setting()], so that each taxon keeps the
#' same colour used in other plots (e.g. [morph_boxplots()], [morph_pca()]).
#' If this attribute is missing, or lacks an entry for a taxon present in a
#' given block, colours are generated automatically from the `viridis`
#' palette for the affected taxa, with a warning. This is independent of
#' `cluster_cols`, which colors the dendrogram *branches* by statistical
#' cluster rather than by taxon.
#'
#' Output files are written to a date-specific directory inside
#' `Figs.dendrogram`.
#'
#' @param analysis_data A numeric data frame containing morphological traits
#' as columns and specimens as rows, with `id` and `taxon` columns and
#' specimen identifiers as row names. Typically the output of
#' [morph_matrix_setting()]. If this object carries a `"taxon_colors"`
#' attribute, those colours are reused for the taxon leaf labels; otherwise a
#' `viridis` palette is generated automatically.
#'
#' @param dendro_blocks Character. Either `"all"`, to run the analysis on
#' every individual block in `base_cols` plus all pairwise and full
#' combinations of these blocks, or a character vector naming one or more
#' blocks in `base_cols` to analyze individually. Default is `"all"`.
#'
#' @param dist_method Character string giving the distance method passed to
#' [vegan::vegdist()] (e.g. `"bray"`, `"jaccard"`, `"euclidean"`,
#' `"manhattan"`, `"gower"`). Default is `"bray"`.
#'
#' @param linkage_candidates Character vector of linkage methods to compare,
#' passed in turn to [stats::hclust()]. Default is
#' `c("complete", "average", "ward.D2", "single")`.
#'
#' @param k_max Integer. Upper bound for the number of clusters considered
#' during the silhouette-width search. Automatically capped at `n - 1` for
#' each trait block. Default is `8`.
#'
#' @param cluster_cols Character vector of colors used for dendrogram
#' branches, one per cluster. If `NULL` (default), a palette is generated
#' automatically for each block, sized exactly to that block's selected
#' number of clusters (`best_k`), so clusters are never forced to repeat
#' colors. If a vector is supplied, it is recycled when the selected number
#' of clusters exceeds its length — a warning is issued if this happens,
#' since repeated colors can make distinct clusters visually
#' indistinguishable. This controls branch colours only — leaf-label
#' colours by taxon come from `analysis_data`'s `"taxon_colors"` attribute
#' (see Details).
#'
#' @param file_formats Character vector specifying which image formats to
#' save the ggplot-based diagnostic plots in. One or more of `"pdf"` and
#' `"jpeg"`. Default is `c("pdf", "jpeg")`. Dendrograms themselves are always
#' saved as both PDF and JPEG.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' including skipped blocks, selected linkage/k, and the output directory,
#' to the console.
#'
#' @return
#' Invisibly returns `NULL`.
#'
#' The function is primarily called for its side effects: diagnostic plots
#' and vertical/horizontal dendrograms are written to disk for each analyzed
#' trait block.
#'
#' @seealso
#' [morph_matrix_setting()],
#' [morph_pca()]
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
#' # Dendrograms for every individual block and every block combination
#' morph_dendrogram(analysis_data = res, dendro_blocks = "all")
#'
#' # Dendrogram for the vegetative block only, using Euclidean distance
#' morph_dendrogram(
#'   analysis_data = res,
#'   dendro_blocks = "veg",
#'   dist_method = "euclidean"
#' )
#' }
#'
#' @importFrom stats hclust cophenetic cutree na.omit prcomp setNames
#' @importFrom ggplot2 aes geom_col geom_line geom_point geom_text geom_vline
#'   ggplot labs theme_bw
#' @importFrom vegan vegdist
#' @importFrom cluster silhouette
#' @importFrom cowplot save_plot
#' @importFrom viridis viridis
#' @importFrom grDevices pdf jpeg dev.off
#' @importFrom graphics par plot legend
#'
#' @export

morph_dendrogram <- function(analysis_data,
                             dendro_blocks = "all",
                             dist_method = "bray",
                             linkage_candidates = c("complete", "average", "ward.D2", "single"),
                             k_max = 8,
                             cluster_cols = NULL,
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

if (!is.character(dendro_blocks) || length(dendro_blocks) < 1) {
stop("`dendro_blocks` must be \"all\" or a character vector of block names.")
}

if (!is.character(dist_method) || length(dist_method) != 1L) {
stop("`dist_method` must be a single character string.")
}

if (!is.character(linkage_candidates) || length(linkage_candidates) < 1) {
stop("`linkage_candidates` must be a character vector of linkage methods.")
}

if (!is.numeric(k_max) || length(k_max) != 1L || k_max < 2) {
stop("`k_max` must be a single integer of at least 2.")
}

file_formats <- match.arg(
file_formats,
choices = c("pdf", "jpeg"),
several.ok = TRUE
)

if (!is.logical(verbose) || length(verbose) != 1L) {
stop("`verbose` must be TRUE or FALSE.")
}

# ---- Species vector -----------------------------------------------------------

species_all <- stats::setNames(
as.character(analysis_data$taxon),
rownames(analysis_data)
)

# ---- Load dendextend quietly ---------------------------------------------------

# Not imported into the namespace; see .quiet_load_namespaces()
.quiet_load_namespaces("dendextend")

# ---- Output folder --------------------------------------------------------------

output_dir <- file.path("Figs.dendrogram", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (verbose) {
message("Running Dendrogram analysis")
}

# ---- Resolve which blocks to run -----------------------------------------------

if (identical(dendro_blocks, "all")) {
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
requested <- intersect(dendro_blocks, names(base_cols))
if (length(requested) == 0) {
  stop(
    "No valid blocks in `dendro_blocks`. Use \"all\" or any of: ",
    paste(names(base_cols), collapse = ", ")
  )
}
run_blocks <- base_cols[requested]
}

if (verbose) {
message("  Blocks to run: ", paste(names(run_blocks), collapse = ", "))
}

# ---- Dendrogram for each trait set -----------------------------------------------

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
  message("Running dendrogram for: ", set_name)
}

# Keeps traits whose contribution to PC1 exceeds the average expected
# contribution if all traits contributed equally (100/p %) -- the same
# reference-line criterion factoextra::fviz_contrib uses to flag
# "significant" variables.
if (ncol(d) > 2) {
  res.pca   <- stats::prcomp(d, scale. = TRUE)
  loadings1 <- res.pca$rotation[, "PC1"]
  contrib1  <- (loadings1^2 / sum(loadings1^2)) * 100
  avg_contrib <- 100 / length(contrib1)

  top_metrics <- names(contrib1)[contrib1 > avg_contrib]
  if (length(top_metrics) < 2) {
    top_metrics <- names(sort(contrib1, decreasing = TRUE))[1:2]
  }

  if (verbose) {
    message(
      "  Traits kept (PC1 contribution > ", round(avg_contrib, 1), "%): ",
      paste(top_metrics, collapse = ", ")
    )
  }
  d_top <- d[, top_metrics, drop = FALSE]
} else {
  d_top <- d
}

# ---- Distance matrix ------------------------------------------------------

dd <- vegan::vegdist(x = d_top, method = dist_method, na.rm = TRUE)

# Builds a tree with each candidate linkage method and measures how well
# its cophenetic distances preserve the original distance matrix (Sokal &
# Rohlf's cophenetic correlation). The linkage with the highest correlation
# is kept.
coph_cor <- sapply(linkage_candidates, function(m) {
  tryCatch({
    hc <- stats::hclust(dd, method = m)
    stats::cor(dd, stats::cophenetic(hc))
  }, error = function(e) NA_real_)
})
names(coph_cor) <- linkage_candidates

if (all(is.na(coph_cor))) {
  warning(
    "Skipping '", set_name,
    "': cophenetic correlation could not be computed for any linkage method."
  )
  next
}

if (verbose && any(is.na(coph_cor))) {
  message(
    "  Linkage method(s) skipped (undefined cophenetic correlation): ",
    paste(names(coph_cor)[is.na(coph_cor)], collapse = ", ")
  )
}

best_linkage <- names(coph_cor)[which.max(coph_cor)]

if (verbose) {
  message(
    "  Best linkage method: ", best_linkage,
    " (cophenetic correlation = ", round(max(coph_cor, na.rm = TRUE), 3), ")"
  )
}

res.den <- stats::hclust(d = dd, method = best_linkage)

# ---- Automated choice of k via average silhouette width --------------------

k_cap     <- min(k_max, nrow(d_top) - 1)
sil_width <- sapply(2:k_cap, function(k) {
  grp <- stats::cutree(res.den, k = k)
  mean(cluster::silhouette(grp, dd)[, "sil_width"])
})
best_k <- (2:k_cap)[which.max(sil_width)]

if (verbose) {
  message(
    "  Best k by average silhouette width: ", best_k,
    " (width = ", round(max(sil_width), 3), ")"
  )
}

# ---- Save diagnostic plots for transparency (linkage choice + silhouette/k) --

if (verbose) {
  message("  Cophenetic correlation by method for '", set_name, "':")
  print(coph_cor)
}

coph_df <- data.frame(
  method = factor(names(coph_cor), levels = linkage_candidates),
  cophenetic_cor = coph_cor
)
p_coph <- ggplot2::ggplot(coph_df, ggplot2::aes(x = .data[["method"]], y = .data[["cophenetic_cor"]])) +
  ggplot2::geom_col(fill = "gray70", na.rm = TRUE) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(is.na(.data[["cophenetic_cor"]]), "NA", round(.data[["cophenetic_cor"]], 3))),
    vjust = -0.4, size = 3
  ) +
  ggplot2::scale_x_discrete(drop = FALSE) +
  ggplot2::labs(
    title = paste0("Cophenetic correlation by linkage - ", set_name),
    x = "Linkage method", y = "Cophenetic correlation"
  ) +
  ggplot2::theme_bw()

sil_df <- data.frame(k = 2:k_cap, avg_sil_width = sil_width)
p_sil <- ggplot2::ggplot(sil_df, ggplot2::aes(x = .data[["k"]], y = .data[["avg_sil_width"]])) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_vline(xintercept = best_k, linetype = "dashed", color = "red") +
  ggplot2::labs(
    title = paste0("Silhouette width by k - ", set_name),
    x = "Number of clusters (k)", y = "Average silhouette width"
  ) +
  ggplot2::theme_bw()

if ("pdf" %in% file_formats) {
  ggplot2::ggsave(
    file.path(output_dir, paste0("linkage_selection_", set_name, ".pdf")),
    plot = p_coph, width = 6 * 1.3, height = 6, units = "in"
  )
  ggplot2::ggsave(
    file.path(output_dir, paste0("nbclust_silhouette_", set_name, ".pdf")),
    plot = p_sil, width = 6 * 1.3, height = 6, units = "in"
  )
}

if ("jpeg" %in% file_formats) {
  ggplot2::ggsave(
    file.path(output_dir, paste0("linkage_selection_", set_name, ".jpeg")),
    plot = p_coph, width = 6 * 1.3, height = 6, units = "in", dpi = 300
  )
  ggplot2::ggsave(
    file.path(output_dir, paste0("nbclust_silhouette_", set_name, ".jpeg")),
    plot = p_sil, width = 6 * 1.3, height = 6, units = "in", dpi = 300
  )
}

# ---- Color vector keyed by species -----------------------------------------

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

# ---- Build dendrogram: branches colored by statistical cluster (best_k),
# leaf labels colored by known taxon, so the two can be visually compared

if (is.null(cluster_cols)) {
  branch_cols <- grDevices::colorRampPalette(
    c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
      "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
  )(best_k)
} else {
  if (best_k > length(cluster_cols)) {
    warning(
      "`cluster_cols` has fewer colors (", length(cluster_cols),
      ") than the selected number of clusters for '", set_name, "' (",
      best_k, "). Colors will be recycled, and some clusters may share ",
      "the same color."
    )
  }
  branch_cols <- rep(cluster_cols, length.out = best_k)
}

dend <- stats::as.dendrogram(res.den)
dend <- dendextend::set(dend, "labels_cex", 0.5)
dend <- dendextend::set(dend, "leaves_pch", 19)
dend <- dendextend::color_branches(dend, k = best_k, col = branch_cols)

leaf_order   <- stats::order.dendrogram(dend)
dend_taxon   <- sp[leaf_order]
leaves_color <- color_vector[dend_taxon]
dend         <- dendextend::set(dend, "leaves_col", leaves_color)

plot_title <- paste0(
  dist_method, " distance, ", best_linkage, " linkage - ", set_name,
  " (k=", best_k, ")"
)

# ---- Vertical dendrogram ----------------------------------------------------

grDevices::pdf(file.path(output_dir, paste0("Dendrogram_", set_name, "_vert.pdf")),
               width = 10, height = 7)
graphics::plot(dend, main = plot_title)
graphics::legend("topright", inset = c(-0.2, 0), legend = names(color_vector),
                 pch = 19, col = color_vector, title = "Taxon", cex = 0.5, text.font = 3)
grDevices::dev.off()

grDevices::jpeg(file.path(output_dir, paste0("Dendrogram_", set_name, "_vert.jpeg")),
                width = 10, height = 7, units = "in", res = 300)
graphics::plot(dend, main = plot_title)
graphics::legend("topright", inset = c(-0.2, 0), legend = names(color_vector),
                 pch = 19, col = color_vector, title = "Taxon", cex = 0.5, text.font = 3)
grDevices::dev.off()

# ---- Horizontal dendrogram ---------------------------------------------------

grDevices::pdf(file.path(output_dir, paste0("Dendrogram_", set_name, "_hor.pdf")),
               width = 10, height = 7)
graphics::par(mar = c(5, 4, 4, 8), xpd = TRUE)
graphics::plot(dend, main = plot_title, horiz = TRUE)
graphics::legend("topright", inset = c(-0.2, 0), legend = names(color_vector),
                 pch = 19, col = color_vector, title = "Taxon", cex = 0.5, text.font = 3)
grDevices::dev.off()

grDevices::jpeg(file.path(output_dir, paste0("Dendrogram_", set_name, "_hor.jpeg")),
                width = 10, height = 7, units = "in", res = 300)
graphics::par(mar = c(5, 4, 4, 8), xpd = TRUE)
graphics::plot(dend, main = plot_title, horiz = TRUE)
graphics::legend("topright", inset = c(-0.2, 0), legend = names(color_vector),
                 pch = 19, col = color_vector, title = "Taxon", cex = 0.5, text.font = 3)
grDevices::dev.off()

if (verbose) {
  message("  Done: ", set_name)
}
}

if (verbose) {
message("Done. Outputs saved in: ", output_dir)
}

invisible(NULL)
}
