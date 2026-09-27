#' K-means clustering on PCA scores of morphological trait blocks
#'
#' @description
#' Performs principal component analysis (PCA) on morphological trait blocks
#' and applies K-means clustering to the resulting principal component scores,
#' generating a standard set of diagnostic figures for each trait block
#' analyzed. The function is designed to work directly on the objects produced
#' by [morph_matrix_setting()].
#'
#' For each requested trait block, the function fits a PCA, determines the
#' number of clusters (either fixed by the user or estimated automatically via
#' the gap statistic), runs K-means on the selected principal components, saves
#' a cluster plot in which point fill encodes taxon identity and dashed
#' ellipses encode cluster membership, saves an equivalent
#' \pkg{factoextra} cluster plot, and — when the number of clusters is
#' estimated — saves an elbow plot of the within-cluster sum of squares.
#'
#' @details
#' `morph_kmeans()` expects `analysis_data` as produced by
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
#' The `kmeans_blocks` argument controls which trait sets are analyzed:
#'
#' \itemize{
#'   \item `"all"`: runs K-means on each individual block in `base_cols`
#'   (typically `"veg"`, `"flo"`, `"fru"`), as well as on every pairwise and
#'   full combination of these blocks (`"vegflo"`, `"vegfru"`, `"flofru"`,
#'   `"vegflofru"`). Combinations whose source blocks are absent from
#'   `base_cols` are silently skipped;
#'   \item a character vector naming one or more blocks in `base_cols`: runs
#'   K-means only on the requested block(s), individually.
#' }
#'
#' For each trait block, rows with missing values in the relevant columns are
#' removed via [stats::na.omit()], and traits with zero variance (constant
#' across the retained specimens) are dropped, since they cannot contribute to
#' a correlation-based PCA. Blocks with fewer than 3 specimens or fewer than 2
#' traits after these filtering steps are skipped with a warning, as are blocks
#' in which any of the components requested in `kmeans_pcs` is unavailable.
#'
#' PCA is fit with [stats::prcomp()] using `scale. = TRUE` (i.e., on the
#' correlation matrix), matching [morph_pca()]. Clustering is then performed
#' with [stats::kmeans()] on the columns named in `kmeans_pcs` only.
#'
#' The number of clusters is resolved as follows:
#'
#' \itemize{
#'   \item if `kmeans_k` is a number, that value is used directly for every
#'   block;
#'   \item if `kmeans_k` is `NULL` (default), the gap statistic is computed
#'   with [cluster::clusGap()] over `1:kmeans_k_max` clusters using
#'   `kmeans_B` bootstrap replicates, and the number of clusters maximizing
#'   the gap is retained, constrained to at least 2. Within each block,
#'   `kmeans_k_max` is additionally capped at `nrow - 1` so that the gap
#'   statistic remains computable in small blocks.
#' }
#'
#' Cluster estimation and fitting are randomized; `kmeans_seed` is passed to
#' [set.seed()] before each block is processed, so that repeated calls with
#' the same data and arguments return identical cluster assignments.
#'
#' Taxon colours are taken from the `"taxon_colors"` attribute attached to
#' `analysis_data` by [morph_matrix_setting()], so that each taxon keeps the
#' same colour used in other plots (e.g. [morph_pca()], [morph_boxplots()]).
#' If this attribute is missing, or lacks an entry for a taxon present in a
#' given block, colours are generated automatically from the `viridis` palette
#' for the affected taxa, with a warning. Cluster ellipses are drawn in
#' greyscale so that taxon identity (fill) and cluster membership (outline)
#' remain visually distinguishable; ellipses are drawn only for clusters with
#' at least 4 specimens, since [ggplot2::stat_ellipse()] cannot be estimated
#' from fewer points.
#'
#' For each analyzed block, the following figures are generated:
#'
#' \itemize{
#'   \item a K-means cluster plot on the first two components listed in
#'   `kmeans_pcs`, with points filled by taxon and dashed t-distribution
#'   ellipses per cluster, and axis labels reporting the proportion of
#'   variance explained by each component;
#'   \item an equivalent cluster plot produced by
#'   [factoextra::fviz_cluster()];
#'   \item an elbow plot of the within-cluster sum of squares against the
#'   number of clusters, with the retained number of clusters marked — only
#'   when `kmeans_k` is `NULL`.
#' }
#'
#' Output files are written to a date-specific directory inside `Figs.KMeans`.
#' Each figure is saved in the formats specified by `file_formats`. File names
#' are prefixed by figure type (e.g. `"KMeans_plot_"`, `"KMeans_fviz_"`,
#' `"KMeans_elbow_"`) and suffixed by the trait block name.
#'
#' @param analysis_data A numeric data frame containing morphological traits
#' as columns and specimens as rows, with `id` and `taxon` columns and
#' specimen identifiers as row names. Typically the object returned by
#' [morph_matrix_setting()]. If it carries a `"taxon_colors"` attribute,
#' those colours are reused; otherwise a `viridis` palette is generated
#' automatically.
#'
#' @param kmeans_blocks Character. Either `"all"`, to cluster every individual
#' block in `base_cols` plus all pairwise and full combinations of these
#' blocks, or a character vector naming one or more blocks in `base_cols` to
#' analyze individually. Default is `"all"`.
#'
#' @param kmeans_k Integer or `NULL`. Fixed number of clusters. If `NULL`
#' (default), the number of clusters is estimated for each block from the gap
#' statistic.
#'
#' @param kmeans_k_max Integer. Maximum number of clusters evaluated when
#' `kmeans_k` is `NULL`. Default is `10`. Capped internally at `nrow - 1`
#' within each block.
#'
#' @param kmeans_B Integer. Number of bootstrap replicates used by
#' [cluster::clusGap()] when estimating the number of clusters. Default is
#' `50`.
#'
#' @param kmeans_nstart Integer. Number of random starting configurations
#' passed to [stats::kmeans()]. Default is `25`.
#'
#' @param kmeans_pcs Character vector naming the principal components used for
#' clustering, e.g. `c("PC1", "PC2")` (default). At least two components must
#' be given; the first two are also the axes of the cluster plots.
#'
#' @param kmeans_seed Integer. Random seed applied before clustering each
#' block, for reproducible cluster assignments. Default is `123`.
#'
#' @param file_formats Character vector specifying which image formats to
#' save. One or more of `"pdf"` and `"jpeg"`. Default is
#' `c("pdf", "jpeg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' per-block cluster summaries, skipped blocks, and the output directory to
#' the console.
#'
#' @return
#' Invisibly returns a data frame of cluster assignments with one row per
#' specimen per analyzed block, or `NULL` if no block could be analyzed. The
#' columns are:
#' \describe{
#'   \item{`block`}{Name of the trait block.}
#'   \item{`id`}{Specimen identifier.}
#'   \item{`taxon`}{Taxon assigned to the specimen.}
#'   \item{`cluster`}{Cluster assigned by K-means, as a factor.}
#'   \item{`optimal_k`}{Number of clusters used for that block.}
#' }
#'
#' The function is also called for its side effects: cluster plots, and
#' (when applicable) elbow plots, are written to disk for each analyzed trait
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
#' # Cluster every individual block and every block combination,
#' # estimating the number of clusters from the gap statistic
#' km <- morph_kmeans(
#'   analysis_data = res,
#'   kmeans_blocks = "all"
#' )
#'
#' # Inspect the cluster assignments
#' head(km)
#' table(km$block, km$cluster)
#'
#' # Cluster the vegetative block only, with a fixed number of clusters
#' morph_kmeans(
#'   analysis_data = res,
#'   kmeans_blocks = "veg",
#'   kmeans_k = 3
#' )
#'
#' # Cluster on the first three components, saving PDF output only
#' morph_kmeans(
#'   analysis_data = res,
#'   kmeans_blocks = c("veg", "flo"),
#'   kmeans_pcs = c("PC1", "PC2", "PC3"),
#'   file_formats = "pdf"
#' )
#' }
#'
#' @importFrom stats prcomp kmeans na.omit setNames var
#' @importFrom ggplot2 aes element_blank element_text geom_hline geom_line
#'   geom_point geom_vline ggplot labs scale_colour_grey scale_fill_manual
#'   scale_x_continuous stat_ellipse theme theme_bw
#' @importFrom cluster clusGap
#' @importFrom cowplot save_plot
#' @importFrom viridis viridis
#' @importFrom rlang .data
#'
#' @export

morph_kmeans <- function(analysis_data,
                         kmeans_blocks = "all",
                         kmeans_k = NULL,
                         kmeans_k_max = 10,
                         kmeans_B = 50,
                         kmeans_nstart = 25,
                         kmeans_pcs = c("PC1", "PC2"),
                         kmeans_seed = 123,
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

if (!is.character(kmeans_blocks) || length(kmeans_blocks) < 1) {
stop("`kmeans_blocks` must be \"all\" or a character vector of block names.")
}

if (!is.null(kmeans_k)) {
if (!is.numeric(kmeans_k) || length(kmeans_k) != 1L || kmeans_k < 2) {
  stop("`kmeans_k` must be NULL or a single integer >= 2.")
}
}

if (!is.numeric(kmeans_k_max) || length(kmeans_k_max) != 1L || kmeans_k_max < 2) {
stop("`kmeans_k_max` must be a single integer >= 2.")
}

if (!is.numeric(kmeans_B) || length(kmeans_B) != 1L || kmeans_B < 1) {
stop("`kmeans_B` must be a positive integer.")
}

if (!is.numeric(kmeans_nstart) || length(kmeans_nstart) != 1L || kmeans_nstart < 1) {
stop("`kmeans_nstart` must be a positive integer.")
}

if (!is.character(kmeans_pcs) || length(kmeans_pcs) < 2) {
stop("`kmeans_pcs` must be a character vector with at least two PC names, ",
     "e.g. c(\"PC1\", \"PC2\").")
}

if (!is.numeric(kmeans_seed) || length(kmeans_seed) != 1L) {
stop("`kmeans_seed` must be a single number.")
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

output_dir <- file.path("Figs.KMeans", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (verbose) {
message("Running K-means clustering")
}

# ---- Resolve blocks ----------------------------------------------------------

if (identical(kmeans_blocks, "all")) {
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
requested <- intersect(kmeans_blocks, names(base_cols))
if (length(requested) == 0) {
  stop(
    "No valid blocks in `kmeans_blocks`. Use \"all\" or any of: ",
    paste(names(base_cols), collapse = ", ")
  )
}
run_blocks <- base_cols[requested]
}

if (verbose) {
message("  Blocks to run: ", paste(names(run_blocks), collapse = ", "))
}

# ---- Helper: taxon colours ---------------------------------------------------

resolve_colors <- function(sp) {

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

color_vector
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

# ---- Helper: K-means for one trait block ------------------------------------

run_kmeans_block <- function(block_name, block_columns) {

valid_cols <- intersect(block_columns, names(analysis_data))
block_data <- stats::na.omit(analysis_data[, valid_cols, drop = FALSE])

# Remove zero-variance traits
keep <- vapply(block_data, function(x) {
  v <- stats::var(x, na.rm = TRUE)
  !is.na(v) && v > 0
}, logical(1))

block_data <- block_data[, keep, drop = FALSE]

if (nrow(block_data) < 3 || ncol(block_data) < 2) {
  warning(
    "Skipping '", block_name, "': not enough data after removing missing ",
    "values and constant variables (n = ", nrow(block_data),
    ", p = ", ncol(block_data), ")."
  )
  return(NULL)
}

sp <- species_all[rownames(block_data)]

if (verbose) {
  message("Running K-means for: ", block_name)
  message("  Specimens: ", nrow(block_data), " | Traits: ", ncol(block_data))
}

# ---- PCA -------------------------------------------------------------------

res_pca <- stats::prcomp(block_data, scale. = TRUE)
explained_var <- summary(res_pca)$importance[2, ]

missing_pcs <- setdiff(kmeans_pcs, colnames(res_pca$x))

if (length(missing_pcs) > 0) {
  warning(
    "Skipping '", block_name, "': component(s) not available: ",
    paste(missing_pcs, collapse = ", ")
  )
  return(NULL)
}

pca_df <- as.data.frame(res_pca$x[, kmeans_pcs, drop = FALSE])
pca_df$taxon <- sp[rownames(pca_df)]

set.seed(kmeans_seed)

# ---- Number of clusters ----------------------------------------------------

if (!is.null(kmeans_k)) {

  if (kmeans_k >= nrow(pca_df)) {
    warning(
      "Skipping '", block_name, "': `kmeans_k` (", kmeans_k,
      ") must be smaller than the number of specimens (", nrow(pca_df), ")."
    )
    return(NULL)
  }

  optimal_k <- kmeans_k

  if (verbose) {
    message("  Using fixed k = ", optimal_k)
  }

} else {

  k_max_eff <- min(kmeans_k_max, nrow(pca_df) - 1)

  if (k_max_eff < 2) {
    warning(
      "Skipping '", block_name, "': too few specimens (n = ", nrow(pca_df),
      ") to estimate the number of clusters."
    )
    return(NULL)
  }

  if (verbose) {
    message("  Estimating optimal k via gap statistic (K.max = ", k_max_eff, ")")
  }

  gap_stat <- cluster::clusGap(
    pca_df[, kmeans_pcs, drop = FALSE],
    FUN = stats::kmeans,
    nstart = kmeans_nstart,
    K.max = k_max_eff,
    B = kmeans_B
  )

  optimal_k <- which.max(gap_stat$Tab[, "gap"])
  optimal_k <- max(2, min(optimal_k, k_max_eff))

  if (verbose) {
    message("  Optimal k = ", optimal_k)
  }

  # ---- Elbow plot ----------------------------------------------------------

  wcss <- vapply(
    seq_len(k_max_eff),
    function(k) {
      stats::kmeans(
        pca_df[, kmeans_pcs, drop = FALSE],
        centers = k,
        nstart = kmeans_nstart
      )$tot.withinss
    },
    numeric(1)
  )

  p_elbow <- ggplot2::ggplot(
    data.frame(k = seq_len(k_max_eff), wcss = wcss),
    ggplot2::aes(x = .data[["k"]], y = .data[["wcss"]])
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(
      xintercept = optimal_k,
      linetype = "dashed",
      color = "red"
    ) +
    ggplot2::scale_x_continuous(breaks = seq_len(k_max_eff)) +
    ggplot2::labs(
      x = "Number of clusters (k)",
      y = "Within-cluster sum of squares",
      title = paste0("Elbow - ", block_name)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  save_outputs(
    paste0("KMeans_elbow_", block_name),
    p_elbow,
    base_height = 5,
    base_aspect_ratio = 1.4
  )
}

# ---- K-means ---------------------------------------------------------------

kmeans_result <- stats::kmeans(
  pca_df[, kmeans_pcs, drop = FALSE],
  centers = optimal_k,
  nstart = kmeans_nstart
)

pca_df$cluster <- as.factor(kmeans_result$cluster)

if (verbose) {
  message("  Cluster sizes: ", paste(table(pca_df$cluster), collapse = " | "))
  for (pc in kmeans_pcs) {
    message(sprintf("  %s: %.1f%%", pc, explained_var[pc] * 100))
  }
}

# ---- Taxon colours ---------------------------------------------------------

color_vector <- resolve_colors(sp)

# ---- Cluster plot ----------------------------------------------------------

pc_x <- kmeans_pcs[1]
pc_y <- kmeans_pcs[2]

# Ellipses require enough points per cluster to be estimated
big_clusters <- names(which(table(pca_df$cluster) >= 4))
ell_df <- pca_df[pca_df$cluster %in% big_clusters, , drop = FALSE]

p_kmeans <- ggplot2::ggplot(
  pca_df,
  ggplot2::aes(x = .data[[pc_x]], y = .data[[pc_y]])
)

if (nrow(ell_df) > 0) {
  p_kmeans <- p_kmeans +
    ggplot2::stat_ellipse(
      data = ell_df,
      ggplot2::aes(
        group = .data[["cluster"]],
        colour = .data[["cluster"]]
      ),
      type = "t",
      linetype = "dashed",
      linewidth = 0.6
    ) +
    ggplot2::scale_colour_grey(start = 0.2, end = 0.7)
}

p_kmeans <- p_kmeans +
  ggplot2::geom_point(
    ggplot2::aes(fill = .data[["taxon"]]),
    shape = 21,
    size = 3,
    colour = "white"
  ) +
  ggplot2::scale_fill_manual(values = color_vector) +
  ggplot2::geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::labs(
    x = paste0(pc_x, " (", round(explained_var[pc_x] * 100, 1), "%)"),
    y = paste0(pc_y, " (", round(explained_var[pc_y] * 100, 1), "%)"),
    fill = "",
    colour = "Cluster",
    title = paste0("K-means (k = ", optimal_k, ") - ", block_name)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(face = "italic")
  )

save_outputs(
  paste0("KMeans_plot_", block_name),
  p_kmeans,
  base_height = 7,
  base_aspect_ratio = 1.3
)

# ---- factoextra cluster plot -----------------------------------------------

.quiet_load_namespaces("factoextra")
p_fviz <- factoextra::fviz_cluster(
  kmeans_result,
  data = pca_df[, kmeans_pcs, drop = FALSE],
  palette = "jco",
  ggtheme = ggplot2::theme_bw(),
  main = paste0("K-means (k = ", optimal_k, ") - ", block_name)
)

save_outputs(
  paste0("KMeans_fviz_", block_name),
  p_fviz,
  base_height = 7,
  base_aspect_ratio = 1.3
)

if (verbose) {
  message("  Done: ", block_name)
}

# ---- Cluster assignments ---------------------------------------------------

data.frame(
  block = block_name,
  id = rownames(pca_df),
  taxon = pca_df$taxon,
  cluster = pca_df$cluster,
  optimal_k = optimal_k,
  row.names = NULL,
  stringsAsFactors = FALSE
)
}

# ---- Run all blocks ----------------------------------------------------------

all_assignments <- Filter(
Negate(is.null),
mapply(
  run_kmeans_block,
  names(run_blocks),
  run_blocks,
  SIMPLIFY = FALSE
)
)

if (length(all_assignments) == 0) {
warning("No trait block could be analyzed.")
return(invisible(NULL))
}

assignments <- do.call(rbind, all_assignments)
rownames(assignments) <- NULL

if (verbose) {
message("K-means clustering done. Outputs saved in: ", output_dir)
}

invisible(assignments)
}
