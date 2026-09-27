#' Mantel test between geographic and morphological distances
#'
#' @description
#' Tests whether morphological dissimilarity among specimens is correlated
#' with the geographic distance between the localities where they were
#' collected (isolation by distance), using a Mantel test.
#'
#' Takes the geographical matrix produced by [geo_matrix_setting()] and the
#' morphometric matrix produced by [morph_matrix_setting()], matches
#' specimens by their identifiers, computes one geographic and one
#' morphological distance matrix, runs [vegan::mantel()], prints a summary
#' of the result, and draws a scatter plot of pairwise morphological versus
#' geographic distances with a linear trend line.
#'
#' @details
#' \strong{Matching specimens.}
#' Specimens are matched by row names. Both [geo_matrix_setting()] and
#' [morph_matrix_setting()] use the specimen identifier (an existing
#' \code{"id"} column, or automatic identifiers \code{"ID0001"},
#' \code{"ID0002"}, ...) as row names, so the two matrices must be built
#' from the same specimens (typically the same spreadsheet, or spreadsheets
#' sharing an \code{"id"} column). Only specimens present in both matrices,
#' with valid coordinates and at least one morphological measurement, are
#' used. Specimens present in only one of the two matrices are ignored.
#'
#' \strong{Morphological matrix.}
#' The \code{"taxon"} column of \code{analysis_data} is \emph{not} part of
#' the morphological distance. Only the numeric trait columns are used, and
#' traits with no observation among the retained specimens are dropped.
#'
#' \strong{Geographic distance.}
#' With \code{geo_dist_method = "haversine"} (default), pairwise
#' great-circle distances are computed from \code{decimalLatitude} and
#' \code{decimalLongitude} on a sphere of radius 6371.0088 km, and are
#' expressed in kilometres. With \code{"euclidean"}, the Euclidean distance
#' between coordinates in decimal degrees is used instead; this is a rough
#' approximation that distorts distances away from the equator and across
#' wide longitudinal ranges.
#'
#' \strong{Morphological distance.}
#' \itemize{
#'   \item \code{"gower"} (default): computed with [cluster::daisy()]. Each
#'   trait is range-standardized and the dissimilarity between two specimens
#'   is averaged over the traits measured in both, so missing values are
#'   handled pairwise: no specimen is dropped and nothing is imputed. An
#'   error is raised if some pair of specimens shares no measured trait.
#'   \item \code{"bray"}, \code{"jaccard"}, \code{"euclidean"},
#'   \code{"manhattan"}: computed with [vegan::vegdist()]. These methods
#'   require complete data, so specimens with any missing trait value are
#'   dropped (with a warning). Traits are used as they are, without any
#'   standardization, so traits measured on larger numeric scales dominate
#'   \code{"euclidean"} and \code{"manhattan"} distances; \code{"bray"} and
#'   \code{"jaccard"} additionally assume non-negative values.
#' }
#'
#' \strong{Mantel test.}
#' The significance of the correlation is assessed by permuting the rows and
#' columns of one distance matrix, so it accounts for the non-independence
#' among pairwise distances. The linear trend line in the plot is purely
#' descriptive and its confidence band should not be interpreted as a test.
#' Because the permutations are random, call [set.seed()] before this
#' function to obtain reproducible p-values.
#'
#' \strong{Outputs.}
#' The plot is saved in a dated subfolder
#' (\code{format(Sys.time(), "\%d\%b\%Y")}) inside \code{"Figs_Mantel"}, in
#' the formats selected by \code{file_formats}. The dated subfolder means
#' repeated runs do not overwrite each other.
#'
#' @param geodata A data frame produced by [geo_matrix_setting()]. It must
#' contain the \code{decimalLatitude} and \code{decimalLongitude} columns,
#' and its row names must be the specimen identifiers.
#'
#' @param analysis_data A data frame produced by [morph_matrix_setting()].
#' Its row names must be the specimen identifiers, and it must contain the
#' numeric trait columns (a \code{"taxon"} column, if present, is ignored).
#'
#' @param cor_method Character string giving the correlation coefficient
#' computed between the two distance matrices. One of \code{"pearson"}
#' (default), \code{"spearman"}, or \code{"kendall"}. Passed to the
#' \code{method} argument of [vegan::mantel()].
#'
#' @param permutations Integer. Number of permutations used to assess
#' significance. Default is \code{999}.
#'
#' @param dist_method Character string giving the morphological distance.
#' One of \code{"gower"} (default), \code{"bray"}, \code{"jaccard"},
#' \code{"euclidean"}, or \code{"manhattan"}. See Details.
#'
#' @param geo_dist_method Character string giving the geographic distance.
#' One of \code{"haversine"} (default; great-circle distance in km) or
#' \code{"euclidean"} (Euclidean distance in decimal degrees). See Details.
#'
#' @param file_formats Character vector specifying the file formats used to
#' save the plot. One or more of \code{"pdf"} and \code{"jpg"}. Default is
#' \code{c("pdf", "jpg")}.
#'
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages
#' and a summary of the Mantel test results.
#'
#' @return
#' A list, returned invisibly, with the elements:
#' \describe{
#'   \item{\code{mantel}}{The [vegan::mantel()] result (class
#'   \code{"mantel"}), including the Mantel statistic (\code{$statistic})
#'   and the permutational p-value (\code{$signif}).}
#'   \item{\code{dist_geo}}{The geographic distance matrix (class
#'   \code{"dist"}).}
#'   \item{\code{dist_morpho}}{The morphological distance matrix (class
#'   \code{"dist"}).}
#'   \item{\code{plot}}{The [ggplot2::ggplot()] object of the scatter plot.}
#'   \item{\code{n}}{Number of specimens used in the test.}
#'   \item{\code{output_dir}}{Path of the folder where the plot was saved.}
#' }
#' The function is also called for its side effects: the plot is written to
#' \code{output_dir} and, when \code{verbose = TRUE}, a summary is printed to
#' the console.
#'
#' @seealso
#' [geo_matrix_setting()], [morph_matrix_setting()],
#' [geo_bioclim_exploratory()]
#'
#' @examples
#' \dontrun{
#' # Both matrices must be built from the same specimens
#' geo <- geo_matrix_setting(
#'   xlsx_path = "specimen_data.xlsx",
#'   taxon_col = "species",
#'   edaphic = FALSE
#' )
#'
#' morph <- morph_matrix_setting(
#'   xlsx_path = "specimen_data.xlsx",
#'   taxon_col = "species",
#'   trait_name_type = "code",
#'   veg_first = "petiole_length/PETIlng",
#'   flo_first = "inflorescence_length/INFLlng",
#'   fru_first = "fruit_stipe_length/FRSTlng"
#' )
#'
#' # Default: Pearson correlation, 999 permutations, Gower distance
#' set.seed(123)
#' res <- geo_morph_mantel(geodata = geo, analysis_data = morph)
#' res$mantel
#'
#' # Rank correlation, more permutations, PDF only
#' res <- geo_morph_mantel(
#'   geodata = geo,
#'   analysis_data = morph,
#'   cor_method = "spearman",
#'   permutations = 9999,
#'   file_formats = "pdf"
#' )
#'
#' # Complete-case Euclidean morphological distance, as in earlier scripts
#' # that used decimal degrees for the geographic distance
#' res <- geo_morph_mantel(
#'   geodata = geo,
#'   analysis_data = morph,
#'   dist_method = "euclidean",
#'   geo_dist_method = "euclidean"
#' )
#' }
#'
#' @importFrom stats as.dist complete.cases dist
#' @importFrom vegan mantel vegdist
#' @importFrom cluster daisy
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth annotate labs
#' @importFrom ggplot2 theme_bw theme element_blank
#' @importFrom cowplot save_plot
#'
#' @export

geo_morph_mantel <- function(geodata,
                             analysis_data,
                             cor_method = c("pearson", "spearman", "kendall"),
                             permutations = 999,
                             dist_method = c("gower", "bray", "jaccard",
                                             "euclidean", "manhattan"),
                             geo_dist_method = c("haversine", "euclidean"),
                             file_formats = c("pdf", "jpg"),
                             verbose = TRUE) {

# ---- Validate input -------------------------------------------------------

if (!is.data.frame(geodata)) {
  stop("`geodata` must be a data frame, typically the output of geo_matrix_setting().")
}

if (!is.data.frame(analysis_data)) {
  stop("`analysis_data` must be a data frame, typically the output of morph_matrix_setting().")
}

if (!all(c("decimalLatitude", "decimalLongitude") %in% names(geodata))) {
  stop("`geodata` must contain `decimalLatitude` and `decimalLongitude` columns.")
}

cor_method <- match.arg(cor_method)
dist_method <- match.arg(dist_method)
geo_dist_method <- match.arg(geo_dist_method)

if (!is.numeric(permutations) || length(permutations) != 1L ||
    is.na(permutations) || permutations < 1 ||
    permutations != round(permutations)) {
  stop("`permutations` must be a single positive integer.")
}

if (!is.character(file_formats) || length(file_formats) < 1L) {
  stop("`file_formats` must be a character vector with at least one format.")
}

file_formats <- match.arg(
  file_formats,
  choices = c("pdf", "jpg"),
  several.ok = TRUE
)

if (!is.logical(verbose) || length(verbose) != 1L) {
  stop("`verbose` must be TRUE or FALSE.")
}

# ---- Match specimens between the two matrices ------------------------------

common_ids <- intersect(rownames(geodata), rownames(analysis_data))

if (length(common_ids) == 0L) {
  stop(
    "No specimen identifiers are shared by `geodata` and `analysis_data`. ",
    "Both must be built from the same specimens (row names are the specimen IDs)."
  )
}

coords <- geodata[
  common_ids,
  c("decimalLongitude", "decimalLatitude"),
  drop = FALSE
]
coords$decimalLongitude <- as.numeric(coords$decimalLongitude)
coords$decimalLatitude <- as.numeric(coords$decimalLatitude)

# Numeric trait columns only: `taxon` (and any other non-numeric column)
# must not enter the morphological distance.
trait_cols <- setdiff(names(analysis_data), "taxon")
trait_cols <- trait_cols[
  vapply(analysis_data[trait_cols], is.numeric, logical(1))
]

if (length(trait_cols) == 0L) {
  stop("`analysis_data` contains no numeric trait columns.")
}

morpho <- analysis_data[common_ids, trait_cols, drop = FALSE]

# Drop specimens without valid coordinates
ok_coords <- stats::complete.cases(coords)
coords <- coords[ok_coords, , drop = FALSE]
morpho <- morpho[ok_coords, , drop = FALSE]

# Drop traits with no observation, then specimens with no observation
morpho <- morpho[, colSums(!is.na(morpho)) > 0, drop = FALSE]
has_obs <- rowSums(!is.na(morpho)) > 0
coords <- coords[has_obs, , drop = FALSE]
morpho <- morpho[has_obs, , drop = FALSE]

# ---- Complete-case filtering for methods that cannot handle NAs ------------

if (dist_method != "gower") {
  incomplete <- !stats::complete.cases(morpho)
  if (any(incomplete)) {
    warning(
      sum(incomplete), " specimens with NAs dropped for method '",
      dist_method, "'. Consider dist_method = \"gower\"."
    )
    morpho <- morpho[!incomplete, , drop = FALSE]
    coords <- coords[!incomplete, , drop = FALSE]
  }
}

if (nrow(coords) < 3L) {
  stop("Fewer than 3 specimens with valid coordinates and morphological data remain.")
}

if (verbose) {
  message("Running Mantel test...")
  message("  Specimens used: ", nrow(coords))
  message("  Morphological traits used: ", ncol(morpho))
}

# ---- Distance matrices -----------------------------------------------------

if (geo_dist_method == "haversine") {
  dist_geo <- .haversine_dist(
    lon = coords$decimalLongitude,
    lat = coords$decimalLatitude,
    labels = rownames(coords)
  )
} else {
  dist_geo <- stats::dist(coords[, c("decimalLatitude", "decimalLongitude")])
}

if (dist_method == "gower") {
  # Gower handles NAs pairwise: no specimens dropped, no imputation needed
  dist_morpho <- stats::as.dist(cluster::daisy(morpho, metric = "gower"))
} else {
  dist_morpho <- vegan::vegdist(morpho, method = dist_method)
}

if (anyNA(dist_morpho)) {
  stop(
    "The morphological distance matrix contains NA values: some pairs of ",
    "specimens share no measured trait. Remove sparsely measured specimens ",
    "or traits before running the Mantel test."
  )
}

# ---- Mantel test -----------------------------------------------------------

mantel_result <- vegan::mantel(
  dist_geo, dist_morpho,
  method = cor_method,
  permutations = permutations
)

if (verbose) {
  message("\n--- Mantel test results ---")
  message("  Specimens          : ", nrow(coords))
  message("  Correlation method : ", cor_method)
  message("  Morpho distance    : ", dist_method)
  message("  Geographic distance: ", geo_dist_method)
  message("  Permutations       : ", permutations)
  message("  Statistic (r)      : ", round(mantel_result$statistic, 4))
  message("  p-value            : ", round(mantel_result$signif, 4), "\n")
}

# ---- Plot ------------------------------------------------------------------

sig_label <- if (mantel_result$signif < 0.001) {
  "p < 0.001"
} else {
  paste0("p = ", round(mantel_result$signif, 3))
}

plot_data <- data.frame(
  geo_dist = as.vector(dist_geo),
  morpho_dist = as.vector(dist_morpho)
)

geo_axis_label <- if (geo_dist_method == "haversine") {
  "Geographic distance (km)"
} else {
  "Geographic distance (decimal degrees)"
}

p_mantel <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = .data$geo_dist, y = .data$morpho_dist)
) +
  ggplot2::geom_point(alpha = 0.4, size = 1.2, color = "steelblue") +
  ggplot2::geom_smooth(
    method = "lm", formula = y ~ x,
    se = TRUE, color = "red", linetype = "dashed", linewidth = 0.8
  ) +
  ggplot2::annotate(
    "text",
    x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
    label = paste0("r = ", round(mantel_result$statistic, 3), "\n", sig_label),
    size = 3.5
  ) +
  ggplot2::labs(
    x = geo_axis_label,
    y = paste0("Morphological distance (", dist_method, ")")
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  )

# ---- Save outputs ----------------------------------------------------------

output_dir <- file.path("Figs_Mantel", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (ext in file_formats) {
  cowplot::save_plot(
    file.path(output_dir, paste0("Mantel_plot.", ext)),
    p_mantel,
    ncol = 1, nrow = 1, base_height = 6, base_aspect_ratio = 1.3
  )
}

if (verbose) {
  message("Done. Outputs saved in: ", output_dir)
}

invisible(list(
  mantel = mantel_result,
  dist_geo = dist_geo,
  dist_morpho = dist_morpho,
  plot = p_mantel,
  n = nrow(coords),
  output_dir = output_dir
))
}

#' Great-circle (haversine) distance matrix
#'
#' @description
#' Internal helper used by [geo_morph_mantel()] to compute pairwise
#' great-circle distances between points given in decimal degrees.
#'
#' @param lon Numeric vector of longitudes (decimal degrees).
#' @param lat Numeric vector of latitudes (decimal degrees).
#' @param labels Optional character vector of point labels, used as the
#' dimension names of the result.
#' @param radius Numeric. Sphere radius in km. Default is the mean Earth
#' radius, 6371.0088 km.
#'
#' @return A `"dist"` object with distances in km.
#'
#' @keywords internal
#'
#' @noRd
.haversine_dist <- function(lon, lat, labels = NULL, radius = 6371.0088) {
rad <- pi / 180
lon <- lon * rad
lat <- lat * rad

dlat <- outer(lat, lat, "-")
dlon <- outer(lon, lon, "-")

a <- sin(dlat / 2)^2 + outer(cos(lat), cos(lat)) * sin(dlon / 2)^2
s <- sqrt(a)
s[s > 1] <- 1   # guard against rounding error above 1
d <- 2 * radius * asin(s)

if (!is.null(labels)) {
  dimnames(d) <- list(labels, labels)
}

stats::as.dist(d)
}
