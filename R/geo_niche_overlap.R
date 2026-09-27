#' Pairwise niche overlap between taxa in a combined environmental space
#'
#' @description
#' Quantifies the environmental niche overlap between every pair of taxa in the
#' geographical/environmental matrix produced by [geo_matrix_setting()]. Climate
#' variables (and, when available, edaphic variables) at the occurrence points
#' are combined in a single principal component analysis (PCA), and niche
#' overlap is measured in the plane of the first two axes with Schoener's D,
#' following the PCA-env framework implemented in \pkg{ecospat}. For each pair,
#' niche equivalency and niche similarity are tested by permutation, and niche
#' dynamics indices (expansion, stability and unfilling) are computed.
#'
#' Environmental variables are taken automatically from the `"env_cols"`
#' attribute of `geodata`, so climate-only and climate + edaphic matrices are
#' both handled without further arguments.
#'
#' @details
#' \strong{Input.} `geodata` must be the data frame returned by
#' [geo_matrix_setting()], carrying the `"env_cols"` attribute. It must contain
#' the column `taxon` and the environmental columns listed in
#' `attr(geodata, "env_cols")`. The climate block (`"clim"`, which includes
#' `"elevation"`) is always used; the edaphic block (`"edaphic"`) is added when
#' present. All taxa in `geodata` are compared; to restrict the analysis, use
#' `species_selected` in [geo_matrix_setting()].
#'
#' \strong{Environmental space.} The workflow is:
#' \enumerate{
#'   \item Environmental variables are converted to numeric and records with
#'   any missing value are removed (at least 20 complete records are
#'   required).
#'   \item Variables with zero variance are dropped.
#'   \item Highly collinear variables (|r| > `cor_cutoff`) are removed with
#'   [caret::findCorrelation()]. The variables removed are documented in
#'   `Removed_correlated_variables.xlsx`.
#'   \item A single PCA (centred and scaled) is fitted to the remaining
#'   variables. The pooled records of all taxa define the environmental space
#'   ("background") in which every pairwise comparison is made, so overlap
#'   values are comparable across pairs.
#' }
#'
#' The niche analyses use the first two principal components (PC1 and PC2).
#' `var_threshold` does not change which axes are used: it defines the
#' cumulative-variance line on the scree plot and the number of axes reported
#' in the progress messages as needed to reach that threshold.
#'
#' \strong{Niche overlap.} For every pair of taxa with at least `min_occ`
#' records each, an occurrence density grid of `grid_size` x `grid_size` cells
#' is built for each taxon with [ecospat::ecospat.grid.clim.dyn()]; Schoener's
#' D is computed with [ecospat::ecospat.niche.overlap()] (D ranges from 0, no
#' overlap, to 1, identical niches); and the following are estimated with
#' `iterations` permutations:
#'
#' \itemize{
#'   \item the \emph{equivalency test}
#'   ([ecospat::ecospat.niche.equivalency.test()]), which asks whether the
#'   niches of the two taxa are identical;
#'   \item the \emph{similarity test}
#'   ([ecospat::ecospat.niche.similarity.test()], `rand.type = 2`), which asks
#'   whether the niches are more (or less) similar than expected by chance
#'   given the environmental space available.
#' }
#'
#' Niche dynamics indices are computed with
#' [ecospat::ecospat.niche.dyn.index()] (`intersection = 0.1`): `expansion`,
#' `stability` and `unfilling`, taken from the weighted index
#' (`dynamic.index.w`). Tests are significant when p < `alpha`. Fewer than 999
#' permutations triggers a warning, as they are generally discouraged for
#' publication. #' The random seed (`seed`) is set at the start of each pair for
#' reproducibility. Permutation results obtained in parallel
#' (`n_cores > 1`) are not guaranteed to be reproducible, because the seed
#' is not propagated to the workers; use `n_cores = 1` for exact
#' reproducibility.
#'
#' A pair that fails at any step produces a warning and is skipped, so one
#' problematic pair does not stop the whole run.
#'
#' \strong{Output.} When `save = TRUE`, files are written to a dated sub-folder
#' of `"Figs_NicheOverlap"` in the working directory (for example
#' `"Figs_NicheOverlap/18Sep2026"`), created if needed:
#'
#' \itemize{
#'   \item `EnvPCA_scree_combined`: cumulative variance
#'   explained by each PCA axis;
#'   \item `Environmental_PCA_loadings.xlsx` and
#'   `Environmental_PCA_variance.xlsx`: PC1-PC2 loadings, and variance
#'   explained by each axis;
#'   \item `Removed_correlated_variables.xlsx`: variables dropped for
#'   collinearity (only if any were dropped);
#'   \item `NicheTest_<taxon1>vs<taxon2>`: null
#'   distributions of the equivalency and similarity tests against the
#'   observed D, for each pair;
#'   \item `NicheOverlap_heatmap`: matrix of pairwise D;
#'   \item `NicheSpace_combined`: PC1-PC2 space with the
#'   density of all records and the occurrences of each taxon;
#'   \item `NicheOverlap_all_results.xlsx`: the full results table.
#' }
#'
#' Characters other than letters, digits, `.`, `_` and `-` in taxon names are
#' replaced by `_` in file names. Plots (scree, null distributions, heatmap and niche space) are saved in
#' the formats given by `file_formats`; tables are always saved as `.xlsx`.
#'
#' @param geodata Data frame returned by [geo_matrix_setting()], carrying the
#' `"env_cols"` attribute.
#'
#' @param var_threshold Numeric between 0 and 1. Cumulative variance threshold
#' displayed on the scree plot and used in the progress message on the number
#' of axes. Default is `0.95`.
#'
#' @param cor_cutoff Numeric between 0 and 1. Environmental variables with
#' pairwise absolute correlation above this value are removed. Default is
#' `0.8`.
#'
#' @param grid_size Integer. Resolution (number of cells per axis) of the
#' environmental density grid. Default is `100`.
#'
#' @param iterations Integer. Number of permutations for the equivalency and
#' similarity tests. Default is `1000`.
#'
#' @param alpha Numeric. Significance threshold for the tests. Default is
#' `0.05`.
#'
#' @param min_occ Integer. Minimum number of records per taxon required for a
#' pair to be analysed. Default is `3`.
#'
#' @param n_cores Integer. Number of cores used by \pkg{ecospat} for the
#' permutation tests. Default is `1`.
#'
#' @param seed Integer or `NULL`. Random seed set before each pair's
#' permutation tests. Default is `123`. If `NULL`, the random state is left
#' untouched.
#'
#' @param save Logical. If `TRUE` (default), figures and tables are written to
#' disk. If `FALSE`, nothing is saved and results are only returned.
#'
#' @param file_formats Character vector specifying which image formats to
#' save. One or more of `"pdf"` and `"jpg"`. Default is
#' `c("pdf", "jpg")`. Ignored when `save = FALSE`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' including the variables removed, the PCA summary and, for each pair, the
#' overlap statistics.
#'
#' @return
#' Invisibly, a named list with the elements:
#' \describe{
#'   \item{`results`}{A data frame with one row per analysed pair:
#'   `taxon1`, `taxon2`, `n1`, `n2` (records), `schoener_D`, `equiv_p`,
#'   `equiv_sig`, `equiv_null_mean`, `equiv_null_q025`, `equiv_null_q975`,
#'   `similarity_p`, `similarity_sig`, `sim_null_mean`, `n_permutations`,
#'   `expansion`, `stability`, `unfilling`, plus `PC1_variance`,
#'   `PC2_variance`, `n_variables` and `variables_used`.}
#'   \item{`D_matrix`}{Symmetric matrix of Schoener's D between taxa (`NA` on
#'   the diagonal and for pairs not analysed).}
#'   \item{`pca`}{The [stats::prcomp()] object of the combined PCA.}
#'   \item{`variables`}{Character vector of the environmental variables used
#'   in the PCA.}
#'   \item{`removed_variables`}{Character vector of the variables removed for
#'   zero variance or collinearity.}
#'   \item{`plots`}{A named list of [ggplot2::ggplot()] objects: `scree`,
#'   `heatmap`, `niche_space`, and `tests`, a list with one null-distribution
#'   plot per pair.}
#'   \item{`output_dir`}{Path of the folder where files were saved, or `NULL`
#'   if `save = FALSE`.}
#' }
#'
#' @seealso
#' [geo_matrix_setting()], [geo_enm()], [geo_bioclim_maps()]
#'
#' @examples
#' \dontrun{
#' geo <- geo_matrix_setting(
#'   xlsx_path = "occurrence_data.xlsx",
#'   taxon_col = "species",
#'   edaphic = TRUE
#' )
#'
#' # Default run: climate + edaphic variables
#' nic <- geo_niche_overlap(geo)
#'
#' # Custom settings
#' nic <- geo_niche_overlap(
#'   geo,
#'   var_threshold = 0.95,
#'   cor_cutoff = 0.8,
#'   grid_size = 100,
#'   iterations = 1000,
#'   alpha = 0.05,
#'   min_occ = 3,
#'   n_cores = 2
#' )
#'
#' # Results table and matrix of Schoener's D
#' nic$results
#' nic$D_matrix
#'
#' # Quick exploratory run without saving files
#' nic <- geo_niche_overlap(geo, iterations = 100, save = FALSE)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_vline
#' @importFrom ggplot2 geom_histogram geom_density geom_tile geom_text
#' @importFrom ggplot2 stat_density_2d after_stat facet_wrap labs theme
#' @importFrom ggplot2 theme_bw element_blank element_text scale_x_continuous
#' @importFrom ggplot2 scale_fill_gradient scale_fill_manual scale_fill_viridis_c
#' @importFrom ggnewscale new_scale_fill
#' @importFrom cowplot save_plot
#' @importFrom ecospat ecospat.grid.clim.dyn ecospat.niche.overlap
#' @importFrom ecospat ecospat.niche.equivalency.test ecospat.niche.similarity.test
#' @importFrom ecospat ecospat.niche.dyn.index
#' @importFrom caret findCorrelation
#' @importFrom openxlsx write.xlsx
#' @importFrom viridis viridis
#' @importFrom rlang .data
#' @importFrom stats complete.cases var cor prcomp na.omit quantile setNames
#' @importFrom utils combn
#'
#' @export

geo_niche_overlap <- function(geodata,
                              var_threshold = 0.95,
                              cor_cutoff = 0.8,
                              grid_size = 100,
                              iterations = 1000,
                              alpha = 0.05,
                              min_occ = 3,
                              n_cores = 1,
                              seed = 123,
                              save = TRUE,
                              file_formats = c("pdf", "jpg"),
                              verbose = TRUE) {

# ---- Validate input ---------------------------------------------------------

if (!is.data.frame(geodata)) {
  stop("`geodata` must be the data frame returned by `geo_matrix_setting()`.")
}

if (!"taxon" %in% names(geodata)) {
  stop("`geodata` must contain a `taxon` column.")
}

env_cols <- attr(geodata, "env_cols")

if (is.null(env_cols) || is.null(env_cols$clim)) {
  stop(
    "`geodata` has no valid \"env_cols\" attribute. ",
    "Build it with `geo_matrix_setting()`."
  )
}

if (!is.numeric(var_threshold) || length(var_threshold) != 1L ||
    var_threshold <= 0 || var_threshold > 1) {
  stop("`var_threshold` must be a single number in (0, 1].")
}

if (!is.numeric(cor_cutoff) || length(cor_cutoff) != 1L ||
    cor_cutoff <= 0 || cor_cutoff > 1) {
  stop("`cor_cutoff` must be a single number in (0, 1].")
}

if (!is.numeric(grid_size) || length(grid_size) != 1L || grid_size < 10) {
  stop("`grid_size` must be a single integer of at least 10.")
}

if (!is.numeric(iterations) || length(iterations) != 1L || iterations < 1) {
  stop("`iterations` must be a single positive integer.")
}

if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
  stop("`alpha` must be a single number between 0 and 1.")
}

if (!is.numeric(min_occ) || length(min_occ) != 1L || min_occ < 2) {
  stop("`min_occ` must be a single integer of at least 2.")
}

if (!is.numeric(n_cores) || length(n_cores) != 1L || n_cores < 1) {
  stop("`n_cores` must be a single positive integer.")
}

if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
  stop("`seed` must be NULL or a single number.")
}

if (!is.logical(save) || length(save) != 1L) {
  stop("`save` must be TRUE or FALSE.")
}

if (!is.character(file_formats) || length(file_formats) < 1L) {
  stop("`file_formats` must be a character vector with at least one format.")
}

file_formats <- match.arg(
  file_formats,
  choices = c("pdf", "jpg"),
  several.ok = TRUE
)

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

if (iterations < 999) {
  warning(
    "Permutation tests with fewer than 999 iterations are generally ",
    "discouraged for publication."
  )
}

# ---- Internal settings (not user arguments) -----------------------------------

base_dir <- "Figs_NicheOverlap" # parent folder for saved outputs
min_complete <- 20 # minimum complete records for the PCA

if (verbose) message("Running Niche Overlap Analysis...")

output_dir <- NULL

if (save) {
  output_dir <- file.path(base_dir, format(Sys.time(), "%d%b%Y"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
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

  if ("jpg" %in% file_formats) {
    cowplot::save_plot(
      file.path(output_dir, paste0(filename_stub, ".jpg")),
      plot_obj,
      ncol = 1,
      nrow = nrow,
      base_height = base_height,
      base_aspect_ratio = base_aspect_ratio,
      base_width = NULL
    )
  }
}

# ENVIRONMENTAL PCA SPACE (climate + edaphic, when available)

env_all_vars <- intersect(unique(c(env_cols$clim, env_cols$edaphic)), names(geodata))

if (length(env_all_vars) < 2) {
  stop("Fewer than 2 environmental variables found in geodata.")
}

env_source <- if (!is.null(env_cols$edaphic)) "climate + edaphic" else "climate"

env_data <- geodata[, env_all_vars, drop = FALSE]
env_data[] <- lapply(env_data, as.numeric)
complete <- stats::complete.cases(env_data)

if (sum(complete) < min_complete) {
  stop(
    "Fewer than ", min_complete,
    " complete environmental records available for the PCA."
  )
}

env_complete <- env_data[complete, , drop = FALSE]
taxon_vec <- geodata$taxon[complete]   # aligned with rows of env_complete

# ---- Drop zero-variance variables ----------------------------------------------
keep_var <- vapply(env_complete, function(x) stats::var(x, na.rm = TRUE) > 0, logical(1))
zero_var <- names(env_complete)[!keep_var]
env_complete <- env_complete[, keep_var, drop = FALSE]

if (ncol(env_complete) < 2) {
  stop("Fewer than 2 variables with variation remain in geodata.")
}

# ---- Drop highly correlated variables (|r| > cor_cutoff) --------------------------
cor_mat <- stats::cor(env_complete)
high_cor <- caret::findCorrelation(cor_mat, cutoff = cor_cutoff, names = TRUE)

removed_variables <- unique(c(zero_var, high_cor))

if (length(high_cor) > 0) {

  if (verbose) {
    message("  Dropping ", length(high_cor), " collinear variable(s): ",
            paste(high_cor, collapse = ", "))
  }

env_complete <- env_complete[, setdiff(names(env_complete), high_cor), drop = FALSE]

# Document exactly which predictors were excluded due to collinearity
  if (save) {
    openxlsx::write.xlsx(
      data.frame(Removed = high_cor),
      file.path(output_dir, "Removed_correlated_variables.xlsx"),
      overwrite = TRUE
    )
  }
}

if (ncol(env_complete) < 2) {
  stop("Fewer than 2 variables remain after collinearity filtering.")
}

# ---- Single combined PCA on the pooled records ---------------------------------------
pca_res <- stats::prcomp(env_complete, scale. = TRUE)

# Per-axis variance, kept separate from the cumulative vector
var_exp <- summary(pca_res)$importance[2, ]
cum_var <- cumsum(var_exp)
n_axes <- max(2, min(which(cum_var >= var_threshold)))

if (verbose) {
  message("  Environmental source: ", env_source,
          " (", ncol(env_complete), " variables)")
  message(sprintf(
    "  Combined PCA: %d axes needed for %.1f%% cumulative variance (%.1f%% reached).",
    n_axes, var_threshold * 100, cum_var[n_axes] * 100
  ))
  message(sprintf(
    "  PC1 = %.1f%% | PC2 = %.1f%% | PC1 + PC2 = %.1f%%",
    var_exp[1] * 100, var_exp[2] * 100, (var_exp[1] + var_exp[2]) * 100
  ))
}

# PCA loadings and variance-explained tables
loadings <- data.frame(
  Variable = rownames(pca_res$rotation),
  PC1 = pca_res$rotation[, 1],
  PC2 = pca_res$rotation[, 2]
)

eig <- data.frame(
  PC = paste0("PC", seq_along(var_exp)),
  Variance = round(var_exp, 4),
  Cumulative = round(cum_var, 4)
)

if (save) {
  openxlsx::write.xlsx(
    loadings, file.path(output_dir, "Environmental_PCA_loadings.xlsx"),
    overwrite = TRUE
  )
  openxlsx::write.xlsx(
    eig, file.path(output_dir, "Environmental_PCA_variance.xlsx"),
    overwrite = TRUE
  )
}

# Scree plot
scree_df <- data.frame(axis = seq_along(cum_var), cum_var = cum_var)

p_scree <- ggplot2::ggplot(scree_df, ggplot2::aes(x = .data$axis, y = .data$cum_var)) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_hline(yintercept = var_threshold, linetype = "dashed", color = "red") +
  ggplot2::geom_vline(xintercept = 2, linetype = "dashed", color = "blue") +
  ggplot2::scale_x_continuous(breaks = seq_along(cum_var)) +
  ggplot2::labs(
    x = "PCA axis", y = "Cumulative variance explained",
    title = paste0("Environmental PCA (", env_source, ") \u2014 cumulative variance"),
    subtitle = "Blue line: axes used for niche overlap (PC1+PC2)"
  ) +
  ggplot2::theme_bw()

if (save) {
  save_outputs(
    "EnvPCA_scree_combined", p_scree,
    base_height = 5, base_aspect_ratio = 1.4
  )
}

# Scores in PC1-PC2 space: the environmental space used by ecospat below
glob_scores <- as.data.frame(pca_res$x[, 1:2, drop = FALSE])
names(glob_scores) <- c("PC1", "PC2")

# RESOLVE TAXA

run_taxa <- sort(unique(stats::na.omit(taxon_vec)))
run_taxa <- run_taxa[nzchar(trimws(run_taxa))]

if (length(run_taxa) < 2) {
  stop("Need at least 2 taxa with complete environmental data to compare niche overlap.")
}

taxon_pairs <- utils::combn(run_taxa, 2, simplify = FALSE)

if (verbose) {
  message("  Taxa: ", paste(run_taxa, collapse = ", "))
  message("  Pairwise comparisons: ", length(taxon_pairs))
}

# NICHE OVERLAP

results_list <- list()
test_plots <- list()

for (pair in taxon_pairs) {

t1 <- pair[1]
t2 <- pair[2]

if (verbose) message("\n  --- ", t1, " vs ", t2, " ---")

if (!is.null(seed)) set.seed(seed)

occ1 <- glob_scores[!is.na(taxon_vec) & taxon_vec == t1, , drop = FALSE]
occ2 <- glob_scores[!is.na(taxon_vec) & taxon_vec == t2, , drop = FALSE]

if (nrow(occ1) < min_occ || nrow(occ2) < min_occ) {
  warning("Skipping ", t1, " vs ", t2, ": fewer than ", min_occ,
          " occurrences (n1 = ", nrow(occ1), ", n2 = ", nrow(occ2), ").",
          call. = FALSE)
  next
}

z1 <- .try_warn(
  ecospat::ecospat.grid.clim.dyn(
    glob = glob_scores, glob1 = glob_scores, sp = occ1, R = grid_size
  ),
  paste0("Grid construction failed for ", t1)
)

z2 <- .try_warn(
  ecospat::ecospat.grid.clim.dyn(
    glob = glob_scores, glob1 = glob_scores, sp = occ2, R = grid_size
  ),
  paste0("Grid construction failed for ", t2)
)

if (is.null(z1) || is.null(z2)) next

ov <- ecospat::ecospat.niche.overlap(z1, z2, cor = TRUE)
D <- unname(ov$D)

eq_test <- .try_warn(
  ecospat::ecospat.niche.equivalency.test(
    z1, z2, rep = iterations, ncores = n_cores
  ),
  "Equivalency test failed"
)

sim_test <- .try_warn(
  ecospat::ecospat.niche.similarity.test(
    z1, z2, rep = iterations, rand.type = 2, ncores = n_cores
  ),
  "Similarity test failed"
)

if (is.null(eq_test) || is.null(sim_test)) next

# Niche dynamics indices (expansion/stability/unfilling), Guisan et al. 2014
dyn <- tryCatch(
  ecospat::ecospat.niche.dyn.index(z1, z2, intersection = 0.1),
  error = function(e) NULL
)

dyn_val <- function(name) {
  if (is.null(dyn)) NA_real_ else round(unname(dyn$dynamic.index.w[name]), 4)
}

if (verbose) {
  message("  D             : ", round(D, 4))
  message("  Equivalency p : ", round(eq_test$p.D, 4))
  message("  Similarity p  : ", round(sim_test$p.D, 4))
  if (!is.null(dyn)) {
    message("  Expansion     : ", dyn_val("expansion"))
    message("  Stability     : ", dyn_val("stability"))
    message("  Unfilling     : ", dyn_val("unfilling"))
  }
}

# ---- Null distribution plot --------------------------------------------------
eq_sim <- eq_test$sim[, "D"]
sim_sim <- sim_test$sim[, "D"]

null_df <- data.frame(
  D = c(eq_sim, sim_sim),
  test = rep(c("Equivalency", "Similarity"), c(length(eq_sim), length(sim_sim)))
)
obs_df <- data.frame(D = c(D, D), test = c("Equivalency", "Similarity"))

p_null <- ggplot2::ggplot(null_df, ggplot2::aes(x = .data$D)) +
ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(.data$density)),
                        bins = 30, fill = "gray80", color = "white") +
ggplot2::geom_density(fill = "gray70", alpha = 0.4) +
ggplot2::geom_vline(data = obs_df, ggplot2::aes(xintercept = .data$D),
                    color = "red", linewidth = 0.8) +
ggplot2::facet_wrap(~ test, scales = "free_y") +
ggplot2::labs(
  x = "Schoener's D", y = "Density",
  title = paste0(t1, " vs ", t2, " \u2014 combined environmental space"),
  subtitle = paste0("Observed D = ", round(D, 4),
                    " | Eq. p = ", round(eq_test$p.D, 4),
                    " | Sim. p = ", round(sim_test$p.D, 4))
) +
ggplot2::theme_bw() +
ggplot2::theme(panel.grid = ggplot2::element_blank())

pair_name <- paste0(t1, "vs", t2)
test_plots[[pair_name]] <- p_null

if (save) {
  file_stub <- gsub("[^[:alnum:]._-]+", "_", pair_name)
  save_outputs(
    paste0("NicheTest_", file_stub), p_null,
    base_height = 4.5, base_aspect_ratio = 1.8
  )
  }

results_list[[pair_name]] <- data.frame(
  taxon1 = t1,
  taxon2 = t2,
  n1 = nrow(occ1),
  n2 = nrow(occ2),
  schoener_D = round(D, 4),
  equiv_p = round(eq_test$p.D, 4),
  equiv_sig = eq_test$p.D < alpha,
  equiv_null_mean = round(mean(eq_sim), 4),
  equiv_null_q025 = round(unname(stats::quantile(eq_sim, 0.025)), 4),
  equiv_null_q975 = round(unname(stats::quantile(eq_sim, 0.975)), 4),
  similarity_p = round(sim_test$p.D, 4),
  similarity_sig = sim_test$p.D < alpha,
  sim_null_mean = round(mean(sim_sim), 4),
  n_permutations = iterations,
  expansion = dyn_val("expansion"),
  stability = dyn_val("stability"),
  unfilling = dyn_val("unfilling")
)
}

if (length(results_list) == 0) {
  stop("No taxon pairs produced valid niche overlap results.")
}

results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL

# PCA/variable metadata added to the results table
results_df$PC1_variance <- round(unname(var_exp[1]), 4)
results_df$PC2_variance <- round(unname(var_exp[2]), 4)
results_df$n_variables <- ncol(env_complete)
results_df$variables_used <- paste(colnames(env_complete), collapse = ", ")

# PAIRWISE D HEATMAP

mat <- matrix(NA_real_, length(run_taxa), length(run_taxa),
              dimnames = list(run_taxa, run_taxa))

for (i in seq_len(nrow(results_df))) {
  mat[results_df$taxon1[i], results_df$taxon2[i]] <- results_df$schoener_D[i]
  mat[results_df$taxon2[i], results_df$taxon1[i]] <- results_df$schoener_D[i]
}

heat_df <- as.data.frame(as.table(mat))
names(heat_df) <- c("taxon1", "taxon2", "D")

p_heat <- ggplot2::ggplot(
  heat_df, ggplot2::aes(x = .data$taxon1, y = .data$taxon2, fill = .data$D)
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(is.na(.data$D), "", sprintf("%.3f", .data$D))),
    size = 3
  ) +
  ggplot2::scale_fill_viridis_c(option = "mako", limits = c(0, 1), na.value = "white") +
  ggplot2::labs(x = NULL, y = NULL, title = "Schoener's D \u2014 combined environmental space") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "italic"),
    axis.text.y = ggplot2::element_text(face = "italic"),
    panel.grid = ggplot2::element_blank()
  )

if (save) {
  save_outputs(
    "NicheOverlap_heatmap", p_heat,
    base_height = 5, base_aspect_ratio = 1
  )
  }

# PCA DENSITY PLOT: pooled space + occurrences by taxon

in_run <- !is.na(taxon_vec) & taxon_vec %in% run_taxa
occ_df <- cbind(glob_scores[in_run, , drop = FALSE], taxon = taxon_vec[in_run])

color_vector <- stats::setNames(viridis::viridis(length(run_taxa)), run_taxa)

p_pca <- ggplot2::ggplot() +
  ggplot2::stat_density_2d(
    data = glob_scores,
    ggplot2::aes(x = .data$PC1, y = .data$PC2, fill = ggplot2::after_stat(.data$level)),
    geom = "polygon", alpha = 0.3, show.legend = FALSE
  ) +
  ggplot2::scale_fill_gradient(low = "white", high = "gray40") +
  ggnewscale::new_scale_fill() +
  ggplot2::geom_point(
    data = occ_df,
    ggplot2::aes(x = .data$PC1, y = .data$PC2, fill = .data$taxon),
    shape = 21, size = 2, colour = "white", alpha = 0.8
  ) +
  ggplot2::scale_fill_manual(values = color_vector) +
  ggplot2::labs(
    x = paste0("PC1 (", round(var_exp[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2] * 100, 1), "%)"),
    fill = "taxon",
    title = paste0("Combined environmental space (", env_source, ")")
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(face = "italic")
  )

if (save) {
  save_outputs(
    "NicheSpace_combined", p_pca,
    base_height = 6, base_aspect_ratio = 1.3
  )
  }

# SAVE RESULTS

if (save) {
  openxlsx::write.xlsx(
    results_df, file.path(output_dir, "NicheOverlap_all_results.xlsx"),
    overwrite = TRUE
  )
}

if (verbose) {
  if (save) {
    message("\nDone. Outputs saved in: ", output_dir)
  } else {
    message("\nDone. Outputs were not saved (`save = FALSE`).")
  }
}

invisible(list(
  results = results_df,
  D_matrix = mat,
  pca = pca_res,
  variables = colnames(env_complete),
  removed_variables = removed_variables,
  plots = list(scree = p_scree, heatmap = p_heat,
               niche_space = p_pca, tests = test_plots),
  output_dir = output_dir
))
}

#' Evaluate an expression, turning errors into warnings
#'
#' @description
#' Internal helper. Evaluates `expr`; if it fails, a warning with `fail_msg`
#' and the original error message is issued and `NULL` is returned.
#'
#' @param expr Expression to evaluate.
#' @param fail_msg Character. Message prefix used in the warning.
#'
#' @return The value of `expr`, or `NULL` on error.
#'
#' @keywords internal
#'
#' @noRd
.try_warn <- function(expr, fail_msg) {
  tryCatch(expr, error = function(e) {
    warning(fail_msg, ": ", conditionMessage(e), call. = FALSE)
    NULL
  })
}
