#' PERMANOVA and PERMDISP on morphological trait blocks
#'
#' @description
#' Tests whether taxa differ in multivariate location (PERMANOVA,
#' [vegan::adonis2()]) and in multivariate dispersion (PERMDISP,
#' [vegan::betadisper()]) on morphological trait blocks derived from a
#' morphometric matrix, and generates a standard set of diagnostic figures
#' and spreadsheets for each trait block analyzed. The function is designed
#' to work directly on the objects produced by [morph_matrix_setting()], and
#' mirrors the structure of [morph_pcoa()] and [morph_dapc()].
#'
#' For each requested trait block, the function computes a specimen-by-specimen
#' distance matrix, tests for a difference in multivariate location among taxa
#' with PERMANOVA, tests for a difference in multivariate dispersion (spread)
#' among taxa with PERMDISP, runs pairwise PERMANOVA and pairwise PERMDISP
#' between every pair of taxa, saves a boxplot of each specimen's distance to
#' its taxon centroid, saves a PCoA ordination of the distance matrix with
#' taxon centroids and confidence ellipses, and saves an Excel workbook with
#' every result table.
#'
#' @details
#' `morph_permanova()` expects `analysis_data` as produced by
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
#' The `permanova_blocks` argument controls which trait sets are analyzed:
#'
#' \itemize{
#'   \item `"all"`: runs the analysis on each individual block in `base_cols`
#'   (typically `"veg"`, `"flo"`, `"fru"`), as well as on every pairwise and
#'   full combination of these blocks (`"vegflo"`, `"vegfru"`, `"flofru"`,
#'   `"vegflofru"`);
#'   \item a character vector naming one or more blocks in `base_cols`: runs
#'   the analysis only on the requested block(s), individually.
#' }
#'
#' For each trait block, rows with missing values in the relevant columns are
#' removed via [stats::na.omit()], and traits with zero variance are dropped.
#' The grouping variable is the `taxon` column of `analysis_data`; any taxon
#' with fewer than `min_n` specimens (after the above filtering) is removed
#' from the block entirely, with a message identifying the removed taxa.
#' Blocks left with fewer than 2 taxa, fewer than 2 traits, or fewer
#' specimens than taxa, are skipped with a warning.
#'
#' \strong{Distance.} The distance between specimens is set by `dist_method`,
#' passed to [vegan::vegdist()]: `"gower"` (default), `"euclidean"` or
#' `"bray"`. For `"euclidean"`, traits are centred and scaled first, so the
#' result matches the distance underlying [morph_pca()]. For `"gower"` and
#' `"bray"`, traits are used on their original scale. `"bray"` requires
#' non-negative trait values and stops with an informative error otherwise.
#'
#' \strong{PERMANOVA.} [vegan::adonis2()] partitions the distance matrix by
#' `taxon`, with `permutations` permutations, testing whether taxa differ in
#' multivariate location (centroid). The reported statistics are the pseudo-F
#' ratio, R\eqn{^2} (proportion of the total sum of squares explained by
#' taxon) and a permutation p-value.
#'
#' \strong{PERMDISP.} [vegan::betadisper()] (`type = "median"`) computes each
#' specimen's distance to its taxon's spatial median in the ordination of the
#' distance matrix, and [vegan::permutest.betadisper()] tests whether these
#' distances differ among taxa, with `permutations` permutations. PERMDISP is
#' the multivariate analogue of a test for homogeneity of variance: a
#' significant PERMANOVA together with a significant PERMDISP means that the
#' taxa may differ in dispersion as well as (or instead of) location, so the
#' PERMANOVA result alone should not be read as evidence of separate
#' centroids.
#'
#' \strong{Pairwise tests.} For every pair of taxa retained in a block,
#' PERMANOVA is repeated on the corresponding sub-matrix of the distance
#' matrix, and PERMDISP pairwise comparisons are obtained from
#' [vegan::permutest.betadisper()] (`pairwise = TRUE`). In both cases,
#' p-values are adjusted across all pairs within a block with
#' [stats::p.adjust()] (`method = p_adjust`).
#'
#' The random seed (`seed`) is set once per block, before the PERMANOVA,
#' PERMDISP, and all pairwise tests of that block, so permutation p-values do
#' not depend on the order in which blocks are analyzed.
#'
#' Taxon colours are taken from the `"taxon_colors"` attribute attached to
#' `analysis_data` by [morph_matrix_setting()], so that each taxon keeps the
#' same colour used in other plots (e.g. [morph_pca()], [morph_pcoa()],
#' [morph_dapc()]). If this attribute is missing, or lacks an entry for a
#' taxon retained in a given block, colours are generated automatically from
#' the `viridis` palette for the affected taxa, with a warning.
#'
#' For each analyzed block, the following figures are generated:
#'
#' \itemize{
#'   \item a boxplot (with individual points) of each specimen's distance to
#'   its taxon's centroid, one box per taxon;
#'   \item a PCoA ordination (axes 1 and 2 of the distance matrix used for
#'   PERMDISP) with a normal confidence ellipse per taxon at level
#'   `ellipse_level`.
#' }
#'
#' Output files are written to a date-specific directory inside
#' `Figs_PERMANOVA`. Each figure is saved in the formats specified by
#' `file_formats`. File names are prefixed by figure type (`"PERMDISP_boxplot_"`,
#' `"PERMANOVA_ordination_"`) and suffixed by the trait block name.
#'
#' @param analysis_data A numeric data frame containing morphological traits
#' as columns and specimens as rows, with `id` and `taxon` columns and
#' specimen identifiers as row names. Typically the `analysis_data` element
#' returned by [morph_matrix_setting()]. If this object carries a
#' `"taxon_colors"` attribute, those colours are reused; otherwise a
#' `viridis` palette is generated automatically.
#'
#' @param permanova_blocks Character. Either `"all"`, to run the analysis on
#' every individual block in `base_cols` plus all pairwise and full
#' combinations of these blocks, or a character vector naming one or more
#' blocks in `base_cols` to analyze individually. Default is `"all"`.
#'
#' @param dist_method Character. Distance between specimens, passed to
#' [vegan::vegdist()]: `"gower"` (default), `"euclidean"` (on scaled traits)
#' or `"bray"` (requires non-negative traits).
#'
#' @param permutations Integer. Number of permutations used by PERMANOVA and
#' PERMDISP (overall and pairwise). Default is `999`.
#'
#' @param min_n Integer. Minimum number of specimens a taxon must have,
#' within a given trait block, to be retained for that block. Taxa with fewer
#' specimens are dropped block-by-block. Default is `3`.
#'
#' @param p_adjust Character. Method used to adjust pairwise p-values across
#' all pairs within a block, passed to [stats::p.adjust()]. Default is
#' `"BH"` (Benjamini-Hochberg).
#'
#' @param ellipse_level Numeric. Confidence level (between 0 and 1) used for
#' the per-taxon normal confidence ellipses drawn on the PCoA ordination
#' plot. Default is `0.95`.
#'
#' @param seed Integer or `NULL`. Random seed set once per block, before its
#' PERMANOVA, PERMDISP and pairwise tests, so a block's permutation p-values
#' do not depend on the other blocks analyzed. Default is `123`. If `NULL`,
#' the random state is left untouched.
#'
#' @param save_xlsx Logical. If `TRUE` (default), an Excel workbook with the
#' PERMANOVA, PERMDISP and pairwise result tables is written for each
#' analyzed block, plus a `"PERMANOVA_summary.xlsx"` collecting one row per
#' block.
#'
#' @param file_formats Character vector specifying which image formats to
#' save. One or more of `"pdf"` and `"jpeg"`. Default is
#' `c("pdf", "jpeg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' including removed taxa, skipped blocks, test statistics and the output
#' directory, to the console.
#'
#' @return
#' Invisibly returns a data frame summarizing the blocks successfully
#' analyzed, with one row per block and columns `block`, `n`, `n_groups`,
#' `permanova_F`, `permanova_R2`, `permanova_p`, `permdisp_F` and
#' `permdisp_p`; or invisibly returns `NULL` if no block could be analyzed.
#'
#' The function is primarily called for its side effects: boxplots, PCoA
#' ordination plots, and per-block and summary Excel workbooks are written to
#' disk for each analyzed trait block.
#'
#' @seealso
#' [morph_matrix_setting()],
#' [morph_pcoa()],
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
#' # Run PERMANOVA + PERMDISP (Gower distance) on every block and combination
#' morph_permanova(
#'   analysis_data = res,
#'   permanova_blocks = "all"
#' )
#'
#' # Vegetative block only, Euclidean distance, 4999 permutations
#' morph_permanova(
#'   analysis_data = res,
#'   permanova_blocks = "veg",
#'   dist_method = "euclidean",
#'   permutations = 4999
#' )
#'
#' # Selected blocks, saving PDF output only
#' perm <- morph_permanova(
#'   analysis_data = res,
#'   permanova_blocks = c("veg", "flo"),
#'   min_n = 5,
#'   file_formats = "pdf"
#' )
#'
#' # Summary table
#' perm
#' }
#'
#' @importFrom stats na.omit setNames var p.adjust as.formula
#' @importFrom vegan vegdist adonis2 betadisper permutest
#' @importFrom ggplot2 aes element_blank element_text geom_boxplot geom_jitter
#'   geom_point geom_hline geom_vline ggplot labs scale_fill_manual
#'   scale_color_manual stat_ellipse theme theme_bw
#' @importFrom rlang .data
#' @importFrom cowplot save_plot
#' @importFrom viridis viridis
#' @importFrom openxlsx write.xlsx
#' @importFrom utils combn
#'
#' @export

morph_permanova <- function(analysis_data,
                            permanova_blocks = "all",
                            dist_method = c("gower", "euclidean", "bray"),
                            permutations = 999,
                            min_n = 3,
                            p_adjust = "BH",
                            ellipse_level = 0.95,
                            seed = 123,
                            save_xlsx = TRUE,
                            file_formats = c("pdf", "jpeg"),
                            verbose = TRUE) {

# ---- Validate input ---------------------------------------------------------

if (!is.data.frame(analysis_data)) {
  stop("`analysis_data` must be a data frame.")
}

if (!all(c("taxon") %in% names(analysis_data))) {
  stop("`analysis_data` must contain taxon` column.")
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

if (!is.character(permanova_blocks) || length(permanova_blocks) < 1) {
  stop("`permanova_blocks` must be \"all\" or a character vector of block names.")
}

dist_method <- match.arg(dist_method)

if (!is.numeric(permutations) || length(permutations) != 1L || permutations < 99) {
  stop("`permutations` must be a single integer of at least 99.")
}

if (!is.numeric(min_n) || length(min_n) != 1L || min_n < 2) {
  stop("`min_n` must be a single integer of at least 2.")
}

if (!is.character(p_adjust) || length(p_adjust) != 1L ||
    !(p_adjust %in% stats::p.adjust.methods)) {
  stop("`p_adjust` must be one of: ", paste(stats::p.adjust.methods, collapse = ", "))
}

if (!is.numeric(ellipse_level) || length(ellipse_level) != 1L ||
    ellipse_level <= 0 || ellipse_level >= 1) {
  stop("`ellipse_level` must be a single numeric value between 0 and 1.")
}

if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
  stop("`seed` must be NULL or a single number.")
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

# ---- Species vector ----------------------------------------------------------

species_all <- stats::setNames(
  as.character(analysis_data$taxon),
  rownames(analysis_data)
)

# ---- Output folder ------------------------------------------------------------

output_dir <- file.path("Figs_PERMANOVA", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (verbose) {
  message("Running PERMANOVA + PERMDISP analysis (", dist_method, " distance)")
}

# ---- Resolve blocks ----------------------------------------------------------

if (identical(permanova_blocks, "all")) {
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
  requested <- intersect(permanova_blocks, names(base_cols))
  if (length(requested) == 0) {
    stop(
      "No valid blocks in `permanova_blocks`. Use \"all\" or any of: ",
      paste(names(base_cols), collapse = ", ")
    )
  }
  run_blocks <- base_cols[requested]
}

if (verbose) {
  message("  Blocks to run: ", paste(names(run_blocks), collapse = ", "))
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

# ---- Helper: pairwise PERMANOVA on a distance matrix --------------------------

pairwise_permanova <- function(dmat, groups, perms) {

  dmat_full <- as.matrix(dmat)
  pairs <- utils::combn(levels(groups), 2, simplify = FALSE)

  out <- lapply(pairs, function(pair) {
    idx <- groups %in% pair
    sub_dist <- stats::as.dist(dmat_full[idx, idx, drop = FALSE])
    sub_grp  <- droplevels(groups[idx])

    res <- tryCatch(
      vegan::adonis2(sub_dist ~ sub_grp, permutations = perms),
      error = function(e) NULL
    )

    if (is.null(res)) {
      return(data.frame(
        taxon1 = pair[1], taxon2 = pair[2],
        F = NA_real_, R2 = NA_real_, p = NA_real_
      ))
    }

    data.frame(
      taxon1 = pair[1], taxon2 = pair[2],
      F  = res[["F"]][1],
      R2 = res[["R2"]][1],
      p  = res[["Pr(>F)"]][1]
    )
  })

  do.call(rbind, out)
}

# ---- Helper: pairwise PERMDISP table from permutest(pairwise = TRUE) ----------

pairwise_permdisp <- function(disp_test, groups) {

pw <- disp_test$pairwise$permuted

if (is.null(pw)) return(NULL)

pair_names <- do.call(rbind, strsplit(names(pw), "-", fixed = TRUE))

data.frame(
  taxon1 = pair_names[, 1],
  taxon2 = pair_names[, 2],
  p = unname(pw)
)
}

# ---- PERMANOVA + PERMDISP for each trait set -----------------------------------

results <- list()

for (set_name in names(run_blocks)) {

if (verbose) message("\n  --- Block: ", set_name, " ---")

valid_cols <- intersect(run_blocks[[set_name]], names(analysis_data))
d <- stats::na.omit(analysis_data[, valid_cols, drop = FALSE])

if (nrow(d) < 2) {
  warning(
    "Skipping '", set_name, "': fewer than 2 specimens with complete data ",
    "(n = ", nrow(d), ")."
  )
  next
}

# Drop traits with no variation
keep_trait <- vapply(d, function(x) {
  v <- stats::var(x)
  !is.na(v) && v > 0
}, logical(1))
d <- d[, keep_trait, drop = FALSE]

sp <- droplevels(factor(species_all[rownames(d)]))
# ---- Remove small groups ---------------------------------------------------
group_counts <- table(sp)
small_groups <- names(group_counts)[group_counts < min_n]

if (length(small_groups) > 0) {
  if (verbose) {
    message("  Removing small groups: ", paste(small_groups, collapse = ", "))
  }
  keep_rows <- !sp %in% small_groups
  d  <- d[keep_rows, , drop = FALSE]
  sp <- droplevels(sp[keep_rows])
}

n_groups <- length(levels(sp))

if (n_groups < 2 || ncol(d) < 2 || nrow(d) <= n_groups) {
  warning(
    "Skipping '", set_name, "': not enough data (n = ", nrow(d),
    ", p = ", ncol(d), ", groups = ", n_groups, ")."
  )
  next
}

if (verbose) {
  message(
    "  Specimens: ", nrow(d),
    " | Traits: ", ncol(d),
    " | Groups: ", n_groups
  )
}

# ---- Distance matrix --------------------------------------------------------

if (dist_method == "bray" && any(d < 0)) {
  warning(
    "Skipping '", set_name, "': Bray-Curtis distance requires non-negative ",
    "trait values. Use `dist_method = \"gower\"` or `\"euclidean\"` instead."
  )
  next
}

d_mat <- if (dist_method == "euclidean") scale(as.matrix(d)) else as.matrix(d)
dmat  <- vegan::vegdist(d_mat, method = dist_method)

if (!is.null(seed)) set.seed(seed)

# ---- PERMANOVA ---------------------------------------------------------------

perm_res <- tryCatch(
  vegan::adonis2(dmat ~ sp, permutations = permutations),
  error = function(e) {
    warning("PERMANOVA failed for '", set_name, "': ", conditionMessage(e))
    NULL
  }
)

if (is.null(perm_res)) next

permanova_F  <- perm_res[["F"]][1]
permanova_R2 <- perm_res[["R2"]][1]
permanova_p  <- perm_res[["Pr(>F)"]][1]

# ---- PERMDISP ------------------------------------------------------------------

betadisp <- tryCatch(
  vegan::betadisper(dmat, sp, type = "median"),
  error = function(e) {
    warning("PERMDISP failed for '", set_name, "': ", conditionMessage(e))
    NULL
  }
)

if (is.null(betadisp)) next

disp_test <- tryCatch(
  vegan::permutest(betadisp, permutations = permutations, pairwise = n_groups > 2),
  error = function(e) {
    warning("PERMDISP permutation test failed for '", set_name, "': ", conditionMessage(e))
    NULL
  }
)

if (is.null(disp_test)) next

permdisp_F <- disp_test$tab[["F"]][1]
permdisp_p <- disp_test$tab[["Pr(>F)"]][1]

if (verbose) {
  message(sprintf("  PERMANOVA : F = %.3f | R2 = %.3f | p = %.4f",
                  permanova_F, permanova_R2, permanova_p))
  message(sprintf("  PERMDISP  : F = %.3f | p = %.4f",
                  permdisp_F, permdisp_p))
}

# ---- Pairwise tests --------------------------------------------------------------

pw_permanova <- pairwise_permanova(dmat, sp, permutations)
pw_permanova$p_adj <- stats::p.adjust(pw_permanova$p, method = p_adjust)

pw_permdisp <- if (n_groups > 2) pairwise_permdisp(disp_test, sp) else NULL
if (!is.null(pw_permdisp)) {
  pw_permdisp$p_adj <- stats::p.adjust(pw_permdisp$p, method = p_adjust)
}

# ---- Taxon colours -----------------------------------------------------------

color_vector <- attr(analysis_data, "taxon_colors")

if (is.null(color_vector)) {

  color_vector <- stats::setNames(
    viridis::viridis(n_groups),
    levels(sp)
  )

} else {

  missing_sp <- setdiff(levels(sp), names(color_vector))

  if (length(missing_sp) > 0) {
    warning(
      "No stored colour found for: ", paste(missing_sp, collapse = ", "),
      ". Assigning automatic colours for these taxa."
    )
    color_vector[missing_sp] <- viridis::viridis(length(missing_sp))
  }

  color_vector <- color_vector[levels(sp)]
}

# ---- Boxplot of distances to centroid ---------------------------------------------

dist_df <- data.frame(
  id = names(betadisp$distances),
  taxon = sp,
  distance = betadisp$distances
)

p_box <- ggplot2::ggplot(
  dist_df,
  ggplot2::aes(x = .data[["taxon"]], y = .data[["distance"]], fill = .data[["taxon"]])
) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.15, size = 1.5, colour = "gray30", alpha = 0.6) +
  ggplot2::scale_fill_manual(values = color_vector) +
  ggplot2::labs(
    x = NULL, y = "Distance to taxon centroid",
    title = paste0("PERMDISP \u2014 ", set_name),
    subtitle = sprintf("F = %.3f | p = %.4f", permdisp_F, permdisp_p)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "italic", angle = 45, hjust = 1),
    panel.grid = ggplot2::element_blank()
  )

save_outputs(paste0("PERMDISP_boxplot_", set_name), p_box,
             base_height = 6, base_aspect_ratio = 1.2)

# ---- PCoA ordination with confidence ellipses -------------------------------------

ord_scores <- as.data.frame(betadisp$vectors[, 1:2, drop = FALSE])
names(ord_scores) <- c("Axis1", "Axis2")
ord_scores$taxon <- sp

eig <- betadisp$eig
rel_eig <- eig / sum(abs(eig))
lab_x <- sprintf("PCoA 1 (%.1f%%)", 100 * rel_eig[1])
lab_y <- sprintf("PCoA 2 (%.1f%%)", 100 * rel_eig[2])

p_ord <- ggplot2::ggplot(
  ord_scores,
  ggplot2::aes(x = .data[["Axis1"]], y = .data[["Axis2"]])
) +
  ggplot2::geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::stat_ellipse(
    ggplot2::aes(color = .data[["taxon"]]),
    type = "norm", level = ellipse_level,
    linetype = "dashed", linewidth = 0.8, show.legend = FALSE
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = .data[["taxon"]]), shape = 21, size = 3,
    colour = "white", alpha = 0.8
  ) +
  ggplot2::scale_fill_manual(values = color_vector) +
  ggplot2::scale_color_manual(values = color_vector) +
  ggplot2::labs(
    x = lab_x, y = lab_y, fill = "taxon",
    title = paste0("PERMANOVA \u2014 ", set_name),
    subtitle = sprintf("F = %.3f | R\u00b2 = %.3f | p = %.4f",
                       permanova_F, permanova_R2, permanova_p)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(face = "italic")
  )

save_outputs(paste0("PERMANOVA_ordination_", set_name), p_ord,
             base_height = 7, base_aspect_ratio = 1.3)

# ---- Save Excel ------------------------------------------------------------------

if (save_xlsx) {
  openxlsx::write.xlsx(
    list(
      permanova = data.frame(
        block = set_name, F = permanova_F, R2 = permanova_R2, p = permanova_p
      ),
      permdisp = data.frame(
        block = set_name, F = permdisp_F, p = permdisp_p
      ),
      pairwise_permanova = pw_permanova,
      pairwise_permdisp = if (!is.null(pw_permdisp)) pw_permdisp else data.frame(),
      distances_to_centroid = dist_df
    ),
    file.path(output_dir, paste0("PERMANOVA_", set_name, ".xlsx"))
  )
}

if (verbose) message("  Done: ", set_name)

results[[set_name]] <- data.frame(
  block = set_name,
  n = nrow(d),
  n_groups = n_groups,
  permanova_F = round(permanova_F, 4),
  permanova_R2 = round(permanova_R2, 4),
  permanova_p = round(permanova_p, 4),
  permdisp_F = round(permdisp_F, 4),
  permdisp_p = round(permdisp_p, 4)
)
}

# ---- Save summary ------------------------------------------------------------

if (length(results) > 0) {
  summary_df <- do.call(rbind, results)
  rownames(summary_df) <- NULL
  if (save_xlsx) {
    openxlsx::write.xlsx(summary_df, file.path(output_dir, "PERMANOVA_summary.xlsx"))
  }
} else {
  summary_df <- NULL
}

if (verbose) {
  message("\nPERMANOVA + PERMDISP analysis done. Outputs saved in: ", output_dir)
}

invisible(summary_df)
}
