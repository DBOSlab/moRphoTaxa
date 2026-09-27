#' Redundancy analysis (RDA) of morphological traits constrained by environment
#'
#' @description
#' Tests how much of the morphological variation among specimens is
#' explained by environmental variables (by default, the WorldClim
#' bioclimatic variables) using a redundancy analysis (RDA), run separately
#' for the vegetative, floral, and fruit trait blocks and for their
#' combinations.
#'
#' Takes the morphometric matrix produced by [morph_matrix_setting()] and
#' the geographical/environmental matrix produced by
#' [geo_matrix_setting()], matches specimens by their identifiers, and, for
#' each trait block, prepares the data, fits the RDA with [vegan::rda()],
#' tests the model with permutations, draws a biplot and a chart of the
#' contribution of each environmental variable, and writes the results to
#' Excel workbooks.
#'
#' @details
#' \strong{Matching specimens.}
#' Specimens are matched by row names. Both [geo_matrix_setting()] and
#' [morph_matrix_setting()] use the specimen identifier (an existing
#' \code{"id"} column, or automatic identifiers \code{"ID0001"},
#' \code{"ID0002"}, ...) as row names, so the two matrices must be built
#' from the same specimens (typically the same spreadsheet, or spreadsheets
#' sharing an \code{"id"} column).
#'
#' \strong{Trait blocks.}
#' The vegetative (\code{"veg"}), floral (\code{"flo"}), and fruit
#' (\code{"fru"}) trait blocks are read from the \code{"base_cols"}
#' attribute of \code{analysis_data}, created by [morph_matrix_setting()].
#' With \code{blocks = "all"}, each available block is run, followed by
#' every combination of two or more available blocks, named by
#' concatenating the block names (\code{"vegflo"}, \code{"vegfru"},
#' \code{"flofru"}, \code{"vegflofru"}). Any single block or combination
#' can also be requested by name.
#'
#' \strong{Environmental variables.}
#' If \code{env_vars = NULL} (default), all bioclimatic columns
#' (\code{bio_1}-\code{bio_19}) recorded in the \code{"env_cols"} attribute
#' of \code{geodata} are used; \code{"elevation"} is not included by
#' default. Any other numeric column of \code{geodata}, such as
#' \code{"elevation"} or a soil variable, can be used by naming it in
#' \code{env_vars}. Specimens with a missing value in any of the selected
#' variables are dropped.
#'
#' \strong{Data preparation (repeated for each block).}
#' \enumerate{
#'   \item Specimens present in both matrices are kept, and specimens with
#'   a missing environmental value or with no measured trait in the block
#'   are dropped.
#'   \item Traits with no observation among the retained specimens are
#'   dropped.
#'   \item Remaining missing trait values are replaced by the mean of the
#'   trait among the retained specimens (mean imputation), because RDA
#'   requires complete data. Blocks with many missing values are therefore
#'   less reliable, and the imputed values shrink the variance of the
#'   affected traits.
#'   \item Traits and environmental variables with zero variance are
#'   dropped.
#'   \item Environmental variables are screened for collinearity: a
#'   variable is dropped when its absolute correlation with any variable
#'   listed after it exceeds \code{env_cor_threshold}.
#'   \item If \code{scale_data = TRUE}, traits and environmental variables
#'   are standardized (mean 0, variance 1).
#' }
#' A block is skipped, with a warning, when fewer than 10 specimens, fewer
#' than 2 traits, or fewer than 2 environmental variables remain (two
#' constrained axes are needed for the biplot), or when the number of
#' environmental variables is not smaller than the number of specimens
#' minus one, in which case the model has no residual degrees of freedom.
#' Reduce \code{env_vars} or use a larger sample in that case.
#'
#' \strong{Model and tests.}
#' The model is \code{rda(traits ~ ., data = environment)}. The overall
#' model is tested with [vegan::anova.cca()], and each environmental
#' variable is tested with \code{by = "terms"}. These are \emph{sequential}
#' tests: each variable is evaluated after those listed before it, so its
#' F statistic and p-value depend on the order of the variables. Both use
#' \code{permutations} permutations; because these are random, call
#' [set.seed()] before this function to obtain reproducible p-values. The
#' percentage of constrained variance is unadjusted and increases with the
#' number of environmental variables relative to the number of specimens;
#' the adjusted R-squared ([vegan::RsquareAdj()]) is reported alongside it.
#' If the tests for individual variables fail, a warning is issued and only
#' the overall test and the biplot are produced for that block.
#'
#' \strong{Biplot.}
#' Specimens are shown in constrained ordination space, coloured by taxon,
#' with the environmental variables as arrows. Taxon colours are taken from
#' the \code{"taxon_colors"} attribute of \code{analysis_data}, so they
#' match those of the other \pkg{moRphoTaxa} plots; if the attribute is
#' missing, colours are drawn from the \code{viridis} palette. Axis labels
#' give the percentage of the total variance explained by each axis. Trait
#' scores are not plotted but are saved in the Excel workbook.
#'
#' \strong{Outputs.}
#' The following files are written to a dated subfolder
#' (\code{format(Sys.time(), "\%d\%b\%Y")}) inside \code{"Figs_RDA"}. The
#' dated subfolder means repeated runs do not overwrite each other.
#' \itemize{
#'   \item \code{RDA_biplot_<block>.<ext>}: biplot;
#'   \item \code{RDA_terms_<block>.<ext>}: F statistic of each
#'   environmental variable, coloured by significance (p < 0.05);
#'   \item \code{RDA_<block>.xlsx}: one workbook per block, with sheets
#'   \code{summary}, \code{by_term}, \code{env_scores}, and
#'   \code{sp_scores};
#'   \item \code{RDA_summary.xlsx}: one row per block.
#' }
#'
#' @param analysis_data A data frame produced by [morph_matrix_setting()],
#' with specimen identifiers as row names, a \code{"taxon"} column, and the
#' \code{"base_cols"} attribute. Do not strip attributes (e.g. by subsetting
#' or rebuilding the data frame) before calling this function.
#'
#' @param geodata A data frame produced by [geo_matrix_setting()], with
#' specimen identifiers as row names. When \code{env_vars = NULL}, it must
#' still carry its \code{"env_cols"} attribute.
#'
#' @param blocks Character vector of blocks to analyze. \code{"all"}
#' (default) runs every available block and every combination of blocks.
#' Otherwise, one or more of \code{"veg"}, \code{"flo"}, \code{"fru"},
#' \code{"vegflo"}, \code{"vegfru"}, \code{"flofru"}, and
#' \code{"vegflofru"}, restricted to the blocks available in
#' \code{analysis_data}. Unavailable names are ignored with a warning.
#'
#' @param env_vars Character vector of column names of \code{geodata} to use
#' as environmental variables. If \code{NULL} (default), all bioclimatic
#' variables (\code{bio_1}-\code{bio_19}) are used.
#'
#' @param permutations Integer. Number of permutations used in the
#' significance tests. Default is \code{999}.
#'
#' @param scale_data Logical. If \code{TRUE} (default), traits and
#' environmental variables are standardized before the analysis.
#'
#' @param env_cor_threshold Numeric between \code{0} and \code{1}. Absolute
#' correlation above which an environmental variable is dropped as
#' collinear. Default is \code{0.95}.
#'
#' @param file_formats Character vector specifying the file formats used to
#' save the figures. One or more of \code{"pdf"} and \code{"jpg"}. Default
#' is \code{c("pdf", "jpg")}.
#'
#' @param verbose Logical. If \code{TRUE} (default), prints progress
#' messages and, for each block, the number of specimens, traits, and
#' environmental variables, the collinear variables dropped, and the
#' variance explained.
#'
#' @return
#' A list, returned invisibly, with the elements:
#' \describe{
#'   \item{\code{results}}{A data frame with one row per analyzed block and
#'   the columns \code{block}, \code{n} (specimens), \code{n_traits},
#'   \code{n_env}, \code{constrained_var} (percentage of constrained
#'   variance), \code{adj_r2} (adjusted R-squared), and \code{p_value}
#'   (overall permutation test), or \code{NULL} if no block was analyzed.}
#'   \item{\code{models}}{A named list of the [vegan::rda()] models, one per
#'   analyzed block.}
#'   \item{\code{plots}}{A named list, one element per analyzed block, each
#'   a list with the [ggplot2::ggplot()] objects \code{biplot} and
#'   \code{terms} (\code{NULL} if the tests for individual variables
#'   failed).}
#'   \item{\code{output_dir}}{Path of the folder where the outputs were
#'   saved.}
#' }
#' The function is also called for its side effects: figures and Excel
#' workbooks are written to \code{output_dir}.
#'
#' @seealso
#' [morph_matrix_setting()], [geo_matrix_setting()],
#' [geo_bioclim_exploratory()], [geo_morph_mantel()]
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
#' # Every block and block combination, all bioclimatic variables
#' set.seed(123)
#' res <- geo_morph_rda(analysis_data = morph, geodata = geo)
#' res$results
#'
#' # Floral traits only, with a chosen set of environmental variables
#' res <- geo_morph_rda(
#'   analysis_data = morph,
#'   geodata = geo,
#'   blocks = "flo",
#'   env_vars = c("bio_1", "bio_4", "bio_12", "bio_15", "elevation")
#' )
#'
#' # Stricter collinearity filter, no standardization, PDF only
#' res <- geo_morph_rda(
#'   analysis_data = morph,
#'   geodata = geo,
#'   env_cor_threshold = 0.8,
#'   scale_data = FALSE,
#'   file_formats = "pdf"
#' )
#' }
#'
#' @importFrom stats complete.cases cor setNames var
#' @importFrom utils combn
#' @importFrom vegan rda anova.cca scores RsquareAdj
#' @importFrom openxlsx write.xlsx
#' @importFrom viridis viridis
#' @importFrom rlang .data
#' @importFrom grid arrow unit
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_text geom_col
#' @importFrom ggplot2 geom_vline geom_hline scale_fill_manual coord_flip labs
#' @importFrom ggplot2 theme_bw theme element_blank element_text
#' @importFrom cowplot save_plot
#'
#' @export

geo_morph_rda <- function(analysis_data,
                          geodata,
                          blocks = "all",
                          env_vars = NULL,
                          permutations = 999,
                          scale_data = TRUE,
                          env_cor_threshold = 0.95,
                          file_formats = c("pdf", "jpg"),
                          verbose = TRUE) {

# ---- Validate input -------------------------------------------------------

if (!is.data.frame(analysis_data)) {
  stop("`analysis_data` must be a data frame, typically the output of morph_matrix_setting().")
}

if (!is.data.frame(geodata)) {
  stop("`geodata` must be a data frame, typically the output of geo_matrix_setting().")
}

if (!"taxon" %in% names(analysis_data)) {
  stop("`analysis_data` must contain a \"taxon\" column.")
}

base_cols <- attr(analysis_data, "base_cols")

if (is.null(base_cols) || length(base_cols) == 0L) {
  stop(
    "`analysis_data` has no \"base_cols\" attribute. ",
    "Pass the data frame returned by morph_matrix_setting() unchanged."
  )
}

if (!is.character(blocks) || length(blocks) < 1L) {
  stop("`blocks` must be a character vector with at least one block name.")
}

if (!is.null(env_vars) && (!is.character(env_vars) || length(env_vars) < 1L)) {
  stop("`env_vars` must be NULL or a character vector with at least one column name.")
}

if (!is.numeric(permutations) || length(permutations) != 1L ||
    is.na(permutations) || permutations < 1 ||
    permutations != round(permutations)) {
  stop("`permutations` must be a single positive integer.")
}

if (!is.logical(scale_data) || length(scale_data) != 1L || is.na(scale_data)) {
  stop("`scale_data` must be TRUE or FALSE.")
}

if (!is.numeric(env_cor_threshold) || length(env_cor_threshold) != 1L ||
    is.na(env_cor_threshold) || env_cor_threshold < 0 || env_cor_threshold > 1) {
  stop("`env_cor_threshold` must be a single number between 0 and 1.")
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

# ---- Resolve blocks -------------------------------------------------------

# Available single blocks, in canonical order
avail <- intersect(c("veg", "flo", "fru"), names(base_cols))

if (length(avail) == 0L) {
  stop("`base_cols` contains none of the blocks \"veg\", \"flo\", or \"fru\".")
}

all_blocks <- base_cols[avail]

# Every combination of two or more available blocks (vegflo, vegfru, ...)
if (length(avail) >= 2L) {
  for (k in 2:length(avail)) {
    for (cmb in utils::combn(avail, k, simplify = FALSE)) {
      all_blocks[[paste(cmb, collapse = "")]] <-
        unlist(base_cols[cmb], use.names = FALSE)
    }
  }
}

if ("all" %in% blocks) {
  run_blocks <- all_blocks
} else {
  requested <- intersect(blocks, names(all_blocks))
  invalid   <- setdiff(blocks, names(all_blocks))

  if (length(requested) == 0L) {
    stop(
      "No valid blocks in `blocks`. Use \"all\" or any of: ",
      paste(names(all_blocks), collapse = ", ")
    )
  }
  if (length(invalid) > 0L) {
    warning("Ignoring unavailable blocks: ", paste(invalid, collapse = ", "))
  }
  run_blocks <- all_blocks[requested]
}

# ---- Resolve environmental variables --------------------------------------

if (is.null(env_vars)) {
  env_cols <- attr(geodata, "env_cols")

  if (is.null(env_cols) || is.null(env_cols$clim)) {
    stop(
      "`geodata` has no \"env_cols\" attribute with a \"clim\" element. ",
      "Pass the data frame returned by geo_matrix_setting() unchanged, ",
      "or set `env_vars` manually."
    )
  }

  env_vars <- setdiff(env_cols$clim, "elevation")

  if (length(env_vars) == 0L) {
    stop("No bioclimatic columns found in `geodata`. Set `env_vars` manually.")
  }

  if (verbose) {
    message(
      "  Auto-selected ", length(env_vars), " bioclimatic variables: ",
      paste(env_vars, collapse = ", ")
    )
  }
}

missing_env <- setdiff(env_vars, names(geodata))

if (length(missing_env) > 0L) {
  stop(
    "Environmental variable(s) not found in `geodata`: ",
    paste(missing_env, collapse = ", ")
  )
}

# ---- Match specimens ------------------------------------------------------

common_ids <- intersect(rownames(analysis_data), rownames(geodata))

if (length(common_ids) == 0L) {
  stop(
    "No specimen identifiers are shared by `analysis_data` and `geodata`. ",
    "Both must be built from the same specimens (row names are the specimen IDs)."
  )
}

taxon_vec <- stats::setNames(
  as.character(analysis_data$taxon),
  rownames(analysis_data)
)
taxon_colors <- attr(analysis_data, "taxon_colors")

# ---- Output folder ----------------------------------------------------------

output_dir <- file.path("Figs_RDA", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (verbose) {
  message("Running Redundancy Analysis (RDA)...")
  message("  Blocks to run: ", paste(names(run_blocks), collapse = ", "))
}

# ---- Run each block -----------------------------------------------------------

results_list <- list()
models <- list()
plots <- list()

for (block_name in names(run_blocks)) {

if (verbose) {
  message("\n  --- Block: ", block_name, " ---")
}

valid_cols <- intersect(run_blocks[[block_name]], names(analysis_data))
valid_cols <- valid_cols[vapply(analysis_data[valid_cols], is.numeric, logical(1))]

if (length(valid_cols) == 0L) {
  warning("Skipping '", block_name, "': no valid columns.")
  next
}

block_data <- analysis_data[common_ids, valid_cols, drop = FALSE]
env_data <- geodata[common_ids, env_vars, drop = FALSE]
env_data[] <- lapply(env_data, as.numeric)

# Drop specimens with missing env or no morphological observation
keep <- stats::complete.cases(env_data) & rowSums(!is.na(block_data)) > 0
block_data <- block_data[keep, , drop = FALSE]
env_data <- env_data[keep, , drop = FALSE]

# Drop traits with no observation
block_data <- block_data[, colSums(!is.na(block_data)) > 0, drop = FALSE]

# Impute column means for the remaining morphological NAs
block_data[] <- lapply(block_data, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

# Drop zero-variance traits and environmental variables
nonconst_trait <- vapply(block_data, function(x) {
  v <- stats::var(x, na.rm = TRUE)
  !is.na(v) && v > 0
}, logical(1))
block_data <- block_data[, nonconst_trait, drop = FALSE]

nonconst_env <- vapply(env_data, function(x) {
  v <- stats::var(x, na.rm = TRUE)
  !is.na(v) && v > 0
}, logical(1))
env_data <- env_data[, nonconst_env, drop = FALSE]

# Drop collinear environmental variables
if (ncol(env_data) > 1L) {
  env_cor <- stats::cor(env_data)
  env_cor[upper.tri(env_cor, diag = TRUE)] <- 0
  drop_env <- which(apply(abs(env_cor) > env_cor_threshold, 2, any))

  if (length(drop_env) > 0L) {
    if (verbose) {
      message(
        "  Dropping collinear env vars: ",
        paste(colnames(env_data)[drop_env], collapse = ", ")
      )
    }
    env_data <- env_data[, -drop_env, drop = FALSE]
  }
}

if (nrow(block_data) < 10L || ncol(block_data) < 2L || ncol(env_data) < 2L) {
  warning("Skipping '", block_name, "': not enough data.")
  next
}

if (ncol(env_data) >= nrow(block_data) - 1L) {
  warning(
    "Skipping '", block_name, "': the number of environmental variables (",
    ncol(env_data), ") must be smaller than the number of specimens minus one (",
    nrow(block_data) - 1L, "). Reduce `env_vars`."
  )
  next
}

if (verbose) {
  message(
    "  Specimens: ", nrow(block_data),
    " | Traits: ", ncol(block_data),
    " | Env vars: ", ncol(env_data)
  )
}

# Standardize if requested
if (scale_data) {
  block_mat <- scale(block_data)
  env_mat <- scale(env_data)
} else {
  block_mat <- as.matrix(block_data)
  env_mat <- as.matrix(env_data)
}
env_df <- as.data.frame(env_mat)

# ---- Fit RDA ----------------------------------------------------------------

rda_model <- tryCatch(
  vegan::rda(block_mat ~ ., data = env_df),
  error = function(e) {
    warning("RDA failed for '", block_name, "': ", conditionMessage(e))
    NULL
  }
)
if (is.null(rda_model)) next

if (length(rda_model$CCA$eig) < 2L) {
  warning("Skipping '", block_name, "': fewer than two constrained axes.")
  next
}

# ---- Permutation tests ----------------------------------------------------------

rda_perm <- vegan::anova.cca(rda_model, permutations = permutations)

rda_by_term <- tryCatch(
  vegan::anova.cca(rda_model, by = "terms", permutations = permutations),
  error = function(e) {
    warning(
      "Tests for individual variables failed for '", block_name, "': ",
      conditionMessage(e)
    )
    NULL
  }
)

p_overall <- rda_perm[1, "Pr(>F)"]

# ---- Variance explained -----------------------------------------------------------

eig_constr <- rda_model$CCA$eig
total_var <- sum(eig_constr) + sum(rda_model$CA$eig)
constrained <- sum(eig_constr) / total_var * 100
unconstrained <- sum(rda_model$CA$eig) / total_var * 100
axis_pct <- eig_constr / total_var * 100
adj_r2 <- vegan::RsquareAdj(rda_model)$adj.r.squared

if (verbose) {
  message("  Constrained var    : ", round(constrained, 2), "%")
  message("  Unconstrained var  : ", round(unconstrained, 2), "%")
  message("  Adjusted R-squared : ", round(adj_r2, 3))
  message("  Overall p-value    : ", format(p_overall, digits = 3))
}

# ---- Biplot -------------------------------------------------------------------------

site_scores <- as.data.frame(vegan::scores(rda_model, display = "sites"))
sp_scores <- as.data.frame(vegan::scores(rda_model, display = "species"))
env_scores  <- as.data.frame(vegan::scores(rda_model, display = "bp"))

site_scores$taxon <- unname(taxon_vec[rownames(site_scores)])
taxa_present <- sort(unique(site_scores$taxon))

if (!is.null(taxon_colors) && all(taxa_present %in% names(taxon_colors))) {
  color_vector <- taxon_colors[taxa_present]
} else {
  color_vector <- stats::setNames(
    viridis::viridis(length(taxa_present)),
    taxa_present
  )
}

sp_scores$trait <- rownames(sp_scores)
sp_scores$label <- gsub("_", " ", gsub("/.*", "", sp_scores$trait))
env_scores$variable <- rownames(env_scores)
env_scores$label <- .short_env_label(env_scores$variable)
rownames(sp_scores)  <- NULL
rownames(env_scores) <- NULL

p_rda <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = site_scores,
    ggplot2::aes(x = .data$RDA1, y = .data$RDA2, fill = .data$taxon),
    shape = 21, size = 2.5, colour = "white", alpha = 0.8
  ) +
  ggplot2::geom_segment(
    data = env_scores,
    ggplot2::aes(x = 0, y = 0, xend = .data$RDA1, yend = .data$RDA2),
    arrow = grid::arrow(length = grid::unit(0.2, "cm")),
    color = "red", linewidth = 0.6
  ) +
  ggplot2::geom_text(
    data = env_scores,
    ggplot2::aes(x = .data$RDA1 * 1.1, y = .data$RDA2 * 1.1, label = .data$label),
    color = "red", size = 3, fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = color_vector) +
  ggplot2::geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0, color = "gray40", linetype = "dashed") +
  ggplot2::labs(
    x = paste0("RDA1 (", round(axis_pct[1], 1), "%)"),
    y = paste0("RDA2 (", round(axis_pct[2], 1), "%)"),
    fill  = "",
    title = paste0(
      "RDA - ", block_name,
      " | Constrained: ", round(constrained, 1), "%",
      " | p = ", format(p_overall, digits = 3)
    )
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid  = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(face = "italic")
  )

for (ext in file_formats) {
  cowplot::save_plot(
    file.path(output_dir, paste0("RDA_biplot_", block_name, ".", ext)),
    p_rda,
    ncol = 1, nrow = 1, base_height = 7, base_aspect_ratio = 1.3
  )
}

# ---- Per-variable contribution plot -----------------------------------------------------

p_terms <- NULL
by_term_out <- NULL

if (!is.null(rda_by_term)) {
  by_term_out <- data.frame(
    variable = rownames(rda_by_term),
    as.data.frame(rda_by_term),
    check.names = FALSE,
    row.names = NULL
  )

  term_df <- by_term_out[!is.na(by_term_out$F), , drop = FALSE]
  term_df$variable <- .short_env_label(term_df$variable)
  term_df <- term_df[order(term_df$F, decreasing = TRUE), , drop = FALSE]
  term_df$variable <- factor(term_df$variable, levels = rev(term_df$variable))
  term_df$sig <- !is.na(term_df[["Pr(>F)"]]) & term_df[["Pr(>F)"]] < 0.05

  p_terms <- ggplot2::ggplot(
    term_df,
    ggplot2::aes(x = .data$variable, y = .data$F, fill = .data$sig)
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#1a9850", "FALSE" = "gray70"),
      labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s.")
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL, y = "F statistic", fill = "",
      title = paste0("RDA - Variable importance: ", block_name)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 8)
    )

  for (ext in file_formats) {
    cowplot::save_plot(
      file.path(output_dir, paste0("RDA_terms_", block_name, ".", ext)),
      p_terms,
      ncol = 1, nrow = 1,
      base_height = max(4, nrow(term_df) * 0.35),
      base_aspect_ratio = 1.5
    )
  }
}

# ---- Save block workbook ------------------------------------------------------------------

sheets <- list(
  summary = data.frame(
    block = block_name,
    n = nrow(block_data),
    n_traits = ncol(block_data),
    n_env = ncol(env_data),
    constrained_var = round(constrained, 4),
    unconstrained_var = round(unconstrained, 4),
    adj_r2 = round(adj_r2, 4),
    p_value = p_overall
  )
)
if (!is.null(by_term_out)) sheets$by_term <- by_term_out
sheets$env_scores <- env_scores
sheets$sp_scores <- sp_scores

openxlsx::write.xlsx(
  sheets,
  file.path(output_dir, paste0("RDA_", block_name, ".xlsx")),
  overwrite = TRUE
)

models[[block_name]] <- rda_model
plots[[block_name]] <- list(biplot = p_rda, terms = p_terms)

results_list[[block_name]] <- data.frame(
  block = block_name,
  n = nrow(block_data),
  n_traits = ncol(block_data),
  n_env = ncol(env_data),
  constrained_var = round(constrained, 2),
  adj_r2 = round(adj_r2, 3),
  p_value = p_overall
)
}

# ---- Save summary ---------------------------------------------------------------------------

if (length(results_list) > 0L) {
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL

  openxlsx::write.xlsx(
    results,
    file.path(output_dir, "RDA_summary.xlsx"),
    overwrite = TRUE
  )
} else {
  results <- NULL
  warning("No block was analyzed; no summary table was written.")
}

if (verbose) {
  message("\nDone. Outputs saved in: ", output_dir)
}

invisible(list(
  results = results,
  models = models,
  plots = plots,
  output_dir = output_dir
))
}

#' Short label for an environmental variable
#'
#' @description
#' Internal helper used by [geo_morph_rda()] to shorten WorldClim variable
#' names for plots (`"bio_12"` becomes `"bio12"`). Other names are returned
#' unchanged.
#'
#' @param x Character vector of variable names.
#'
#' @return A character vector of the same length as `x`.
#'
#' @keywords internal
#'
#' @noRd
.short_env_label <- function(x) {
  gsub("^bio_", "bio", x)
}
