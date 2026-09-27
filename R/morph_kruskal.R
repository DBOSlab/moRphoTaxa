#' Kruskal-Wallis test of morphological traits across a grouping variable
#'
#' @description
#' Runs a Kruskal-Wallis test for every morphological trait in each requested
#' trait block, testing whether trait distributions differ across the levels
#' of one or more grouping variables, and generates a summary bar plot and
#' spreadsheet for each block x group combination. When a trait is
#' significant, Dunn's post-hoc test can be run automatically for every pair
#' of groups, with a pairwise p-value heatmap per trait. The function is
#' designed to work directly on the objects produced by
#' [morph_matrix_setting()].
#'
#' @details
#' `morph_kruskal()` expects `analysis_data` as produced by
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
#' Grouping variables are read from columns of `analysis_data` itself, named
#' in `kw_groups` — by default just `"taxon"`, the grouping column that
#' [morph_matrix_setting()] always attaches. Passing any other name only
#' works if that column has been added to `analysis_data` beforehand;
#' `morph_kruskal()` does not read from any other object. Every name in
#' `kw_groups` is run, not only the first.
#'
#' The `kw_blocks` argument controls which trait sets are analyzed:
#'
#' \itemize{
#'   \item `"all"`: runs the test on each individual block in `base_cols`
#'   (typically `"veg"`, `"flo"`, `"fru"`), as well as on every pairwise and
#'   full combination of these blocks (`"vegflo"`, `"vegfru"`, `"flofru"`,
#'   `"vegflofru"`). Combinations whose source blocks are absent from
#'   `base_cols` are silently skipped;
#'   \item a character vector naming one or more blocks in `base_cols`: runs
#'   the test only on the requested block(s), individually.
#' }
#'
#' Within each block x group combination, specimens with a missing value for
#' the grouping variable are removed; the combination is skipped with a
#' warning if fewer than 3 specimens or fewer than 2 group levels remain.
#' Each trait is then tested individually with [stats::kruskal.test()]: rows
#' with a missing value for that trait are dropped for that test only, and a
#' trait is skipped (silently, since this is expected for traits scored only
#' in part of a block) if it has fewer than 3 non-missing values, or if the
#' test fails to fit (for instance because it is constant within every
#' group).
#'
#' P-values are adjusted with [stats::p.adjust()] using `kw_p_adjust`,
#' separately within each block x group combination — that is, the
#' correction accounts for the number of traits tested in that block, not for
#' the total number of tests across all blocks and groups in the call. A
#' trait is flagged `significant` when its adjusted p-value is below
#' `kw_alpha`.
#'
#' When `kw_posthoc = TRUE` (default), every trait flagged significant in a
#' given block x group combination is followed up with Dunn's test
#' ([dunn.test::dunn.test()]) for all pairwise group comparisons, using the
#' one-sided-equivalent p-values (`altp = TRUE`) and the same
#' `kw_p_adjust` method. Because `dunn.test` uses different method names
#' than [stats::p.adjust()], `kw_p_adjust` is translated internally
#' (`"BH"` to `"bh"`, `"BY"` to `"by"`, and so on); `"hommel"` has no
#' equivalent in `dunn.test` and cannot be combined with `kw_posthoc = TRUE`.
#' A pairwise comparison whose Dunn test fails to fit — typically because one
#' of the two groups has too few specimens with a non-missing value for that
#' trait — is skipped with a warning rather than aborting the block. For each
#' trait with post-hoc results, a group x group heatmap of adjusted p-values
#' is saved, with asterisks marking conventional significance bands
#' (`*** p < 0.001`, `** p < 0.01`, `* p < 0.05`) independently of
#' `kw_alpha`.
#'
#' Output files are written to a date-specific directory inside
#' `Figs.KruskalWallis`. For each block x group combination, a bar plot of H
#' statistics (traits sorted from highest to lowest, coloured by
#' significance) is saved in the formats specified by `file_formats`, named
#' `"KW_<block>_<group>"`. When `save_xlsx = TRUE`, a spreadsheet with that
#' combination's results is written alongside it, and a single
#' `"KW_all_results.xlsx"` combining every block x group combination is
#' written once all combinations have been processed. Post-hoc output for a
#' given block x group combination is written to a `"posthoc_<block>_<group>"`
#' subdirectory: one heatmap per significant trait, a spreadsheet of that
#' combination's pairwise comparisons, and, when `save_xlsx = TRUE`, a single
#' `"KW_posthoc_all_results.xlsx"` combining post-hoc results from every
#' combination is written to the top-level output directory once all
#' combinations have been processed.
#'
#' @param analysis_data A numeric data frame containing morphological traits
#' as columns and specimens as rows, with `id` and `taxon` columns and
#' specimen identifiers as row names. Typically the object returned by
#' [morph_matrix_setting()].
#'
#' @param kw_blocks Character. Either `"all"`, to test every individual block
#' in `base_cols` plus all pairwise and full combinations of these blocks, or
#' a character vector naming one or more blocks in `base_cols` to test
#' individually. Default is `"all"`.
#'
#' @param kw_groups Character vector of column names in `analysis_data` to
#' use as grouping variables, one test run per name. Default is `"taxon"`.
#'
#' @param kw_p_adjust Character. Method passed to [stats::p.adjust()] for
#' correcting p-values across the traits tested within each block x group
#' combination, and (translated) to [dunn.test::dunn.test()] for the
#' post-hoc pairwise comparisons. One of `"bonferroni"`, `"holm"`,
#' `"hochberg"`, `"BH"` (default), `"BY"`, or `"none"`; `"hommel"` is
#' accepted only when `kw_posthoc = FALSE`, since `dunn.test` has no
#' equivalent method.
#'
#' @param kw_alpha Numeric between 0 and 1. Adjusted-p-value threshold below
#' which a trait (and, in the post-hoc step, a pairwise comparison) is
#' flagged as significant. Default is `0.05`.
#'
#' @param kw_posthoc Logical. If `TRUE` (default), runs Dunn's post-hoc test
#' and saves a pairwise heatmap for every trait flagged significant.
#'
#' @param save_xlsx Logical. If `TRUE` (default), writes spreadsheets of
#' results — per block x group combination, per post-hoc combination, and
#' combined across the whole run — to the output directory.
#'
#' @param file_formats Character vector specifying which image formats to
#' save. One or more of `"pdf"` and `"jpeg"`. Default is
#' `c("pdf", "jpeg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' per-combination counts of traits tested and found significant, skipped
#' combinations, and the output directory to the console.
#'
#' @return
#' Invisibly returns a list with two elements, or `NULL` if no block x group
#' combination could be analyzed:
#' \describe{
#'   \item{`kw`}{A data frame combining the Kruskal-Wallis results of every
#'   analyzed block x group combination, with columns `trait`, `H_stat`,
#'   `df`, `p_value`, `p_adjusted`, `significant`, `block`, and `group`.}
#'   \item{`posthoc`}{A data frame combining the Dunn post-hoc results of
#'   every significant trait across every analyzed combination, with columns
#'   `comparison`, `Z`, `p_value`, `p_adjusted`, `significant`, `trait`,
#'   `block`, and `group`; `NULL` if `kw_posthoc = FALSE` or no trait was
#'   significant.}
#' }
#'
#' The function is also called for its side effects: bar plots, post-hoc
#' heatmaps, and (when `save_xlsx = TRUE`) spreadsheets of results are
#' written to disk.
#'
#' @seealso
#' [morph_matrix_setting()],
#' [morph_pca()],
#' [morph_kmeans()],
#' [morph_nmm()],
#' [morph_anova()]
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
#' # Test every individual block and every block combination against taxon,
#' # with automatic Dunn post-hoc for significant traits
#' kw <- morph_kruskal(
#'   analysis_data = res,
#'   kw_blocks = "all"
#' )
#'
#' # Traits that differ significantly among taxa
#' subset(kw$kw, significant)
#'
#' # Pairwise comparisons behind those traits
#' kw$posthoc
#'
#' # Vegetative block only, without post-hoc testing
#' morph_kruskal(
#'   analysis_data = res,
#'   kw_blocks = "veg",
#'   kw_posthoc = FALSE
#' )
#' }
#'
#' @importFrom stats kruskal.test p.adjust as.formula
#' @importFrom ggplot2 aes element_blank element_text geom_col geom_text
#'   geom_tile ggplot coord_flip labs scale_fill_gradient2 scale_fill_manual
#'   theme theme_bw
#' @importFrom cowplot save_plot
#' @importFrom openxlsx write.xlsx
#' @importFrom dunn.test dunn.test
#' @importFrom rlang .data
#'
#' @export

morph_kruskal <- function(analysis_data,
                          kw_blocks = "all",
                          kw_groups = "taxon",
                          kw_p_adjust = "BH",
                          kw_alpha = 0.05,
                          kw_posthoc = TRUE,
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

if (!is.character(kw_blocks) || length(kw_blocks) < 1) {
  stop("`kw_blocks` must be \"all\" or a character vector of block names.")
}

if (!is.character(kw_groups) || length(kw_groups) < 1) {
  stop("`kw_groups` must be a character vector of column names in `analysis_data`.")
}

missing_groups <- setdiff(kw_groups, names(analysis_data))

if (length(missing_groups) > 0) {
  stop(
    "Group column(s) not found in `analysis_data`: ",
    paste(missing_groups, collapse = ", ")
  )
}

kw_p_adjust <- match.arg(kw_p_adjust, choices = stats::p.adjust.methods)

if (!is.logical(kw_posthoc) || length(kw_posthoc) != 1L) {
  stop("`kw_posthoc` must be TRUE or FALSE.")
}

# `dunn.test` has no "hommel" adjustment; fail fast rather than partway
# through the run if post-hoc testing is requested with it.
dunn_method_map <- c(
  holm = "holm",
  hochberg = "hochberg",
  bonferroni = "bonferroni",
  BH = "bh",
  BY = "by",
  none = "none"
)

if (kw_posthoc && !(kw_p_adjust %in% names(dunn_method_map))) {
  stop(
    "`kw_p_adjust = \"", kw_p_adjust, "\"` has no equivalent method in ",
    "`dunn.test`, so it cannot be used with `kw_posthoc = TRUE`. Choose ",
    "one of ", paste(shQuote(names(dunn_method_map)), collapse = ", "),
    ", or set `kw_posthoc = FALSE`."
  )
}

if (!is.numeric(kw_alpha) ||
    length(kw_alpha) != 1L ||
    kw_alpha <= 0 ||
    kw_alpha >= 1) {
  stop("`kw_alpha` must be a single number between 0 and 1.")
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

# ---- Output folder ------------------------------------------------------------

output_dir <- file.path("Figs.KruskalWallis", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (verbose) {
  message("Running Kruskal-Wallis test")
}

# ---- Resolve blocks ----------------------------------------------------------

if (identical(kw_blocks, "all")) {
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
  requested <- intersect(kw_blocks, names(base_cols))
  if (length(requested) == 0) {
    stop(
      "No valid blocks in `kw_blocks`. Use \"all\" or any of: ",
      paste(names(base_cols), collapse = ", ")
    )
  }
  run_blocks <- base_cols[requested]
}

if (verbose) {
  message("  Blocks to run: ", paste(names(run_blocks), collapse = ", "))
}

# ---- Helper: trait display label ---------------------------------------------

trait_label <- function(x) {
  gsub("_", " ", gsub("/.*", "", x))
}

# ---- Helper: save a plot in the requested formats ----------------------------

save_outputs <- function(dir, filename_stub, plot_obj, base_height, base_aspect_ratio) {
  if ("pdf" %in% file_formats) {
    cowplot::save_plot(
      file.path(dir, paste0(filename_stub, ".pdf")),
      plot_obj,
      ncol = 1,
      nrow = 1,
      base_height = base_height,
      base_aspect_ratio = base_aspect_ratio,
      base_width = NULL
    )
  }

  if ("jpeg" %in% file_formats) {
    cowplot::save_plot(
      file.path(dir, paste0(filename_stub, ".jpeg")),
      plot_obj,
      ncol = 1,
      nrow = 1,
      base_height = base_height,
      base_aspect_ratio = base_aspect_ratio,
      base_width = NULL
    )
  }
}

# ---- Helper: Dunn post-hoc for one trait -------------------------------------

run_dunn_trait <- function(trait, block_data, dunn_dir, block_name, group_var) {

  y <- block_data[[trait]]
  g <- block_data$group__

  complete <- !is.na(y)
  y <- y[complete]
  g <- droplevels(g[complete])

  dunn_res <- tryCatch(
    dunn.test::dunn.test(
      y, g,
      method = dunn_method_map[[kw_p_adjust]],
      altp = TRUE,
      kw = FALSE,
      table = FALSE
    ),
    error = function(e) {
      warning(
        "Dunn post-hoc failed for trait '", trait, "' in '", block_name,
        "' x '", group_var, "': ", conditionMessage(e)
      )
      NULL
    }
  )

if (is.null(dunn_res)) return(NULL)

dunn_df <- data.frame(
  comparison = dunn_res$comparisons,
  Z = round(dunn_res$Z, 4),
  p_value = round(dunn_res$altP, 4),
  p_adjusted = round(dunn_res$altP.adjusted, 4),
  significant = dunn_res$altP.adjusted < kw_alpha,
  stringsAsFactors = FALSE
)
dunn_df$trait <- trait

# ---- Pairwise heatmap -------------------------------------------------------

groups_present <- levels(g)
mat <- matrix(
  NA_real_,
  length(groups_present),
  length(groups_present),
  dimnames = list(groups_present, groups_present)
)

for (i in seq_len(nrow(dunn_df))) {
  parts <- strsplit(dunn_df$comparison[i], " - ")[[1]]
  if (length(parts) == 2 &&
      parts[1] %in% groups_present &&
      parts[2] %in% groups_present) {
    mat[parts[1], parts[2]] <- dunn_df$p_adjusted[i]
    mat[parts[2], parts[1]] <- dunn_df$p_adjusted[i]
  }
}
diag(mat) <- 1

heat_df <- as.data.frame(as.table(mat))
names(heat_df) <- c("group1", "group2", "p_adj")
heat_df$p_adj <- as.numeric(heat_df$p_adj)
heat_df$sig <- ifelse(
  is.na(heat_df$p_adj), "",
  ifelse(heat_df$p_adj < 0.001, "***",
         ifelse(heat_df$p_adj < 0.01, "**",
                ifelse(heat_df$p_adj < 0.05, "*", "n.s.")
         )
  )
)

p_heat <- ggplot2::ggplot(
  heat_df,
  ggplot2::aes(
    x = .data[["group1"]],
    y = .data[["group2"]],
    fill = .data[["p_adj"]]
  )
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(ggplot2::aes(label = .data[["sig"]]), size = 4) +
  ggplot2::scale_fill_gradient2(
    low = "#d73027", mid = "white", high = "gray90",
    midpoint = 0.05, limits = c(0, 1),
    name = "adj. p", na.value = "gray95"
  ) +
  ggplot2::labs(
    x = NULL, y = NULL,
    title = paste0("Dunn post-hoc - ", trait_label(trait))
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "italic"),
    axis.text.y = ggplot2::element_text(face = "italic"),
    panel.grid = ggplot2::element_blank()
  )

tname <- gsub("[^a-zA-Z0-9]", "_", trait)

save_outputs(
  dunn_dir,
  paste0("Dunn_", tname),
  p_heat,
  base_height = 5,
  base_aspect_ratio = 1
)

dunn_df
}

# ---- Helper: KW (+ optional Dunn) for one block x group combination ----------

run_kw_block <- function(block_name, block_columns, group_var) {

valid_cols <- intersect(block_columns, names(analysis_data))
block_data <- analysis_data[, valid_cols, drop = FALSE]

# Grouping variable, read directly from analysis_data so it stays aligned
# with block_data by specimen ID (row name) regardless of row order
block_data$group__ <- as.factor(analysis_data[rownames(block_data), group_var])

# Drop specimens with a missing group
block_data <- block_data[!is.na(block_data$group__), , drop = FALSE]

if (nrow(block_data) < 3 || length(unique(block_data$group__)) < 2) {
  warning(
    "Skipping '", block_name, "' x '", group_var,
    "': fewer than 3 specimens or fewer than 2 group levels after removing ",
    "missing groups."
  )
  return(NULL)
}

if (verbose) {
  message("Running Kruskal-Wallis for: ", block_name, " x ", group_var)
}

trait_cols <- setdiff(names(block_data), "group__")

# ---- Kruskal-Wallis per trait -----------------------------------------------

results <- lapply(trait_cols, function(trait) {
  y <- block_data[[trait]]
  if (sum(!is.na(y)) < 3) return(NULL)

  tryCatch({
    fit <- stats::kruskal.test(
      stats::as.formula(paste0("`", trait, "` ~ group__")),
      data = block_data
    )
    data.frame(
      trait = trait,
      H_stat = round(fit$statistic, 4),
      df = fit$parameter,
      p_value = fit$p.value,
      row.names = NULL
    )
  }, error = function(e) NULL)
})

results <- do.call(rbind, Filter(Negate(is.null), results))

if (is.null(results) || nrow(results) == 0) {
  warning(
    "Skipping '", block_name, "' x '", group_var,
    "': no trait had enough data for a valid test."
  )
  return(NULL)
}

# ---- Adjust p-values ---------------------------------------------------------

results$p_adjusted <- stats::p.adjust(results$p_value, method = kw_p_adjust)
results$significant <- results$p_adjusted < kw_alpha
results$block <- block_name
results$group <- group_var

n_sig <- sum(results$significant, na.rm = TRUE)

if (verbose) {
  message(
    "  Traits tested: ", nrow(results),
    " | Significant (adj. p < ", kw_alpha, "): ", n_sig
  )
}

# ---- H-statistic bar chart ---------------------------------------------------

plot_df <- results[order(results$H_stat, decreasing = TRUE), ]
plot_df$trait_lab <- trait_label(plot_df$trait)
plot_df$trait_lab <- factor(plot_df$trait_lab, levels = rev(plot_df$trait_lab))
plot_df$sig_label <- ifelse(
  is.na(plot_df$significant), "n.s.",
  ifelse(plot_df$significant, "sig", "n.s.")
)

p_bar <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = .data[["trait_lab"]],
    y = .data[["H_stat"]],
    fill = .data[["sig_label"]]
  )
) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(
    values = c("n.s." = "gray70", "sig" = "#d6604d"),
    labels = c("n.s." = "n.s.", "sig" = paste0("adj. p < ", kw_alpha))
  ) +
  ggplot2::coord_flip(clip = "off") +
  ggplot2::labs(
    x = NULL, y = "H statistic", fill = "",
    title = paste0("Kruskal-Wallis - ", block_name, " | ", group_var)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 7)
  )

fname <- paste0("KW_", block_name, "_", group_var)

save_outputs(
  output_dir,
  fname,
  p_bar,
  base_height = max(4, nrow(plot_df) * 0.3),
  base_aspect_ratio = 1.5
)

if (save_xlsx) {
  openxlsx::write.xlsx(
    results[, c("trait", "H_stat", "df", "p_value", "p_adjusted", "significant")],
    file.path(output_dir, paste0(fname, ".xlsx"))
  )
}

# ---- Dunn post-hoc for significant traits ------------------------------------

posthoc_df <- NULL

if (kw_posthoc && n_sig > 0) {

  sig_traits <- results$trait[results$significant]
  dunn_dir <- file.path(output_dir, paste0("posthoc_", block_name, "_", group_var))
  dir.create(dunn_dir, recursive = TRUE, showWarnings = FALSE)

  all_dunn <- lapply(sig_traits, run_dunn_trait,
                     block_data = block_data,
                     dunn_dir = dunn_dir,
                     block_name = block_name,
                     group_var = group_var
  )

  all_dunn <- Filter(Negate(is.null), all_dunn)

  if (length(all_dunn) > 0) {
    posthoc_df <- do.call(rbind, all_dunn)
    posthoc_df$block <- block_name
    posthoc_df$group <- group_var
    rownames(posthoc_df) <- NULL

    if (save_xlsx) {
      openxlsx::write.xlsx(
        posthoc_df,
        file.path(dunn_dir, "Dunn_all_results.xlsx")
      )
    }
  }
}

if (verbose) {
  message("  Done: ", block_name, " x ", group_var)
}

list(kw = results, posthoc = posthoc_df)
}

# ---- Run all block x group combinations ---------------------------------------

all_kw <- list()
all_posthoc <- list()

for (group_var in kw_groups) {
  for (block_name in names(run_blocks)) {
    key <- paste0(block_name, "_", group_var)
    block_result <- run_kw_block(block_name, run_blocks[[block_name]], group_var)

    if (!is.null(block_result)) {
      all_kw[[key]] <- block_result$kw
      if (!is.null(block_result$posthoc)) {
        all_posthoc[[key]] <- block_result$posthoc
      }
    }
  }
}

combined_kw <- do.call(rbind, all_kw)

if (is.null(combined_kw) || nrow(combined_kw) == 0) {
  warning("No block x group combination could be analyzed.")
  return(invisible(NULL))
}

rownames(combined_kw) <- NULL

combined_posthoc <- do.call(rbind, all_posthoc)
if (!is.null(combined_posthoc)) {
  rownames(combined_posthoc) <- NULL
}

if (save_xlsx) {
  openxlsx::write.xlsx(combined_kw, file.path(output_dir, "KW_all_results.xlsx"))

  if (!is.null(combined_posthoc)) {
    openxlsx::write.xlsx(
      combined_posthoc,
      file.path(output_dir, "KW_posthoc_all_results.xlsx")
    )
  }
}

if (verbose) {
  message("Kruskal-Wallis analysis done. Outputs saved in: ", output_dir)
}

invisible(list(kw = combined_kw, posthoc = combined_posthoc))
}
