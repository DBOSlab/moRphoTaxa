#' Lag-1 autocorrelation of morphological traits along a geographic gradient
#'
#' @description
#' Measures, for each morphological trait, how similar the values of
#' neighbouring specimens are when specimens are ordered along a
#' geographic axis (latitude or longitude) or kept in their original
#' spreadsheet order. A positive lag-1 autocorrelation indicates that
#' specimens close to each other in the ordering tend to have similar
#' values (a gradient or geographic structure in the trait), a value near
#' zero indicates no ordering-related structure, and a negative value
#' indicates that neighbouring specimens tend to alternate.
#'
#' Takes the morphometric matrix produced by [morph_matrix_setting()] and,
#' when specimens are ordered by geography, the geographical matrix
#' produced by [geo_matrix_setting()]. For each requested trait block, the
#' function computes the autocorrelation of every trait, draws a ranked
#' bar chart, and writes a combined summary table.
#'
#' @details
#' \strong{Ordering specimens.}
#' Specimens are matched by row names (the specimen identifiers shared by
#' [morph_matrix_setting()] and [geo_matrix_setting()]) and sorted once,
#' before any block is analyzed:
#' \itemize{
#'   \item \code{"latitude"} (default): ascending \code{decimalLatitude};
#'   \item \code{"longitude"}: ascending \code{decimalLongitude};
#'   \item \code{"none"}: the row order of \code{analysis_data} is kept
#'   (the spreadsheet order, e.g. collection number). \code{geodata} is
#'   not needed in this case.
#' }
#' When ordering by latitude or longitude, specimens of \code{analysis_data}
#' that are absent from \code{geodata} (or have a missing coordinate) are
#' dropped. Specimens with identical coordinates keep their relative
#' spreadsheet order.
#'
#' \strong{Lag-1 autocorrelation.}
#' For each trait, the specimens are taken in the chosen order and the
#' trait value of each specimen is correlated with the value of the
#' preceding specimen, using the coefficient given by \code{cor_method}.
#' A pair of consecutive specimens is used only if both have a measured
#' value for that trait; the autocorrelation is \code{NA}, and the trait is
#' left out of the results, when fewer than 3 such pairs exist or when
#' the trait is constant among them. The coefficient is descriptive: no
#' significance test is performed. It only detects structure along the
#' chosen ordering axis.
#'
#' \strong{Pooled taxa.}
#' All specimens in \code{analysis_data} are ordered together. If several
#' taxa are included and they occupy different parts of the ordering axis,
#' differences between taxa contribute to the autocorrelation. To examine
#' geographic structure within a single taxon, build \code{analysis_data}
#' with \code{species_selected} set to that taxon in
#' [morph_matrix_setting()].
#'
#' \strong{Trait blocks.}
#' The vegetative (\code{"veg"}), floral (\code{"flo"}), and fruit
#' (\code{"fru"}) trait blocks are read from the \code{"base_cols"}
#' attribute of \code{analysis_data}, created by [morph_matrix_setting()].
#' With \code{blocks = "all"}, each available block is run, followed by
#' every combination of two or more available blocks, named by
#' concatenating the block names (\code{"vegflo"}, \code{"vegfru"},
#' \code{"flofru"}, \code{"vegflofru"}). Any single block or combination
#' can also be requested by name. Because the autocorrelation is computed
#' trait by trait, a trait has the same value in every block that contains
#' it; combined blocks only group traits in a single chart and table.
#'
#' \strong{Outputs.}
#' One chart per block (\code{Autocorrelation_<block>.<ext>}) and a single
#' \code{Autocorrelation_summary.csv} with all blocks are written to a
#' dated subfolder (\code{format(Sys.time(), "\%d\%b\%Y")}) inside
#' \code{"Figs_Autocorrelation"}. The dated subfolder means repeated runs
#' do not overwrite each other. Trait labels in the charts and the summary
#' are the descriptive part of the column name, with underscores replaced
#' by spaces (\code{"petiole_length/PETIlng"} becomes
#' \code{"petiole length"}); columns already renamed to codes are shown as
#' they are.
#'
#' @param analysis_data A data frame produced by [morph_matrix_setting()],
#' with specimen identifiers as row names and the \code{"base_cols"}
#' attribute. Do not strip attributes (e.g. by subsetting or rebuilding the
#' data frame) before calling this function.
#'
#' @param geodata A data frame produced by [geo_matrix_setting()], with
#' specimen identifiers as row names and the \code{decimalLatitude} and
#' \code{decimalLongitude} columns. Required when \code{order_by} is
#' \code{"latitude"} or \code{"longitude"}; ignored when \code{order_by =
#' "none"}. Default is \code{NULL}.
#'
#' @param blocks Character vector of blocks to analyze. \code{"all"}
#' (default) runs every available block and every combination of blocks.
#' Otherwise, one or more of \code{"veg"}, \code{"flo"}, \code{"fru"},
#' \code{"vegflo"}, \code{"vegfru"}, \code{"flofru"}, and
#' \code{"vegflofru"}, restricted to the blocks available in
#' \code{analysis_data}. Unavailable names are ignored with a warning.
#'
#' @param cor_method Character string giving the correlation coefficient
#' used for the lag-1 autocorrelation. One of \code{"pearson"} (default),
#' \code{"spearman"}, or \code{"kendall"}. Passed to [stats::cor()].
#'
#' @param order_by Character string giving how specimens are ordered before
#' computing the autocorrelation. One of \code{"latitude"} (default),
#' \code{"longitude"}, or \code{"none"} (keep the row order of
#' \code{analysis_data}).
#'
#' @param file_formats Character vector specifying the file formats used to
#' save the charts. One or more of \code{"pdf"} and \code{"jpg"}. Default
#' is \code{c("pdf", "jpg")}.
#'
#' @param verbose Logical. If \code{TRUE} (default), prints progress
#' messages and, for each block, the number of traits computed and the
#' traits with the highest and lowest autocorrelation.
#'
#' @return
#' A list, returned invisibly, with the elements:
#' \describe{
#'   \item{\code{results}}{A data frame with one row per trait and block,
#'   with columns \code{block}, \code{trait}, and \code{autocor}, or
#'   \code{NULL} if no block produced results.}
#'   \item{\code{plots}}{A named list of [ggplot2::ggplot()] objects, one per
#'   block.}
#'   \item{\code{specimen_order}}{Character vector of specimen identifiers in
#'   the order used.}
#'   \item{\code{output_dir}}{Path of the folder where the charts and the
#'   summary table were saved.}
#' }
#' The function is also called for its side effects: charts and the summary
#' table are written to \code{output_dir}.
#'
#' @seealso
#' [morph_matrix_setting()], [geo_matrix_setting()],
#' [geo_morph_mantel()]
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
#' # Every block and block combination, specimens ordered by latitude
#' res <- morph_autocorrelation(analysis_data = morph, geodata = geo)
#' head(res$results)
#'
#' # Floral traits only, ordered by longitude, Spearman correlation
#' res <- morph_autocorrelation(
#'   analysis_data = morph,
#'   geodata = geo,
#'   blocks = "flo",
#'   cor_method = "spearman",
#'   order_by = "longitude"
#' )
#'
#' # Keep the original spreadsheet order (geodata not needed)
#' res <- morph_autocorrelation(
#'   analysis_data = morph,
#'   blocks = c("veg", "fru"),
#'   order_by = "none",
#'   file_formats = "pdf"
#' )
#' }
#'
#' @importFrom stats cor reorder
#' @importFrom utils combn write.csv
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_col geom_text geom_hline coord_flip
#' @importFrom ggplot2 scale_y_continuous expansion labs theme_bw theme
#' @importFrom ggplot2 element_blank element_text
#' @importFrom viridis scale_fill_viridis
#' @importFrom cowplot save_plot
#'
#' @export

morph_autocorrelation <- function(analysis_data,
                                  geodata = NULL,
                                  blocks = "all",
                                  cor_method = c("pearson", "spearman", "kendall"),
                                  order_by = c("latitude", "longitude", "none"),
                                  file_formats = c("pdf", "jpg"),
                                  verbose = TRUE) {

# ---- Validate input -------------------------------------------------------

if (!is.data.frame(analysis_data)) {
  stop("`analysis_data` must be a data frame, typically the output of morph_matrix_setting().")
}

base_cols <- attr(analysis_data, "base_cols")

if (is.null(base_cols) || length(base_cols) == 0L) {
  stop(
    "`analysis_data` has no \"base_cols\" attribute. ",
    "Pass the data frame returned by morph_matrix_setting() unchanged."
  )
}

cor_method <- match.arg(cor_method)
order_by   <- match.arg(order_by)

if (order_by != "none") {
  if (!is.data.frame(geodata)) {
    stop(
      "`geodata` (the output of geo_matrix_setting()) is required when ",
      "`order_by` is \"", order_by, "\". Use order_by = \"none\" to keep the ",
      "row order of `analysis_data`."
    )
  }
  if (!all(c("decimalLatitude", "decimalLongitude") %in% names(geodata))) {
    stop("`geodata` must contain `decimalLatitude` and `decimalLongitude` columns.")
  }
}

if (!is.character(blocks) || length(blocks) < 1L) {
  stop("`blocks` must be a character vector with at least one block name.")
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

# ---- Order specimens (once, for all blocks) ---------------------------------

specimen_ids <- rownames(analysis_data)

if (order_by != "none") {
  geo_col    <- if (order_by == "latitude") "decimalLatitude" else "decimalLongitude"
  common_ids <- intersect(specimen_ids, rownames(geodata))

  if (length(common_ids) == 0L) {
    stop(
      "No specimen identifiers are shared by `analysis_data` and `geodata`. ",
      "Both must be built from the same specimens (row names are the specimen IDs)."
    )
  }

  coord_vals <- as.numeric(geodata[common_ids, geo_col])
  has_coord  <- !is.na(coord_vals)
  common_ids <- common_ids[has_coord]
  coord_vals <- coord_vals[has_coord]

  if (verbose && length(common_ids) < length(specimen_ids)) {
    message(
      "  Specimens without coordinates in `geodata` dropped: ",
      length(specimen_ids) - length(common_ids)
    )
  }

  specimen_ids <- common_ids[order(coord_vals)]
}

if (length(specimen_ids) < 4L) {
  stop("At least 4 specimens are required to compute lag-1 autocorrelation.")
}

ordered_data <- analysis_data[specimen_ids, , drop = FALSE]

# ---- Output folder ------------------------------------------------------------

output_dir <- file.path("Figs_Autocorrelation", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (verbose) {
  message("Running autocorrelation analysis...")
  message("  Specimens used: ", length(specimen_ids), " (ordered by ", order_by, ")")
  message("  Blocks to run: ", paste(names(run_blocks), collapse = ", "))
}

# ---- Run each block -----------------------------------------------------------

results_list <- list()
plots        <- list()

for (block_name in names(run_blocks)) {

  if (verbose) {
    message("\n  --- Block: ", block_name, " ---")
  }

  valid_cols <- intersect(run_blocks[[block_name]], names(ordered_data))
  valid_cols <- valid_cols[vapply(ordered_data[valid_cols], is.numeric, logical(1))]

  if (length(valid_cols) == 0L) {
    warning("Skipping '", block_name, "': no valid columns.")
    next
  }

  autocor <- vapply(
    valid_cols,
    function(col) .lag1_autocor(ordered_data[[col]], method = cor_method),
    numeric(1)
  )

  df_plot <- data.frame(
    trait   = gsub("_", " ", gsub("/.*", "", valid_cols)),
    autocor = unname(autocor),
    stringsAsFactors = FALSE
  )
  df_plot <- df_plot[!is.na(df_plot$autocor), , drop = FALSE]
  df_plot <- df_plot[order(df_plot$autocor, decreasing = TRUE), , drop = FALSE]

  if (nrow(df_plot) == 0L) {
    warning("Skipping '", block_name, "': no traits with enough data.")
    next
  }

  if (verbose) {
    message("  Traits computed: ", nrow(df_plot))
    message(
      "  Highest autocorrelation: ", round(max(df_plot$autocor), 3),
      " (", df_plot$trait[1], ")"
    )
    message(
      "  Lowest autocorrelation : ", round(min(df_plot$autocor), 3),
      " (", df_plot$trait[nrow(df_plot)], ")"
    )
  }

  # ---- Plot -------------------------------------------------------------------

  label_hjust <- ifelse(df_plot$autocor >= 0, -0.1, 1.1)

  p_autocor <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x    = stats::reorder(.data$trait, .data$autocor),
      y    = .data$autocor,
      fill = .data$autocor
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      ggplot2::aes(label = round(.data$autocor, 3)),
      hjust = label_hjust, size = 3
    ) +
    viridis::scale_fill_viridis(option = "mako", discrete = FALSE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.12)) +
    ggplot2::labs(
      x     = NULL,
      y     = "Lag-1 Autocorrelation",
      fill  = "Autocorrelation",
      title = paste0("Autocorrelation - ", block_name)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(size = 7)
    )

  for (ext in file_formats) {
    cowplot::save_plot(
      file.path(output_dir, paste0("Autocorrelation_", block_name, ".", ext)),
      p_autocor,
      ncol = 1, nrow = 1,
      base_height = max(4, nrow(df_plot) * 0.3),
      base_aspect_ratio = 1.5
    )
  }

  plots[[block_name]] <- p_autocor

  results_list[[block_name]] <- data.frame(
    block   = block_name,
    trait   = df_plot$trait,
    autocor = df_plot$autocor,
    stringsAsFactors = FALSE
  )
}

# ---- Save combined summary ------------------------------------------------------

if (length(results_list) > 0L) {
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL

  utils::write.csv(
    results,
    file.path(output_dir, "Autocorrelation_summary.csv"),
    row.names = FALSE
  )
} else {
  results <- NULL
  warning("No block produced results; no summary table was written.")
}

if (verbose) {
  message("\nDone. Outputs saved in: ", output_dir)
}

invisible(list(
  results        = results,
  plots          = plots,
  specimen_order = specimen_ids,
  output_dir     = output_dir
))
}

#' Lag-1 autocorrelation of a single trait
#'
#' @description
#' Internal helper used by [morph_autocorrelation()]. Correlates each value
#' of `x` with the value that precedes it, using only pairs in which both
#' values are observed.
#'
#' @param x Numeric vector of trait values, already in the desired order.
#' @param method Character. Correlation method passed to [stats::cor()].
#'
#' @return A single number, or `NA_real_` when fewer than 3 complete pairs
#' exist or the correlation is undefined (e.g. constant values).
#'
#' @keywords internal
#'
#' @noRd
.lag1_autocor <- function(x, method) {
n <- length(x)
if (n < 2L) return(NA_real_)

x_now <- x[-1L]
x_lag <- x[-n]

valid <- !is.na(x_now) & !is.na(x_lag)
if (sum(valid) < 3L) return(NA_real_)

suppressWarnings(
  stats::cor(x_now[valid], x_lag[valid], method = method)
)
}
