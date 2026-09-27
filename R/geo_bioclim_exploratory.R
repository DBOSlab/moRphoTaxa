#' Rank bioclimatic variables by variance, redundancy, and PCA contribution
#'
#' @description
#' Takes the occurrence + climate data frame produced by
#' [geo_matrix_setting()] and evaluates which of the 19 WorldClim
#' bioclimatic variables (`bio_1`-`bio_19`) are most informative for
#' downstream ecogeographic/niche analyses. It computes per-variable
#' variance, a pairwise correlation matrix, and a PCA on the standardized
#' bioclimatic variables, combines these into a single ranked table of
#' variable importance, and draws a bioclimatic-space scatter plot of the
#' occurrence records on two user-chosen bioclimatic axes.
#'
#' @details
#' Only the 19 bioclimatic variables (`bio_1`-`bio_19`) are evaluated here.
#' `"elevation"`, which [geo_matrix_setting()] extracts and stores as a
#' separate, non-bioclimatic column, is deliberately excluded from the
#' variance/correlation/PCA analysis. (An earlier, ad hoc version of this
#' script located bioclimatic columns with the pattern `"bio_0"`.."bio_19"`,
#' as if elevation were a twentieth bioclimatic variable folded into the
#' stack under a spurious `"bio_0"` label. Since [geo_matrix_setting()] now
#' extracts elevation as its own `"elevation"` column instead, that pattern
#' no longer matches anything and is not used here.) The bioclimatic
#' columns to analyze are instead taken directly from
#' `attr(geodata, "env_cols")$clim`, with `"elevation"` dropped, rather
#' than being re-derived from a column-name regular expression. This means
#' `geo_bioclim_exploratory()` always matches whatever climate columns
#' [geo_matrix_setting()] actually extracted, instead of relying on a
#' naming convention that can drift out of sync with it.
#'
#' Four blocks of output are produced:
#' \enumerate{
#'   \item **Variance + correlation** - a horizontal bar chart of
#'   per-variable variance and a correlation heatmap of all retained
#'   bioclimatic variables, combined into one two-panel figure
#'   (`"BIOCLIM_variance_correlation"`).
#'   \item **PCA** - a scree plot (percentage and cumulative variance
#'   explained per principal component), a variable biplot (colored by
#'   contribution to PC1 + PC2), and loading bar charts for PC1 and PC2,
#'   combined into one four-panel figure (`"BIOCLIM_PCA"`).
#'   \item **Bioclimatic scatter plot** - occurrence records plotted on the
#'   `bio_x` and `bio_y` bioclimatic axes, colored by `taxon`, to visualize
#'   how taxa separate along two chosen variables
#'   (`"BIOCLIM_scatter_bixbioy"`).
#'   \item **Ranked variable table** - one row per bioclimatic variable,
#'   ordered by an importance score (the Euclidean norm of its PC1 and PC2
#'   loadings), together with variance, PC1/PC2 loadings, percentage
#'   contribution to PC1 and PC2, and a companion worksheet flagging
#'   variable pairs with `|r| > cor_threshold` as redundant.
#' }
#'
#' Records with a missing value in any retained bioclimatic column are
#' dropped before the variance/correlation/PCA calculations (via
#' `complete.cases()`). The scatter plot instead uses every record with
#' non-missing `bio_x` and `bio_y` values, since it does not depend on the
#' full bioclimatic matrix being complete.
#'
#' If `geodata` carries a `"taxon_colors"` attribute (as produced by
#' [morph_matrix_setting()]), that taxon -> color mapping is reused for the
#' scatter plot so taxon colors stay consistent with other plots in the
#' package; otherwise colors are assigned automatically from the `viridis`
#' palette.
#'
#' All figures are written to a dated subfolder
#' (`format(Sys.time(), "%d%b%Y")`) inside `"Figs.bioclimPCA"`, using the
#' formats selected by `file_formats`. The ranked table and redundant-pairs
#' sheet are written to the same folder as a single `.xlsx` workbook.
#' The dated subfolder means repeated runs do not overwrite each other.
#'
#' @param geodata A data frame produced by [geo_matrix_setting()], i.e. one
#' containing a `"taxon"` column, the bioclimatic columns (`bio_1`-`bio_19`),
#' and the `"env_cols"` attribute produced by that function. `geodata` must
#' still carry its `"env_cols"` attribute (do not strip attributes, e.g. by
#' rebuilding a plain copy with `as.data.frame()`/`data.frame()`, before
#' calling this function).
#'
#' @param bio_x Integer between `1` and `19`. The bioclimatic variable
#' (`bio<bio_x>`) plotted on the x-axis of the bioclimatic scatter plot.
#' Default is `15` (Precipitation Seasonality).
#'
#' @param bio_y Integer between `1` and `19`. The bioclimatic variable
#' (`bio<bio_y>`) plotted on the y-axis of the bioclimatic scatter plot.
#' Default is `12` (Annual Precipitation).
#'
#' @param cor_threshold Numeric between `0` and `1`. Absolute correlation
#' above which a pair of bioclimatic variables is flagged as redundant in
#' the `"redundant_pairs"` worksheet. Default is `0.8`.
#'
#' @param file_formats Character vector specifying the file formats used to
#' save figures. One or more of `"pdf"` and `"jpg"`. Default is
#' `c("pdf", "jpg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages
#' and a console summary of the top 5 most informative variables.
#'
#' @return
#' `NULL`, invisibly. This function is called for its side effects: the
#' variance/correlation figure, the PCA figure, the bioclimatic scatter
#' plot, and the ranked-variable workbook are written to a dated subfolder
#' of `"Figs.bioclimPCA"` (see Details), and (when `verbose = TRUE`) a
#' summary of the top 5 most informative variables is printed to the
#' console. The ranked table
#' and other intermediate objects (`bio_matrix`, `cor_mat`, `pca_res`,
#' `summary_table`) stay internal to the function.
#'
#' @seealso
#' [geo_matrix_setting()], [morph_matrix_setting()]
#'
#' @examples
#' \dontrun{
#' geo <- geo_matrix_setting(
#'   xlsx_path = "occurrence_data.xlsx",
#'   taxon_col = "species",
#'   edaphic = FALSE
#' )
#'
#' geo_bioclim_exploratory(
#'   geodata = geo,
#'   bio_x = 15,
#'   bio_y = 12
#' )
#' # Figures + ranked-variable workbook are written to Figs.bioclimPCA/<date>/
#'
#' # Save PDF only
#' geo_bioclim_exploratory(
#'   geodata = geo,
#'   bio_x = 15,
#'   bio_y = 12,
#'   file_formats = "pdf"
#' )
#'
#' # Save JPEG only
#' geo_bioclim_exploratory(
#'   geodata = geo,
#'   bio_x = 15,
#'   bio_y = 12,
#'   file_formats = "jpg"
#' )
#' }
#'
#' @importFrom stats complete.cases cor prcomp var setNames
#' @importFrom utils head
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_col geom_tile geom_text geom_line
#' @importFrom ggplot2 geom_point coord_flip labs theme_bw theme
#' @importFrom ggplot2 element_text element_blank scale_fill_gradient2
#' @importFrom viridis scale_fill_viridis viridis
#' @importFrom reshape2 melt
#' @importFrom cowplot plot_grid save_plot
#' @importFrom openxlsx createWorkbook addWorksheet writeData addStyle
#' @importFrom openxlsx createStyle saveWorkbook
#'
#' @export

geo_bioclim_exploratory <- function(geodata,
                                    bio_x = 15,
                                    bio_y = 12,
                                    cor_threshold = 0.8,
                                    file_formats = c("pdf", "jpg"),
                                    verbose = TRUE) {

# ---- Validate input ---------------------------------------------------------

if (!is.data.frame(geodata)) {
  stop("`geodata` must be a data frame, typically the output of geo_matrix_setting().")
}

env_cols <- attr(geodata, "env_cols")

if (is.null(env_cols) || is.null(env_cols$clim)) {
  stop(
    "`geodata` has no \"env_cols\" attribute with a \"clim\" element. ",
    "Pass the data frame returned by geo_matrix_setting() unchanged."
  )
}

if (!"taxon" %in% names(geodata)) {
  stop("`geodata` must contain a \"taxon\" column.")
}

if (!is.numeric(bio_x) || length(bio_x) != 1L || !(bio_x %in% 1:19)) {
  stop("`bio_x` must be a single integer between 1 and 19.")
}

if (!is.numeric(bio_y) || length(bio_y) != 1L || !(bio_y %in% 1:19)) {
  stop("`bio_y` must be a single integer between 1 and 19.")
}

if (!is.numeric(cor_threshold) || length(cor_threshold) != 1L ||
    cor_threshold < 0 || cor_threshold > 1) {
  stop("`cor_threshold` must be a single number between 0 and 1.")
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

# ---- Output folder ------------------------------------------------------------

folder_name <- file.path("Figs.bioclimPCA", format(Sys.time(), "%d%b%Y"))
dir.create(folder_name, recursive = TRUE, showWarnings = FALSE)

#_______________________________________________________________________________
# Extract BIO columns ####

# Bioclimatic columns come straight from the "env_cols" attribute attached by
# geo_matrix_setting(), rather than a column-name regular expression, so this
# function always matches whatever climate columns were actually extracted.
# "elevation" is a separate, non-bioclimatic layer and is excluded here.
bio_cols <- setdiff(env_cols$clim, "elevation")

if (length(bio_cols) == 0) {
  stop("No bioclimatic (bio_1-bio_19) columns found in `geodata`.")
}

bio_x_col <- paste0("bio_", bio_x)
bio_y_col <- paste0("bio_", bio_y)

if (!(bio_x_col %in% names(geodata))) {
  stop("`bio_x` = ", bio_x, " refers to column \"", bio_x_col,
       "\", which is not present in `geodata`.")
}

if (!(bio_y_col %in% names(geodata))) {
  stop("`bio_y` = ", bio_y, " refers to column \"", bio_y_col,
       "\", which is not present in `geodata`.")
}

bio_matrix <- geodata[, bio_cols, drop = FALSE]
bio_matrix <- bio_matrix[stats::complete.cases(bio_matrix), ]

if (nrow(bio_matrix) < 3) {
  stop("Fewer than 3 complete records remain across the bioclimatic variables; ",
       "cannot compute variance/correlation/PCA.")
}

if (verbose) {
  message("Ranking bioclimatic variables...")
  message("  Bioclimatic variables: ", length(bio_cols))
  message("  Complete records used for variance/correlation/PCA: ", nrow(bio_matrix))
}

# Short display names: BIO1 ... BIO19 (elevation excluded, see Details)
colnames(bio_matrix) <- paste0(
  "BIO",
  gsub("[^0-9]", "", bio_cols)
)

bio_desc <- c(
  BIO1 = "Annual Mean Temperature",        BIO2 = "Mean Diurnal Range",
  BIO3 = "Isothermality",                  BIO4 = "Temperature Seasonality",
  BIO5 = "Max Temp Warmest Month",         BIO6 = "Min Temp Coldest Month",
  BIO7 = "Temperature Annual Range",       BIO8 = "Mean Temp Wettest Quarter",
  BIO9 = "Mean Temp Driest Quarter",       BIO10 = "Mean Temp Warmest Quarter",
  BIO11 = "Mean Temp Coldest Quarter",     BIO12 = "Annual Precipitation",
  BIO13 = "Precipitation Wettest Month",   BIO14 = "Precipitation Driest Month",
  BIO15 = "Precipitation Seasonality",     BIO16 = "Precipitation Wettest Quarter",
  BIO17 = "Precipitation Driest Quarter",  BIO18 = "Precipitation Warmest Quarter",
  BIO19 = "Precipitation Coldest Quarter"
)

#_______________________________________________________________________________
# PART 1 - Variance bar chart + correlation heatmap ####

var_df <- data.frame(
  variable = colnames(bio_matrix),
  variance = apply(bio_matrix, 2, stats::var, na.rm = TRUE)
)

p_variance <- ggplot2::ggplot(
  var_df,
  ggplot2::aes(
    x = stats::reorder(.data$variable, .data$variance),
    y = .data$variance,
    fill = .data$variance
  )
) +
  ggplot2::geom_col(width = 0.7) +
  viridis::scale_fill_viridis(option = "plasma", direction = -1, guide = "none") +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Variance per bioclimatic variable",
    subtitle = "Higher = more discriminating power across specimens",
    x = NULL, y = "Variance"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9),
    panel.grid.major.y = ggplot2::element_blank()
  )

cor_mat <- stats::cor(bio_matrix, use = "pairwise.complete.obs")
cor_melt <- reshape2::melt(cor_mat)
names(cor_melt) <- c("Var1", "Var2", "r")

p_cor <- ggplot2::ggplot(
  cor_melt,
  ggplot2::aes(x = .data$Var1, y = .data$Var2, fill = .data$r)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = round(.data$r, 1)),
                     size = 2.2, colour = "white") +
  viridis::scale_fill_viridis(option = "inferno", limits = c(-1, 1), name = "r") +
  ggplot2::labs(
    title = "Bioclimatic variable correlation matrix",
    subtitle = paste0("Variables with |r| > ", cor_threshold, " are largely redundant"),
    x = NULL, y = NULL
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = ggplot2::element_text(size = 8),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9)
  )

fig_part1 <- cowplot::plot_grid(
  p_variance, p_cor,
  labels = c("a)", "b)"), label_size = 16,
  ncol = 2, rel_widths = c(0.4, 0.6)
)

#_______________________________________________________________________________
# PART 2 - PCA on BIO variables ####

pca_res <- stats::prcomp(bio_matrix, scale. = TRUE)

# Scree plot
eig_df <- data.frame(
  PC = factor(paste0("PC", seq_along(pca_res$sdev)),
              levels = paste0("PC", seq_along(pca_res$sdev))),
  var_pct = (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
)
eig_df$cumulative <- cumsum(eig_df$var_pct)

p_scree <- ggplot2::ggplot(eig_df, ggplot2::aes(x = .data$PC)) +
  ggplot2::geom_col(ggplot2::aes(y = .data$var_pct, fill = .data$var_pct), width = 0.6) +
  ggplot2::geom_line(ggplot2::aes(y = .data$cumulative, group = 1),
                     colour = "red", linewidth = 0.7) +
  ggplot2::geom_point(ggplot2::aes(y = .data$cumulative), colour = "red", size = 2) +
  viridis::scale_fill_viridis(option = "plasma", direction = -1, guide = "none") +
  ggplot2::labs(
    title = "Scree plot - PCA on bioclimatic variables",
    subtitle = "Bars = % variance per PC  |  Red = cumulative",
    x = NULL, y = "Variance explained (%)"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9)
  )

# Variable biplot
.quiet_load_namespaces("factoextra")
p_biplot <- factoextra::fviz_pca_var(
  pca_res,
  col.var = "contrib",
  gradient.cols = c("#440154", "#21908C", "#FDE725"),
  repel = TRUE,
  labelsize = 3
) +
  ggplot2::labs(
    title = "PCA variable biplot",
    subtitle = "Colour = contribution to PC1 + PC2",
    color = "Contrib (%)"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9),
    panel.grid = ggplot2::element_blank()
  )

# PC1 and PC2 loading bars
make_loading_bar <- function(pc_col, pc_label) {
  df <- data.frame(
    variable = rownames(pca_res$rotation),
    loading = pca_res$rotation[, pc_col]
  )
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = stats::reorder(.data$variable, abs(.data$loading)),
      y = .data$loading,
      fill = .data$loading
    )
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "gray90", high = "#d6604d",
                                  midpoint = 0, guide = "none") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = paste(pc_label, "loadings"), x = NULL, y = "Loading") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid.major.y = ggplot2::element_blank()
    )
}

fig_part2 <- cowplot::plot_grid(
  cowplot::plot_grid(p_scree, p_biplot,
                     labels = c("a)", "b)"), label_size = 16, ncol = 2),
  cowplot::plot_grid(make_loading_bar("PC1", "PC1"),
                     make_loading_bar("PC2", "PC2"),
                     labels = c("c)", "d)"), label_size = 16, ncol = 2),
  ncol = 1
)

#_______________________________________________________________________________
# PART 3 - Bioclimatic scatter plot (bio_x vs bio_y) ####

scatter_data <- geodata[
  stats::complete.cases(geodata[, c(bio_x_col, bio_y_col)]),
  c("taxon", bio_x_col, bio_y_col)
]

taxon_colors <- attr(geodata, "taxon_colors")
taxa_present <- sort(unique(scatter_data$taxon))

if (is.null(taxon_colors)) {
  taxon_colors <- stats::setNames(viridis::viridis(length(taxa_present)), taxa_present)
} else {
  taxon_colors <- taxon_colors[taxa_present]
}

bio_x_label <- paste0(
  "BIO", bio_x, " - ",
  bio_desc[[paste0("BIO", bio_x)]]
)

bio_y_label <- paste0(
  "BIO", bio_y, " - ",
  bio_desc[[paste0("BIO", bio_y)]]
)

p_scatter <- ggplot2::ggplot(
  scatter_data,
  ggplot2::aes(x = .data[[bio_x_col]], y = .data[[bio_y_col]], colour = .data$taxon)
) +
  ggplot2::geom_point(alpha = 0.8, size = 2) +
  ggplot2::scale_colour_manual(values = taxon_colors, name = "Taxon") +
  ggplot2::labs(
    title = paste0("Bioclimatic space: BIO", bio_x, " vs BIO", bio_y),
    x = bio_x_label, y = bio_y_label
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
  )

#_______________________________________________________________________________
# Ranked summary table ####

contrib_pc1 <- (pca_res$rotation[, 1]^2 / sum(pca_res$rotation[, 1]^2)) * 100
contrib_pc2 <- (pca_res$rotation[, 2]^2 / sum(pca_res$rotation[, 2]^2)) * 100
importance <- sqrt(pca_res$rotation[, 1]^2 + pca_res$rotation[, 2]^2)

vars <- colnames(bio_matrix)
summary_table <- data.frame(
  rank = NA_integer_,
  variable = vars,
  description = bio_desc[vars],
  variance = round(var_df$variance[match(vars, var_df$variable)], 3),
  PC1_loading = round(pca_res$rotation[vars, 1], 3),
  PC2_loading = round(pca_res$rotation[vars, 2], 3),
  contrib_PC1_pct = round(contrib_pc1[vars], 1),
  contrib_PC2_pct = round(contrib_pc2[vars], 1),
  importance_score = round(importance[vars], 3),
  row.names = NULL
)
summary_table <- summary_table[order(summary_table$importance_score, decreasing = TRUE), ]
summary_table$rank <- seq_len(nrow(summary_table))

# Flag redundant pairs
high_cor <- which(abs(cor_mat) > cor_threshold & upper.tri(cor_mat), arr.ind = TRUE)

cor_notes <- if (nrow(high_cor) > 0) {
  apply(high_cor, 1, function(i) {
    paste0(
      rownames(cor_mat)[i[1]], "~",
      colnames(cor_mat)[i[2]],
      " (r=", round(cor_mat[i[1], i[2]], 2), ")"
    )
  })
} else {
  "None"
}

#_______________________________________________________________________________
# Save outputs ####

save_outputs <- function(name, plot, ncol = 1, nrow = 1,
                         base_height = 10, base_aspect_ratio = 1.1) {
  for (ext in paste0(".", file_formats)) {
    cowplot::save_plot(
      file.path(folder_name, paste0(name, ext)), plot,
      ncol = ncol, nrow = nrow, base_height = base_height,
      base_aspect_ratio = base_aspect_ratio, base_width = NULL
    )
  }
}

save_outputs("BIOCLIM_variance_correlation", fig_part1,
             ncol = 2, nrow = 1, base_height = 9)

save_outputs("BIOCLIM_PCA", fig_part2,
             ncol = 2, nrow = 2,
             base_height = 8, base_aspect_ratio = 1.0)

save_outputs(paste0("BIOCLIM_scatter_bio", bio_x, "bio", bio_y),
             p_scatter,
             ncol = 1, nrow = 1,
             base_height = 7, base_aspect_ratio = 1.3)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "ranked_variables")
openxlsx::writeData(wb, "ranked_variables", summary_table)
openxlsx::addStyle(
  wb, "ranked_variables",
  openxlsx::createStyle(fgFill = "#C6EFCE", textDecoration = "bold"),
  rows = 2:6, cols = 1:ncol(summary_table), gridExpand = TRUE
)
openxlsx::addWorksheet(wb, "redundant_pairs")
openxlsx::writeData(
  wb, "redundant_pairs",
  data.frame(highly_correlated_pairs = cor_notes)
)
openxlsx::saveWorkbook(
  wb, file.path(folder_name, "BIOCLIM_variable_ranking.xlsx"),
  overwrite = TRUE
)

if (verbose) {
  message(
    "\nDone. Top 5 most informative variables: ",
    paste(utils::head(summary_table$variable, 5), collapse = ", ")
  )
  message("Outputs saved in: ", folder_name)
}

invisible(NULL)
}
