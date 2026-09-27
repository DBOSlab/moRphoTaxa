#' Plot morphological trait distributions among taxa
#'
#' @description
#' Generates boxplots, violin plots, or both for morphological traits measured
#' across taxa. The function is designed for exploratory visualization of
#' morphometric variation using the objects produced by
#' [morph_matrix_setting()].
#'
#' For each morphological trait with a sufficient number of non-missing
#' observations, the function creates a separate plot showing the distribution
#' of trait values among taxa and saves the resulting figures to disk.
#'
#' @details
#' `morph_boxplots()` combines a numeric morphometric matrix with the taxonomic
#' identity of each specimen and generates one plot per trait.
#'
#' Three plotting modes are available through `plot_type`:
#'
#' \itemize{
#'   \item `"boxplot"`: generates standard boxplots;
#'   \item `"violin"`: generates violin plots with an embedded boxplot;
#'   \item `"both"`: generates both plot types for each eligible trait.
#' }
#'
#' When `show_jitter = TRUE`, individual specimen observations are displayed
#' as jittered points over the corresponding distributions.
#'
#' Traits with fewer than `min_n_total` non-missing observations are skipped.
#' This threshold is applied across all retained taxa combined.
#'
#' Taxon colours are taken from the `"taxon_colors"` attribute attached to
#' `analysis_data` by [morph_matrix_setting()]. If this attribute is missing,
#' colours are generated automatically from the `viridis` palette.
#' Taxon names are displayed on the x-axis and morphological trait values on
#' the y-axis.
#'
#' Output files are written to a date-specific directory inside
#' `Figs.boxplot`. Each eligible trait is saved in the formats specified by
#' `file_formats` for the selected plot type or types.
#'
#' File names are generated from the trait names after replacing characters
#' other than letters, numbers, and underscores with underscores.
#'
#' @param analysis_data A numeric data frame or matrix containing morphological
#' traits as columns and specimens as rows. Row names must correspond to the
#' specimen identifiers used in `species_all`. Typically obtained from the
#' `analysis_data` element returned by [morph_matrix_setting()]. If this
#' object carries a `"taxon_colors"` attribute, those colours are reused;
#' otherwise a `viridis` palette is generated automatically.
#'
#' @param plot_type Character string specifying the type of distribution plot
#' to generate. One of `"boxplot"`, `"violin"`, or `"both"`. Default is
#' `"violin"`.
#'
#' @param show_jitter Logical. If `TRUE` (default), individual specimen
#' observations are displayed as jittered points over boxplots or violin plots.
#'
#' @param min_n_total Integer. Minimum total number of non-missing observations
#' required for a morphological trait to be plotted. Traits with fewer
#' observations are skipped. Default is `10`.
#'
#' @param file_formats Character vector specifying which image formats to
#' save. One or more of `"pdf"` and `"jpeg"`. Default is
#' `c("pdf", "jpeg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' including skipped traits and the output directory, to the console.
#'
#' @return
#' Invisibly returns `NULL`.
#'
#' The function is primarily called for its side effects: one or more PDF and/or
#' JPEG figures are written to the output directory for each eligible
#' morphological trait.
#'
#' @seealso
#' [morph_matrix_setting()],
#' [morph_normality()]
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
#' # Generate violin plots for all traits
#' morph_boxplots(
#'   analysis_data = res$analysis_data,
#'   species_all = res$species_all,
#'   plot_type = "violin",
#'   show_jitter = TRUE,
#'   min_n_total = 10
#' )
#'
#' # Generate both boxplots and violin plots for selected taxa
#' morph_boxplots(
#'   analysis_data = res$analysis_data,
#'   species_all = res$species_all,
#'   species_selected = c("Species A", "Species B"),
#'   plot_type = "both",
#'   show_jitter = TRUE,
#'   min_n_total = 10
#' )
#'
#' # Generate violin plots as PDF only
#' morph_boxplots(
#'   analysis_data = res$analysis_data,
#'   plot_type = "violin",
#'   show_jitter = TRUE,
#'   min_n_total = 10,
#'   file_formats = "pdf"
#' )
#'
#' # Generate both boxplots and violin plots as PDF and JPEG
#' morph_boxplots(
#'   analysis_data = res$analysis_data,
#'   plot_type = "both",
#'   show_jitter = TRUE,
#'   min_n_total = 10,
#'   file_formats = c("pdf", "jpeg")
#' )
#' }
#'
#' @importFrom cowplot save_plot
#' @importFrom ggplot2 aes element_blank element_text geom_boxplot geom_jitter
#'   geom_violin ggplot scale_color_manual scale_fill_manual theme theme_bw
#'   xlab ylab
#' @importFrom grid unit
#' @importFrom viridis viridis
#'
#' @export
#'

morph_boxplots <- function(analysis_data = NULL,
                           plot_type = "violin",
                           show_jitter = TRUE,
                           min_n_total = 10,
                           file_formats = c("pdf", "jpeg"),
                           verbose = TRUE) {


# ---- Validate input ------------------------------------------------------

plot_type <- match.arg(
  plot_type,
  choices = c("boxplot", "violin", "both")
)

file_formats <- match.arg(
  file_formats,
  choices = c("pdf", "jpeg"),
  several.ok = TRUE
)

if (!is.logical(show_jitter) || length(show_jitter) != 1L) {
  stop("`show_jitter` must be TRUE or FALSE.")
}

if (!is.numeric(min_n_total) ||
    length(min_n_total) != 1L ||
    min_n_total < 1) {
  stop("`min_n_total` must be a positive integer.")
}

if (!is.data.frame(analysis_data)) {
  stop("`analysis_data` must be a data frame.")
}

if (!"taxon" %in% names(analysis_data)) {
  stop("`analysis_data` must contain a `taxon` column")
}

if (is.null(rownames(analysis_data))) {
  stop("`analysis_data` must have specimen IDs as row names.")
}

if (verbose) {
  message("Running Boxplot analysis")
}

# ---- Output folder -------------------------------------------------------

if (!exists("output_dir")) {
  output_dir <- file.path(
    "Figs.boxplot",
    format(Sys.time(), "%d%b%Y")
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ---- Build plotting dataframe -------------------------------------------

bp_data <- analysis_data

bp_data$taxon <- factor(bp_data$taxon)

trait_vars <- setdiff(
  names(bp_data),
  c("id", "taxon")
)

all_species <- sort(unique(bp_data$taxon))

# ---- Taxon colours ---------------------------------------------------------
color_palette <- attr(analysis_data, "taxon_colors")

if (is.null(color_palette)) {
  color_palette <- stats::setNames(
    viridis::viridis(length(all_species)),
    as.character(all_species)
  )
} else {
  missing_sp <- setdiff(as.character(all_species), names(color_palette))
  if (length(missing_sp) > 0) {
    warning("No stored colour found for: ", paste(missing_sp, collapse = ", "),
            ". Assigning automatic colours for these taxa.")
    color_palette[missing_sp] <- viridis::viridis(length(missing_sp))
  }
}

# ---- Helper: clean y-axis label ------------------------------------------

clean_label <- function(x) {
  gsub("[_/]", " ", x)
}

# ---- Helper: build one boxplot -------------------------------------------

make_boxplot <- function(df, traitname, color_palette, show_jitter) {

  df <- df[!is.na(df[[traitname]]), ]

  if (nrow(df) < 3) {
    return(NULL)
  }

  present_sp <- sort(unique(df$taxon))

  color_vector <- color_palette[
    names(color_palette) %in% present_sp
  ]

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["taxon"]],
      y = .data[[traitname]],
      fill = .data[["taxon"]]
    )
  ) +
    ggplot2::geom_boxplot(
      lwd = 0.1,
      alpha = 0.5,
      outlier.shape = NA
    ) +
    {
      if (show_jitter)
        ggplot2::geom_jitter(
          ggplot2::aes(color = .data[["taxon"]]),
          width = 0.2,
          alpha = 0.4,
          size = 1.2
        )
    } +
    ggplot2::scale_fill_manual(values = color_vector) +
    ggplot2::scale_color_manual(values = color_vector) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none",
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        face = "italic"
      ),
      plot.margin = grid::unit(
        c(1, 1, 1, 1),
        "cm"
      )
    ) +
    ggplot2::ylab(
      paste0(clean_label(traitname), " (cm)")
    ) +
    ggplot2::xlab("")

  p
}

# ---- One plot per trait --------------------------------------------------

p_boxplots <- list()

for (traitname in trait_vars) {

df <- bp_data[
  ,
  c("taxon", traitname),
  drop = FALSE
]

if (sum(!is.na(df[[traitname]])) < min_n_total) {

  if (verbose) {
    message(
      "  Skipping ",
      traitname,
      ": too few observations."
    )
  }

  next
}

safe_name <- gsub(
  "[^A-Za-z0-9_]",
  "_",
  traitname
)

if (plot_type %in% c("boxplot", "both")) {

  p_box <- make_boxplot(
    df,
    traitname,
    color_palette,
    show_jitter
  )

  if (!is.null(p_box)) {

    if ("pdf" %in% file_formats) {
      cowplot::save_plot(
        file.path(
          output_dir,
          paste0("boxplot_", safe_name, ".pdf")
        ),
        p_box,
        ncol = 1,
        nrow = 1,
        base_height = 8.5,
        base_aspect_ratio = 1.3
      )
    }

    if ("jpeg" %in% file_formats) {
      cowplot::save_plot(
        file.path(
          output_dir,
          paste0("boxplot_", safe_name, ".jpeg")
        ),
        p_box,
        ncol = 1,
        nrow = 1,
        base_height = 8.5,
        base_aspect_ratio = 1.3
      )
    }
  }
}

if (plot_type %in% c("violin", "both")) {

  p_vio <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["taxon"]],
      y = .data[[traitname]],
      fill = .data[["taxon"]]
    )
  ) +
    ggplot2::geom_violin(
      alpha = 0.6,
      trim = FALSE,
      lwd = 0.3
    ) +
    ggplot2::geom_boxplot(
      width = 0.1,
      fill = "white",
      alpha = 0.8,
      outlier.shape = NA,
      lwd = 0.3
    ) +
    {
      if (show_jitter)
        ggplot2::geom_jitter(
          ggplot2::aes(color = .data[["taxon"]]),
          width = 0.1,
          alpha = 0.3,
          size = 1
        )
    } +
    ggplot2::scale_fill_manual(
      values = color_palette[
        names(color_palette) %in% sort(unique(df$taxon))
      ]
    ) +
    ggplot2::scale_color_manual(
      values = color_palette[
        names(color_palette) %in% sort(unique(df$taxon))
      ]
    ) +
    ggplot2::ylab(
      paste0(clean_label(traitname), " (cm)")
    ) +
    ggplot2::xlab("") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none",
      axis.text.x = ggplot2::element_text(
        angle = 45,
        vjust = 1,
        hjust = 1,
        face = "italic"
      )
    )

  if ("pdf" %in% file_formats) {
    cowplot::save_plot(
      file.path(
        output_dir,
        paste0("violin_", safe_name, ".pdf")
      ),
      p_vio,
      ncol = 1,
      nrow = 1,
      base_height = 8.5,
      base_aspect_ratio = 1.3
    )
  }

  if ("jpeg" %in% file_formats) {
    cowplot::save_plot(
      file.path(
        output_dir,
        paste0("violin_", safe_name, ".jpeg")
      ),
      p_vio,
      ncol = 1,
      nrow = 1,
      base_height = 8.5,
      base_aspect_ratio = 1.3
    )
  }
}

if (verbose) {
  message("  Done: ", traitname)
}
}

if (verbose) {
  message(
    "Boxplot analysis done. Outputs saved in: ",
    output_dir
  )
}
}
