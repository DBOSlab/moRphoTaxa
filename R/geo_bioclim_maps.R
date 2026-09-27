#' Map occurrences over a bioclimatic layer and plot them in climate space
#'
#' @description
#' Builds a two-panel bioclimatic figure from the geographical/environmental
#' matrix produced by [geo_matrix_setting()]: (a) a distribution map in which
#' occurrence points are drawn over a WorldClim raster background, optionally
#' with sub-national borders, rivers and a locator inset; and (b) a scatter
#' plot of the occurrences in the climate space defined by two WorldClim
#' variables.
#'
#' The map extent, the raster cropping and the locator inset are all
#' determined automatically from the geographical distribution of the
#' occurrence records. The function reuses the climate raster stored in the
#' `"clim_stack"` attribute of the `geodata` object, so no WorldClim data are
#' downloaded again. Figures are optionally saved as PDF and JPEG, and the
#' \pkg{ggplot2} objects are returned so they can be further customised.
#'
#' @details
#' \strong{Input.} `geodata` must be the data frame returned by
#' [geo_matrix_setting()]. It must contain the columns `taxon`,
#' `decimalLatitude` and `decimalLongitude`, the climate columns listed in
#' `attr(geodata, "env_cols")$clim`, and the `"clim_stack"` attribute
#' (a [terra::SpatRaster]).
#'
#' \strong{Choosing climate variables.} `bio_map_layer`, `bio_x` and `bio_y`
#' accept an integer from `1` to `19` (WorldClim bioclimatic variables
#' BIO1-BIO19) or the character string `"elevation"`. Elevation is a separate,
#' non-bioclimatic layer: it is \emph{not} "BIO0". Bioclimatic layers are
#' matched against the layer names present in the data (both the `"bio_1"`
#' and `"bio1"` naming styles are recognised). The reference for the
#' variables is:
#'
#' \preformatted{
#' BIO1 = Annual Mean Temperature        BIO11 = Mean Temp Coldest Quarter
#' BIO2 = Mean Diurnal Range             BIO12 = Annual Precipitation
#' BIO3 = Isothermality                  BIO13 = Precipitation Wettest Month
#' BIO4 = Temperature Seasonality        BIO14 = Precipitation Driest Month
#' BIO5 = Max Temp Warmest Month         BIO15 = Precipitation Seasonality
#' BIO6 = Min Temp Coldest Month         BIO16 = Precipitation Wettest Quarter
#' BIO7 = Temperature Annual Range       BIO17 = Precipitation Driest Quarter
#' BIO8 = Mean Temp Wettest Quarter      BIO18 = Precipitation Warmest Quarter
#' BIO9 = Mean Temp Driest Quarter       BIO19 = Precipitation Coldest Quarter
#' BIO10 = Mean Temp Warmest Quarter     elevation = Elevation (m)
#' }
#'
#' \strong{Automatic map extent.} The visible area of the map is the bounding
#' box of the occurrence points plus a margin of 15\% of the longitudinal and
#' latitudinal range of the data (at least 1 degree). The margin on the
#' right-hand side is 50\% larger, leaving room for the legend, scale bar and
#' north arrow. The climate raster and rivers are cropped to this area plus an
#' additional 25\% buffer, so the raster always fills the map panel. The
#' locator inset (when `inset = TRUE`) shows the mapped area, marked by a red
#' rectangle, in its wider geographical context: it extends the map area by
#' its own width and height in each direction (at least 15 degrees). All
#' extents are clamped to valid longitudes and latitudes. Because the extent
#' is computed from the range of the coordinates, datasets that cross the
#' antimeridian (180 degrees) are not handled specially.
#'
#' \strong{Point shapes.} Each taxon is drawn with a distinct point shape,
#' assigned automatically (in alphabetical order of taxa) from a pool of 24
#' shapes. If there are more than 24 taxa, supply your own `pt_shapes`.
#'
#' \strong{Optional layers.} Country and state boundaries come from
#' \pkg{rnaturalearth} (`ne_countries()` at medium scale; sub-national borders
#' via `ne_states()` when `show_country_states` is not `NULL`, which requires
#' the \pkg{rnaturalearthhires} package). Rivers are downloaded with
#' `ne_download()` and therefore require an internet connection.
#'
#' \strong{Output.} When `save = TRUE`, figures are written to a dated
#' sub-folder of `"Figs.bioclimaps"` in the working directory (for example
#' `"Figs.bioclimaps/18Sep2026"`), created if needed. Three figures are saved
#' in the formats selected by `file_formats`: `MAP_RASTER` (panel a),
#' `PLOT_RASTER` (panel b) and `MAP_RASTER_PLOT` (both panels combined and
#' labelled "a)" and "b)"). The folder path is returned in the
#' `output_dir` element of the result.
#'
#' @param geodata Data frame returned by [geo_matrix_setting()], carrying the
#' `"env_cols"` and `"clim_stack"` attributes.
#'
#' @param bio_map_layer Integer from `1` to `19`, or `"elevation"`. Climate
#' variable shown as the raster background of the map. Default is `12`
#' (annual precipitation).
#'
#' @param bio_x,bio_y Integer from `1` to `19`, or `"elevation"`. Climate
#' variables plotted on the x and y axes of the scatter plot. Default is
#' `1` (annual mean temperature) for `bio_x` and `12` (annual precipitation)
#' for `bio_y`.
#'
#' @param show_country_states Character string with a country name (e.g.
#' `"brazil"`) whose sub-national borders are drawn on the map, passed to
#' `rnaturalearth::ne_states()`. If `NULL` (default), no sub-national borders
#' are drawn.
#'
#' @param show_rivers Logical. If `TRUE`, river centre lines are downloaded
#' and drawn on the map. Default is `FALSE`.
#'
#' @param river_scale Numeric. Natural Earth resolution for rivers, one of
#' `10` (fine detail), `50` or `110` (coarser and faster). Only used when
#' `show_rivers = TRUE`. Default is `10`.
#'
#' @param inset Logical. If `TRUE` (default), a small locator map is added to
#' the lower-left corner of the main map, with a red rectangle marking the
#' mapped region.
#'
#' @param pt_shapes Named numeric vector mapping each taxon to a point shape
#' (see [graphics::points()]), or `NULL` (default) to assign shapes
#' automatically. When supplied, it must contain an entry for every taxon in
#' `geodata`.
#'
#' @param save Logical. If `TRUE` (default), figures are written to disk. If
#' `FALSE`, nothing is saved and the plots are only returned.
#'
#' @param file_formats Character vector specifying the file formats used to
#' save figures. One or more of `"pdf"` and `"jpg"`. Default is
#' `c("pdf", "jpg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' including the variables plotted and the output folder.
#'
#' @return
#' Invisibly, a named list with the elements:
#' \describe{
#'   \item{`map`}{The main distribution map, including the inset when
#'   `inset = TRUE` (a \pkg{cowplot} composite).}
#'   \item{`scatter`}{The climate scatter plot (a [ggplot2::ggplot()]
#'   object).}
#'   \item{`combined`}{The two-panel figure, labelled "a)" and "b)".}
#'   \item{`inset`}{The locator inset ([ggplot2::ggplot()]), or `NULL` if
#'   `inset = FALSE`.}
#'   \item{`output_dir`}{Path of the folder where figures were saved, or
#'   `NULL` if `save = FALSE`.}
#' }
#'
#' @seealso
#' [geo_matrix_setting()]
#'
#' @examples
#' \dontrun{
#' geo <- geo_matrix_setting(
#'   xlsx_path = "occurrence_data.xlsx",
#'   taxon_col = "species",
#'   edaphic = FALSE
#' )
#'
#' # Default figure: precipitation background, temperature x precipitation
#' figs <- geo_bioclim_maps(geo)
#'
#' # Brazil with state borders and rivers, precipitation map,
#' # temperature seasonality x precipitation seasonality scatter plot
#' figs <- geo_bioclim_maps(
#'   geo,
#'   bio_map_layer = 12,
#'   bio_x = 4,
#'   bio_y = 15,
#'   show_country_states = "brazil",
#'   show_rivers = TRUE
#' )
#'
#' # Elevation as the map background, without saving files
#' figs <- geo_bioclim_maps(geo, bio_map_layer = "elevation", save = FALSE)
#'
#' # Further customise a returned plot
#' figs$scatter + ggplot2::theme(legend.position = "bottom")
#'
#' # Save PDF only
#' figs <- geo_bioclim_maps(
#'   geo,
#'   file_formats = "pdf"
#' )
#'
#' # Save JPEG only
#' figs <- geo_bioclim_maps(
#'   geo,
#'   file_formats = "jpg"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_sf geom_tile geom_point geom_rug
#' @importFrom ggplot2 annotate coord_sf labs theme theme_bw element_blank
#' @importFrom ggplot2 element_rect element_text scale_fill_viridis_c
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom cowplot ggdraw draw_plot plot_grid save_plot
#' @importFrom ggspatial annotation_scale annotation_north_arrow
#' @importFrom ggspatial north_arrow_orienteering
#' @importFrom rnaturalearth ne_countries ne_states ne_download
#' @importFrom sf st_crop
#' @importFrom terra ext crop
#' @importFrom grid unit
#' @importFrom rlang .data
#' @importFrom stats setNames
#'
#' @export

geo_bioclim_maps <- function(geodata,
                             bio_map_layer = 12,
                             bio_x = 1,
                             bio_y = 12,
                             show_country_states = NULL,
                             show_rivers = FALSE,
                             river_scale = 10,
                             inset = TRUE,
                             pt_shapes = NULL,
                             save = TRUE,
                             file_formats = c("pdf", "jpg"),
                             verbose = TRUE) {

# ---- Validate input ---------------------------------------------------------

if (!is.data.frame(geodata)) {
  stop("`geodata` must be the data frame returned by `geo_matrix_setting()`.")
}

if (!all(c("taxon", "decimalLatitude", "decimalLongitude") %in% names(geodata))) {
  stop("`geodata` must contain `taxon`, `decimalLatitude` and `decimalLongitude` columns.")
}

clim_stack <- attr(geodata, "clim_stack")

if (is.null(clim_stack) || !inherits(clim_stack, "SpatRaster")) {
  stop(
    "`geodata` has no valid \"clim_stack\" attribute. ",
    "Build it with `geo_matrix_setting()`."
  )
}

if (!is.logical(show_rivers) || length(show_rivers) != 1L) {
  stop("`show_rivers` must be TRUE or FALSE.")
}

if (show_rivers && !(length(river_scale) == 1L && river_scale %in% c(10, 50, 110))) {
  stop("`river_scale` must be one of 10, 50, or 110.")
}

if (!is.null(show_country_states) &&
    (!is.character(show_country_states) || length(show_country_states) != 1L)) {
  stop("`show_country_states` must be a single country name or NULL.")
}

if (!is.logical(inset) || length(inset) != 1L) {
  stop("`inset` must be TRUE or FALSE.")
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

if (!is.logical(verbose) || length(verbose) != 1L) {
  stop("`verbose` must be TRUE or FALSE.")
}

# ---- Internal settings (not user arguments) -----------------------------------

base_dir <- "Figs.bioclimaps" # parent folder for saved figures
pt_size <- 4 # point aesthetics
pt_stroke <- 0.1
pt_fill <- "gray30"
pt_colour <- "gray30"
pt_alpha <- 0.8
palette <- "inferno" # viridis palette for the raster background

# ---- Resolve climate variables -----------------------------------------------

env_cols <- attr(geodata, "env_cols")
clim_cols <- if (!is.null(env_cols$clim)) env_cols$clim else names(clim_stack)
clim_cols <- intersect(clim_cols, names(geodata))

col_map <- .resolve_clim_var(bio_map_layer, names(clim_stack), "bio_map_layer")
col_x <- .resolve_clim_var(bio_x, clim_cols, "bio_x")
col_y <- .resolve_clim_var(bio_y, clim_cols, "bio_y")

lab_map <- .clim_labels(col_map)
lab_x <- .clim_labels(col_x)
lab_y <- .clim_labels(col_y)

if (verbose) {
  message("Building bioclimatic maps...")
  message("  Map background: ", lab_map$axis)
  message("  Scatter plot:   ", lab_x$axis, "  x  ", lab_y$axis)
}

# ---- Point shapes by taxon ---------------------------------------------------

taxa_sorted <- sort(unique(geodata$taxon))

if (is.null(pt_shapes)) {
  shape_pool <- c(21, 22, 24, 23, 25, 8, 4, 3, 7, 9, 10, 11, 12, 13, 14,
                  15, 16, 17, 18, 0, 1, 2, 5, 6)
  if (length(taxa_sorted) > length(shape_pool)) {
    stop(
      "There are ", length(taxa_sorted), " taxa but only ", length(shape_pool),
      " automatic point shapes. Supply `pt_shapes`."
    )
  }
  pt_shapes <- stats::setNames(shape_pool[seq_along(taxa_sorted)], taxa_sorted)
} else {
  if (is.null(names(pt_shapes))) {
    stop("`pt_shapes` must be a named vector (taxon name = shape).")
  }
  missing_taxa <- setdiff(taxa_sorted, names(pt_shapes))
  if (length(missing_taxa) > 0) {
    stop("`pt_shapes` is missing an entry for: ", paste(missing_taxa, collapse = ", "))
  }
}

# ---- Automatic extents from the occurrence distribution ------------------------

lon_rng <- range(geodata$decimalLongitude, na.rm = TRUE)
lat_rng <- range(geodata$decimalLatitude,  na.rm = TRUE)

ext_auto <- .auto_extent(lon_rng, lat_rng)
map_lim <- ext_auto$map
crop_lim <- ext_auto$crop
inset_lim <- ext_auto$inset

# ---- Spatial reference layers --------------------------------------------------

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

states <- if (!is.null(show_country_states)) {
  rnaturalearth::ne_states(country = show_country_states, returnclass = "sf")
} else {
  NULL
}

rivers <- if (show_rivers) {
  if (verbose) message("  Downloading rivers (scale = ", river_scale, ")...")
  rnaturalearth::ne_download(
    scale = river_scale, type = "rivers_lake_centerlines",
    category = "physical", returnclass = "sf"
  )
} else {
  NULL
}

# ---- Inset: locator map --------------------------------------------------------

neo <- NULL

if (inset) {
  neo <- ggplot2::ggplot(data = world) +
    ggplot2::geom_sf(fill = "gray95", color = "gray70", linewidth = 0.05) +
    ggplot2::coord_sf(
      xlim = c(inset_lim[["xmin"]], inset_lim[["xmax"]]),
      ylim = c(inset_lim[["ymin"]], inset_lim[["ymax"]]),
      expand = FALSE
    ) +
    ggplot2::annotate(
      "rect",
      xmin = map_lim[["xmin"]], xmax = map_lim[["xmax"]],
      ymin = map_lim[["ymin"]], ymax = map_lim[["ymax"]],
      fill = NA, colour = "red", linewidth = 0.5
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}

# ---- Crop raster and build tile data frame --------------------------------------

ext_crop <- terra::ext(
  crop_lim[["xmin"]], crop_lim[["xmax"]],
  crop_lim[["ymin"]], crop_lim[["ymax"]]
)

rivers_crop <- if (!is.null(rivers)) {
  suppressWarnings(
    sf::st_crop(
      rivers,
      xmin = crop_lim[["xmin"]], xmax = crop_lim[["xmax"]],
      ymin = crop_lim[["ymin"]], ymax = crop_lim[["ymax"]]
    )
  )
} else {
  NULL
}

raster_bio <- terra::crop(clim_stack[[col_map]], ext_crop)
pointsdf_bio <- as.data.frame(raster_bio, xy = TRUE, na.rm = TRUE)
names(pointsdf_bio)[3] <- "bio_map"

# ---- Main raster map ----------------------------------------------------------------

map_raster <- ggplot2::ggplot(data = world) +
  ggplot2::geom_tile(
    data = pointsdf_bio,
    ggplot2::aes(x = .data$x, y = .data$y, fill = .data$bio_map),
    alpha = 0.9
  ) +
  ggplot2::scale_fill_viridis_c(
    name = lab_map$short, option = palette, direction = -1
  ) +
  ggplot2::geom_sf(colour = "gray70", fill = NA, linewidth = 0.2) +
  { if (!is.null(states))
    ggplot2::geom_sf(data = states, colour = "gray70", fill = NA,
                     linewidth = 0.2, alpha = 0.4) } +
  { if (!is.null(rivers_crop))
    ggplot2::geom_sf(data = rivers_crop, colour = "deepskyblue3",
                     linewidth = 0.25, alpha = 0.8) } +
  ggplot2::coord_sf(
    xlim = c(map_lim[["xmin"]], map_lim[["xmax"]]),
    ylim = c(map_lim[["ymin"]], map_lim[["ymax"]]),
    expand = FALSE
  ) +
  ggplot2::geom_point(
    data = geodata,
    ggplot2::aes(x = .data$decimalLongitude, y = .data$decimalLatitude,
                 shape = .data$taxon),
    alpha = pt_alpha, size = pt_size, stroke = pt_stroke,
    fill = pt_fill, colour = pt_colour
  ) +
  ggplot2::scale_shape_manual(values = pt_shapes) +
  ggspatial::annotation_scale(location = "tr", width_hint = 0.2) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = grid::unit(0.75, "in"), pad_y = grid::unit(0.5, "in"),
    style = ggspatial::north_arrow_orienteering
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    legend.position = "inside",
    legend.position.inside = c(1, 0.4),
    legend.key = ggplot2::element_rect(fill = NA, colour = NA),
    legend.text = ggplot2::element_text(size = 8, face = "italic"),
    legend.title = ggplot2::element_blank(),
    plot.margin = grid::unit(c(1, 1, 1, 1), "cm")
  )

fullmap_raster <- if (inset) {
  cowplot::ggdraw() +
    cowplot::draw_plot(map_raster) +
    cowplot::draw_plot(neo, x = 0.05, y = 0.05, width = 0.25, height = 0.25)
} else {
  map_raster
}

# ---- Climate scatter plot -------------------------------------------------------------

pt_raster <- ggplot2::ggplot(
  geodata,
  ggplot2::aes(x = .data[[col_x]], y = .data[[col_y]])
) +
  ggplot2::geom_rug(position = ggplot2::position_jitter(seed = 123), linewidth = 0.1,
                    colour = "red", alpha = 0.6) +
  ggplot2::geom_point(
    ggplot2::aes(shape = .data$taxon),
    alpha = pt_alpha, size = pt_size, stroke = pt_stroke,
    fill = pt_fill, colour = pt_colour
  ) +
  ggplot2::scale_shape_manual(values = pt_shapes) +
  ggplot2::labs(
    x = lab_x$axis, y = lab_y$axis,
    title = paste0(nrow(geodata), " specimens \u2014 WorldClim climate models")
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.8),
    legend.key = ggplot2::element_rect(fill = NA, colour = NA),
    legend.text = ggplot2::element_text(size = 10, face = "italic"),
    legend.title = ggplot2::element_blank(),
    plot.margin = grid::unit(c(1, 1, 1, 1), "cm")
  )

# ---- Combine ---------------------------------------------------------------------------

fig_raster <- cowplot::plot_grid(
  fullmap_raster, pt_raster,
  labels = c("a)", "b)"), label_size = 25,
  ncol = 2, nrow = 1, align = "hv"
)

# ---- Save ------------------------------------------------------------------------------

folder_name <- NULL

if (save) {

  folder_name <- file.path(base_dir, format(Sys.time(), "%d%b%Y"))
  dir.create(folder_name, recursive = TRUE, showWarnings = FALSE)

  save_outputs <- function(name, plot, ncol = 1, nrow = 1,
                           base_height = 8.5, base_aspect_ratio = 1.3) {
    for (ext in paste0(".", file_formats)) {
      cowplot::save_plot(
        file.path(folder_name, paste0(name, ext)), plot,
        ncol = ncol, nrow = nrow, base_height = base_height,
        base_aspect_ratio = base_aspect_ratio, base_width = NULL
      )
    }
  }

  save_outputs("MAP_RASTER",      fullmap_raster)
  save_outputs("PLOT_RASTER",     pt_raster)
  save_outputs("MAP_RASTER_PLOT", fig_raster,
  ncol = 2, nrow = 1, base_height = 10, base_aspect_ratio = 1.0)

  if (verbose) {
    message("Done. Outputs saved in: ", folder_name)
  }
} else if (verbose) {
  message("Done. Figures were not saved (`save = FALSE`).")
}

invisible(list(
  map = fullmap_raster,
  scatter = pt_raster,
  combined = fig_raster,
  inset = neo,
  output_dir = folder_name
))
}

#' Compute automatic map, raster-crop and inset extents
#'
#' @description
#' Internal helper deriving all map extents from the range of the occurrence
#' coordinates. The map area is the bounding box plus a 15\% margin (at least
#' 1 degree; 50\% larger on the right to leave room for the legend and north
#' arrow). The raster-crop area adds a further 25\% buffer around the map area.
#' The inset area extends the map area by its own width and height (at least
#' 15 degrees) in each direction. All limits are clamped to valid coordinates.
#'
#' @param lon_rng Numeric vector of length 2. Range of longitudes.
#' @param lat_rng Numeric vector of length 2. Range of latitudes.
#'
#' @return A list with elements `map`, `crop` and `inset`, each a named
#' numeric vector with `xmin`, `xmax`, `ymin` and `ymax`.
#'
#' @keywords internal
#'
#' @noRd
.auto_extent <- function(lon_rng, lat_rng) {

clamp <- function(e) {
  c(xmin = max(e[["xmin"]], -180), xmax = min(e[["xmax"]], 180),
    ymin = max(e[["ymin"]],  -90), ymax = min(e[["ymax"]],  90))
}

pad_x <- max(0.15 * diff(lon_rng), 1)
pad_y <- max(0.15 * diff(lat_rng), 1)

map <- clamp(c(
  xmin = lon_rng[1] - pad_x,
  xmax = lon_rng[2] + 1.5 * pad_x,
  ymin = lat_rng[1] - pad_y,
  ymax = lat_rng[2] + pad_y
))

w <- map[["xmax"]] - map[["xmin"]]
h <- map[["ymax"]] - map[["ymin"]]

crop <- clamp(c(
  xmin = map[["xmin"]] - 0.25 * w,
  xmax = map[["xmax"]] + 0.25 * w,
  ymin = map[["ymin"]] - 0.25 * h,
  ymax = map[["ymax"]] + 0.25 * h
))

inset_x <- max(w, 15)
inset_y <- max(h, 15)

inset <- clamp(c(
  xmin = map[["xmin"]] - inset_x,
  xmax = map[["xmax"]] + inset_x,
  ymin = map[["ymin"]] - inset_y,
  ymax = map[["ymax"]] + inset_y
))

list(map = map, crop = crop, inset = inset)
}

#' Resolve a climate variable to a column/layer name
#'
#' @description
#' Internal helper translating a user-supplied climate variable (an integer
#' 1-19 for BIO1-BIO19, `"elevation"`, or an exact layer name) into the name
#' actually present in the data or raster, accepting both the `"bio_N"` and
#' `"bioN"` naming styles.
#'
#' @param x Integer, or character string.
#' @param available Character vector of available column/layer names.
#' @param arg Character. Name of the argument, used in the error message.
#'
#' @return A single character string, present in `available`.
#'
#' @keywords internal
#'
#' @noRd
.resolve_clim_var <- function(x, available, arg) {

if (length(x) != 1L || is.na(x)) {
  stop("`", arg, "` must be a single value.")
}

if (is.numeric(x)) {
  if (!(x %in% 1:19)) {
    stop(
      "`", arg, "` must be an integer from 1 to 19 or \"elevation\". ",
      "Elevation is not a bioclimatic variable (there is no BIO0)."
    )
  }
  candidates <- c(paste0("bio_", x), paste0("bio", x))
} else {
  candidates <- as.character(x)
}

hit <- candidates[candidates %in% available]

if (length(hit) == 0L) {
  stop(
    "`", arg, "` (", x, ") was not found. Available variables: ",
    paste(available, collapse = ", "), "."
  )
}

hit[1]
}

#' Build human-readable labels for a climate variable
#'
#' @description
#' Internal helper returning a short label (used for the map legend) and an
#' axis label (used for the scatter plot) for a WorldClim variable.
#'
#' @param col Character. Column/layer name, e.g. `"bio_12"` or `"elevation"`.
#'
#' @return A list with elements `short` and `axis`.
#'
#' @keywords internal
#'
#' @noRd
.clim_labels <- function(col) {

if (col == "elevation") {
  return(list(short = "Elevation (m)", axis = "Elevation (m)"))
}

lookup <- c(
  "1" = "Annual Mean Temperature",
  "2" = "Mean Diurnal Range",
  "3" = "Isothermality",
  "4" = "Temperature Seasonality",
  "5" = "Max Temp Warmest Month",
  "6" = "Min Temp Coldest Month",
  "7" = "Temperature Annual Range",
  "8" = "Mean Temp Wettest Quarter",
  "9" = "Mean Temp Driest Quarter",
  "10" = "Mean Temp Warmest Quarter",
  "11" = "Mean Temp Coldest Quarter",
  "12" = "Annual Precipitation (mm)",
  "13" = "Precipitation Wettest Month (mm)",
  "14" = "Precipitation Driest Month (mm)",
  "15" = "Precipitation Seasonality \u2013 CV",
  "16" = "Precipitation Wettest Quarter (mm)",
  "17" = "Precipitation Driest Quarter (mm)",
  "18" = "Precipitation Warmest Quarter (mm)",
  "19" = "Precipitation Coldest Quarter (mm)"
)

n <- gsub("\\D", "", col)
short <- unname(lookup[n])

if (is.na(short)) {
  return(list(short = col, axis = col))
}

list(short = short, axis = paste0("BIO", n, " \u2013 ", short))
}
