#' Prepare a geographical/environmental matrix for ecogeographic analyses
#'
#' @description
#' Reads a structured Excel spreadsheet of georeferenced specimen
#' occurrences and builds a combined geographical and environmental
#' matrix for downstream ecogeographic and niche analyses in
#' \pkg{moRphoTaxa}.
#'
#' The function standardizes the taxon column, optionally filters taxa,
#' removes records with missing coordinates, and extracts WorldClim
#' bioclimatic variables and elevation at each occurrence point.
#' Optionally, SoilGrids edaphic variables can also be extracted.
#'
#' By default, the output retains only essential occurrence columns
#' (taxon, decimalLatitude, and decimalLongitude), together with
#' the extracted environmental variables. Original occurrence metadata
#' can be retained using \code{keep_original = TRUE}.
#'
#' Climate and, when requested, soil raster layers used for extraction
#' are kept and returned as attributes so they can be reused by
#' downstream mapping or niche-modelling functions without being
#' fetched again.
#'
#' @details
#' The input spreadsheet is expected to contain one row per occurrence
#' record, a column identifying the taxon of each record, and Darwin
#' Core-style \code{decimalLatitude} and \code{decimalLongitude}
#' coordinate columns.
#'
#' The taxon column is internally standardized to \code{"taxon"}.
#' Records can optionally be filtered using \code{species_selected}.
#' Records with missing taxon assignments or missing geographic
#' coordinates are removed.
#'
#' By default, only the essential occurrence columns are retained:
#' \code{taxon}, \code{decimalLatitude}, and
#' \code{decimalLongitude}. If an \code{"id"} column is present,
#' it is also retained. All other original spreadsheet columns are
#' excluded from the returned analysis matrix.
#'
#' When \code{keep_original = TRUE}, all original spreadsheet columns
#' are retained before the environmental variables are appended.
#'
#' Climate data are fetched from WorldClim
#' ([geodata::worldclim_global()]) at the resolution given by
#' \code{wc_res}: the 19 bioclimatic variables
#' (\code{bio1}-\code{bio19}) plus elevation.
#'
#' Elevation is extracted as a separate environmental variable named
#' \code{"elevation"}. It is not a twentieth bioclimatic variable
#' and is not named \code{"bio0"}.
#'
#' If the original spreadsheet contains an elevation column, it is
#' excluded from the default analysis matrix to avoid duplication
#' with the WorldClim-derived elevation variable.
#'
#' Edaphic (soil) data are extracted only when \code{edaphic = TRUE}.
#' Each variable in \code{edaphic_vars} is fetched from SoilGrids
#' ([geodata::soil_world()]) at the depth given by
#' \code{edaphic_depth}.
#'
#' A variable that fails to download is skipped with a warning rather
#' than stopping the whole function. Only successfully fetched
#' variables are added to \code{geodata} and recorded in the
#' \code{"env_cols"} attribute.
#'
#' When edaphic layers are extracted, the resulting soil raster stack
#' is resampled to match the climate raster's resolution before being
#' attached as an attribute.
#'
#' Occurrence points are converted once to a [terra::vect()] object
#' (WGS84, \code{"EPSG:4326"}) and reused for climate and edaphic
#' extraction. This object is also returned as the \code{"pts"}
#' attribute.

#' @param xlsx_path Character. Path to the `.xlsx` file containing the
#' occurrence data.
#'
#' @param sheet Character string or integer specifying the Excel worksheet to
#' read. Passed directly to [openxlsx::read.xlsx()]. Default is `1`.
#'
#' @param taxon_col Character string giving the name of the column containing
#' taxon assignments. This column is internally renamed to `"taxon"`. If
#' `NULL`, no column is explicitly renamed, and the spreadsheet must already
#' contain a `"taxon"` column.
#'
#' @param species_selected Character vector specifying taxa to retain in the
#' analysis. Values must match entries in the taxon column. If `NULL`
#' (default), all taxa are retained.
#'
#' @param wc_res Numeric. WorldClim resolution, in minutes of arc. One of
#' `0.5`, `2.5`, `5`, or `10`. Default is `2.5`.
#'
#' @param edaphic Logical. If `TRUE` (default), also extracts edaphic (soil)
#' variables at each occurrence point. If `FALSE`, only occurrence and
#' climate data are returned, and `edaphic_vars`/`edaphic_depth` are
#' ignored.
#'
#' @param edaphic_vars Character vector of SoilGrids variable codes to
#' extract when `edaphic = TRUE`. Default is `c("clay", "sand", "silt",
#' "soc", "phh2o", "cec", "nitrogen", "bdod", "cfvo", "ocd", "ocs")`.
#'
#' @param edaphic_depth Numeric. Soil depth, in centimeters, at which
#' edaphic variables are extracted when `edaphic = TRUE`. One of `5`, `15`,
#' `30`, `60`, `100`, or `200`. Default is `5`.
#'
#' @param data_path Character string giving the directory where downloaded
#' WorldClim and SoilGrids raster files are cached. Created if it does not
#' already exist. Default is `"variables_data/"`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' including the number of occurrence records retained, the climate and
#' edaphic variables added, and any variables that failed to download.
#'
#' @return
#' A data frame (`geodata`) with one row per retained occurrence record:
#'
#' \itemize{
#'   \item `taxon`: taxonomic assignment.
#'   \item `decimalLatitude`: geographic latitude.
#'   \item `decimalLongitude`: geographic longitude.
#'   \item `id`: specimen or occurrence identifier, if present in the
#'   original data.
#'   \item `bio1`-`bio19`: WorldClim bioclimatic variables.
#'   \item `elevation`: WorldClim elevation.
#'   \item Soil variables: when `edaphic = TRUE` and successfully fetched.
#' }
#'
#' Four attributes are attached:
#' \describe{
#'   \item{`env_cols`}{A named list with a `"clim"` element
#'   (character vector of climate column names, including
#'   `"elevation"`) and, when applicable, an `"edaphic"` element
#'   (character vector of successfully extracted edaphic columns).}
#'   \item{`clim_stack`}{The [terra::SpatRaster] of climate layers
#'   (`bio1`-`bio19` and `elevation`) used for extraction.}
#'   \item{`soil_stack`}{The [terra::SpatRaster] of edaphic layers used
#'   for extraction, resampled to the climate raster's resolution;
#'   `NULL` if `edaphic = FALSE` or no soil variable was successfully
#'   fetched.}
#'   \item{`pts`}{The [terra::SpatVector] of occurrence points
#'   (`"EPSG:4326"`) built from the retained records.}
#' }
#'
#' @seealso
#' [morph_matrix_setting()]
#'
#' @examples
#' \dontrun{
#' # Build a geographical/environmental matrix with climate only
#' geo <- geo_matrix_setting(
#'   xlsx_path = "occurrence_data.xlsx",
#'   taxon_col = "species",
#'   edaphic = FALSE
#' )
#'
#' # Occurrence + climate + edaphic matrix
#' geo <- geo_matrix_setting(
#'   xlsx_path = "occurrence_data.xlsx",
#'   taxon_col = "species",
#'   species_selected = c("Species A", "Species B"),
#'   wc_res = 2.5,
#'   edaphic = TRUE,
#'   edaphic_depth = 5
#' )
#'
#' # Occurrence + environmental data
#' geo
#'
#' # Which columns belong to which environmental block
#' attr(geo, "env_cols")
#'
#' # Reuse the cached climate raster and occurrence points downstream
#' attr(geo, "clim_stack")
#' attr(geo, "pts")
#' }
#'
#' @importFrom openxlsx read.xlsx
#' @importFrom stats complete.cases
#' @importFrom terra vect extract rast resample
#' @importFrom geodata worldclim_global soil_world
#'
#' @export

geo_matrix_setting <- function(xlsx_path,
                               sheet = 1,
                               taxon_col = NULL,
                               species_selected = NULL,
                               wc_res = 2.5,
                               edaphic = TRUE,
                               edaphic_vars = c("clay", "sand", "silt", "soc",
                                                "phh2o", "cec", "nitrogen",
                                                "bdod", "cfvo", "ocd", "ocs"),
                               edaphic_depth = 5,
                               data_path = "variables_data/",
                               verbose = TRUE) {

# ---- Validate input ---------------------------------------------------------

if (!is.character(xlsx_path) || length(xlsx_path) != 1L) {
  stop("`xlsx_path` must be a single character string.")
}

if (!is.numeric(wc_res) || length(wc_res) != 1L || !(wc_res %in% c(0.5, 2.5, 5, 10))) {
  stop("`wc_res` must be one of 0.5, 2.5, 5, or 10 (minutes of arc).")
}

if (!is.logical(edaphic) || length(edaphic) != 1L) {
  stop("`edaphic` must be TRUE or FALSE.")
}

if (edaphic) {
  if (!is.character(edaphic_vars) || length(edaphic_vars) < 1L) {
    stop("`edaphic_vars` must be a character vector with at least one variable code.")
  }
  if (!is.numeric(edaphic_depth) || length(edaphic_depth) != 1L ||
      !(edaphic_depth %in% c(5, 15, 30, 60, 100, 200))) {
    stop("`edaphic_depth` must be one of 5, 15, 30, 60, 100, or 200 (cm).")
  }
}

if (!is.character(data_path) || length(data_path) != 1L) {
  stop("`data_path` must be a single character string.")
}

if (!is.logical(verbose) || length(verbose) != 1L) {
  stop("`verbose` must be TRUE or FALSE.")
}

# ---- Read occurrence data ----------------------------------------------------

  raw_data <- openxlsx::read.xlsx(
    xlsxFile = xlsx_path,
    sheet = sheet
  )

  # Create stable specimen IDs before any filtering
  if (!"id" %in% names(raw_data)) {
    raw_data$id <- sprintf("ID%04d", seq_len(nrow(raw_data)))
  }

if (!is.null(taxon_col)) {
  tf <- names(raw_data) %in% taxon_col
  names(raw_data)[tf] <- "taxon"
}

if (!"taxon" %in% names(raw_data)) {
  stop(
    "No `taxon` column found. Supply `taxon_col`, or make sure the ",
    "spreadsheet already has a column named \"taxon\"."
  )
}

if (!all(c("decimalLatitude", "decimalLongitude") %in% names(raw_data))) {
  stop("`xlsx_path` must contain `decimalLatitude` and `decimalLongitude` columns.")
}

# ---- Filter occurrence records -----------------------------------------------

geodata <- raw_data

# Convert coordinates to numeric
geodata$decimalLatitude <- as.numeric(geodata$decimalLatitude)
geodata$decimalLongitude <- as.numeric(geodata$decimalLongitude)

# Filter selected taxa
if (!is.null(species_selected)) {
  geodata <- geodata[
    geodata$taxon %in% species_selected,
    ,
    drop = FALSE
  ]
}

# Remove records without taxon
geodata <- geodata[
  !is.na(geodata$taxon),
  ,
  drop = FALSE
]

# Remove records without coordinates
geodata <- geodata[
  stats::complete.cases(
    geodata[, c("decimalLongitude", "decimalLatitude")]
  ),
  ,
  drop = FALSE
]

# ---- Keep only essential occurrence columns ---------------------------------

keep_cols <- c(
  "id",
  "taxon",
  "decimalLatitude",
  "decimalLongitude"
)

# Keep only columns that actually exist
keep_cols <- intersect(keep_cols, names(geodata))

# Remove all other original columns, including original elevation
geodata <- geodata[
  ,
  keep_cols,
  drop = FALSE
]

if (nrow(geodata) == 0) {
  stop("No records remain after filtering. Check `xlsx_path` and `species_selected`.")
}

if (verbose) {
  message("Building geographical/environmental matrix...")
  message("  Occurrence records: ", nrow(geodata))
}

# ---- Occurrence points (shared by climate + edaphic extraction) --------------

pts <- terra::vect(
  geodata,
  geom = c("decimalLongitude", "decimalLatitude"),
  crs = "EPSG:4326"
)

# ---- Climate: WorldClim bioclimatic variables + elevation --------------------

dir.create(data_path, recursive = TRUE, showWarnings = FALSE)

if (verbose) {
  message("  Fetching WorldClim rasters (res = ", wc_res, ")...")
}

bio_stack <- geodata::worldclim_global(var = "bio",  res = wc_res, path = data_path)
elev_stack <- geodata::worldclim_global(var = "elev", res = wc_res, path = data_path)

names(bio_stack) <- gsub("^wc2\\.1_[0-9.]+m_", "", names(bio_stack))

# `elev_stack` is the elevation layer, NOT a 20th bioclimatic variable, so it
# is named "elevation" here rather than folded into the "bio*" naming scheme
# (avoiding the previous, incorrect "bio_0" label).
names(elev_stack) <- "elevation"

clim_stack <- c(bio_stack, elev_stack)

clim_vals <- terra::extract(clim_stack, pts)[, -1, drop = FALSE]
geodata <- cbind(geodata, clim_vals)

if (verbose) {
  message(
    "  Climate variables added: ", ncol(clim_vals),
    " (", paste(names(clim_vals), collapse = ", "), ")"
  )
}

env_cols <- list(clim = names(clim_vals))
soil_stack <- NULL

# ---- Edaphic (soil) variables, optional ---------------------------------------

if (edaphic) {

  if (verbose) {
    message("Extracting edaphic (soil) data...")
  }

  edaphic_path <- file.path(data_path, "soil")
  dir.create(edaphic_path, recursive = TRUE, showWarnings = FALSE)

  soil_list <- list()
  soil_rasters <- list()

  for (v in edaphic_vars) {

    if (verbose) {
      message("  Fetching: ", v, " (depth = ", edaphic_depth, "cm)")
    }

    r <- tryCatch(
      geodata::soil_world(
        var = v, depth = edaphic_depth,
        stat = "mean", path = edaphic_path
      ),
      error = function(e) {
        warning("Could not fetch '", v, "': ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(r)) next

    names(r) <- v
    soil_rasters[[v]] <- r

    vals <- terra::extract(r, pts)[, -1, drop = FALSE]
    names(vals) <- v
    soil_list[[v]] <- vals
  }

  if (length(soil_list) == 0) {
    warning("No edaphic variables were successfully extracted; edaphic data were skipped.")
  } else {
    soil_vals <- do.call(cbind, soil_list)
    geodata <- cbind(geodata, soil_vals)
    env_cols$edaphic <- names(soil_vals)

    soil_stack <- terra::rast(soil_rasters)
    soil_stack <- terra::resample(soil_stack, clim_stack[[1]])

    if (verbose) {
      message(
        "Done. Edaphic variables added: ",
        paste(env_cols$edaphic, collapse = ", ")
      )
    }
  }
}

# ---- Attach attributes ---------------------------------------------------------

attr(geodata, "env_cols") <- env_cols
attr(geodata, "clim_stack") <- clim_stack
attr(geodata, "soil_stack") <- soil_stack
attr(geodata, "pts") <- pts

if (verbose) {
  message(
    "\nDone. `geodata` built with ", nrow(geodata), " records and ",
    ncol(geodata), " columns."
  )
}

rownames(geodata) <- geodata$id
geodata$id <- NULL

return(geodata)
}
