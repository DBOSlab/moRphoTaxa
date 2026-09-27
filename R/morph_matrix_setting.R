#' Prepare a morphological matrix for downstream analyses
#'
#' @description
#' Reads a structured Excel spreadsheet containing morphological measurements
#' and prepares the data for downstream morphometric analyses in
#' \pkg{moRphoTaxa}.
#'
#' The function standardizes the taxon column, optionally filters taxa,
#' converts morphological variables to numeric format, removes empty rows and
#' columns, assigns unique specimen identifiers, stores the taxonomic identity
#' associated with each retained specimen, organizes traits into vegetative,
#' floral, and fruit blocks, and defines a colour for each retained taxon to be
#' reused consistently by downstream plotting functions.
#'
#' @details
#' The input spreadsheet is expected to contain one row per specimen and one
#' column identifying the taxon of each specimen. Morphological trait columns
#' should occur after the descriptive or specimen metadata columns.
#'
#' Trait columns may follow the \pkg{moRphoTaxa} naming convention:
#'
#' \preformatted{
#' full_trait_name/TRAITcode
#' }
#'
#' For example:
#'
#' \preformatted{
#' petiole_length/PETIlng
#' leaf_length/LEAFlng
#' inflorescence_length/INFLlng
#' }
#'
#' The `trait_name_type` argument controls how these compound names are handled:
#'
#' \itemize{
#'   \item `"full"` retains only the full descriptive name, such as
#'   `"petiole_length"`;
#'   \item `"code"` retains only the abbreviated trait code, such as
#'   `"PETIlng"`;
#'   \item `NULL` leaves the original column names unchanged.
#' }
#'
#' The morphological matrix begins at the column specified by `veg_first`.
#' All columns from this point to the end of the spreadsheet are treated as
#' morphological variables.
#'
#' Three optional arguments define the first trait in major morphological
#' blocks:
#'
#' \itemize{
#'   \item `veg_first`: first vegetative trait;
#'   \item `flo_first`: first floral trait;
#'   \item `fru_first`: first fruit trait.
#' }
#'
#' These positions are used to create lists of trait names corresponding to
#' the vegetative (`"veg"`), floral (`"flo"`), and fruit (`"fru"`) blocks.
#'
#' Morphological values are converted to numeric format. During this process,
#' leading and trailing whitespace is removed, commas used as decimal
#' separators are converted to periods, empty strings are converted to `NA`,
#' and infinite values are replaced with `NA`.
#'
#' Rows and columns containing no morphological observations after conversion
#' are removed.
#'
#' Each retained specimen receives a unique identifier in the form
#' `"ID0001"`, `"ID0002"`, and so forth. These identifiers are used as row
#' names in the returned morphometric matrix.
#'
#' A colour is assigned to each taxon retained in the final matrix (i.e. after
#' all filtering steps). If `colors = NULL` (default), colours are drawn
#' automatically from the `viridis` palette. Alternatively, the user may
#' supply a named character vector (taxon name = colour) via `colors` to
#' define custom colours. This taxon-colour mapping is stored as the
#' `"taxon_colors"` attribute of the returned `analysis_data` object.
#'
#' @param xlsx_path Character. Path to the `.xlsx` file containing the
#' morphological data.
#'
#' @param sheet Character string or integer specifying the Excel worksheet to
#' read. Passed directly to [readxl::read_excel()]. Default is `1`.
#'
#' @param taxon_col Character string giving the name of the column containing
#' taxon or morphological-group assignments. This column is internally renamed
#' to `"taxon"`. If `NULL`, no column is explicitly renamed.
#'
#' @param species_selected Character vector specifying taxa to retain in the
#' analysis. Values must match entries in the taxon column. If `NULL`
#' (default), all taxa are retained.
#'
#' @param trait_name_type Character string specifying how morphological trait
#' column names should be represented. One of `"full"`, `"code"`, or `NULL`.
#' `"full"` retains the descriptive trait name, `"code"` retains the abbreviated
#' trait code, and `NULL` leaves column names unchanged.
#'
#' @param veg_first Character string giving the name of the first vegetative
#' trait column. This also defines the first column included in the
#' morphometric analysis matrix.
#'
#' @param flo_first Character string giving the name of the first floral trait
#' column. Used to define the beginning of the floral trait block.
#'
#' @param fru_first Character string giving the name of the first fruit trait
#' column. Used to define the beginning of the fruit trait block.
#'
#' @param colors Named character vector of colours to assign to each taxon
#' (e.g. `c("Species A" = "#1B9E77", "Species B" = "#D95F02")`), or `NULL`
#' (default) to assign colours automatically from the `viridis` palette. When
#' supplied, `colors` must contain an entry for every taxon retained in the
#' final `analysis_data`; an error is raised otherwise. The resulting mapping
#' is stored as the `"taxon_colors"` attribute of the returned object.
#'
#' @return
#' A numeric data frame (`analysis_data`) with two additional attributes:
#' \describe{
#'   \item{`base_cols`}{A named list of trait columns per block ("veg","flo","fru").}
#'   \item{`taxon_colors`}{A named character vector mapping each taxon to a colour.}
#' }
#'
#' @seealso
#' [morph_add_traits()]
#'
#' @examples
#' \dontrun{
#' # Prepare a matrix using abbreviated trait codes
#' morph <- morph_matrix_setting(
#'   xlsx_path = "morphological_data.xlsx",
#'   sheet = 1,
#'   taxon_col = "species",
#'   trait_name_type = "code",
#'   veg_first = "petiole_length/PETIlng",
#'   flo_first = "inflorescence_length/INFLlng",
#'   fru_first = "fruit_stipe_length/FRSTlng"
#' )
#'
#' # Morphological matrix
#' morph$analysis_data
#'
#' # Trait blocks
#' morph$base_cols
#'
#' # Restrict the analysis to selected taxa
#' morph_selected <- morph_matrix_setting(
#'   xlsx_path = "morphological_data.xlsx",
#'   taxon_col = "species",
#'   species_selected = c("Species A", "Species B"),
#'   trait_name_type = "full",
#'   veg_first = "petiole_length/PETIlng",
#'   flo_first = "inflorescence_length/INFLlng",
#'   fru_first = "fruit_stipe_length/FRSTlng"
#' )
#' }
#'
#' @importFrom openxlsx read.xlsx
#' @importFrom stats setNames
#' @importFrom viridis viridis
#'
#' @export

morph_matrix_setting <- function(xlsx_path,
                                 sheet = 1,
                                 taxon_col = NULL,
                                 species_selected = NULL,
                                 trait_name_type = NULL,
                                 veg_first = NULL,
                                 flo_first = NULL,
                                 fru_first = NULL,
                                 colors = NULL) {


# ---- Resolve trait column naming (full name vs. code) ----------------------

# Read Excel structured spreadsheet
raw_data <- openxlsx::read.xlsx(xlsxFile = xlsx_path, sheet = sheet)

tf <- names(raw_data) %in% taxon_col
names(raw_data)[tf] <- "taxon"

names(raw_data) <- .pick_name(names(raw_data), trait_name_type)
veg_first <- .pick_name(veg_first, trait_name_type)
flo_first <- .pick_name(flo_first, trait_name_type)
fru_first <- .pick_name(fru_first, trait_name_type)

# ---- Specimen ID -----------------------------------------------------------

# Create IDs before filtering so they remain stable across functions
if (!"id" %in% names(raw_data)) {
  raw_data$id <- sprintf("ID%04d", seq_len(nrow(raw_data)))
}

# ---- Keep the selected species ---------------------------------------------

if (!is.null(species_selected)) {
  raw_data <- raw_data[raw_data$taxon %in% species_selected, , drop = FALSE]
}

raw_data <- raw_data[!is.na(raw_data$taxon), , drop = FALSE]
raw_data <- raw_data[rowSums(!is.na(raw_data)) > 0, , drop = FALSE]

# ---- Morphometric matrix ---------------------------------------------------

start_idx <- match(veg_first, names(raw_data))

if (is.na(start_idx)) {
  stop("veg_first column not found: ", veg_first)
}

# Keep specimen ID and taxon with morphological variables
trait_range <- setdiff(names(raw_data)[start_idx:ncol(raw_data)], "id")

analysis_data <- raw_data[
  ,
  c("id", "taxon", trait_range),
  drop = FALSE
]

# ---- Numeric conversion ----------------------------------------------------

trait_cols <- names(analysis_data)[
  !(names(analysis_data) %in% c("id", "taxon"))
]

analysis_data[trait_cols] <- lapply(
  analysis_data[trait_cols],
  function(x) {
    x <- trimws(as.character(x))
    x <- gsub(",", ".", x, fixed = TRUE)
    x[x == ""] <- NA
    x <- suppressWarnings(as.numeric(x))
    x[is.infinite(x)] <- NA
    x
  }
)

# ---- Remove empty rows and columns -----------------------------------------

# Remove specimens with no morphological observations
analysis_data <- analysis_data[
  rowSums(!is.na(analysis_data[trait_cols])) > 0,
  ,
  drop = FALSE
]

# Remove empty morphological variables only
trait_cols <- names(analysis_data)[
  !(names(analysis_data) %in% c("id", "taxon"))
]

analysis_data <- analysis_data[
  ,
  c(
    "id",
    "taxon",
    trait_cols[
      colSums(!is.na(analysis_data[trait_cols])) > 0
    ]
  ),
  drop = FALSE
]

# Specimen IDs as row names
rownames(analysis_data) <- analysis_data$id
analysis_data$id <- NULL

# ---- Block column lists -----------------------------------------------------

block_starts <- c(
  veg = veg_first,
  flo = flo_first,
  fru = fru_first
)

col_names <- names(analysis_data)[
  !(names(analysis_data) %in% c("id", "taxon"))
]

idx <- match(block_starts, col_names)

if (anyNA(idx)) {
  missing_blocks <- names(block_starts)[is.na(idx)]

  warning(
    "The following morphological blocks could not be identified: ",
    paste(missing_blocks, collapse = ", ")
  )

  block_starts <- block_starts[!is.na(idx)]
  idx <- idx[!is.na(idx)]
}

ord <- order(idx)

ends <- c(idx[ord][-1] - 1, length(col_names))

base_cols <- stats::setNames(
  lapply(seq_along(ord), function(i) {
    col_names[idx[ord][i]:ends[i]]
  }),
  names(block_starts)[ord]
)

attr(analysis_data, "base_cols") <- base_cols

# ---- Taxon colours -----------------------------------------------------------
taxa_present <- sort(unique(analysis_data$taxon))

if (is.null(colors)) {
  taxon_colors <- stats::setNames(
    viridis::viridis(length(taxa_present)),
    taxa_present
  )
} else {
  if (is.null(names(colors)) || any(names(colors) == "")) {
    stop("`colors` must be a named character vector, e.g. ...")
  }
  missing_taxa <- setdiff(taxa_present, names(colors))
  if (length(missing_taxa) > 0) {
    stop("`colors` is missing an entry for: ", paste(missing_taxa, collapse = ", "))
  }
  taxon_colors <- colors[taxa_present]
}

attr(analysis_data, "taxon_colors") <- taxon_colors

return(analysis_data)
}

#' Select a morphological trait naming format
#'
#' @description
#' Internal helper used to extract either the descriptive name or abbreviated
#' code from \pkg{moRphoTaxa} trait names formatted as
#' `"full_name/TRAITcode"`.
#'
#' @param x Character vector of trait names.
#' @param type Character. Either `"full"`, `"code"`, or `NULL`.
#'
#' @return A character vector with trait names converted according to `type`.
#'
#' @keywords internal
#'
#' @noRd
.pick_name <- function(x, type) {
  if (is.null(type)) return(x)  # no change requested — keep "full/CODE" as-is

  parts <- strsplit(x, "/", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) == 2) {
      if (type == "full") p[1] else p[2]
    } else {
      p[1]  # no "/" in this column (e.g. "taxon", "id") -> leave untouched
    }
  }, character(1))
}
