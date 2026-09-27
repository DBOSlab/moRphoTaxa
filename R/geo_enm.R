#' Ensemble species distribution models (ENM) for each taxon
#'
#' @description
#' Fits ensemble species distribution models (ecological niche models, ENM)
#' for every taxon in the geographical/environmental matrix produced by
#' [geo_matrix_setting()], using the \pkg{flexsdm} framework. For each taxon,
#' the function delimits a calibration area, samples pseudo-absences, fits up
#' to five algorithms (GLM, GAM, GBM, neural network and random forest),
#' combines them into a weighted ensemble, predicts habitat suitability across
#' the study region, estimates predictor importance, and produces a
#' suitability map with a locator inset.
#'
#' Predictors are taken automatically from the raster stacks stored as
#' attributes of `geodata`: the climate stack (`"clim_stack"`) and, when
#' available, the edaphic stack (`"soil_stack"`). The study region and the
#' locator inset are derived automatically from the geographical distribution
#' of the records.
#'
#' @details
#' \strong{Input.} `geodata` must be the data frame returned by
#' [geo_matrix_setting()], carrying the `"clim_stack"` attribute and,
#' optionally, the `"soil_stack"` attribute. It must contain the columns
#' `taxon`, `decimalLatitude` and `decimalLongitude`. All taxa present in
#' `geodata` are modelled; to restrict the analysis to some taxa, use
#' `species_selected` in [geo_matrix_setting()].
#'
#' \strong{Predictors.} If `"soil_stack"` is present, climate and edaphic
#' layers are combined (`env_source = "clim_edaphic"`); otherwise only climate
#' layers are used (`env_source = "clim"`). The soil stack is resampled
#' (bilinear) to the climate raster grid only if their geometries differ.
#' Highly collinear predictors are then removed: pairwise Pearson correlations
#' are computed on a random sample of 5000 raster cells, and predictors are
#' dropped with [caret::findCorrelation()] using a cutoff of |r| > 0.8.
#' Layers with no variation in the sample are also dropped. The predictors
#' retained are used for all taxa.
#'
#' \strong{Study region.} The predictor stack is cropped to the bounding box of
#' all occurrence records, expanded by `buffer_km` (converted to degrees at
#' about 111 km per degree) plus a 2-degree safety margin. This avoids
#' predicting over the whole globe and reduces memory use. Memory fraction of
#' \pkg{terra} is lowered to 0.2 during the run and restored on exit.
#'
#' \strong{Workflow per taxon.}
#' \enumerate{
#'   \item Occurrence records are selected. Taxa with fewer than `min_occ`
#'   records are skipped with a warning.
#'   \item A calibration area is built with [flexsdm::calib_area()]
#'   (`"bmcp"` method: buffered minimum convex polygon of width `buffer_km`),
#'   and predictors are cropped and masked to it.
#'   \item Pseudo-absences are sampled with [flexsdm::sample_pseudoabs()]
#'   (`"env_const"` method: environmentally constrained), at `pa_ratio`
#'   pseudo-absences per presence.
#'   \item Predictor values are extracted at presences and pseudo-absences,
#'   incomplete rows are removed, and data are split into `folds`
#'   cross-validation partitions ([flexsdm::part_random()], `"kfold"`).
#'   \item The five algorithms are fitted, each tuned over a small
#'   hyperparameter grid and thresholded at the maximum sensitivity + specificity
#'   (`"max_sens_spec"`). An algorithm that fails to converge is dropped; if
#'   fewer than two algorithms converge, the taxon is skipped.
#'   \item Models are combined with [flexsdm::fit_ensemble()], weighting by the
#'   Sorensen index, and habitat suitability (0-1) is predicted across the
#'   study region with [flexsdm::sdm_predict()].
#'   \item Predictor importance is estimated by permutation with
#'   [flexsdm::sdm_varimp()], repeated `importance_perm` times, and averaged
#'   over the Sorensen metric.
#'   \item A suitability map with the occurrence points, optional state
#'   borders and a locator inset is built.
#' }
#'
#' A failure at any step for one taxon produces a warning and moves on to the
#' next taxon, so a single problematic taxon does not stop the whole run. The
#' number of records must be at least `folds` for cross-validation partitions
#' to be created.
#'
#' \strong{Output.} When `save = TRUE`, files are written to a dated sub-folder
#' of `"Figs_ENM"` in the working directory (for example
#' `"Figs_ENM/18Sep2026"`), created if needed. For each successfully modelled
#' taxon:
#'
#' \itemize{
#'   \item `ENM_map_<taxon>.pdf` and/or `.jpg`: suitability map;
#'   \item `ENM_suitability_<taxon>.tif`: suitability raster (GeoTIFF), for
#'   use in downstream niche-overlap analyses;
#'   \item `ENM_importance_<taxon>.pdf` and/or `.jpg`: predictor importance
#'   plot;
#'   \item `ENM_<taxon>.xlsx`: ensemble performance and predictor importance.
#' }
#'
#' Plus `ENM_summary.xlsx`, one row per modelled taxon. Characters other than
#' letters, digits, `.`, `_` and `-` in taxon names are replaced by `_` in file
#' names.
#'
#' @param geodata Data frame returned by [geo_matrix_setting()], carrying the
#' `"clim_stack"` attribute and, optionally, the `"soil_stack"` attribute.
#'
#' @param buffer_km Numeric. Width, in kilometres, of the buffer used to
#' build each taxon's calibration area and to expand the study region.
#' Default is `100`.
#'
#' @param pa_ratio Numeric. Number of pseudo-absences sampled per presence.
#' Default is `2`.
#'
#' @param folds Integer. Number of cross-validation folds. Default is `5`.
#'
#' @param min_occ Integer. Minimum number of occurrence records required to
#' model a taxon. Default is `3`. Note that `folds` also needs to be
#' satisfiable by the number of records.
#'
#' @param importance_perm Integer. Number of permutations used to estimate
#' predictor importance. Default is `10`.
#'
#' @param seed Integer or `NULL`. Random seed set at the start of each
#' taxon's model (pseudo-absence sampling, cross-validation partitions,
#' algorithm fitting and permutation importance), so results for a taxon do
#' not depend on the other taxa in `geodata`. Default is `123`. If `NULL`,
#' the random state is left untouched.
#'
#' @param show_country_states Character string with a country name (e.g.
#' `"brazil"`) whose sub-national borders are drawn on the maps, passed to
#' `rnaturalearth::ne_states()` (requires the \pkg{rnaturalearthhires}
#' package). If `NULL` (default), no sub-national borders are drawn.
#'
#' @param inset Logical. If `TRUE` (default), a small locator map is added to
#' each suitability map, with a red rectangle marking the study region.
#'
#' @param save Logical. If `TRUE` (default), maps, rasters, plots and tables
#' are written to disk. If `FALSE`, nothing is saved and results are only
#' returned.
#'
#' @param file_formats Character vector specifying which image formats to
#' save. One or more of `"pdf"` and `"jpg"`. Default is
#' `c("pdf", "jpg")`. Ignored when `save = FALSE`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#' including retained predictors, taxa modelled, algorithms fitted and
#' ensemble performance.
#'
#' @return
#' Invisibly, a named list with the elements:
#' \describe{
#'   \item{`summary`}{A data frame with one row per modelled taxon:
#'   `taxon`, `n_occ` (records used), `n_algos` (algorithms in the ensemble),
#'   `algos` (algorithm names separated by `"|"`), `env_source`
#'   (`"clim"` or `"clim_edaphic"`) and `sorensen` (ensemble Sorensen index).
#'   `NULL` if no taxon could be modelled.}
#'   \item{`suitability`}{A named list of [terra::SpatRaster] objects with
#'   the predicted habitat suitability (0-1) of each modelled taxon.}
#'   \item{`maps`}{A named list of suitability maps ([ggplot2::ggplot()] or
#'   \pkg{cowplot} objects), one per modelled taxon.}
#'   \item{`predictors`}{Character vector with the names of the predictors
#'   retained after the collinearity filter.}
#'   \item{`output_dir`}{Path of the folder where files were saved, or `NULL`
#'   if `save = FALSE`.}
#' }
#'
#' @seealso
#' [geo_matrix_setting()], [geo_bioclim_maps()]
#'
#' @examples
#' \dontrun{
#' geo <- geo_matrix_setting(
#'   xlsx_path = "occurrence_data.xlsx",
#'   taxon_col = "species",
#'   edaphic = TRUE
#' )
#'
#' # Default run: climate + edaphic predictors
#' enm <- geo_enm(geo)
#'
#' # Custom settings, with state borders
#' enm <- geo_enm(
#'   geo,
#'   buffer_km = 100,
#'   pa_ratio = 2,
#'   folds = 5,
#'   min_occ = 3,
#'   importance_perm = 10,
#'   show_country_states = "brazil"
#' )
#'
#' # Summary table and a suitability raster
#' enm$summary
#' terra::plot(enm$suitability[["Species A"]])
#'
#' # Climate-only models, without saving files
#' geo_clim <- geo_matrix_setting("occurrence_data.xlsx",
#'                                taxon_col = "species", edaphic = FALSE)
#' enm <- geo_enm(geo_clim, save = FALSE)
#'
#' # Save PDF maps only
#' enm <- geo_enm_ensemble(geo, file_formats = "pdf")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_raster geom_sf geom_point geom_col
#' @importFrom ggplot2 coord_sf coord_flip labs theme theme_bw theme_void
#' @importFrom ggplot2 element_blank element_rect scale_fill_viridis_c
#' @importFrom cowplot ggdraw draw_plot save_plot
#' @importFrom terra terraOptions spatSample nlyr ext crop mask extract
#' @importFrom terra resample compareGeom writeRaster
#' @importFrom caret findCorrelation
#' @importFrom flexsdm calib_area sample_pseudoabs part_random fit_glm fit_gam
#' @importFrom flexsdm fit_gbm fit_net fit_raf fit_ensemble sdm_predict sdm_varimp
#' @importFrom rnaturalearth ne_countries ne_states
#' @importFrom sf st_sfc st_polygon
#' @importFrom dplyr bind_rows
#' @importFrom openxlsx write.xlsx
#' @importFrom rlang .data
#' @importFrom stats complete.cases cor na.omit aggregate as.formula sd
#'
#' @export

geo_enm <- function(geodata,
                     buffer_km = 100,
                     pa_ratio = 2,
                     folds = 5,
                     min_occ = 3,
                     importance_perm = 10,
                     seed = 123,
                     show_country_states = NULL,
                     inset = TRUE,
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
soil_stack <- attr(geodata, "soil_stack")

if (is.null(clim_stack) || !inherits(clim_stack, "SpatRaster")) {
  stop(
    "`geodata` has no valid \"clim_stack\" attribute. ",
    "Build it with `geo_matrix_setting()`."
  )
}

if (!is.numeric(buffer_km) || length(buffer_km) != 1L || buffer_km <= 0) {
  stop("`buffer_km` must be a single positive number.")
}

if (!is.numeric(pa_ratio) || length(pa_ratio) != 1L || pa_ratio <= 0) {
  stop("`pa_ratio` must be a single positive number.")
}

if (!is.numeric(folds) || length(folds) != 1L || folds < 2) {
  stop("`folds` must be a single integer of at least 2.")
}

if (!is.numeric(min_occ) || length(min_occ) != 1L || min_occ < 2) {
  stop("`min_occ` must be a single integer of at least 2.")
}

if (!is.numeric(importance_perm) || length(importance_perm) != 1L || importance_perm < 1) {
  stop("`importance_perm` must be a single positive integer.")
}

if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
  stop("`seed` must be NULL or a single number.")
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

base_dir <- "Figs_ENM"   # parent folder for saved outputs
cor_cutoff <- 0.8        # |r| above which predictors are dropped
n_sample <- 5000       # cells sampled for the correlation matrix

if (verbose) message("Running Ensemble Species Distribution Models...")

# Lower memory usage to avoid allocation failures during prediction;
# the previous setting is restored on exit
old_opts <- tryCatch(terra::terraOptions(print = FALSE), error = function(e) NULL)
terra::terraOptions(memfrac = 0.2)
on.exit(
  if (is.list(old_opts) && !is.null(old_opts$memfrac)) {
    terra::terraOptions(memfrac = old_opts$memfrac)
  },
  add = TRUE
)

# ---- Output folder -------------------------------------------------------------

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

# ---- Predictor raster stack (climate + edaphic if available) --------------------

use_soil <- !is.null(soil_stack) && inherits(soil_stack, "SpatRaster")
env_source <- if (use_soil) "clim_edaphic" else "clim"

if (use_soil) {
  if (!terra::compareGeom(soil_stack, clim_stack, stopOnError = FALSE)) {
    soil_stack <- terra::resample(soil_stack, clim_stack, method = "bilinear")
  }
  pred_stack_full <- c(clim_stack, soil_stack)
} else {
  pred_stack_full <- clim_stack
}

# Drop predictors with no variation and highly correlated predictors
# (|r| > cutoff) to avoid overfitting/singular models
sample_vals <- terra::spatSample(pred_stack_full, n_sample, na.rm = TRUE)

const_vars <- names(sample_vals)[vapply(sample_vals, function(x) stats::sd(x) == 0, logical(1))]
const_vars <- const_vars[!is.na(const_vars)]

if (length(const_vars) > 0) {
  sample_vals <- sample_vals[, setdiff(names(sample_vals), const_vars), drop = FALSE]
}

cor_mat <- stats::cor(sample_vals, use = "pairwise.complete.obs")
high_cor <- caret::findCorrelation(cor_mat, cutoff = cor_cutoff, names = TRUE)
drop_vars <- unique(c(const_vars, high_cor))

if (length(drop_vars) > 0) {
  if (verbose) {
    message("  Dropping ", length(drop_vars), " constant/collinear predictors: ",
            paste(drop_vars, collapse = ", "))
  }
  pred_stack <- pred_stack_full[[setdiff(names(pred_stack_full), drop_vars)]]
} else {
  pred_stack <- pred_stack_full
}

if (terra::nlyr(pred_stack) < 2) {
  stop("Fewer than 2 predictors remain after removing collinear variables.")
}

if (verbose) {
  message("  Environmental source: ", env_source)
  message("  Predictors (", terra::nlyr(pred_stack), "): ",
          paste(names(pred_stack), collapse = ", "))
}

# ---- Crop predictors to the study region (automatic) ----------------------------

margin_deg <- (buffer_km / 111) + 2   # ~111 km per degree + 2 degrees safety margin

lon_rng <- range(geodata$decimalLongitude, na.rm = TRUE)
lat_rng <- range(geodata$decimalLatitude,  na.rm = TRUE)

study_bbox <- c(
  xmin = max(lon_rng[1] - margin_deg, -180),
  xmax = min(lon_rng[2] + margin_deg,  180),
  ymin = max(lat_rng[1] - margin_deg,  -90),
  ymax = min(lat_rng[2] + margin_deg,   90)
)

study_extent <- terra::ext(
  study_bbox[["xmin"]], study_bbox[["xmax"]],
  study_bbox[["ymin"]], study_bbox[["ymax"]]
)

pred_stack <- terra::crop(pred_stack, study_extent)

if (verbose) {
  message("  Study extent (margin = ", round(margin_deg, 1), " deg): ",
          paste(round(study_bbox, 2), collapse = ", "))
}

# ---- Reusable locator inset (same for every taxon) --------------------------------

p_inset <- NULL

if (inset) {

inset_lim <- .enm_inset_limits(study_bbox)

world_sf <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")

bbox_ring <- matrix(
  c(study_bbox[["xmin"]], study_bbox[["ymin"]],
    study_bbox[["xmax"]], study_bbox[["ymin"]],
    study_bbox[["xmax"]], study_bbox[["ymax"]],
    study_bbox[["xmin"]], study_bbox[["ymax"]],
    study_bbox[["xmin"]], study_bbox[["ymin"]]),   # close ring
  ncol = 2, byrow = TRUE
)

bbox_poly <- sf::st_sfc(sf::st_polygon(list(bbox_ring)), crs = 4326)

p_inset <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = world_sf, fill = "grey85", colour = "grey60",
                   linewidth = 0.15) +
  ggplot2::geom_sf(data = bbox_poly, fill = NA, colour = "red",
                   linewidth = 0.9) +
  ggplot2::coord_sf(
    xlim = c(inset_lim[["xmin"]], inset_lim[["xmax"]]),
    ylim = c(inset_lim[["ymin"]], inset_lim[["ymax"]]),
    expand = FALSE
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white", colour = "grey40",
                                             linewidth = 0.4)
  )
}

# ---- Sub-national borders (downloaded once) ------------------------------------------

borders <- if (!is.null(show_country_states)) {
  rnaturalearth::ne_states(country = show_country_states, returnclass = "sf")
} else {
  NULL
}

# ---- Taxa to model (all taxa present in geodata) -------------------------------------

run_taxa <- sort(unique(stats::na.omit(geodata$taxon)))
run_taxa <- run_taxa[nzchar(trimws(run_taxa))]

if (length(run_taxa) == 0) {
  stop("No valid taxa found in geodata$taxon.")
}

if (verbose) {
  message("  Taxa to model: ", paste(run_taxa, collapse = ", "))
}

# ---- Helper: run ENM for one taxon ------------------------------------------------------

run_enm_taxon <- function(taxon_name) {

if (verbose) message("\n  --- Taxon: ", taxon_name, " ---")

if (!is.null(seed)) set.seed(seed)

file_stub <- gsub("[^[:alnum:]._-]+", "_", taxon_name)

occ <- geodata[geodata$taxon == taxon_name, ]
occ <- occ[stats::complete.cases(occ[, c("decimalLongitude", "decimalLatitude")]), ]

if (nrow(occ) < min_occ) {
  warning("Skipping '", taxon_name, "': fewer than ", min_occ,
          " occurrence records (n = ", nrow(occ), ").")
  return(NULL)
}

occ_pts <- data.frame(
  id = seq_len(nrow(occ)),
  x = occ$decimalLongitude,
  y = occ$decimalLatitude,
  pr_ab = 1
)

# -- Calibration area
ca <- .safe_run(
  flexsdm::calib_area(
    data = occ_pts,
    x = "x", y = "y",
    method = c("bmcp", width = buffer_km * 1000),
    crs = "EPSG:4326"
  ),
  paste0("Calibration area failed for '", taxon_name, "'")
)

if (is.null(ca)) return(NULL)

pred_stack_ca <- terra::crop(pred_stack, ca)
pred_stack_ca <- terra::mask(pred_stack_ca, ca)

# -- Pseudo-absences
pa <- .safe_run(
  flexsdm::sample_pseudoabs(
    data = occ_pts,
    x = "x", y = "y",
    n = nrow(occ_pts) * pa_ratio,
    method = c("env_const", env = pred_stack_ca),
    rlayer = pred_stack_ca[[1]],
    calibarea = ca
  ),
  paste0("Pseudo-absence sampling failed for '", taxon_name, "'")
)

if (is.null(pa)) return(NULL)
pa$pr_ab <- 0

# -- Modelling data
full_data <- dplyr::bind_rows(occ_pts, pa)

env_vals <- terra::extract(pred_stack, full_data[, c("x", "y")])[, -1, drop = FALSE]
full_data <- cbind(full_data, env_vals)
full_data <- full_data[stats::complete.cases(full_data[, names(env_vals)]), ]

if (sum(full_data$pr_ab == 1) < min_occ) {
  warning("Skipping '", taxon_name, "': fewer than ", min_occ,
          " presences remain after extracting predictor values.")
  return(NULL)
}

full_data <- .safe_run(
  flexsdm::part_random(
    data = full_data,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = folds)
  ),
  paste0("Data partitioning failed for '", taxon_name, "' (too few records for ",
         folds, " folds?)")
)

if (is.null(full_data)) return(NULL)

preds <- names(env_vals)

# -- Algorithms
fits <- list(
  glm = .try_fit(
    flexsdm::fit_glm(data = full_data, response = "pr_ab",
                     predictors = preds, partition = ".part",
                     thr = "max_sens_spec", poly = 1)
  ),
  gam = .try_fit(
    flexsdm::fit_gam(data = full_data, response = "pr_ab",
                     predictors = preds, partition = ".part",
                     thr = "max_sens_spec", k = 3)
  ),
  gbm = .try_fit(
    flexsdm::fit_gbm(data = full_data, response = "pr_ab",
                     predictors = preds, partition = ".part",
                     thr = "max_sens_spec",
                     grid = expand.grid(n.trees = c(50, 100, 250, 500),
                                        shrinkage = c(0.1, 0.5, 1),
                                        n.minobsinnode = c(1, 3, 5, 7, 9)))
  ),
  net = .try_fit(
    flexsdm::fit_net(data = full_data, response = "pr_ab",
                     predictors = preds, partition = ".part",
                     thr = "max_sens_spec",
                     grid = expand.grid(size = c(2, 4, 6, 8, 10),
                                        decay = c(0.001, 0.05, 0.1, 0.5, 1, 3)))
  ),
  raf = .try_fit(
    flexsdm::fit_raf(data = full_data, response = "pr_ab",
                     predictors = preds, partition = ".part",
                     thr = "max_sens_spec",
                     grid = expand.grid(mtry = c(1, 3, 5)))
  )
)

failed <- names(fits)[vapply(fits, is.null, logical(1))]
fits <- Filter(Negate(is.null), fits)

if (length(fits) < 2) {
  warning("Fewer than 2 algorithms converged for '", taxon_name, "'. Skipping ensemble.")
  return(NULL)
}

if (verbose) {
  message("  Algorithms fitted: ", paste(names(fits), collapse = ", "))
  if (length(failed) > 0) {
    message("  Algorithms that failed: ", paste(failed, collapse = ", "))
  }
}

# -- Ensemble
ens <- .safe_run(
  flexsdm::fit_ensemble(
    models = fits,
    thr = "max_sens_spec",
    thr_model = "max_sens_spec",
    metric = "SORENSEN"
  ),
  paste0("Ensemble fitting failed for '", taxon_name, "'")
)

if (is.null(ens)) return(NULL)

sorensen_col <- grep("sorensen", names(ens$performance),
                     ignore.case = TRUE, value = TRUE)[1]

if (is.na(sorensen_col)) {
  warning("No SORENSEN-like column found in ens$performance for '", taxon_name,
          "'. Columns available: ", paste(names(ens$performance), collapse = ", "))
  sorensen_val <- NA_real_
} else {
  sorensen_val <- ens$performance[[sorensen_col]][1]
}

if (verbose) message("  Ensemble Sorensen : ", round(sorensen_val, 3))

gc(verbose = FALSE)

# -- Prediction
suit_map <- .safe_run(
  flexsdm::sdm_predict(models = ens, pred = pred_stack, con_thr = FALSE)[[1]],
  paste0("Prediction failed for '", taxon_name, "'")
)

if (is.null(suit_map)) return(NULL)

suit_df <- as.data.frame(suit_map, xy = TRUE, na.rm = TRUE)
names(suit_df)[3] <- "suitability"

# -- Suitability map
p_map <- ggplot2::ggplot() +
  ggplot2::geom_raster(
    data = suit_df,
    ggplot2::aes(x = .data$x, y = .data$y, fill = .data$suitability)
  )

if (!is.null(borders)) {
  p_map <- p_map +
    ggplot2::geom_sf(data = borders, fill = NA, colour = "white",
                     linewidth = 0.3, inherit.aes = FALSE)
}

p_map <- p_map +
  ggplot2::geom_point(
    data = occ,
    ggplot2::aes(x = .data$decimalLongitude, y = .data$decimalLatitude),
    shape = 21, size = 1.5, fill = "white", color = "black", stroke = 0.3
  ) +
  ggplot2::scale_fill_viridis_c(option = "turbo", limits = c(0, 1),
                                name = "Habitat\nsuitability") +
  ggplot2::coord_sf(xlim = range(suit_df$x), ylim = range(suit_df$y)) +
  ggplot2::labs(
    title = paste0("ENM \u2014 ", taxon_name, " (", env_source, ")"),
    x = "Longitude", y = "Latitude"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

if (!is.null(p_inset)) {
  p_map <- cowplot::ggdraw(p_map) +
    cowplot::draw_plot(p_inset, x = 0.68, y = 0.02, width = 0.28, height = 0.28)
}

if (save) {
  save_outputs(
    paste0("ENM_map_", file_stub), p_map,
    base_height = 7, base_aspect_ratio = 1.1
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  .safe_run(
    terra::writeRaster(
      suit_map,
      file.path(output_dir, paste0("ENM_suitability_", file_stub, ".tif")),
      overwrite = TRUE
    ),
    paste0("Saving suitability raster failed for '", taxon_name, "'")
  )
}

# -- Predictor importance
importance_df <- .compute_variable_importance(fits, full_data, preds, importance_perm)

if (!is.null(importance_df) && nrow(importance_df) > 0) {

  # SORENSEN is the importance metric (matches the ensemble weighting metric)
  imp_metric <- "SORENSEN"

  if (!"predictors" %in% names(importance_df) || !imp_metric %in% names(importance_df)) {
    warning("Expected columns not found in sdm_varimp() output for '", taxon_name,
            "'. Columns found: ", paste(names(importance_df), collapse = ", "))
    importance_df <- NULL
  } else {
    importance_df <- stats::aggregate(
      stats::as.formula(paste0(imp_metric, " ~ predictors")),
      data = importance_df, FUN = mean, na.rm = TRUE
    )
    names(importance_df) <- c("predictor", "importance")
  }
}

if (!is.null(importance_df) && nrow(importance_df) > 0) {

  importance_df$taxon <- taxon_name
  importance_df <- importance_df[order(importance_df$importance, decreasing = TRUE), ]
  importance_df$predictor <- factor(importance_df$predictor,
                                    levels = rev(importance_df$predictor))

  p_imp <- ggplot2::ggplot(
    importance_df,
    ggplot2::aes(x = .data$predictor, y = .data$importance)
  ) +
    ggplot2::geom_col(fill = "#2166ac") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL, y = "Permutation importance",
      title = paste0("Variable importance \u2014 ", taxon_name)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  if (save) {
    save_outputs(
      paste0("ENM_importance_", file_stub), p_imp,
      base_height = max(4, nrow(importance_df) * 0.3), base_aspect_ratio = 1.5
    )
  }
} else {
  warning("Permutation importance could not be computed for '", taxon_name, "'.")
}

if (save) {
  openxlsx::write.xlsx(
    list(
      ensemble_performance = ens$performance,
      predictor_importance = if (!is.null(importance_df)) importance_df else data.frame()
    ),
    file.path(output_dir, paste0("ENM_", file_stub, ".xlsx"))
  )
}

result <- data.frame(
  taxon = taxon_name,
  n_occ = nrow(occ),
  n_algos = length(fits),
  algos = paste(names(fits), collapse = "|"),
  env_source = env_source,
  sorensen = round(sorensen_val, 4)
)

rm(fits, ens, full_data, env_vals, importance_df)
gc(verbose = FALSE)

list(row = result, suit_map = suit_map, map = p_map)
}

# ---- Run all taxa -------------------------------------------------------------------------

all_results <- Filter(Negate(is.null), lapply(run_taxa, run_enm_taxon))

if (length(all_results) == 0) {
  warning("No taxon could be modelled.")
  summary_df <- NULL
  suitability <- list()
  maps <- list()
} else {
  summary_df <- do.call(rbind, lapply(all_results, `[[`, "row"))
  suitability <- stats::setNames(lapply(all_results, `[[`, "suit_map"), summary_df$taxon)
  maps <- stats::setNames(lapply(all_results, `[[`, "map"),      summary_df$taxon)

  if (save) {
    openxlsx::write.xlsx(summary_df, file.path(output_dir, "ENM_summary.xlsx"))
  }
}

if (verbose) {
  if (save) {
    message("\nDone. Outputs saved in: ", output_dir)
    message("Suitability rasters (ENM_suitability_*.tif) can be reused in niche-overlap analyses.")
  } else {
    message("\nDone. Outputs were not saved (`save = FALSE`).")
  }
}

invisible(list(
  summary = summary_df,
  suitability = suitability,
  maps = maps,
  predictors = names(pred_stack),
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
.safe_run <- function(expr, fail_msg) {
  tryCatch(expr, error = function(e) {
    warning(fail_msg, ": ", conditionMessage(e), call. = FALSE)
    NULL
  })
}

#' Fit an algorithm silently
#'
#' @description
#' Internal helper. Evaluates a model-fitting expression and returns `NULL`
#' if it fails, so that a non-converging algorithm is simply left out of the
#' ensemble.
#'
#' @param expr Expression fitting a \pkg{flexsdm} model.
#'
#' @return The fitted model, or `NULL` on error.
#'
#' @keywords internal
#'
#' @noRd
.try_fit <- function(expr) {
  tryCatch(expr, error = function(e) NULL)
}

#' Variable importance via flexsdm::sdm_varimp()
#'
#' @description
#' Internal helper computing permutation-based variable importance for the
#' fitted models of one taxon. Returns `NULL` (with a warning) if the
#' computation fails.
#'
#' @param fits Named list of fitted \pkg{flexsdm} models.
#' @param full_data Data frame with the response, predictors and partitions.
#' @param preds Character vector of predictor names.
#' @param n_sim Integer. Number of permutations.
#'
#' @return A data frame with importance values, or `NULL`.
#'
#' @keywords internal
#'
#' @noRd
.compute_variable_importance <- function(fits, full_data, preds, n_sim) {
vi <- tryCatch(
  flexsdm::sdm_varimp(
    models = unname(fits),
    data = full_data,
    response = "pr_ab",
    predictors = preds,
    thr = "max_sens_spec",
    n_sim = n_sim,
    n_cores = 1
  ),
  error = function(e) {
    warning("sdm_varimp failed: ", conditionMessage(e), call. = FALSE)
    NULL
  }
)
if (is.null(vi) || NROW(vi) == 0) return(NULL)
as.data.frame(vi)
}

#' Compute the locator-inset extent from the study region
#'
#' @description
#' Internal helper. Extends the study-region bounding box by its own width and
#' height (at least 15 degrees) in each direction, clamped to valid
#' coordinates, so the inset shows the study region in its wider context.
#'
#' @param bbox Named numeric vector with `xmin`, `xmax`, `ymin` and `ymax`.
#'
#' @return A named numeric vector with `xmin`, `xmax`, `ymin` and `ymax`.
#'
#' @keywords internal
#'
#' @noRd
.enm_inset_limits <- function(bbox) {

w <- bbox[["xmax"]] - bbox[["xmin"]]
h <- bbox[["ymax"]] - bbox[["ymin"]]

ix <- max(w, 15)
iy <- max(h, 15)

c(
  xmin = max(bbox[["xmin"]] - ix, -180),
  xmax = min(bbox[["xmax"]] + ix,  180),
  ymin = max(bbox[["ymin"]] - iy,  -90),
  ymax = min(bbox[["ymax"]] + iy,   90)
)
}
