# Tests for environmental data set-up and exploration:
# geo_matrix_setting(), geo_bioclim_exploratory(), geo_bioclim_maps()
# WorldClim/SoilGrids/Natural Earth downloads are mocked.

# Mocks for geodata downloads, returning synthetic rasters
mock_geodata_downloads <- function(env = parent.frame(), soil_fail = character(0)) {
  local_mocked_bindings(
    worldclim_global = function(var, res, path, ...) {
      if (var == "bio") {
        make_clim_stack(names = paste0("wc2.1_", res, "m_bio_", 1:19))
      } else {
        make_clim_stack(names = paste0("wc2.1_", res, "m_elev"), seed = 5)
      }
    },
    soil_world = function(var, depth, stat, path, ...) {
      if (var %in% soil_fail) stop("server unavailable")
      make_clim_stack(names = var, seed = nchar(var), res = 1)
    },
    .package = "geodata",
    .env = env
  )
}

#_______________________________________________________________________________
# geo_matrix_setting() ####

test_that("geo_matrix_setting extracts climate and soil values for each record", {
  local_test_dir()
  write_raw_xlsx("raw.xlsx")
  mock_geodata_downloads()

  expect_message(
    geo <- geo_matrix_setting("raw.xlsx", taxon_col = "species",
                              edaphic_vars = c("clay", "sand")),
    "geodata` built with"
  )

  # Records without taxon or coordinates are dropped
  expect_equal(nrow(geo), 28)
  expect_true(all(grepl("^ID[0-9]{4}$", rownames(geo))))
  expect_false("id" %in% names(geo))
  expect_type(geo$decimalLatitude, "double")
  expect_true(all(c(paste0("bio_", 1:19), "elevation", "clay", "sand") %in% names(geo)))

  env <- attr(geo, "env_cols")
  expect_equal(env$clim, c(paste0("bio_", 1:19), "elevation"))
  expect_equal(env$edaphic, c("clay", "sand"))
  expect_s4_class(attr(geo, "clim_stack"), "SpatRaster")
  expect_s4_class(attr(geo, "soil_stack"), "SpatRaster")
  expect_true(terra::compareGeom(attr(geo, "soil_stack"), attr(geo, "clim_stack")))
  expect_s4_class(attr(geo, "pts"), "SpatVector")
})

test_that("geo_matrix_setting skips edaphic data on request or failure", {
  local_test_dir()
  write_raw_xlsx("raw.xlsx", with_id = TRUE)
  mock_geodata_downloads(soil_fail = c("clay", "silt"))

  geo <- geo_matrix_setting("raw.xlsx", taxon_col = "species", edaphic = FALSE,
                            species_selected = "alpha", verbose = FALSE)
  expect_null(attr(geo, "env_cols")$edaphic)
  expect_null(attr(geo, "soil_stack"))
  expect_true(all(geo$taxon == "alpha"))
  expect_true(all(grepl("^SPEC", rownames(geo))))

  expect_warning(
    geo2 <- geo_matrix_setting("raw.xlsx", taxon_col = "species",
                               edaphic_vars = c("clay", "sand"), verbose = FALSE),
    "Could not fetch 'clay': server unavailable"
  )
  expect_equal(attr(geo2, "env_cols")$edaphic, "sand")

  w <- testthat::capture_warnings(
    geo3 <- geo_matrix_setting("raw.xlsx", taxon_col = "species",
                               edaphic_vars = c("clay", "silt"), verbose = FALSE)
  )
  expect_true(any(grepl("No edaphic variables were successfully extracted", w)))
  expect_null(attr(geo3, "env_cols")$edaphic)
})

test_that("geo_matrix_setting validates its inputs", {
  local_test_dir()
  write_raw_xlsx("raw.xlsx")
  mock_geodata_downloads()

  expect_error(geo_matrix_setting(1), "xlsx_path")
  expect_error(geo_matrix_setting("raw.xlsx", wc_res = 3), "wc_res")
  expect_error(geo_matrix_setting("raw.xlsx", edaphic = "yes"), "edaphic")
  expect_error(geo_matrix_setting("raw.xlsx", edaphic_vars = character(0)),
               "edaphic_vars")
  expect_error(geo_matrix_setting("raw.xlsx", edaphic_depth = 7), "edaphic_depth")
  expect_error(geo_matrix_setting("raw.xlsx", data_path = 1), "data_path")
  expect_error(geo_matrix_setting("raw.xlsx", verbose = "yes"), "verbose")
  expect_error(geo_matrix_setting("raw.xlsx"), "No `taxon` column found")
  expect_error(geo_matrix_setting("raw.xlsx", taxon_col = "species",
                                  species_selected = "nothing", verbose = FALSE),
               "No records remain")

  raw <- openxlsx::read.xlsx("raw.xlsx")
  raw$decimalLatitude <- NULL
  openxlsx::write.xlsx(raw, "no_coords.xlsx")
  expect_error(geo_matrix_setting("no_coords.xlsx", taxon_col = "species"),
               "decimalLatitude")
})

#_______________________________________________________________________________
# geo_bioclim_exploratory() ####

test_that("geo_bioclim_exploratory ranks bioclimatic variables", {
  local_test_dir()
  geo <- make_geodata()
  expect_message(
    out <- geo_bioclim_exploratory(geo, bio_x = 1, bio_y = 12,
                                   file_formats = c("pdf", "jpg")),
    "Top 5 most informative variables"
  )
  expect_null(out)

  od <- out_dir("Figs.bioclimPCA")
  files <- list.files(od)
  expect_true(all(c("BIOCLIM_variance_correlation.pdf",
                    "BIOCLIM_PCA.jpg",
                    "BIOCLIM_scatter_bio1bio12.pdf",
                    "BIOCLIM_variable_ranking.xlsx") %in% files))
  ranking <- openxlsx::read.xlsx(file.path(od, "BIOCLIM_variable_ranking.xlsx"),
                                 sheet = "ranked_variables")
  expect_equal(nrow(ranking), 19)
  expect_equal(ranking$rank, 1:19)
  expect_equal(openxlsx::getSheetNames(file.path(od, "BIOCLIM_variable_ranking.xlsx")),
               c("ranked_variables", "redundant_pairs"))
})

test_that("geo_bioclim_exploratory uses stored taxon colours and reports no redundancy", {
  local_test_dir()
  geo <- make_geodata()
  attr(geo, "taxon_colors") <- c(alpha = "red", beta = "blue", gamma = "green")
  geo_bioclim_exploratory(geo, cor_threshold = 1, file_formats = "pdf",
                          verbose = FALSE)
  pairs <- openxlsx::read.xlsx(
    file.path(out_dir("Figs.bioclimPCA"), "BIOCLIM_variable_ranking.xlsx"),
    sheet = "redundant_pairs"
  )
  expect_equal(pairs$highly_correlated_pairs, "None")
})

test_that("geo_bioclim_exploratory validates its inputs", {
  local_test_dir()
  geo <- make_geodata()
  expect_error(geo_bioclim_exploratory("x"), "must be a data frame")
  no_attr <- geo
  attr(no_attr, "env_cols") <- NULL
  expect_error(geo_bioclim_exploratory(no_attr), "env_cols")
  no_taxon <- geo
  no_taxon$taxon <- NULL
  expect_error(geo_bioclim_exploratory(no_taxon), "taxon")
  expect_error(geo_bioclim_exploratory(geo, bio_x = 20), "bio_x")
  expect_error(geo_bioclim_exploratory(geo, bio_y = 0), "bio_y")
  expect_error(geo_bioclim_exploratory(geo, cor_threshold = 2), "cor_threshold")
  expect_error(geo_bioclim_exploratory(geo, file_formats = character(0)), "file_formats")
  expect_error(geo_bioclim_exploratory(geo, file_formats = "png"))
  expect_error(geo_bioclim_exploratory(geo, verbose = "a"), "verbose")

  only_elev <- geo
  attr(only_elev, "env_cols") <- list(clim = "elevation")
  expect_error(geo_bioclim_exploratory(only_elev), "No bioclimatic")

  no_bio1 <- geo
  no_bio1$bio_1 <- NULL
  expect_error(geo_bioclim_exploratory(no_bio1, bio_x = 1, bio_y = 12),
               "bio_x` = 1 refers to column")
  expect_error(geo_bioclim_exploratory(no_bio1, bio_x = 12, bio_y = 1),
               "bio_y` = 1 refers to column")

  few <- geo
  few[3:nrow(few), "bio_5"] <- NA
  expect_error(geo_bioclim_exploratory(few, verbose = FALSE),
               "Fewer than 3 complete records")
})

#_______________________________________________________________________________
# geo_bioclim_maps() ####

mock_natural_earth <- function(env = parent.frame()) {
  local_mocked_bindings(
    ne_countries = fake_world_sf,
    ne_states = fake_world_sf,
    ne_download = fake_rivers_sf,
    .package = "rnaturalearth",
    .env = env
  )
}

test_that("geo_bioclim_maps draws the map, scatter and combined figure", {
  local_test_dir()
  mock_natural_earth()
  geo <- make_geodata()

  expect_message(
    res <- geo_bioclim_maps(geo, bio_map_layer = 12, bio_x = 1, bio_y = "elevation",
                            show_country_states = "Brazil", show_rivers = TRUE,
                            file_formats = c("pdf", "jpg")),
    "Outputs saved in"
  )
  expect_named(res, c("map", "scatter", "combined", "inset", "output_dir"))
  expect_s3_class(res$scatter, "ggplot")
  expect_s3_class(res$inset, "ggplot")
  expect_true(all(c("MAP_RASTER.pdf", "PLOT_RASTER.jpg", "MAP_RASTER_PLOT.pdf") %in%
                    list.files(res$output_dir)))
})

test_that("geo_bioclim_maps works without inset, saving or custom shapes", {
  local_test_dir()
  mock_natural_earth()
  geo <- make_geodata()

  expect_message(
    res <- geo_bioclim_maps(geo, bio_map_layer = "elevation", inset = FALSE,
                            pt_shapes = c(alpha = 21, beta = 22, gamma = 24),
                            save = FALSE),
    "not saved"
  )
  expect_null(res$inset)
  expect_null(res$output_dir)
  expect_s3_class(res$map, "ggplot")
  expect_false(dir.exists("Figs.bioclimaps"))
})

test_that("geo_bioclim_maps validates its inputs", {
  local_test_dir()
  mock_natural_earth()
  geo <- make_geodata()

  expect_error(geo_bioclim_maps("x"), "must be the data frame")
  expect_error(geo_bioclim_maps(geo[, c("taxon", "bio_1")]), "decimalLatitude")
  no_stack <- geo
  attr(no_stack, "clim_stack") <- NULL
  expect_error(geo_bioclim_maps(no_stack), "clim_stack")
  expect_error(geo_bioclim_maps(geo, show_rivers = "yes"), "show_rivers")
  expect_error(geo_bioclim_maps(geo, show_rivers = TRUE, river_scale = 20), "river_scale")
  expect_error(geo_bioclim_maps(geo, show_country_states = c("A", "B")),
               "show_country_states")
  expect_error(geo_bioclim_maps(geo, inset = 1), "inset")
  expect_error(geo_bioclim_maps(geo, save = "no"), "save")
  expect_error(geo_bioclim_maps(geo, file_formats = character(0)), "file_formats")
  expect_error(geo_bioclim_maps(geo, verbose = "a"), "verbose")
  expect_error(geo_bioclim_maps(geo, pt_shapes = c(21, 22, 24)), "named vector")
  expect_error(geo_bioclim_maps(geo, pt_shapes = c(alpha = 21)),
               "missing an entry for: beta, gamma")
  expect_error(geo_bioclim_maps(geo, bio_map_layer = 0), "Elevation is not a bioclimatic")

  many <- geo
  many$taxon <- paste0("t", seq_len(nrow(many)))
  expect_error(geo_bioclim_maps(many), "automatic point shapes")
})

test_that("bioclim map helpers resolve variables, labels and extents", {
  avail <- c(paste0("bio_", 1:19), "elevation")
  expect_equal(.resolve_clim_var(12, avail, "x"), "bio_12")
  expect_equal(.resolve_clim_var(3, "bio3", "x"), "bio3")
  expect_equal(.resolve_clim_var("elevation", avail, "x"), "elevation")
  expect_error(.resolve_clim_var(c(1, 2), avail, "x"), "single value")
  expect_error(.resolve_clim_var(NA, avail, "x"), "single value")
  expect_error(.resolve_clim_var(25, avail, "x"), "from 1 to 19")
  expect_error(.resolve_clim_var("soil", avail, "x"), "was not found")

  expect_equal(.clim_labels("elevation")$short, "Elevation (m)")
  expect_equal(.clim_labels("bio_12")$short, "Annual Precipitation (mm)")
  expect_match(.clim_labels("bio_1")$axis, "^BIO1 ")
  expect_equal(.clim_labels("clay"), list(short = "clay", axis = "clay"))

  ext <- .auto_extent(c(-50, -40), c(-20, -10))
  expect_named(ext, c("map", "crop", "inset"))
  expect_equal(unname(ext$map["xmin"]), -51.5)
  expect_true(ext$crop[["xmin"]] < ext$map[["xmin"]])
  world <- .auto_extent(c(-179, 179), c(-89, 89))
  expect_equal(unname(world$inset), c(-180, 180, -90, 90))
})
