# Shared fixtures for the moRphoTaxa test suite
#
# Every analysis function writes figures/tables into the working directory,
# so tests run inside a temporary directory created with `local_test_dir()`.
# All data are synthetic and every download (WorldClim, SoilGrids, Natural
# Earth) is mocked, so the suite runs offline.

# Several plot titles contain an em dash, which the default Latin-1 encoding
# of pdf() cannot represent (one warning per glyph). WinAnsi (CP1252) has it.
grDevices::pdf.options(encoding = "WinAnsi.enc")

#_______________________________________________________________________________
# Temporary working directory ####

local_test_dir <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(pattern = "moRphoTaxa_test_", .local_envir = env)
  withr::local_dir(dir, .local_envir = env)
  invisible(dir)
}

# Run an expression silencing messages and printed output (warnings still
# propagate, so they can be tested with expect_warning()).
quiet <- function(expr) {
  utils::capture.output(res <- suppressMessages(expr))
  invisible(res)
}

# Today's dated output folder, as built by the analysis functions
out_dir <- function(base) {
  file.path(base, format(Sys.time(), "%d%b%Y"))
}

#_______________________________________________________________________________
# Morphometric matrix (as returned by morph_matrix_setting()) ####

make_morph_matrix <- function(n_per = 12,
                              taxa = c("alpha", "beta", "gamma"),
                              seed = 1,
                              with_na = TRUE,
                              colors = TRUE) {
  set.seed(seed)
  n <- n_per * length(taxa)
  taxon <- rep(taxa, each = n_per)
  shift <- as.numeric(factor(taxon, levels = taxa))

  veg <- c("petiole_length/PETIlng", "petiole_width/PETIwid",
           "leaf_length/LEAFlng", "leaf_width/LEAFwid")
  flo <- c("inflorescence_length/INFLlng", "pedicel_length/PDCLlng",
           "calyx_length/CLYXlng")
  fru <- c("fruit_stipe_length/FRSTlng", "fruit_length/FRUTlng",
           "seed_length/SEEDlng")

  traits <- c(veg, flo, fru)
  df <- data.frame(taxon = taxon, stringsAsFactors = FALSE)
  for (i in seq_along(traits)) {
    df[[traits[i]]] <- 10 + i + 3 * shift * (i %% 3 + 1) +
      stats::rnorm(n, sd = 1 + i / 10)
  }

  if (with_na) {
    df[c(2, 15), fru[1]] <- NA
    df[c(5, 30), flo[2]] <- NA
  }

  rownames(df) <- sprintf("ID%04d", seq_len(n))
  attr(df, "base_cols") <- list(veg = veg, flo = flo, fru = fru)
  if (colors) {
    attr(df, "taxon_colors") <- stats::setNames(
      c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")[seq_along(taxa)],
      taxa
    )
  }
  df
}

# Keep attributes when subsetting rows of a morph matrix
subset_rows <- function(df, keep) {
  out <- df[keep, , drop = FALSE]
  attr(out, "base_cols") <- attr(df, "base_cols")
  attr(out, "taxon_colors") <- attr(df, "taxon_colors")
  out
}

#_______________________________________________________________________________
# Raw spreadsheet (input of morph_matrix_setting / geo_matrix_setting) ####

write_raw_xlsx <- function(path = "raw_data.xlsx", with_id = FALSE) {
  set.seed(42)
  n <- 30
  sp <- rep(c("alpha", "beta", "gamma"), each = 10)
  shift <- rep(1:3, each = 10)
  raw <- data.frame(
    occurrenceID = paste0("occ", seq_len(n)),
    species = sp,
    decimalLatitude = as.character(round(-20 + shift * 2 + stats::rnorm(n), 3)),
    decimalLongitude = round(-50 + shift + stats::rnorm(n), 3),
    `petiole_length/PETIlng` = round(10 + shift + stats::rnorm(n), 2),
    `petiole_width/PETIwid` = as.character(round(2 + shift + stats::rnorm(n), 2)),
    `leaf_length/LEAFlng` = round(50 + 5 * shift + stats::rnorm(n), 2),
    `empty_trait/EMPTlng` = NA_real_,
    `inflorescence_length/INFLlng` = round(30 + shift + stats::rnorm(n), 2),
    `calyx_length/CLYXlng` = round(5 + shift + stats::rnorm(n), 2),
    `fruit_stipe_length/FRSTlng` = round(3 + shift + stats::rnorm(n), 2),
    `fruit_length/FRUTlng` = round(20 + 2 * shift + stats::rnorm(n), 2),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  # Decimal commas and blank strings must be coerced to numbers/NA
  raw$`petiole_width/PETIwid`[1] <- "3,5"
  raw$`petiole_width/PETIwid`[2] <- ""
  # A specimen without taxon and one without any trait
  raw$species[3] <- NA
  trait_cols <- grep("/", names(raw), value = TRUE)
  raw[4, trait_cols] <- NA
  # Missing coordinates
  raw$decimalLongitude[5] <- NA
  if (with_id) raw$id <- sprintf("SPEC%03d", seq_len(n))
  openxlsx::write.xlsx(raw, path)
  invisible(path)
}

#_______________________________________________________________________________
# Environmental data (as returned by geo_matrix_setting()) ####

clim_layer_names <- function() c(paste0("bio_", 1:19), "elevation")

make_clim_stack <- function(names = clim_layer_names(), seed = 99,
                            xmin = -70, xmax = -30, ymin = -35, ymax = 5,
                            res = 0.5) {
  set.seed(seed)
  r <- terra::rast(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                   resolution = res, crs = "EPSG:4326",
                   nlyrs = length(names))
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  vals <- vapply(seq_along(names), function(i) {
    a <- stats::runif(1, -1, 1)
    b <- stats::runif(1, -1, 1)
    100 + i * 10 + a * xy[, 1] * i + b * xy[, 2] * (20 - i) +
      stats::rnorm(nrow(xy), sd = 5 * i)
  }, numeric(nrow(xy)))
  terra::values(r) <- vals
  names(r) <- names
  r
}

make_soil_stack <- function(vars = c("clay", "sand"), seed = 7) {
  make_clim_stack(names = vars, seed = seed)
}

make_geodata <- function(analysis_data = make_morph_matrix(),
                         edaphic = FALSE, seed = 3) {
  set.seed(seed)
  ids <- rownames(analysis_data)
  taxon <- analysis_data$taxon
  shift <- as.numeric(factor(taxon))
  geo <- data.frame(
    taxon = taxon,
    decimalLatitude = -25 + 5 * shift + stats::rnorm(length(ids), sd = 1.5),
    decimalLongitude = -55 + 4 * shift + stats::rnorm(length(ids), sd = 1.5),
    stringsAsFactors = FALSE
  )
  clim <- make_clim_stack()
  pts <- terra::vect(geo, geom = c("decimalLongitude", "decimalLatitude"),
                     crs = "EPSG:4326")
  clim_vals <- terra::extract(clim, pts)[, -1, drop = FALSE]
  geo <- cbind(geo, clim_vals)
  env_cols <- list(clim = names(clim_vals))
  soil <- NULL
  if (edaphic) {
    soil <- make_soil_stack()
    soil_vals <- terra::extract(soil, pts)[, -1, drop = FALSE]
    geo <- cbind(geo, soil_vals)
    env_cols$edaphic <- names(soil_vals)
  }
  rownames(geo) <- ids
  attr(geo, "env_cols") <- env_cols
  attr(geo, "clim_stack") <- clim
  attr(geo, "soil_stack") <- soil
  attr(geo, "pts") <- pts
  geo
}

# Minimal sf world used to mock rnaturalearth downloads
fake_world_sf <- function(...) {
  poly <- sf::st_polygon(list(rbind(
    c(-80, -40), c(-30, -40), c(-30, 10), c(-80, 10), c(-80, -40)
  )))
  sf::st_sf(name = "Fakeland", geometry = sf::st_sfc(poly, crs = 4326))
}

fake_rivers_sf <- function(...) {
  line <- sf::st_linestring(rbind(c(-60, -20), c(-45, -10), c(-40, -5)))
  sf::st_sf(name = "Fake river", geometry = sf::st_sfc(line, crs = 4326))
}
