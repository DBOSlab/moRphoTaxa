# Tests for analyses combining morphology and geography:
# geo_morph_mantel(), morph_autocorrelation(), geo_morph_rda()

#_______________________________________________________________________________
# geo_morph_mantel() ####

test_that("geo_morph_mantel correlates morphological and geographic distances", {
  local_test_dir()
  m <- make_morph_matrix()
  geo <- make_geodata(m)

  expect_message(
    res <- geo_morph_mantel(geo, m, permutations = 99,
                            file_formats = c("pdf", "jpg")),
    "Mantel test results"
  )
  expect_named(res, c("mantel", "dist_geo", "dist_morpho", "plot", "n", "output_dir"))
  expect_s3_class(res$mantel, "mantel")
  expect_equal(res$n, 36)
  expect_s3_class(res$dist_geo, "dist")
  expect_s3_class(res$plot, "ggplot")
  # Taxa are geographically and morphologically structured
  expect_gt(res$mantel$statistic, 0)
  expect_true(all(file.exists(file.path(res$output_dir,
                                        c("Mantel_plot.pdf", "Mantel_plot.jpg")))))
})

test_that("geo_morph_mantel supports other distances and drops incomplete specimens", {
  local_test_dir()
  m <- make_morph_matrix()
  geo <- make_geodata(m)
  geo[1, "decimalLatitude"] <- NA

  expect_warning(
    res <- geo_morph_mantel(geo, m, cor_method = "spearman", permutations = 99,
                            dist_method = "euclidean", geo_dist_method = "euclidean",
                            file_formats = "pdf", verbose = FALSE),
    "specimens with NAs dropped for method 'euclidean'"
  )
  # 1 without coordinates + 4 with missing traits
  expect_equal(res$n, 31)
  expect_equal(res$plot$labels$x, "Geographic distance (decimal degrees)")
})

test_that("geo_morph_mantel validates its inputs", {
  local_test_dir()
  m <- make_morph_matrix()
  geo <- make_geodata(m)
  expect_error(geo_morph_mantel("x", m), "`geodata` must be a data frame")
  expect_error(geo_morph_mantel(geo, "x"), "`analysis_data` must be a data frame")
  expect_error(geo_morph_mantel(geo[, "taxon", drop = FALSE], m), "decimalLatitude")
  expect_error(geo_morph_mantel(geo, m, cor_method = "x"))
  expect_error(geo_morph_mantel(geo, m, permutations = 1.5), "permutations")
  expect_error(geo_morph_mantel(geo, m, file_formats = character(0)), "file_formats")
  expect_error(geo_morph_mantel(geo, m, verbose = 1), "verbose")

  other <- geo
  rownames(other) <- paste0("X", seq_len(nrow(other)))
  expect_error(geo_morph_mantel(other, m), "No specimen identifiers are shared")

  no_num <- m[, "taxon", drop = FALSE]
  expect_error(geo_morph_mantel(geo, no_num), "no numeric trait columns")

  few <- geo
  few$decimalLatitude[3:nrow(few)] <- NA
  expect_error(geo_morph_mantel(few, m, verbose = FALSE), "Fewer than 3 specimens")

  gappy <- m
  gappy[1, -1] <- NA
  gappy[1, "petiole_length/PETIlng"] <- 5
  gappy[2, -1] <- NA
  gappy[2, "seed_length/SEEDlng"] <- 5
  expect_error(geo_morph_mantel(geo, gappy, verbose = FALSE),
               "distance matrix contains NA values")
})

test_that(".haversine_dist returns great-circle distances in km", {
  d <- .haversine_dist(lon = c(0, 0, 90), lat = c(0, 1, 0), labels = c("a", "b", "c"))
  expect_s3_class(d, "dist")
  dm <- as.matrix(d)
  expect_equal(dm["a", "b"], 111.19, tolerance = 0.01)
  expect_equal(dm["a", "c"], pi / 2 * 6371.0088, tolerance = 1e-6)
  expect_null(attr(.haversine_dist(0, 0), "Labels"))
})

#_______________________________________________________________________________
# morph_autocorrelation() ####

test_that("morph_autocorrelation computes lag-1 autocorrelation along latitude", {
  local_test_dir()
  m <- make_morph_matrix()
  geo <- make_geodata(m)
  geo[1, "decimalLatitude"] <- NA

  msgs <- testthat::capture_messages(
    res <- morph_autocorrelation(m, geo, file_formats = c("pdf", "jpg"))
  )
  expect_true(any(grepl("without coordinates in `geodata` dropped: 1", msgs)))
  expect_named(res, c("results", "plots", "specimen_order", "output_dir"))
  expect_setequal(unique(res$results$block),
                  c("veg", "flo", "fru", "vegflo", "vegfru", "flofru", "vegflofru"))
  expect_true(all(abs(res$results$autocor) <= 1))
  expect_length(res$specimen_order, 35)
  lat <- geo[res$specimen_order, "decimalLatitude"]
  expect_false(is.unsorted(lat))
  expect_true(file.exists(file.path(res$output_dir, "Autocorrelation_summary.csv")))
  expect_true(file.exists(file.path(res$output_dir, "Autocorrelation_veg.jpg")))
})

test_that("morph_autocorrelation supports longitude, no ordering and block subsets", {
  local_test_dir()
  m <- make_morph_matrix()
  geo <- make_geodata(m)

  lon <- morph_autocorrelation(m, geo, blocks = "vegflo", order_by = "longitude",
                               cor_method = "spearman", file_formats = "pdf",
                               verbose = FALSE)
  expect_equal(unique(lon$results$block), "vegflo")
  expect_false(is.unsorted(geo[lon$specimen_order, "decimalLongitude"]))

  expect_warning(
    none <- morph_autocorrelation(m, blocks = c("veg", "zzz"), order_by = "none",
                                  file_formats = "pdf", verbose = FALSE),
    "Ignoring unavailable blocks: zzz"
  )
  expect_equal(none$specimen_order, rownames(m))
})

test_that("morph_autocorrelation skips blocks without usable traits", {
  local_test_dir()
  m <- make_morph_matrix()
  m$`leaf_length/LEAFlng` <- as.character(m$`leaf_length/LEAFlng`)
  attr(m, "base_cols") <- list(veg = "leaf_length/LEAFlng",
                               flo = attr(m, "base_cols")$flo)
  expect_warning(
    res <- morph_autocorrelation(m, blocks = c("veg", "flo"), order_by = "none",
                                 file_formats = "pdf", verbose = FALSE),
    "Skipping 'veg': no valid columns"
  )
  expect_equal(unique(res$results$block), "flo")

  sparse <- make_morph_matrix()
  sparse[seq(1, 36, by = 2), attr(sparse, "base_cols")$fru] <- NA
  attr(sparse, "base_cols") <- list(fru = attr(sparse, "base_cols")$fru)
  w <- testthat::capture_warnings(
    res2 <- morph_autocorrelation(sparse, order_by = "none", file_formats = "pdf",
                                  verbose = FALSE)
  )
  expect_true(any(grepl("no traits with enough data", w)))
  expect_true(any(grepl("No block produced results", w)))
  expect_null(res2$results)
})

test_that("morph_autocorrelation validates its inputs", {
  local_test_dir()
  m <- make_morph_matrix()
  geo <- make_geodata(m)

  expect_error(morph_autocorrelation(m), "`geodata` .* is required")
  expect_error(morph_autocorrelation(m, geo[, "taxon", drop = FALSE]),
               "decimalLatitude")
  expect_error(morph_autocorrelation(m, geo, cor_method = "x"))
  expect_error(morph_autocorrelation(m, geo, order_by = "x"))

  no_blocks <- m
  attr(no_blocks, "base_cols") <- list(other = "leaf_length/LEAFlng")
  expect_error(morph_autocorrelation(no_blocks, order_by = "none"),
               "none of the blocks")

  other <- geo
  rownames(other) <- paste0("X", seq_len(nrow(other)))
  expect_error(morph_autocorrelation(m, other), "No specimen identifiers are shared")

  small <- subset_rows(m, 1:3)
  expect_error(morph_autocorrelation(small, order_by = "none"), "At least 4 specimens")
})

test_that(".lag1_autocor handles short and sparse series", {
  expect_true(is.na(.lag1_autocor(1, "pearson")))
  expect_true(is.na(.lag1_autocor(c(1, NA, 3, NA), "pearson")))
  expect_equal(.lag1_autocor(1:10, "pearson"), 1)
})

#_______________________________________________________________________________
# geo_morph_rda() ####

test_that("geo_morph_rda constrains trait variation by climate", {
  local_test_dir()
  m <- make_morph_matrix()
  geo <- make_geodata(m)

  msgs <- testthat::capture_messages(
    res <- geo_morph_rda(m, geo, blocks = c("veg", "vegflofru"),
                         permutations = 99, file_formats = c("pdf", "jpg"))
  )
  expect_true(any(grepl("Auto-selected 19 bioclimatic variables", msgs)))
  expect_named(res, c("results", "models", "plots", "output_dir"))
  expect_equal(res$results$block, c("veg", "vegflofru"))
  expect_s3_class(res$models$veg, "rda")
  expect_s3_class(res$plots$veg$biplot, "ggplot")
  expect_true(all(res$results$constrained_var > 0 & res$results$constrained_var < 100))

  files <- list.files(res$output_dir)
  expect_true(all(c("RDA_summary.xlsx", "RDA_veg.xlsx", "RDA_biplot_veg.pdf",
                    "RDA_biplot_veg.jpg", "RDA_terms_veg.pdf") %in% files))
  expect_true("by_term" %in% openxlsx::getSheetNames(file.path(res$output_dir,
                                                               "RDA_veg.xlsx")))
})

test_that("geo_morph_rda accepts user variables, unscaled data and collinearity filtering", {
  local_test_dir()
  m <- make_morph_matrix(colors = FALSE)
  geo <- make_geodata(m)
  geo$bio_1_copy <- geo$bio_1 * 2 + 1

  msgs <- testthat::capture_messages(
    expect_warning(
      res <- geo_morph_rda(m, geo, blocks = c("flo", "nope"),
                           env_vars = c("bio_1", "bio_1_copy", "bio_12", "bio_15"),
                           permutations = 99, scale_data = FALSE,
                           file_formats = "pdf"),
      "Ignoring unavailable blocks: nope"
    )
  )
  # The earlier variable of a collinear pair is the one dropped
  expect_true(any(grepl("Dropping collinear env vars: bio_1\\s*$", msgs)))
  expect_equal(res$results$n_env, 3)
})

test_that("geo_morph_rda skips blocks it cannot analyse", {
  local_test_dir()
  m <- make_morph_matrix()
  geo <- make_geodata(m)

  small <- subset_rows(m, 1:8)
  w <- testthat::capture_warnings(
    res <- geo_morph_rda(small, geo, blocks = "veg", permutations = 99,
                         env_vars = c("bio_1", "bio_12"),
                         file_formats = "pdf", verbose = FALSE)
  )
  expect_true(any(grepl("Skipping 'veg': not enough data", w)))
  expect_true(any(grepl("No block was analyzed", w)))
  expect_null(res$results)

  # 12 specimens cannot support 19 environmental predictors
  mid <- subset_rows(m, 1:12)
  w_mid <- testthat::capture_warnings(
    geo_morph_rda(mid, geo, blocks = "veg", permutations = 99,
                  file_formats = "pdf", verbose = FALSE)
  )
  expect_true(any(grepl("must be smaller than the number of specimens minus one",
                        w_mid)))

  chr <- m
  chr$`leaf_length/LEAFlng` <- as.character(chr$`leaf_length/LEAFlng`)
  attr(chr, "base_cols") <- list(veg = "leaf_length/LEAFlng")
  w_chr <- testthat::capture_warnings(
    geo_morph_rda(chr, geo, permutations = 99, env_vars = c("bio_1", "bio_12"),
                  file_formats = "pdf", verbose = FALSE)
  )
  expect_true(any(grepl("Skipping 'veg': no valid columns", w_chr)))
})

test_that("geo_morph_rda reports failures from vegan", {
  local_test_dir()
  m <- make_morph_matrix()
  geo <- make_geodata(m)

  local_mocked_bindings(rda = function(...) stop("rda broke"), .package = "vegan")
  w <- testthat::capture_warnings(
    geo_morph_rda(m, geo, blocks = "veg", permutations = 99,
                  env_vars = c("bio_1", "bio_12", "bio_15"),
                  file_formats = "pdf", verbose = FALSE)
  )
  expect_true(any(grepl("RDA failed for 'veg': rda broke", w)))
})

test_that("geo_morph_rda validates its inputs", {
  local_test_dir()
  m <- make_morph_matrix()
  geo <- make_geodata(m)

  expect_error(geo_morph_rda("x", geo), "`analysis_data` must be a data frame")
  expect_error(geo_morph_rda(m, "x"), "`geodata` must be a data frame")
  no_taxon <- m
  no_taxon$taxon <- NULL
  expect_error(geo_morph_rda(no_taxon, geo), "taxon")
  no_bc <- m
  attr(no_bc, "base_cols") <- NULL
  expect_error(geo_morph_rda(no_bc, geo), "base_cols")
  expect_error(geo_morph_rda(m, geo, blocks = 1), "blocks")
  expect_error(geo_morph_rda(m, geo, env_vars = 1), "env_vars")
  expect_error(geo_morph_rda(m, geo, permutations = 0), "permutations")
  expect_error(geo_morph_rda(m, geo, scale_data = NA), "scale_data")
  expect_error(geo_morph_rda(m, geo, env_cor_threshold = 2), "env_cor_threshold")
  expect_error(geo_morph_rda(m, geo, file_formats = character(0)), "file_formats")
  expect_error(geo_morph_rda(m, geo, verbose = "a"), "verbose")
  expect_error(geo_morph_rda(m, geo, blocks = "zzz"), "No valid blocks")
  expect_error(geo_morph_rda(m, geo, env_vars = "nope"), "not found in `geodata`: nope")

  other_blocks <- m
  attr(other_blocks, "base_cols") <- list(other = "leaf_length/LEAFlng")
  expect_error(geo_morph_rda(other_blocks, geo), "none of the blocks")

  no_env <- geo
  attr(no_env, "env_cols") <- NULL
  expect_error(geo_morph_rda(m, no_env), "env_cols")
  only_elev <- geo
  attr(only_elev, "env_cols") <- list(clim = "elevation")
  expect_error(geo_morph_rda(m, only_elev, verbose = FALSE), "No bioclimatic columns")

  other_ids <- geo
  rownames(other_ids) <- paste0("X", seq_len(nrow(other_ids)))
  expect_error(geo_morph_rda(m, other_ids, verbose = FALSE),
               "No specimen identifiers are shared")
})

test_that(".short_env_label shortens bioclim names", {
  expect_equal(.short_env_label(c("bio_12", "clay")), c("bio12", "clay"))
})
