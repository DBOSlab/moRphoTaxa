# Tests for niche modelling functions: geo_enm() and geo_niche_overlap()
#
# flexsdm model fitting is mocked (a full ensemble takes minutes per taxon),
# so these tests exercise the geo_enm() workflow around the models. Natural
# Earth downloads are mocked as well. geo_niche_overlap() runs ecospat for
# real on a small grid with few permutations.

mock_natural_earth <- function(env = parent.frame()) {
  local_mocked_bindings(
    ne_countries = fake_world_sf,
    ne_states = fake_world_sf,
    .package = "rnaturalearth",
    .env = env
  )
}

# `fail`: names of flexsdm functions that should throw an error.
# `perf_col`: name of the performance column returned by fit_ensemble().
# `varimp`: "ok", "bad_cols" or "empty" output of sdm_varimp().
mock_flexsdm <- function(env = parent.frame(), fail = character(0),
                         perf_col = "SORENSEN_mean", varimp = "ok") {
  maybe_fail <- function(name) {
    if (name %in% fail) stop(name, " exploded", call. = FALSE)
  }
  fit <- function(name) {
    force(name)
    function(...) {
      maybe_fail(paste0("fit_", name))
      list(model = name)
    }
  }
  local_mocked_bindings(
    calib_area = function(data, x, y, method, crs) {
      maybe_fail("calib_area")
      e <- terra::ext(min(data[[x]]) - 2, max(data[[x]]) + 2,
                      min(data[[y]]) - 2, max(data[[y]]) + 2)
      terra::as.polygons(e, crs = crs)
    },
    sample_pseudoabs = function(data, x, y, n, method, rlayer, calibarea) {
      maybe_fail("sample_pseudoabs")
      s <- terra::spatSample(rlayer, n, na.rm = TRUE, xy = TRUE, method = "random")
      data.frame(x = s$x, y = s$y, pr_ab = 0)
    },
    part_random = function(data, pr_ab, method) {
      maybe_fail("part_random")
      data$.part <- rep_len(seq_len(as.numeric(method[["folds"]])), nrow(data))
      data
    },
    fit_glm = fit("glm"),
    fit_gam = fit("gam"),
    fit_gbm = fit("gbm"),
    fit_net = fit("net"),
    fit_raf = fit("raf"),
    fit_ensemble = function(models, ...) {
      maybe_fail("fit_ensemble")
      perf <- data.frame(model = "mean", value = 0.87)
      names(perf)[2] <- perf_col
      list(models = models, performance = perf)
    },
    sdm_predict = function(models, pred, con_thr) {
      maybe_fail("sdm_predict")
      r <- pred[[1]]
      mm <- terra::minmax(r)
      list(mean = (r - mm[1]) / (mm[2] - mm[1]))
    },
    sdm_varimp = function(models, data, response, predictors, ...) {
      maybe_fail("sdm_varimp")
      switch(varimp,
             ok = data.frame(model = rep(c("a", "b"), each = length(predictors)),
                             predictors = rep(predictors, 2),
                             SORENSEN = seq_len(2 * length(predictors)) / 10),
             bad_cols = data.frame(var = predictors, imp = 1),
             empty = data.frame())
    },
    .package = "flexsdm",
    .env = env
  )
}

#_______________________________________________________________________________
# geo_enm() ####

test_that("geo_enm models every taxon and saves maps, rasters and tables", {
  local_test_dir()
  mock_natural_earth()
  mock_flexsdm()
  geo <- make_geodata(edaphic = TRUE)

  msgs <- testthat::capture_messages(
    res <- geo_enm(geo, show_country_states = "Brazil",
                   file_formats = c("pdf", "jpg"))
  )
  expect_true(any(grepl("Environmental source: clim_edaphic", msgs)))
  expect_true(any(grepl("Dropping .* constant/collinear predictors", msgs)))
  expect_named(res, c("summary", "suitability", "maps", "predictors", "output_dir"))
  expect_equal(res$summary$taxon, c("alpha", "beta", "gamma"))
  expect_true(all(res$summary$n_algos == 5))
  expect_equal(res$summary$sorensen, rep(0.87, 3))
  expect_s4_class(res$suitability$alpha, "SpatRaster")
  expect_gte(length(res$predictors), 2)

  files <- list.files(res$output_dir)
  expect_true(all(c("ENM_summary.xlsx", "ENM_alpha.xlsx", "ENM_map_alpha.pdf",
                    "ENM_map_alpha.jpg", "ENM_importance_alpha.pdf",
                    "ENM_suitability_alpha.tif") %in% files))
})

test_that("geo_enm runs without saving, inset or soil data", {
  local_test_dir()
  mock_natural_earth()
  mock_flexsdm(fail = "fit_gbm")
  geo <- make_geodata()
  # A constant layer is removed before the collinearity check
  clim <- attr(geo, "clim_stack")
  const <- clim[[1]]
  terra::values(const) <- 1
  names(const) <- "constant"
  attr(geo, "clim_stack") <- c(clim, const)

  msgs <- testthat::capture_messages(
    res <- geo_enm(geo, inset = FALSE, save = FALSE, seed = NULL)
  )
  expect_true(any(grepl("Environmental source: clim$", trimws(msgs))))
  expect_true(any(grepl("Algorithms that failed: gbm", msgs)))
  expect_true(any(grepl("Outputs were not saved", msgs)))
  expect_false("constant" %in% res$predictors)
  expect_null(res$output_dir)
  expect_true(all(res$summary$n_algos == 4))
  expect_false(dir.exists("Figs_ENM"))
})

test_that("geo_enm resamples soil rasters with a different geometry", {
  local_test_dir()
  mock_natural_earth()
  mock_flexsdm()
  geo <- make_geodata(edaphic = TRUE)
  attr(geo, "soil_stack") <- make_clim_stack(names = c("clay", "sand"), res = 1)

  res <- geo_enm(geo, save = FALSE, verbose = FALSE)
  expect_equal(nrow(res$summary), 3)
})

test_that("geo_enm skips taxa when a workflow step fails", {
  local_test_dir()
  mock_natural_earth()
  geo <- make_geodata()

  steps <- c(calib_area = "Calibration area failed",
             sample_pseudoabs = "Pseudo-absence sampling failed",
             part_random = "Data partitioning failed",
             fit_ensemble = "Ensemble fitting failed",
             sdm_predict = "Prediction failed")
  for (step in names(steps)) {
    local({
      mock_flexsdm(fail = step)
      w <- testthat::capture_warnings(
        res <- geo_enm(geo, save = FALSE, inset = FALSE, verbose = FALSE)
      )
      expect_true(any(grepl(steps[[step]], w)), info = step)
      expect_true(any(grepl("No taxon could be modelled", w)), info = step)
      expect_null(res$summary)
      expect_length(res$maps, 0)
    })
  }
})

test_that("geo_enm needs at least two fitted algorithms", {
  local_test_dir()
  mock_natural_earth()
  mock_flexsdm(fail = c("fit_glm", "fit_gam", "fit_gbm", "fit_net"))
  w <- testthat::capture_warnings(
    res <- geo_enm(make_geodata(), save = FALSE, inset = FALSE, verbose = FALSE)
  )
  expect_true(any(grepl("Fewer than 2 algorithms converged", w)))
  expect_null(res$summary)
})

test_that("geo_enm tolerates missing Sorensen and importance outputs", {
  local_test_dir()
  mock_natural_earth()
  mock_flexsdm(perf_col = "TSS_mean", varimp = "bad_cols")
  w <- testthat::capture_warnings(
    res <- geo_enm(make_geodata(), inset = FALSE, file_formats = "pdf",
                   verbose = FALSE)
  )
  expect_true(any(grepl("No SORENSEN-like column", w)))
  expect_true(any(grepl("Expected columns not found in sdm_varimp", w)))
  expect_true(any(grepl("Permutation importance could not be computed", w)))
  expect_true(all(is.na(res$summary$sorensen)))
  expect_true(file.exists(file.path(res$output_dir, "ENM_alpha.xlsx")))
})

test_that("geo_enm reports sdm_varimp failures and empty importance tables", {
  local_test_dir()
  mock_natural_earth()
  local({
    mock_flexsdm(fail = "sdm_varimp")
    w <- testthat::capture_warnings(
      geo_enm(make_geodata(), save = FALSE, inset = FALSE, verbose = FALSE)
    )
    expect_true(any(grepl("sdm_varimp failed", w)))
  })
  local({
    mock_flexsdm(varimp = "empty")
    w <- testthat::capture_warnings(
      geo_enm(make_geodata(), save = FALSE, inset = FALSE, verbose = FALSE)
    )
    expect_true(any(grepl("Permutation importance could not be computed", w)))
  })
})

test_that("geo_enm skips taxa with too few records", {
  local_test_dir()
  mock_natural_earth()
  mock_flexsdm()
  w <- testthat::capture_warnings(
    res <- geo_enm(make_geodata(), min_occ = 13, save = FALSE, inset = FALSE,
                   verbose = FALSE)
  )
  expect_true(any(grepl("fewer than 13 occurrence records", w)))
  expect_null(res$summary)
})

test_that("geo_enm validates its inputs", {
  local_test_dir()
  mock_natural_earth()
  mock_flexsdm()
  geo <- make_geodata()

  expect_error(geo_enm("x"), "must be the data frame")
  expect_error(geo_enm(geo[, c("taxon", "bio_1")]), "decimalLatitude")
  no_stack <- geo
  attr(no_stack, "clim_stack") <- NULL
  expect_error(geo_enm(no_stack), "clim_stack")
  expect_error(geo_enm(geo, buffer_km = 0), "buffer_km")
  expect_error(geo_enm(geo, pa_ratio = -1), "pa_ratio")
  expect_error(geo_enm(geo, folds = 1), "folds")
  expect_error(geo_enm(geo, min_occ = 1), "min_occ")
  expect_error(geo_enm(geo, importance_perm = 0), "importance_perm")
  expect_error(geo_enm(geo, seed = "a"), "seed")
  expect_error(geo_enm(geo, show_country_states = 1), "show_country_states")
  expect_error(geo_enm(geo, inset = "a"), "inset")
  expect_error(geo_enm(geo, save = "a"), "save")
  expect_error(geo_enm(geo, file_formats = character(0)), "file_formats")
  expect_error(geo_enm(geo, file_formats = "png"))
  expect_error(geo_enm(geo, verbose = "a"), "verbose")

  collinear <- geo
  r <- attr(geo, "clim_stack")[[1]]
  r2 <- r * 2
  names(r2) <- "bio_1_copy"
  attr(collinear, "clim_stack") <- c(r, r2)
  expect_error(geo_enm(collinear, save = FALSE, verbose = FALSE),
               "Fewer than 2 predictors remain")

  no_taxa <- geo
  no_taxa$taxon <- " "
  expect_error(geo_enm(no_taxa, save = FALSE, inset = FALSE, verbose = FALSE),
               "No valid taxa")
})

test_that(".enm_inset_limits pads and clamps the study box", {
  lim <- .enm_inset_limits(c(xmin = -50, xmax = -40, ymin = -20, ymax = -10))
  expect_equal(unname(lim), c(-65, -25, -35, 5))
  big <- .enm_inset_limits(c(xmin = -170, xmax = 170, ymin = -80, ymax = 80))
  expect_equal(unname(big), c(-180, 180, -90, 90))
})

test_that("ENM wrappers turn errors into NULL", {
  expect_null(.try_fit(stop("x")))
  expect_equal(.try_fit(1), 1)
  expect_warning(out <- .safe_run(stop("bad"), "Step failed"), "Step failed: bad")
  expect_null(out)
  expect_equal(.safe_run(2, "never"), 2)
})

#_______________________________________________________________________________
# geo_niche_overlap() ####

test_that("geo_niche_overlap compares every pair of taxa", {
  local_test_dir()
  geo <- make_geodata()

  w <- testthat::capture_warnings(
    msgs <- testthat::capture_messages(
      res <- geo_niche_overlap(geo, grid_size = 15, iterations = 9,
                               file_formats = c("pdf", "jpg"))
    )
  )
  expect_true(any(grepl("fewer than 999 iterations", w)))
  expect_true(any(grepl("Pairwise comparisons: 3", msgs)))

  expect_named(res, c("results", "D_matrix", "pca", "variables",
                      "removed_variables", "plots", "output_dir"))
  expect_equal(nrow(res$results), 3)
  expect_true(all(res$results$schoener_D >= 0 & res$results$schoener_D <= 1))
  expect_true(all(c("equiv_p", "similarity_p", "expansion", "stability",
                    "unfilling", "variables_used") %in% names(res$results)))
  expect_equal(dim(res$D_matrix), c(3, 3))
  expect_true(isSymmetric(unname(res$D_matrix[!is.na(res$D_matrix[, 1]),
                                              !is.na(res$D_matrix[1, ])])))
  expect_s3_class(res$plots$heatmap, "ggplot")
  expect_length(res$plots$tests, 3)

  files <- list.files(res$output_dir)
  expect_true(all(c("NicheOverlap_all_results.xlsx", "NicheOverlap_heatmap.pdf",
                    "NicheSpace_combined.jpg", "EnvPCA_scree_combined.pdf",
                    "Environmental_PCA_loadings.xlsx",
                    "Environmental_PCA_variance.xlsx",
                    "NicheTest_alphavsbeta.pdf") %in% files))
})

test_that("geo_niche_overlap handles edaphic data, collinearity and no saving", {
  local_test_dir()
  geo <- make_geodata(edaphic = TRUE)
  geo$bio_1_copy <- geo$bio_1 * 3
  geo$constant <- 7
  attr(geo, "env_cols")$clim <- c(attr(geo, "env_cols")$clim, "bio_1_copy", "constant")

  msgs <- testthat::capture_messages(
    res <- suppressWarnings(
      geo_niche_overlap(geo, cor_cutoff = 0.99, grid_size = 10, iterations = 5,
                        seed = NULL, save = FALSE)
    )
  )
  expect_true(any(grepl("climate \\+ edaphic", msgs)))
  expect_true(any(grepl("Dropping .* collinear variable", msgs)))
  expect_true(any(grepl("Outputs were not saved", msgs)))
  expect_true("constant" %in% res$removed_variables)
  expect_null(res$output_dir)
  expect_false(dir.exists("Figs_NicheOverlap"))
})

test_that("geo_niche_overlap saves the removed-variable table when saving", {
  local_test_dir()
  geo <- make_geodata()
  geo$bio_1_copy <- geo$bio_1 * 3
  attr(geo, "env_cols")$clim <- c(attr(geo, "env_cols")$clim, "bio_1_copy")
  res <- suppressWarnings(
    geo_niche_overlap(geo, cor_cutoff = 0.99, grid_size = 10, iterations = 5,
                      file_formats = "pdf", verbose = FALSE)
  )
  expect_true(file.exists(file.path(res$output_dir,
                                    "Removed_correlated_variables.xlsx")))
})

test_that("geo_niche_overlap skips pairs and reports ecospat failures", {
  local_test_dir()
  geo <- make_geodata()

  expect_error(
    suppressWarnings(
      geo_niche_overlap(geo, min_occ = 13, grid_size = 10, iterations = 5,
                        save = FALSE, verbose = FALSE)
    ),
    "No taxon pairs produced valid niche overlap results"
  )

  local({
    local_mocked_bindings(ecospat.grid.clim.dyn = function(...) stop("grid"),
                          .package = "ecospat")
    w <- testthat::capture_warnings(
      expect_error(
        geo_niche_overlap(geo, grid_size = 10, iterations = 5, save = FALSE,
                          verbose = FALSE),
        "No taxon pairs"
      )
    )
    expect_true(any(grepl("Grid construction failed for alpha: grid", w)))
  })

  local({
    local_mocked_bindings(ecospat.niche.equivalency.test = function(...) stop("eq"),
                          .package = "ecospat")
    w <- testthat::capture_warnings(
      expect_error(
        geo_niche_overlap(geo, grid_size = 10, iterations = 5, save = FALSE,
                          verbose = FALSE),
        "No taxon pairs"
      )
    )
    expect_true(any(grepl("Equivalency test failed: eq", w)))
  })

  local({
    # ecospat's own permutation tests call ecospat.niche.dyn.index() with
    # intersection = 0; only the direct call from geo_niche_overlap() fails.
    orig_dyn <- ecospat::ecospat.niche.dyn.index
    local_mocked_bindings(
      ecospat.niche.dyn.index = function(z1, z2, intersection = 0) {
        if (identical(intersection, 0.1)) stop("dyn")
        orig_dyn(z1, z2, intersection)
      },
      .package = "ecospat"
    )
    res <- suppressWarnings(
      geo_niche_overlap(geo, grid_size = 10, iterations = 5, save = FALSE,
                        verbose = TRUE) |>
        suppressMessages()
    )
    expect_true(all(is.na(res$results$expansion)))
  })
})

test_that("geo_niche_overlap validates its inputs", {
  local_test_dir()
  geo <- make_geodata()

  expect_error(geo_niche_overlap("x"), "must be the data frame")
  no_taxon <- geo
  no_taxon$taxon <- NULL
  expect_error(geo_niche_overlap(no_taxon), "taxon")
  no_env <- geo
  attr(no_env, "env_cols") <- NULL
  expect_error(geo_niche_overlap(no_env), "env_cols")
  expect_error(geo_niche_overlap(geo, var_threshold = 0), "var_threshold")
  expect_error(geo_niche_overlap(geo, cor_cutoff = 2), "cor_cutoff")
  expect_error(geo_niche_overlap(geo, grid_size = 5), "grid_size")
  expect_error(geo_niche_overlap(geo, iterations = 0), "iterations")
  expect_error(geo_niche_overlap(geo, alpha = 1), "alpha")
  expect_error(geo_niche_overlap(geo, min_occ = 1), "min_occ")
  expect_error(geo_niche_overlap(geo, n_cores = 0), "n_cores")
  expect_error(geo_niche_overlap(geo, seed = "a"), "seed")
  expect_error(geo_niche_overlap(geo, save = "a"), "save")
  expect_error(geo_niche_overlap(geo, file_formats = character(0)), "file_formats")
  expect_error(geo_niche_overlap(geo, file_formats = "png"))
  expect_error(geo_niche_overlap(geo, verbose = "a"), "verbose")

  one_var <- geo
  attr(one_var, "env_cols") <- list(clim = "bio_1")
  expect_error(suppressWarnings(geo_niche_overlap(one_var, verbose = FALSE)),
               "Fewer than 2 environmental variables")

  few <- geo
  few[5:nrow(few), "bio_1"] <- NA
  expect_error(suppressWarnings(geo_niche_overlap(few, verbose = FALSE)),
               "Fewer than 20 complete environmental records")

  flat <- geo
  attr(flat, "env_cols") <- list(clim = c("bio_1", "bio_2"))
  flat$bio_2 <- 1
  expect_error(suppressWarnings(geo_niche_overlap(flat, verbose = FALSE)),
               "Fewer than 2 variables with variation")

  dup <- geo
  dup$bio_1_copy <- dup$bio_1 * 2
  attr(dup, "env_cols") <- list(clim = c("bio_1", "bio_1_copy"))
  expect_error(suppressWarnings(geo_niche_overlap(dup, save = FALSE, verbose = FALSE)),
               "Fewer than 2 variables remain after collinearity filtering")

  one_taxon <- geo
  one_taxon$taxon <- "alpha"
  expect_error(suppressWarnings(geo_niche_overlap(one_taxon, save = FALSE,
                                                  verbose = FALSE)),
               "Need at least 2 taxa")
})
