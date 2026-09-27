# Tests for ordination functions: morph_pca(), morph_pcoa(), morph_nmds()

#_______________________________________________________________________________
# morph_pca() ####

test_that("morph_pca runs every block and trait combination", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_message(
    out <- morph_pca(m, top_contrib = 5, file_formats = c("pdf", "jpeg")),
    "PCA analysis done"
  )
  expect_null(out)

  files <- list.files(out_dir("Figs.PCA"))
  for (b in c("veg", "flo", "fru", "vegflo", "vegfru", "flofru", "vegflofru")) {
    expect_true(paste0("PCA_", b, ".pdf") %in% files, info = b)
    expect_true(paste0("PCA_", b, ".jpeg") %in% files, info = b)
    expect_true(paste0("Dim1_contrib_variables_", b, ".pdf") %in% files, info = b)
    expect_true(paste0("Dim2_contrib_variables_", b, ".pdf") %in% files, info = b)
    expect_true(paste0("PCA_panels_", b, ".pdf") %in% files, info = b)
  }
})

test_that("morph_pca handles selected blocks, missing colours and sparse blocks", {
  local_test_dir()
  m <- make_morph_matrix(colors = FALSE)
  morph_pca(m, pca_blocks = "veg", file_formats = "pdf", verbose = FALSE)
  files <- list.files(out_dir("Figs.PCA"))
  expect_true("PCA_veg.pdf" %in% files)
  expect_false("PCA_flo.pdf" %in% files)

  m2 <- make_morph_matrix()
  attr(m2, "taxon_colors") <- c(alpha = "red", beta = "blue")
  expect_warning(
    morph_pca(m2, pca_blocks = "fru", file_formats = "pdf", verbose = FALSE),
    "No stored colour found for: gamma"
  )

  # Two traits: no PC1-PC3 panels; one trait: skipped
  m3 <- make_morph_matrix()
  attr(m3, "base_cols") <- list(veg = attr(m3, "base_cols")$veg[1:2],
                                flo = attr(m3, "base_cols")$flo[1])
  unlink("Figs.PCA", recursive = TRUE)
  expect_warning(
    morph_pca(m3, pca_blocks = c("veg", "flo"), file_formats = "pdf",
              verbose = FALSE),
    "Skipping 'flo': not enough data"
  )
  files <- list.files(out_dir("Figs.PCA"))
  expect_true("PCA_veg.pdf" %in% files)
  expect_false("PCA_panels_veg.pdf" %in% files)
})

test_that("morph_pca validates top_contrib", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_pca(m, top_contrib = 0), "top_contrib")
  expect_error(morph_pca(m, top_contrib = c(1, 2)), "top_contrib")
})

#_______________________________________________________________________________
# morph_pcoa() ####

test_that("morph_pcoa returns scores and eigenvalues for each block", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_message(
    res <- morph_pcoa(m, top_vars = 3, file_formats = c("pdf", "jpeg")),
    "PCoA analysis done"
  )
  expect_setequal(names(res),
                  c("veg", "flo", "fru", "vegflo", "vegfru", "flofru", "vegflofru"))
  veg <- res$veg
  expect_named(veg, c("scores", "eigenvalues", "n_specimens", "n_traits"))
  expect_true(all(c("PCo1", "PCo2", "taxon") %in% names(veg$scores)))
  expect_equal(veg$n_traits, 4)
  expect_equal(veg$n_specimens, 36)
  expect_equal(nrow(veg$scores), 36)

  files <- list.files(out_dir("Figs.PCoA"))
  expect_true(all(c("PCoA_veg.pdf", "PCoA_veg.jpeg", "PCoA_eigenvalues_veg.pdf",
                    "Dim1_corr_variables_veg.pdf", "Dim2_corr_variables_veg.pdf",
                    "PCoA_panels_veg.pdf") %in% files))
})

test_that("morph_pcoa supports euclidean distances and eigenvalue corrections", {
  local_test_dir()
  m <- make_morph_matrix(colors = FALSE)
  eu <- morph_pcoa(m, pcoa_blocks = "veg", dist_method = "euclidean",
                   file_formats = "pdf", verbose = FALSE)
  expect_named(eu, "veg")

  li <- morph_pcoa(m, pcoa_blocks = "flo", correction = "lingoes",
                   file_formats = "pdf", verbose = FALSE)
  expect_true(all(li$flo$eigenvalues$eigenvalue >= -1e-8))

  ca <- morph_pcoa(m, pcoa_blocks = "fru", correction = "cailliez",
                   file_formats = "pdf", verbose = FALSE)
  expect_named(ca, "fru")

  expect_error(morph_pcoa(m, dist_method = "bray"))
  expect_error(morph_pcoa(m, correction = "foo"))
  expect_error(morph_pcoa(m, top_vars = 0), "top_vars")
})

test_that("morph_pcoa reports negative eigenvalues for Gower without correction", {
  local_test_dir()
  m <- make_morph_matrix()
  msgs <- testthat::capture_messages(
    morph_pcoa(m, pcoa_blocks = "veg", file_formats = "pdf", verbose = TRUE)
  )
  expect_true(any(grepl("negative eigenvalue", msgs)))
})

test_that("combined blocks are only reachable through \"all\"", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_pcoa(m, pcoa_blocks = "vegflofru", verbose = FALSE),
               "No valid blocks")
})

test_that("morph_pcoa skips blocks without enough specimens or variable traits", {
  local_test_dir()
  m <- make_morph_matrix()
  m[, "seed_length/SEEDlng"] <- 5 # constant trait
  attr(m, "base_cols") <- list(
    veg = attr(m, "base_cols")$veg,
    flo = c("fruit_length/FRUTlng", "seed_length/SEEDlng"),
    fru = attr(m, "base_cols")$fru
  )
  m[3:36, "fruit_stipe_length/FRSTlng"] <- NA

  w <- testthat::capture_warnings(
    res <- morph_pcoa(m, pcoa_blocks = c("veg", "flo", "fru"),
                      file_formats = "pdf", verbose = FALSE)
  )
  expect_true(any(grepl("Skipping 'flo': not enough data", w)))
  expect_true(any(grepl("Skipping 'fru': fewer than 3 specimens", w)))
  expect_named(res, "veg")

  # Two-axis solution (no PCo1-PCo3 panels)
  attr(m, "base_cols") <- list(veg = attr(m, "base_cols")$veg[1:2])
  unlink("Figs.PCoA", recursive = TRUE)
  res2 <- morph_pcoa(m, pcoa_blocks = "veg", dist_method = "euclidean",
                     file_formats = "pdf", verbose = FALSE)
  expect_named(res2, "veg")

  # Missing colours are filled in
  attr(m, "taxon_colors") <- c(alpha = "red")
  expect_warning(
    morph_pcoa(m, pcoa_blocks = "veg", file_formats = "pdf", verbose = FALSE),
    "No stored colour found"
  )
})

#_______________________________________________________________________________
# morph_nmds() ####

test_that("morph_nmds returns a stress summary for each block", {
  local_test_dir()
  m <- make_morph_matrix()
  # metaMDS may warn about near-zero stress on well-separated synthetic taxa
  expect_output(
    res <- suppressWarnings(suppressMessages(
      morph_nmds(m, nmds_trymax = 20, file_formats = c("pdf", "jpeg"))
    )),
    "Stress"
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("block", "n", "distance", "stress", "quality"))
  expect_equal(nrow(res), 7)
  expect_true(all(res$quality %in% c("excellent", "good", "acceptable", "poor")))

  files <- list.files(out_dir("Figs_NMDS"))
  expect_true(all(c("NMDS_plot_veg.pdf", "NMDS_plot_veg.jpeg") %in% files))
})

test_that("morph_nmds supports vegan distances and handles incomplete data", {
  local_test_dir()
  m <- make_morph_matrix(colors = FALSE)
  res <- suppressWarnings(
    morph_nmds(m, nmds_blocks = "veg", nmds_dist = "euclidean",
               nmds_trymax = 10, file_formats = "pdf", verbose = FALSE)
  )
  expect_equal(res$distance, "euclidean")

  m[3:36, "fruit_stipe_length/FRSTlng"] <- NA
  expect_warning(
    res2 <- morph_nmds(m, nmds_blocks = "fru", nmds_dist = "bray",
                       nmds_trymax = 10, file_formats = "pdf", verbose = FALSE),
    "not enough complete specimens"
  )
  expect_equal(nrow(res2), 0)
})

test_that("morph_nmds skips blocks with too few specimens", {
  local_test_dir()
  m <- make_morph_matrix()
  m[3:36, attr(m, "base_cols")$flo] <- NA
  attr(m, "taxon_colors") <- c(alpha = "red")
  w <- testthat::capture_warnings(
    res <- morph_nmds(m, nmds_blocks = c("veg", "flo"), nmds_trymax = 10,
                      file_formats = "pdf", verbose = FALSE)
  )
  expect_true(any(grepl("Skipping 'flo': not enough specimens", w)))
  expect_true(any(grepl("No stored colour found", w)))
  expect_equal(res$block, "veg")
})

test_that("morph_nmds errors on distances unknown to vegan", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(
    morph_nmds(m, nmds_blocks = "veg", nmds_dist = "notadistance",
               file_formats = "pdf", verbose = FALSE),
    "invalid distance method"
  )
})

test_that("morph_nmds reports failures from metaMDS", {
  local_test_dir()
  m <- make_morph_matrix()
  local_mocked_bindings(
    metaMDS = function(...) stop("did not converge"),
    .package = "vegan"
  )
  expect_warning(
    res <- morph_nmds(m, nmds_blocks = "veg", file_formats = "pdf",
                      verbose = FALSE),
    "NMDS failed for block 'veg': did not converge"
  )
  expect_equal(nrow(res), 0)
})

test_that("morph_nmds validates its arguments", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_nmds(m, nmds_dist = c("a", "b")), "nmds_dist")
  expect_error(morph_nmds(m, nmds_k = 0), "nmds_k")
  expect_error(morph_nmds(m, nmds_trymax = 0), "nmds_trymax")
  expect_error(morph_nmds(m, nmds_seed = "a"), "nmds_seed")
})
