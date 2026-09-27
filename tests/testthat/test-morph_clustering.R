# Tests for clustering functions: morph_dendrogram(), morph_kmeans(),
# morph_nmm()

#_______________________________________________________________________________
# morph_dendrogram() ####

test_that("morph_dendrogram draws diagnostics and dendrograms for each block", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_message(
    expect_output(
      out <- morph_dendrogram(m, k_max = 4, file_formats = c("pdf", "jpeg"))
    ),
    "Best linkage method"
  )
  expect_null(out)

  files <- list.files(out_dir("Figs.dendrogram"))
  for (b in c("veg", "flo", "fru", "vegflofru")) {
    expect_true(paste0("Dendrogram_", b, "_vert.pdf") %in% files, info = b)
    expect_true(paste0("Dendrogram_", b, "_hor.jpeg") %in% files, info = b)
    expect_true(paste0("linkage_selection_", b, ".pdf") %in% files, info = b)
    expect_true(paste0("nbclust_silhouette_", b, ".jpeg") %in% files, info = b)
  }
})

test_that("morph_dendrogram handles custom colours, linkages and two-trait blocks", {
  local_test_dir()
  m <- make_morph_matrix(colors = FALSE)
  attr(m, "base_cols") <- list(veg = attr(m, "base_cols")$veg[1:2])

  expect_warning(
    msgs <- testthat::capture_messages(
      morph_dendrogram(m, dist_method = "euclidean",
                       linkage_candidates = c("average", "not_a_method"),
                       k_max = 6, cluster_cols = "black",
                       file_formats = "pdf", verbose = TRUE) |>
        utils::capture.output()
    ),
    "`cluster_cols` has fewer colors"
  )
  expect_true(any(grepl("Linkage method\\(s\\) skipped.*not_a_method", msgs)))
  expect_true(file.exists(file.path(out_dir("Figs.dendrogram"),
                                    "Dendrogram_veg_vert.pdf")))
})

test_that("morph_dendrogram skips blocks it cannot cluster", {
  local_test_dir()
  m <- make_morph_matrix()
  attr(m, "base_cols") <- list(veg = attr(m, "base_cols")$veg,
                               flo = attr(m, "base_cols")$flo[1])
  attr(m, "taxon_colors") <- c(alpha = "red")

  w <- testthat::capture_warnings(
    morph_dendrogram(m, dendro_blocks = c("veg", "flo"),
                     linkage_candidates = "complete", cluster_cols = rep("red", 8),
                     file_formats = "pdf", verbose = FALSE)
  )
  expect_true(any(grepl("Skipping 'flo': not enough data", w)))
  expect_true(any(grepl("No stored colour found", w)))

  expect_warning(
    morph_dendrogram(m, dendro_blocks = "veg", linkage_candidates = "nope",
                     file_formats = "pdf", verbose = FALSE),
    "cophenetic correlation could not be computed"
  )
})

test_that("morph_dendrogram validates its arguments", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_dendrogram(m, dist_method = c("a", "b")), "dist_method")
  expect_error(morph_dendrogram(m, linkage_candidates = 1), "linkage_candidates")
  expect_error(morph_dendrogram(m, k_max = 1), "k_max")
})

#_______________________________________________________________________________
# morph_kmeans() ####

test_that("morph_kmeans estimates k with the gap statistic", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_message(
    res <- morph_kmeans(m, kmeans_blocks = c("veg", "flo"), kmeans_k_max = 4,
                        kmeans_B = 5, kmeans_nstart = 5,
                        file_formats = c("pdf", "jpeg")),
    "K-means clustering done"
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("block", "id", "taxon", "cluster", "optimal_k"))
  expect_setequal(unique(res$block), c("veg", "flo"))
  expect_true(all(res$optimal_k >= 2 & res$optimal_k <= 4))

  files <- list.files(out_dir("Figs.KMeans"))
  expect_true(all(c("KMeans_elbow_veg.pdf", "KMeans_plot_veg.jpeg",
                    "KMeans_fviz_flo.pdf") %in% files))
})

test_that("morph_kmeans accepts a fixed k and runs every block", {
  local_test_dir()
  m <- make_morph_matrix(colors = FALSE)
  res <- morph_kmeans(m, kmeans_k = 3, kmeans_nstart = 5,
                      file_formats = "pdf", verbose = TRUE) |>
    suppressMessages()
  expect_setequal(unique(res$block),
                  c("veg", "flo", "fru", "vegflo", "vegfru", "flofru", "vegflofru"))
  expect_true(all(res$optimal_k == 3))
  expect_false(any(grepl("elbow", list.files(out_dir("Figs.KMeans")))))
})

test_that("morph_kmeans skips blocks it cannot cluster", {
  local_test_dir()
  m <- make_morph_matrix()
  attr(m, "taxon_colors") <- c(alpha = "red")

  expect_warning(
    expect_warning(
      res <- morph_kmeans(m, kmeans_blocks = "veg", kmeans_k = 100,
                          file_formats = "pdf", verbose = FALSE),
      "must be smaller than the number of specimens"
    ),
    "No trait block could be analyzed"
  )
  expect_null(res)

  expect_warning(
    expect_warning(
      morph_kmeans(m, kmeans_blocks = "veg", kmeans_pcs = c("PC1", "PC9"),
                   kmeans_k = 2, file_formats = "pdf", verbose = FALSE),
      "component\\(s\\) not available: PC9"
    ),
    "No trait block"
  )

  m1 <- make_morph_matrix()
  attr(m1, "base_cols") <- list(veg = attr(m1, "base_cols")$veg[1])
  expect_warning(
    expect_warning(
      morph_kmeans(m1, kmeans_k = 2, file_formats = "pdf", verbose = FALSE),
      "not enough data after removing missing values"
    ),
    "No trait block"
  )

  w <- testthat::capture_warnings(
    morph_kmeans(m, kmeans_blocks = "fru", kmeans_k = 2,
                 file_formats = "pdf", verbose = FALSE)
  )
  expect_true(any(grepl("No stored colour found", w)))
})

test_that("morph_kmeans validates its arguments", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_kmeans(m, kmeans_k = 1), "kmeans_k")
  expect_error(morph_kmeans(m, kmeans_k_max = 1), "kmeans_k_max")
  expect_error(morph_kmeans(m, kmeans_B = 0), "kmeans_B")
  expect_error(morph_kmeans(m, kmeans_nstart = 0), "kmeans_nstart")
  expect_error(morph_kmeans(m, kmeans_pcs = "PC1"), "kmeans_pcs")
  expect_error(morph_kmeans(m, kmeans_seed = "a"), "kmeans_seed")
})

#_______________________________________________________________________________
# morph_nmm() ####

test_that("morph_nmm fits normal mixture models and saves group tables", {
  local_test_dir()
  m <- make_morph_matrix(n_per = 15, with_na = FALSE)
  expect_message(
    res <- morph_nmm(m, nmm_blocks = c("veg", "flo"), nmm_G = 1:3,
                     file_formats = c("pdf", "jpeg")),
    "NMM analysis done"
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("block", "n", "model", "k", "bic", "top_traits"))
  expect_equal(res$block, c("veg", "flo"))

  od <- out_dir("Figs.NMM")
  expect_true(file.exists(file.path(od, "NMM_summary.xlsx")))
  groups <- openxlsx::read.xlsx(file.path(od, "NMM_groups_veg.xlsx"))
  expect_true(all(c("id", "taxon", "nmm_group", "uncertainty") %in% names(groups)))
  files <- list.files(od)
  expect_true(all(c("NMM_BIC_veg.pdf", "NMM_classification_veg.jpeg",
                    "NMM_uncertainty_flo.pdf") %in% files))
})

test_that("morph_nmm runs a single-group model without ellipses or Excel", {
  local_test_dir()
  m <- make_morph_matrix(n_per = 15, with_na = FALSE)
  msgs <- testthat::capture_messages(
    res <- morph_nmm(m, nmm_blocks = "veg", nmm_G = 1, nmm_top_traits = 2,
                     save_xlsx = FALSE, file_formats = "pdf", verbose = TRUE)
  )
  expect_true(any(grepl("skipping confidence ellipses", msgs)))
  expect_equal(res$k, 1)
  expect_equal(lengths(strsplit(res$top_traits, "|", fixed = TRUE)), 2)
  expect_false(file.exists(file.path(out_dir("Figs.NMM"), "NMM_summary.xlsx")))
})

test_that("morph_nmm skips blocks it cannot model", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_warning(
    expect_warning(
      res <- morph_nmm(m, nmm_blocks = "veg", nmm_min_n = 500,
                       file_formats = "pdf", verbose = FALSE),
      "Skipping 'veg': only 36 of 36 specimens"
    ),
    "No trait block could be analyzed"
  )
  expect_null(res)

  m[, attr(m, "base_cols")$flo[-1]] <- 1
  expect_warning(
    expect_warning(
      morph_nmm(m, nmm_blocks = "flo", file_formats = "pdf", verbose = FALSE),
      "fewer than 2 traits vary"
    ),
    "No trait block"
  )

  expect_warning(
    expect_warning(
      morph_nmm(make_morph_matrix(), nmm_blocks = "veg",
                nmm_model_names = "NOTAMODEL",
                file_formats = "pdf", verbose = FALSE),
      "NMM failed for block 'veg'"
    ),
    "No trait block"
  )
})

test_that("morph_nmm validates its arguments", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_nmm(m, nmm_top_traits = 1), "nmm_top_traits")
  expect_error(morph_nmm(m, nmm_min_n = 2), "nmm_min_n")
  expect_error(morph_nmm(m, nmm_G = 0), "nmm_G")
  expect_error(morph_nmm(m, nmm_model_names = 1), "nmm_model_names")
  expect_error(morph_nmm(m, ellipse_level = 1), "ellipse_level")
  expect_error(morph_nmm(m, nmm_seed = "a"), "nmm_seed")
  expect_error(morph_nmm(m, save_xlsx = "yes"), "save_xlsx")
})
