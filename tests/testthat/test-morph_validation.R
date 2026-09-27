# Input checks shared by the morph_* analysis functions

analysis_fns <- list(
  morph_normality  = function(x) morph_normality(x, verbose = FALSE),
  morph_similarity = function(x) morph_similarity(x, verbose = FALSE),
  morph_pca        = function(x) morph_pca(x, verbose = FALSE),
  morph_pcoa       = function(x) morph_pcoa(x, verbose = FALSE),
  morph_nmds       = function(x) morph_nmds(x, verbose = FALSE),
  morph_dendrogram = function(x) morph_dendrogram(x, verbose = FALSE),
  morph_kmeans     = function(x) morph_kmeans(x, verbose = FALSE),
  morph_nmm        = function(x) morph_nmm(x, verbose = FALSE),
  morph_anova      = function(x) morph_anova(x, verbose = FALSE),
  morph_kruskal    = function(x) morph_kruskal(x, verbose = FALSE),
  morph_permanova  = function(x) morph_permanova(x, verbose = FALSE),
  morph_lda        = function(x) morph_lda(x, verbose = FALSE),
  morph_cva        = function(x) morph_cva(x, verbose = FALSE),
  morph_dapc       = function(x) morph_dapc(x, verbose = FALSE)
)

test_that("morph_* functions reject non-data-frame input", {
  local_test_dir()
  for (nm in names(analysis_fns)) {
    expect_error(analysis_fns[[nm]](list(a = 1)), "must be a data frame", info = nm)
  }
  expect_error(morph_boxplots(list(a = 1), verbose = FALSE), "must be a data frame")
  expect_error(morph_autocorrelation(list(a = 1)), "must be a data frame")
})

test_that("morph_* functions require a taxon column", {
  local_test_dir()
  m <- make_morph_matrix()
  no_taxon <- m
  no_taxon$taxon <- NULL
  attr(no_taxon, "base_cols") <- attr(m, "base_cols")
  for (nm in names(analysis_fns)) {
    expect_error(analysis_fns[[nm]](no_taxon), "taxon", info = nm)
  }
  expect_error(morph_boxplots(no_taxon, verbose = FALSE), "taxon")
})

test_that("morph_* functions require the base_cols attribute", {
  local_test_dir()
  m <- make_morph_matrix()
  attr(m, "base_cols") <- NULL
  for (nm in names(analysis_fns)) {
    expect_error(analysis_fns[[nm]](m), "base_cols", info = nm)
  }
  expect_error(morph_autocorrelation(m, order_by = "none"), "base_cols")
})

test_that("morph_* functions reject unknown block names", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_similarity(m, sim_blocks = "zzz", verbose = FALSE), "No valid blocks")
  expect_error(morph_pca(m, pca_blocks = "zzz", verbose = FALSE), "No valid blocks")
  expect_error(morph_pcoa(m, pcoa_blocks = "zzz", verbose = FALSE), "No valid blocks")
  expect_error(morph_nmds(m, nmds_blocks = "zzz", verbose = FALSE), "No valid blocks")
  expect_error(morph_dendrogram(m, dendro_blocks = "zzz", verbose = FALSE), "No valid blocks")
  expect_error(morph_kmeans(m, kmeans_blocks = "zzz", verbose = FALSE), "No valid blocks")
  expect_error(morph_nmm(m, nmm_blocks = "zzz", verbose = FALSE), "No valid blocks")
  expect_error(morph_anova(m, anova_blocks = "zzz", verbose = FALSE), "No valid blocks")
  expect_error(morph_kruskal(m, kw_blocks = "zzz", verbose = FALSE), "No valid blocks")
  expect_error(morph_permanova(m, permanova_blocks = "zzz", verbose = FALSE), "No valid blocks")
  expect_error(morph_lda(m, lda_blocks = "zzz", verbose = FALSE), "No valid blocks")
  expect_error(morph_cva(m, cva_blocks = "zzz", verbose = FALSE), "No valid blocks")
  expect_error(morph_dapc(m, dapc_blocks = "zzz", verbose = FALSE), "No valid blocks")
  expect_error(morph_autocorrelation(m, blocks = "zzz", order_by = "none",
                                     verbose = FALSE), "No valid blocks")
})

test_that("morph_* functions reject non-character block arguments", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_pca(m, pca_blocks = 1), "pca_blocks")
  expect_error(morph_pcoa(m, pcoa_blocks = 1), "pcoa_blocks")
  expect_error(morph_nmds(m, nmds_blocks = 1), "nmds_blocks")
  expect_error(morph_dendrogram(m, dendro_blocks = 1), "dendro_blocks")
  expect_error(morph_kmeans(m, kmeans_blocks = 1), "kmeans_blocks")
  expect_error(morph_nmm(m, nmm_blocks = 1), "nmm_blocks")
  expect_error(morph_anova(m, anova_blocks = 1), "anova_blocks")
  expect_error(morph_kruskal(m, kw_blocks = 1), "kw_blocks")
  expect_error(morph_permanova(m, permanova_blocks = 1), "permanova_blocks")
  expect_error(morph_lda(m, lda_blocks = 1), "lda_blocks")
  expect_error(morph_cva(m, cva_blocks = 1), "cva_blocks")
  expect_error(morph_dapc(m, dapc_blocks = 1), "dapc_blocks")
  expect_error(morph_autocorrelation(m, blocks = 1, order_by = "none"), "blocks")
})

test_that("morph_* functions validate `verbose`", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_similarity(m, verbose = "yes"), "verbose")
  expect_error(morph_pca(m, verbose = NA_character_), "verbose")
  expect_error(morph_pcoa(m, verbose = "yes"), "verbose")
  expect_error(morph_nmds(m, verbose = "yes"), "verbose")
  expect_error(morph_dendrogram(m, verbose = "yes"), "verbose")
  expect_error(morph_kmeans(m, verbose = "yes"), "verbose")
  expect_error(morph_nmm(m, verbose = "yes"), "verbose")
  expect_error(morph_anova(m, verbose = "yes"), "verbose")
  expect_error(morph_kruskal(m, verbose = "yes"), "verbose")
  expect_error(morph_permanova(m, verbose = "yes"), "verbose")
  expect_error(morph_lda(m, verbose = "yes"), "verbose")
  expect_error(morph_cva(m, verbose = "yes"), "verbose")
  expect_error(morph_dapc(m, verbose = "yes"), "verbose")
  expect_error(morph_autocorrelation(m, order_by = "none", verbose = "yes"), "verbose")
})

test_that("morph_* functions validate `file_formats`", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_normality(m, file_formats = "png", verbose = FALSE))
  expect_error(morph_boxplots(m, file_formats = "png", verbose = FALSE))
  expect_error(morph_pca(m, file_formats = "png", verbose = FALSE))
  expect_error(morph_pcoa(m, file_formats = "png", verbose = FALSE))
  expect_error(morph_nmds(m, file_formats = "png", verbose = FALSE))
  expect_error(morph_dendrogram(m, file_formats = "png", verbose = FALSE))
  expect_error(morph_kmeans(m, file_formats = "png", verbose = FALSE))
  expect_error(morph_nmm(m, file_formats = "png", verbose = FALSE))
  expect_error(morph_anova(m, file_formats = "png", verbose = FALSE))
  expect_error(morph_kruskal(m, file_formats = "png", verbose = FALSE))
  expect_error(morph_permanova(m, file_formats = "png", verbose = FALSE))
  expect_error(morph_lda(m, file_formats = "png", verbose = FALSE))
  expect_error(morph_cva(m, file_formats = "png", verbose = FALSE))
  expect_error(morph_dapc(m, file_formats = "png", verbose = FALSE))
  expect_error(morph_autocorrelation(m, order_by = "none", file_formats = "png"))
  expect_error(morph_autocorrelation(m, order_by = "none", file_formats = character(0)),
               "file_formats")
})
