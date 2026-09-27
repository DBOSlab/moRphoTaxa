# Tests for confirmatory and discriminant functions:
# morph_anova(), morph_kruskal(), morph_permanova(), morph_lda(),
# morph_cva(), morph_dapc()

# Matrix with an extra grouping column (e.g. a region) besides `taxon`
with_region <- function(m = make_morph_matrix()) {
  m$region <- rep(c("north", "south"), length.out = nrow(m))
  m
}

#_______________________________________________________________________________
# morph_anova() ####

test_that("morph_anova tests every trait by group and block", {
  local_test_dir()
  m <- with_region()
  expect_message(
    res <- morph_anova(m, anova_groups = c("taxon", "region"),
                       file_formats = c("pdf", "jpeg")),
    "ANOVA analysis done"
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("trait", "F_statistic", "p_value", "p_adjusted",
                    "significant", "block", "group") %in% names(res)))
  expect_setequal(unique(res$group), c("taxon", "region"))
  expect_equal(length(unique(res$block)), 7)
  # Taxa are well separated in the synthetic data
  expect_true(all(res$significant[res$group == "taxon" & res$block == "veg"]))

  od <- out_dir("Figs.ANOVA")
  expect_true(file.exists(file.path(od, "ANOVA_all_results.xlsx")))
  expect_true(file.exists(file.path(od, "ANOVA_veg_taxon.pdf")))
  expect_true(file.exists(file.path(od, "ANOVA_veg_region.jpeg")))
  expect_true(file.exists(file.path(od, "ANOVA_veg_taxon.xlsx")))
})

test_that("morph_anova skips unusable blocks and can skip Excel output", {
  local_test_dir()
  m <- make_morph_matrix()
  m$single <- "one"
  expect_warning(
    expect_warning(
      res <- morph_anova(m, anova_blocks = "veg", anova_groups = "single",
                         file_formats = "pdf", verbose = FALSE),
      "fewer than 2 group levels"
    ),
    "No block x group combination could be analyzed"
  )
  expect_null(res)

  m2 <- make_morph_matrix()
  m2[, attr(m2, "base_cols")$flo] <- NA
  m2[1:2, attr(m2, "base_cols")$flo] <- 1
  w <- testthat::capture_warnings(
    res2 <- morph_anova(m2, anova_blocks = c("veg", "flo"), save_xlsx = FALSE,
                        anova_p_adjust = "bonferroni", file_formats = "pdf",
                        verbose = FALSE)
  )
  expect_true(any(grepl("no trait had enough data", w)))
  expect_equal(unique(res2$block), "veg")
  expect_false(file.exists(file.path(out_dir("Figs.ANOVA"), "ANOVA_all_results.xlsx")))
})

test_that("morph_anova validates its arguments", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_anova(m, anova_groups = 1), "anova_groups")
  expect_error(morph_anova(m, anova_groups = "nope"), "not found in `analysis_data`: nope")
  expect_error(morph_anova(m, anova_p_adjust = "nope"))
  expect_error(morph_anova(m, anova_alpha = 1), "anova_alpha")
  expect_error(morph_anova(m, save_xlsx = "yes"), "save_xlsx")
})

#_______________________________________________________________________________
# morph_kruskal() ####

test_that("morph_kruskal runs Kruskal-Wallis with Dunn post-hoc tests", {
  local_test_dir()
  m <- with_region()
  expect_message(
    res <- morph_kruskal(m, kw_groups = c("taxon", "region"),
                         file_formats = c("pdf", "jpeg")),
    "Kruskal-Wallis analysis done"
  )
  expect_named(res, c("kw", "posthoc"))
  expect_true(all(c("trait", "H_stat", "df", "p_value", "p_adjusted",
                    "significant", "block", "group") %in% names(res$kw)))
  expect_true(all(c("comparison", "Z", "p_adjusted", "trait") %in%
                    names(res$posthoc)))
  # Three taxa -> three pairwise comparisons per significant trait
  veg_ph <- res$posthoc[res$posthoc$block == "veg" & res$posthoc$group == "taxon", ]
  expect_equal(nrow(veg_ph), 3 * 4)

  od <- out_dir("Figs.KruskalWallis")
  expect_true(file.exists(file.path(od, "KW_all_results.xlsx")))
  expect_true(file.exists(file.path(od, "KW_posthoc_all_results.xlsx")))
  ph_dir <- file.path(od, "posthoc_veg_taxon")
  expect_true(file.exists(file.path(ph_dir, "Dunn_all_results.xlsx")))
  expect_true(any(grepl("^Dunn_.*\\.jpeg$", list.files(ph_dir))))
})

test_that("morph_kruskal can skip post-hoc tests and Excel output", {
  local_test_dir()
  m <- make_morph_matrix()
  res <- morph_kruskal(m, kw_blocks = "veg", kw_posthoc = FALSE,
                       kw_p_adjust = "hommel", save_xlsx = FALSE,
                       file_formats = "pdf", verbose = FALSE)
  expect_null(res$posthoc)
  expect_false(dir.exists(file.path(out_dir("Figs.KruskalWallis"),
                                    "posthoc_veg_taxon")))
  expect_false(file.exists(file.path(out_dir("Figs.KruskalWallis"),
                                     "KW_all_results.xlsx")))
})

test_that("morph_kruskal reports Dunn failures and unusable blocks", {
  local_test_dir()
  m <- make_morph_matrix()
  local_mocked_bindings(
    dunn.test = function(...) stop("dunn exploded"),
    .package = "dunn.test"
  )
  w <- testthat::capture_warnings(
    res <- morph_kruskal(m, kw_blocks = "veg", file_formats = "pdf",
                         verbose = FALSE)
  )
  expect_true(any(grepl("Dunn post-hoc failed.*dunn exploded", w)))
  expect_null(res$posthoc)

  m$single <- "one"
  expect_warning(
    expect_warning(
      res2 <- morph_kruskal(m, kw_blocks = "veg", kw_groups = "single",
                            file_formats = "pdf", verbose = FALSE),
      "fewer than 2 group levels"
    ),
    "No block x group combination"
  )
  expect_null(res2)

  m3 <- make_morph_matrix()
  m3[, attr(m3, "base_cols")$flo] <- NA
  w3 <- testthat::capture_warnings(
    morph_kruskal(m3, kw_blocks = c("veg", "flo"), kw_posthoc = FALSE,
                  file_formats = "pdf", verbose = FALSE)
  )
  expect_true(any(grepl("no trait had enough data", w3)))
})

test_that("morph_kruskal validates its arguments", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_kruskal(m, kw_groups = 1), "kw_groups")
  expect_error(morph_kruskal(m, kw_groups = "nope"), "not found")
  expect_error(morph_kruskal(m, kw_p_adjust = "nope"))
  expect_error(morph_kruskal(m, kw_posthoc = "yes"), "kw_posthoc")
  expect_error(morph_kruskal(m, kw_p_adjust = "hommel"), "no equivalent method")
  expect_error(morph_kruskal(m, kw_alpha = 0), "kw_alpha")
  expect_error(morph_kruskal(m, save_xlsx = "yes"), "save_xlsx")
})

#_______________________________________________________________________________
# morph_permanova() ####

test_that("morph_permanova runs PERMANOVA and PERMDISP for each block", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_message(
    res <- morph_permanova(m, permutations = 99, file_formats = c("pdf", "jpeg")),
    "PERMANOVA \\+ PERMDISP analysis done"
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("block", "n", "n_groups", "permanova_F", "permanova_R2",
                      "permanova_p", "permdisp_F", "permdisp_p"))
  expect_equal(nrow(res), 7)
  expect_true(all(res$n_groups == 3))
  expect_true(all(res$permanova_p <= 0.05))

  od <- out_dir("Figs_PERMANOVA")
  expect_true(file.exists(file.path(od, "PERMANOVA_summary.xlsx")))
  sheets <- openxlsx::getSheetNames(file.path(od, "PERMANOVA_veg.xlsx"))
  expect_equal(sheets, c("permanova", "permdisp", "pairwise_permanova",
                         "pairwise_permdisp", "distances_to_centroid"))
  expect_true(file.exists(file.path(od, "PERMDISP_boxplot_veg.jpeg")))
  expect_true(file.exists(file.path(od, "PERMANOVA_ordination_veg.pdf")))
})

test_that("morph_permanova supports other distances, two groups and small groups", {
  local_test_dir()
  m <- make_morph_matrix(colors = FALSE)
  eu <- morph_permanova(m, permanova_blocks = "veg", dist_method = "euclidean",
                        permutations = 99, seed = NULL, save_xlsx = FALSE,
                        file_formats = "pdf", verbose = FALSE)
  expect_equal(eu$block, "veg")
  expect_false(file.exists(file.path(out_dir("Figs_PERMANOVA"),
                                     "PERMANOVA_summary.xlsx")))

  br <- morph_permanova(m, permanova_blocks = "flo", dist_method = "bray",
                        permutations = 99, file_formats = "pdf", verbose = FALSE)
  expect_equal(br$block, "flo")

  # Drop most "gamma" specimens: they fall below min_n and are removed
  m2 <- subset_rows(make_morph_matrix(), 1:26)
  attr(m2, "taxon_colors") <- c(alpha = "red")
  msgs <- NULL
  expect_warning(
    msgs <- testthat::capture_messages(
      res <- morph_permanova(m2, permanova_blocks = "veg", min_n = 5,
                             permutations = 99, file_formats = "pdf")
    ),
    "No stored colour found for: beta"
  )
  expect_true(any(grepl("Removing small groups: gamma", msgs)))
  expect_equal(res$n_groups, 2)
})

test_that("morph_permanova skips blocks it cannot analyse", {
  local_test_dir()
  m <- make_morph_matrix()
  neg <- m
  neg[1, "petiole_length/PETIlng"] <- -5
  expect_warning(
    res <- morph_permanova(neg, permanova_blocks = "veg", dist_method = "bray",
                           permutations = 99, file_formats = "pdf",
                           verbose = FALSE),
    "Bray-Curtis distance requires non-negative"
  )
  expect_null(res)

  sparse <- make_morph_matrix()
  sparse[2:36, "fruit_stipe_length/FRSTlng"] <- NA
  expect_warning(
    morph_permanova(sparse, permanova_blocks = "fru", permutations = 99,
                    file_formats = "pdf", verbose = FALSE),
    "fewer than 2 specimens with complete data"
  )

  expect_warning(
    morph_permanova(make_morph_matrix(), permanova_blocks = "veg", min_n = 20,
                    permutations = 99, file_formats = "pdf", verbose = FALSE),
    "not enough data"
  )
})

test_that("morph_permanova reports failures from vegan", {
  local_test_dir()
  m <- make_morph_matrix()
  local_mocked_bindings(adonis2 = function(...) stop("boom"), .package = "vegan")
  expect_warning(
    res <- morph_permanova(m, permanova_blocks = "veg", permutations = 99,
                           file_formats = "pdf", verbose = FALSE),
    "PERMANOVA failed for 'veg': boom"
  )
  expect_null(res)
})

test_that("morph_permanova reports PERMDISP failures", {
  local_test_dir()
  m <- make_morph_matrix()
  local_mocked_bindings(betadisper = function(...) stop("no disp"), .package = "vegan")
  expect_warning(
    morph_permanova(m, permanova_blocks = "veg", permutations = 99,
                    file_formats = "pdf", verbose = FALSE),
    "PERMDISP failed for 'veg': no disp"
  )
})

test_that("morph_permanova reports PERMDISP permutation-test failures", {
  local_test_dir()
  m <- make_morph_matrix()
  local_mocked_bindings(permutest = function(...) stop("no perm"), .package = "vegan")
  expect_warning(
    morph_permanova(m, permanova_blocks = "veg", permutations = 99,
                    file_formats = "pdf", verbose = FALSE),
    "PERMDISP permutation test failed for 'veg': no perm"
  )
})

test_that("morph_permanova validates its arguments", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_permanova(m, dist_method = "jaccard"))
  expect_error(morph_permanova(m, permutations = 10), "permutations")
  expect_error(morph_permanova(m, min_n = 1), "min_n")
  expect_error(morph_permanova(m, p_adjust = "nope"), "p_adjust")
  expect_error(morph_permanova(m, ellipse_level = 0), "ellipse_level")
  expect_error(morph_permanova(m, seed = "a"), "seed")
  expect_error(morph_permanova(m, save_xlsx = "a"), "save_xlsx")
})

#_______________________________________________________________________________
# morph_lda() ####

test_that("morph_lda fits discriminant functions and saves scores", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_message(
    res <- morph_lda(m, file_formats = c("pdf", "jpeg")),
    "LDA analysis done"
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("block", "group", "n", "n_groups", "n_lds", "ld1_var"))
  expect_equal(nrow(res), 7)
  expect_true(all(res$n_lds == 2))

  od <- out_dir("Figs.LDA")
  expect_true(file.exists(file.path(od, "LDA_summary.xlsx")))
  scores <- openxlsx::read.xlsx(file.path(od, "LDA_scores_veg_taxon.xlsx"))
  expect_true(all(c("id", "taxon", "LD1", "LD2") %in% names(scores)))
  expect_true(file.exists(file.path(od, "LDA_plot_veg_taxon.jpeg")))
  expect_true(file.exists(file.path(od, "LDA_loadings_veg_taxon.pdf")))
})

test_that("morph_lda handles two groups, other grouping columns and group-mean imputation", {
  local_test_dir()
  m <- with_region()
  msgs <- testthat::capture_messages(
    res <- morph_lda(m, lda_blocks = "veg", lda_group = "region",
                     lda_impute = "group_mean", save_xlsx = FALSE,
                     file_formats = "pdf")
  )
  # Only one discriminant axis exists for two groups
  expect_true(any(grepl("Requested LD\\(s\\) not available: LD2", msgs)))
  expect_equal(res$n_lds, 1)
  expect_true(file.exists(file.path(out_dir("Figs.LDA"), "LDA_plot_veg_region.pdf")))
  expect_false(file.exists(file.path(out_dir("Figs.LDA"), "LDA_summary.xlsx")))

  res1 <- morph_lda(make_morph_matrix(), lda_blocks = "flo", lda_lds = "LD1",
                    file_formats = "pdf", verbose = FALSE)
  expect_equal(res1$block, "flo")
})

test_that("morph_lda falls back to the overall mean for groups without data", {
  local_test_dir()
  m <- make_morph_matrix()
  m[m$taxon == "gamma", "leaf_length/LEAFlng"] <- NA
  expect_warning(
    morph_lda(m, lda_blocks = "veg", lda_impute = "group_mean",
              file_formats = "pdf", verbose = FALSE),
    "Group 'gamma' has no non-missing values"
  )
})

test_that("morph_lda removes small groups and skips unusable blocks", {
  local_test_dir()
  m <- subset_rows(make_morph_matrix(), 1:26)
  attr(m, "taxon_colors") <- c(alpha = "red")
  w <- NULL
  msgs <- testthat::capture_messages(
    w <- testthat::capture_warnings(
      res <- morph_lda(m, lda_blocks = "veg", lda_min_n = 5, file_formats = "pdf")
    )
  )
  expect_true(any(grepl("Removing groups with fewer than 5 specimens: gamma", msgs)))
  expect_true(any(grepl("No stored colour found for: beta", w)))
  expect_equal(res$n_groups, 2)

  expect_warning(
    expect_warning(
      morph_lda(make_morph_matrix(), lda_blocks = "veg", lda_min_n = 50,
                file_formats = "pdf", verbose = FALSE),
      "fewer than 2 groups remaining"
    ),
    "No trait block could be analyzed"
  )

  const <- make_morph_matrix()
  const[, attr(const, "base_cols")$flo] <- 1
  expect_warning(
    expect_warning(
      morph_lda(const, lda_blocks = "flo", file_formats = "pdf", verbose = FALSE),
      "no variable traits after filtering"
    ),
    "No trait block"
  )
})

test_that("morph_lda reports MASS::lda failures", {
  local_test_dir()
  m <- make_morph_matrix()
  local_mocked_bindings(lda = function(...) stop("singular"), .package = "MASS")
  expect_warning(
    expect_warning(
      morph_lda(m, lda_blocks = "veg", file_formats = "pdf", verbose = FALSE),
      "LDA failed for 'veg': singular"
    ),
    "No trait block"
  )
})

test_that("morph_lda validates its arguments", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_lda(m, lda_group = c("a", "b")), "lda_group")
  expect_error(morph_lda(m, lda_group = "nope"), "not found")
  expect_error(morph_lda(m, lda_lds = c("LD1", "LD2", "LD3")), "lda_lds")
  expect_error(morph_lda(m, lda_ellipse = 2), "lda_ellipse")
  expect_error(morph_lda(m, lda_min_n = 1), "lda_min_n")
  expect_error(morph_lda(m, lda_impute = "median"))
  expect_error(morph_lda(m, save_xlsx = 1), "save_xlsx")
})

#_______________________________________________________________________________
# morph_cva() ####

test_that("morph_cva computes canonical variates for each block", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_message(
    expect_output(
      res <- morph_cva(m, file_formats = c("pdf", "jpeg")),
      "Canonical axes"
    ),
    "CVA analysis done"
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("block", "n", "n_groups", "can1_pct", "can1_r"))
  expect_equal(nrow(res), 7)
  # Percentages of explained variance (not x100)
  expect_true(all(res$can1_pct > 0 & res$can1_pct <= 100))
  expect_true(all(res$can1_r > 0 & res$can1_r <= 1))

  od <- out_dir("Figs_CVA")
  expect_true(file.exists(file.path(od, "CVA_summary.xlsx")))
  expect_equal(openxlsx::getSheetNames(file.path(od, "CVA_veg.xlsx")),
               c("summary", "scores", "std_coefficients", "structure_coefs"))
  expect_true(all(c("CVA_scatter_veg.jpeg", "CVA_coefficients_veg.pdf",
                    "CVA_structure_veg.pdf") %in% list.files(od)))
})

test_that("morph_cva handles two groups, small groups and missing colours", {
  local_test_dir()
  m <- subset_rows(make_morph_matrix(), 1:26)
  attr(m, "taxon_colors") <- c(alpha = "red")
  w <- NULL
  msgs <- testthat::capture_messages(
    w <- testthat::capture_warnings(
      res <- morph_cva(m, cva_blocks = "veg", cva_min_n = 5, file_formats = "pdf",
                       verbose = TRUE) |>
        utils::capture.output()
    )
  )
  expect_true(any(grepl("Removing small groups: gamma", msgs)))
  expect_true(any(grepl("No stored colour found for: beta", w)))

  res2 <- morph_cva(make_morph_matrix(colors = FALSE), cva_blocks = "flo",
                    file_formats = "pdf", verbose = FALSE)
  expect_equal(res2$n_groups, 3)
})

test_that("morph_cva skips blocks it cannot analyse", {
  local_test_dir()
  expect_warning(
    res <- morph_cva(make_morph_matrix(), cva_blocks = "veg", cva_min_n = 50,
                     file_formats = "pdf", verbose = FALSE),
    "fewer than 2 groups"
  )
  expect_null(res)

  m <- make_morph_matrix()
  attr(m, "base_cols") <- list(veg = attr(m, "base_cols")$veg[1])
  expect_warning(
    morph_cva(m, file_formats = "pdf", verbose = FALSE),
    "Skipping 'veg': not enough data"
  )

  local_mocked_bindings(.candisc_fit = function(...) stop("no cva"))
  expect_warning(
    morph_cva(make_morph_matrix(), cva_blocks = "veg", file_formats = "pdf",
              verbose = FALSE),
    "CVA failed for 'veg': no cva"
  )
})

test_that("morph_cva reports MANOVA failures", {
  local_test_dir()
  local_mocked_bindings(manova = function(...) stop("no manova"), .package = "stats")
  expect_warning(
    morph_cva(make_morph_matrix(), cva_blocks = "veg", file_formats = "pdf",
              verbose = FALSE),
    "MANOVA failed for 'veg': no manova"
  )
})

test_that(".candisc_fit returns canonical variates of a one-way MANOVA", {
  m <- make_morph_matrix(with_na = FALSE)
  X <- as.data.frame(scale(m[, -1]))
  X$group <- factor(m$taxon)
  fit <- stats::manova(as.matrix(X[, 1:10]) ~ group, data = X)
  cv <- .candisc_fit(fit)

  # 3 groups -> 2 canonical axes
  expect_equal(cv$rank, 2)
  expect_equal(colnames(cv$coeffs.std), c("Can1", "Can2"))
  expect_equal(sum(cv$pct), 100)
  expect_true(all(cv$canrsq > 0 & cv$canrsq < 1))
  expect_equal(dim(cv$structure), c(10, 2))
  expect_true(all(abs(cv$structure) <= 1))
  expect_equal(names(cv$scores), c("group", "Can1", "Can2"))
  expect_equal(rownames(cv$scores), rownames(m))
  # Canonical scores are uncorrelated and centred
  expect_equal(unname(colMeans(cv$scores[, -1])), c(0, 0), tolerance = 1e-10)
  expect_lt(abs(stats::cor(cv$scores$Can1, cv$scores$Can2)), 1e-8)
  # Can1 separates the taxa: squared canonical correlation = eta^2 of Can1
  a <- summary(stats::aov(cv$scores$Can1 ~ X$group))[[1]]
  expect_equal(cv$canrsq[1], a[["Sum Sq"]][1] / sum(a[["Sum Sq"]]), tolerance = 1e-8)

  # Two groups -> a single axis
  two <- X[X$group != "gamma", ]
  two$group <- droplevels(two$group)
  fit2 <- stats::manova(as.matrix(two[, 1:10]) ~ group, data = two)
  expect_equal(.candisc_fit(fit2)$rank, 1)
})

test_that("morph_cva validates its arguments", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_cva(m, cva_ellipse = 1), "cva_ellipse")
  expect_error(morph_cva(m, cva_min_n = 0), "cva_min_n")
})

#_______________________________________________________________________________
# morph_dapc() ####

test_that("morph_dapc runs DAPC and membership plots for each block", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_message(
    expect_output(
      res <- morph_dapc(m, file_formats = c("pdf", "jpeg")),
      "PCs used"
    ),
    "DAPC analysis done"
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("block", "n", "n_groups", "n_pca", "n_da"))
  expect_equal(nrow(res), 7)
  expect_true(all(res$n_da == 2))

  od <- out_dir("Figs_DAPC")
  expect_true(file.exists(file.path(od, "DAPC_summary.xlsx")))
  expect_equal(openxlsx::getSheetNames(file.path(od, "DAPC_veg.xlsx")),
               c("membership", "da_scores"))
  expect_true(all(c("DAPC_scatter_veg.jpeg", "DAPC_membership_veg.pdf") %in%
                    list.files(od)))
})

test_that("morph_dapc honours user-set PCs/DAs and handles small groups", {
  local_test_dir()
  m <- make_morph_matrix(colors = FALSE)
  res <- morph_dapc(m, dapc_blocks = "veg", dapc_n_pca = 3, dapc_n_da = 1,
                    file_formats = "pdf", verbose = FALSE)
  expect_equal(res$n_pca, 3)
  expect_equal(res$n_da, 1)

  m2 <- subset_rows(make_morph_matrix(), 1:26)
  attr(m2, "taxon_colors") <- c(alpha = "red")
  w <- NULL
  msgs <- testthat::capture_messages(
    w <- testthat::capture_warnings(
      res2 <- morph_dapc(m2, dapc_blocks = "flo", dapc_min_n = 5,
                         file_formats = "pdf") |>
        utils::capture.output()
    )
  )
  expect_true(any(grepl("Removing small groups: gamma", msgs)))
  expect_true(any(grepl("No stored colour found for: beta", w)))
})

test_that("morph_dapc skips blocks it cannot analyse", {
  local_test_dir()
  expect_warning(
    res <- morph_dapc(make_morph_matrix(), dapc_blocks = "veg", dapc_min_n = 50,
                      file_formats = "pdf", verbose = FALSE),
    "fewer than 2 groups"
  )
  expect_null(res)

  m <- make_morph_matrix()
  attr(m, "base_cols") <- list(veg = attr(m, "base_cols")$veg[1])
  expect_warning(
    morph_dapc(m, file_formats = "pdf", verbose = FALSE),
    "Skipping 'veg': not enough data"
  )

  local_mocked_bindings(dapc = function(...) stop("no dapc"), .package = "adegenet")
  expect_warning(
    morph_dapc(make_morph_matrix(), dapc_blocks = "veg", file_formats = "pdf",
               verbose = FALSE),
    "DAPC failed for 'veg': no dapc"
  )
})

test_that("morph_dapc validates its arguments", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_dapc(m, dapc_n_pca = 0), "dapc_n_pca")
  expect_error(morph_dapc(m, dapc_n_da = 0), "dapc_n_da")
  expect_error(morph_dapc(m, dapc_min_n = 0), "dapc_min_n")
  expect_error(morph_dapc(m, dapc_ellipse = 1), "dapc_ellipse")
})
