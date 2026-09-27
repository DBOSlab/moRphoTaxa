# Tests for data-quality and exploratory functions:
# morph_normality(), morph_boxplots(), morph_similarity()

#_______________________________________________________________________________
# morph_normality() ####

test_that("morph_normality writes summary, Q-Q plots and Excel tables", {
  local_test_dir()
  m <- make_morph_matrix()

  expect_message(
    expect_output(
      morph_normality(m, file_formats = c("pdf", "jpeg"), verbose = TRUE),
      "Overall Recommendation"
    ),
    "Running Normality Tests"
  )

  od <- out_dir("Figs_Normality")
  expect_true(file.exists(file.path(od, "Normality_all traits.pdf")))
  expect_true(file.exists(file.path(od, "Normality_all traits.jpeg")))
  expect_true(file.exists(file.path(od, "Normality_all_results.xlsx")))

  res <- openxlsx::read.xlsx(file.path(od, "Normality_all_results.xlsx"))
  expect_equal(nrow(res), 10)
  expect_true(all(res$recommended %in% c("ANOVA", "Kruskal-Wallis")))

  qq <- list.files(file.path(od, "QQ_plots"))
  expect_length(grep("\\.pdf$", qq), 10)
  expect_length(grep("\\.jpeg$", qq), 10)
})

test_that("morph_normality can skip Q-Q plots and too-short traits", {
  local_test_dir()
  m <- make_morph_matrix()
  m[, "seed_length/SEEDlng"] <- NA
  m[1:2, "seed_length/SEEDlng"] <- c(1, 2)

  quiet(morph_normality(m, norm_qq = FALSE, file_formats = "pdf", verbose = FALSE))

  od <- out_dir("Figs_Normality")
  expect_false(dir.exists(file.path(od, "QQ_plots")))
  res <- openxlsx::read.xlsx(file.path(od, "Normality_all_results.xlsx"))
  expect_equal(nrow(res), 9)
  expect_false("seed_length/SEEDlng" %in% res$trait)
})

test_that("morph_normality returns quietly when no trait can be tested", {
  local_test_dir()
  m <- make_morph_matrix(n_per = 1, taxa = c("alpha", "beta"), with_na = FALSE)
  quiet(morph_normality(m, file_formats = "pdf", verbose = FALSE))
  expect_false(file.exists(file.path(out_dir("Figs_Normality"),
                                     "Normality_all_results.xlsx")))
})

test_that("morph_normality validates its arguments", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_normality(m, norm_qq = "yes", verbose = FALSE), "norm_qq")
})

#_______________________________________________________________________________
# morph_boxplots() ####

test_that("morph_boxplots draws violin plots (default) for each trait", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_message(
    morph_boxplots(m, file_formats = c("pdf", "jpeg"), verbose = TRUE),
    "Boxplot analysis done"
  )
  files <- list.files(out_dir("Figs.boxplot"))
  expect_length(grep("^violin_.*\\.pdf$", files), 10)
  expect_length(grep("^violin_.*\\.jpeg$", files), 10)
  expect_length(grep("^boxplot_", files), 0)
})

test_that("morph_boxplots draws boxplots, both types, and without jitter", {
  local_test_dir()
  m <- make_morph_matrix()
  morph_boxplots(m, plot_type = "boxplot", show_jitter = FALSE,
                 file_formats = c("pdf", "jpeg"), verbose = FALSE)
  files <- list.files(out_dir("Figs.boxplot"))
  expect_length(grep("^boxplot_.*\\.pdf$", files), 10)
  expect_length(grep("^boxplot_.*\\.jpeg$", files), 10)

  unlink("Figs.boxplot", recursive = TRUE)
  morph_boxplots(m, plot_type = "both", file_formats = "pdf", verbose = FALSE)
  files <- list.files(out_dir("Figs.boxplot"))
  expect_length(grep("^boxplot_", files), 10)
  expect_length(grep("^violin_", files), 10)
})

test_that("morph_boxplots skips sparse traits and fills missing colours", {
  local_test_dir()
  m <- make_morph_matrix()
  m[, "seed_length/SEEDlng"] <- NA
  m[1:2, "seed_length/SEEDlng"] <- c(1, 2)
  attr(m, "taxon_colors") <- c(alpha = "red")

  expect_warning(
    expect_message(
      morph_boxplots(m, plot_type = "both", min_n_total = 1,
                     file_formats = "pdf", verbose = TRUE),
      "Done"
    ),
    "No stored colour found for: beta, gamma"
  )
  files <- list.files(out_dir("Figs.boxplot"))
  # Too few observations for a boxplot, but a violin is still drawn
  expect_false(any(grepl("^boxplot_seed_length", files)))
  expect_true(any(grepl("^violin_seed_length", files)))

  unlink("Figs.boxplot", recursive = TRUE)
  expect_message(
    morph_boxplots(m, min_n_total = 10, file_formats = "pdf", verbose = TRUE) |>
      suppressWarnings(),
    "Skipping seed_length/SEEDlng: too few observations"
  )
})

test_that("morph_boxplots uses viridis colours when none are stored", {
  local_test_dir()
  m <- make_morph_matrix(colors = FALSE)
  morph_boxplots(m, file_formats = "pdf", verbose = FALSE)
  expect_length(list.files(out_dir("Figs.boxplot")), 10)
})

test_that("morph_boxplots validates its arguments", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_error(morph_boxplots(m, plot_type = "bars"))
  expect_error(morph_boxplots(m, show_jitter = "no"), "show_jitter")
  expect_error(morph_boxplots(m, min_n_total = 0), "min_n_total")
  expect_error(morph_boxplots(m, min_n_total = "a"), "min_n_total")
})

#_______________________________________________________________________________
# morph_similarity() ####

test_that("morph_similarity correlates trait pairs across all blocks", {
  local_test_dir()
  m <- make_morph_matrix()
  res <- quiet(morph_similarity(m, file_formats = c("pdf", "jpeg")))

  expect_type(res, "list")
  expect_named(res, c("by_block", "combined", "output_dir"))
  expect_setequal(names(res$by_block),
                  c("veg", "flo", "fru", "vegflo", "vegfru", "flofru", "vegflofru"))
  # 10 traits -> 45 pairs in the full combination
  expect_equal(nrow(res$by_block$vegflofru), 45)
  expect_true(all(abs(res$combined$cor) <= 1))
  expect_true(all(c("p_adjusted", "significant") %in% names(res$combined)))
  # Labels are cleaned in the combined table
  expect_false(any(grepl("/", res$combined$trait1)))

  od <- res$output_dir
  expect_true(file.exists(file.path(od, "Similarity_all_results.xlsx")))
  expect_true(file.exists(file.path(od, "Similarity_heatmap_veg.pdf")))
  expect_true(file.exists(file.path(od, "Similarity_heatmap_veg.jpeg")))
  expect_true(file.exists(file.path(od, "Similarity_veg.xlsx")))
})

test_that("morph_similarity supports single blocks, methods and custom folders", {
  local_test_dir()
  m <- make_morph_matrix()
  res <- quiet(morph_similarity(m, sim_blocks = "veg", sim_method = "spearman",
                                output_dir = "custom_out", file_formats = "pdf"))
  expect_named(res$by_block, "veg")
  expect_equal(res$output_dir, "custom_out")
  expect_true(file.exists("custom_out/Similarity_veg.xlsx"))

  res_k <- quiet(morph_similarity(m, sim_blocks = "flo", sim_method = "kendall",
                                  file_formats = "pdf", verbose = FALSE))
  expect_named(res_k$by_block, "flo")
  expect_error(morph_similarity(m, sim_method = "cosine"))
})

test_that("morph_similarity draws scatter plots for requested trait pairs", {
  local_test_dir()
  m <- make_morph_matrix()
  pairs <- list(c("petiole_length/PETIlng", "leaf_length/LEAFlng"))
  expect_output(
    res <- suppressMessages(
      morph_similarity(m, sim_blocks = "veg", sim_pairs = pairs,
                       file_formats = c("pdf", "jpeg"))
    ),
    "Scatter saved"
  )
  files <- list.files(res$output_dir, pattern = "^Scatter_")
  expect_length(files, 2)
})

test_that("morph_similarity warns about unusable scatter pairs", {
  local_test_dir()
  m <- make_morph_matrix()
  expect_warning(
    quiet(morph_similarity(m, sim_blocks = "veg",
                           sim_pairs = list(c("petiole_length/PETIlng", "nope")),
                           file_formats = "pdf")),
    "column\\(s\\) not found: nope"
  )
  expect_warning(
    quiet(morph_similarity(m, sim_blocks = "veg", sim_min_n = 1000,
                           sim_pairs = list(c("petiole_length/PETIlng",
                                              "leaf_length/LEAFlng")),
                           file_formats = "pdf")),
    "fewer than 1000 complete observations"
  )
})

test_that("morph_similarity skips blocks with too few traits or observations", {
  local_test_dir()
  m <- make_morph_matrix()
  attr(m, "base_cols") <- list(veg = "petiole_length/PETIlng",
                               flo = attr(m, "base_cols")$flo)
  expect_warning(
    res <- quiet(morph_similarity(m, sim_blocks = c("veg", "flo"),
                                  file_formats = "pdf")),
    "fewer than 2 traits"
  )
  expect_named(res$by_block, "flo")

  res2 <- quiet(morph_similarity(m, sim_blocks = "flo", sim_min_n = 1000,
                                 file_formats = "pdf"))
  expect_length(res2$by_block, 0)
  expect_null(res2$combined)
})

test_that(".clean_label strips codes and underscores", {
  expect_equal(.clean_label(c("leaf_length/LEAFlng", "taxon")),
               c("leaf length", "taxon"))
})
