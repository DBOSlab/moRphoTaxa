# Tests for data preparation: traits_database, morph_add_traits(),
# fill_missing_traits(), morph_matrix_setting() and .pick_name()

#_______________________________________________________________________________
# traits_database ####

test_that("traits_database has the documented structure", {
  data("traits_database", package = "moRphoTaxa", envir = environment())
  expect_s3_class(traits_database, "data.frame")
  expect_named(traits_database,
               c("trait_type", "trait_scope", "trait_name", "trait_code"))
  expect_gt(nrow(traits_database), 0)
  expect_true(all(traits_database$trait_type %in% c("qualitative", "quantitative")))
  expect_true("default" %in% traits_database$trait_scope)
  expect_false(anyNA(traits_database$trait_code))
})

#_______________________________________________________________________________
# morph_add_traits() ####

test_that("morph_add_traits builds an empty trait template by default", {
  out <- morph_add_traits(save = FALSE, verbose = FALSE)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0)
  expect_true(all(c("leaf_length", "petiole_indumentum") %in% names(out)))
  expect_false(anyDuplicated(names(out)) > 0)
})

test_that("morph_add_traits supports trait codes and full/code labels", {
  codes <- morph_add_traits(trait_name = "trait code", save = FALSE, verbose = FALSE)
  both <- morph_add_traits(trait_name = "both", save = FALSE, verbose = FALSE)
  full <- morph_add_traits(trait_name = "trait full name", save = FALSE, verbose = FALSE)

  expect_true("LEAFtxt" %in% names(codes))
  expect_true("leaf_texture/LEAFtxt" %in% names(both))
  expect_equal(ncol(codes), ncol(full))
  expect_equal(ncol(both), ncol(full))
  expect_error(morph_add_traits(trait_name = "nonsense", save = FALSE))
})

test_that("morph_add_traits can drop default qualitative/quantitative traits", {
  quali <- morph_add_traits(quanti_default = FALSE, save = FALSE, verbose = FALSE)
  quanti <- morph_add_traits(quali_default = FALSE, save = FALSE, verbose = FALSE)
  none <- morph_add_traits(quali_default = FALSE, quanti_default = FALSE,
                           save = FALSE, verbose = FALSE)
  full <- morph_add_traits(save = FALSE, verbose = FALSE)

  expect_equal(ncol(quali) + ncol(quanti), ncol(full))
  expect_equal(ncol(none), 0)
})

test_that("morph_add_traits adds taxon-specific presets (case-insensitive)", {
  base <- morph_add_traits(save = FALSE, verbose = FALSE)
  pass <- morph_add_traits(quali_specific = "Passifloraceae",
                           quanti_specific = "PASSIFLORACEAE",
                           save = FALSE, verbose = FALSE)
  expect_gt(ncol(pass), ncol(base))

  expect_error(morph_add_traits(quali_specific = "notafamily", save = FALSE))
  expect_error(morph_add_traits(quanti_specific = "notafamily", save = FALSE))
})

test_that("every taxon-specific preset can be built", {
  presets <- c("araceae", "asteraceae", "bromeliaceae", "bryophyta",
               "cactaceae", "caesalpinioideae", "eriocaulaceae",
               "euphorbiaceae", "lecythidaceae", "lycopodiopsida", "moraceae",
               "orchidaceae", "papilionoideae", "passifloraceae", "pinophyta",
               "piperaceae", "poaceae", "polypodiopsida", "rubiaceae")
  for (p in presets) {
    out <- morph_add_traits(quali_specific = p, quanti_specific = p,
                            save = FALSE, verbose = FALSE)
    expect_gt(ncol(out), 0)
    expect_false(anyDuplicated(names(out)) > 0)
  }
})

test_that("non-angiosperm presets drop the angiosperm default traits", {
  expect_message(
    expect_message(
      moss <- morph_add_traits(quali_specific = "bryophyta",
                               quanti_specific = "bryophyta",
                               save = FALSE, verbose = TRUE),
      "quanti_specific"
    ),
    "quali_specific"
  )
  expect_false("petal_color" %in% names(moss))
  expect_false("leaf_length" %in% names(moss))

  # No message when defaults are switched off explicitly
  expect_no_message(
    morph_add_traits(quali_specific = "pinophyta", quanti_specific = "pinophyta",
                     quali_default = FALSE, quanti_default = FALSE,
                     save = FALSE, verbose = TRUE)
  )
})

test_that("morph_add_traits reports trait prefixes missing from the order list", {
  expect_message(
    morph_add_traits(quali_specific = "eriocaulaceae",
                     quanti_specific = "eriocaulaceae",
                     save = FALSE, verbose = TRUE),
    "character_order"
  )
})

test_that("morph_add_traits keeps base_df columns first and fills new traits with NA", {
  base_df <- data.frame(taxon = c("a", "b"), leaf_length = c(1, 2))
  out <- morph_add_traits(base_df = base_df, save = FALSE, verbose = FALSE)
  expect_equal(names(out)[1], "taxon")
  expect_equal(nrow(out), 2)
  expect_equal(out$leaf_length, c(1, 2))
  expect_true(all(is.na(out$petiole_indumentum)))
})

test_that("morph_add_traits writes the template to a user-given xlsx file", {
  local_test_dir()
  expect_message(
    morph_add_traits(save = TRUE, file_name = "my_traits.xlsx", verbose = TRUE),
    "exported"
  )
  expect_true(file.exists("my_traits.xlsx"))
  back <- openxlsx::read.xlsx("my_traits.xlsx")
  expect_true("leaf_length" %in% names(back))
})

test_that("morph_add_traits creates output_data/ when no file name is given", {
  local_test_dir()
  morph_add_traits(save = TRUE, verbose = FALSE)
  files <- list.files("output_data", pattern = "^trait_data_.*\\.xlsx$")
  expect_length(files, 1)
})

#_______________________________________________________________________________
# fill_missing_traits() ####

test_that("fill_missing_traits imputes numeric, integer and character columns", {
  set.seed(1)
  df <- data.frame(
    taxon = c("a", "b", "c", "d", "e"),
    num = c(1.5, NA, 3.2, NA, 2.2),
    int = c(1L, 2L, NA, 4L, 5L),
    chr = c("red", "", NA, "blue", "red"),
    fac = factor(c("x", NA, "y", "x", "y")),
    stringsAsFactors = FALSE
  )
  out <- fill_missing_traits(df, start_col = 2)

  expect_false(anyNA(out$num))
  expect_false(anyNA(out$int))
  expect_type(out$int, "integer")
  expect_false(anyNA(out$chr))
  expect_true(all(out$chr %in% c("red", "blue")))
  expect_false(anyNA(out$fac))
  # Observed values are left untouched
  expect_equal(out$num[c(1, 3, 5)], df$num[c(1, 3, 5)])
  expect_equal(out$taxon, df$taxon)
})

test_that("fill_missing_traits handles all-NA, complete and constant columns", {
  set.seed(2)
  df <- data.frame(
    all_na = c(NA_real_, NA_real_, NA_real_),
    complete = c(1, 2, 3),
    constant = c(5, 5, NA)
  )
  out <- fill_missing_traits(df, start_col = 1)
  expect_true(all(is.na(out$all_na)))
  expect_equal(out$complete, df$complete)
  # sd = 0 falls back to 10% of the mean instead of producing NaN
  expect_false(is.na(out$constant[3]))
})

test_that("fill_missing_traits skips columns before start_col", {
  df <- data.frame(id = c(NA, "b"), x = c(NA, 2))
  out <- fill_missing_traits(df, start_col = 2)
  expect_true(is.na(out$id[1]))
})

#_______________________________________________________________________________
# .pick_name() ####

test_that(".pick_name resolves full names and codes", {
  x <- c("leaf_length/LEAFlng", "taxon", "id")
  expect_equal(.pick_name(x, NULL), x)
  expect_equal(.pick_name(x, "full"), c("leaf_length", "taxon", "id"))
  expect_equal(.pick_name(x, "code"), c("LEAFlng", "taxon", "id"))
})

#_______________________________________________________________________________
# morph_matrix_setting() ####

test_that("morph_matrix_setting builds a numeric matrix with block attributes", {
  local_test_dir()
  write_raw_xlsx("raw.xlsx")

  m <- morph_matrix_setting(
    xlsx_path = "raw.xlsx",
    taxon_col = "species",
    veg_first = "petiole_length/PETIlng",
    flo_first = "inflorescence_length/INFLlng",
    fru_first = "fruit_stipe_length/FRSTlng"
  )

  expect_s3_class(m, "data.frame")
  expect_equal(names(m)[1], "taxon")
  # Rows without taxon or without any trait are dropped
  expect_equal(nrow(m), 28)
  expect_true(all(grepl("^ID[0-9]{4}$", rownames(m))))
  # Empty trait columns are dropped
  expect_false("empty_trait/EMPTlng" %in% names(m))
  # Decimal commas and blank strings are coerced
  expect_type(m$`petiole_width/PETIwid`, "double")
  expect_equal(m["ID0001", "petiole_width/PETIwid"], 3.5)
  expect_true(is.na(m["ID0002", "petiole_width/PETIwid"]))

  bc <- attr(m, "base_cols")
  expect_named(bc, c("veg", "flo", "fru"))
  expect_equal(bc$veg, c("petiole_length/PETIlng", "petiole_width/PETIwid",
                         "leaf_length/LEAFlng"))
  expect_equal(bc$fru, c("fruit_stipe_length/FRSTlng", "fruit_length/FRUTlng"))

  cols <- attr(m, "taxon_colors")
  expect_named(cols, c("alpha", "beta", "gamma"))
})

test_that("morph_matrix_setting renames traits to codes or full names", {
  local_test_dir()
  write_raw_xlsx("raw.xlsx")

  m_code <- morph_matrix_setting(
    "raw.xlsx", taxon_col = "species", trait_name_type = "code",
    veg_first = "petiole_length/PETIlng",
    flo_first = "inflorescence_length/INFLlng",
    fru_first = "fruit_stipe_length/FRSTlng"
  )
  expect_true(all(c("PETIlng", "INFLlng", "FRSTlng") %in% names(m_code)))
  expect_equal(attr(m_code, "base_cols")$flo[1], "INFLlng")

  m_full <- morph_matrix_setting(
    "raw.xlsx", taxon_col = "species", trait_name_type = "full",
    veg_first = "petiole_length/PETIlng",
    flo_first = "inflorescence_length/INFLlng",
    fru_first = "fruit_stipe_length/FRSTlng"
  )
  expect_true(all(c("petiole_length", "fruit_length") %in% names(m_full)))
})

test_that("morph_matrix_setting keeps selected species and existing ids", {
  local_test_dir()
  write_raw_xlsx("raw.xlsx", with_id = TRUE)

  m <- morph_matrix_setting(
    "raw.xlsx", taxon_col = "species", species_selected = c("alpha", "gamma"),
    veg_first = "petiole_length/PETIlng",
    flo_first = "inflorescence_length/INFLlng",
    fru_first = "fruit_stipe_length/FRSTlng"
  )
  expect_setequal(unique(m$taxon), c("alpha", "gamma"))
  expect_true(all(grepl("^SPEC", rownames(m))))
  expect_false("id" %in% names(m))
})

test_that("morph_matrix_setting validates block starts and colours", {
  local_test_dir()
  write_raw_xlsx("raw.xlsx")

  expect_error(
    morph_matrix_setting("raw.xlsx", taxon_col = "species",
                         veg_first = "nope/NOPE"),
    "veg_first column not found"
  )

  expect_warning(
    m <- morph_matrix_setting("raw.xlsx", taxon_col = "species",
                              veg_first = "petiole_length/PETIlng",
                              flo_first = "missing/MISS",
                              fru_first = "fruit_stipe_length/FRSTlng"),
    "could not be identified: flo"
  )
  expect_named(attr(m, "base_cols"), c("veg", "fru"))

  cols <- c(alpha = "red", beta = "blue", gamma = "green")
  m_col <- morph_matrix_setting("raw.xlsx", taxon_col = "species",
                                veg_first = "petiole_length/PETIlng",
                                flo_first = "inflorescence_length/INFLlng",
                                fru_first = "fruit_stipe_length/FRSTlng",
                                colors = cols)
  expect_equal(attr(m_col, "taxon_colors"), cols)

  expect_error(
    morph_matrix_setting("raw.xlsx", taxon_col = "species",
                         veg_first = "petiole_length/PETIlng",
                         colors = c("red", "blue", "green")),
    "named character vector"
  )
  expect_error(
    morph_matrix_setting("raw.xlsx", taxon_col = "species",
                         veg_first = "petiole_length/PETIlng",
                         colors = c(alpha = "red")),
    "missing an entry for: beta, gamma"
  )
})
