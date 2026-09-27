#' Assess normality of morphological traits
#'
#' @description
#' Performs trait-wise Shapiro-Wilk normality tests on a morphometric data
#' matrix and generates tabular and graphical summaries for downstream
#' statistical analyses.
#'
#' The function is designed to work directly with the `analysis_data` object
#' returned by [morph_matrix_setting()], reading the trait blocks from its
#' `"base_cols"` attribute automatically.
#'
#' @details
#' For each morphological trait, `normality()` removes missing values and
#' performs a Shapiro-Wilk test of normality.
#'
#' Traits with fewer than three non-missing observations are skipped. Because
#' [stats::shapiro.test()] accepts at most 5000 observations, traits with more
#' than 5000 non-missing values are tested using a random sample of 5000
#' observations.
#'
#' A trait is classified as approximately normally distributed when its
#' Shapiro-Wilk p-value is greater than or equal to `norm_alpha`. Based on this
#' classification, the output includes a simple downstream recommendation:
#'
#' \itemize{
#'   \item `"ANOVA"` for traits classified as normal;
#'   \item `"Kruskal-Wallis"` for traits classified as non-normal.
#' }
#'
#' These recommendations are intended as exploratory guidance only. Selection
#' of an inferential test should also consider the study design, independence
#' of observations, homogeneity of variance, sample size, and other model
#' assumptions.
#'
#' The function generates an output directory named `Figs_Normality` with a
#' date-specific subdirectory. The following files are created:
#'
#' \itemize{
#'   \item a summary bar plot of Shapiro-Wilk p-values in the formats selected
#'   by `file_formats`;
#'   \item an Excel table containing the normality results;
#'   \item an Excel file combining all normality results;
#'   \item optionally, one Q-Q plot per trait in a `QQ_plots` subdirectory,
#'   using the formats selected by `file_formats`.
#' }
#'
#' At present, all traits listed in `base_cols` are combined into a single
#' analysis block named `"all traits"`. Trait blocks are read from `attr(analysis_data, "base_cols")`, which
#' [morph_matrix_setting()] attaches automatically. If `analysis_data` lacks
#' this attribute, the function stops with an informative error.
#'
#' @param analysis_data A numeric data frame or matrix containing morphological
#' traits as columns and specimens as rows. Typically obtained from the
#' `analysis_data` element returned by [morph_matrix_setting()].
#'
#' @param norm_alpha Numeric. Significance threshold used to classify traits
#' according to the Shapiro-Wilk test. Traits with p-values greater than or
#' equal to this value are classified as normal. Default is `0.05`.
#'
#' @param norm_qq Logical. Controls whether individual Q-Q plots are generated
#' for each tested trait. Default is `TRUE`.
#'
#' @param file_formats Character vector specifying which image formats to save.
#' One or more of `"pdf"` and `"jpeg"`. Default is `c("pdf", "jpeg")`.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages and
#' output information to the console.
#'
#' @return
#' Invisibly returns `NULL`.
#'
#' The function is primarily called for its side effects: it performs
#' normality tests, prints summaries to the console, and writes figures and
#' Excel files to disk.
#'
#' @seealso
#' [morph_matrix_setting()],
#' [stats::shapiro.test()],
#' [stats::qqnorm()]
#'
#' @examples
#' \dontrun{
#' # Build the morphometric matrix
#' res <- morph_matrix_setting(
#'   xlsx_path = "output_data/all_data.xlsx",
#'   taxon_col = "taxon",
#'   species_selected = NULL,
#'   trait_name_type = "code",
#'   veg_first = "petiole_length/PETIlng",
#'   flo_first = "inflorescence_length/INFLlng",
#'   fru_first = "fruit_stipe_length/FRSTlng"
#' )
#'
#' # Assess normality of all morphological traits
#' morph_normality(
#'   analysis_data = res,
#'   norm_alpha = 0.05,
#'   norm_qq = TRUE,
#'   verbose = TRUE
#' )
#' }
#'
#' @importFrom ggplot2 aes coord_flip element_blank element_text geom_col
#'   geom_hline geom_point geom_qq_line ggplot labs scale_fill_manual theme
#'   theme_bw
#' @importFrom cowplot save_plot
#' @importFrom openxlsx write.xlsx
#' @importFrom stats qqnorm shapiro.test
#'
#' @export
#'

morph_normality <- function(analysis_data,
                            norm_alpha = 0.05,
                            norm_qq = TRUE,
                            file_formats = c("pdf", "jpeg"),
                            verbose = TRUE) {

# ---- Validate input ------------------------------------------------------

if (!is.data.frame(analysis_data)) {
  stop("`analysis_data` must be a data frame.")
}

if (!"taxon" %in% names(analysis_data)) {
  stop("`analysis_data` must contain a `taxon` column")
}

if (is.null(rownames(analysis_data))) {
  stop("`analysis_data` must have specimen IDs as row names.")
}

base_cols <- attr(analysis_data, "base_cols")

if (!is.list(base_cols) || is.null(names(base_cols)) || any(names(base_cols) == "")) {
  stop(
    "`analysis_data` has no valid \"base_cols\" attribute. ",
    "Build `analysis_data` with `morph_matrix_setting()`, which attaches ",
    "this attribute automatically."
  )
}

file_formats <- match.arg(
  file_formats,
  choices = c("pdf", "jpeg"),
  several.ok = TRUE
)

if (!is.logical(norm_qq) || length(norm_qq) != 1L) {
  stop("`norm_qq` must be TRUE or FALSE.")
}

if (verbose) {
message("Running Normality Tests...")
}

# ---- Output folder -----------------------------------------------------------
output_dir <- file.path("Figs_Normality", format(Sys.time(), "%d%b%Y"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (norm_qq) dir.create(file.path(output_dir, "QQ_plots"),
                        recursive = TRUE, showWarnings = FALSE)

# ---- Resolve blocks ----------------------------------------------------------
run_blocks <- list("all traits" = unique(unlist(base_cols)))

# ---- Helper: test one block --------------------------------------------------
run_normality_block <- function(block_name, block_columns) {

if (verbose) {
  message("\n  --- Block: ", block_name, " ---")
}

valid_cols <- intersect(block_columns, names(analysis_data))
block_data <- analysis_data[, valid_cols, drop = FALSE]

results <- lapply(valid_cols, function(trait) {
  x <- block_data[[trait]]
  x <- x[!is.na(x)]

  if (length(x) < 3) return(NULL)

  # Shapiro-Wilk (max n = 5000)
  sw <- tryCatch(
    stats::shapiro.test(if (length(x) > 5000) sample(x, 5000) else x),
    error = function(e) NULL
  )
  if (is.null(sw)) return(NULL)

  data.frame(
    block = block_name,
    trait = trait,
    trait_label = gsub("_", " ", gsub("/.*", "", trait)),
    n = length(x),
    W_stat = round(sw$statistic, 4),
    p_value = round(sw$p.value, 4),
    normal = sw$p.value >= norm_alpha,
    recommended = ifelse(sw$p.value >= norm_alpha, "ANOVA", "Kruskal-Wallis")
  )
})

results <- do.call(rbind, Filter(Negate(is.null), results))
if (is.null(results) || nrow(results) == 0) return(invisible(NULL))

n_normal <- sum(results$normal)
n_nonnormal <- sum(!results$normal)

cat("  Traits tested :", nrow(results), "\n")
cat("  Normal (ANOVA):", n_normal, "\n")
cat("  Non-normal (KW):", n_nonnormal, "\n\n")

# ---- Summary bar chart -----------------------------------------------------
plot_df <- results[order(results$p_value), ]
plot_df$trait_label <- factor(plot_df$trait_label,
                              levels = rev(plot_df$trait_label))

p_summary <- ggplot2::ggplot(plot_df,
                             ggplot2::aes(x = .data[["trait_label"]], y = .data[["p_value"]], fill = .data[["normal"]])) +
  ggplot2::geom_col() +
  ggplot2::geom_hline(yintercept = norm_alpha, linetype = "dashed",
                      color = "red", linewidth = 0.8) +
  ggplot2::scale_fill_manual(values = c("TRUE"  = "#2166ac",
                                        "FALSE" = "#d6604d"),
                             labels = c("TRUE"  = "Normal (ANOVA)",
                                        "FALSE" = "Non-normal (KW)")) +
  ggplot2::coord_flip(clip = "off") +
  ggplot2::labs(x = NULL, y = "Shapiro-Wilk p-value", fill = "",
                title = paste0("Normality Test:", block_name),
                subtitle = paste0("Red line = p < ", norm_alpha,
                                  " | Normal: ", n_normal,
                                  " | Non-normal: ", n_nonnormal)) +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major.y = element_blank(),
                 axis.text.y        = element_text(size = 7))

if ("pdf" %in% file_formats) {
  cowplot::save_plot(
    file.path(
      output_dir,
      paste0("Normality_", block_name, ".pdf")
    ),
    p_summary,
    ncol = 1,
    nrow = 1,
    base_height = max(4, nrow(plot_df) * 0.3),
    base_aspect_ratio = 1.5
  )
}

if ("jpeg" %in% file_formats) {
  cowplot::save_plot(
    file.path(
      output_dir,
      paste0("Normality_", block_name, ".jpeg")
    ),
    p_summary,
    ncol = 1,
    nrow = 1,
    base_height = max(4, nrow(plot_df) * 0.3),
    base_aspect_ratio = 1.5
  )
}

# ---- Q-Q plots -------------------------------------------------------------
if (norm_qq) {
  qq_dir <- file.path(output_dir, "QQ_plots")

  for (i in seq_len(nrow(results))) {
    trait <- results$trait[i]
    x <- block_data[[trait]]
    x <- x[!is.na(x)]
    normal <- results$normal[i]
    trait_lab <- results$trait_label[i]

    qq <- stats::qqnorm(x, plot.it = FALSE)

    qq_df <- data.frame(
      theoretical = qq$x,
      sample = qq$y
    )

    p_qq <- ggplot2::ggplot(qq_df, ggplot2::aes(x = .data[["theoretical"]], y = .data[["sample"]])) +
      ggplot2::geom_point(size = 1.5, alpha = 0.6,
                          color = ifelse(normal, "#2166ac", "#d6604d")) +
      ggplot2::geom_qq_line(aes(sample = .data[["sample"]]),
                            color = "black", linetype = "dashed", linewidth = 0.6) +
      ggplot2::labs(x = "Theoretical Quantiles",
                    y = "Sample Quantiles",
                    title = trait_lab,
                    subtitle = paste0("W = ", results$W_stat[i],
                                      " | p = ", results$p_value[i],
                                      " | ", results$recommended[i])) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid = element_blank())

    safe_name <- gsub("[^a-zA-Z0-9]", "_", trait)

    if ("pdf" %in% file_formats) {
      cowplot::save_plot(
        file.path(
          qq_dir,
          paste0("QQ_", block_name, "_", safe_name, ".pdf")
        ),
        p_qq,
        ncol = 1,
        nrow = 1,
        base_height = 4,
        base_aspect_ratio = 1.2
      )
    }

    if ("jpeg" %in% file_formats) {
      cowplot::save_plot(
        file.path(
          qq_dir,
          paste0("QQ_", block_name, "_", safe_name, ".jpeg")
        ),
        p_qq,
        ncol = 1,
        nrow = 1,
        base_height = 4,
        base_aspect_ratio = 1.2
      )
    }
  }
}

# ---- Save Excel ------------------------------------------------------------
write.xlsx(results,
           file.path(output_dir, paste0("Normality_", block_name, ".xlsx")))

results
}

# ---- Run all blocks ----------------------------------------------------------
all_results <- Filter(Negate(is.null),
                      mapply(run_normality_block,
                             names(run_blocks), run_blocks,
                             SIMPLIFY = FALSE))

# ---- Save combined summary ---------------------------------------------------
if (length(all_results) > 0) {
  combined <- do.call(rbind, all_results)
  write.xlsx(combined, file.path(output_dir, "Normality_all_results.xlsx"))

  # Print overall recommendation
  cat("\n=== Overall Recommendation ===\n")
  cat("Normal traits    :", sum(combined$normal[!duplicated(combined$trait)]), "\n")
  cat("Non-normal traits:", sum(!combined$normal[!duplicated(combined$trait)]), "\n")
  cat("Suggested test   :",
      ifelse(mean(combined$normal) >= 0.5, "ANOVA", "Kruskal-Wallis"),
      "(majority of traits)\n\n")
}

if (verbose) {
  message("Done. Outputs saved in: ", output_dir)
}

}
