#' Add random measures for preliminar morphometric analysis
#'
#' @author
#' João Dornelas & Domingos Cardoso
#'
#' @description
#' Creates and associates a data frame of morphological traits for flowering plants
#' with species occurrence matrices, including traits specific to selected families.
#' The resulting data frame can be exported to a specified directory as a .xlsx file.
#'
#' @usage
#' fill_missing_traits(
#'   df = NULL,
#'   start_col = 1
#' )
#'
#' @param df A data frame containing species occurrence records, previously
#' loaded from the working directory. If NULL, only the morphological traits
#' data frame is created.
#'
#' @param start_col Character. Choose either "trait full name" to load full
#' trait names or "trait code" to load abbreviated trait codes.
#'
#' @return A data frame combining occurrence data with the selected morphological traits,
#' exported as a .xlsx file.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' fill_missing_traits(
#'   base_df = my_occurrence_data,
#'   trait_name = "trait full name",
#'   quali_specific = "passifloraceae")
#' }
#'
#' @importFrom stats rnorm sd
#'
#' @export
#'

fill_missing_traits <- function(df = NULL,
                                start_col = 1) {

  for (i in start_col:ncol(df)) {
    x <- df[[i]]

    # Treat empty strings as NA for character/factor
    if (is.character(x) || is.factor(x)) x[x == ""] <- NA

    miss <- which(is.na(x))
    if (!length(miss)) next   # skip if no missing

    vals <- x[!is.na(x)]
    if (!length(vals)) next   # skip if column all NA

    if (is.numeric(x)) {
      # Numeric case → generate values around mean ± sd of existing
      mu <- mean(vals, na.rm = TRUE)
      sigma <- stats::sd(vals, na.rm = TRUE)
      if (is.na(sigma) || sigma == 0) sigma <- abs(mu) * 0.1  # fallback if no variation

      newvals <- stats::rnorm(length(miss), mean = mu, sd = sigma)

      # Keep integers if original was integer
      if (is.integer(df[[i]])) newvals <- as.integer(round(newvals))
    } else {
      # Character/factor case → resample existing values
      newvals <- sample(vals, length(miss), TRUE)
    }

    df[[i]][miss] <- newvals
  }
  df
}
