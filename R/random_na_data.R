
#' Create random example data for imputation
#'
#' This function generates some independent normal variables with missing values
#' according to MCAR (Missing Completely at Random) mechanism.
#'
#' @importFrom stats rnorm
#' @importFrom stats runif
#'
#' @param n number of samples
#' @param p number of variables
#' @param ratio fraction of missing values to generate
#'
#' @keywords internal
#'

random_mcar_data <- function(n, p, ratio = 0.2) {
  X <- matrix(stats::rnorm(n * p), nrow = n)
  X[stats::runif(n * p) <= ratio] <- NA
  data.frame(X)
}


#' Standard exponential imputation
#'
#' This is an example function that imputes from univariate Exp(1) distribution.
#'
#' @importFrom stats rexp
#'
#' @param X_miss a dataset with missing values
#'
#' @keywords internal
#'

random_imputation <- function(X_miss) {
  X_miss[is.na(X_miss)] <- stats::rexp(sum(is.na(X_miss)))

  X_miss
}


#' Standard normal imputation
#'
#' This is an example function that imputes from N(0, 1)
#'
#' @importFrom stats rnorm
#'
#' @param X_miss a dataset with missing values
#'
#' @keywords internal
#'

random_imputation <- function(X_miss) {
  X_miss[is.na(X_miss)] <- stats::rnorm(sum(is.na(X_miss)))

  X_miss
}

