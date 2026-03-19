
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

random_mcar_mixed_data <- function(n, p, n_fac = 1, ratio = 0.2) {
  X <- matrix(stats::rnorm(n * p), nrow = n)

  factors <- apply(
    matrix(sample(1:4, n * n_fac, replace = TRUE), nrow = n),
    2,
    function(x) factor(x)
  )

  factors[stats::runif(n * n_fac) <= ratio] <- NA
  X[stats::runif(n * p) <= ratio] <- NA

  X <- data.frame(X)
  factors <- data.frame(factors)

  X <- cbind(X, factors)
  colnames(X) <- paste0("col", 1:ncol(X))

  for (i in (ncol(X) - n_fac + 1):ncol(X)) {
    X[[i]] <- factor(X[[i]])
  }

  X
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

norm_imputation <- function(X_miss) {
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

exp_imputation <- function(X_miss) {
  X_miss[is.na(X_miss)] <- stats::rnorm(sum(is.na(X_miss)))

  X_miss
}

#' Standard normal imputation
#'
#' This is an example function that imputes from N(0, 1)
#'
#' @importFrom stats median
#'
#' @param X_miss a dataset with missing values
#'
#' @keywords internal
#'

median_mode_imputation <- function(X_miss) {
  for (col in names(X_miss)) {
    if (is.numeric(X_miss[[col]])) {
      med <- stats::median(X_miss[[col]], na.rm = TRUE)
      X_miss[[col]][is.na(X_miss[[col]])] <- med
    }
    else if (is.factor(X_miss[[col]])) {
      mode_val <- names(sort(table(X_miss[[col]]), decreasing = TRUE))[1]
      X_miss[[col]][is.na(X_miss[[col]])] <- mode_val
      X_miss[[col]] <- factor(X_miss[[col]])
    }
  }

  X_miss
}


