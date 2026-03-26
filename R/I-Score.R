
#' @title Calculates Imputation Score for imputation function
#'
#' @inheritParams energy_Iscore_num
#'
#' @param ... other arguments for the functions \link[Iscores]{energy_Iscore_num}
#' and \link[Iscores]{energy_Iscore_cat}.
#'
#' @details This function relies on functions \link[Iscores]{energy_Iscore_num}
#' and \link[Iscores]{energy_Iscore_cat}. Depending on the presence of factor-type
#' data, these functions compute a score either for purely numerical data or for
#' mixed data types.
#'
#' If you want to compute the score for numerical data, make sure that the
#' dataset does not contain any factor-type variables.
#'
#' If you want to compute the score for categorical data, ensure that all
#' categorical variables are preserved as factors.
#'
#' If your imputation method does not support categorical variables represented
#' as factors, implement a wrapper function that handles the appropriate data
#' type conversions before and after imputation.
#'
#' @return a numerical value denoting weighted Imputation Score obtained for
#' provided imputation function and a table with scores and weights calculated
#' for particular columns.
#'
#' @examples
#' set.seed(111)
#' X <- Iscores:::random_mcar_data(100, 4)
#' imputation_func <- Iscores:::exp_imputation
#' energy_IScore(X, imputation_func)
#'
#' X <-  Iscores:::random_mcar_mixed_data(100, 4, 2)
#' imputation_func <- Iscores:::median_mode_imputation
#' energy_IScore(X, imputation_func)
#'
#' @export
#'


energy_IScore <- function(X, imputation_func, X_imp = NULL, ...) {

  X <- as.data.frame(X, check.names = FALSE)

  categoricals <- sapply(X, function(i) is.factor(i))

  if(any(categoricals)) {
    message("There are some factor variables in the dataset.
    The energy-I-Score for mixed datasets will be calculated.")
    mixed <- TRUE
  } else {
    mixed <- FALSE
  }

  if(!is.function(imputation_func))
    stop("Imputation_func must be a function!")

  if(is.null(X_imp)) {
    X_imp <- try({imputation_func(X)})

    if(inherits(X_imp, "try-error") | any(is.na(X_imp)))
      stop("Errored imputing X using provided imputation_func!")
  }

  X_imp <- as.data.frame(X_imp, check.names = FALSE)

  if(mixed) {
    score <- energy_Iscore_cat(X, imputation_func, X_imp,...)
  } else {
    score <- energy_Iscore_num(X, imputation_func, X_imp, ...)
  }

  score
}


#' @title Calculates IScores for multiple imputation functions
#'
#' @inheritParams energy_Iscore_num
#'
#' @param score a vector of names of scores to calculate. It can be
#' \code{"energy_IScore"} and \code{"DR_IScore"}.
#' @param methods_list a named list of imputing functions.
#' @param ... other arguments to be passed to  \link[Iscores]{energy_IScore} or
#'  \link[Iscores]{DR_IScore}
#'
#'
#' @return a vector of IScores for provided methods
#'
#' @examples
#' set.seed(111)
#' X <- Iscores:::random_mcar_data(100, 3, 0.2)
#' methods_list <- list(exp = Iscores:::exp_imputation,
#'                        norm = Iscores:::norm_imputation)
#' compare_Iscores(X, methods_list = methods_list)
#'
#' @export
#'

compare_Iscores <- function(X,
                            methods_list,
                            score = c("energy_IScore", "DR_IScore"),
                            ...) {

  score <- match.arg(score, c("energy_IScore", "DR_IScore"), several.ok = TRUE)

  lapply(score, function(ith_score) {

    score_fun <- get(ith_score)
    methods <- names(methods_list)

    lapply(seq_along(methods_list), function(ith_method) {

      print(sprintf("Calculating %s for method: %s",
                    ith_score,
                    methods[ith_method]))

      imputation_func <- methods_list[[ith_method]]

      add_args <- list(...)

      args <- c(list(X = X, imputation_func = imputation_func),
                add_args[names(add_args) %in% names(formals(score_fun))])

      score <- do.call(score_fun, args)

      data.frame(score = as.numeric(score),
                 score_name = ith_score,
                 method = methods[ith_method])
    }) |>
      do.call(args = _, rbind)

  }) |>
    do.call(args = _, rbind)


}

