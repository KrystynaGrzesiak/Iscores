
#' @title Calculates Imputation Score for imputation function
#'
#' @inheritParams energy_Iscore
#'
#' @param ... other arguments for the functions \link[Iscores]{energy_Iscore}
#' and \link[Iscores]{energy_Iscore_cat}.
#'
#' @details This function relies on functions \link[Iscores]{energy_Iscore} and
#' \link[Iscores]{energy_Iscore_cat}. Depending on the presence of factor-type
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
#' imputation_func <- Iscores:::random_exp_imputation
#' energy_IScore(X, imputation_func)
#'
#' X <-  Iscores:::random_mcar_mixed_data(100, 4, 2)
#' imputation_func <- Iscores:::median_mode_imputation
#' energy_IScore(X, imputation_func)
#'
#' @export
#'


energy_IScore <- function(X, imputation_func, X_imp = NULL, ...) {

  mixed <- FALSE

  X <- as.data.frame(X, check.names = FALSE)

  categoricals <- sapply(X, function(i) is.factor(i))

  if(any(categoricals)) {
    message("There are some factor variables in the dataset.
    The energy-I-Score for mixed datasets will be calculated.")
    mixed <- TRUE
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
    score <- energy_Iscore(X, imputation_func, X_imp, ...)
  }

  score
}


#' @title Calculates IScores for multiple imputation functions
#'
#' @inheritParams energy_Iscore
#'
#' @param imputation_list a named list of imputing functions.
#'
#' @return a vector of IScores for provided methods
#'
#' @examples
#' set.seed(111)
#' X <- Iscores:::random_mcar_data(100, 4, 0.4)
#' imputation_list <- list(exp = Iscores:::exp_imputation,
#'                        norm = Iscores:::norm_imputation)
#' compare_Iscores(X, imputation_list)
#'
#' @export
#'

compare_Iscores <- function(X, imputation_list, ...) {

  if(length(names(imputation_list)) < length(imputation_list))
    stop("You must provide a named list of imutation methods!")

  methods <- names(imputation_list)
  iscores <- sapply(seq_along(methods), function(ith_method) {

    print(paste0("Calculating score for method: ", methods[ith_method]))

    X_imp <- imputation_list[[ith_method]](X)

    score <- energy_IScore(X, imputation_list[[ith_method]], ...)

    as.numeric(score)
  })

  names(iscores) <- methods
  sort(iscores)
}





