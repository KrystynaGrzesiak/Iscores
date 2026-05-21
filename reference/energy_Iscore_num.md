# Calculates score for one imputation function

Calculates score for one imputation function

## Usage

``` r
energy_Iscore_num(
  X,
  imputation_func,
  X_imp = imputation_func(X),
  multiple = TRUE,
  N = 50,
  max_length = NULL,
  skip_if_needed = TRUE,
  scale = FALSE,
  n_cores = 1,
  silent = TRUE
)
```

## Arguments

- X:

  data containing missing values denoted with NA's.

- imputation_func:

  a function that imputes data.

- X_imp:

  imputed dataset of the same size as `X`. It's `NULL` by default
  meaning that it will be obtained by imputation of `X` using the
  `imputation_func`.

- multiple:

  a logical indicating whether provided imputation method is a multiple
  imputation approach (i.e. it generates different values to impute for
  each call). Default to TRUE. Note that if multiple equals to FALSE, N
  is automatically set to 1.

- N:

  a numeric value. Number of samples from imputation distribution H.
  Default to 50.

- max_length:

  Maximum number of variables \\X_j\\ to consider, can speed up the
  code. Default to `NULL` meaning that all the columns will be taken
  under consideration.

- skip_if_needed:

  logical, indicating whether some observations should be skipped to
  obtain complete columns for scoring. If FALSE, NA will be returned for
  column with no observed variable for training.

- scale:

  a logical value. If TRUE, each variable is scaled in the score.

- n_cores:

  a number of cores for parallelization.

- silent:

  logical indicating whether warnings and messages should be printed.

## Value

a numerical value denoting weighted Imputation Score obtained for
provided imputation function and a table with scores and weights
calculated for particular columns.

## References

This method is described in detail in:

Näf, J., Grzesiak, K., and Scornet, E. (2025a). How to rank imputation
methods? arXiv preprint arXiv:2507.11297
(<https://doi.org/10.48550/arXiv.2507.11297>).

## Examples

``` r
set.seed(123)

X <- Iscores:::random_mcar_data(n = 100, p = 4, ratio = 0.2)

imp_fun <- Iscores:::norm_imputation

sc <- Iscores::energy_Iscore_num(X = X, imputation_func = imp_fun, N = 5)
#> Error: 'energy_Iscore_num' is not an exported object from 'namespace:Iscores'

sc
#> Error: object 'sc' not found
```
