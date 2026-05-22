# energy-I-Score for imputation of mixed data (categorical and numerical)

energy-I-Score for imputation of mixed data (categorical and numerical)

## Usage

``` r
energy_Iscore_cat(
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

## Details

The categorical variables should be stored as factors. If you need
additional conversion of the data (for example one-hot encoding) for
imputation, please, implement everything within `imputation_func`
parameter. You can use `miceDRF:::onehot_to_factor` and
`miceDRF:::factor_to_onehot` functions.

## References

This method is described in detail in:

Näf, J., Grzesiak, K., and Scornet, E. (2025). How to rank imputation
methods? arXiv preprint.
[doi:10.48550/arXiv.2507.11297](https://doi.org/10.48550/arXiv.2507.11297)
.

## Examples

``` r
set.seed(123)

X <- Iscores:::random_mcar_mixed_data(n = 100, p = 3, n_fac = 1, ratio = 0.2)

imp_fun <- Iscores:::median_mode_imputation

sc <- Iscores:::energy_Iscore_cat(X = X, imputation_func = imp_fun, N = 5)

sc
#> [1] 0.7699801
#> attr(,"dat")
#>      column_id weight     score n_columns_used
#> col1         1 0.1824 0.7096268              1
#> col2         2 0.1659 0.7301908              1
#> col4         4 0.1600 0.9355567              1
#> col3         3 0.1539 0.7122622              1
```
