# Calculates Imputation Score for imputation function

Calculates Imputation Score for imputation function

## Usage

``` r
energy_IScore(
  X,
  imputation_func,
  X_imp = NULL,
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

This function relies on functions
[energy_Iscore_num](https://krystynagrzesiak.github.io/Iscores/reference/energy_Iscore_num.md)
and
[energy_Iscore_cat](https://krystynagrzesiak.github.io/Iscores/reference/energy_Iscore_cat.md).
Depending on the presence of factor-type data, these functions compute a
score either for purely numerical data or for mixed data types.

If you want to compute the score for numerical data, make sure that the
dataset does not contain any factor-type variables.

If you want to compute the score for categorical data, ensure that all
categorical variables are preserved as factors.

If your imputation method does not support categorical variables
represented as factors, implement a wrapper function that handles the
appropriate data type conversions before and after imputation.

## References

Näf, J., Grzesiak, K., and Scornet, E. (2025). How to rank imputation
methods? arXiv preprint.
[doi:10.48550/arXiv.2507.11297](https://doi.org/10.48550/arXiv.2507.11297)
.

## Examples

``` r
set.seed(111)
X <- random_mcar_data(100, 4)
imputation_func <- exp_imputation
energy_IScore(X, imputation_func)
#> [1] 0.8307913
#> attr(,"dat")
#>    column_id weight     score n_columns_used
#> X1         1 0.1875 0.9536624              1
#> X2         2 0.1824 0.7968742              1
#> X4         4 0.1600 0.7883675              1
#> X3         3 0.1476 0.7626066              1

X <-  random_mcar_mixed_data(100, 4, 2)
imputation_func <- median_mode_imputation
energy_IScore(X, imputation_func)
#> Factor variables detected.
#>              Calculating the energy-I-Score for mixed data.
#> [1] 0.8551675
#> attr(,"dat")
#>      column_id weight     score n_columns_used
#> col3         3 0.1771 0.8678183              1
#> col6         6 0.1716 1.0606602              1
#> col1         1 0.1476 0.7664758              1
#> col2         2 0.1411 0.7829227              1
#> col5         5 0.1275 0.8599947              1
#> col4         4 0.1131 0.7240102              1
```
