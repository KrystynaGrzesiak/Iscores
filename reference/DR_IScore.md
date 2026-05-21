# Compute the imputation KL-based scoring rules

Compute the imputation KL-based scoring rules

## Usage

``` r
DR_IScore(
  X,
  imputation_func = NULL,
  X_imp = NULL,
  m = 5,
  n_proj = 100,
  n_trees_per_proj = 5,
  min_node_size = 10,
  n_cores = 1,
  projection_function = NULL,
  ...
)
```

## Arguments

- X:

  data containing missing values denoted with NA's.

- imputation_func:

  an imputing function. If `NULL`, please provide imputed datasets
  `X_imp` and `m`.

- X_imp:

  a list of imputed datasets. If `NULL` it will be obtained using
  `imputation_func`.

- m:

  the number of multiple imputations to consider, default to 5.

- n_proj:

  an integer specifying the number of projections to consider for the
  score.

- n_trees_per_proj:

  an integer, the number of trees per projection.

- min_node_size:

  the minimum number of nodes in a tree.

- n_cores:

  an integer, the number of cores to use.

- projection_function:

  a function providing the user-specific projections.

- ...:

  used for compatibility

## Value

a vector made of the scores for each imputation method.

## Examples

``` r
set.seed(111)
X <- Iscores:::random_mcar_data(100, 3, 0.2)
imputation_func <- Iscores:::exp_imputation
DR_IScore(X, imputation_func)
#> [1] 2.985856
```
