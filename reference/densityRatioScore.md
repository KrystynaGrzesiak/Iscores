# Computation of the density ratio score

Computes the density ratio score using a random forest model based on
random projections.

## Usage

``` r
densityRatioScore(
  X,
  X_imp,
  pattern = NULL,
  n_proj = 10,
  n_trees_per_proj = 1,
  projection_function = NULL,
  min_node_size = 1,
  normal_proj = TRUE
)
```

## Arguments

- X:

  A numeric matrix of observed data that may contain missing values
  denoted by `NA`.

- X_imp:

  A numeric matrix of imputed values with the same dimensions as `X`.

- pattern:

  A vector or pattern indicating the missingness structure.

- n_proj:

  An integer specifying the number of random projections.

- n_trees_per_proj:

  An integer specifying the number of trees grown per projection.

- projection_function:

  A function that generates user-defined projections.

- min_node_size:

  An integer specifying the minimum number of observations in a terminal
  node (leaf) of each tree.

- normal_proj:

  Logical. If `TRUE`, sampling is performed from both missing (NA) and
  observed values. If `FALSE`, sampling is performed only from missing
  (NA) values.

## Value

An object representing a fitted random forest model based on random
projections.

## Details

The method builds multiple random forests on projected versions of the
data to estimate the density ratio between observed and imputed
distributions.
