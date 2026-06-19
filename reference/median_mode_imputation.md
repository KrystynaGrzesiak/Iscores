# Median/mode imputation

Imputes numerical variables using their median and categorical variables
using their most frequent observed category.

## Usage

``` r
median_mode_imputation(X_miss)
```

## Arguments

- X_miss:

  A data set containing missing values.

## Value

A completed data set with all missing values imputed.

## Examples

``` r
X <- random_mcar_mixed_data(100, 3, n_fac = 1)
X_imp <- median_mode_imputation(X)
```
