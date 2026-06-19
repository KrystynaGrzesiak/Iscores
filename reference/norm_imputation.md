# Standard normal imputation

Imputes all missing values by independent draws from a standard normal
distribution.

## Usage

``` r
norm_imputation(X_miss)
```

## Arguments

- X_miss:

  A data set containing missing values.

## Value

A completed data set with all missing values replaced by draws from a
\\N(0,1)\\ distribution.

## Examples

``` r
X <- random_mcar_data(100, 3)
X_imp <- norm_imputation(X)
```
