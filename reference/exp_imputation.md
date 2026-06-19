# Standard exponential imputation

Imputes all missing values by independent draws from an exponential
distribution with rate 1.

## Usage

``` r
exp_imputation(X_miss)
```

## Arguments

- X_miss:

  A data set containing missing values.

## Value

A completed data set with all missing values replaced by draws from an
`Exp(1)` distribution.

## Examples

``` r
X <- random_mcar_data(100, 3)
X_imp <- exp_imputation(X)
```
