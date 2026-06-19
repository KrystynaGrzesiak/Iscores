# Generate random data with MCAR missing values

Generates a numerical dataset consisting of independent standard normal
variables and introduces missing values according to a Missing
Completely at Random (MCAR) mechanism.

## Usage

``` r
random_mcar_data(n, p, ratio = 0.2)
```

## Arguments

- n:

  Number of observations.

- p:

  Number of numerical variables.

- ratio:

  Proportion of entries to replace with missing values.

## Value

A data frame with `n` rows and `p` numerical variables containing
missing values.

## Examples

``` r
X <- random_mcar_data(100, 3, ratio = 0.2)
head(X)
#>           X1          X2        X3
#> 1         NA          NA 0.2358444
#> 2         NA -1.60646457        NA
#> 3  0.7598488          NA 1.5902015
#> 4 -0.3961216 -0.90981348 1.4488694
#> 5 -1.5779950          NA 2.1717999
#> 6  0.2172440 -0.02315266 0.2438218
```
