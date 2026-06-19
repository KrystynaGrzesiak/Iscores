# Generate random mixed data with MCAR missing values

Generates a mixed dataset containing independent standard normal
variables and categorical variables, then introduces missing values
according to a Missing Completely at Random (MCAR) mechanism.

## Usage

``` r
random_mcar_mixed_data(n, p, n_fac = 1, ratio = 0.2)
```

## Arguments

- n:

  Number of observations.

- p:

  Number of numerical variables.

- n_fac:

  Number of categorical variables.

- ratio:

  Proportion of entries to replace with missing values.

## Value

A data frame containing `p` numerical variables and `n_fac` factor
variables with missing values.

## Examples

``` r
X <- random_mcar_mixed_data(100, 3, n_fac = 2, ratio = 0.2)
str(X)
#> 'data.frame':    100 obs. of  5 variables:
#>  $ col1: num  0.882 -0.815 0.59 0.64 -0.22 ...
#>  $ col2: num  0.954 0.0439 1.0666 NA 0.1515 ...
#>  $ col3: num  0.17 0.874 0.706 1.141 0.249 ...
#>  $ col4: Factor w/ 4 levels "1","2","3","4": 3 2 3 3 4 4 2 NA 1 NA ...
#>  $ col5: Factor w/ 4 levels "1","2","3","4": 2 1 NA 4 NA 4 1 4 NA 1 ...
```
