# Create random example data for imputation

This function generates some independent normal variables with missing
values according to MCAR (Missing Completely at Random) mechanism.

## Usage

``` r
random_mcar_mixed_data(n, p, n_fac = 1, ratio = 0.2)
```

## Arguments

- n:

  number of samples

- p:

  number of variables

- ratio:

  fraction of missing values to generate
