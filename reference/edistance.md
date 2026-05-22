# Energy distance

Calculating energy distance/statistic.

## Usage

``` r
edistance(X, X_imp, rescale = FALSE)
```

## Arguments

- X:

  a complete original dataset

- X_imp:

  an imputed dataset

- rescale:

  a logical, indicating whether the returned value should be rescaled.
  Default to `FALSE`. See "details" section for more information.

## Details

This function uses the
[eqdist.e](https://rdrr.io/pkg/energy/man/eqdist.etest.html) function.
According to this implementation, by default, the function returns the
energy statistic which is given by \$\$E(X, Y) = \frac{nm}{n + m}
\hat{\varepsilon}{(X, Y)},\$\$ where \\\hat{\varepsilon}{(X, Y)}\\ is
the raw energy distance. To obtain raw energy distance use
`rescale = TRUE`.

## Examples

``` r
X <- matrix(rnorm(100), nrow = 25)
X_imp <- matrix(rnorm(100), nrow = 25)
edistance(X, X_imp)
#> E-statistic 
#>    2.905211 
```
