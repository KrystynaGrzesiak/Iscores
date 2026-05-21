# Balancing of Classes

Balancing of Classes

## Usage

``` r
class.balancing(X_proj_complete, Y.proj, drawA, X_imp, ids.with.missing, vars)
```

## Arguments

- X_proj_complete:

  matrix with complete projected observations.

- Y.proj:

  matrix with projected imputed observations.

- drawA:

  vector of indices corresponding to current missingness pattern.

- X_imp:

  matrix of full imputed observations.

- ids.with.missing:

  vector of indices of observations with missing values.

- vars:

  vectors of variables in projection.

## Value

a list of new X_proj_complete and Y.proj.
