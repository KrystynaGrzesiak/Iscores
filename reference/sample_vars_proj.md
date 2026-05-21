# Sampling of Projections

Sampling of Projections

## Usage

``` r
sample_vars_proj(ids_x_na, X, projection_function = NULL, normal_proj = TRUE)
```

## Arguments

- ids_x_na:

  a vector of indices corresponding to NA in the given missingness
  pattern.

- X:

  a matrix of the observed data containing missing values.

- projection_function:

  a function providing the user-specific projections.

- normal_proj:

  a boolean, if TRUE, sample from the NA of the pattern and additionally
  from the non-NA. If FALSE, sample only from the NA of the pattern.

## Value

a vector of variables corresponding to the projection.
