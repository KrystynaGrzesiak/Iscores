# Compute the density ratio score

Compute the density ratio score

## Usage

``` r
compute_drScore(object, Z = Z, n_trees_per_proj, n_proj)
```

## Arguments

- object:

  a crf object.

- Z:

  a matrix of candidate points.

- n_trees_per_proj:

  an integer, the number of trees per projection.

- n_proj:

  an integer specifying the number of projections.

## Value

a numeric value, the DR I-Score.
