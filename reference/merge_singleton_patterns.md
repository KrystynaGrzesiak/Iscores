# Merge singleton missingness patterns

Merges missingness patterns that occur only once (singleton patterns)
into a single pattern. If the merged pattern already exists among the
current patterns, the corresponding groups of observations are combined.
Otherwise, a new pattern is created and appended.

## Usage

``` r
merge_singleton_patterns(patterns, groups, ind_singletons)
```

## Arguments

- patterns:

  A numeric matrix where each row represents a unique missingness
  pattern.

- groups:

  A list of integer vectors. Each element contains the indices of
  observations corresponding to a given pattern in `patterns`.

- ind_singletons:

  An integer vector indicating indices of patterns in `patterns` that
  occur only once.

## Value

A list with two elements:

- patterns:

  Updated matrix of unique missingness patterns.

- groups:

  Updated list of observation indices grouped by pattern.
