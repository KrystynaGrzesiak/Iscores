# Extract and group missing-data patterns

Identifies unique missingness patterns in a data matrix and groups
observations according to these patterns. If more than one pattern
occurs only once, such singleton patterns are merged into a single
group.

## Usage

``` r
get_pattern_data(X)
```

## Arguments

- X:

  A matrix or data frame that may contain missing values.

## Value

A list with three elements:

- patterns:

  A matrix of unique missingness patterns.

- groups:

  A list of integer vectors giving row indices for each pattern.

- average_diff:

  A logical indicating whether singleton patterns were merged.

## Details

Missingness patterns are represented by a logical matrix obtained from
`is.na(X)`. Only rows containing at least one missing value are used to
define the unique patterns.

If more than one pattern is represented by a single observation, these
singleton patterns are merged using
[`merge_singleton_patterns()`](https://krystynagrzesiak.github.io/Iscores/reference/merge_singleton_patterns.md).
