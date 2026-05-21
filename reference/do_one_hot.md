# Convert a factor vector to one-hot encoding

Converts a factor vector into a one-hot encoded matrix with one column
per factor level.

## Usage

``` r
do_one_hot(vec)
```

## Arguments

- vec:

  A factor vector to be encoded.

## Value

A numeric matrix with one row per element of \`vec\` and one column per
factor level. Column names are prefixed with \`"level\_"\`.

## Details

Missing values in \`vec\` are preserved as rows containing \`NA\`
values.
