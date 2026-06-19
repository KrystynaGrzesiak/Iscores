# Calculates IScores for multiple imputation functions

Calculates IScores for multiple imputation functions

## Usage

``` r
compare_Iscores(X, methods_list, score = c("energy_IScore", "DR_IScore"), ...)
```

## Arguments

- X:

  data containing missing values denoted with NA's.

- methods_list:

  a named list of imputing functions.

- score:

  a vector of names of scores to calculate. It can be `"energy_IScore"`
  and `"DR_IScore"`.

- ...:

  other arguments to be passed to
  [energy_IScore](https://krystynagrzesiak.github.io/Iscores/reference/energy_Iscore.md)
  or
  [DR_IScore](https://krystynagrzesiak.github.io/Iscores/reference/DR_IScore.md)

## Value

a vector of IScores for provided methods

## Examples

``` r
set.seed(111)
X <- random_mcar_data(100, 3, 0.2)
methods_list <- list(exp = exp_imputation,
                       norm = norm_imputation)
compare_Iscores(X, methods_list = methods_list, m = 2,
                n_proj = 10, n_trees_per_proj = 2 )
#> Calculating the energy_IScore for method exp ...
#> Calculating the energy_IScore for method norm ...
#> Calculating the DR_IScore for method exp ...
#> Calculating the DR_IScore for method norm ...
#>       score    score_name method
#> 1 0.7804594 energy_IScore    exp
#> 2 0.5733052 energy_IScore   norm
#> 3 2.9562649     DR_IScore    exp
#> 4 4.4875566     DR_IScore   norm
```
