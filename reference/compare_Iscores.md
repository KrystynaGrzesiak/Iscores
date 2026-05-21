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
X <- Iscores:::random_mcar_data(100, 3, 0.2)
methods_list <- list(exp = Iscores:::exp_imputation,
                       norm = Iscores:::norm_imputation)
compare_Iscores(X, methods_list = methods_list)
#> Calculating the energy_IScore for method exp ...
#> Calculating the energy_IScore for method norm ...
#> Calculating the DR_IScore for method exp ...
#> Calculating the DR_IScore for method norm ...
#>       score    score_name method
#> 1 0.5982370 energy_IScore    exp
#> 2 0.7892513 energy_IScore   norm
#> 3 3.1261742     DR_IScore    exp
#> 4 3.1731275     DR_IScore   norm
```
