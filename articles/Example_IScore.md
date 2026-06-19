# Energy-I-Score: First Steps

### Introduction

The `Iscores` package provides tools for evaluating imputation methods
using Imputation Scores (IScores). In particular, the package implements
the DR I-Score and energy-I-Score, which measure the quality of imputed
datasets by comparing the relationships between observed and imputed
values. The methodology is described in detail in Näf et al. (2022) and
Näf et al. (2025).

The package supports:

- numerical datasets,
- mixed datasets containing both numerical and categorical variables,
- evaluation of a single imputation method,
- comparison of multiple imputation methods.

The main functions are:

- [`energy_IScore()`](https://krystynagrzesiak.github.io/Iscores/reference/energy_Iscore.md)
  — calculates the Energy-I-Score for a single imputation method,
- [`compare_Iscores()`](https://krystynagrzesiak.github.io/Iscores/reference/compare_Iscores.md)
  — compares several imputation methods using selected IScores.

This vignette presents the basic workflow for computing Imputation
Scores and comparing imputation approaches.

### Installation

The stable version of the package can be installed from CRAN (soon):

``` r

install.packages("Iscores")
```

The development version can be installed from GitHub:

``` r

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("KrystynaGrzesiak/Iscores")
```

After installation, load the package and set a random seed to ensure
reproducibility:

``` r

library(Iscores)
```

### Preparing data

The package expects a `data.frame` containing missing values represented
as `NA`.

Numerical variables should be stored as numeric vectors. Categorical
variables should be stored as factors if they are intended to be treated
as categorical during score calculation.

For demonstration purposes, we use some randomly generated data with
MCAR missings:

``` r

set.seed(10)

X <- random_mcar_data(100, 4)

head(X)
#>            X1         X2         X3         X4
#> 1  0.01874617 -0.7618043  1.2155138  1.5025446
#> 2 -0.18425254  0.4193754  0.3308765         NA
#> 3 -1.37133055         NA  1.3902751 -0.6306855
#> 4          NA  0.7115740  0.8720470  0.7923495
#> 5  0.29454513         NA -1.0808170         NA
#> 6  0.38979430  0.5631747  0.4958216  0.3227550
```

### Imputation Function

Before computing an Imputation Score, we first need to define the
**imputation method** that will be applied to the incomplete dataset.
The `Iscores` package is flexible and allows the user to evaluate any
imputation approach, provided that the imputation function satisfies a
few simple requirements.

Your imputation function:

1.  **Must accept** a dataset with missing values as its first argument.
2.  **Must return** a completed dataset with the same dimensions as the
    input data.
3.  **Should return a dataset without missing values**.
4.  Can represent:

- a simple custom imputation strategy,
- a wrapper around an external function, package or programming
  language.

For example, let’s define *zero imputation* below:

``` r

impute_zero <- function(X) { 
  
  X[is.na(X)] <- 0
  
  return(X) 
}
```

The function can now be passed directly to
[`energy_IScore()`](https://krystynagrzesiak.github.io/Iscores/reference/energy_Iscore.md),
[`DR_IScore()`](https://krystynagrzesiak.github.io/Iscores/reference/DR_IScore.md),
or
[`compare_Iscores()`](https://krystynagrzesiak.github.io/Iscores/reference/compare_Iscores.md).

### Calculating the energy-I-Score

Once an imputation function has been defined, we can evaluate it using
[`energy_IScore()`](https://krystynagrzesiak.github.io/Iscores/reference/energy_Iscore.md).

In the example below, we calculate the energy-I-Score for the
`impute_zero()` method.

``` r

sc <- energy_IScore(X = X, imputation_func = impute_zero)

sc
#> [1] 0.799712
#> attr(,"dat")
#>    column_id weight     score n_columns_used
#> X1         1 0.1771 0.7738732              1
#> X3         3 0.1771 0.7737987              1
#> X2         2 0.1539 0.7941482              1
#> X4         4 0.1411 0.8707366              1
```

The result is a single weighted score summarizing the imputation
performance across all variables containing missing values.

In addition, detailed information for each variable is returned as an
attribute of the result. This table contains:

- the variable-specific scores,
- aggregation weights,
- the number of variables used during training.

The table can be accessed using
[`attr()`](https://rdrr.io/r/base/attr.html):

``` r

attr(sc, "dat")
#>    column_id weight     score n_columns_used
#> X1         1 0.1771 0.7738732              1
#> X3         3 0.1771 0.7737987              1
#> X2         2 0.1539 0.7941482              1
#> X4         4 0.1411 0.8707366              1
```

Note that the weighted score is simply a weighted mean of scores for
particular variables:

``` r

sum(attr(sc, "dat")[["score"]] * attr(sc, "dat")[["weight"]]) / sum(attr(sc, "dat")[["weight"]])
#> [1] 0.799712
```

#### Important parameters

The
[`energy_IScore()`](https://krystynagrzesiak.github.io/Iscores/reference/energy_Iscore.md)
function exposes several parameters that control the scoring procedure.

##### Number of imputations: `N`

The parameter `N` controls how many times the missing part is re-imputed
during score estimation.

This parameter is mainly relevant for stochastic or multiple imputation
methods. For deterministic methods, setting `multiple = FALSE`
automatically forces `N = 1`.

``` r

energy_IScore(X = X, imputation_func = impute_zero, N = 5)
#> [1] 0.799712
#> attr(,"dat")
#>    column_id weight     score n_columns_used
#> X1         1 0.1771 0.7738732              1
#> X3         3 0.1771 0.7737987              1
#> X2         2 0.1539 0.7941482              1
#> X4         4 0.1411 0.8707366              1
```

> Note that hen `N = 1`, or when the imputation method always returns
> the same completed dataset, the predictive distribution is effectively
> represented by a single point estimate. In such cases, the
> energy-I-Score is computed from a degenerate empirical distribution,
> which may provide a less reliable approximation of uncertainty.
> Consequently, the score naturally favors imputation methods that
> generate realistic variability and sample well from the conditional
> distribution of the missing values, rather than methods that always
> return fixed imputations.

##### Limiting the number of scored variables: `max_length`

For datasets with many incomplete variables, score computation may
become time-consuming.

The `max_length` argument allows the user to limit the number of
variables used during score calculation. Variables with the largest
number of missing values are selected first.

By default, `max_length = NULL`, meaning that all incomplete variables
are included.

``` r

energy_IScore(X = X, imputation_func = impute_zero, max_length = 2)
#> [1] 0.773836
#> attr(,"dat")
#>    column_id weight     score n_columns_used
#> X1         1 0.1771 0.7738732              1
#> X3         3 0.1771 0.7737987              1
```

##### Handling incomplete training sets: `skip_if_needed`

Some variables may not have enough fully observed predictors available
for training.

In such situations:

- `skip_if_needed = TRUE` (default) removes a minimal number of
  observations to construct a valid training set,

- `skip_if_needed = FALSE` returns `NA` for variables where no complete
  predictors can be identified.

##### Scaling variables: `scale`

Setting `scale = TRUE` standardizes variables internally before score
calculation.

This can be useful when variables have very different numerical ranges,
preventing large-scale variables from dominating the score.

#### Energy-I-Score for mixed data

As mentioned earlier, categorical variables must be stored as factors in
order to be handled correctly by
[`energy_IScore()`](https://krystynagrzesiak.github.io/Iscores/reference/energy_Iscore.md).
If the input data contains at least one factor variable,
[`energy_IScore()`](https://krystynagrzesiak.github.io/Iscores/reference/energy_Iscore.md)
automatically switches to the mixed-data version of the score.
Therefore, users do not need to call a separate function for mixed
datasets.

Below we construct a simple mixed dataset containing both numerical and
categorical variables. Missing values are generated according to the
MCAR mechanism.

``` r

set.seed(10)

X_cat <- random_mcar_mixed_data(100, 4)

head(X_cat)
#>          col1       col2      col3       col4 col5
#> 1  0.01874617 -0.7618043        NA         NA    1
#> 2 -0.18425254         NA 0.3308765  0.5904095    3
#> 3 -1.37133055 -1.0399434 1.3902751 -0.6306855 <NA>
#> 4 -0.59916772  0.7115740 0.8720470  0.7923495    4
#> 5  0.29454513         NA        NA  0.1253846 <NA>
#> 6  0.38979430  0.5631747 0.4958216         NA    3
```

We use a simple median/mode imputation function. Numerical variables are
imputed with the median, while factor variables are imputed with the
most frequent category.

``` r

impute_mean_mode <- median_mode_imputation
```

The score can be calculated with the same public function as for
numerical data:

``` r

energy_IScore(X = X_cat, imputation_func = impute_mean_mode)
#> Factor variables detected.
#>              Calculating the energy-I-Score for mixed data.
#> [1] 0.8155521
#> attr(,"dat")
#>      column_id weight     score n_columns_used
#> col3         3 0.1924 0.7758805              1
#> col1         1 0.1771 0.7280336              1
#> col5         5 0.1539 0.9428090              1
#> col2         2 0.1411 0.8319434              1
#> col4         4 0.1411 0.8243028              1
```

Internally, categorical variables are transformed using one-hot encoding
and the score is then computed through multivariate energy distances.

For additional implementation details and methodological discussion, see
the vignette *“Energy-I-Score: Implementation Details”*.

### Calculating the DR-I-Score

The package also provides the
[`DR_IScore()`](https://krystynagrzesiak.github.io/Iscores/reference/DR_IScore.md)
function based on density-ratio estimation and random projection
forests.

``` r

sc_dr <- DR_IScore(X = X,
                   imputation_func = impute_zero,
                   m = 3,
                   n_proj = 10,
                   n_trees_per_proj = 2,
                   n_cores = 1)

sc_dr
#> [1] -8.327919
```

Unlike
[`energy_IScore()`](https://krystynagrzesiak.github.io/Iscores/reference/energy_Iscore.md),
which evaluates predictive distributions through scoring rules,
[`DR_IScore()`](https://krystynagrzesiak.github.io/Iscores/reference/DR_IScore.md)
compares the distributions of observed and imputed data using projected
random forests and density-ratio estimation.

The parameter `m` controls the number of imputed datasets generated by
the imputation method. Increasing `m` may improve score stability for
stochastic imputation procedures.

#### Parameters

The parameters `n_proj` and `n_trees_per_proj` control the complexity of
the random projection forests used internally:

- larger `n_proj` increases the number of random projections considered,
- larger `n_trees_per_proj` increases the number of trees grown for each
  projection.

Increasing these parameters may improve stability and precision of the
score, but also increases computational cost.

### Comparing multiple imputation methods

The
[`compare_Iscores()`](https://krystynagrzesiak.github.io/Iscores/reference/compare_Iscores.md)
function can be used to compare several imputation methods
simultaneously.

Below we define two additional imputation strategies from `mice`
package.

``` r

library(mice)


impute_mice_norm <- function(X) {
  imp <- mice(X, m = 1, method = "norm", maxit = 5, printFlag = FALSE)
  
  complete(imp)
}

impute_mice_rf <- function(X) {
  imp <- mice(X, m = 1, method = "rf", maxit = 5, printFlag = FALSE)
  
  complete(imp)
}
```

Now we place the methods in a named list:

``` r

methods_list <- list(zero = impute_zero,
                     mice_norm = impute_mice_norm,
                     mice_rf = impute_mice_rf)
```

We can now compare the methods using the energy-I-Score:

``` r

sc_comparison <- compare_Iscores(X = X,
                                 methods_list = methods_list,
                                 score = "energy_IScore",
                                 N = 10,
                                 silent = TRUE)
#> Calculating the energy_IScore for method zero ...
#> Calculating the energy_IScore for method mice_norm ...
#> Calculating the energy_IScore for method mice_rf ...

sc_comparison
#>       score    score_name    method
#> 1 0.7997120 energy_IScore      zero
#> 2 0.6247271 energy_IScore mice_norm
#> 3 0.6839281 energy_IScore   mice_rf
```

The resulting data frame contains one row per imputation method.

### Comparing methods using multiple scores

The package also allows simultaneous comparison using multiple scoring
rules.

``` r

comparison_all <- compare_Iscores(X = X,
                                  methods_list = methods_list,
                                  score = c("energy_IScore", "DR_IScore"),
                                  N = 10,
                                  m = 3,
                                  n_proj = 10,
                                  n_trees_per_proj = 2,
                                  silent = TRUE)
#> Calculating the energy_IScore for method zero ...
#> Calculating the energy_IScore for method mice_norm ...
#> Calculating the energy_IScore for method mice_rf ...
#> Calculating the DR_IScore for method zero ...
#> Calculating the DR_IScore for method mice_norm ...
#> Calculating the DR_IScore for method mice_rf ...

comparison_all
#>        score    score_name    method
#> 1  0.7997120 energy_IScore      zero
#> 2  0.6473443 energy_IScore mice_norm
#> 3  0.6753519 energy_IScore   mice_rf
#> 4 -9.2208219     DR_IScore      zero
#> 5  1.4280980     DR_IScore mice_norm
#> 6  2.4883888     DR_IScore   mice_rf
```

When multiple scores are requested, additional arguments passed to
[`compare_Iscores()`](https://krystynagrzesiak.github.io/Iscores/reference/compare_Iscores.md)
are automatically forwarded to the corresponding scoring functions.

In the example above:

- `N` and `silent` are arguments used by
  [`energy_IScore()`](https://krystynagrzesiak.github.io/Iscores/reference/energy_Iscore.md),
- `m`, `n_proj`, and `n_trees_per_proj` are arguments used by
  [`DR_IScore()`](https://krystynagrzesiak.github.io/Iscores/reference/DR_IScore.md).

Therefore, when combining several scoring rules, users should provide
the parameters required for each selected score.

### Summary of best practices

- Use `multiple = TRUE` with a genuinely multiple imputers and determine
  `N` for stable estimates.

- Consider `scale = TRUE` when mixing variables on different scales.

- Use `max_length` for quick experiments; remove it for final runs.

- Keep `skip_if_needed = TRUE` unless you explicitly want to flag
  unscorable columns with NA.

### Energy score

If you have access to the original dataset before imputation, you can
also use the energy distance as an additional evaluation metric. Our
package provides an easy-to-use wrapper
[`edistance()`](https://krystynagrzesiak.github.io/Iscores/reference/edistance.md)
around the
[`energy::eqdist.e`](https://rdrr.io/pkg/energy/man/eqdist.etest.html)
function from the energy package.

``` r


X_observed <- matrix(rnorm(2000), ncol = 4)  

X_miss <- X_observed
X_miss[runif(nrow(X_miss) * ncol(X_miss)) < 0.2] <- NA

edistance(X_observed, impute_zero(X_miss))
#> E-statistic 
#>    4.360287

edistance(X_observed, impute_mice_norm(X_miss))
#> E-statistic 
#>   0.5447785
```

## References

Näf, Jeffrey, Krystyna Grzesiak, and Erwan Scornet. 2025. *How to Rank
Imputation Methods?* <https://arxiv.org/abs/2507.11297>.

Näf, Jeffrey, Meta-Lina Spohn, Loris Michel, and Nicolai Meinshausen.
2022. *Imputation Scores*. <https://arxiv.org/abs/2106.03742>.
