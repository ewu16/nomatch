# Estimate covariate means for population represented in `nomatch()` estimand

Estimate covariate means for the weighted principal strata population
represented in
[`nomatch()`](https://ewu16.github.io/nomatch/reference/nomatch.md)
estimand. Since the principal strata populations are not directly
observed, simple summary statistics cannot be used. Instead, the means
are estimated as weighted means where the weights depend on the both the
original estimand weights and on predicted probabilities of belonging to
different principal strata. Currently only implemented for the case when
the estimand weights represent the observed distribution of exposure
days and covariates among the exposed.

## Usage

``` r
summarize_nomatch_population(data, fit)
```

## Arguments

- data:

  A data frame. The same data passed to
  [`nomatch()`](https://ewu16.github.io/nomatch/reference/nomatch.md).

- fit:

  A `nomatchfit` object created by
  [`nomatch()`](https://ewu16.github.io/nomatch/reference/nomatch.md).
  Must have been fit with `keep_models = TRUE` (the default) and
  `weights_source = "observed"`.

## Value

A data frame of the summary statistics with columns:

- variable:

  Covariate name

- label:

  Factor level label (empty string for numeric variables)

- mean:

  Mean (or weighted proportion for factor levels)

## Examples

``` r
fit <- nomatch(
  data = simdata,
  outcome_time = "Y",
  outcome_status = "event",
  exposure = "V",
  exposure_time = "D_obs",
  covariates = c("x1", "x2"),
  timepoints = seq(30, 180, by = 30),
  immune_lag = 14
)
summarize_nomatch_population(simdata, fit)
#>   variable label      mean
#> 1       x1       0.5035969
#> 2       x2       8.0659550
```
