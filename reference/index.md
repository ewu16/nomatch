# Package index

## Main functions

Primary user-facing functions for proposed method

- [`nomatch()`](https://ewu16.github.io/nomatch/reference/nomatch.md) :
  Main function to estimate marginal cumulative incidences and derived
  effect measures without matching

- [`add_simultaneous_ci()`](https://ewu16.github.io/nomatch/reference/add_simultaneous_ci.md)
  :

  Add simultaneous confidence intervals to `nomatchfit` object

- [`plot(`*`<nomatchfit>`*`)`](https://ewu16.github.io/nomatch/reference/plot.nomatchfit.md)
  : Plot method for nomatchfit objects

## Supplementary functions for nomatch

Secondary functions for proposed method

- [`get_observed_weights()`](https://ewu16.github.io/nomatch/reference/get_observed_weights.md)
  : Estimate observed distributions of exposure times and covariates

- [`summarize_nomatch_population()`](https://ewu16.github.io/nomatch/reference/summarize_nomatch_population.md)
  :

  Estimate covariate means for population represented in
  [`nomatch()`](https://ewu16.github.io/nomatch/reference/nomatch.md)
  estimand

## Matching-based estimation

Functions for matching cohort design and estimation

- [`match_rolling_cohort()`](https://ewu16.github.io/nomatch/reference/match_rolling_cohort.md)
  : Match exposed and unexposed individuals using rolling cohort design
- [`matching()`](https://ewu16.github.io/nomatch/reference/matching.md)
  : Compute marginal cumulative incidence estimates for a matched cohort

## Utilities

Helper functions for data manipulation and formatting

- [`print(`*`<nomatchfit>`*`)`](https://ewu16.github.io/nomatch/reference/print.nomatchfit.md)
  :

  Print method for `nomatchfit` objects

- [`summary(`*`<nomatchfit>`*`)`](https://ewu16.github.io/nomatch/reference/summary.nomatchfit.md)
  :

  Summary method for `nomatchfit` objects

- [`estimates_to_df()`](https://ewu16.github.io/nomatch/reference/estimates_to_df.md)
  : Convert estimates from a fitted object into a tidy data frame

## Data

- [`simdata`](https://ewu16.github.io/nomatch/reference/simdata.md) :
  Simulated dataset
