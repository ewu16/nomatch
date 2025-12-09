# Package index

## Main functions

Primary user-facing functions

- [`nomatch()`](https://ewu16.github.io/nomatch/reference/nomatch.md) :
  Main function to estimate marginal cumulative incidences and derived
  effect measures without matching
- [`add_simultaneous_ci()`](https://ewu16.github.io/nomatch/reference/add_simultaneous_ci.md)
  : Add simultaneous confidence intervals to effectiveness fit
- [`plot(`*`<nomatchfit>`*`)`](https://ewu16.github.io/nomatch/reference/plot.nomatchfit.md)
  : Plot method for nomatchfit objects

## Proposed method estimation

Core estimation functions for the G-computation approach

- [`fit_model_0()`](https://ewu16.github.io/nomatch/reference/fit_model_0.md)
  [`fit_model_1()`](https://ewu16.github.io/nomatch/reference/fit_model_0.md)
  : Fit survival models to estimate exposure-specific hazards for
  G-computation approach
- [`get_observed_weights()`](https://ewu16.github.io/nomatch/reference/get_observed_weights.md)
  : Estimate observed distributions of exposure times and covariates
- [`compute_psi_dx_t0()`](https://ewu16.github.io/nomatch/reference/compute_psi_dx_t0.md)
  : Compute time- and covariate- specific cumulative incidences
- [`marginalize_psi_dx_t0()`](https://ewu16.github.io/nomatch/reference/marginalize_psi_dx_t0.md)
  : Marginalize conditional cumulative incidences
- [`compute_psi_bar_times()`](https://ewu16.github.io/nomatch/reference/compute_psi_bar_times.md)
  : Compute overall cumulative incidences at multiple timepoints

## Matching-based estimation

Functions for matching cohort design and estimation

- [`match_rolling_cohort()`](https://ewu16.github.io/nomatch/reference/match_rolling_cohort.md)
  : Match exposed and unexposed individuals using rolling cohort design
- [`matching()`](https://ewu16.github.io/nomatch/reference/matching.md)
  : Compute marginal cumulative incidence estimates for a matched cohort

## Utilities

Helper functions for data manipulation and formatting

- [`print(`*`<nomatchfit>`*`)`](https://ewu16.github.io/nomatch/reference/print.nomatchfit.md)
  : Print method for effectiveness fits
- [`summary(`*`<nomatchfit>`*`)`](https://ewu16.github.io/nomatch/reference/summary.nomatchfit.md)
  : Summary method for effectiveness fits
- [`estimates_to_df()`](https://ewu16.github.io/nomatch/reference/estimates_to_df.md)
  : Convert estimates from a fitted object into a tidy data frame
- [`compare_ve_fits()`](https://ewu16.github.io/nomatch/reference/compare_ve_fits.md)
  : Compare two effectiveness fit objects

## Data

- [`simdata`](https://ewu16.github.io/nomatch/reference/simdata.md) :
  Simulated dataset
