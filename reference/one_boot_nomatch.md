# Compute one bootstrap replicate for G-computation method

This function computes the point estimate for a bootstrap replicate.

## Usage

``` r
one_boot_nomatch(
  data,
  outcome_time,
  outcome_status,
  exposure,
  exposure_time,
  covariates,
  tau,
  timepoints,
  formula_0,
  formula_1,
  custom_gp_list
)
```

## Arguments

- data:

  A data frame with one row per individual containing the columns named
  in `outcome_time`, `outcome_status`, `exposure`, `exposure_time`, and
  `covariates`. Missing values for all columns except `exposure_time`
  are not allowed.

- outcome_time:

  Name of the follow-up time for the outcome of interest, i.e. time to
  either the event or right-censoring, whichever occurs first. Time
  should be measured from a chosen time origin (e.g. study start,
  enrollment, or age).

- outcome_status:

  Name of the event indicator for the outcome. The underlying column
  should be numeric (`1` = event, `0` = censored).

- exposure:

  Name of the exposure indicator. The underlying column should be
  numeric (`1` = exposed during follow-up, `0` = never exposed during
  follow-up).

- exposure_time:

  Name of the time to exposure, measured on the same time scale as that
  used for `outcome_time`. Must be a non-missing numeric value exposed
  individuals and must be set to `NA` for unexposed individuals.

- covariates:

  Character vector of covariates to adjust for when fitting the hazard
  models. Include all known, measured confounders of exposure and
  censoring. Covariates must be measured or defined at the chosen time
  origin.

- timepoints:

  Numeric vector specifying the timepoints at which to compute
  cumulative incidence and the derived effect measures. The timepoints
  should be expressed in terms of time since exposure and use the same
  units as that used for `outcome_time` and `exposure_time` (e.g. days).
  All values must be greater than ` immune_lag` and and should
  correspond to clinically meaningful follow-up durations, such as 30,
  60, or 90 days after exposure. A fine grid of timepoints (e.g.,
  `timepoints = (immune_lag + 1):100`) can be provided if cumulative
  incidence curves over time are desired.

## Value

The `pt_estimates` component returned by
[`get_one_nomatch()`](https://ewu16.github.io/nomatch/reference/get_one_nomatch.md),
for a bootstrap sample.
