# Compute one bootstrap replicate of matching analysis method

Compute one bootstrap replicate of matching analysis method

## Usage

``` r
one_boot_matching(
  matched_data,
  outcome_time,
  outcome_status,
  exposure,
  exposure_time,
  tau,
  timepoints,
  limit_type = "fixed",
  data = NULL,
  id_name = NULL,
  matching_vars = NULL,
  replace = NULL,
  keep_boot_samples = TRUE
)
```

## Arguments

- matched_data:

  A data frame for the matched cohort created using
  [`match_rolling_cohort()`](https://ewu16.github.io/nomatch/reference/match_rolling_cohort.md).

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
  used for `outcome_time`. Must be a non-missing numeric value for
  exposed individuals and must be set to `NA` for unexposed individuals.

- timepoints:

  Numeric vector specifying the timepoints at which to compute
  cumulative incidence and the derived effect measures. The timepoints
  should be expressed in terms of time since exposure and use the same
  units as that used for `outcome_time` and `exposure_time` (e.g. days).
  All values must be greater than ` immune_lag` and and should
  correspond to clinically meaningful follow-up durations, such as 30,
  60, or 90 days after exposure. A fine grid of timepoints (e.g.,
  `timepoints = (immune_lag + 1):100`) can be provided if cumulative
  incidence curves over time are desired. By default, the sequence from
  `immune_lag + 1` to the maximum event time in the exposed group, by
  units of 1, is used.

## Value

The `pt_estimates` component returned by
[`get_one_matching()`](https://ewu16.github.io/nomatch/reference/get_one_matching.md),
for a bootstrap sample.
