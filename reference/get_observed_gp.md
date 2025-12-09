# Estimate the marginalizing distribution in the observed or matched analysis-eligible populations.

`get_observed_gp` returns marginalizing distributions based on the
original observed data `get_matching_gp` returns marginalizing
distributions based on the matched-analysis data

## Usage

``` r
get_observed_gp(data, outcome_time, exposure, exposure_time, covariates, tau)
```

## Arguments

- outcome_time:

  Name of the follow-up time for the outcome of interest, i.e. time to
  either the event or right-censoring, whichever occurs first. Time
  should be measured from a chosen time origin (e.g. study start,
  enrollment, or age).

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

- tau:

  Numeric. Only vaccinated who are at-risk tau days after vaccination
  are included.

## Value

A list of two data frames containing the marginalizing distributions
returned from calls to
[`get_gp()`](https://ewu16.github.io/nomatch/reference/get_gp.md)
