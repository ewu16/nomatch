# Get empirical probability distributions g(d\|x) and p(x) based on input data

Get empirical probability distributions g(d\|x) and p(x) based on input
data

## Usage

``` r
get_gp(df, outcome_time, exposure, exposure_time, covariates, tau)
```

## Arguments

- df:

  A data frame representing the target population for the estimated
  distributions

  - typically, population involves the "analysis-eligible" subset of a
    given data source

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

A list of the following:

- g_weights:

  A data frame containing the covariate-conditional day probabilities.
  Columns include `group_id` identifying the covariate group,
  `<exposure_time> specifying the day, `prob`specifying the covariate-conditional day probability, and each variables in`covariates`} \item{p_weights}{A data frame containing the covariate probabilities. Columns include `group_id` identifying the covariate group,`prob`specifying the covariate probability, and each variable in`covariates\`
