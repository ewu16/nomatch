# Get analysis matched dataset from matched cohort

This function modifies the original matched data, preparing it for use
in analysis. Namely, this includes

- (if requested) censoring matched pairs at the time of the unvaccinated
  individual's vaccination

- creating the outcome time from matching index date to endpoint

- restricting to matched pairs in which both individuals are eligible
  `tau` days after the matching index date

## Usage

``` r
clean_matched_data(
  matched_data,
  outcome_time,
  outcome_status,
  exposure,
  exposure_time,
  tau
)
```

## Arguments

- matched_data:

  A data frame representing a matched cohort

- outcome_time:

  Character string representing the original outcome variable in
  `matched_data`

- outcome_status:

  Character string representing the original event variable in
  `matched_data`

- exposure:

  Character string representing the original vaccination status variable
  in `matched_data`

- exposure_time:

  Character string representing the original vaccination time variable
  in `matched_data`

- tau:

  The time excluded after vaccination to allow building up of immunity

## Value

A data frame representing the analysis-eligible matched dataset (a
subset of `matched data`). Contains all variables in `matched_data`,
plus the following variables

- match\_\<outcome_status\>:

  Event variable after adjusting for pair_censoring

- match\_\<outcome_time\>:

  Outcome variable after adjusting for pair_censoring

- match_T:

  Time to endpoint from matched index date after adjusting for
  pair_censoring
