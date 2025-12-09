# Estimate observed distributions of exposure times and covariates

Computes two empirical probability distributions which are used to
marginalize time- and covariate- specific cumulative incidences:

- \\g(d \mid x)\\: the distribution of exposure times within each
  covariate subgroup of exposed individuals who remain at risk
  `immune_lag` days after exposure.

- \\p(x)\\: the distribution of covariates among exposed individuals who
  remain at risk `immune_lag` days after exposure.

Called internally by
[`nomatch()`](https://ewu16.github.io/nomatch/reference/nomatch.md).
Provides an example of the structure for user-specified
`custom_weights`.

## Usage

``` r
get_observed_weights(
  data,
  outcome_time,
  exposure,
  exposure_time,
  covariates,
  immune_lag
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

- immune_lag:

  Non-negative numeric value specifying the time after exposure
  (sometimes denoted by `tau`) that should be excluded from the risk
  evaluation period. This argument is primarily intended for vaccination
  exposures, where it is common to exclude the time after vaccination
  when immunity is still building. Time must be measured in the same
  units as that used for `outcome_time` and `exposure_time` (e.g. days)
  and should reflect the biological understanding of when
  vaccine-induced immunity develops (usually 1-2 weeks). For non-vaccine
  exposures, ` immune_lag` can be set to 0 (no delay period for
  evaluating risk).

## Value

A list with two data frames:

- `g_weights`: covariate-conditional exposure time probabilities (\\g(d
  \mid x)\\)

- `p_weights`: covariate probabilities (\\p(x)\\)

  Each data frame includes all variables in `covariates`, and a `prob`
  column with empirical probabilities. `g_weights` additionally includes
  a `<exposure_time>` column for exposure times.

## Examples

``` r
weights <- get_observed_weights(simdata, "Y", "V", "D_obs",
                                   c("x1", "x2"), immune_lag = 14)
str(weights)
#> List of 2
#>  $ g_weights:'data.frame':   1010 obs. of  4 variables:
#>   ..$ x1   : int [1:1010] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ x2   : int [1:1010] 5 5 5 5 5 5 5 5 5 5 ...
#>   ..$ D_obs: num [1:1010] 2 3 4 5 6 7 8 9 10 11 ...
#>   ..$ prob : num [1:1010] 0.0256 0.022 0.011 0.0586 0.0183 ...
#>  $ p_weights:'data.frame':   14 obs. of  3 variables:
#>   ..$ x1  : int [1:14] 0 0 0 0 0 0 0 1 1 1 ...
#>   ..$ x2  : int [1:14] 5 6 7 8 9 10 11 5 6 7 ...
#>   ..$ prob: num [1:14] 0.0675 0.0672 0.0705 0.064 0.0796 ...
```
