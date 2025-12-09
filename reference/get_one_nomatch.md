# Internal function to compute cumulative incidence and effect measures

This is an internal function that performs the actual estimation of
cumulative incidences and derived effect measures using the
G-computation style approach. It is called by
[`nomatch()`](https://ewu16.github.io/nomatch/reference/nomatch.md) and
[`one_boot_nomatch()`](https://ewu16.github.io/nomatch/reference/one_boot_nomatch.md).
Users should typically call these functions rather than calling this
function directly. For historical reasons, this function handles more
complex inputs than exposed in the
[`nomatch()`](https://ewu16.github.io/nomatch/reference/nomatch.md)
interface. In particular, it includes an options to

- set censor time in hazard model for the exposed

- weight by the weights from a matched dataset and

- provide hazard model formulas or prefit objects

## Usage

``` r
get_one_nomatch(
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
  censor_time = max(timepoints),
  custom_gp_list = NULL,
  keep_models = TRUE,
  return_gp_list = TRUE
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

- censor_time:

  Time after exposure at which exposed individuals are administratively
  censored during model fitting. Default: `max(timepoints)`. This limits
  estimation to the observed period of interest and prevents
  extrapolation beyond the largest evaluation time. Typically, users
  should leave this at the default.

- keep_models:

  Logical; return the two fitted hazard models used to compute
  cumulative incidences? Default: `TRUE`.

- return_gp_list:

  Logical; return marginalizing weights? Default is TRUE.

- weights_source:

  Character string specifying how to construct marginalizing weights.
  One of:

  - `"observed"` (default): estimate weights from the observed data;

  - `"custom"`: use user-specified weights provided via `custom_weights`
    argument;

- return_matching:

  Logical; return matched datasets? Default is TRUE if
  `weights_source = "matched"`. When `weights_source != "matched"`, this
  is automatically set to `FALSE`.

## Value

List with components:

- pt_estimates:

  Matrix of point estimates with one row per evaluation time and one
  column per measure (`cuminc_0`, `cuminc_1`,
  `risk_ratio`,`relative_risk_reduction`).

- model_0, model_1:

  Fitted Cox models (if `keep_models = TRUE`)

- gp_list:

  Marginalizing distributions (if `return_gp_list = TRUE`)
