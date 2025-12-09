# Fit survival models to estimate exposure-specific hazards for G-computation approach

`fit_model_0()` fits a Cox model to estimate risk for unexposed
individuals on the original time scale. Includes all individuals,
censoring exposed individuals at their time of exposure. By default,
model is adjusted for by `<covariates>`, included as simple linear
terms.

`fit_model_1()` fits a Cox model to estimate risk for exposed
individuals on the time scale of time since exposure. Includes exposed
individuals who remain at risk `tau` days after exposure. Individuals
are additionally censored at `censor_time` days after exposure to avoid
extrapolation beyond the time period of interest. By default, model is
adjusted for `<covariates>`, included as simple linear terms, and
exposure time is included as a natural cubic spline with 4 degrees of
freedom.

## Usage

``` r
fit_model_0(
  data,
  outcome_time,
  outcome_status,
  exposure,
  exposure_time,
  formula_0 = NULL
)

fit_model_1(
  data,
  outcome_time,
  outcome_status,
  exposure,
  exposure_time,
  tau,
  formula_1 = NULL,
  censor_time = NULL
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

- formula_0:

  Optional right hand side of the formula for model 0. By default, uses
  `covariates`.

- tau:

  Delay period so that events before tau are not included.

- formula_1:

  Optional right hand side of the formula for model 1. By default, uses
  `covariates` plus natural spline of vaccination day (4 df).

- censor_time:

  Time after exposure at which exposed individuals are censored during
  model fitting to prevent extrapolation. By default, no censoring is
  applied.

## Value

A fitted `coxph` object with additional component `$data` containing the
analysis dataset used for fitting:

- For `fit_model_0()`: includes the survival tuple `(Y`, `event`) and
  covariates adjusted for in model, where `Y` is the time from time
  origin to first of endpoint, censoring or exposure time (for exposed
  individuals).

- For `fit_model_1()`: includes the survival tuple `(T1`, `event`),
  `<exposure_time>`, and covariates adjusted for in model, where `T1` is
  the time from exposure to endpoint or censoring, with additional
  censoring by `censor_time`. Only includes exposed individuals at risk
  `tau` days after exposure.
