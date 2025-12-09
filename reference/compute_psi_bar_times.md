# Compute overall cumulative incidences at multiple timepoints

Wrapper that computes cumulative incidences at multiple timepoints.
Calls
[`compute_psi_bar_t0()`](https://ewu16.github.io/nomatch/reference/compute_psi_bar_t0.md)
internally for each timepoint. Models are fitted once before calling
this function to allow efficient evaluation at multiple timepoints
without refitting.

## Usage

``` r
compute_psi_bar_times(
  fit_0,
  fit_1,
  exposure_time,
  timepoints,
  tau,
  newdata,
  gp_list
)
```

## Arguments

- fit_0:

  A fitted model returned from
  [`fit_model_0()`](https://ewu16.github.io/nomatch/reference/fit_model_0.md)

- fit_1:

  A fitted model returned from
  [`fit_model_1()`](https://ewu16.github.io/nomatch/reference/fit_model_0.md)

- exposure_time:

  Name of the time-to-exposure variable in `newdata`. Used to compute
  \\\psi_0(t_0; d,x)\\ where \\d + \tau\\ and \\d + t_0\\ are needed.

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

- tau:

  Delay period. Internally used argument equivalent to user facing
  argument `immune_lag`

- newdata:

  New data at which to do predictions.

- gp_list:

  A list with two data frames:

  - **g_weights** Data frame of covariate-conditional exposure-day
    probabilities \\g(d \mid x)\\. Must include:

    - `group_id`: covariate group identifier

    - `<exposure_time>`: exposure time variable

    - `prob_g`: probability of exposure time given the covariates

    - all variables in `covariates`

  - **p_weights** Data frame of covariate probabilities \\p(x)\\. Must
    include:\\

    - `group_id`: covariate group identifier

    - `prob_p`: marginal probability of each covariate group

    - all variables in `covariates`

  Default is `NULL` in which case each row of `psi_dx` gets equal
  weight.

## Value

A matrix of estimates where the columns are the terms `cuminc_0` and
`cuminc_1`, and the rows are the time points of interest.

## Details

This function expects models already fitted via
[`fit_model_0()`](https://ewu16.github.io/nomatch/reference/fit_model_0.md)
and
[`fit_model_1()`](https://ewu16.github.io/nomatch/reference/fit_model_0.md).
It performs G-computation by predicting and marginalizing over
conditional risks.
