# Compute overall cumulative incidences at a single timepoint

Compute overall cumulative incidences at a single timepoint

## Usage

``` r
compute_psi_bar_t0(fit_0, fit_1, exposure_time, t0, tau, newdata, gp_list)
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

- t0:

  Time since exposure at which to evaluate cumulative incidence.

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

Named numeric vector with `cuminc_0`, `cuminc_1`
