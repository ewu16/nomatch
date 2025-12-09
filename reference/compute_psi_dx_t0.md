# Compute time- and covariate- specific cumulative incidences

Wrapper that calls internal functions for predicting cumulative
incidences from fitted hazard models. Returns the predicted
exposure-specific cumulative incidences side by side for each time- and
covariate- pair in `newdata`.

## Usage

``` r
compute_psi_dx_t0(fit_0, fit_1, exposure_time, t0, tau, newdata)
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

## Value

A data frame with one row per row of `newdata` with the predicted
cumulative incidences and component survival probabilities:

- `psi_0_dx`, `surv_0_d_plus_tau`, `surv_0_d_plus_t0`

- `psi_1_dx`, `surv_1_tau`, `surv_1_t0`

## Details

Definitions of the cumulative incidences returned: \$\$\psi_0(t_0; d,x)
= 1 - S_0(d+t_0; x)\\/\\S_0(d+\tau; x)\$\$ \$\$\psi_1(t_0; d,x) = 1 -
S_1(t_0; d,x)\\/\\S_1(\tau; d,x)\$\$

where \\d\\ represents exposure time, \\x\\ represents covariates, and
\\S_v\\ represents the survival probability from hazard model for
exposure \\v\\.

## See also

[`predict_from_model_0()`](https://ewu16.github.io/nomatch/reference/predict_from_model_0.md),
[`predict_from_model_1()`](https://ewu16.github.io/nomatch/reference/predict_from_model_0.md)
