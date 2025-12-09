# Compute conditional cumulative incidences from fitted Cox models

These functions compute the conditional time- and covariate-specific
cumulative incidence for each row of `newdata`.

- `predict_from_model_0()` computes \\\psi_0(t_0; d, x) = 1 - S_0(d+t_0
  \mid x)/S_0(d+\tau \mid x)\\ by calling
  `predict(fit_0, newdata, type = "survival")` at times `d + t0` and
  `d + tau`.

- `predict_from_model_1()` computes \\\psi_1(t_0; d, x) = 1 - S_1(t_0
  \mid d, x)/S_1(\tau \mid d, x)\\ by calling
  `predict(fit_1, newdata, type = "survival")` at times `t0` and `tau`.

## Usage

``` r
predict_from_model_0(fit_0, exposure_time, t0, tau, newdata)

predict_from_model_1(fit_1, exposure_time, t0, tau, newdata)
```

## Value

The `newdata` data frame with three additional columns:

- `predict_from_model_0()`: `surv_0_d_plus_tau`, `surv_0_d_plus_t0`, and
  `psi_0_dx`

- `predict_from_model_1()`: `surv_1_tau`, `surv_1_t0`, and `psi_1_dx`
