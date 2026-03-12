# Internal plotting function for effect panel plots

Internal plotting function for effect panel plots

## Usage

``` r
plot_effect_panel(
  plot_data,
  ci_type,
  alpha,
  effect = c("risk_ratio", "risk_difference", "relative_risk_reduction"),
  trt_0_label = "Unexposed",
  trt_1_label = "Exposed",
  colors = NULL
)
```

## Arguments

- plot_data:

  A data frame containing estimates to plot. Must include columns: `t0`
  (time), `estimate`, `term` (one of "cuminc_0", "cuminc_1",
  "relative_risk_reduction"), `<ci_type>_lower`, `<ci_type>_upper`, and
  `method`

- ci_type:

  Character string specifying the type of confidence interval band to
  plot, one of `"wald", "percentile", "simul"`, or `"none"`. Must choose
  a `ci_type` whose lower and upper bounds are already computed in
  `estimates` component of `x`. Uses `"wald"` by default if available.

- effect:

  The effect measure to plot next to the cumulative incidence plots.
  Either `"risk_ratio"`(default), `"relative_risk_reduction"` or
  `"risk_difference"`.
