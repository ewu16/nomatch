# Compare two effectiveness fit objects

Plot cumulative incidence and effectiveness estimates for two different
methods using colors to distinguish methods

## Usage

``` r
compare_ve_fits(
  fit1,
  fit2,
  effect = c("risk_ratio", "risk_difference", "relative_risk_reduction"),
  labels = c("Method 1", "Method 2"),
  ci_type = NULL,
  colors = c("#F8766D", "#00BFC4")
)
```

## Arguments

- fit1:

  A nomatchfit object (typically from
  [`matching()`](https://ewu16.github.io/nomatch/reference/matching.md))

- fit2:

  A nomatchfit object (typically from
  [`nomatch()`](https://ewu16.github.io/nomatch/reference/nomatch.md))

- effect:

  The effect measure to plot next to the cumulative incidence plots.
  Either `"risk_ratio"`(default), `"relative_risk_reduction"` or
  `"risk_difference"`.

- labels:

  Character vector of length 2 providing labels for the two methods.
  Default is `c("Method 1", "Method 2")`.

- ci_type:

  Character string specifying the type of confidence interval to plot.
  One of "wald", "percentile", or "simul". If `NULL` (default), uses the
  CI type from `fit1`. If the object has `ci_type = "both"`, defaults to
  "wald".

- colors:

  Character vector of length 2 providing colors for the two methods.
  Default is `c("#F8766D", "#00BFC4")` (ggplot2's default red and cyan).

## Value

A ggplot2 object with three faceted panels showing cumulative incidence
and VE estimates for both methods.

## Details

Both `fit1` and `fit2` must have the same significance level (alpha).
The function will stop with an error if alphas differ.

For cumulative incidence panels, y-axis limits are shared across methods
to facilitate comparison. The effect measure panel uses free y-axis
scaling.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit both methods
fit_nomatch <- nomatch(data = simdata, ...)
fit_match <- matching(matched_data = matched_data, ...)
# Compare with custom labels
compare_ve_fits(
  fit_match,
  fit_nomatch,
  labels = c("Matching", "G-computation"))
} # }
```
