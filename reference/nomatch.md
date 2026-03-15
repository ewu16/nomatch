# Main function to estimate marginal cumulative incidences and derived effect measures without matching

`nomatch()` estimates marginal cumulative incidences under exposure and
no exposure using a G-computation approach. The method fits two
conditional hazard models- one for the exposed group and one for the
unexposed group- and uses these models to predict time- and
covariate-specific cumulative incidences. These time- and
covariate-specific cumulative incidences are then marginalized to
compute overall (marginal) cumulative incidences. By default, the
cumulative incidences are marginalized over the observed distribution of
exposure times and covariates among the exposed. The resulting
cumulative incidences can be summarized as risk ratios (RR = 1 -
risk_exposed/risk_unexposed), relative risk reduction (1 - RR), or risk
differences (RD = risk_exposed - risk_unexposed).

## Usage

``` r
nomatch(
  data,
  outcome_time,
  outcome_status,
  exposure,
  exposure_time,
  covariates,
  immune_lag = 0,
  timepoints = NULL,
  weights_source = c("observed", "custom"),
  custom_weights = NULL,
  ci_type = c("wald", "percentile", "both"),
  boot_reps = 0,
  alpha = 0.05,
  keep_models = TRUE,
  keep_boot_samples = TRUE,
  seed = NULL,
  formula_unexposed = NULL,
  formula_exposed = NULL
)
```

## Arguments

- data:

  A data frame with one row per individual containing the columns named
  in `outcome_time`, `outcome_status`, `exposure`, `exposure_time`, and
  `covariates`. Missing values in any of these columns except
  `exposure_time` are not allowed.

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
  used for `outcome_time`. Must be a non-missing numeric value for
  exposed individuals and must be set to `NA` for unexposed individuals.

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
  evaluating risk). Default: 0,

- timepoints:

  Numeric vector specifying the timepoints at which to compute
  cumulative incidence and the derived effect measures. The timepoints
  should be expressed in terms of time since exposure and use the same
  units as that used for `outcome_time` and `exposure_time` (e.g. days).
  All values must be greater than ` immune_lag` and and should
  correspond to clinically meaningful follow-up durations, such as 30,
  60, or 90 days after exposure. A fine grid of timepoints (e.g.,
  `timepoints = (immune_lag + 1):100`) can be provided if cumulative
  incidence curves over time are desired. By default, the sequence from
  `immune_lag + 1` to the maximum event time in the exposed group, by
  units of 1, is used.

- weights_source:

  Character string specifying the type of marginalizing weights to use
  for marginalization. Either:

  - `"observed"` (default): use the empirical distribution of exposure
    times and covariates among the exposed as the marginalizing weights.
    If `immune_lag > 0`, this is the distribution of exposure times and
    covariates among the exposed who remain at-risk `immune_lag` days
    after exposure. These weights provide close alignment with the
    weights implicitly used in matching.

  - `"custom"`: use the user-specified weights provided in the
    `custom_weights` argument. See **Marginalizing weights** section for
    more details.

- custom_weights:

  a named `list(g_weights, p_weights)` providing weights for
  marginalizing the time- and covariate-specific cumulative incidences.
  Must have the following format:

  - `g_weights`: data frame with the following columns:

    - columns named after each variable in `covariates`

    - `exposure_time` (time of exposure),

    - `prob` (probability of exposure at the given time within the
      covariate-group; should sum to 1 within each covariate-group)

  - `p_weights`: data frame with the following columns:

    - columns named after each variable in `covariates`

    - `prob` (probability of covariate-group; should sum to 1 over all
      covariate groups.)

- ci_type:

  Method for constructing pointwise bootstrap confidence intervals. One
  of `"wald"`, `"percentile"`, or `"both"`.

  - `"wald"` (default): Computes Wald-style intervals using bootstrap
    standard errors. See **Confidence intervals** section for details.

  - `"percentile"`: Computes percentile bootstrap intervals.

  - `"both"`: Computes and returns both Wald and percentile intervals.

- boot_reps:

  Number of bootstrap replicates for confidence intervals. Recommended
  to use at least 1000 for publication-quality results. Use smaller
  values (e.g., 10-100) for initial exploration. Default: `0` (no
  bootstrapping). Bootstrap procedure can be parallelizaed- see
  **Parallelization** section for details.

- alpha:

  Significance level for confidence intervals (Confidence level =
  100\*(1-alpha)%). Default: `0.05`.

- keep_models:

  Logical; return the two fitted hazard models used to compute
  cumulative incidences? Default: `TRUE`.

- keep_boot_samples:

  Logical; return bootstrap samples? Default: `TRUE`. Must be set to
  `TRUE` if user plans to use
  [`add_simultaneous_ci()`](https://ewu16.github.io/nomatch/reference/add_simultaneous_ci.md)
  to obtain simultaneous confidence intervals.

- seed:

  Integer seed for reproducible bootstrap results. Default: `NULL`.

- formula_unexposed, formula_exposed:

  Optional specification of the right-hand side of the formula to use
  for fitting the hazard model for the unexposed and exposed groups,
  respectively. Each accepts one of the following:

  - a one-sided formula object (e.g. `~ x1 + x2`)

  - a string representation of a formula (e.g. `"x1 + x2"` or
    `"~ x1 + x2"`), or a character vector of term names (e.g.
    `c("x1", "x2")`).

  The set of variables included in these formulas must be included in
  `covariates`. In `formula_exposed`, it is strongly recommended to
  model `exposure_time` flexibly (e.g. with a spline).

  By default, `formula_unexposed` is a main-effects formula using all
  `covariates` and `formula_exposed` is a main-effects formula using all
  `covariates` and a natural cubic spline of `exposure_time` (4 df).

  See **Modeling** section for details.

## Value

An object of class `nomatchfit` containing:

- estimates:

  Named list of matrices containing the cumulative incidence and effect
  estimates:

  `cuminc_0`

  :   marginal cumulative incidence under no exposure

  `cuminc_1`

  :   marginal cumulative incidence under exposure

  `risk_difference`

  :   `cuminc_1 - cuminc_0`

  `risk_ratio`

  :   `cuminc_1/cuminc_0`

  `relative_risk_reduction`

  :   `1 - risk_ratio`

  Each matrix has one row per value in `timepoints` and columns for the
  point estimate (`estimate`) and confidence limits
  (`{wald/percentile}_lower`, `{wald/percentile}_upper`), when
  applicable.

- model_0:

  (If `keep_models = TRUE`) Fitted hazard model for the unexposed group.
  See **Modeling** section for details.

- model_1:

  (If `keep_models = TRUE`) Fitted hazard model for the exposed group.
  See **Modeling** section for details.

- weights:

  List with dataframes `g_weights`, `p_weights` specifying the
  marginalizing weights used for averaging over exposure times and
  covariates.

- n_success_boot:

  Integer vector indicating the number of successful bootstrap
  replications per timepoint.

- boot_samples:

  (If `keep_boot_samples = TRUE`) Named list containing the `original`
  bootstrap estimates and, if Wald confidence intervals are used, the
  `transformed` bootstrap draws, which are the same bootstrap estimates
  on the scales used for inference. Both `original` and `transformed`
  are named lists containing the bootstrap estimates for each term. The
  bootstrap estimates for each term are stored as matrices, with rows
  indexing bootstrap replicates and columns indexing `timepoints`.

The `nomatchfit` object has methods for
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html). Use
[`add_simultaneous_ci()`](https://ewu16.github.io/nomatch/reference/add_simultaneous_ci.md)
to add simultaneous confidence intervals.

## Details

**Modeling.** Two Cox proportional hazards models are fit to estimate
exposure-specific cumulative incidences. For the unexposed model, the
outcome is modeled on the original time scale and includes all
individuals, with exposed individuals censored at their time of
exposure. For the exposed model, the outcome is modeled on the time
scale of time since exposure and includes individuals who were exposed
and who remained at risk `immune_lag` days after exposure. By default,
both models adjust for the specified covariates as linear, main-effect
terms to help control for confounding. The second model also flexibly
adjusts for exposure time (by default, as a natural cubic spline with 4
degrees of freedom) to capture time-varying background risk. Custom
formulas specifiying the right hand side of the formulas to be used in
these models can be provided through the optional arguments
`formula_unexposed` and `formula_exposed`. Predicted risks from both
models are then marginalized over the specified exposure-time and
covariate weights to obtain G-computation style cumulative incidence
estimates.

**Marginalizing weights.** When `weights_source = "observed"`, the
marginalizing weights are the empirical distributions of exposure times
and covariates among the exposed who remain at-risk `immune_lag` days
after exposure. These weights are returned in the `nomatchfit` object
under `weights`. They can also be obtained prior to the call to
`nomatch()` by calling
[`get_observed_weights()`](https://ewu16.github.io/nomatch/reference/get_observed_weights.md).

**Confidence intervals.** Wald and percentile confidence intervals are
constructed for cumulative incidence and effectiveness parameters at
each timepoint. The Wald pointwise confidence intervals are constructed
on transformed scales: \\\text{logit}\\ for cumulative incidence,
identity transformation for risk differences, \\\log{RR}\\ for risk
ratios, and \\\log{1 - RR}\\ for relative risk reduction using bootstrap
standard errors of the transformed estimates. These confidence intervals
are then back-transformed to the original scale. Wald p-values test the
null of no effect (RD = 0, RR = 1, 1-RR = 0) on the same transformed
scales used to construct the Wald confidence intervals. To obtain
simultaneous confidence intervals, use
[`add_simultaneous_ci()`](https://ewu16.github.io/nomatch/reference/add_simultaneous_ci.md)
after saving the original fit.

**Parallelization.** Bootstraps can be parallelized using the `future`
framework. Set a parallel plan before calling `nomatch()`: e.g.

    future::plan(future::multisession, workers = 4)
    fit <- nomatch(..., boot_reps = 1000, seed = 42)
    future::plan(future::sequential)  # reset when done

If no plan is set, bootstraps run sequentially. `multisession` works on
all operating systems and is recommended for most users. See the [future
package documentation](https://future.futureverse.org) for additional
plans and details on setup.

## Examples

``` r
# Fit nomatch using simulated data

fit <- nomatch(
  data = simdata,
  outcome_time = "Y",
  outcome_status = "event",
  exposure = "V",
  exposure_time = "D_obs",
  covariates = c("x1", "x2"),
  timepoints = seq(30, 90, by = 30),
  immune_lag = 14,
  boot_reps = 5,
  seed = 123
)
#> Bootstrapping 5 samples...
#> Bootstrap completed in 1.05 secs

# View basic results
fit$estimates
#> $cuminc_0
#>      estimate  wald_lower wald_upper wald_pval wald_n
#> 30 0.01163891 0.009815874 0.01379582        NA      5
#> 60 0.03790915 0.032558421 0.04409914        NA      5
#> 90 0.05857293 0.052444528 0.06536807        NA      5
#> 
#> $cuminc_1
#>       estimate  wald_lower wald_upper wald_pval wald_n
#> 30 0.006211499 0.003822016 0.01007975        NA      5
#> 60 0.022820770 0.019600908 0.02655524        NA      5
#> 90 0.035207076 0.028754881 0.04304289        NA      5
#> 
#> $risk_difference
#>        estimate   wald_lower   wald_upper    wald_pval wald_n
#> 30 -0.005427414 -0.009417768 -0.001437060 7.680255e-03      5
#> 60 -0.015088377 -0.022777832 -0.007398923 1.201233e-04      5
#> 90 -0.023365858 -0.031704713 -0.015027002 3.976111e-08      5
#> 
#> $risk_ratio
#>     estimate wald_lower wald_upper    wald_pval wald_n
#> 30 0.5336838  0.2976856  0.9567758 3.500360e-02      5
#> 60 0.6019858  0.4669828  0.7760178 8.960706e-05      5
#> 90 0.6010810  0.4851894  0.7446542 3.194904e-06      5
#> 
#> $relative_risk_reduction
#>     estimate wald_lower wald_upper    wald_pval wald_n
#> 30 0.4663162 0.04322417  0.7023144 3.500360e-02      5
#> 60 0.3980142 0.22398220  0.5330172 8.960706e-05      5
#> 90 0.3989190 0.25534575  0.5148106 3.194904e-06      5
#> 
```
