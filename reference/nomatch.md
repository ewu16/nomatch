# Main function to estimate marginal cumulative incidences and derived effect measures without matching

`nomatch()` estimates marginal cumulative incidences under exposure and
no exposure using a G-computation approach. The method fits two
conditional hazard models- one for each exposure type- and uses these
models to predict time- and covariate-specific cumulative incidences.
These cumulative incidences are then marginalized to compute overall
(marginal) cumulative incidences. By default, the cumulative incidences
are marginalized over the observed distribution of exposure times and
covariates among the exposed. The resulting cumulative incidences can be
summarized as risk ratios (RR = 1 - risk_exposed/risk_unexposed),
relative risk reduction (1 - RR), or risk differences (RD =
risk_exposed - risk_unexposed).

## Usage

``` r
nomatch(
  data,
  outcome_time,
  outcome_status,
  exposure,
  exposure_time,
  covariates,
  immune_lag,
  timepoints,
  weights_source = c("observed", "custom"),
  custom_weights = NULL,
  ci_type = c("wald", "percentile", "both"),
  boot_reps = 0,
  alpha = 0.05,
  keep_models = TRUE,
  keep_boot_samples = TRUE,
  n_cores = 1,
  formula_unexposed = NULL,
  formula_exposed = NULL
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

- weights_source:

  Character string specifying the type of marginalizing weights to use.
  Either:

  - `"observed"` (default): use the empirical distribution of exposure
    times and covariates among the exposed as the marginalizing weights.
    This provides close alignment with the weights implicitly used in
    matching.

  - `"custom"`: use the user-specified weights provided in the
    `custom_weights` argument.

- custom_weights:

  a `list(g_weights, p_weights)` providing weights for marginalizing the
  time- and covariate-specific cumulative incidences. Must have the
  following format:

  - `g_weights`: data frame with columns

    - all variables in `covariates`

    - `exposure_time` (time of exposure),

    - `prob` (probability of exposure at the given time within the
      covariate-group; should sum to 1 within each covariate-group)

  - `p_weights`: data frame with columns

    - all variables in `covariates`

    - `prob` (probability of covariate-group; should sum to 1 over all
      covariate groups.)

- ci_type:

  Method for constructing bootstrap confidence intervals. One of
  `"wald"`, `"percentile"`, or `"both"`.

  - `"wald"` (default): Computes Wald-style intervals using bootstrap
    standard errors. See **Confidence intervals** section for details.

  - `"percentile"`: Computes percentile bootstrap intervals.

  - `"both"`: Computes and returns both Wald and percentile intervals.

- boot_reps:

  Number of bootstrap replicates for confidence intervals. Recommended
  to use at least 1000 for publication-quality results. Use smaller
  values (e.g., 10-100) for initial exploration. Default: `0` (no
  bootstrapping).

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

- n_cores:

  Integer; parallel cores for bootstrapping passed to
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) as
  `mc.cores`. On Unix-like OS only; not available on Windows. Default:
  `1`.

- formula_unexposed:

  A character specification of the right-hand side of the formula to use
  for fitting the hazard model for the unexposed. The set of variables
  in the formula must be identical to the set of variables in
  `covariates`.

- formula_exposed:

  A character specification of the right-hand side of the formula to use
  for fitting the hazard model for the exposed. The set of variables in
  the formula must contain the variables in `covariates` and
  `exposure_time`. It is recommended that `exposure_time` be modeled
  flexibly (e.g. using splines).

## Value

An object of class `nomatchfit` containing:

- estimates:

  Named list of matrices containing the cumulative incidence and effect
  estimates.

  `cuminc_0`

  :   marginal cumulative incidence under no exposure

  `cuminc_1`

  :   marginal cumulative incidence under exposure

  `risk_difference`

  :   `cuminc_1 - cuminc_0`

  `risk_ratio`

  :   `cuminc_1/cuminc_0`

  `vaccine_effectivess`

  :   `1 - risk_ratio`

  Each matrix has one row per value in `timepoints` and columns
  including the point estimate (`estimate`) and, when requested,
  confidence limits of the form (`{wald/percentile}_lower`,
  `{wald/percentile}_upper`).

- weights:

  List with dataframes `g_weights`, `p_weights` specifying the
  marginalizing weights used for averaging over exposure times and
  covariates.

- model_0:

  Fitted hazard model for the unexposed group. See **Modeling** section
  for details.

- model_1:

  Fitted hazard model for the exposed group. See **Modeling** section
  for details.

- n_success_boot:

  Integer vector indicating the number of successful bootstrap
  replications per timepoint.

- boot_samples:

  (If `keep_boot_samples = TRUE`) Named list of bootstrap draws (stored
  as matrices) for each term. Rows index bootstrap replicates and
  columns index `timepoints`.

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
exposure. For the exposed model, the outcome is modeled in terms of time
since exposure among individuals who were exposed and who remained at
risk `tau` days after exposure. Both models adjust for the specified
covariates to help control for confounding. The second model also
flexibly adjusts for exposure time (by default, as a natural cubic
spline with 4 degrees of freedom) to capture time-varying background
risk. Predicted risks from both models are then marginalized over the
specified exposure-time and covariate distributions to obtain
G-computation style cumulative incidence estimates.

**Marginalizing weights.** When `weights_source = "observed"`, the
marginalizing weights are the empirical distributions of exposure times
and covariates among the exposed who remain at-risk `tau` days after
exposure. These weights are returned in the `nomatchfit` object under
`weights`. They can also be obtained prior to the call to `nomatch()` by
calling
[`get_observed_weights()`](https://ewu16.github.io/nomatch/reference/get_observed_weights.md).

**Confidence intervals.** Wald CIs are constructed on transformed
scales: \\\text{logit}\\ for cumulative incidence; \\\log{RR}\\ for risk
ratios/relative risk reduction, using bootstrap SEs. These are then
back-transformed to the original scale. No transformation is used for
risk differences.

**Parallelization.** Bootstraps can be parallelized on Unix via
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) by
providing `n_cores` argument.

## Examples

``` r
# Fit effectiveness model using simulated data

fit <- nomatch(
  data = simdata,
  outcome_time = "Y",
  outcome_status = "event",
  exposure = "V",
  exposure_time = "D_obs",
  covariates = c("x1", "x2"),
  timepoints = seq(30, 180, by = 30),
  immune_lag = 14,
  boot_reps = 5,
  n_cores = 1
)
#> Bootstrapping 5 samples...
#> Time difference of 1.800435 secs

# View basic results
fit$estimates
#> $cuminc_0
#>       estimate wald_lower wald_upper    wald_se wald_pval wald_n
#> 30  0.01163891 0.01007963 0.01343613 0.07419198        NA      5
#> 60  0.03790915 0.03331788 0.04310488 0.06829663        NA      5
#> 90  0.05857293 0.05262835 0.06514281 0.05781346        NA      5
#> 120 0.06691947 0.06032600 0.07417670 0.05651542        NA      5
#> 150 0.07548793 0.06861566 0.08298717 0.05247947        NA      5
#> 180 0.09576440 0.08852914 0.10352382 0.04414821        NA      5
#> 
#> $cuminc_1
#>        estimate wald_lower wald_upper    wald_se wald_pval wald_n
#> 30  0.006219935 0.00490254 0.00788853 0.12210980        NA      5
#> 60  0.022925618 0.02174377 0.02417012 0.02762124        NA      5
#> 90  0.035346676 0.03197190 0.03906329 0.05298012        NA      5
#> 120 0.044244807 0.04130257 0.04738628 0.03667769        NA      5
#> 150 0.055146409 0.04937548 0.06154816 0.05950454        NA      5
#> 180 0.079293545 0.07547382 0.08328917 0.02730204        NA      5
#> 
#> $risk_difference
#>         estimate   wald_lower   wald_upper     wald_se    wald_pval wald_n
#> 30  -0.005418978 -0.007191528 -0.003646428 0.000904379 2.073618e-09      5
#> 60  -0.014983530 -0.019314424 -0.010652635 0.002209681 1.194644e-11      5
#> 90  -0.023226258 -0.031540712 -0.014911803 0.004242146 4.372170e-08      5
#> 120 -0.022674663 -0.030836500 -0.014512826 0.004164279 5.179400e-08      5
#> 150 -0.020341520 -0.029851959 -0.010831081 0.004852354 2.763923e-05      5
#> 180 -0.016470857 -0.025719682 -0.007222032 0.004718875 4.822615e-04      5
#> 
#> $risk_ratio
#>      estimate wald_lower wald_upper    wald_se    wald_pval wald_n
#> 30  0.5344086  0.4211850  0.6780692 0.12147601 1.267008e-04      5
#> 60  0.6047516  0.5378709  0.6799484 0.05979658 3.846257e-11      5
#> 90  0.6034643  0.5060940  0.7195683 0.08977948 1.001848e-05      5
#> 120 0.6611649  0.5753649  0.7597597 0.07091909 1.772505e-06      5
#> 150 0.7305328  0.6292230  0.8481543 0.07616898 4.035446e-04      5
#> 180 0.8280065  0.7474679  0.9172230 0.05220997 9.867902e-04      5
#> 
#> $relative_risk_reduction
#>      estimate wald_lower wald_upper    wald_se    wald_pval wald_n
#> 30  0.4655914  0.3219308  0.5788150 0.12147601 1.267008e-04      5
#> 60  0.3952484  0.3200516  0.4621291 0.05979658 3.846257e-11      5
#> 90  0.3965357  0.2804317  0.4939060 0.08977948 1.001848e-05      5
#> 120 0.3388351  0.2402403  0.4246351 0.07091909 1.772505e-06      5
#> 150 0.2694672  0.1518457  0.3707770 0.07616898 4.035446e-04      5
#> 180 0.1719935  0.0827770  0.2525321 0.05220997 9.867902e-04      5
#> 
```
