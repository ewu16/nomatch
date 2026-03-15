# Compute marginal cumulative incidence estimates for a matched cohort

This function is the main function for estimating cumulative incidence
in a matched dataset based on Kaplan Meier estimation. It creates the
analysis dataset and performs the analysis for the matched cohort.

## Usage

``` r
matching(
  matched_data,
  outcome_time,
  outcome_status,
  exposure,
  exposure_time,
  immune_lag = 0,
  timepoints = NULL,
  ci_type = c("wald", "percentile", "both"),
  boot_reps = 0,
  alpha = 0.05,
  keep_models = TRUE,
  keep_boot_samples = TRUE,
  seed = NULL
)
```

## Arguments

- matched_data:

  A data frame for the matched cohort created using
  [`match_rolling_cohort()`](https://ewu16.github.io/nomatch/reference/match_rolling_cohort.md).

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

  Logical; return the fitted Kaplan Meier object used to compute
  cumulative incidences? Default: `TRUE`

- keep_boot_samples:

  Logical; return bootstrap samples? Default: `TRUE`. Must be set to
  `TRUE` if user plans to use
  [`add_simultaneous_ci()`](https://ewu16.github.io/nomatch/reference/add_simultaneous_ci.md)
  to obtain simultaneous confidence intervals.

- seed:

  Integer seed for reproducible bootstrap results. Default: `NULL`.

## Value

An object of class `nomatchfit` containing:

- matched_adata:

  The analytic dataset used for the matching analysis, often a subset of
  the original matched dataset when `immune_lag > 0`. The data frame
  consists of the original variables in `matched_data` plus:

  `match_<outcome_status>`

  :   Event indicator for the outcome in the matched analysis

  `match_<outcome_time>`

  :   Follow-up time for outcome of interest in the matched analysis,
      relative to original time-origin

  `match_T`

  :   Follow-up time for outcome of interest in the matched analysis,
      relative to time of matching

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

- models:

  Fitted Kaplan Meier object

- n_success_boot:

  Integer vector indicating the number of successful bootstrap
  replications per timepoint.

- boot_samples:

  (If `keep_boot_samples = TRUE`) Named list of bootstrap draws (stored
  as matrices) for each term. Rows index bootstrap replicates and
  columns index `timepoints`.

## Details

**Matching analysis**: To create the matched analysis dataset, both
individuals in the matched pair are censored when the control receives
the exposure, and matched pairs in which either individual experiences
the event or censoring with `immune_lag` days of matching is excluded.
Cumulative incidence in the matched analysis dataset is estimated using
Kaplan Meier stratified only on exposure group.

**Confidence intervals.** Wald and percentile confidence intervals are
constructed for cumulative incidence and effectiveness parameters at
each timepoint. The Wald pointwise confidence intervals are constructed
on transformed scales: \\\text{logit}\\ for cumulative incidence;
\\\log{RR}\\ for risk ratios, and \\\log{1 - RR}\\ for relative risk
reduction, using bootstrap standard errors. These confidence intervals
are then back-transformed to the original scale. Identity transformation
is used for risk differences. To obtain simultaneous confidence
intervals, use
[`add_simultaneous_ci()`](https://ewu16.github.io/nomatch/reference/add_simultaneous_ci.md)
after saving the original fit.

**Parallelization.** Bootstraps can be parallelized using the `future`
framework. Set a parallel plan before calling `matching()`: e.g.

    future::plan(future::multisession, workers = 4)
    fit <- matching(..., boot_reps = 1000, seed = 42)
    future::plan(future::sequential)  # reset when done

If no plan is set, bootstraps run sequentially. `multisession` works on
all operating systems and is recommended for most users. See the [future
package documentation](https://future.futureverse.org) for additional
plans and details on setup.

## Examples

``` r
# Perform matching analysis on simulated data

# Create matched data
matched_cohort <- match_rolling_cohort(
  data = simdata,
  outcome_time =  "Y",
  exposure = "V",
  exposure_time = "D_obs",
  matching_vars = c("x1", "x2"),
  id_name = "ID",
  seed = 5678
)

# Analyze matched data
fit_matching <- matching(
  matched_data = matched_cohort$matched_data,
  outcome_time = "Y",
  outcome_status = "event",
  exposure = "V",
  exposure_time = "D_obs",
  timepoints = seq(30, 90, by = 30),
  immune_lag = 14,
  boot_reps = 5,
  seed = 123
)
#> Bootstrapping 5 samples...
#> Bootstrap completed in 1.02 secs

# View basic results
fit_matching
#> 
#>  Risk Ratio Estimates 
#> ================================================== 
#> Call: matching(matched_data = matched_cohort$matched_data, outcome_time = "Y", 
#>     outcome_status = "event", exposure = "V", exposure_time = "D_obs", 
#>     immune_lag = 14, timepoints = seq(30, 90, by = 30), boot_reps = 5, 
#>     seed = 123) 
#> 
#> Result:
#>   Timepoint Estimate 95% Wald CI: Lower 95% Wald CI: Upper Wald p-value
#> 1        30    0.497              0.196              1.260     0.140809
#> 2        60    0.517              0.443              0.604     0.000000
#> 3        90    0.611              0.457              0.818     0.000934
#> 
#> Use summary() for more details
#> Use plot() to visualize results
```
