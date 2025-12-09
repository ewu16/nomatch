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
  immune_lag,
  timepoints,
  ci_type = c("wald", "percentile", "both"),
  boot_reps = 0,
  alpha = 0.05,
  keep_models = TRUE,
  keep_boot_samples = TRUE,
  n_cores = 1
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
  used for `outcome_time`. Must be a non-missing numeric value exposed
  individuals and must be set to `NA` for unexposed individuals.

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

- models:

  Fitted Kaplan Meier

- n_success_boot:

  Integer vector indicating the number of successful bootstrap
  replications per timepoint.

- boot_samples:

  (If `keep_boot_samples = TRUE`) Named list of bootstrap draws (stored
  as matrices) for each term. Rows index bootstrap replicates and
  columns index `timepoints`.
