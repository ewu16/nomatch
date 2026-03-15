# Add simultaneous confidence intervals to `nomatchfit` object

Computes simultaneous confidence intervals, which maintain the specified
coverage level for the parameter of interest across all evaluation
timepoints jointly. This is useful for making inferences about the
entire cumulative incidence or effectiveness curve. The simultaneous
confidence intervals returned contrast with the pointwise confidence
intervals automatically returned with a `nomatchfit` object when
`boot_reps > 0`.

## Usage

``` r
add_simultaneous_ci(object, seed = NULL)
```

## Arguments

- object:

  An object of class `nomatchfit` created by
  [`nomatch()`](https://ewu16.github.io/nomatch/reference/nomatch.md) or
  [`matching()`](https://ewu16.github.io/nomatch/reference/matching.md).
  Must

  - contain evaluations at at least 2 timepoints
    (`length(object$timepoints) > = 2`),

  - contain more than one bootstrap replicate (`boot_reps` \> 1)

  - contain non-null bootstrap samples (`keep_boot_samples = TRUE` when
    fitting).

- seed:

  Integer seed to ensure reproducible critical values for simultaneous
  confidence intervals (relied on random number generation). Default is
  `NULL` (no seed set).

## Value

The original `nomatchfit` object with the following modifications:

- estimates:

  The matrix for each term in the `estimates` list gets additional
  columns describing the simultaneous confidence interval bounds and
  construction: `simul_lower`, `simul_upper`, `simul_n`

- simul_z_star:

  Named vector containing the critical values used for each term

- simul_excluded_timepoints:

  Named list containing a vector of timepoints excluded from the
  simultaneous bands for each term due to insufficient bootstrap samples

## Details

Critical values are computed using the bootstrap covariance structure
across timepoints. If any timepoint has more than 5% missing bootstrap
samples, it is excluded from the simultaneous band and a warning is
issued.

## Examples

``` r
# Fit model with bootstrap samples
fit <- nomatch(
  data = simdata,
  outcome_time = "Y",
  outcome_status = "event",
  exposure = "V",
  exposure_time = "D_obs",
  covariates = c("x1", "x2"),
  timepoints = seq(30, 180, by = 30),
  immune_lag = 14,
  boot_reps = 100,         # must be > 1
  keep_boot_samples = TRUE # must be set to TRUE (default)
)
#> Bootstrapping 100 samples...
#> Bootstrap completed in 32.65 secs

# Add simultaneous CIs
fit_simul <- add_simultaneous_ci(fit, seed = 123)

# Look at results
fit_simul
#> 
#>  Risk Ratio Estimates 
#> ================================================== 
#> Call: nomatch(data = simdata, outcome_time = "Y", outcome_status = "event", 
#>     exposure = "V", exposure_time = "D_obs", covariates = c("x1", 
#>         "x2"), immune_lag = 14, timepoints = seq(30, 180, by = 30), 
#>     boot_reps = 100, keep_boot_samples = TRUE) 
#> 
#> Result:
#>   Timepoint Estimate 95% Wald CI: Lower 95% Wald CI: Upper Wald p-value
#> 1        30    0.534              0.346              0.826     4.81e-03
#> 2        60    0.605              0.490              0.746     2.74e-06
#> 3        90    0.603              0.511              0.712     2.28e-09
#> 4       120    0.661              0.561              0.779     7.28e-07
#> 5       150    0.731              0.629              0.849     4.28e-05
#> 6       180    0.828              0.715              0.959     1.18e-02
#>   95% Simul CI: Lower 95% Simul CI: Upper
#> 1               0.308               0.928
#> 2               0.463               0.789
#> 3               0.489               0.744
#> 4               0.537               0.813
#> 5               0.604               0.884
#> 6               0.687               0.997
#> 
#> Use summary() for more details
#> Use plot() to visualize results

# Visualize
plot(fit_simul, ci_type = "simul")
```
