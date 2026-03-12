# Helper for constructing bootstrapped confidence intervals

This function computes Wald and percentile bootstrapped confidence
intervals.

## Usage

``` r
estimate_bootstrap_ci(
  one_boot_function,
  one_boot_args,
  ci_type,
  boot_reps,
  pt_est = NULL,
  alpha = 0.05,
  keep_boot_samples = TRUE,
  seed = NULL
)
```

## Arguments

- one_boot_function:

  Function that computes one bootstrap iteration and returns the
  bootstrap estimates

- one_boot_args:

  List of arguments to pass to `one_boot_function`

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

- pt_est:

  A matrix of point estimates, with columns corresponding to cumulative
  incidence and effect measures, and rows representing evaluation
  timepoints. This argument must be provided when `ci_type = "wald"`.
  Ignored when `ci_type = "percentile"`.

- alpha:

  Significance level for confidence intervals (Confidence level =
  100\*(1-alpha)%). Default: `0.05`.

- keep_boot_samples:

  Logical; return bootstrap samples? Default: `TRUE`. Must be set to
  `TRUE` if user plans to use
  [`add_simultaneous_ci()`](https://ewu16.github.io/nomatch/reference/add_simultaneous_ci.md)
  to obtain simultaneous confidence intervals.

- seed:

  Integer seed for reproducible bootstrap results. Default: `NULL`.

## Value

A list containing the following:

- ci_estimates:

  A list of matrices, one for each term, containing the lower and upper
  confidence interval bounds

- n_success_boot:

  A numeric vector of the number of successful bootstrap samples for
  each time point.(Success bootstrap samples are those that result in
  non-missing valid point estimates.)

- boot_samples:

  If `keep_boot_samples = TRUE`, a list of matrices, one for each term,
  containing the bootstrap estimates. Rows are the bootstrap iterations
  and columns are the time points.

When `boot_reps = 0`, returns `NULL`.
