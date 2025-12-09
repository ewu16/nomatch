# Calculate critical value for simultaneous confidence bands

Computes the critical value for simultaneous confidence intervals by
simulating from the bootstrap covariance structure (10,000 draws from
multivariate normal).

## Usage

``` r
get_z_star(mat, alpha, seed = NULL)
```

## Arguments

- mat:

  Matrix of bootstrap estimates (rows = iterations, columns =
  timepoints)

- alpha:

  Significance level (e.g., 0.05 for 95% confidence)

- seed:

  Integer seed for reproducibility. Default is `NULL`

## Value

List containing:

- z_star:

  Critical value (scalar)

- excluded_timepoints:

  Timepoints excluded due to \>5% missing bootstrap estimates
