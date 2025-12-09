# Convert estimates from a fitted object into a tidy data frame

A convenience function for extracting estimates from a `nomatchfit`
object (or its `estimates` component) into a tidy data frame for
plotting and summary tables.

## Usage

``` r
estimates_to_df(x, collapse = TRUE)
```

## Arguments

- x:

  A fitted `nomatchfit` object or its `estimates` component (a named
  list of matrices).

- collapse:

  Logical; if `TRUE` (default), returns a single long-format data frame.
  If `FALSE`, returns a list of data frames (one per term).

## Value

Either:

- A long-format data frame (if `collapse = TRUE`), or

- A list of term-specific data frames (if `collapse = FALSE`).

The data frame(s) have columns:

- `t0`: evaluation timepoint

- `term`: estimated quantity

- other columns from the original matrices (e.g. estimate, confidence
  intervals)
