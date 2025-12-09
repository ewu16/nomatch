# Internal transformations

These define how estimates are transformed before applying the Wald CI
(e.g., taking log or logit). These define how to back-transform point
estimates Functions take arguments `(eta, sd, z)` and return values on
the original scale.

## Usage

``` r
.transformations
```

## Format

An object of class `list` of length 4.
