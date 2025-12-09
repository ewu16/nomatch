# Print method for effectiveness fits

Prints a concise summary of effectiveness estimates from a fitted model.

## Usage

``` r
# S3 method for class 'nomatchfit'
print(
  x,
  digits = 3,
  effect = c("risk_ratio", "risk_difference", "relative_risk_reduction"),
  ...
)
```

## Arguments

- x:

  An object of class `nomatchfit` created by
  [`nomatch()`](https://ewu16.github.io/nomatch/reference/nomatch.md) or
  [`matching()`](https://ewu16.github.io/nomatch/reference/matching.md).

- digits:

  Integer indicating the number of decimal places to display. Default is
  3.

- effect:

  The effect measure to output. Either `"risk_ratio"`(default),
  `"relative_risk_reduction"` or `"risk_difference"`.

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns the input object `x`. Called primarily for its side
effect of printing a summary to the console.
