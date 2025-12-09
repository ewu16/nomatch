# Summary method for effectiveness fits

Summarizes how effectiveness estimates were obtained

## Usage

``` r
# S3 method for class 'nomatchfit'
summary(object, digits = 4, show_models = FALSE, ...)
```

## Arguments

- object:

  An object of class `nomatchfit` created by
  [`nomatch()`](https://ewu16.github.io/nomatch/reference/nomatch.md) or
  [`matching()`](https://ewu16.github.io/nomatch/reference/matching.md).

- digits:

  Integer indicating the number of decimal places to display. Default is
  4.

- show_models:

  Logical; print model details? Default is FALSE.

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns the input object `object`. Called primarily for its
side effect of printing a detailed summary to the console.
