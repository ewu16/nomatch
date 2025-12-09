# Simulated dataset

Minimal dataset mimicking data from an observational vaccine study. Used
to apply/test methods.

## Usage

``` r
simdata
```

## Format

### `simdata`

A data frame with 10,000 rows and 7 columns:

- ID:

  Subject identifier

- x1,x2:

  Baseline covariates to adjust for

- V:

  Indicator of ever being vaccinated while under follow-up

- D_obs:

  Day of observed vaccination (relative to study start)

- Y:

  Time to event in days (relative to study start)

- event:

  Event indicator, event = 0 if censored

## Source

Created from code
