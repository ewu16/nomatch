# Get the marginal cumulative incidence in the treated and untreated groups based on Kaplan Meier estimation.

Get the marginal cumulative incidence in the treated and untreated
groups based on Kaplan Meier estimation.

## Usage

``` r
compute_km_ve(
  adata,
  adata_outcome_name,
  adata_event_name,
  adata_trt_name,
  timepoints,
  keep_models = TRUE
)
```

## Arguments

- adata:

  A data frame that represents the analysis data set of a clinical
  trial.

- adata_outcome_name:

  Character string specifying the time to event variable in `adata`. The
  time should be the time to event from vaccination/matched index date

- adata_event_name:

  Character string specifying the event variable in `adata`

- adata_trt_name:

  Character string specifying the treatment variable in `adata`

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

- keep_models:

  Logical; return the two fitted hazard models used to compute
  cumulative incidences? Default: `TRUE`.

## Value

A list containing the following:

- estimates:

  A matrix of estimates. The columns of the matrix are the cumulative
  incidence/VE terms and the rows are the requested time points for
  evaluation.

- models:

  If `keep_models = TRUE`, the models used to compute risk and VE
