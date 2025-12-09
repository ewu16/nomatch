# Marginalize conditional cumulative incidences

Averages day- and covariate-specific cumulative incidences from
[`compute_psi_dx_t0()`](https://ewu16.github.io/nomatch/reference/compute_psi_dx_t0.md)
to produce overall cumulative incidences under each exposure. By
default, performs a simple average. If weights are provided for in
`gp_list`, performs a weighted average.

## Usage

``` r
marginalize_psi_dx_t0(psi_dx, gp_list = NULL, show_intermediates = FALSE)
```

## Arguments

- psi_dx:

  A data frame with conditional risk predictions containing:

  - `<covariates>`: all covariate columns

  - `<exposure_time>`: vaccination day column

  - `psi_0_dx`: cumulative incidence under no vaccine,\\\psi_0(t_0; d,
    x)\\

  - `psi_1_dx`: cumulative incidence under vaccine,\\\psi_1(t_0; d, x)\\

  Typically output from
  [`compute_psi_dx_t0()`](https://ewu16.github.io/nomatch/reference/compute_psi_dx_t0.md).

- gp_list:

  A list with two data frames:

  - **g_weights** Data frame of covariate-conditional exposure-day
    probabilities \\g(d \mid x)\\. Must include:

    - `group_id`: covariate group identifier

    - `<exposure_time>`: exposure time variable

    - `prob_g`: probability of exposure time given the covariates

    - all variables in `covariates`

  - **p_weights** Data frame of covariate probabilities \\p(x)\\. Must
    include:\\

    - `group_id`: covariate group identifier

    - `prob_p`: marginal probability of each covariate group

    - all variables in `covariates`

  Default is `NULL` in which case each row of `psi_dx` gets equal
  weight.

- show_intermediates:

  Logical that only applies when `gp_list` is not `NULL`.

  - When `FALSE` (default), the function performs marginalization in a
    single step and returns only the overall cumulative incidences
    \\\bar{\psi}\_0(t_0)\\ and \\\bar{\psi}\_1(t_0)\\.

  - When `TRUE`, performs the two-step marginalization (first over days
    \\d\\, then covariate groups \\x\\) and returns intermediate
    group-level results \\\psi_v(t_0; x)\\.

## Value

When `show_intermediates = FALSE`, a named numeric vector with two
elements:

- `cuminc_0`: Overall cumulative incidence under no exposure

- `cuminc_1`: Overall cumulative incidence under exposure

When `show_intermediates = TRUE`, a list containing both overall results
(as a vector) and group-specific results (as a data frame).

## See also

- [`compute_psi_dx_t0()`](https://ewu16.github.io/nomatch/reference/compute_psi_dx_t0.md)
  for generating `psi_dx` input

- `canonicalize_weights()` for creating `gp_list`
