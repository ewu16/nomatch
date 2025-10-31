#' Validate inputs to main functions
#'
#' @description
#' Internal validation function to check that all inputs to `nomatch()/matching()` are
#' properly formatted and logically consistent.
#'
#' @inheritParams nomatch
#'
#' @return Invisibly returns `NULL` if all checks pass. Throws an error with
#'   descriptive message if any validation fails.
#'
#' @keywords internal
#' @noRd
validate_nomatch_inputs <- function(data, outcome_time, outcome_status, exposure, exposure_time, covariates, immune_lag, timepoints) {
    validate_data(data)
    required_args <- list(
        outcome_time = outcome_time,
        outcome_status = outcome_status,
        exposure = exposure,
        exposure_time = exposure_time,
        covariates = covariates,
        immune_lag = immune_lag,
        timepoints = timepoints)

    validate_args_not_null(required_args)

    cols <- c(outcome_time, outcome_status, exposure, exposure_time, covariates)
    validate_cols_exist(data, cols)

    validate_outcome_time(data, outcome_time)
    validate_outcome_status(data, outcome_status)
    validate_exposure_args(data, exposure, exposure_time)
    validate_covariates(data, covariates)
    validate_time_args(immune_lag, timepoints)
    invisible(NULL)
}

validate_match_rolling_cohort_inputs <- function(data, outcome_time, exposure, exposure_time, covariates) {
    validate_data(data)
    required_args <- list(
        outcome_time = outcome_time,
        exposure = exposure,
        exposure_time = exposure_time,
        covariates = covariates)

    validate_args_not_null(required_args)

    cols <- c(outcome_time, exposure, exposure_time, covariates)
    validate_cols_exist(data, cols)

    validate_outcome_time(data, outcome_time)
    validate_exposure_args(data, exposure, exposure_time)
    validate_covariates(data, covariates)
    invisible(NULL)
}

validate_matching_inputs <- function(data, outcome_time, outcome_status, exposure, exposure_time,  immune_lag, timepoints) {
    validate_data(data)
    required_args <- list(
        outcome_time = outcome_time,
        outcome_status = outcome_status,
        exposure = exposure,
        exposure_time = exposure_time,
        immune_lag = immune_lag,
        timepoints = timepoints)

    validate_args_not_null(required_args)

    cols <- c(outcome_time, outcome_status, exposure, exposure_time)
    validate_cols_exist(data, cols)

    validate_outcome_time(data, outcome_time)
    validate_outcome_status(data, outcome_status)
    validate_exposure_args(data, exposure, exposure_time)
    validate_time_args(immune_lag, timepoints)
    invisible(NULL)
}



