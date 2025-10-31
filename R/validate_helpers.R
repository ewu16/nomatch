#' Internal functions for helping to validate input
#'
#' @keywords internal
#' @noRd
NULL

validate_data <- function(data){
    if(!is.data.frame(data)) stop("'data' must be a data.frame, not ", class(data)[1])
    if(nrow(data) == 0)  stop("'data' has no rows")
    invisible(NULL)
}

validate_args_not_null <- function(args) {
    if(is.null(names(args)) || any(names(args) == "")) stop("All elements of args must be named.")

    missing <- names(args)[vapply(args, is.null, logical(1))]
    if (length(missing)) {
        stop("NULL values are not allowed for the following arguments: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    invisible(NULL)
}


validate_cols_exist <- function(data, cols){
    missing_cols <- setdiff(cols, names(data))
    if(length(missing_cols) > 0) {
        stop("Missing required column(s) in data: ", paste(missing_cols, collapse = ", "))
    }
    invisible(NULL)
}

validate_outcome_time <- function(data, outcome_time){
    if (anyNA(data[[outcome_time]])) {
        stop("Missing values in `outcome_time` are not supported. Please remove these observations.")
    }
    if (any(data[[outcome_time]] < 0)) {
        stop("Outcome time variable '", outcome_time, "' must be non-negative.")
    }
    invisible(NULL)
}


validate_outcome_status <- function(data, outcome_status){
    if (anyNA(data[[outcome_status]])) {
        stop("Missing values in `outcome_status` are not supported. Please remove these observations.")
    }
    if (!all(data[[outcome_status]] %in% c(0,1))) {
        stop("Outcome status variable '", outcome_status, "' must be coded 0/1.")
    }
    invisible(NULL)
}

validate_exposure_args <- function(data, exposure, exposure_time){
    # Check exposure
    if (anyNA(data[[exposure]])) {
        stop("Missing values in <exposure> are not supported. Please remove these observations.")
    }
    if (!all(data[[exposure]] %in% c(0,1))) {
        stop("Exposure variable '", exposure, "' must be coded 0/1.")
    }

    # Check exposure_time
    et <- data[[exposure_time]]
    if (any(et[!is.na(et)] < 0)) {
        stop("Exposure time variable '", exposure_time, "' must be non-negative.")
    }
    exposed <- data[[exposure]] == 1
    if (any(is.na(et[exposed])))       stop("All exposed individuals must have non-missing ", exposure_time, ".")
    if (any(!is.na(et[!exposed])))     stop("All unexposed individuals must have NA ", exposure_time, ".")

    invisible(NULL)
}


validate_covariates <- function(data, covariates){
    if(!is.null(covariates) && sum(is.na(data[, covariates])) > 0){
        stop("Missing values in covariates '",
             paste(covariates, collapse = ","), "' are not supported.",
             "Please remove these observations")
    }
    invisible(NULL)
}

validate_time_args <- function(immune_lag, timepoints){
    if (!is.numeric(immune_lag) || length(immune_lag) != 1L || is.na(immune_lag))
        stop("`immune_lag` must be a single non-missing numeric value.", call. = FALSE)
    if (immune_lag < 0) stop("`immune_lag` must be >= 0.", call. = FALSE)
    if (!is.numeric(timepoints) || length(timepoints) < 1 || anyNA(timepoints))
        stop("`timepoints` must be a numeric vector with length > 0 and no missing values.", call. = FALSE)
    if (any(timepoints <= immune_lag))
        stop("All `timepoints` must be > `immune_lag` (", immune_lag, ").", call. = FALSE)
    if (immune_lag > 28) warning("immune_lag = ", immune_lag, " days is unusually long. Please verify this is intentional.")
    invisible(NULL)
}


check_reserved_vars <- function(vars, reserved_vars, vars_label) {
    conflicts <- intersect(reserved_vars, vars)
    if (length(conflicts) > 0) {
        stop(
            paste(vars_label, " contain reserved variable names: "),
            paste(conflicts, collapse = ", "),
            "\nPlease rename these variables before fitting the model.",
            call. = FALSE
        )
    }
    invisible(NULL)
}

