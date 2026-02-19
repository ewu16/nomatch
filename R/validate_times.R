#' Internal functions for helping to validate input
#'
#' @keywords internal
#' @noRd




# Parameter checks --------------------------------------------------------
validate_immune_lag <- function(immune_lag){
    if (!is.numeric(immune_lag) || length(immune_lag) != 1L || is.na(immune_lag))
        stop("`immune_lag` must be a single non-missing numeric value.", call. = FALSE)
    if (immune_lag < 0) stop("`immune_lag` must be >= 0.", call. = FALSE)
    if (immune_lag > 28) warning("immune_lag = ", immune_lag, " days is unusually long. Please verify this is intentional.")
    invisible(NULL)
}

resolve_timepoints <- function(data, outcome_time, outcome_status, exposure, timepoints, immune_lag){
    if(is.null(timepoints)){
        max_follow_up <- max(data[[outcome_time]][data[[exposure]] == 1 & data[[outcome_status]] == 1])
        if(max_follow_up <= immune_lag){
            stop("Maximum follow-up time among exposed events (", max_follow_up, ") is <= `immune_lag` (", immune_lag, "). ",
                 "The default `timepoints` is therefore not valid. Please choose a smaller `immune_lag` or manually specify
                 `timepoints` so that the maximum timepoint > `immune_lag` .")
        }
        timepoints <- seq_len(max_follow_up)
    }else{
        if(anyDuplicated(timepoints) > 0){
            warning("Non-unique timepoints were provided. Using only unique timepoint values.")
        }
        timepoints <- sort(unique(timepoints))
    }
    timepoints
}


validate_timepoints <- function(timepoints, immune_lag){
    if (!is.numeric(timepoints) || length(timepoints) < 1 || anyNA(timepoints))
        stop("`timepoints` must be a numeric vector with length > 0 and no missing values.", call. = FALSE)
    if (max(timepoints) <= immune_lag)
        stop("The maximum timepoint must be > `immune_lag` (", immune_lag, ") because the exposed hazard is fit using data censored at the maximum timepoint.", call. = FALSE)
    invisible(NULL)
}



