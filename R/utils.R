# Format estimates --------------------------------------------------------

#' Convert estimates from a fitted object into a tidy data frame
#'
#' @description
#' A convenience function for extracting estimates from a
#' [`vefit`] object (or its `estimates` component) into a tidy data frame
#' for plotting and summary tables.
#'
#'  @param x A fitted [`vefit`] object or its `estimates` component (a named list
#'   of matrices).
#' @param collapse Logical; if `TRUE` (default), returns a single long-format
#'   data frame. If `FALSE`, returns a list of data frames (one per term).
#'
#' @return Either:
#'   - A long-format data frame (if `collapse = TRUE`), or
#'   - A list of term-specific data frames (if `collapse = FALSE`).
#'   The data frame(s) have columns:
#'    \itemize{
#'     \item `t0`: evaluation timepoint
#'     \item `term`: estimated quantity
#'     \item other columns from the original matrices (e.g. estimate, confidence intervals)
#'   }
#'
#' @export
#'
estimates_to_df <- function(x, collapse = TRUE){
    # detect if full vefit object
    if (inherits(x, "vefit")) {
        estimates <- x$estimates
    } else {
        estimates <- x
    }

    # validate
    if (!is.list(estimates) || !all(vapply(estimates, is.matrix, logical(1)))) {
        stop("Input must be a 'vefit' object or a list of matrices.", call. = FALSE)
    }

    df_list <- lapply(seq_along(estimates), \(i){
        df <-data.frame(estimates[[i]])
        df$term <- names(estimates)[i]
        df$t0 <- as.numeric(rownames(df))
        df <- df[, c("t0", "term", colnames(estimates[[i]]))]
    })

    if(collapse){
        df <- do.call("rbind", df_list)
        df <-  df[order(df$t0),]
        rownames(df) <- NULL
        df
    }else{
        df_list
    }
}


#' Helper to save warning calls when a function is called
#'
#' @noRd
capture_warnings <- function(expr){
    w <- NULL
    result <- withCallingHandlers(
        expr,
        warning = function(m) {
            w <<- c(w, conditionMessage(m))
            invokeRestart("muffleWarning")
        }
    )
    list(result = result, warnings = w)
}


get_basic_descriptives_nomatch <- function(data, outcome_time, outcome_status, exposure, exposure_time, tau){
    n <- nrow(data)
    n_events <- sum(data[[outcome_status]])

    exposed <- data[[exposure]] == 1 & (data[[outcome_time]] - data[[exposure_time]] > 0)
    n_exposed <- sum(exposed)

    exposed_at_tau <- data[[exposure]] == 1 & (data[[outcome_time]] - data[[exposure_time]] > tau)
    n_exposed_at_tau <- sum(exposed_at_tau)

    exposure_times_at_tau <- data[[exposure_time]][exposed_at_tau]

    list(n = n,
         n_events = n_events,
         n_exposed = n_exposed,
         n_exposed_at_tau = n_exposed_at_tau,
         exposure_times_at_tau = exposure_times_at_tau)
}

get_basic_descriptives_matching <- function(matched_data, matched_adata,  outcome_status, exposure){
    n_matched <- nrow(matched_data)
    n_matched_analysis <- nrow(matched_adata)

    n_events <- sum(matched_adata[[paste0("match_", outcome_status)]])

    exposure_times <- matched_adata$match_index_time[matched_adata[[paste0("match_", exposure)]] == 1]

    list(n_matched = n_matched,
         n_matched_analysis = n_matched_analysis,
         n_events = n_events,
         exposure_times = exposure_times)
}



check_pt_estimates <- function(pt_estimates){

    both_zero <- which(pt_estimates[, "cuminc_0"] == 0 & pt_estimates[, "cuminc_1"] == 0 &
                           is.nan(pt_estimates[, "risk_ratio"]))
    if(length(both_zero) > 0){
        warning("Cumulative incidences in unexposed and exposed groups are both zero, ",
                "resulting in undefined (NaN) \n risk ratios/VE estimates at the following timepoints: ",
                paste(names(both_zero), collapse = ", "))
    }

    unexposed_zero <- which(pt_estimates[, "cuminc_0"] == 0 & pt_estimates[, "cuminc_1"] > 0 &
                                is.nan(pt_estimates[, "risk_ratio"]))

    if(length(unexposed_zero) > 0){
        warning("Cumulative incidence in unexposed group is zero, ",
                "resulting in undefined (NaN) \n risk ratios/VE estimates at the following timepoints: ",
                paste(names(both_zero), collapse = ", "))
    }
    invisible(NULL)
}


