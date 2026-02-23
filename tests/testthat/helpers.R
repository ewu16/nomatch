make_nomatch_fit <- function(...) {
    defaults <- list(
        data           = simdata,
        outcome_time   = "Y",
        outcome_status = "event",
        exposure       = "V",
        exposure_time  = "D_obs",
        covariates     = c("x1", "x2"),
        immune_lag     = 14,
        timepoints     = c(30, 60, 90)
    )
    do.call(nomatch, modifyList(defaults, list(...)))
}


make_matching_fit <- function(...){
    matched_data <- match_rolling_cohort(simdata, "Y", "V", "D_obs", 
                                         c("x1","x2"), "ID", seed = 123)$matched_data
    
    defaults <- list(
        matched_data   = matched_data,
        outcome_time   = "Y",
        outcome_status = "event",
        exposure       = "V",
        exposure_time  = "D_obs",
        immune_lag     = 14,
        timepoints     = c(30, 60, 90)
    )
    do.call(matching, modifyList(defaults, list(...)))
}