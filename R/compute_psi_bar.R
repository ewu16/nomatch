#' Compute overall cumulative incidences at multiple timepoints
#'
#' @description
#' Wrapper that computes marginal cumulative incidences at multiple
#' timepoints. Calls `compute_psi_bar_t0()` internally for each timepoint.
#' Models are fitted once before calling this function to allow efficient evaluation
#' at multiple timepoints without refitting.
#'
#' @inheritParams compute_psi_bar_t0
#' @inheritParams nomatch
#'
#' @return A matrix of estimates where the columns are the terms `cuminc_0` and `cuminc_1`,
#' and the rows are the time points of interest.
#' @details
#'  This function expects models already fitted via [fit_model_0()] and
#' [fit_model_1()]. It performs G-computation by predicting and marginalizing over
#' conditional risks.
#'
#' @keywords internal
compute_psi_bar_times <- function(fit_0, fit_1, exposure_time, timepoints, tau, newdata, gp_list){
    out <- sapply(timepoints, \(t){
        compute_psi_bar_t0(fit_0, fit_1, exposure_time, t0 = t, tau, newdata, gp_list)
    })
    colnames(out) <- timepoints
    t(out)
}

#' Compute overall cumulative incidences at a single timepoint. Internally calls
#'  [compute_psi_dx_t0()] to predict time and covariate-specific cumulative incidences 
#'  and [marginalize_psi_dx_t0] to marginalize these cumulative incidences over 
#'  the selected weight functions for time and covariates. 
#'
#' @inheritParams compute_psi_dx_t0
#' @inheritParams  marginalize_psi_dx_t0
#' @return Named numeric vector with `cuminc_0`, `cuminc_1`
#' @keywords internal
compute_psi_bar_t0 <- function(fit_0, fit_1, exposure_time, t0, tau, newdata, gp_list){
    if (t0 <= tau) {
        # By definition under the immune-lag/principal-stratum estimand
        return(c(cuminc_0 = 0, cuminc_1 = 0))
    }
    
    # Otherwise compute psi_bar_estimate
    # Step 1: Get conditional risks for all (d, x)
    psi_dx <- compute_psi_dx_t0(fit_0, fit_1, exposure_time, t0, tau, newdata)

    # Step 2: Marginalize to population level
    psi_bar <- marginalize_psi_dx_t0(psi_dx, gp_list)

    psi_bar
}
