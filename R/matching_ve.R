#' Compute a matching-based estimator of VE with confidence intervals
#'
#' @description This function is the main function for computing a matching-based estimator
#' based on Kaplan Meier estimation.
#'
#' @inheritParams nomatchVE
#' @inheritParams clean_matched_data
#' @param matched_data A data frame for the matched cohort

#' @return A list containing the following:
#' \describe{
#'  \item{estimates}{A list of matrices of the estimates at each timepoint. Rows of
#'  each matrix are the terms "cuminc_0", "cuminc_1", "vaccine_effectiveness". Columns of each matrix
#'  gives the point estimate and confidence intervals at the specified time point.}
#'  \item{eval_times}{The timepoints at which VE was evaluated}
#'  \item{n_success_boot}{A numeric vector of the number of successful bootstrap samples for each time point.(Success bootstrap samples are
#'  those that result in non-missing valid point estimates.)}
#'  \item{boot_samples}{If `keep_boot_samples = TRUE`, a list of matrices for each term that contain the bootstrap estimates where the rows are the bootstrap iterations and
#'  the columns are the time points.}
#' }
#' @export
#'
matching_ve <- function(matched_data,
                        outcome_time,
                        outcome_status,
                        exposure,
                        exposure_time,
                        immune_lag,
                        eval_times,
                        effect = c("vaccine_effectiveness", "risk_ratio"),
                        ci_type = c("wald", "percentile", "both"),
                        boot_reps = 0,
                        alpha = 0.05,
                        keep_models = TRUE,
                        keep_boot_samples = TRUE,
                        n_cores = 1){

    call <- match.call()
    effect = match.arg(effect)
    ci_type <- match.arg(ci_type)
    tau <- immune_lag


    # Check data/inputs
    validate_ve_inputs(
        data = matched_data,
        outcome_time = outcome_time,
        outcome_status = outcome_status,
        exposure = exposure,
        exposure_time = exposure_time,
        covariates = NULL,
        tau = tau,
        eval_times = eval_times
    )

    # --------------------------------------------------------------------------
    # 1 - Get original estimate
    # --------------------------------------------------------------------------
    estimation_args <- list(matched_data = matched_data,
                            outcome_time = outcome_time,
                            outcome_status = outcome_status,
                            exposure = exposure,
                            exposure_time = exposure_time,
                            tau = tau,
                            eval_times = eval_times)

    original <- do.call(get_one_matching_ve, estimation_args)

    descrip <- get_basic_descriptives_matching(
        matched_data,
        original$matched_adata,
        outcome_status = outcome_status,
        exposure = exposure)


    # --------------------------------------------------------------------------
    # 2 - Get bootstrap CI
    # --------------------------------------------------------------------------

    # Helper returns NULL if boot_reps = 0
    boot <- estimate_bootstrap_ci(
        one_boot_function  = one_boot_matching,
        one_boot_args      = estimation_args,
        ci_type            = ci_type,
        boot_reps          = boot_reps,
        pt_est             = original$pt_estimates,
        alpha              = alpha,
        keep_boot_samples  = keep_boot_samples,
        n_cores            = n_cores
    )
    ci_est <- boot$ci_estimates
    boot_samples <- boot$boot_samples

    # --------------------------------------------------------------------------
    # 3 - Add p-values and format
    # --------------------------------------------------------------------------
    terms_keep <- names(ci_est)

    est <- lapply(names(ci_est), \(term){
        #Add p-values
        wald <- ci_type %in% c("wald", "both")
        percentile <- ci_type %in% c("percentile", "both")

        wald_pval <-       if(wald) compute_wald_pval(term, ci_est[[term]]) else NULL
        percentile_pval <- if(percentile) compute_percentile_pval(term, boot_samples[[term]]) else NULL

        x <- cbind(ci_est[[term]],
                   wald_pval = wald_pval,
                   percentile_pval = percentile_pval)

        #Format column order
        col_order <- c("estimate",
                       paste0("wald_", c("lower", "upper", "se", "pval", "n")),
                       paste0("percentile_", c("lower", "upper", "pval", "n")))

        x[, intersect(col_order, colnames(x)), drop = FALSE]

    })
    names(est) <- terms_keep


    # --------------------------------------------------------------------------
    # 4 - Return
    # --------------------------------------------------------------------------

    # Build return object
    out <- list(
        # Core output
        estimates = est,
        models = original$models,

        # Bootstrap information if available
        n_success_boot   = boot$n_success_boot,
        boot_errors  = boot$boot_errors,
        boot_nas     = boot$boot_nas,
        boot_samples     = if (keep_boot_samples) boot_samples[terms_keep] else NULL,

        # User-provided or default specifications
        outcome_time = outcome_time,
        outcome_status = outcome_status,
        exposure = exposure,
        exposure_time = exposure_time,
        immune_lag = tau,
        eval_times = eval_times,
        effect = effect,
        ci_type = ci_type,
        boot_reps = boot_reps,
        alpha = alpha,

        # Meta information
        call = call,
        descrip = descrip,
        method =  "matching"
    )

    class(out) <- "vefit"

    return(out)
}




