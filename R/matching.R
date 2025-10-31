#' Compute marginal cumulative incidence estimates for a matched cohort
#'
#' @description This function is the main function for estimating cumulative incidence
#' in a matched dataset based on Kaplan Meier estimation. It creates the
#' analysis dataset and performs the analysis for the matched cohort.
#'
#' @inheritParams nomatch
#' @inheritParams clean_matched_data
#' @param matched_data A data frame for the matched cohort created using [match_rolling_cohort()].

#'@return An object of class `nomatchfit` containing:
#' \describe{
#'   \item{estimates}{Named list of matrices containing the cumulative incidence and
#'   effect estimates.
#'   \describe{
#'      \item{`cuminc_0`}{ marginal cumulative incidence under no exposure}
#'      \item{`cuminc_1`}{ marginal cumulative incidence under exposure}
#'      \item{`risk_difference`}{ `cuminc_1 - cuminc_0`}
#'      \item{`risk_ratio`}{ `cuminc_1/cuminc_0`}
#'      \item{`vaccine_effectivess`}{ `1 - risk_ratio`}
#'   }
#'      Each matrix has one row per value in `timepoints` and columns including the
#'     point estimate (`estimate`) and, when requested, confidence limits of the form
#'     (`{wald/percentile}_lower`, `{wald/percentile}_upper`). }
#'   \item{models}{Fitted Kaplan Meier}
#'   \item{n_success_boot}{Integer vector indicating the
#'   number of successful bootstrap replications per timepoint.}
#'   \item{boot_samples}{(If `keep_boot_samples = TRUE`) Named list of bootstrap draws
#'   (stored as matrices) for each term. Rows index bootstrap replicates and columns index `timepoints`.}
#' }
#'
#'
#' @export
#'
matching <- function(matched_data,
                        outcome_time,
                        outcome_status,
                        exposure,
                        exposure_time,
                        immune_lag,
                        timepoints,
                        ci_type = c("wald", "percentile", "both"),
                        boot_reps = 0,
                        alpha = 0.05,
                        keep_models = TRUE,
                        keep_boot_samples = TRUE,
                        n_cores = 1){

    call <- match.call()
    ci_type <- match.arg(ci_type)
    tau <- immune_lag


    # Check data/inputs
    validate_matching_inputs(
        data = matched_data,
        outcome_time = outcome_time,
        outcome_status = outcome_status,
        exposure = exposure,
        exposure_time = exposure_time,
        immune_lag = tau,
        timepoints = timepoints
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
                            timepoints = timepoints)

    original <- do.call(get_one_matching, estimation_args)

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
    if(is.null(ci_est)){
        pt_est  <- original$pt_estimates
        est <- lapply(colnames(pt_est), \(term) cbind(estimate = pt_est[,term]))
        names(est) <- colnames(pt_est)

    }else{
        est <- lapply(names(ci_est), \(term){
            #Add p-values
            wald <- ci_type %in% c("wald", "both")
            wald_pval <-       if(wald) compute_wald_pval(term, ci_est[[term]]) else NULL

            x <- cbind(ci_est[[term]],
                       wald_pval = wald_pval)

            #Format column order
            col_order <- c("estimate",
                           paste0("wald_", c("lower", "upper", "se", "pval", "n")),
                           paste0("percentile_", c("lower", "upper", "pval", "n")))

            x[, intersect(col_order, colnames(x)), drop = FALSE]

        })
        names(est) <- names(ci_est)
    }


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
        boot_samples     = if (keep_boot_samples) boot_samples else NULL,

        # User-provided or default specifications
        outcome_time = outcome_time,
        outcome_status = outcome_status,
        exposure = exposure,
        exposure_time = exposure_time,
        immune_lag = tau,
        timepoints = timepoints,
        ci_type = ci_type,
        boot_reps = boot_reps,
        alpha = alpha,

        # Meta information
        call = call,
        descrip = descrip,
        method =  "matching"
    )

    class(out) <- "nomatchfit"

    return(out)
}




