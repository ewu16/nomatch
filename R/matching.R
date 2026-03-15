#' Compute marginal cumulative incidence estimates for a matched cohort
#'
#' @description This function is the main function for estimating cumulative
#'   incidence in a matched dataset based on Kaplan Meier estimation. It creates
#'   the analysis dataset and performs the analysis for the matched cohort.
#'
#' @inheritParams nomatch
#' @inheritParams make_matched_adata
#' @param matched_data A data frame for the matched cohort created using
#'   [match_rolling_cohort()].
#' @param keep_models Logical; return the fitted Kaplan Meier object used to
#'   compute cumulative incidences? Default: `TRUE`

#'@return An object of class `nomatchfit` containing:
#' \describe{
#'   \item{matched_adata}{The analytic dataset used for the matching analysis, often
#'      a subset of the original matched dataset when `immune_lag > 0`. The data frame
#'      consists of the original variables in `matched_data` plus:
#'      \describe{
#'        \item{`match_<outcome_status>`}{Event indicator for the outcome in the matched analysis}
#'        \item{`match_<outcome_time>`}{Follow-up time for outcome of interest in the matched analysis, relative to original time-origin}
#'        \item{`match_T`}{Follow-up time for outcome of interest in the matched analysis, relative to time of matching}
#'       }
#'   }
#'  \item{estimates}{Named list of matrices containing the cumulative incidence and
#'   effect estimates:
#'   \describe{
#'      \item{`cuminc_0`}{ marginal cumulative incidence under no exposure}
#'      \item{`cuminc_1`}{ marginal cumulative incidence under exposure}
#'      \item{`risk_difference`}{ `cuminc_1 - cuminc_0`}
#'      \item{`risk_ratio`}{ `cuminc_1/cuminc_0`}
#'      \item{`relative_risk_reduction`}{ `1 - risk_ratio`}
#'   }
#'      Each matrix has one row per value in `timepoints` and columns for the
#'     point estimate (`estimate`) and confidence limits
#'     (`{wald/percentile}_lower`, `{wald/percentile}_upper`), when applicable.}
#'  \item{models}{Fitted Kaplan Meier object}
#'   \item{n_success_boot}{Integer vector indicating the
#'   number of successful bootstrap replications per timepoint.}
#'   \item{boot_samples}{(If `keep_boot_samples = TRUE`) Named list of bootstrap draws
#'   (stored as matrices) for each term. Rows index bootstrap replicates and columns index `timepoints`.}
#'}
#'
#'@details
#'  **Matching analysis**: To create
#'the matched analysis dataset, both individuals in the matched pair are
#'censored when the control receives the exposure, and matched pairs in which
#'either individual experiences the event or censoring with  `immune_lag` days
#'of matching is excluded. Cumulative
#'incidence in the matched analysis dataset is estimated using Kaplan Meier
#'stratified only on exposure group.
#'
#' **Confidence intervals.** Wald and percentile confidence intervals are constructed
#' for cumulative incidence and effectiveness parameters at each timepoint. 
#'  The Wald pointwise confidence intervals are constructed on transformed scales:
#' \eqn{\text{logit}} for cumulative incidence; \eqn{\log{RR}} for risk ratios, and \eqn{\log{1 - RR}} for relative risk 
#' reduction, using bootstrap standard errors. These confidence intervals are  
#' then back-transformed to the original scale. Identity transformation is used for risk differences.
#' To obtain simultaneous confidence intervals, use [add_simultaneous_ci()] after
#' saving the original fit. 
#'
#' **Parallelization.** Bootstraps can be parallelized using the `future`
#'framework. Set a parallel plan before calling `matching()`: e.g.
#'
#' ```r
#' future::plan(future::multisession, workers = 4)
#' fit <- matching(..., boot_reps = 1000, seed = 42)
#' future::plan(future::sequential)  # reset when done
#' ```
#'
#'If no plan is set, bootstraps run sequentially. `multisession` works on all
#'operating systems and is recommended for most users. See the [future package
#'documentation](https://future.futureverse.org) for additional plans and
#'details on setup.
#'
#'
#' @export
#' @examples
#' # Perform matching analysis on simulated data
#'
#' # Create matched data
#' matched_cohort <- match_rolling_cohort(
#'   data = simdata,
#'   outcome_time =  "Y",
#'   exposure = "V",
#'   exposure_time = "D_obs",
#'   matching_vars = c("x1", "x2"),
#'   id_name = "ID",
#'   seed = 5678
#' )
#'
#' # Analyze matched data
#' fit_matching <- matching(
#'   matched_data = matched_cohort$matched_data,
#'   outcome_time = "Y",
#'   outcome_status = "event",
#'   exposure = "V",
#'   exposure_time = "D_obs",
#'   timepoints = seq(30, 90, by = 30),
#'   immune_lag = 14,
#'   boot_reps = 5,
#'   seed = 123
#' )
#'
#' # View basic results
#' fit_matching
matching <- function(matched_data,
                        outcome_time,
                        outcome_status,
                        exposure,
                        exposure_time,
                        immune_lag = 0,
                        timepoints = NULL,
                        ci_type = c("wald", "percentile", "both"),
                        boot_reps = 0,
                        alpha = 0.05,
                        keep_models = TRUE,
                        keep_boot_samples = TRUE,
                        seed = NULL){

    call <- match.call()
    ci_type <- match.arg(ci_type)
    tau <- immune_lag


    # Check data/inputs
    validate_data_inputs(
        data = matched_data,
        core_args = list(outcome_time = outcome_time,
                         outcome_status = outcome_status,
                         exposure = exposure,
                         exposure_time = exposure_time)
    )

    validate_immune_lag(immune_lag)
    timepoints <- resolve_timepoints(matched_data, outcome_time, outcome_status, 
                                     exposure, timepoints, immune_lag)
    validate_timepoints(timepoints, immune_lag)
    
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
    boot <- run_bootstrap_inference(
        one_boot_function  = one_boot_matching,
        one_boot_args      = estimation_args,
        ci_type            = ci_type,
        boot_reps          = boot_reps,
        pt_est             = original$pt_estimates,
        alpha              = alpha,
        seed               = seed
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
            x <- cbind(ci_est[[term]])
            #Format column order
            col_order <- c("estimate",
                           paste0("wald_", c("lower", "upper", "pval", "n")),
                           paste0("percentile_", c("lower", "upper", "n")))
            
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
        matched_adata =  original$matched_adata,
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




