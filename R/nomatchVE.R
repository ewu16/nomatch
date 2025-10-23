#'Main function to estimate marginal cumulative incidences and derived effect
#'measures without matching
#'
#'@description `nomatchVE()` estimates marginal cumulative incidences under
#'  exposure and no exposure using a G-computation approach. The method fits two
#'  conditional hazard models- one for each exposure type- and uses
#'  these models to predict time- and covariate-specific cumulative incidences.
#'  The conditional cumulative incidences are then marginalized to compute overall
#'  (marginal) cumulative incidences. By default, the cumulative incidences
#'  are marginalized over the observed distribution of exposure times and covariates
#'  among the exposed. The resulting cumulative incidences can be
#'   summarized as risk ratios (RR = 1 - risk_exposed/risk_unexposed) or
#'   vaccine effectiveness  (VE = 1 - RR).
#'
#'@param data A data frame with one row per individual containing the columns
#'  named in `outcome_time`, `outcome_status`, `exposure`, `exposure_time`, and
#'  `covariates`.
#'@param outcome_time Name of the follow-up time for the outcome of interest, i.e.
#'  time to either the event or right-censoring, whichever occurs first. Time should
#'  be measured from a given time origin (e.g. study start, enrollment, or age)
#'  for all individuals.
#'@param outcome_status Name of the event indicator. The underlying column should be
#'  numeric (`1` = event, `0` = censored).
#'@param exposure Name of the exposure indicator. The underlying column should
#'  be numeric (`1` = exposed during follow-up, `0` = never exposed during
#'  follow-up).
#'@param exposure_time Name of the time to exposure, measured on the same time scale
#'  as that used for `outcome time`. Set `exposure_time` of all unexposed
#'  individuals to `NA`.
#'@param covariates Character vector of covariates to adjust for when fitting
#'  the hazard models; should include all known confounders of
#'  exposure and censoring.  Covariates must be measured or defined at the chosen time origin.
#'@param immune_lag Non-negative numeric value specifying the time after exposure (tau) that
#'  should be excluded from the risk evaluation period. This argument is
#'  primarily intended for vaccination exposures, where it is common to exclude
#'  the time after vaccination when immunity is still building. Time must be
#'  measured in the same units as that used for `outcome_time` and `exposure_time`
#'  and should reflect the biological understanding of when vaccine-induced
#'  immunity develops (usually 1-2 weeks). For non-vaccine exposures, ` immune_lag` can
#'  be set to 0 (no delay period).
#'@param eval_times Numeric vector specifying the timepoints at which to compute
#'  cumulative incidence and the derived effect measures. The timepoints should
#'  be expressed in terms of time since exposure. All values must be greater
#'  than ` immune_lag` and and should correspond to clinically meaningful follow-up
#'  durations, such as 30, 60, or 90 days after exposure. A fine grid of
#'  timepoints (e.g., `eval_times = (immune_lag + 1):100`) can be provided if cumulative
#'  incidence curves over time are desired.
#'@param effect Character. Type of effect measure to return (a contrast of
#'  the marginal cumulative incidences). Either
#'  `"vaccine_effectiveness"` (default) or `"risk_ratio"`.
#'@param weights_source Character string specifying the type of marginalizing weights
#'  to use. Either:
#'   - `"observed"` (default): set the marginalizing weights to the empirical
#'  distribution of exposure eval_times and covariates among the exposed. This
#'  provides close alignment with the weights implicitly used in matching.
#'   - `"custom"`: use the user-specified weights provided in the `custom_weights` argument.
#'@param custom_weights a `list(g_weights, p_weights)` providing weights for
#'  marginalizing the time- and covariate-specific cumulative incidences. Must
#'  have the following format:
#'   - `g_weights`: data frame with columns
#'      *  all variables in `covariates`
#'      * `exposure_time` (time of exposure),
#'      * `prob` (probability of exposure at the given time within the covariate-group;
#'      should sum to 1 within each covariate-group)
#'   - `p_weights`: data frame with columns
#'      *  all variables in `covariates`
#'      * `prob` (probability of covariate-group; should sum to 1 over all covariate groups.)
#'@param ci_type Method for constructing bootstrap confidence intervals. One of
#'  `"wald"`, `"percentile"`, or `"both"`.
#'   - `"wald"` (default): Computes Wald-style intervals using bootstrap standard errors.
#'   See **Confidence intervals** section for details.
#'   - `"percentile"`: Computes percentile bootstrap intervals.
#'   - `"both"`: Computes and returns both sets of intervals.
#'
#'@param boot_reps Number of bootstrap replicates for confidence intervals.
#'  Recommended to use at least 1000 for publication-quality results. Use
#'  smaller values (e.g., 10-100) for initial exploration. Default: `0` (no
#'  bootstrapping).
#'@param alpha Significance level for confidence intervals (Confidence level =
#'  100*(1-`alpha`)%). Default: `0.05`.
#'@param keep_models Logical; return the two fitted hazard models used to compute
#' cumulative incidences?
#'  Default: `TRUE`.
#'@param keep_boot_samples Logical; return bootstrap samples? Default:
#'  `TRUE`. Must be set to `TRUE` if user plans to use [add_simultaneous_ci()]
#'  to obtain simultaneous confidence intervals.
#'@param n_cores Integer; parallel cores for bootstrapping. Passed to
#'  `parallel::mclapply` as `mc.cores`. On Unix-like OS only; not available on
#'  Windows. Default: `1`.
#'
#'
#'@return An object of class `vefit` containing:
#' \describe{
#'   \item{estimates}{List of matrices:
#'   \describe{
#'      \item{`cuminc_0`}{ marginal cumulative incidence under no exposure}
#'      \item{`cuminc_1`}{ marginal cumulative incidence under exposure}
#'      \item{`<effect>`}{the selected effect measure}
#'   }
#'      Each matrix has one row per value in `eval_times` and columns including the
#'     point estimate (`estimate`) and, when requested, confidence limits of the form
#'     (`{wald/percentile}_lower`, `{wald/percentile}_upper`). }
#'   \item{weights}{List with dataframes `g_weights`, `p_weights` specifying
#'   the marginalizing weights used for averaging over exposure times and covariates.}
#'   \item{model_0}{Fitted hazard model for the unexposed group.
#'   See **Modeling** section for details.}
#'   \item{model_1}{Fitted hazard model for the exposed group.
#'   See **Modeling** section for details.}
#'   \item{n_success_boot}{Integer vector indicating the
#'   number of successful bootstrap replications per timepoint.}
#'   \item{boot_samples}{(If `keep_boot_samples = TRUE`) List of bootstrap draws
#'   (stored as matrices) for each
#'   returned quantity with names mirroring those in `estimates` (i.e. `cuminc_0`, `cuminc_1`, `<effect>`).
#'   Rows index bootstrap replicates and columns index `eval_times`.}
#' }
#'
#' The `vefit` object has methods for [print()], [summary()], and [plot()].
#' Use [add_simultaneous_ci()] to add simultaneous confidence intervals.
#'
#'
#'@details
#'
#' **Modeling.** Two Cox proportional hazards models are fit to estimate
#' exposure-specific cumulative incidences. The first models the hazard of the
#' outcome over the chosen time scale among individuals who have not yet been
#' exposed, with exposed individuals censored at their time of exposure. The
#' second models the hazard over time since exposure among individuals who were
#' exposed and remain at risk `tau` days after exposure. Both models adjust for
#' the specified covariates to help control for confounding. The second model
#' also flexibly adjusts for exposure time (as a natural cubic spline with 4
#' degrees of freedom) to capture time-varying background risk. Predicted risks
#' from both models are then marginalized over the specified covariate and
#' exposure-time distributions to obtain G-computation style cumulative
#' incidence estimates.
#'
#'
#'**Marginalizing weights.** When `weights_source = "observed"`, the marginalizing weights
#'are the empirical distributions of exposure eval_times and covariates among the
#'exposed who remain at-risk `tau` days after exposure. These weights are
#'returned in the `vefit` object under `gp_list`. They can also be obtained
#'prior to the call to `nomatchVE()` by calling `get_observed_weights()`.
#'
#'
#' **Confidence intervals.** Wald CIs are constructed on transformed scales:
#'\eqn{\text{logit}} for cumulative incidence; \eqn{\log{RR}} for risk ratios,
#'\eqn{\log{1 - VE}} for vaccine effectiveness, using bootstrap SEs. These are
#'then back-transformed to the original scale.
#'
#' **Parallelization.** Bootstraps can be parallelized on Unix via [parallel::mclapply()]
#'by providing `n_cores` argument.
#'
#'@export
#'
#' @examples
#' # Fit vaccine effectiveness model using simulated data
#' data(simdata)
#'
#' fit <- nomatchVE(
#'   data = simdata,
#'   outcome_time = "Y",
#'   outcome_status = "event",
#'   exposure = "V",
#'   exposure_time = "D_obs",
#'   covariates = c("x1", "x2"),
#'   eval_times = seq(30, 180, by = 30),
#'   immune_lag = 14,
#'   boot_reps = 5,
#'   n_cores = 2
#' )
#'
#' # View basic results
#' fit$estimates

nomatchVE <- function(data,
                  outcome_time,
                  outcome_status,
                  exposure,
                  exposure_time,
                  covariates,
                  immune_lag,
                  eval_times,
                  effect = c("vaccine_effectiveness", "risk_ratio"),
                  weights_source = c("observed", "custom"),
                  custom_weights = NULL,
                  ci_type = c("wald", "percentile", "both"),
                  boot_reps = 0,
                  alpha = 0.05,
                  keep_models = TRUE,
                  keep_boot_samples = TRUE,
                  n_cores = 1
                  ){

    # --------------------------------------------------------------------------
    # 0 - Prep
    # --------------------------------------------------------------------------

    # Normalize user choices
    call <- match.call()

    effect <- match.arg(effect)
    weights_source      <- match.arg(weights_source)
    ci_type   <- match.arg(ci_type)
    tau <- immune_lag

    # Validate inputs
    validate_ve_inputs(
        data = data,
        outcome_time = outcome_time,
        outcome_status = outcome_status,
        exposure = exposure,
        exposure_time = exposure_time,
        covariates = covariates,
        tau = tau,
        eval_times = eval_times
    )

     if(identical(weights_source, "custom") ){
         validate_marginalizing_weights(
             custom_weights = custom_weights,
             exposure_time      = exposure_time,
             covariates    = covariates
         )

         # Format weights to improve efficiency of internal calls
         custom_gp_list <- canonicalize_weights(
             weights = weights,
             exposure_time  = exposure_time,
             covariates    = covariates
         )
     }else{
        custom_gp_list <- NULL
     }

     descrip <- get_basic_descriptives_nomatch(data,
                                       outcome_time = outcome_time,
                                       outcome_status = outcome_status,
                                       exposure = exposure,
                                       exposure_time = exposure_time,
                                       tau = tau)

     # --------------------------------------------------------------------------
     # 1 - Get original estimate
     # --------------------------------------------------------------------------
     estimation_args <- list(data = data,
                             outcome_time = outcome_time,
                             outcome_status = outcome_status,
                             exposure = exposure,
                             exposure_time = exposure_time,
                             covariates = covariates,
                             tau = tau,
                             eval_times = eval_times,
                             custom_gp_list = custom_gp_list
                            )

     original <- do.call(get_one_nomatch_ve, estimation_args)

    # --------------------------------------------------------------------------
    # 2 - Add bootstrap confidence intervals to point estimates
    # --------------------------------------------------------------------------
    # Helper returns NULL if boot_reps = 0
     boot <- estimate_bootstrap_ci(
         one_boot_function  = one_boot_nomatch,
         one_boot_args      = estimation_args,
         ci_type            = ci_type,
         boot_reps          = boot_reps,
         pt_est             = original$pt_estimates,
         alpha              = alpha,
         keep_boot_samples  = keep_boot_samples,
         n_cores            = n_cores
     )

     ci_est <-boot$ci_estimates
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

     # Define the weights to return
     if(identical(weights_source, "custom")){
         weights <- custom_weights
     }else{
         weights <- gp_to_weights(original$gp_list)
     }

    # Build return object
     out <- list(
         # Core output
         estimates = est,
         model_0 = original$model_0,
         model_1 = original$model_1,
         weights = weights,

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
         covariates = covariates,
         immune_lag = tau,
         eval_times = eval_times,
         effect = effect,
         ci_type = ci_type,
         boot_reps = boot_reps,
         alpha = alpha,

         # Meta information
         call = call,
         descrip = descrip,
         method =  "nomatchVE (G-computation)"
         )

     class(out) <- "vefit"

    return(out)
}


