#'Main function to estimate marginal cumulative incidences and derived effect
#'measures without matching
#'
#'@description `nomatch()` estimates marginal cumulative incidences under
#'  exposure and no exposure using a G-computation approach. The method fits two
#'  conditional hazard models- one for each exposure type- and uses
#'  these models to predict time- and covariate-specific cumulative incidences.
#'  These cumulative incidences are then marginalized to compute overall
#'  (marginal) cumulative incidences. By default, the cumulative incidences
#'  are marginalized over the observed distribution of exposure times and covariates
#'  among the exposed. The resulting cumulative incidences can be
#'   summarized as risk ratios (RR = 1 - risk_exposed/risk_unexposed),
#'   vaccine effectiveness  (VE = 1 - RR), or risk differences (RD = risk_exposed - risk_unexposed).
#'
#'@param data A data frame with one row per individual containing the columns
#'  named in `outcome_time`, `outcome_status`, `exposure`, `exposure_time`, and
#'  `covariates`. Missing values for all columns except `exposure_time` are not allowed.
#'@param outcome_time Name of the follow-up time for the outcome of interest, i.e.
#'  time to either the event or right-censoring, whichever occurs first. Time should
#'  be measured from a chosen time origin (e.g. study start, enrollment, or age).
#'@param outcome_status Name of the event indicator for the outcome. The underlying column should be
#'  numeric (`1` = event, `0` = censored).
#'@param exposure Name of the exposure indicator. The underlying column should
#'  be numeric (`1` = exposed during follow-up, `0` = never exposed during
#'  follow-up).
#'@param exposure_time Name of the time to exposure, measured on the same time scale
#'  as that used for `outcome_time`. Must be a non-missing numeric value exposed individuals
#'  and must be set to `NA` for unexposed individuals.
#'@param covariates Character vector of covariates to adjust for when fitting
#'  the hazard models. Include all known, measured confounders of
#'  exposure and censoring.  Covariates must be measured or defined at the chosen time origin.
#'@param immune_lag Non-negative numeric value specifying the time after exposure (sometimes denoted by `tau`) that
#'  should be excluded from the risk evaluation period. This argument is
#'  primarily intended for vaccination exposures, where it is common to exclude
#'  the time after vaccination when immunity is still building. Time must be
#'  measured in the same units as that used for `outcome_time` and `exposure_time`  (e.g. days)
#'  and should reflect the biological understanding of when vaccine-induced
#'  immunity develops (usually 1-2 weeks). For non-vaccine exposures, ` immune_lag` can
#'  be set to 0 (no delay period for evaluating risk).
#'@param eval_times Numeric vector specifying the timepoints at which to compute
#'  cumulative incidence and the derived effect measures. The timepoints should
#'  be expressed in terms of time since exposure and use the same units
#'  as that used for `outcome_time` and `exposure_time` (e.g. days). All values must be greater
#'  than ` immune_lag` and and should correspond to clinically meaningful follow-up
#'  durations, such as 30, 60, or 90 days after exposure. A fine grid of
#'  timepoints (e.g., `eval_times = (immune_lag + 1):100`) can be provided if cumulative
#'  incidence curves over time are desired.
#'@param effect Character. Primary effect measure of interest (a contrast of
#'  the marginal cumulative incidences). Either
#'  `"risk_ratio"`(default), `"vaccine_effectiveness"`  or `"risk_difference"`.
#'@param weights_source Character string specifying the type of marginalizing weights
#'  to use. Either:
#'   - `"observed"` (default): use  the empirical
#'  distribution of exposure times and covariates among the exposed as the marginalizing weights. This
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
#'   - `"both"`: Computes and returns both Wald and percentile intervals.
#'
#'@param boot_reps Number of bootstrap replicates for confidence intervals.
#'  Recommended to use at least 1000 for publication-quality results. Use
#'  smaller values (e.g., 10-100) for initial exploration. Default: `0` (no
#'  bootstrapping).
#'@param alpha Significance level for confidence intervals (Confidence level =
#'  100*(1-alpha)%). Default: `0.05`.
#'@param keep_models Logical; return the two fitted hazard models used to compute
#' cumulative incidences?
#'  Default: `TRUE`.
#'@param keep_boot_samples Logical; return bootstrap samples? Default:
#'  `TRUE`. Must be set to `TRUE` if user plans to use [add_simultaneous_ci()]
#'  to obtain simultaneous confidence intervals.
#'@param n_cores Integer; parallel cores for bootstrapping passed to
#'  `parallel::mclapply` as `mc.cores`. On Unix-like OS only; not available on
#'  Windows. Default: `1`.
#'  
#'@param formula_unexposed A character specification of the right-hand side of the formula 
#' to use for fitting the hazard model for the unexposed. The set of variables in the formula
#' must be identical to the set of variables in `covariates`. 
#' 
#' @param formula_exposed A character specification of the right-hand side of the formula 
#' to use for fitting the hazard model for the exposed. The set of variables in the formula
#' must contain the variables in `covariates` and `exposure_time`. It is recommended
#' that `exposure_time` be modeled flexibly (e.g. using splines). 
#' 
#'@param n_cores Integer; parallel cores for bootstrapping passed to
#'  `parallel::mclapply` as `mc.cores`. On Unix-like OS only; not available on
#'  Windows. Default: `1`.
#'
#'
#'@return An object of class `vefit` containing:
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
#'   \item{boot_samples}{(If `keep_boot_samples = TRUE`) Named list of bootstrap draws
#'   (stored as matrices) for each term. Rows index bootstrap replicates and columns index `eval_times`.}
#' }
#'
#' The `vefit` object has methods for [print()], [summary()], and [plot()].
#' Use [add_simultaneous_ci()] to add simultaneous confidence intervals.
#'
#'
#'@details
#'
#' **Modeling.** Two Cox proportional hazards models are fit to estimate
#' exposure-specific cumulative incidences. For the unexposed model, the outcome is modeled
#' on the original time scale and includes all individuals, with exposed individuals
#' censored at their time of exposure. For the exposed model, the outcome is modeled in
#' terms of time since exposure among individuals who were
#' exposed and who remained at risk `tau` days after exposure. Both models adjust for
#' the specified covariates to help control for confounding. The second model
#' also flexibly adjusts for exposure time (by default, as a natural cubic spline with 4
#' degrees of freedom) to capture time-varying background risk. Predicted risks
#' from both models are then marginalized over the specified
#' exposure-time and covariate distributions to obtain G-computation style cumulative
#' incidence estimates.
#'
#'
#'**Marginalizing weights.** When `weights_source = "observed"`, the marginalizing weights
#'are the empirical distributions of exposure times and covariates among the
#'exposed who remain at-risk `tau` days after exposure. These weights are
#'returned in the `vefit` object under `weights`. They can also be obtained
#'prior to the call to `nomatch()` by calling `get_observed_weights()`.
#'
#'
#' **Confidence intervals.** Wald CIs are constructed on transformed scales:
#'\eqn{\text{logit}} for cumulative incidence; \eqn{\log{RR}} for risk ratios,
#'\eqn{\log{1 - VE}} for vaccine effectiveness, using bootstrap SEs. These are
#'then back-transformed to the original scale. No transformation is used for risk differences.
#'
#' **Parallelization.** Bootstraps can be parallelized on Unix via [parallel::mclapply()]
#'by providing `n_cores` argument.
#'
#'@export
#'
#' @examples
#' # Fit vaccine effectiveness model using simulated data
#'
#' fit <- nomatch(
#'   data = simdata,
#'   outcome_time = "Y",
#'   outcome_status = "event",
#'   exposure = "V",
#'   exposure_time = "D_obs",
#'   covariates = c("x1", "x2"),
#'   eval_times = seq(30, 180, by = 30),
#'   immune_lag = 14,
#'   boot_reps = 5,
#'   n_cores = 1
#' )
#'
#' # View basic results
#' fit$estimates

nomatch <- function(data,
                  outcome_time,
                  outcome_status,
                  exposure,
                  exposure_time,
                  covariates,
                  immune_lag,
                  eval_times,
                  effect = c("risk_ratio", "vaccine_effectiveness", "risk_difference"),
                  weights_source = c("observed", "custom"),
                  custom_weights = NULL,
                  ci_type = c("wald", "percentile", "both"),
                  boot_reps = 0,
                  alpha = 0.05,
                  keep_models = TRUE,
                  keep_boot_samples = TRUE,
                  n_cores = 1,
                  formula_unexposed = NULL,
                  formula_exposed = NULL
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
    validate_nomatch_inputs(
        data = data,
        outcome_time = outcome_time,
        outcome_status = outcome_status,
        exposure = exposure,
        exposure_time = exposure_time,
        covariates = covariates,
        immune_lag = tau,
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

     descrip <- get_basic_descriptives_nomatch(
         data = data,
         outcome_time = outcome_time,
         outcome_status = outcome_status,
         exposure = exposure,
         exposure_time = exposure_time,
         tau = tau
         )
     
     hazard_formulas <- resolve_hazard_formulas(
         formula_unexposed = formula_unexposed,
         formula_exposed = formula_exposed,
         covariates = covariates,
         exposure_time = exposure_time
         )
     
     validate_formulas(
         formula_unexposed = hazard_formulas$formula_unexposed,
         formula_exposed = hazard_formulas$formula_exposed,
         covariates = covariates,
         exposure_time = exposure_time
        )

     

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
                             custom_gp_list = custom_gp_list,
                             formula_0 = hazard_formulas$formula_unexposed,
                             formula_1 = hazard_formulas$formula_exposed
                            )

     original <- do.call(get_one_nomatch, estimation_args)

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
         model_0 = if(keep_models) original$model_0 else NULL,
         model_1 = if(keep_models) original$model_1 else NULL,
         weights = weights,

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
         method =  "nomatch (G-computation)"
         )

     class(out) <- "vefit"

    return(out)
}


