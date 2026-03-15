#'Main function to estimate marginal cumulative incidences and derived effect
#'measures without matching
#'
#'@description `nomatch()` estimates marginal cumulative incidences under
#'  exposure and no exposure using a G-computation approach. The method fits two
#'  conditional hazard models- one for the exposed group and one for the unexposed group- and uses
#'  these models to predict time- and covariate-specific cumulative incidences.
#'  These time- and covariate-specific cumulative incidences are then marginalized to compute overall
#'  (marginal) cumulative incidences. By default, the cumulative incidences
#'  are marginalized over the observed distribution of exposure times and covariates
#'  among the exposed. The resulting cumulative incidences can be
#'   summarized as risk ratios (RR = 1 - risk_exposed/risk_unexposed),
#'   relative risk reduction (1 - RR), or risk differences (RD = risk_exposed - risk_unexposed).
#'
#'@param data A data frame with one row per individual containing the columns
#'  named in `outcome_time`, `outcome_status`, `exposure`, `exposure_time`, and
#'  `covariates`. Missing values in any of these columns except `exposure_time` are not allowed.
#'@param outcome_time Name of the follow-up time for the outcome of interest, i.e.
#'  time to either the event or right-censoring, whichever occurs first. Time should
#'  be measured from a chosen time origin (e.g. study start, enrollment, or age).
#'@param outcome_status Name of the event indicator for the outcome. The underlying column should be
#'  numeric (`1` = event, `0` = censored).
#'@param exposure Name of the exposure indicator. The underlying column should
#'  be numeric (`1` = exposed during follow-up, `0` = never exposed during
#'  follow-up).
#'@param exposure_time Name of the time to exposure, measured on the same time scale
#'  as that used for `outcome_time`. Must be a non-missing numeric value for exposed individuals
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
#'  be set to 0 (no delay period for evaluating risk). Default: 0, 
#'@param timepoints Numeric vector specifying the timepoints at which to compute
#'  cumulative incidence and the derived effect measures. The timepoints should
#'  be expressed in terms of time since exposure and use the same units
#'  as that used for `outcome_time` and `exposure_time` (e.g. days). All values must be greater
#'  than ` immune_lag` and and should correspond to clinically meaningful follow-up
#'  durations, such as 30, 60, or 90 days after exposure. A fine grid of
#'  timepoints (e.g., `timepoints = (immune_lag + 1):100`) can be provided if cumulative
#'  incidence curves over time are desired. By default, the sequence from `immune_lag + 1` to
#'  the maximum event time in the exposed group, by units of 1, is used. 
#'@param weights_source Character string specifying the type of marginalizing weights
#'  to use for marginalization. Either:
#'   - `"observed"` (default): use  the empirical
#'  distribution of exposure times and covariates among the exposed as the marginalizing weights. If `immune_lag > 0`, 
#'  this is the distribution of exposure times and covariates among the exposed who remain at-risk `immune_lag` days after exposure. 
#'  These weights provide close alignment with the weights implicitly used in matching.
#'   - `"custom"`: use the user-specified weights provided in the `custom_weights` argument.
#'   See **Marginalizing weights** section for more details. 
#'@param custom_weights a named `list(g_weights, p_weights)` providing weights for
#'  marginalizing the time- and covariate-specific cumulative incidences. Must
#'  have the following format:
#'   - `g_weights`: data frame with the following columns: 
#'      *  columns named after each variable in `covariates` 
#'      * `exposure_time` (time of exposure),
#'      * `prob` (probability of exposure at the given time within the covariate-group;
#'      should sum to 1 within each covariate-group)
#'   - `p_weights`: data frame with the following columns: 
#'      *  columns named after each variable in `covariates` 
#'      * `prob` (probability of covariate-group; should sum to 1 over all covariate groups.)
#'@param ci_type Method for constructing pointwise bootstrap confidence intervals. One of
#'  `"wald"`, `"percentile"`, or `"both"`.
#'   - `"wald"` (default): Computes Wald-style intervals using bootstrap standard errors.
#'   See **Confidence intervals** section for details.
#'   - `"percentile"`: Computes percentile bootstrap intervals.
#'   - `"both"`: Computes and returns both Wald and percentile intervals.
#'
#'@param boot_reps Number of bootstrap replicates for confidence intervals.
#'  Recommended to use at least 1000 for publication-quality results. Use
#'  smaller values (e.g., 10-100) for initial exploration. Default: `0` (no
#'  bootstrapping). Bootstrap procedure can be parallelizaed- see **Parallelization** section for details. 
#'  
#'@param alpha Significance level for confidence intervals (Confidence level =
#'  100*(1-alpha)%). Default: `0.05`.
#'@param keep_models Logical; return the two fitted hazard models used to compute
#' cumulative incidences?
#'  Default: `TRUE`.
#'@param keep_boot_samples Logical; return bootstrap samples? Default:
#'  `TRUE`. Must be set to `TRUE` if user plans to use [add_simultaneous_ci()]
#'  to obtain simultaneous confidence intervals.
#'  
#' @param seed Integer seed for reproducible bootstrap results. Default: `NULL`.
#'   
#'@param formula_unexposed,formula_exposed  Optional specification of the right-hand side of the formula 
#' to use for fitting the hazard model for the unexposed and exposed groups, respectively. 
#' Each accepts one of the following:
#' -  a one-sided formula object (e.g. `~ x1 + x2`)
#' -  a string representation of a formula (e.g. `"x1 + x2"` or `"~ x1 + x2"`),
#' or a character vector of term names (e.g. `c("x1", "x2")`). 
#' 
#' The set of variables included in these formulas must be included in `covariates`. 
#' In `formula_exposed`, it is strongly recommended to model `exposure_time` flexibly (e.g. with a spline).
#' 
#' By default, `formula_unexposed` is a main-effects formula using all `covariates` 
#' and `formula_exposed` is a main-effects formula using all `covariates` and a natural cubic spline of `exposure_time` (4 df).
#'
#' See **Modeling** section for details. 
#'
#'@return An object of class `nomatchfit` containing:
#' \describe{
#'   \item{estimates}{Named list of matrices containing the cumulative incidence and
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
#'     
#'   \item{model_0}{(If `keep_models = TRUE`) Fitted hazard model for the unexposed group.
#'   See **Modeling** section for details.}
#'   \item{model_1}{(If `keep_models = TRUE`) Fitted hazard model for the exposed group.
#'   See **Modeling** section for details.}
#'   
#'    \item{weights}{List with dataframes `g_weights`, `p_weights` specifying
#'   the marginalizing weights used for averaging over exposure times and covariates.}

#'   \item{n_success_boot}{Integer vector indicating the
#'   number of successful bootstrap replications per timepoint.}
#'   \item{boot_samples}{(If `keep_boot_samples = TRUE`) Named list containing the `original` bootstrap
#'   estimates and, if Wald confidence intervals are used, the `transformed` bootstrap draws, which are the same
#'   bootstrap estimates on the scales used for inference. Both `original` and `transformed` are named lists
#'   containing the bootstrap estimates for each term. The bootstrap estimates for each term are stored
#'   as matrices, with rows indexing bootstrap replicates and columns indexing `timepoints`.}
#'}
#'
#'The `nomatchfit` object has methods for [print()], [summary()], and [plot()].
#'Use [add_simultaneous_ci()] to add simultaneous confidence intervals.
#'
#'
#'@details
#'
#' **Modeling.** Two Cox proportional hazards models are fit to estimate
#'exposure-specific cumulative incidences. For the unexposed model, the outcome
#'is modeled on the original time scale and includes all individuals, with
#'exposed individuals censored at their time of exposure. For the exposed model,
#'the outcome is modeled on the time scale of time since exposure and includes
#'individuals who were exposed and who remained at risk `immune_lag` days after
#'exposure. By default, both models adjust for the specified covariates as
#'linear, main-effect terms to help control for confounding. The second model
#'also flexibly adjusts for exposure time (by default, as a natural cubic spline
#'with 4 degrees of freedom) to capture time-varying background risk. Custom
#'formulas specifiying the right hand side of the formulas to be used in these
#'models can be provided through the optional arguments `formula_unexposed` and
#'`formula_exposed`. Predicted risks from both models are then marginalized over
#'the specified exposure-time and covariate weights to obtain G-computation
#'style cumulative incidence estimates.
#'
#'
#'**Marginalizing weights.** When `weights_source = "observed"`, the marginalizing weights
#'are the empirical distributions of exposure times and covariates among the
#'exposed who remain at-risk `immune_lag` days after exposure. These weights are
#'returned in the `nomatchfit` object under `weights`. They can also be obtained
#'prior to the call to `nomatch()` by calling `get_observed_weights()`.
#'
#'
#' **Confidence intervals.** Wald and percentile confidence intervals are constructed
#'for cumulative incidence and effectiveness parameters at each timepoint. The
#'Wald pointwise confidence intervals are constructed on transformed scales:
#'\eqn{\text{logit}} for cumulative incidence, identity transformation for risk
#'differences, \eqn{\log{RR}} for risk ratios, and \eqn{\log{1 - RR}} for
#'relative risk reduction using bootstrap standard errors of the transformed
#'estimates. These confidence intervals are then back-transformed to the
#'original scale. Wald p-values test the null of no effect (RD = 0, RR = 1, 1-RR
#'= 0) on the same transformed scales used to construct the Wald confidence
#'intervals. To obtain simultaneous confidence intervals, use
#'[add_simultaneous_ci()] after saving the original fit.
#'
#'
#' **Parallelization.** Bootstraps can be parallelized using the `future`
#'framework. Set a parallel plan before calling `nomatch()`: e.g.
#'
#' ```r
#' future::plan(future::multisession, workers = 4)
#' fit <- nomatch(..., boot_reps = 1000, seed = 42)
#' future::plan(future::sequential)  # reset when done
#' ```
#'
#'If no plan is set, bootstraps run sequentially. `multisession` works on all
#'operating systems and is recommended for most users. See the [future package
#'documentation](https://future.futureverse.org) for additional plans and
#'details on setup.
#'
#'@export
#'
#' @examples
#' # Fit nomatch using simulated data
#'
#' fit <- nomatch(
#'   data = simdata,
#'   outcome_time = "Y",
#'   outcome_status = "event",
#'   exposure = "V",
#'   exposure_time = "D_obs",
#'   covariates = c("x1", "x2"),
#'   timepoints = seq(30, 90, by = 30),
#'   immune_lag = 14,
#'   boot_reps = 5,
#'   seed = 123
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
                  immune_lag = 0,
                  timepoints = NULL,
                  weights_source = c("observed", "custom"),
                  custom_weights = NULL,
                  ci_type = c("wald", "percentile", "both"),
                  boot_reps = 0,
                  alpha = 0.05,
                  keep_models = TRUE,
                  keep_boot_samples = TRUE,
                  seed = NULL,
                  formula_unexposed = NULL,
                  formula_exposed = NULL
                  ){

    # --------------------------------------------------------------------------
    # 0 - Prep
    # --------------------------------------------------------------------------
    # Normalize user choices
    call <- match.call()

    weights_source      <- match.arg(weights_source)
    ci_type   <- match.arg(ci_type)
    tau <- immune_lag

    # Validate inputs
    validate_data_inputs(
        data = data,
        core_args = list(outcome_time = outcome_time,
                         outcome_status = outcome_status,
                         exposure = exposure,
                         exposure_time = exposure_time),
        covariates = covariates
    )
    
    validate_immune_lag(immune_lag)
    timepoints <- resolve_timepoints(data, outcome_time, outcome_status, 
                                     exposure, timepoints, immune_lag)
    validate_timepoints(timepoints, immune_lag)
    
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

     if(identical(weights_source, "custom") ){
         validate_marginalizing_weights(
             custom_weights = custom_weights,
             exposure_time      = exposure_time,
             covariates    = covariates
         )

         # Format weights to improve efficiency of internal calls
         custom_gp_list <- canonicalize_weights(
             weights = custom_weights,
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
                             timepoints = timepoints,
                             custom_gp_list = custom_gp_list,
                             formula_0 = hazard_formulas$formula_unexposed,
                             formula_1 = hazard_formulas$formula_exposed
                            )

     original <- do.call(get_one_nomatch, estimation_args)

    # --------------------------------------------------------------------------
    # 2 - Add bootstrap confidence intervals to point estimates
    # --------------------------------------------------------------------------
    # Helper returns NULL if boot_reps = 0
     boot <- run_bootstrap_inference(
         one_boot_function  = one_boot_nomatch,
         one_boot_args      = estimation_args,
         ci_type            = ci_type,
         boot_reps          = boot_reps,
         pt_est             = original$pt_estimates,
         alpha              = alpha,
         seed               = seed
     )

     ci_est <-boot$ci_estimates
     boot_samples <- boot$boot_samples

     # --------------------------------------------------------------------------
     # 3 - Format results
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
         weights_source = weights_source, 
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
         timepoints = timepoints,
         ci_type = ci_type,
         boot_reps = boot_reps,
         alpha = alpha,

         # Meta information
         call = call,
         descrip = descrip,
         method =  "nomatch (G-computation)"
         )

     class(out) <- "nomatchfit"

    return(out)
}


