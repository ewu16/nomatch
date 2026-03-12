#' Estimate covariate means for population represented in `nomatch()` estimand
#'
#' @description
#' Estimate covariate means for the weighted principal strata population 
#' represented in `nomatch()` estimand. Since the principal strata populations
#' are not directly observed, simple summary statistics cannot be used. Instead,
#' the means are estimated as weighted means where the weights depend on 
#' the both the original estimand weights and on 
#' predicted probabilities of belonging to different principal strata.
#' Currently only implemented for the case when the estimand weights 
#' represent the observed distribution of exposure days and covariates among the exposed.
#'
#' @param data A data frame. The same data passed to [nomatch()].
#' @param fit A `nomatchfit` object created by [nomatch()]. Must have been
#'   fit with `keep_models = TRUE` (the default) and `weights_source = "observed"`. 
#'
#' @return A data frame of the summary statistics with columns:
#' \describe{
#'   \item{variable}{Covariate name}
#'   \item{label}{Factor level label (empty string for numeric variables)}
#'   \item{mean}{Mean (or weighted proportion for factor levels)}
#' }
#'
#' @export
#'
#' @examples
#' fit <- nomatch(
#'   data = simdata,
#'   outcome_time = "Y",
#'   outcome_status = "event",
#'   exposure = "V",
#'   exposure_time = "D_obs",
#'   covariates = c("x1", "x2"),
#'   timepoints = seq(30, 180, by = 30),
#'   immune_lag = 14
#' )
#' summarize_nomatch_population(simdata, fit)
#' 
summarize_nomatch_population <- function(data, fit){
    # Input checks
    if (!is.nomatchfit(fit)) {
        stop("`fit` must be a nomatchfit object.", call. = FALSE)
    }
    if (!grepl("nomatch", fit$method)) {
        stop("`make_table1()` is only implemented for nomatch() fits.", call. = FALSE)
    }
    if (is.null(fit$model_0)) {
        stop("`fit` must be run with `keep_models = TRUE`.", call. = FALSE)
    }
    if (!is.data.frame(data) || nrow(data) == 0) {
        stop("`data` must be a non-empty data frame.", call. = FALSE)
    }
    if(fit$weights_source != "observed"){
        stop("Currently only implemented for when observed weights are used")
    }
    
    # Extract info from fit 
    covariates <- fit$covariates
    outcome_time <- fit$outcome_time
    exposure <- fit$exposure
    exposure_time <- fit$exposure_time
    tau <- fit$immune_lag
    
    model_0 <- fit$model_0
    surv_vars <- all.vars(stats::formula(model_0)[[2]])
    time_var <- surv_vars[1]
    event_var <- surv_vars[2]
    
    # Restrict to exposed and at-risk past tau
    data_pop <- data[
        data[[exposure]] == 1 & 
        data[[outcome_time]] - data[[exposure_time]] > tau, ,
        drop = FALSE
    ]
    
    # Build newdata for survival prediction at D* + tau
    newdata <- data_pop[,  c(covariates, exposure_time)]
    newdata[[time_var]] <- newdata[[exposure_time]] + tau
    newdata[[event_var]] <- 0
    newdata$surv <- stats::predict(model_0, newdata, type = "survival")
    
    # Expand covariates to dummy variables
    dummy_formula <- stats::as.formula(paste("~", paste(covariates, collapse = " + ")))
    newdata_mat   <- stats::model.matrix(dummy_formula, newdata)
    newdata_mat   <- newdata_mat[, colnames(newdata_mat) != "(Intercept)", drop = FALSE]
    dummy_vars    <- colnames(newdata_mat)
    
    newdata_dummy <- cbind(
        as.data.frame(newdata_mat),
        newdata[, !names(newdata) %in% covariates, drop = FALSE]
    )
    
    # Compute exposure-time specific weights 
    groups <- split(newdata_dummy, newdata_dummy[[exposure_time]])
    add_weights <- function(g) {
        g$n      <- nrow(g)
        g$sum    <- sum(g$surv)
        g$denom  <- g$sum / g$n
        g$weight <- g$surv / g$denom
        g
    }
    newdata_dummy <- do.call(rbind, lapply(groups, add_weights))
    
    # Weighted means for each dummy variable
    mean_est <- sapply(newdata_dummy[dummy_vars], function(x) mean(x * newdata_dummy$weight))
    
    # Fill in the dropped reference level for each factor covariate
    factor_covariates <- covariates[!sapply(data[covariates], is.numeric)]
    missing_levels    <- rep(NA_real_, length(factor_covariates))
    
    for (i in seq_along(factor_covariates)) {
        var           <- factor_covariates[i]
        all_levels    <- sort(paste0(var, unique(data[[var]])))
        present       <- grep(var, names(mean_est), value = TRUE)
        missing_levels[i] <- 1 - sum(mean_est[present])
        names(missing_levels)[i] <- setdiff(all_levels, present)
    }
    
    final <- c(mean_est, missing_levels)
    final <- final[sort(names(final))]
    
    # Format output
    mean_df          <- data.frame(name = names(final), mean = unname(final),
                                   stringsAsFactors = FALSE)
    pattern_vars     <- paste0(paste0("^", covariates), collapse = "|")
    mean_df$label    <- gsub(pattern_vars, "", mean_df$name)
    pattern_labels   <- paste0(paste0(mean_df$label, "$"), collapse = "|")
    mean_df$variable <- gsub(pattern_labels, "", mean_df$name)
    
    mean_df[, c("variable", "label", "mean")]
}

