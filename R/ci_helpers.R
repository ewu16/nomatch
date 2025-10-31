#' Compute Wald or percentile bootstrapped confidence interval
#'
#' @description
#' Internally calls [compute_wald_ci()] and/or [compute_percentile_ci()]
#' depending on the value of `ci_type`
#'
#'
#' @param x A numeric vector of point estimates
#' @param boot_x A matrix of bootstrapped estimates where the rows are the
#' bootstrap iterations and the columns are the time points of interest.
#' @param ci_type Character string indicating which type of confidence interval
#'   to return ("wald", "percentile", "both")
#' @param alpha  Significance level used to compute confidence intervals.
#'   Confidence intervals have nominal level `1 - alpha`.
#' @param transform If `ci_type = "wald"`, a character string indicating the scale on which
#' to compute Wald confidence intervals that are then transformed back to the original scale.
#' Options are `logit` for the transformation log(x/(1-x)) or `log_ve` for the transformation
#' log(1-x).
#' @param z_star If `ci_type = "wald"`, a specific critical value used to
#' to compute Wald confidence intervals (assumed to be positive). If used, `alpha` argument is ignored.
#'
#' @return A matrix containing the lower and upper confidence intervals
#'
#' @keywords internal
#' @noRd
#'
compute_boot_ci <- function(x, boot_x, ci_type, alpha = .05, transform = NULL, z_star = NULL){
    if(ci_type == "wald" | ci_type == "both"){
         wald <- compute_wald_ci(x, boot_x, alpha, transform, z_star)
         percentile <- NULL
    }
    if(ci_type == "percentile" | ci_type == "both"){
        percentile <- compute_percentile_ci(boot_x, alpha)
        wald <- if(ci_type == "percentile") NULL else wald

    }

    list(ci = cbind(wald$ci, percentile),
         transform = wald$transform)
}


#' Internal transformations
#'
#' These define how estimates are transformed before applying
#' the Wald CI (e.g., taking log or logit).
#' These define how to back-transform point estimates
#' Functions take arguments `(eta, sd, z)` and return values
#' on the original scale.
#' @keywords internal
.transformations <- list(
    "logit" = list(
        fwd = stats::qlogis,
        lower = function(eta, sd, z) stats::plogis(eta - z * sd),
        upper = function(eta, sd, z) stats::plogis(eta + z * sd)
    ),
    "log_ve" = list(
        fwd = function(y) log(1-y),
        lower = function(eta, sd, z) 1 - exp(eta + z * sd),
        upper = function(eta, sd, z) 1 - exp(eta - z * sd)
    ),
    "log_rr" = list(
        fwd = log,
        lower = function(eta, sd, z) exp(eta - z * sd),
        upper = function(eta, sd, z) exp(eta + z * sd)
    ),
    "identity" = list(fwd = identity,
                  lower = function(eta, sd, z) eta - z * sd,
                  upper = function(eta, sd, z) eta + z * sd
    )
)

#' Internal map of term → transformation type
#' @keywords internal
.term_transformations <- list(
    cuminc_0              = .transformations[["logit"]],
    cuminc_1              = .transformations[["logit"]],
    risk_difference       = .transformations[["identity"]],
    risk_ratio            = .transformations[["log_rr"]],
    relative_risk_reduction = .transformations[["log_ve"]]
)



#' @rdname compute_boot_ci
#' @keywords internal
#' @noRd
compute_wald_ci <- function(x, boot_x,  alpha = .05, transform, z_star = NULL){

    z <- if (!is.null(z_star)) z_star else stats::qnorm(1 - alpha/2)

    # transform
    eta_x    <-  transform$fwd(x)
    eta_boot <-  transform$fwd(boot_x)

    # bootstrap SDs/counts ignoring non-finite draws
    col_sd  <- apply(eta_boot, 2, function(col) stats::sd(col[is.finite(col)], na.rm = TRUE))
    boot_n  <- apply(eta_boot, 2, function(col) sum(is.finite(col)))

    # back-transform CI limits
    lower <- transform$lower(eta_x, col_sd, z)
    upper <- transform$upper(eta_x, col_sd, z)

    transform$eta_boot <- eta_boot

    ci_out <- cbind(wald_lower = lower, wald_upper = upper, wald_se = col_sd, wald_n = boot_n)

    list(ci =  ci_out,
         transform = transform)
}


#' @rdname compute_boot_ci
#' @keywords internal
#' @noRd
compute_percentile_ci <- function(boot_x, alpha = .05){
    lower <-  apply(boot_x, 2, \(x) stats::quantile(x, alpha/2, na.rm = TRUE))
    upper <-  apply(boot_x, 2, \(x) stats::quantile(x, 1 - alpha/2, na.rm = TRUE))
    boot_n <- apply(boot_x, 2, \(x) sum(!is.na(x)))

    cbind(percentile_lower = lower, percentile_upper = upper, percentile_n = boot_n)
}
