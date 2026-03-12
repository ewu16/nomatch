#' Compute Wald-style p-value based on bootstrap standard errors 
#'
#' @param term String specifying the effect measure of interest. 
#'  Either `"risk_ratio"`, `"relative_risk_reduction"`  or
#'   `"risk_difference"`.
#' @param term_ci_estimate A matrix containing the point estimate (`estimate`) and 
#' and bootstrap standard error (`wald_se`) for the specified term. 
#'
#' @return p-value 
#' @keywords internal 
#' @noRd
compute_wald_pval <- function(term, term_ci_estimate){
    null_val <- switch(term,
                       "risk_difference" = 0,
                       "risk_ratio" = 1,
                       "relative_risk_reduction" = 0,
                        NULL)

    if(is.null(null_val)){
        return(NA)
    }

    estimate <- term_ci_estimate[, "estimate"]
    se <- term_ci_estimate[, "wald_se"]

    z <- (estimate - null_val)/se

    2*(1 - stats::pnorm(abs(z)))
}

