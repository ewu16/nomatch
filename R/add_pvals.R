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

# compute_percentile_pval <- function(term, term_boot){
#     null_val <- switch(term,
#                        "risk_difference" = 0,
#                        "risk_ratio" = 1,
#                        "relative_risk_reduction" = 0,
#                        NULL)
#
#     if(is.null(null_val)){
#         return(NA)
#     }
#
#    apply(term_boot, 2, \(x) mean(abs(x) >= null_val))
# }
