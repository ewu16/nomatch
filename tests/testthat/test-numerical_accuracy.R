check_probability_bounds <- function(fit){
    cuminc_0 <- fit$estimates$cuminc_0[, "estimate"]
    cuminc_1 <- fit$estimates$cuminc_1[, "estimate"]
    risk_difference <- fit$estimates$risk_difference[, "estimate"]
    risk_ratio <- fit$estimates$risk_ratio[, "estimate"]
    relative_risk_reduction <- fit$estimates$relative_risk_reduction[, "estimate"]
    
    expect_true(all(cuminc_0 >= 0 & cuminc_0 <= 1))
    expect_true(all(cuminc_1 >= 0 & cuminc_1 <= 1))
    expect_true(all(risk_difference >= -1 & risk_difference <= 1))
    expect_true(all(risk_ratio >= 0))
    expect_true(all(relative_risk_reduction <= 1))
}

check_ci_ordering <- function(fit){
    ci_types <- c("wald")
    
    for (term in names(fit$estimates)) {
        cols <- colnames(fit$estimates[[term]])
        
        for (ci in ci_types) {
            lower_col <- paste0(ci, "_lower")
            upper_col <- paste0(ci, "_upper")
            if (!all(c(lower_col, upper_col) %in% cols)) next
            
            lower <- fit$estimates[[term]][, lower_col]
            upper <- fit$estimates[[term]][, upper_col]
            
            expect_true(all(lower <= upper, na.rm = TRUE),
                        label = paste(term, ci, "lower <= upper"))
        }
    }
}

test_that("nomatch estimates are numerically well-behaved", {
    fit <- make_nomatch_fit()  
    check_probability_bounds(fit)
    check_ci_ordering(fit)
    
})


test_that("matching estimates are numerically well-behaved", {
    fit <- make_matching_fit()  
    check_probability_bounds(fit)
    check_ci_ordering(fit)
    
})


# More specific numerical checks ------------------------------------------
test_that("cumulative incidence is 0 at t0 <= immune_lag", {
    fit <- make_nomatch_fit()  
    gp_list <- canonicalize_weights(fit$weights, "D_obs", c("x1","x2"))
    result <- compute_psi_bar_t0(fit$model_0, fit$model_1, 
                                 exposure_time = "D_obs", t0=5, tau=14, 
                                 newdata=gp_list$g_weights,
                                 gp_list = gp_list)
    expect_equal(unname(result), c(0,0))
})


test_that("marginalization with uniform weights equals simple mean", {
    psi_dx <- data.frame(
        group_id = c(1, 1, 1),
        D_obs    = c(1, 2, 3),
        psi_0_dx = c(0.1, 0.2, 0.3),
        psi_1_dx = c(0.05, 0.15, 0.25)
    )
    
    gp_list <- list(
        g_weights = data.frame(
            group_id = c(1, 1, 1),
            D_obs    = c(1, 2, 3),
            prob_g   = c(1/3, 1/3, 1/3)
        ),
        p_weights = data.frame(
            group_id = 1,
            prob_p   = 1
        )
    )
    
    weighted   <- marginalize_psi_dx_t0(psi_dx, gp_list)
    unweighted <- marginalize_psi_dx_t0(psi_dx, gp_list = NULL)
    
    expect_equal(weighted, unweighted)
})

