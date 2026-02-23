matched_data <- match_rolling_cohort(simdata, "Y", "V", "D_obs", 
                                     c("x1","x2"), "ID", seed = 123)$matched_data


# Correct output ----------------------------------------------------------
test_that("matching() returns a nomatchfit with correct structure", {
    fit <- matching(matched_data, "Y", "event", "V", "D_obs",
                    immune_lag = 14, timepoints = seq(30, 180, 30), boot_reps = 0)
    
    expect_s3_class(fit, "nomatchfit")
    expect_setequal(names(fit$estimates),
                    c("cuminc_0", "cuminc_1", "risk_difference", "risk_ratio", "relative_risk_reduction"))
})

test_that("matching() results match gold standard", {
    set.seed(1234)
    fit <-  matching(matched_data, "Y", "event", "V", "D_obs",
                     immune_lag = 14, timepoints = seq(30, 180, 30), 
                     ci_type = "both", boot_reps = 5)
    gold_standard <- readRDS(test_path("fixtures", "gold_matching.rds"))
    expect_equal(fit$estimates, gold_standard$estimates)
})



# Confidence intervals ----------------------------------------------------
test_that("nomatch() adds CI columns when boot_reps > 0", {
    matching_with_ci <- function(ci_type){
        matching(matched_data, "Y", "event", "V", "D_obs",
                 immune_lag = 14, timepoints = seq(30, 90, 30), boot_reps = 5,
                 ci_type = ci_type)
    }
    fit <- matching_with_ci("wald")
    expect_true(all(c("wald_lower", "wald_upper") %in% 
                        colnames(fit$estimates$cuminc_0)))
    
    fit <- matching_with_ci("percentile")
    expect_true(all(c("percentile_lower", "percentile_upper") %in% 
                        colnames(fit$estimates$cuminc_0)))
    
    fit <- matching_with_ci("both")
    expect_true(all(c("wald_lower", "wald_upper", "percentile_lower", "percentile_upper") %in% 
                        colnames(fit$estimates$cuminc_0)))
    
})



# Defaults ----------------------------------------------------------------
test_that("matching() resolves timepoints = NULL without error", {
    fit <- matching(matched_data, "Y", "event", "V", "D_obs",
                    immune_lag = 14, timepoints = NULL, boot_reps = 0)
    expect_true(length(fit$timepoints) > 0)
    expect_true(all(fit$timepoints > 14))
})



# Failure modes -----------------------------------------------------------
test_that("matching() fails gracefully when outcomes are missing", {
    matched_data <- match_rolling_cohort(simdata, "Y", "V", "D_obs", c("x1","x2"), "ID")[[1]]
    matched_data$Y[1] <- NA

    expect_error(
        matching(matched_data, "Y", "event", "V", "D_obs", immune_lag = 14, timepoints = c(30,60)),
        "Missing values in `outcome_time`", ignore.case = TRUE
    )
})
