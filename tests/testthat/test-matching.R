test_that("matching() fails gracefully when outcomes are missing", {
    matched_data <- match_rolling_cohort(simdata, "Y", "V", "D_obs", c("x1","x2"), "ID")[[1]]
    matched_data$Y[1] <- NA

    expect_error(
        matching(matched_data, "Y", "event", "V", "D_obs", immune_lag = 14, timepoints = c(30,60)),
        "Missing values in `outcome_time`", ignore.case = TRUE
    )
})
