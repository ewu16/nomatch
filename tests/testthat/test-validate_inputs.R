test_that("validate_nomatch_inputs() passes with good input", {
    expect_no_error(validate_nomatch_inputs(
        data = simdata,
        outcome_time = "Y",
        outcome_status = "event",
        exposure = "V",
        exposure_time = "D_obs",
        covariates = c("x1","x2"),
        immune_lag = 14,
        eval_times = c(30,60)
    ))
})

test_that("validate_nomatch_inputs() surfaces sub-validator errors", {
    # missing column
    dat <- simdata[, -which(names(simdata) == "x2")]
    expect_error(
        validate_nomatch_inputs(dat, "Y", "event", "V", "D_obs", c("x1","x2"), 14, c(30,60)),
        "Missing required column", ignore.case = TRUE
    )

    # invalid eval_times < immune_lag
    expect_error(
        validate_nomatch_inputs(simdata, "Y", "event", "V", "D_obs", c("x1","x2"), 14, c(7,14)),
        "must be > `immune_lag`", ignore.case = TRUE
    )
})

test_that("validate_match_rolling_cohort_inputs() passes good input", {
    expect_no_error(
        validate_match_rolling_cohort_inputs(
            data = simdata,
            outcome_time = "Y",
            exposure = "V",
            exposure_time = "D_obs",
            covariates = c("x1","x2")
        )
    )
})

test_that("validate_match_rolling_cohort_inputs() surfaces obvious errors", {
    dat <- simdata; dat$V[1] <- NA
    expect_error(
        validate_match_rolling_cohort_inputs(dat, "Y", "V", "D_obs", c("x1","x2")),
        "Missing values in <exposure>", fixed = TRUE
    )
})
