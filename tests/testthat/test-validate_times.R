# Tests for validate_immune_lag() -----------------------------------------
test_that("validate_immune_lag() enforces types, ranges, relationships, and warns on long lag", {
    expect_no_error(validate_immune_lag(immune_lag = 14))

    # immune_lag issues
    expect_error(validate_immune_lag(immune_lag = NULL), "single non-missing", ignore.case = TRUE)
    expect_error(validate_immune_lag(immune_lag = NA_real_), "single non-missing", ignore.case = TRUE)
    expect_error(validate_immune_lag(immune_lag = "a"), "numeric", ignore.case = TRUE)
    expect_error(validate_immune_lag(immune_lag = -1), "must be >= 0", ignore.case = TRUE)
    expect_warning(validate_immune_lag(immune_lag = 29), "unusually long", ignore.case = TRUE)
})


# Tests for resolve_timepoints() ------------------------------------------
test_that("resolve_timepoints() auto-resolves from data when NULL", {
    result <- resolve_timepoints(simdata, "Y", "event", "V", timepoints = NULL, immune_lag = 14)
    expect_true(all(result > 14))
    expect_equal(result, seq(15, max(simdata$Y[simdata$V == 1 & simdata$event == 1])))
})

test_that("resolve_timepoints() errors when max event time <= immune_lag", {
    dat <- simdata
    dat$Y[dat$V == 1 & dat$event == 1] <- 10  # force max below immune_lag
    expect_error(
        resolve_timepoints(dat, "Y", "event", "V", timepoints = NULL, immune_lag = 14),
        "immune_lag", ignore.case = TRUE
    )
})
test_that("resolve_timepoints() sorts provided timepoints", {
    result <- resolve_timepoints(simdata, "Y", "event", "V", timepoints = c(90, 30, 60), immune_lag = 14)
    expect_equal(result, c(30, 60, 90))
})

test_that("resolve_timepoints() deduplicates and warns on non-unique timepoints", {
    expect_warning(
        result <- resolve_timepoints(simdata, "Y", "event", "V", timepoints = c(30, 30, 60), immune_lag = 14),
        "Non-unique", ignore.case = TRUE
    )
    expect_equal(result, c(30, 60))
})


# Tests for validate_timepoints() -----------------------------------------
test_that("validate_timepoints() enforces structure, relationships", {
    expect_no_error(validate_timepoints(15:20, immune_lag = 14))
    expect_no_error(validate_timepoints(1:20, immune_lag = 14))
    
    # timepoints issues
    expect_error(validate_timepoints(timepoints = numeric(0), immune_lag = 14), "length > 0", ignore.case = TRUE)
    expect_error(validate_timepoints(timepoints = c(NA, 60), immune_lag = 14), "no missing values", ignore.case = TRUE)
    expect_error(validate_timepoints(timepoints = c(7, 14), immune_lag = 14), "must be > `immune_lag`", ignore.case = TRUE)
})


