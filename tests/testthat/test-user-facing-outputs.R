# tests/test-nomatch-basic.R
test_that("basic nomatch run works", {
    set.seed(123)
    fit <- nomatch(simdata, "Y","event","V","D_obs", c("x1","x2"),
                   immune_lag = 14, timepoints = seq(30,180,30), boot_reps = 5)
    expect_s3_class(fit, "nomatchfit")
    expect_no_error(summary(fit))
    p <- plot(fit)
    expect_s3_class(p, "ggplot")
    
    df <- estimates_to_df(fit)
    expect_true(all(c("term","t0","estimate") %in% names(df)))
})

# Shared fixture
make_fit <- function(boot_reps = 5) {
    set.seed(123)
    nomatch(simdata, "Y", "event", "V", "D_obs", c("x1", "x2"),
            immune_lag = 14, timepoints = seq(30, 180, 30), boot_reps = boot_reps)
}

# nomatch() with bootstrap ------------------------------------------------
test_that("nomatch() print/summary/plot work without error with bootstrap", {
    fit <- make_fit()
    expect_no_error(print(fit))
    expect_output(print(fit), regexp = "Estimate")
    expect_no_error(summary(fit))
    expect_output(summary(fit), regexp = "Summary")
    expect_s3_class(plot(fit), "ggplot")
})

test_that("estimates_to_df() returns expected columns", {
    df <- estimates_to_df(make_fit())
    expect_true(all(c("term", "t0", "estimate") %in% names(df)))
})

# nomatch() without bootstrap ---------------------------------------------
test_that("nomatch() print/summary/plot work without error when boot_reps = 0", {
    fit <- make_fit(boot_reps = 0)
    expect_equal(colnames(fit$estimates$cuminc_0), "estimate")
    expect_no_error(print(fit))
    expect_no_error(summary(fit))
    expect_s3_class(plot(fit), "ggplot")
})

# matching() --------------------------------------------------------------
test_that("matching() print/summary/plot work without error", {
    matched <- match_rolling_cohort(simdata, "Y", "V", "D_obs",
                                    c("x1", "x2"), "ID")[[1]]
    fit <- matching(matched, "Y", "event", "V", "D_obs",
                    immune_lag = 14, timepoints = c(30, 60, 90), boot_reps = 5)
    expect_no_error(print(fit))
    expect_no_error(summary(fit))
    expect_s3_class(plot(fit), "ggplot")
})