test_that("add_simultaneous_ci() succeeds with >=2 timepoints and boot samples", {
    fit <- nomatch(simdata, "Y","event","V","D_obs", c("x1","x2"),
                     immune_lag = 14, timepoints = seq(30,180,30), boot_reps = 5)
    fit2 <- add_simultaneous_ci(fit)
    expect_s3_class(fit2, "nomatchfit")

    # Check that simult columns/flags appear (tailor to your object)
    df2 <- estimates_to_df(fit2)
    expect_true(any(grepl("simul", names(df2), ignore.case = TRUE)))
})

test_that("add_simultaneous_ci() succeeds with percentile confidence intervals", {
    fit <- nomatch(simdata, "Y","event","V","D_obs", c("x1","x2"),
                     immune_lag = 14, timepoints = seq(30,180,30),
                     ci_type = "percentile", boot_reps = 5)
    fit2 <- add_simultaneous_ci(fit)
    expect_s3_class(fit2, "nomatchfit")

    # Check that simult columns/flags appear (tailor to your object)
    df2 <- estimates_to_df(fit2)
    expect_true(any(grepl("simul", names(df2), ignore.case = TRUE)))
})

test_that("add_simultaneous_ci() errors clearly when only one timepoint", {
    fit_single <- nomatch(simdata, "Y","event","V","D_obs", c("x1","x2"),
                     immune_lag = 14, timepoints = 50, boot_reps = 5)
    expect_error(add_simultaneous_ci(fit_single),
                 "more than 1 timepoint|>=\\s*2", ignore.case = TRUE)
})

test_that("add_simultaneous_ci() errors clearly when no boot samples", {
    fit_no_boot <- nomatch(simdata, "Y","event","V","D_obs", c("x1","x2"),
                             immune_lag = 14, timepoints = seq(30,180,30), boot_reps = 0)
    expect_error(add_simultaneous_ci(fit_no_boot),
                 "bootstrap|boot", ignore.case = TRUE)
})
