fit_nomatch <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1", "x2"),
                    immune_lag = 14, timepoints = seq(30, 180, 30), boot_reps = 5, seed = 456)

fit_matching <- local({
    m <- match_rolling_cohort(simdata, "Y", "V", "D_obs", c("x1","x2"), "ID", seed = 123)
    matching(m$matched_data, "Y", "event", "V", "D_obs",
             immune_lag = 14, timepoints = seq(30, 180, 30), boot_reps = 5, seed = 456)
})

# Correct output -------------------------------------------------------------
test_that("add_simultaneous_ci()  returns nomatchfit with simul columns and metadata", {
    fit2 <- add_simultaneous_ci(fit_nomatch, seed = 1)
    
    expect_s3_class(fit2, "nomatchfit")
    expect_true(all(c("simul_lower", "simul_upper") %in% colnames(fit2$estimates$cuminc_0)))
    
    expect_false(is.null(fit2$simul_z_star))
    expect_false(any(is.null(fit2$simul_excluded_timepoints)))
})

# test_that("nomatch simul z_star matches gold standard", {
#     fit2 <- add_simultaneous_ci(fit_nomatch, seed = 1)
#     expect_snapshot(fit2$simul_z_star)
# })
# 
# test_that("matching simul z_star matches gold standard", {
#     fit2 <- add_simultaneous_ci(fit_matching, seed = 1)
#     expect_snapshot(fit2$simul_z_star)
# })



# Failure modes (input validation)  -----------------------------------------
test_that("add_simultaneous_ci() errors on non-nomatchfit", {
    expect_error(add_simultaneous_ci(list()), "nomatchfit", ignore.case = TRUE)
})

test_that("add_simultaneous_ci() errors when only one timepoint", {
    fit_single <- suppressWarnings(nomatch(simdata, "Y","event","V","D_obs", c("x1","x2"),
                     immune_lag = 14, timepoints = 50, boot_reps = 5))
    expect_error(add_simultaneous_ci(fit_single),
                 "more than 1 timepoint|>=\\s*2", ignore.case = TRUE)
})

test_that("add_simultaneous_ci() errors when no or too few boot samples", {
    fit_no_boot <- suppressWarnings(nomatch(simdata, "Y","event","V","D_obs", c("x1","x2"),
                             immune_lag = 14, timepoints = seq(30,180,30), boot_reps = 0))
    expect_error(add_simultaneous_ci(fit_no_boot),
                 "bootstrap|boot", ignore.case = TRUE)
    
    fit_no_boot <- suppressWarnings(nomatch(simdata, "Y","event","V","D_obs", c("x1","x2"),
                                            immune_lag = 14, timepoints = seq(30,180,30), boot_reps = 1))
    expect_error(add_simultaneous_ci(fit_no_boot),
                 "bootstrap|boot", ignore.case = TRUE)
})


test_that("add_simultaneous_ci() errors when  boot samples not returned", {
    fit_no_boot <- suppressWarnings(nomatch(simdata, "Y","event","V","D_obs", c("x1","x2"),
                                            immune_lag = 14, timepoints = seq(30,180,30), boot_reps = 10,
                                            keep_boot_samples = FALSE))
    expect_error(add_simultaneous_ci(fit_no_boot),
                 "bootstrap|boot", ignore.case = TRUE)
})

# --- Reproducibility ---
test_that("seed controls z_star reproducibility", {
    z1 <- add_simultaneous_ci(fit_nomatch, seed = 42)$simul_z_star
    z2 <- add_simultaneous_ci(fit_nomatch, seed = 42)$simul_z_star
    z3 <- add_simultaneous_ci(fit_nomatch, seed = 99)$simul_z_star
    
    expect_equal(z1, z2)
    expect_false(identical(z1, z3))
})


# Edge Cases --------------------------------------------------------------
test_that("add_simultaneous_ci() succeeds with percentile confidence intervals", {
    fit_pct <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1", "x2"),
                       immune_lag = 14, timepoints = seq(30, 180, 30),
                       ci_type = "percentile", boot_reps = 5)
    expect_null(fit_pct$boot_samples$transformed)  # confirm precondition
    expect_no_error(add_simultaneous_ci(fit_pct, seed = 1))
})


