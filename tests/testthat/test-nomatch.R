# Correct output ----------------------------------------------------------
test_that("nomatch() results match gold standard", {
    fit <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                   immune_lag = 14, timepoints = seq(30,180, by = 30), 
                   ci_type = "both",
                   boot_reps = 5,
                   seed = 1234)
    gold_standard <- readRDS(test_path("fixtures", "gold_nomatch.rds"))
    expect_equal(fit$estimates, gold_standard$estimates)
})



test_that("nomatch() returns a nomatchfit with correct structure", {
    fit <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                   immune_lag = 14, timepoints = seq(30, 180, 30), boot_reps = 0)
    expect_s3_class(fit, "nomatchfit")
    expect_setequal(names(fit$estimates),
                    c("cuminc_0", "cuminc_1", "risk_difference", "risk_ratio", "relative_risk_reduction"))
    
    fit <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                   immune_lag = 0, timepoints = seq(30, 180, 30), boot_reps = 0)
    expect_s3_class(fit, "nomatchfit")
    expect_setequal(names(fit$estimates),
                    c("cuminc_0", "cuminc_1", "risk_difference", "risk_ratio", "relative_risk_reduction"))
})


test_that("nomatch() is reproducible under fixed seed", {
    fit1 <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                    immune_lag = 14, timepoints = seq(30,180,30),
                    boot_reps = 5, seed = 42)
    fit2 <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                    immune_lag = 14, timepoints = seq(30,180,30),
                    boot_reps = 5, seed = 42)
    expect_identical(fit1$estimates, fit2$estimates)
    
    
    future::plan(future::multisession, workers = 2)
    fit3 <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                    timepoints = seq(30,180,30), immune_lag = 14,
                    boot_reps = 5, seed = 42)
    future::plan(future::sequential)
    expect_identical(fit1$estimates, fit3$estimates)
})

test_that("nomatch() gives different results with different seeds", {
    fit1 <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                    immune_lag = 14, timepoints = seq(30,180,30),
                    boot_reps = 5, seed = 42)
    fit2 <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                    immune_lag = 14, timepoints = seq(30,180,30),
                    boot_reps = 5, seed = 99)
    expect_false(identical(fit1$estimates, fit2$estimates))
})


# Confidence intervals ----------------------------------------------------
test_that("nomatch() adds CI columns when boot_reps > 0", {
    nomatch_with_ci <- function(ci_type){
        nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                immune_lag = 14, timepoints = seq(30, 90, 30),
                boot_reps = 5, ci_type = ci_type)
    }
    fit <- nomatch_with_ci("wald")
    expect_true(all(c("wald_lower", "wald_upper") %in% 
                        colnames(fit$estimates$cuminc_0)))
    
    fit <- nomatch_with_ci("percentile")
    expect_true(all(c("percentile_lower", "percentile_upper") %in% 
                        colnames(fit$estimates$cuminc_0)))
    
    #ci_type = "boot"
    fit <- nomatch_with_ci("both")
    expect_true(all(c("wald_lower", "wald_upper", "percentile_lower", "percentile_upper") %in% 
                        colnames(fit$estimates$cuminc_0)))
    
})


# Defaults ----------------------------------------------------------------
test_that("nomatch() resolves timepoints = NULL without error", {
    fit <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                   immune_lag = 14, timepoints = NULL, boot_reps = 0)
    expect_true(length(fit$timepoints) > 0)
    expect_true(all(fit$timepoints > 14))
})



# Custom formulas ---------------------------------------------------------
test_that("nomatch() accepts custom formulas", {
    #formulas as strings without ~ 
    fit <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                   immune_lag = 14, timepoints = seq(30, 90, 30),
                   boot_reps = 0,
                   formula_unexposed = "x1 + x2",
                   formula_exposed = "x1 + x2 + splines::ns(D_obs, df = 3)")
    expect_s3_class(fit, "nomatchfit")
    
    #formulas as strings with ~ 
    fit <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                   immune_lag = 14, timepoints = seq(30, 90, 30),
                   boot_reps = 0,
                   formula_unexposed = "~ x1 + x2",
                   formula_exposed = "~ x1 + x2 + splines::ns(D_obs, df = 3)")
    expect_s3_class(fit, "nomatchfit")
    
    #formulas as vectors
    fit <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                   immune_lag = 14, timepoints = seq(30, 90, 30),
                   boot_reps = 0,
                   formula_unexposed = c("x1", "x2"),
                   formula_exposed = c("x1", "x2", "splines::ns(D_obs, df = 3)"))
    expect_s3_class(fit, "nomatchfit")
    
    #actual formulas
    fit <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                   immune_lag = 14, timepoints = seq(30, 90, 30),
                   boot_reps = 0,
                   formula_unexposed = ~ x1 + x2,
                   formula_exposed = ~ x1 + x2 + splines::ns(D_obs, df = 3))
    expect_s3_class(fit, "nomatchfit")
})


# Custom weights ----------------------------------------------------------
test_that("nomatch() accepts custom weights", {
    w <- get_observed_weights(simdata, "Y", "V", "D_obs", c("x1","x2"), immune_lag = 14)
    fit <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1","x2"),
                   immune_lag = 14, timepoints = seq(30, 90, 30),
                   weights_source = "custom", custom_weights = w,
                   boot_reps = 0)
    expect_s3_class(fit, "nomatchfit")
})


# Edge cases --------------------------------------------------------------
test_that("nomatch() works with covariates = NULL", {
    fit <- nomatch(simdata, "Y", "event", "V", "D_obs",
                   covariates = NULL, immune_lag = 14,
                   timepoints = seq(30, 90, 30), boot_reps = 0)
    expect_s3_class(fit, "nomatchfit")
})


test_that("individuals with exposure_time > outcome_time are treated as unexposed", {
    # Setup: take the simulated data and corrupt ~50% of exposed individuals
    # by shifting their exposure time to AFTER their event time
    dat <- simdata
    set.seed(123)
    late_exposure <- rbinom(n = sum(dat$V == 1), size = 1, prob = 0.5)
    dat$D_obs[dat$V == 1] <- ifelse(late_exposure,  
                                    dat$D_obs[dat$V == 1] + dat$Y[dat$V == 1],  # push exposure past event
                                    dat$D_obs[dat$V == 1])
    tau <- 14
    fit <- nomatch(dat, "Y", "event", "V", "D_obs",
                   covariates = c("x1"), immune_lag = tau,
                   timepoints = seq(30, 90, 30), boot_reps = 0)
    
    # Verify that model_1 and weights only include individuals 
    # who were exposed before their event 
    n_valid_exposed <- sum(dat$V == 1 & dat$Y > dat$D_obs + tau)
    expect_equal(n_valid_exposed, fit$model_1$n)
    expect_equal(mean(dat$x1[dat$V == 1 & dat$Y > dat$D_obs + tau]),
                 fit$weights$p_weights$prob[fit$weights$p_weights$x1 == 1])
    
    # Reclassifying individuals with exposure_time > outcome time as unexposed
    # should produce identical results since only interested in exposures
    # that happen prior to first event 
    dat2 <- dat
    dat2$V[dat$V == 1]    <- ifelse(late_exposure, 0, 1)          # recode as unexposed
    dat2$D_obs[dat$V == 1] <- ifelse(late_exposure, NA, dat2$D_obs[dat$V == 1])  # clear their exposure time
    fit2 <- nomatch(dat2, "Y", "event", "V", "D_obs",
                    covariates = c("x1"), immune_lag = tau,
                    timepoints = seq(30, 90, 30), boot_reps = 0)
    
    fit$call <- fit2$call <- NULL
    expect_equal(fit, fit2)
})


test_that("nomatch() works with timepoint = 0 but not timepoint < 0", {
    expect_no_error(
        suppressWarnings(
        fit <- nomatch(simdata, "Y", "event", "V", "D_obs", c("x1", "x2"),
            immune_lag = 14, timepoints = 0:20)
        )
    )
    df <- estimates_to_df(fit)
    at_tau <- df[df$t0 <= 14 & df$term %in% c("cuminc_0", "cuminc_1"), ]
    expect_true(all(at_tau$estimate == 0))
    
    
    expect_error(
        nomatch(simdata, "Y", "event", "V", "D_obs", c("x1", "x2"),
                   immune_lag = 14, timepoints = -5:20)
    )
})






# Key failure modes -------------------------------------------------------
test_that("nomatch() fails loudly on bad inputs", {
    # Inconsistent exposure status and exposure time
    dat <- simdata
    dat$V[is.na(dat$D_obs)][1] <- 1
    expect_error(
        nomatch(dat, "Y", "event", "V", "D_obs", c("x1", "x2"),
                immune_lag = 14, timepoints = 30),
        "non-missing", ignore.case = TRUE
    )
})

test_that("nomatch() fails for reserved vars", {
    dat <- simdata
    names(dat)[names(dat) %in% c("x1",  "Y", "event")]  <- c("event",  "outcome_time", "outcome_status")
    
    expect_error(
        nomatch(dat, "outcome_time", "outcome_status", "V", "D_obs", covariates = c("event", "x2"),
                immune_lag = 14, timepoints = 30),
        "conflict with internal", ignore.case = TRUE
    )
    
    dat <- simdata
    names(dat)[names(dat) %in% c("x2",  "Y", "event")]  <- c("T1",  "outcome_time", "outcome_status")
    
    expect_error(
        nomatch(dat, "outcome_time", "outcome_status", "V", "D_obs", covariates = c("x1", "T1"),
                immune_lag = 14, timepoints = 30),
        "conflict with internal", ignore.case = TRUE
    )
})

test_that("nomatch() fails when no events", {
    dat <- simdata
    dat$event <- 0
    expect_error(
        nomatch(dat, "Y", "event", "V", "D_obs", c("x1", "x2"),
                immune_lag = 14, timepoints = 30),
        "no events", ignore.case = TRUE
    )
    
    dat <- simdata
    dat$event[dat$V == 1] <- 0
    expect_error(
        nomatch(dat, "Y", "event", "V", "D_obs", c("x1", "x2"),
                immune_lag = 14, timepoints = 30),
        "no events", ignore.case = TRUE
    )
})

test_that("nomatch() fails when no individuals to fit the exposed model ", {
    dat <- simdata
    dat$Y[dat$V == 1] <- dat$D_obs[dat$V==1] + 1
    expect_error(
        nomatch(dat, "Y", "event", "V", "D_obs", c("x1", "x2"),
                immune_lag = 14, timepoints = 30),
        "no eligible individuals", ignore.case = TRUE
    )
})




