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

