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
    # matched_data_gold <- readRDS(test_path("fixtures", "gold_matched_data.rds"))
    # matched_data <- match_rolling_cohort(simdata, "Y", "V", "D_obs", 
    #                                      c("x1","x2"), "ID", seed = 123)$matched_data
    fit <-  matching(matched_data, "Y", "event", "V", "D_obs",
                     immune_lag = 14, timepoints = seq(30, 180, 30), 
                     ci_type = "both", boot_reps = 5, seed = 1234)
    gold_standard <- readRDS(test_path("fixtures", "gold_matching.rds"))
    expect_equal(fit$estimates, gold_standard$estimates)
})

test_that("matching() is reproducible under fixed seed", {
    fit1 <- matching(matched_data, "Y", "event", "V", "D_obs",
                     immune_lag = 14, timepoints = seq(30, 180, 30),
                     boot_reps = 5, seed = 42)
    fit2 <- matching(matched_data, "Y", "event", "V", "D_obs",
                     immune_lag = 14, timepoints = seq(30, 180, 30),
                     boot_reps = 5, seed = 42)
    expect_identical(fit1$estimates, fit2$estimates)
    
    future::plan(future::multisession, workers = 2)
    fit3 <- matching(matched_data, "Y", "event", "V", "D_obs",
                     immune_lag = 14, timepoints = seq(30, 180, 30),
                     boot_reps = 5, seed = 42)
    future::plan(future::sequential)
    expect_identical(fit1$estimates, fit3$estimates)
})

test_that("matching() gives different results with different seeds", {
    fit1 <- matching(matched_data, "Y", "event", "V", "D_obs",
                     immune_lag = 14, timepoints = seq(30, 180, 30),
                     boot_reps = 5, seed = 42)
    fit2 <- matching(matched_data, "Y", "event", "V", "D_obs",
                     immune_lag = 14, timepoints = seq(30, 180, 30),
                     boot_reps = 5, seed = 99)
    expect_false(identical(fit1$estimates, fit2$estimates))
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


# Edge cases --------------------------------------------------------------
test_that("individuals with exposure_time > outcome_time are treated as unexposed", {
    # Setup: take the simulated data and corrupt ~50% of exposed individuals
    # by shifting their exposure time to AFTER their event time
    dat <- simdata
    set.seed(123)
    late_exposure <- rbinom(n = sum(dat$V == 1), size = 1, prob = 0.5)
    dat$D_obs[dat$V == 1] <- ifelse(late_exposure,  
                                    dat$D_obs[dat$V == 1] + dat$Y[dat$V == 1],  # push exposure past event
                                    dat$D_obs[dat$V == 1])
    
    out <- match_rolling_cohort(dat, "Y", "V", "D_obs", c("x1", "x2"), "ID", seed = 123)
    
    
    # Reclassifying individuals with exposure_time > outcome time as unexposed
    # should produce identical results since only interested in exposures
    # that happen prior to first event 
    dat2 <- dat
    dat2$V[dat$V == 1]    <- ifelse(late_exposure, 0, 1)          # recode as unexposed
    dat2$D_obs[dat$V == 1] <- ifelse(late_exposure, NA, dat2$D_obs[dat$V == 1])  # clear their exposure time
    out2 <- match_rolling_cohort(dat2, "Y", "V", "D_obs", c("x1", "x2"), "ID", seed = 123)
    
    vars <- c("ID", "Y", "event", "match_index_time", "match_type", "match_V", "match_id")
    expect_equal(out[[1]][,vars],
                 out2[[1]][, vars])
    
    
    #Downstream matching analysis also handles individuals with exposure_time >
    # outcome time as unexposed
    fit <- matching(out[[1]], "Y", "event", "V", "D_obs",
                    immune_lag = 14, timepoints = seq(30,90, by = 30), boot_reps = 0)
    fit2 <- matching(out2[[1]], "Y", "event", "V", "D_obs",
                     immune_lag = 14, timepoints = seq(30,90, by = 30), boot_reps = 0)
    
    fit$call <- fit2$call <- NULL
    fit$matched_adata <- fit2$matched_adata <- NULL
    expect_equal(fit, fit2)
    
})

test_that("matching() works with timepoint = 0 but not timepoint < 0", {
    expect_no_error(
        suppressWarnings(
        fit <- matching(matched_data, "Y", "event", "V", "D_obs",
                 immune_lag = 14, timepoints = 0:20)
        )
    )
    
    df <- estimates_to_df(fit)
    at_tau <- df[df$t0 <= 14 & df$term %in% c("cuminc_0", "cuminc_1"), ]
    expect_true(all(at_tau$estimate == 0))
    
    expect_error(
        fit <- matching(matched_data, "Y", "event", "V", "D_obs",
                        immune_lag = 14, timepoints = -5:20)
    )
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


test_that("matching() warns when no events", {
    dat <- matched_data 
    dat$event <- 0
    suppressWarnings(
    expect_warning(
        fit <- matching(dat, "Y", "event", "V", "D_obs", immune_lag = 14, timepoints = c(30,60)),
        "No events", ignore.case = TRUE
    )
    )
})
