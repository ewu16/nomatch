# Tests for validate_data_inputs() -----------------------------------------------
test_that("validate_data_inputs() passes on clean input", {
    core_args = list(outcome_time = "Y",
                     outcome_status = "event",
                     exposure = "V",
                     exposure_time = "D_obs")
    
    expect_no_error(validate_data_inputs(simdata, core_args))
    expect_no_error(validate_data_inputs(simdata, core_args, covariates = "x1"))
    expect_no_error(validate_data_inputs(simdata, core_args, covariates = c("x1", "x2")))
    
    core_args_with_null <- core_args
    core_args_with_null$outcome_status <- NULL
    expect_no_error(validate_data_inputs(simdata, core_args_with_null))
})

test_that("validate_data_inputs() catches required args", {
    core_args = list(outcome_time = "Y",
                     outcome_status = NULL,
                     exposure = "V",
                     exposure_time = "D_obs")
    
    expect_error(validate_data_inputs(simdata, core_args), "non-null")
})

test_that("validate_data_inputs() checks columns", {
    core_args  = list(outcome_time = "Y",
                              outcome_status = "event",
                              exposure = "V",
                              exposure_time = "D_obs")
    
    core_args_not_present <- core_args; core_args_not_present$outcome_time <- "not_here"
    expect_error(validate_data_inputs(simdata, core_args_not_present), "not_here")
    
    data <- simdata; data$Y <- -data$Y
    expect_error(validate_data_inputs(data, core_args), "non-negative")
    
    data <- simdata; data$event <- data$Y
    expect_error(validate_data_inputs(data, core_args), "0/1")
    
    data <- simdata; data$x1 <- data$D_obs
    expect_error(validate_data_inputs(data, core_args, covariates = "x1"), "Missing values")
})


# Tests for validate_data_str() -------------------------------------------
test_that("validate_data_str() passes on clean input", {
    expect_no_error(validate_data_str(simdata))
})

test_that("validate_data_str() fails on non-data.frame and empty data", {
    expect_error(validate_data_str(NULL), "must be a data.frame")
    expect_error(validate_data_str(as.matrix(simdata)), "must be a data.frame")
    expect_error(validate_data_str(simdata[0,]), "no rows")
})

# Tests for validate_column_args() ----------------------------------------
test_that("validate_column_args() passes on clean input", {
    expect_no_error(validate_column_args(list(a = "a", b = "b")))
    expect_no_error(validate_column_args(list(a = "a", b = "b"), covariates = "c"))
    expect_no_error(validate_column_args(list(a = "a", b = "b"), covariates = c("c", "d")))
})

test_that("validate_column_args() fails on unnamed list", {
    expect_error(validate_column_args(list("a",  "b")), "named list")
})

test_that("validate_column_args() fails when core_args is not a list of single, non-null character strings", {
    expect_error(validate_column_args(list(a = 1, b = 2)), "single, non-null, character strings: a, b")
    expect_error(validate_column_args(list(a = "a", b = 1)), "single, non-null, character strings: b")
    expect_error(validate_column_args(list(a = "a", b = NULL)), "single, non-null, character strings: b")
    expect_error(validate_column_args(list(a = "a", b = TRUE)), "single, non-null, character strings: b")
    expect_error(validate_column_args(list(a = "a", b = c("b", "c"))), "single, non-null, character strings: b")
})

test_that("validate_column_args() fails when covariates is not a character vector", {
    expect_error(validate_column_args(list(a = "a", b = "b"), covariates = 1), "`covariates`")
    expect_error(validate_column_args(list(a = "a", b = "b"), covariates = 1:2), "`covariates`")
    expect_error(validate_column_args(list(a = "a", b = "b"), covariates = TRUE), "`covariates`")
})


# Tests for validate_cols_exist -------------------------------------------
test_that("validate_cols_exist() passes when all variable present", {
    core_args <-  list(outcome_time = "Y", outcome_status = "event")
    expect_no_error(validate_cols_exist(simdata, core_args))
    expect_no_error(validate_cols_exist(simdata, core_args, covariates = c("x1")))
    expect_no_error(validate_cols_exist(simdata, core_args, covariates = c("x1", "x2")))
})

test_that("validate_cols_exist() errors listing missing columns", {
    expect_error(
        validate_cols_exist(simdata, list(outcome_time = "not_here", outcome_status = "event")), 
        "Missing `outcome_time`"
    )
    expect_error(
        validate_cols_exist(simdata, list(outcome_time = "Y", outcome_status = "not_here")), 
        "Missing `outcome_status`"
    )
    
    expect_error(
        validate_cols_exist(simdata, list(outcome_time = "Y", outcome_status = "event"), covariates = c("x1", "not_here")), 
        "Missing `covariates` .* not_here"
    )
})


# Tests for column-level checks -------------------------------------------
test_that("validate_outcome_time() catches NA, non-numerics, and negatives", {
    expect_no_error(validate_outcome_time(simdata, "Y"))
    
    dat <- simdata; dat$Y[1] <- NA
    expect_error(validate_outcome_time(dat, "Y"), "Missing values.*outcome_time", ignore.case = TRUE)
    
    dat <- simdata; dat$Y <- as.character(dat$Y)
    expect_error(validate_outcome_time(dat, "Y"), "must be numeric", ignore.case = TRUE)
    
    dat <- simdata; dat$Y[2] <- -1
    expect_error(validate_outcome_time(dat, "Y"), "must be non-negative", ignore.case = TRUE)

})


test_that("validate_outcome_status() catches NA and non-binary coding", {
    expect_no_error(validate_outcome_status(simdata, "event"))
    
    dat <- simdata; dat$event[1] <- NA
    expect_error(validate_outcome_status(dat, "event"), "Missing values.*outcome_status", ignore.case = TRUE)
    
    dat <- simdata; dat$event[2] <- 2
    expect_error(validate_outcome_status(dat, "event"), "0/1", ignore.case = TRUE)
})


test_that("validate_exposure_args() checks exposure coding, time sign, and consistency", {
    expect_no_error(validate_exposure_args(simdata, "V", "D_obs"))
    
    # NA in exposure
    dat <- simdata; dat$V[1] <- NA
    expect_error(validate_exposure_args(dat, "V", "D_obs"), "Missing values in <exposure>", fixed = TRUE)
    
    # non-binary exposure
    dat <- simdata; dat$V[1] <- 2
    expect_error(validate_exposure_args(dat, "V", "D_obs"), "0/1", ignore.case = TRUE)
    
    # non-numeric exposure time
    dat <- simdata; dat$D_obs <- as.character(dat$D_obs)
    expect_error(validate_exposure_args(dat, "V", "D_obs"), "must be numeric", ignore.case = TRUE)
    
    
    # negative exposure time
    dat <- simdata; dat$D_obs[2] <- -3
    expect_error(validate_exposure_args(dat, "V", "D_obs"), "must be non-negative", ignore.case = TRUE)
    
    # missing time for exposed
    dat <- simdata; dat$D_obs[dat$V == 1][1] <- NA
    expect_error(validate_exposure_args(dat, "V", "D_obs"), "All exposed individuals must have non-missing", ignore.case = TRUE)
    
    # non-NA time for unexposed
    dat <- simdata; dat$D_obs[dat$V == 0][1] <- 9
    expect_error(validate_exposure_args(dat, "V", "D_obs"), "All unexposed individuals must have .* NA", ignore.case = TRUE)
})


test_that("validate_covariates() errors on missingness in any covariate", {
    expect_no_error(validate_covariates(simdata, NULL))
    expect_no_error(validate_covariates(simdata, c("x1", "x2")))
    
    dat <- simdata; dat$x1[1] <- NA
    expect_error(
        validate_covariates(dat, c("x1", "x2")),
        "Missing values in covariates", ignore.case = TRUE
    )
    
    dat <- simdata; dat$x1[1] <- NA; dat$x2[2] <- NA
    expect_error(
        validate_covariates(dat, c("x1", "x2")),
        "Missing values in covariates", ignore.case = TRUE
    )
})
