test_that("validate_data() passes on clean input", {
    expect_no_error(validate_data(simdata))
})

test_that("validate_data() fails on non-data.frame and empty data", {
    expect_error(validate_data(NULL), "must be a data.frame")
    expect_error(validate_data(as.matrix(simdata)), "must be a data.frame")
    expect_error(validate_data(simdata[0,]), "no rows")
})


test_that("validate_args_not_null() requires names and flags NULLs", {
    expect_no_error(validate_args_not_null(list(a = 1, b = 2)))
    expect_error(validate_args_not_null(list(1, b = 2)), "must be named", ignore.case = TRUE)
    expect_error(
        validate_args_not_null(list(a = 1, b = NULL, c = NULL)),
        "NULL values are not allowed.*b, c", ignore.case = TRUE
    )
})


test_that("validate_cols_exist() passes when all present", {
    expect_no_error(validate_cols_exist(simdata, c("Y", "event", "V")))
})

test_that("validate_cols_exist() errors listing missing columns", {
    expect_error(
        validate_cols_exist(simdata, c("Y", "event", "not_here")),
        "Missing required column\\(s\\).*not_here", ignore.case = TRUE
    )
})


test_that("validate_outcome_time() catches NA and negatives", {
    expect_no_error(validate_outcome_time(simdata, "Y"))

    dat <- simdata; dat$Y[1] <- NA
    expect_error(validate_outcome_time(dat, "Y"), "Missing values.*outcome_time", ignore.case = TRUE)

    dat <- simdata; dat$Y[2] <- -1
    expect_error(validate_outcome_time(dat, "Y"), "must be non-negative", ignore.case = TRUE)
})


test_that("validate_outcome_status() enforces 0/1 and no NA", {
    expect_no_error(validate_outcome_status(simdata, "event"))

    dat <- simdata; dat$event[1] <- NA
    expect_error(validate_outcome_status(dat, "event"), "Missing values.*outcome_status", ignore.case = TRUE)

    dat <- simdata; dat$event[2] <- 2
    expect_error(validate_outcome_status(dat, "event"), "coded 0/1", ignore.case = TRUE)
})


test_that("validate_exposure_args() checks exposure coding, time sign, and consistency", {
    expect_no_error(validate_exposure_args(simdata, "V", "D_obs"))

    # NA in exposure
    dat <- simdata; dat$V[1] <- NA
    expect_error(validate_exposure_args(dat, "V", "D_obs"), "Missing values in <exposure>", fixed = TRUE)

    # non-binary exposure
    dat <- simdata; dat$V[1] <- 2
    expect_error(validate_exposure_args(dat, "V", "D_obs"), "must be coded 0/1", ignore.case = TRUE)

    # negative exposure time
    dat <- simdata; dat$D_obs[2] <- -3
    expect_error(validate_exposure_args(dat, "V", "D_obs"), "must be non-negative", ignore.case = TRUE)

    # missing time for exposed
    dat <- simdata; dat$D_obs[dat$V == 1][1] <- NA
    expect_error(validate_exposure_args(dat, "V", "D_obs"), "All exposed individuals must have non-missing", ignore.case = TRUE)

    # non-NA time for unexposed
    dat <- simdata; dat$D_obs[dat$V == 0][1] <- 9
    expect_error(validate_exposure_args(dat, "V", "D_obs"), "All unexposed individuals must have NA", ignore.case = TRUE)
})

test_that("validate_covariates() catches NA", {
    expect_no_error(validate_covariates(simdata, NULL))
    expect_no_error(validate_covariates(simdata, c("x1", "x2")))
})

test_that("validate_covariates() errors on missingness in any covariate", {
    expect_no_error(validate_covariates(simdata, c("x1", "x2")))

    dat <- simdata; dat$x1[1] <- NA
    expect_error(
        validate_covariates(dat, c("x1", "x2")),
        "Missing values in covariates", ignore.case = TRUE
    )
})


test_that("validate_time_args() enforces types, ranges, relationships, and warns on long lag", {
    # happy path
    expect_no_error(validate_time_args(immune_lag = 14, timepoints = c(30, 60)))

    # immune_lag issues
    expect_error(validate_time_args(immune_lag = NA_real_, timepoints = 30), "single non-missing", ignore.case = TRUE)
    expect_error(validate_time_args(immune_lag = -1, timepoints = 30), "must be >= 0", ignore.case = TRUE)

    # timepoints issues
    expect_error(validate_time_args(14, timepoints = numeric(0)), "length > 0", ignore.case = TRUE)
    expect_error(validate_time_args(14, timepoints = c(NA, 60)), "no missing values", ignore.case = TRUE)
    expect_error(validate_time_args(14, timepoints = c(7, 14)), "must be > `immune_lag`", ignore.case = TRUE)

    # warning on long lag
    expect_warning(validate_time_args(immune_lag = 29, timepoints = c(60, 90)), "unusually long", ignore.case = TRUE)
})



test_that("check_reserved_vars() passes when there are no conflicts", {
    expect_no_error(check_reserved_vars(vars = c("x1", "x2"), reserved_vars = c(".t0", ".estimate"), vars_label = "covariates"))
})

test_that("check_reserved_vars() errors when vars collide with reserved names", {
    reserved <- c(".t0", ".estimate", "group_id")
    expect_error(
        check_reserved_vars(vars = c("x1", "group_id"), reserved_vars = reserved, vars_label = "covariates"),
        "contain reserved variable names.*group_id", ignore.case = TRUE
    )
})

