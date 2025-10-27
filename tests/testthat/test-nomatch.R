test_that("nomatch() fails gracefully when covariates missing", {
    bad <- simdata[, -which(names(simdata) == "x2")]
    expect_error(
        nomatch(bad, "Y", "event", "V", "D_obs", c("x1","x2"), immune_lag = 14, eval_times = c(30,60)),
        "Missing required column", ignore.case = TRUE
    )
})

# tests/test-nomatch-basic.R
test_that("basic nomatch run works", {
    set.seed(123)
    fit <- nomatch(simdata, "Y","event","V","D_obs", c("x1","x2"),
                     immune_lag = 14, eval_times = seq(30,180,30), boot_reps = 5)
    expect_s3_class(fit, "vefit")
    expect_no_error(summary(fit))
    p <- plot(fit)
    expect_s3_class(p, "ggplot")

    df <- estimates_to_df(fit)
    expect_true(all(c("term","t0","estimate") %in% names(df)))
})

test_that("nomatch with risk_difference works", {
    set.seed(123)
    fit <- nomatch(simdata, "Y","event","V","D_obs", c("x1","x2"),
                   immune_lag = 14, eval_times = seq(30,180,30),
                   effect = "risk_difference",
                   boot_reps = 5)
    expect_s3_class(fit, "vefit")
    expect_no_error(print(fit))
    expect_no_error(summary(fit))
    p <- plot(fit)
    expect_s3_class(p, "ggplot")

    df <- estimates_to_df(fit)
    expect_true(all(c("term","t0","estimate") %in% names(df)))
})

test_that("nomatch with vaccine_effectiveness works", {
    set.seed(123)
    fit <- nomatch(simdata, "Y","event","V","D_obs", c("x1","x2"),
                   immune_lag = 14, eval_times = seq(30,180,30),
                   effect = "vaccine_effectiveness",
                   boot_reps = 5)
    expect_s3_class(fit, "vefit")
    expect_no_error(print(fit))
    expect_no_error(summary(fit))
    p <- plot(fit)
    expect_s3_class(p, "ggplot")

    df <- estimates_to_df(fit)
    expect_true(all(c("term","t0","estimate") %in% names(df)))
})
