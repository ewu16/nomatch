# Shared fixture — fit once, reuse across all tests
make_fit <- function(...) {
    args <- list(
        data         = simdata,
        outcome_time = "Y",
        outcome_status = "event",
        exposure     = "V",
        exposure_time = "D_obs",
        covariates   = c("x1", "x2"),
        immune_lag   = 14,
        timepoints   = seq(30, 180, 30)
    )
    do.call(nomatch, modifyList(args, list(...)))
}

fit <- make_fit() 


# ── Input validation ──────────────────────────────────────────────────────────

test_that("summarize_nomatch_population() errors on bad `fit` argument", {
    expect_error(summarize_nomatch_population(simdata, fit = list()),
                 "nomatchfit", ignore.case = TRUE)
})

test_that("summarize_nomatch_population() errors when fit is from matching()", {
    fit_m <- fit
    fit_m$method <- "matching"
    expect_error(summarize_nomatch_population(simdata, fit_m),
                 "nomatch", ignore.case = TRUE)
})

test_that("summarize_nomatch_population() errors when models were not kept", {
    fit_no_models <- make_fit(keep_models = FALSE)
    expect_error(summarize_nomatch_population(simdata, fit_no_models),
                 "keep_models", ignore.case = TRUE)
})

test_that("summarize_nomatch_population() errors on empty data", {
    expect_error(summarize_nomatch_population(simdata[0, ], fit),
                 "non-empty", ignore.case = TRUE)
})


# ── Output structure ──────────────────────────────────────────────────────────

test_that("summarize_nomatch_population() matches gold standard", {
    out <- summarize_nomatch_population(simdata, fit)
    
    gold <- data.frame(
        variable = c("x1", "x2"),
        label    = c("", ""),
        mean     = c(0.5035969, 8.0659550)
    )
    expect_equal(out, gold)
})



# ── Numerical correctness ─────────────────────────────────────────────────────

test_that("summarize_nomatch_population() factor level proportions sum to 1 per variable", {
    # Use a version of simdata with a factor covariate to test this path
    dat <- simdata
    dat$x1 <- as.character(dat$x1)  # force non-numeric
    fit_factor <- local(make_fit(data = dat))
    out <- summarize_nomatch_population(dat, fit_factor)
    
    x1_rows <- out[out$variable == "x1", ]
    expect_equal(sum(x1_rows$mean), 1, tolerance = 1e-6)
})

