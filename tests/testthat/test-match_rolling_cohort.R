# Helpers -----------------------------------------------------------------
make_matchable_data <- function() {
    data.frame(
        ID        = 1:10,
        Y         = c(100, 100, 100, 100, 100, 50, 50, 50, 50, 50),
        event     = rep(0, 10),
        V         = c(1,   1,   1,   1,   1,  0,  0,  0,  0,  0),
        D_obs     = c(10,  10,  20,  20,  30, NA, NA, NA, NA, NA),
        x1        = c(1,   1,   2,   2,   1,  1,  1,  2,  2,  1)
    )
}
make_unmatchable_data <- function() {
    data.frame(
        ID    = 1:4,
        Y     = c(100, 100,  5,  5),   # controls have outcome before any exposure
        event = c(0, 0, 1, 1),
        V     = c(1, 1, 0, 0),
        D_obs = c(10, 10, NA, NA),
        x1    = c(1, 1, 1, 1)
    )
}


# Validate inputs ---------------------------------------------------------
test_that("passes with valid input", {
    dat <- make_matchable_data()
    expect_no_error(
        match_rolling_cohort(dat, "Y", "V", "D_obs", "x1", "ID")
    )
})

test_that("errors when a required column is missing", {
    dat <- make_matchable_data()[, -which(names(make_matchable_data()) == "x1")]
    expect_error(
        match_rolling_cohort(dat, "Y", "V", "D_obs", "x1", "ID"),
        "missing.*covariates|x1", ignore.case = TRUE
    )
})

test_that("errors when id_name is not unique", {
    dat <- make_matchable_data()
    dat$ID[2] <- dat$ID[1]   # duplicate
    expect_error(
        match_rolling_cohort(dat, "Y", "V", "D_obs", "x1", "ID"),
        "unique identifier", ignore.case = TRUE
    )
})

test_that("errors when reserved variable name already exists in data", {
    dat <- make_matchable_data()
    dat$match_id <- 999
    dat$match_index_time <- 999
    dat$match_type <- 999
    dat$match_V <- 999
    expect_error(
        match_rolling_cohort(dat, "Y", "V", "D_obs", "x1", "ID"),
        "conflict with internal", ignore.case = TRUE
    )
})



# Check output ------------------------------------------------------------
test_that("returns a list with the three expected components", {
    dat <- make_matchable_data()
    out <- match_rolling_cohort(dat, "Y", "V", "D_obs", "x1", "ID")
    expect_named(out, c("matched_data", "n_unmatched_cases", "discarded"))
})

test_that("matched_data has exactly 2 rows per match_id", {
    dat <- make_matchable_data()
    out <- match_rolling_cohort(dat, "Y", "V", "D_obs", "x1", "ID")
    rows_per_pair <- table(out$matched_data$match_id)
    expect_true(all(rows_per_pair == 2))
})

test_that("required new columns are present in matched_data", {
    dat        <- make_matchable_data()
    out        <- match_rolling_cohort(dat, "Y", "V", "D_obs", "x1", "ID")
    md         <- out$matched_data
    expect_true(all(c("match_index_time", "match_type", "match_V", "match_id") %in% names(md)))
})

test_that("discarded is a logical vector with length == nrow(data)", {
    dat <- make_matchable_data()
    out <- match_rolling_cohort(dat, "Y", "V", "D_obs", "x1", "ID")
    expect_type(out$discarded, "logical")
    expect_length(out$discarded, nrow(dat))
})



# Matching logic ----------------------------------------------------------

test_that("controls are unexposed at their match_index_time", {
    out <- match_rolling_cohort(simdata, "Y", "V", "D_obs", c("x1", "x2"), "ID")
    md  <- out$matched_data
    
    controls <- md[md$match_type == "control", ]
    # A control must be unvaccinated OR vaccinated after the index date
    unexposed_at_index <- is.na(controls$D_obs) | controls$D_obs > controls$match_index_time
    expect_true(all(unexposed_at_index))
})

test_that("controls are event-free at their match_index_time", {
    out <- match_rolling_cohort(simdata, "Y", "V", "D_obs", c("x1", "x2"), "ID")
    md  <- out$matched_data
    
    controls <- md[md$match_type == "control", ]
    expect_true(all(controls$Y > controls$match_index_time))
})

test_that("exact matching on matching_vars holds within each pair", {
    out <- match_rolling_cohort(simdata, "Y", "V", "D_obs", c("x1", "x2"), "ID")
    md  <- out$matched_data
    
    cases    <- md[md$match_type == "case",    ]
    controls <- md[md$match_type == "control", ]
    
    # order both by match_id so rows align
    cases    <- cases[order(cases$match_id), ]
    controls <- controls[order(controls$match_id), ]
    
    expect_identical(cases$x1, controls$x1)
    expect_identical(cases$x2, controls$x2)
})

test_that("with replace = FALSE, each control ID appears at most once", {
    out <- match_rolling_cohort(simdata, "Y", "V", "D_obs", c("x1", "x2"), "ID",
                                replace = FALSE)
    controls <- out$matched_data[out$matched_data$match_type == "control", ]
    expect_false(any(duplicated(controls$ID)))
})

test_that("with replace = TRUE, control IDs may repeat", {
    # Use a dataset with many cases but few controls to force reuse
    set.seed(1)
    n  <- 200
    dat <- data.frame(
        ID    = 1:n,
        Y     = rep(200, n),
        event = rep(0, n),
        V     = c(rep(1, 190), rep(0, 10)),
        D_obs = c( rep(10:28, each = 10), rep(NA, 10)),
        x1    = rep(1, n)
    )
    out <- match_rolling_cohort(dat, "Y", "V", "D_obs", "x1", "ID", replace = TRUE)
    controls <- out$matched_data[out$matched_data$match_type == "control", ]
    expect_true(any(duplicated(controls$ID)))
})





# Edge cases --------------------------------------------------------------
test_that("warns and returns empty matched_data when no matches possible", {
    dat <- make_unmatchable_data()
    expect_warning(
        out <- match_rolling_cohort(dat, "Y", "V", "D_obs", "x1", "ID"),
        "No matches found", ignore.case = TRUE
    )
    expect_equal(nrow(out$matched_data), 0)
})






# Reproducability ---------------------------------------------------------
test_that("seed produces reproducible results", {
    out1 <- match_rolling_cohort(simdata, "Y", "V", "D_obs", c("x1", "x2"), "ID", seed = 42)
    out2 <- match_rolling_cohort(simdata, "Y", "V", "D_obs", c("x1", "x2"), "ID", seed = 42)
    expect_identical(out1$matched_data, out2$matched_data)
})

test_that("different seeds can produce different results", {
    out1 <- match_rolling_cohort(simdata, "Y", "V", "D_obs", c("x1", "x2"), "ID", seed = 1)
    out2 <- match_rolling_cohort(simdata, "Y", "V", "D_obs", c("x1", "x2"), "ID", seed = 2)
    # It's highly unlikely both produce identical control assignments on a large dataset
    expect_false(identical(out1$matched_data$ID, out2$matched_data$ID))
})