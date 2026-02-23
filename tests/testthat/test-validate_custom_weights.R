make_valid_weights <- function() {
    get_observed_weights(simdata, "Y", "V", "D_obs", c("x1", "x2"), immune_lag = 14)
}


test_that("validate_custom_weights() passes on clean input", {
    expect_no_error(
        validate_marginalizing_weights(make_valid_weights(), "D_obs", c("x1", "x2"))
    )
})

test_that("validate_custom_weights() passes on edge case (no covariates)", {
    w <-  get_observed_weights(simdata, "Y", "V", "D_obs", NULL, immune_lag = 14)
    expect_no_error(
        validate_marginalizing_weights(w, "D_obs", NULL)
    )
})


test_that("validate_custom_weights() errors for incorrect structure", {
    expect_error(
        validate_marginalizing_weights(NULL, "D_obs", c("x1", "x2")),
        "must be a list"
    )
    expect_error(
        validate_marginalizing_weights(list(a = 1, b =2), "D_obs", c("x1", "x2")),
        "must include components"
    )
    expect_error(
        validate_marginalizing_weights(list(g_weights = 1, p_weights =2), "D_obs", c("x1", "x2")),
        "must be data frames"
    )
})

test_that("validate_custom_weights() errors for missing columns", {
    w <- make_valid_weights()
    w$g_weights$D_obs <- NULL
    w$g_weights$prob <- NULL
    w$g_weights$x1 <- NULL
    
    expect_error(
        validate_marginalizing_weights(w, "D_obs", c("x1", "x2")),
        "missing required columns: D_obs, prob, x1"
    )
    
    w <- make_valid_weights()
    w$p_weights$prob <- NULL
    w$p_weights$x1 <- NULL
    
    expect_error(
        validate_marginalizing_weights(w, "D_obs", c("x1", "x2")),
        "missing required columns: prob, x1"
    )
})

test_that("validate_custom_weights() requires proper probabilties", {
    w <- make_valid_weights()
    w$g_weights$prob <- 2

    expect_error(
        validate_marginalizing_weights(w, "D_obs", c("x1", "x2")),
        "g_weights"
    )
    
    w <- make_valid_weights()
    w$p_weights$prob <- NA
    expect_error(
        validate_marginalizing_weights(w, "D_obs", c("x1", "x2")),
        "p_weights"
    )
})

test_that("validate_custom_weights() checks probabilities sum to 1", {
    w <-  get_observed_weights(simdata[1:10,], "Y", "V", "D_obs", "x1", immune_lag = 14)
    w$g_weights$prob[1] <- 0.2
    expect_error(
        validate_marginalizing_weights(w, "D_obs", "x1"),
        "must sum to 1 within each group"
    )
    
    w <-  get_observed_weights(simdata[1:10,], "Y", "V", "D_obs", "x1", immune_lag = 14)
    w$p_weights$prob[1] <- 0.2
    expect_error(
        validate_marginalizing_weights(w, "D_obs", "x1"),
        "must sum to 1"
    )
})



test_that("validate_prob_column() requires proper probabilities", {
    w <- make_valid_weights()
    prob <- w$g_weights$prob 
    expect_error( validate_prob_column(as.character(prob), "prob"), "numeric")
    
    prob <- w$g_weights$prob; prob[1] <- NA
    expect_error( validate_prob_column(prob, "prob"), "missing")
    
    prob <- w$g_weights$prob; prob[1] <- -1
    expect_error(validate_prob_column(prob, "prob"), "0 and 1")
    
    prob <- w$g_weights$prob; prob[1] <- 5
    expect_error(validate_prob_column(prob, "prob"), "0 and 1")
})