#' Validate custom marginalizing weights
#'
#' @description
#' Internal validation function that checks whether a user-supplied
#' list of marginalizing weights (`g_weights` and `p_weights`) has the required
#' structure and valid probability values.
#'
#' @param custom_weights List containing `g_weights` and `p_weights` data frames.
#' @param exposure_time Character; name of the vaccination time variable.
#' @param covariates Character vector of covariate names.
#'
#' @return Invisibly returns `NULL` if all checks pass. Throws an error with a
#' descriptive message otherwise.
#' @keywords internal
#' @noRd
validate_marginalizing_weights <- function(custom_weights, exposure_time, covariates) {
    # Check structure
    if (!is.list(custom_weights))
        stop("`custom_weights` must be a list containing 'g_weights' and 'p_weights'.", call. = FALSE)
    if (!all(c("g_weights", "p_weights") %in% names(custom_weights)))
        stop("List must include components 'g_weights' and 'p_weights'.", call. = FALSE)

    g <- custom_weights$g_weights
    p <- custom_weights$p_weights

    if (!is.data.frame(g) || !is.data.frame(p))
        stop("'g_weights' and 'p_weights' must be data frames.", call. = FALSE)

    # Required columns
    required_g <- c(exposure_time, "prob", covariates)
    required_p <- c("prob", covariates)

    check_missing_cols <- function(df, required, name) {
        missing <- setdiff(required, names(df))
        if (length(missing))
            stop(sprintf("`%s` missing required columns: %s",
                         name, paste(missing, collapse = ", ")), call. = FALSE)
    }

    check_missing_cols(g, required_g, "g_weights")
    check_missing_cols(p, required_p, "p_weights")

    # Validate probabilities
    validate_prob_column(g$prob, "g_weights")
    validate_prob_column(p$prob, "p_weights")

    # Check sums
    if (any(abs(tapply(g$prob, g[covariates], sum) - 1) > 1e-6))
        stop("In 'g_weights', probabilities must sum to 1 within each group.", call. = FALSE)

    if (abs(sum(p$prob) - 1) > 1e-6)
        stop("In 'p_weights', probabilities must sum to 1.", call. = FALSE)

    invisible(NULL)
}

#' Validate numeric probability vector
#' @keywords internal
#' @noRd
validate_prob_column <- function(prob, name) {
    if (!is.numeric(prob))
        stop(sprintf("'prob' column in '%s' must be numeric.", name), call. = FALSE)
    if (any(is.na(prob)))
        stop(sprintf("'prob' column in '%s' contains missing values.", name), call. = FALSE)
    if (any(prob < 0 | prob > 1))
        stop(sprintf("'prob' column in '%s' must be between 0 and 1.", name), call. = FALSE)
    invisible(NULL)
}

# g <- weights$g_weights
# g$prob[g$x2 == 5] <- 0
# p_weights <- weights$p_weights
# p_weights$prob[p_weights$x2 == 5] <- NA
# p <- p_weights


#' Convert weights into standard format for estimation procedure
#'
#' @keywords internal
#' @noRd
canonicalize_weights <- function(weights, exposure_time, covariates) {
    g <- weights$g_weights
    p <- weights$p_weights

    # Rename prob columns if needed
    if (!"prob_g" %in% names(g)) names(g)[names(g)=="prob"] <- "prob_g"
    if (!"prob_p" %in% names(p)) names(p)[names(p)=="prob"] <- "prob_p"

    # Drop rows with non-positive weight
    gp <- merge(g, p, all = TRUE)
    prob_g_zero <- gp$prob_g == 0 & !is.na(gp$prob_g)
    prob_p_zero <- gp$prob_p == 0 & !is.na(gp$prob_p)

    gp <- gp[!prob_g_zero & !prob_p_zero,]

    if(anyNA(gp$prob_g * gp$prob_p)){
        stop("Covariate groups differ between g_weights and p_weights.", call. = FALSE)
    }
    # Add group_id
    groups <- unique(gp[covariates])
    groups <- groups[do.call(order, groups[covariates]), , drop = FALSE]
    groups$group_id <- seq_len(nrow(groups))

    g <- merge(groups, g, all.x = TRUE)
    p <- merge(groups, p, all.x = TRUE)

    # Return canonical structure
    list(
        g_weights = g[, c("group_id", covariates, exposure_time, "prob_g"), drop = FALSE],
        p_weights = p[, c("group_id", covariates, "prob_p"), drop = FALSE]
    )
}




