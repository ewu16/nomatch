#' Print method for effectiveness fits
#'
#' @description
#' Prints a concise summary of effectiveness estimates from a fitted model.
#'
#' @param x An object of class `nomatchfit` created by [nomatch()] or [matching()].
#' @param digits Integer indicating the number of decimal places to display. Default is 3.
#' @param effect The effect measure to output. By default,
#' this is the `effect` measure chosen in the main function call. Must be one of `risk_difference`, `risk_ratio`,
#' or `relative_risk_reduction`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' Invisibly returns the input object `x`. Called primarily for its side effect
#' of printing a summary to the console.
#'
#'
#'
#' @export
print.nomatchfit <- function(x, digits = 3, effect = NULL,...) {
    if(is.null(effect)){
        effect <- x$effect
    }

    if(x$effect == "risk_difference"){
        title <- "Risk Difference Estimates"
    }else if(x$effect == "risk_ratio"){
        title <- "Risk Ratio Estimates"
    }else if(x$effect == "relative_risk_reduction"){
        title <- "Vaccine Effectiveness Estimates"
    }else{
        stop("Effect must be one of 'risk_difference', 'risk_ratio', or 'relative_risk_reduction'")
    }

    cat("\n", title, "\n")
    cat(strrep("=", 50), "\n")

    # Key identifying info
    #cat("Method:", x$method, "\n")
    cat("Call:", paste(deparse(x$call), collapse = "\n"), "\n")

    # The main result
    cat("\nResult:\n")
    effect <- stats::setNames(list(x$estimates[[effect]]), effect)
    effect_df <- estimates_to_df(effect)

    ci_level <- (1-x$alpha)*100
    name_map <- c(
        t0                = "Timepoint",
        estimate          = "Estimate",
        wald_lower        = paste0(ci_level, "% Wald CI: Lower"),
        wald_upper        = paste0(ci_level, "% Wald CI: Upper"),
        wald_pval         = "Wald p-value",
        percentile_lower  = paste0(ci_level, "% Percentile CI: Lower"),
        percentile_upper  = paste0(ci_level, "% Percentile CI: Upper"),
        percentile_pval   = "Percentile p-value",
        simul_lower  = paste0(ci_level, "% Simul CI: Lower"),
        simul_upper  = paste0(ci_level, "% Simul CI: Upper")
    )

    display_cols <- intersect(names(effect_df), names(name_map))

    display_labels <- name_map[display_cols]
    names(effect_df)[names(effect_df) %in% display_cols] <- display_labels

    print(utils::head(effect_df[display_labels], 10), digits = digits)

    # Print how many rows remain
    remaining <- nrow(effect_df) - 10
    if(remaining > 0){
        cat("\nRemaining rows:", remaining, "\n")
    }

    # Hint for more info
    cat("\nUse summary() for more details\n")
    cat("Use plot() to visualize results\n")


    invisible(x)
}

#' Summary method for effectiveness fits
#' @description
#' Summarizes how effectiveness estimates were obtained
#'
#' @param object An object of class `nomatchfit` created by [nomatch()] or [matching()].
#' @param digits Integer indicating the number of decimal places to display. Default is 4.
#' @param show_models Logical; print model details? Default is FALSE.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' Invisibly returns the input object `object`. Called primarily for its side effect
#' of printing a detailed summary to the console.
#'
#' @export
summary.nomatchfit <- function(object, digits = 4, show_models = FALSE,...) {
    if(object$effect == "risk_difference"){
        title <- "Risk Difference"
    }else if(object$effect == "risk_ratio"){
        title <- "Risk Ratio"
    }else if(object$effect == "relative_risk_reduction"){
        title <- "Vaccine Effectiveness"
    }

    cat("\n")
    cat(strrep("=", 70), "\n")
    cat(title, "Analysis Summary\n")
    cat(strrep("=", 70), "\n\n")

    # ---- Key Parameters ----
    cat("Method:             ", object$method, "\n")
    cat("Evaluation times:   ", paste(utils::head(object$timepoints), collapse = ", "),
        ifelse(length(object$timepoints) > 6, ", ...", ""), "\n")
    cat("Immune lag:         ", object$immune_lag, "\n")

    if(object$method ==  "nomatch (G-computation)"){
        if (length(object$covariates) > 0) {
            cat("Adjusted for:       ", paste(object$covariates, collapse = ", "), "\n")
        }
    }



    # ---- Bootstrap Info (if applicable) ----
    if (object$boot_reps > 0) {
        cat("\nBootstrap:          ", object$boot_reps, "replicates\n")
        cat("Confidence level:   ", (1 - object$alpha) * 100, "%\n")
        cat("Successful samples: ",
            paste(range(object$n_success_boot), collapse = "-"),
            " (range across timepoints)\n")
    } else {
        cat("\nNo bootstrap performed (boot_reps = 0)\n")
    }


    if(object$method ==  "nomatch (G-computation)"){
        cat("\n")
        cat(strrep("-", 70), "\n")
        cat("Sample:\n")
        cat(strrep("-", 70), "\n")

        descrip <- object$descrip
        cat("N total:", descrip$n, "\n")
        cat("Number of events:", descrip$n_events, "\n")
        cat("\n")
        cat("N exposed:", descrip$n_exposed, "\n")
        cat("N exposed at-risk <immune_lag> days after exposure:", descrip$n_exposed_at_tau, "\n")

        cat("\n")
        cat("Distribution of exposure times among at-risk <immune_lag> days after exposure:\n")
        cat(" Range: ", paste(range(descrip$exposure_times_at_tau), collapse = " - "), "| ")
        cat(" Median (IQR): ",
            stats::median(descrip$exposure_times_at_tau),
            paste0("(", stats::quantile(descrip$exposure_times_at_tau, 0.25), " - ",
                   stats::quantile(descrip$exposure_times_at_tau, 0.75), ")"), "| ")
        cat(" Mean: ",  round(mean(descrip$exposure_times_at_tau), 1))

        # ---- Main Results ----
        if(!is.null(object$model_0) & !is.null(object$model_1)){
            cat("\n\n")
            cat(strrep("-", 70), "\n")
            cat("Model for unexposed:\n")
            cat(strrep("-", 70), "\n")
            cat("N =", object$model_0$n, "| Number of events =", object$model_0$nevent, "\n\n")
            if(show_models){
                print(round(stats::coef(summary(object$model_0)), 3))
            }else{
                cat("Use '$model_0' to see model details.\n")
            }

            cat("\n")
            cat(strrep("-", 70), "\n")
            cat("Model for exposed:\n")
            cat(strrep("-", 70), "\n")
            cat("N =", object$model_1$n, "| Number of events =", object$model_1$nevent, "\n\n")
            if(show_models){
                print(round(stats::coef(summary(object$model_1)), 3))

            }else{
                cat("Use '$model_1' to see model details.\n")
            }

        }

    }else{
        cat("\n")
        cat(strrep("-", 70), "\n")
        cat("Sample:\n")
        cat(strrep("-", 70), "\n")

        descrip <- object$descrip
        cat("N matched:", descrip$n_matched, "\n")
        cat("N matched analysis:",  descrip$n_matched_analysis, "\n",
            "(excludes pairs in which either individual is not at risk\n `immune_lag` days after matching index day)\n\n")
        cat("Number of events in matched analysis:", descrip$n_events, "\n")
        cat("\n")
        cat("Distribution of exposure times in matched analysis:\n")
        cat(" Range: ", paste(range(descrip$exposure_times), collapse = " - "), "| ")
        cat(" Median (IQR): ",
            stats::median(descrip$exposure_times),
            paste0("(", stats::quantile(descrip$exposure_times, 0.25), " - ",
                        stats::quantile(descrip$exposure_times, 0.75), ")"), "| ")
        cat(" Mean: ",  round(mean(descrip$exposure_times), 1))

        cat("\n\n")
        cat(strrep("-", 70), "\n")
        cat("Kaplan Meier for matched analysis:\n")
        cat(strrep("-", 70), "\n\n")
        object$models$call <- NULL
        print(object$models)

    }


    cat("\n")
    cat(strrep("=", 70), "\n\n")

    invisible(object)
}

#' Plot method for nomatchfit objects
#'
#' @description
#' Create a panel plot of cumulative incidence and effectiveness estimates across all evaluation time points.
#'
#'
#' @param x An object of class `nomatchfit` created by [nomatch()] or [matching()].
#' @param effect The effect measure to plot next to the cumulative incidence plots. By default,
#' this is the `effect` measure chosen in the main function call. Must be one of `risk_difference`, `risk_ratio`,
#' or `relative_risk_reduction`.
#' @param ci_type Character string specifying the type of confidence interval band to plot. By default,
#'   `"wald"` if available, otherwise set to `"percentile"` or `none`.
#'   One of `"wald", "percentile", "simul"`, or `"none"`. Must choose a `ci_type` whose lower
#'   and upper bounds are already computed in `estimates` component of `x`.
#' @param color Aesthetic value to map data values to. Default: `"#0072B2"` (blue)
#' @param ... Additional arguments (currently ignored).
#' @return a ggplot2 object with three faceted panels (for cumulative incidences and the chosen effect measure)
#'
#' @details
#' For cumulative incidence panels, y-axis limits are shared across methods to
#' facilitate comparison. The VE panel uses free y-axis scaling.
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' fit <- nomatch(
#'  data = simdata,
#'  outcome_time = "Y",
#'  outcome_status = "event",
#'  exposure = "V",
#'  exposure_time = "D_obs",
#'  covariates = c("x1", "x2"),
#'  timepoints = seq(30, 180, by = 30),
#'  immune_lag = 14,
#'  boot_reps = 5,
#'  n_cores = 2
#' )
#' plot(fit)
#'
plot.nomatchfit <- function(x, effect = NULL, ci_type = NULL, color = "#0072B2", ...) {

    if(is.null(effect)){
        effect <- x$effect
    }

    if(length(x$timepoints) < 2){
        message("Only one evaluation time supplied; returning point estimate(s) instead of a plot")
        est_df <- estimates_to_df(x$estimates)
        print(est_df)
        return(invisible(est_df))
    }

    ci_type <- validate_confint_type(x, ci_type)
    plot_data <- estimates_to_df(x$estimates)
    plot_data$method <- "dummy"
    alpha <- x$alpha
    plot_ve_panel(plot_data, ci_type, alpha, colors = color,
                  effect = effect,
                  trt_0_label = paste(x$exposure, "= 0"),
                  trt_1_label = paste(x$exposure, "= 1"))
}

#' Check if x is a nomatchfit
#' @keywords internal
#' @noRd
is.nomatchfit <- function(x) {
    inherits(x, "nomatchfit")
}

