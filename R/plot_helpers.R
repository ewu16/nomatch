#' Internal plotting function for VE panel plots
#'
#' @description
#' Creates a three-panel plot showing cumulative incidence (no vaccine),
#' cumulative incidence (vaccine), and vaccine effectiveness over time.
#' This function handles both single method and comparison plots.
#'
#' @param plot_data A data frame containing estimates to plot. Must include
#'   columns: `t0` (time), `estimate`, `term` (one of "cuminc_0", "cuminc_1", "vaccine_effectiveness"),
#'   `<ci_type>_lower`, `<ci_type>_upper`, and `method`
#' @param ci_type Character string specifying the type of confidence interval.
#'   One of "wald", "percentile", "simul", or "none"
#' @param alpha Numeric significance level for confidence intervals (e.g., 0.05 for 95% CIs).
#' @param colors Character vector of colors. Length should match the number of
#'   unique methods in `plot_data`. If `NULL` or length doesn't match, ggplot2's
#'   default colors are used.
#'
#' @return A ggplot2 object with three faceted panels.
#'
#' @details
#' This is an internal helper function used by `plot.vefit()` and `compare_ve_fits()`.
#' The function automatically detects whether it's plotting a single method or
#' comparing multiple methods based on the `method` column in `plot_data`.
#'
#' For cumulative incidence panels, y-axis limits are shared to facilitate comparison.
#' The VE panel uses free y-axis scaling.
#'
#' @keywords internal
#' @noRd

plot_ve_panel <- function(plot_data,
                           ci_type,
                           alpha,
                           effect,
                           trt_0_label = "Unexposed",
                           trt_1_label = "Exposed",
                           colors = NULL) {

  n_methods <- length(unique(plot_data$method))
  has_ci <- (ci_type != "none")
  terms <- c("cuminc_0", "cuminc_1", effect)
  plot_data <- plot_data[plot_data$term %in% terms,]

  #Validate inputs
  if(!has_ci){
    lower <- "estimate"
    upper <- "estimate"
  }else{
    lower <- paste0(ci_type, "_lower")
    upper <- paste0(ci_type, "_upper")

    if(!all(c(lower, upper) %in% names(plot_data))){
      stop(paste0("Confidence intervals for ci_type = ", ci_type,
                  "were not found in plot_data"))
    }
  }



  #Set defaults for plot based on effect
  if(effect == "risk_difference"){
    d <- list(label = "Risk Difference",
              x = "Time since exposure",
              label_function = scales::label_percent())
  }else if(effect == "risk_ratio"){
    d <- list(label = "Risk Ratio",
              x = "Time since exposure",
              label_function =  ggplot2::waiver())

  }else if(effect == "vaccine_effectiveness"){
    d <- list(label = "Vaccine Effectiveness",
              x = "Time since vaccination",
              label_function = scales::label_percent())
  }

  #Clean up term labels before plotting
  cuminc_0_label <-  paste("Cumulative Incidence:\n", trt_0_label)
  cuminc_1_label <-  paste("Cumulative Incidence:\n", trt_1_label)
  term_labels <-  c(cuminc_0_label, cuminc_1_label, d$label)

  plot_data$term_label <- factor(plot_data$term, terms, term_labels)


  #Make main plot
  p <- ggplot2::ggplot(plot_data,
          ggplot2::aes(x = .data$t0,
                       y = .data$estimate,
                       color = .data$method,
                       fill = .data$method)) +
    ggplot2::geom_line()


  # For panel plot
  ## Calculate shared y-limits for cumulative incidence
  cuminc_data <- plot_data[plot_data$term %in% c("cuminc_0", "cuminc_1"), ]
  y_min <- min(cuminc_data[[lower]], na.rm = TRUE)
  y_max <- max(cuminc_data[[upper]], na.rm = TRUE)

  p <- p + ggplot2::facet_wrap(~term_label, scales = "free") +
    ggh4x::facetted_pos_scales(
      y = list(
        ggplot2::scale_y_continuous(labels = d$label_function,
                                    limits = c(y_min, y_max)),
        ggplot2::scale_y_continuous(labels = d$label_function,
                                    limits = c(y_min, y_max)),
        ggplot2::scale_y_continuous(labels = d$label_function)
      )
    )


  # Add confidence interval bands
  if (has_ci) {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[lower]], ymax = .data[[upper]]),
                           alpha = if(n_methods > 1) 0.2 else 0.3,
                           linewidth = if(n_methods > 1) 0.1 else 0, show.legend = FALSE)
  }

  # Add color scales if appropriate, otherwise use default scales
  if (n_methods == length(colors)) {
    p <- p +
      ggplot2::scale_color_manual(values = colors, name = "Method") +
      ggplot2::scale_fill_manual(values = colors, name = "Method")
  }

  # Make captions about confidence intervals
  ci_labels <- c(wald = "pointwise Wald",
                 percentile = "pointwise percentile",
                 simul = "simultaneous")

  caption <- ifelse(has_ci,
                    paste("\nShaded areas represent",  paste0((1 - alpha) * 100, "%"),
                          ci_labels[ci_type], "confidence intervals."),
                    "\nPoint estimates only (no confidence intervals computed/requested).")

  #Add plot labels
    p <- p +
    ggplot2::labs(
      x = d$x,
      y = "Estimate",
      caption = caption
    )

  # Add theme
  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0),
      strip.background = ggplot2::element_rect(fill = "white", color = "white"),
      strip.text = ggplot2::element_text(size = 11, colour = "black"),
      axis.text = ggplot2::element_text(size = 9, colour = "black"),
      axis.title = ggplot2::element_text(size = 11, colour = "black")
    )


  # Add legend positioning if using method colors
  if (n_methods == 1) {
    p <- p + ggplot2::theme(
      legend.position = "none"
    )
  }else{
    p <- p + ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 10)
    )
  }

  return(p)
}


#' Validate and normalize CI type for plotting
#'
#' @param object A vefit object
#' @param ci_type Character string or NULL
#'
#' @return A validated ci_type string for plotting (one of "wald", "percentile", "simul")
#' @keywords internal
#' @noRd
validate_confint_type <- function(object, ci_type = NULL) {

  if(object$boot_reps == 0){
    return("none")
  }

  # Use object's ci_type if not specified
  if (is.null(ci_type)) {
    ci_type <- object$ci_type
  }

  if(ci_type  == "none"){
    return("none")
  }
  # Handle "both" - default to "wald"
  if (ci_type == "both") {
    ci_type <- "wald"
  }

  # Validate it's one of the allowed types
  if (!ci_type %in% c("wald", "percentile", "simul", "none")) {
    stop("ci_type must be one of 'wald', 'percentile', 'simul', or 'none'")
  }

  # Check if simultaneous CIs are requested but not computed
  if (ci_type == "simul") {
    if (is.null(object$simul_z_star)) {
      stop("Simultaneous CIs not yet computed. Call add_simultaneous_ci()")
    }
  }


  #Check if requested CI type exists in estimates
  lower <- paste0(ci_type, "_lower")
  upper <- paste0(ci_type, "_upper")
  if (!all(c(lower, upper) %in% colnames(object$estimates[[tolower(object$effect)]]))) {
    stop("CI type '", ci_type, "' not found in object.\n")
  }


  return(ci_type)
}



