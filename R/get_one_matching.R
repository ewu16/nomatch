#' Internal function to compute cumulative incidence and effectiveness measures
#' for a matched dataset
#'
#' @description First, creates the analysis matched data set based on
#' `matched_data`. Then computes Kaplan Meier estimates of cumulative incidence
#' in the matched analysis data set.
#'
#' @inheritParams matching
#' @param keep_adata Logical; return the matched analysis dataset? Default: TRUE
#' @return The object returned by [compute_km_ve()] with the matched analysis
#'   dataset, `matched_adata`, attached, if requested.
#'   
#' @keywords internal
#' 
get_one_matching <- function(matched_data,
                                outcome_time,
                                outcome_status,
                                exposure,
                                exposure_time,
                                tau,
                                timepoints,
                                keep_models = TRUE,
                                keep_adata = TRUE){


    matched_adata <- make_matched_adata(matched_data = matched_data,
                                        outcome_time = outcome_time,
                                        outcome_status = outcome_status,
                                        exposure = exposure,
                                        exposure_time = exposure_time,
                                        tau = tau)

   out <- compute_km_ve(adata = matched_adata,
                        adata_outcome_name = "match_T",
                        adata_event_name =  paste0("match_", outcome_status),
                        adata_trt_name = paste0("match_", exposure),
                        timepoints = timepoints,
                        keep_models = keep_models)

   check_pt_estimates(out$pt_estimates)

   #add some additional info
   if(keep_adata){
       out$matched_adata <- matched_adata
   }

   return(out)

}


#' Get analysis matched dataset from matched cohort
#'
#' @description This function modifies the original matched data, preparing it for use in
#' analysis. Namely, this includes
#' *  (if requested) censoring matched pairs at the time
#' of the unvaccinated individual's vaccination
#' * creating the outcome time from matching index date to endpoint
#' * restricting to matched  pairs in which both individuals are eligible `tau` days after the matching
#' index date
#'

#' @param matched_data A data frame representing a matched cohort
#' @param outcome_time Character string representing the original outcome variable in  `matched_data`
#' @param outcome_status Character string representing the original event variable in `matched_data`
#' @param exposure Character string representing the original vaccination status variable in `matched_data`
#' @param exposure_time Character string representing the original vaccination time variable in `matched_data`
#' @param tau The time excluded after vaccination to allow building up of
#'   immunity
#' @return A data frame representing the analysis-eligible matched dataset (a subset of `matched data`).
#' Contains all variables in `matched_data`, plus the following variables
#' \describe{
#' \item{match_<outcome_status>}{Event variable after adjusting for pair_censoring}
#' \item{match_<outcome_time>}{Outcome variable after adjusting for pair_censoring}
#' \item{match_T}{Time to endpoint from matched index date after adjusting for pair_censoring}
#' }
#' @keywords internal
#'
make_matched_adata <- function(matched_data, outcome_time, outcome_status, exposure, exposure_time,
                               tau){
    
    # Check for name conflicts for variables that will be added
    reserved_vars <- c("match_T", paste0("match_", outcome_time), paste0("match_", outcome_status))
    conflicts <- intersect(reserved_vars, names(matched_data))
    
    if(length(conflicts) > 0) {
        stop(
            "Data contains reserved variable names: ",
            paste(conflicts, collapse = ", "), "\n",
            "Please rename these variables before matching.",
            call. = FALSE
        )
    }
    
    #Censor cases and controls if control later gets vaccinated
    ordered_data <- matched_data[order(matched_data$match_id),]
    #new variables
    new_outcome_name <- paste0("match_", outcome_time)
    new_event_name <- paste0("match_", outcome_status)
    ordered_data[[new_outcome_name]] <-ordered_data[[new_event_name]] <-  NA
    
    cases <- ordered_data[ordered_data$match_type == "case",]
    controls <-  ordered_data[ordered_data$match_type == "control",]
    stopifnot(all.equal(cases$match_id, controls$match_id))
    
    #controls who later got vaccinated
    trt_later <- controls[[exposure]] == 1 & 
        controls[[exposure_time]] < controls[[outcome_time]]
    
    #censor controls at their time of vaccination
    controls[[new_outcome_name]] <- ifelse(trt_later, controls[[exposure_time]], controls[[outcome_time]])
    controls[[new_event_name]]   <- ifelse(trt_later, 0, controls[[outcome_status]])
    
    #censor cases at control's time of vaccination if still at risk
    at_risk <-   cases[[outcome_time]] >  controls[[exposure_time]]
    cases[[new_outcome_name]] <- ifelse(trt_later & at_risk, controls[[exposure_time]], cases[[outcome_time]])
    cases[[new_event_name]]   <- ifelse(trt_later & at_risk, 0, cases[[outcome_status]])
    
    tmp_data <- rbind(cases, controls)
    
    #Create matching analysis dataset
    ##compute survival time from matching index date
    tmp_data$match_T <- tmp_data[[new_outcome_name]] - tmp_data$match_index_time
    ##exclude pairs where at least one individual is no longer at risk by tau
    excluded_match_ids <- unique(tmp_data$match_id[!tmp_data$match_T > tau])
    adata <- subset(tmp_data, !tmp_data$match_id %in% excluded_match_ids)
    
    adata[order(adata$match_id),]
    
}


#' Get the marginal cumulative incidence in the exposed and unexposed
#' groups based on Kaplan Meier estimation.
#'
#' @inheritParams matching
#' @param adata A data frame that represents the analysis data set
#' @param adata_outcome_name Character string specifying the time to event variable in `adata`. The
#' time should be the time to event from exposure or matching 
#' @param adata_event_name Character string specifying the event variable in `adata`
#' @param adata_trt_name Character string specifying the treatment variable in `adata`
#'
#' @return A list containing the following:
#' \describe{
#' \item{estimates}{A matrix of estimates. The columns of the matrix are the cumulative
#'  incidence/effectiveness terms and the rows are the requested timepoints for evaluation.}
#' \item{models}{If `keep_models = TRUE`, the fitted Kaplan Meier object}
#' }
#'
#' @keywords internal

compute_km_ve <- function(adata,
                          adata_outcome_name,
                          adata_event_name,
                          adata_trt_name,
                          timepoints,
                          keep_models = TRUE){

    # Build formula and fit KM
    outcome <- paste0("survival::Surv(", adata_outcome_name, ",", adata_event_name, ")")
    km_fit <- survival::survfit(stats::reformulate(adata_trt_name, response = outcome) , data = adata)
    
    if(all(km_fit$n.event == 0)){
        warning("No events observed in matched analysis dataset. All cumulative incidence estimates will be 0.",
                call. = FALSE)
    }

    # Summarize at requested times
    surv_probs <- summary(km_fit, times = timepoints)

    ## survival probabilities fo each group (assumed to have 2 groups_
    surv_0 <- surv_probs$surv[surv_probs$strata == levels(surv_probs$strata)[1]]
    surv_1 <- surv_probs$surv[surv_probs$strata == levels(surv_probs$strata)[2]]

    #Compute cumulative incidences and derived effect measures
    cuminc_0 <- 1 - surv_0
    cuminc_1 <- 1 - surv_1
    rd <- cuminc_1 - cuminc_0
    rr <- cuminc_1/cuminc_0
    ve <- 1 - rr

    est  <- cbind(
        cuminc_0 = cuminc_0,
        cuminc_1 = cuminc_1,
        risk_difference = rd,
        risk_ratio = rr,
        relative_risk_reduction = ve
    )

    rownames(est) <- timepoints

    out <- list(
        pt_estimates = est,
        models = if(keep_models) km_fit else NULL
    )
}





