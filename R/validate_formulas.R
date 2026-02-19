resolve_hazard_formulas <- function(formula_unexposed, formula_exposed, covariates, exposure_time) {
    if(is.null(formula_unexposed)){
        formula_unexposed <- if(is.null(covariates)) (~ 1) else stats::reformulate(covariates)

    }else if(!inherits(formula_unexposed, "formula")){
        formula_unexposed <- string_to_formula(formula_unexposed, "formula_unexposed")
    }

    if(is.null(formula_exposed)){
        d_term <- paste0("splines::ns(", exposure_time, ", df = 4)")
        formula_exposed <- stats::reformulate(c(covariates, d_term))
    }else if(!inherits(formula_exposed, "formula")){
        formula_exposed <- string_to_formula(formula_exposed, "formula_exposed")
    }
    
    formula_unexposed <- check_lhs(formula_unexposed, "formula_unexposed")
    formula_exposed <- check_lhs(formula_exposed, "formula_exposed")
    
    list(formula_unexposed = formula_unexposed, formula_exposed = formula_exposed)
}



validate_formulas <- function(formula_unexposed, formula_exposed, covariates, exposure_time){
    unexposed_vars <- get_rhs_vars(formula_unexposed)
    exposed_vars <- get_rhs_vars(formula_exposed)
    
    extra_unexposed <- setdiff(unexposed_vars, covariates)
    extra_exposed <- setdiff(exposed_vars, c(covariates, exposure_time))
    
    if(length(extra_unexposed) != 0){
        stop("`formula_unexposed` contains terms that are not included in `covariates`: ", paste(extra_unexposed, collapse = ", "), call. = FALSE)
    }
    if(length(extra_exposed) != 0){
        stop("`formula_exposed` contains terms that are not included in `covariates`: ", paste(extra_exposed, collapse = ", "), call. = FALSE)
    }
    
    missing_unexposed <-  setdiff(covariates, unexposed_vars)
    missing_exposed   <-  setdiff(covariates, exposed_vars)

    if(length(missing_unexposed) != 0){
        warning("`formula_unexposed` does not include a term for the following `covariates`: ", paste(missing_unexposed, collapse = ","), call.= FALSE)
    }
    
    if(length(missing_exposed) != 0){
        warning("`formula_exposed` does not include a term for the following `covariates`: ", paste(missing_exposed, collapse = ", "), call.= FALSE)
    }
    
    if(!exposure_time %in% get_rhs_vars(formula_exposed)){
        warning("`formula_exposed` does not include a term for `exposure_time`",  call.= FALSE)
    }
    
    invisible(NULL)
}


