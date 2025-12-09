resolve_hazard_formulas <- function(formula_unexposed, formula_exposed, covariates, exposure_time) {
    
    if(is.null(formula_unexposed)){
        formula_unexposed <- stats::reformulate(covariates)
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
 

    if(any(!covariates %in% get_rhs_vars(formula_unexposed))){
        warning("Formula for unexposed does not include all variables in `covariates`")
    }
    
    if(any(!covariates %in% get_rhs_vars(formula_exposed))){
        warning("Formula for exposed does not include all variables in `covariates`")
    }
    
    if(!exposure_time %in% get_rhs_vars(formula_exposed)){
        warning("Formula for exposed does not include a term for `exposure_time`")
    }
    
   
    invisible(NULL)
}


