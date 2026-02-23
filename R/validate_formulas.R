resolve_hazard_formulas <- function(formula_unexposed, formula_exposed, covariates, exposure_time) {
    resolve_formula <- function(f, default, label) {
        if (is.null(f)) return(default)
        if (inherits(f, "formula")) return(f)
        if (is.character(f) && length(f) > 1) return(stats::reformulate(f))
        string_to_formula(f, label)
    }
    
    # default formulas
    default_unexposed <- if (is.null(covariates)) (~ 1) else stats::reformulate(covariates)
    d_term <- paste0("splines::ns(", exposure_time, ", df = 4)")
    default_exposed <- stats::reformulate(c(covariates, d_term))
    
    # resolve
    formula_unexposed <- resolve_formula(formula_unexposed, default_unexposed, "formula_unexposed")
    formula_exposed   <- resolve_formula(formula_exposed,   default_exposed,   "formula_exposed")
    
    # drop LHS of formulas
    formula_unexposed <- force_one_sided(formula_unexposed, "formula_unexposed")
    formula_exposed   <- force_one_sided(formula_exposed,   "formula_exposed")
    
    list(formula_unexposed = formula_unexposed,
         formula_exposed   = formula_exposed)
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



# Helpers for working with formulas ---------------------------------------
string_to_formula <- function(rhs, label){
    form_label <- paste0("`", label, "`")
    
    if(!(is.character(rhs) && length(rhs) == 1)){
        stop(form_label, " must be a single string")
    }
    
    rhs_trim <- trimws(rhs)
    form <- ifelse(grepl("~", rhs_trim, fixed = TRUE), rhs_trim, paste("~", rhs_trim))
    
    tryCatch({
        stats::as.formula(form)
    }, error = function(e) {stop(form_label, " appears to be invalid.", call. = FALSE)}
    )
}

force_one_sided <- function(formula, label){
    form_label <- paste0("`", label, "`")
    
    if(attr(stats::terms(formula), "response") == 1){
        warning("Left hand side of ", form_label, " is not allowed. Only the right hand side will be used")
    }
    stats::update(formula, NULL ~ . )
}

get_rhs_vars <- function(f) {
    if (length(f) == 3){
        all.vars(f[[3]])
    }else{
        all.vars(f[[2]])
    }
}


