string_to_formula <- function(rhs, label){
    form_label <- paste0("`", label, "`")
    
    if(!(is.character(rhs) && length(rhs) == 1)){
        stop(form_label, " must be a single string")
    }
    if(grepl("~", rhs, fixed = TRUE)){
        stop("Do not include '~' in ", form_label, "; provide only the right-hand side.")
    }
    
    rhs_trim <- trimws(rhs)
    tryCatch({
        stats::as.formula(paste("~", rhs_trim))
    }, error = function(e) {stop(form_label, " appears to be invalid.")}
    )
}

check_lhs <- function(formula, label){
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
