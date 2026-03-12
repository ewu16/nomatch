#' Validate inputs to main functions
#'
#' @description
#' Internal validation function to check that all inputs to `nomatch()/matching()` are
#' properly formatted and logically consistent.
#'
#' @inheritParams nomatch
#'
#' @return Invisibly returns `NULL` if all checks pass. Throws an error with
#'   descriptive message if any validation fails.
#'
#' @keywords internal
#' @noRd
validate_data_inputs <- function(data, core_args, covariates = NULL) {
    validate_data_str(data)
    validate_column_args(core_args, covariates)
    validate_cols_exist(data, core_args, covariates)
    
    outcome_time <- core_args$outcome_time
    exposure <- core_args$exposure
    exposure_time <- core_args$exposure_time
    outcome_status <- core_args$outcome_status
    
    validate_outcome_time(data, outcome_time)
    validate_exposure_args(data, exposure, exposure_time)
    if (!is.null(outcome_status)) validate_outcome_status(data, outcome_status)
    if (!is.null(covariates)) validate_covariates(data, covariates)
    
    invisible(NULL)
}

## Structural data checks -------------
validate_data_str <- function(data){
    if(!is.data.frame(data)) stop("'data' must be a data.frame, not ", class(data)[1])
    if(nrow(data) == 0)  stop("'data' has no rows")
    invisible(NULL)
}


validate_column_args <- function(core_args, covariates = NULL){
    if(!is.list(core_args) || is.null(names(core_args))){
        stop("`core_args` must be a named list")
    }
    # Check that elements of core_args are a single string
    not_single_string <- !vapply(core_args, \(x) is.character(x) && length(x) == 1, logical(1))
    if (any(not_single_string)) {
        stop("The following must be single, non-null, character strings: ",
             paste(names(core_args)[not_single_string], collapse = ", "), call. = FALSE)
    }
    # Check covariates is character vector or NULL
    not_vector <- !(is.null(covariates) || (is.character(covariates) && length(covariates) >= 1))
    if (not_vector) {
        stop("`covariates` must be be a character vector or NULL")
    }
}

validate_cols_exist <- function(data, core_args, covariates = NULL){
    # core_args should be a named list 
    #Check that columns exist in data
    all_args <- c(core_args, list(covariates = covariates))
    lapply(seq_along(all_args), \(i){
        cols <- all_args[[i]]
        missing_cols <- setdiff(cols, names(data))
        if(length(missing_cols) > 0) {
            stop("Missing `", paste(names(all_args)[i]), "` column(s) in data: ", paste(missing_cols, collapse = ", "), call. = FALSE)
        }
    })
}


## Column-level data checks --------
validate_outcome_time <- function(data, outcome_time){
    if (anyNA(data[[outcome_time]])) {
        stop("Missing values in `outcome_time` are not supported. Please remove these observations.")
    }
    
    if (!is.numeric(data[[outcome_time]])) {
        stop("Outcome time variable '", outcome_time, "' must be numeric.")
    }
    
    if (any(data[[outcome_time]] < 0)) {
        stop("Outcome time variable '", outcome_time, "' must be non-negative.")
    }
    invisible(NULL)
}


validate_outcome_status <- function(data, outcome_status){
    if (anyNA(data[[outcome_status]])) {
        stop("Missing values in `outcome_status` are not supported. Please remove these observations.")
    }
    if (!all(data[[outcome_status]] %in% c(0,1))) {
        stop("Outcome status variable '", outcome_status, "' must be coded as 0/1 numeric.")
    }
    invisible(NULL)
}

validate_exposure_args <- function(data, exposure, exposure_time){
    # Check exposure
    if (anyNA(data[[exposure]])) {
        stop("Missing values in <exposure> are not supported. Please remove these observations.")
    }
    if (!all(data[[exposure]] %in% c(0,1))) {
        stop("Exposure variable '", exposure, "' must be coded as 0/1 numeric.")
    }
    
    # Check exposure_time
    et <- data[[exposure_time]]
    if (!is.numeric(et)) {
        stop("Exposure time variable '", exposure_time, "' must be numeric.")
    }
    
    if (any(et[!is.na(et)] < 0)) {
        stop("Exposure time variable '", exposure_time, "' must be non-negative.")
    }
    exposed <- data[[exposure]] == 1
    if (any(is.na(et[exposed])))       stop("All exposed individuals must have non-missing ", exposure_time, ".")
    if (any(!is.na(et[!exposed])))     stop("All unexposed individuals must have ", exposure_time, " set to NA.")
    
    invisible(NULL)
}


validate_covariates <- function(data, covariates){
    if(!is.null(covariates) && sum(is.na(data[, covariates])) > 0){
        stop("Missing values in covariates ",
             paste(paste0("'", covariates, "'"), collapse = ", "), "' are not supported. ",
             "Please remove these observations.")
    }
    invisible(NULL)
}

