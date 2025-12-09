# Tests for resolve_hazard_formulas ---------------------------------------
test_that("defaults are used when both inputs are NULL", {
    covars <- c("x1", "x2")
    expo   <- "D_obs"
    
    fm <- resolve_hazard_formulas(
        formula_unexposed = NULL,
        formula_exposed   = NULL,
        covariates = covars,
        exposure_time = expo
    )
    
    # both should be formulas and one-sided (length 2)
    expect_s3_class(fm$formula_unexposed, "formula")
    expect_s3_class(fm$formula_exposed, "formula")
    
    # no response in either
    expect_identical(attr(terms(fm$formula_unexposed), "response"), 0L)
    expect_identical(attr(terms(fm$formula_exposed),   "response"), 0L)
    
    # term labels equal to defaults
    tu <- terms(fm$formula_unexposed)
    te <- terms(fm$formula_exposed)
    
    expect_setequal(attr(tu, "term.labels"), covars)
    expect_true("x1" %in% attr(te, "term.labels"))
    expect_true("x2" %in% attr(te, "term.labels"))
    expect_true(paste0("splines::ns(", expo, ", df = 4)") %in% attr(te, "term.labels"))
})

test_that("formula are returned for valid inputs", {
    covars <- c("x1", "x2")
    expo   <- "D_obs"
    fm <- resolve_hazard_formulas(
        formula_unexposed =  "x1 + I(x2^2)",
        formula_exposed   = "x1 + I(x2^2) + splines::ns(D_obs, df = 4)",
        covariates = covars,
        exposure_time = expo
    )
    
    # both should be formulas and one-sided (length 2)
    expect_s3_class(fm$formula_unexposed, "formula")
    expect_s3_class(fm$formula_exposed, "formula")
})

test_that("errors from string_to_formula are propoagated", {
    covars <- c("x1", "x2")
    expo   <- "D_obs"
   
    expect_error(
        resolve_hazard_formulas(" ~ x1 + x2", NULL, covars, expo),
        regexp = "Do not include '~'.*`formula_unexposed`", fixed = FALSE
    )
    expect_error(
        resolve_hazard_formulas(NULL, "", covars, expo),
        regexp = "invalid", fixed = FALSE
    )
})

test_that("formulas with LHS are returned without LHS", {
    covars <- c("x1", "x2")
    expo   <- "D_obs"
    
    expect_warning(
    fm <- resolve_hazard_formulas(
        formula_unexposed =  as.formula("y1 ~ x1 + I(x2^2)"),
        formula_exposed   = as.formula(" ~ x1 + I(x2^2) + splines::ns(D_obs, df = 4)"),
        covariates = covars,
        exposure_time = expo
    ), 
    )
    
    expect_warning(
        fm <- resolve_hazard_formulas(
            formula_unexposed =  as.formula(" ~ x1 + I(x2^2)"),
            formula_exposed   = as.formula(" y2 ~ x1 + I(x2^2) + splines::ns(D_obs, df = 4)"),
            covariates = covars,
            exposure_time = expo
        ), 
    )
    
    # both should be formulas and one-sided (length 2)
    expect_s3_class(fm$formula_unexposed, "formula")
    expect_s3_class(fm$formula_exposed, "formula")
    
    # no response in either
    expect_identical(attr(terms(fm$formula_unexposed), "response"), 0L)
    expect_identical(attr(terms(fm$formula_exposed),   "response"), 0L)
})



# Tests for validate_formulas ---------------------------------------------

test_that("valid formulas do not error", {
    covars <- c("x1", "x2")
    expo   <- "D_obs"
    expect_no_error(
        validate_formulas(
            formula_unexposed =  as.formula("~ x1 + I(x2^2)"),
            formula_exposed   = as.formula(" ~ x1 + I(x2^2) + splines::ns(D_obs, df = 4)"),
            covariates = covars,
            exposure_time = expo)
    )
})

test_that("Formula for unexposed throws warning for missing covariates", {
    covars <- c("x1", "x2")
    expo   <- "D_obs"
    expect_warning(
        validate_formulas(
            formula_unexposed =  as.formula("~  I(x2^2)"),
            formula_exposed   = as.formula("~ x1 + I(x2^2) + splines::ns(D_obs, df = 4)"),
            covariates = covars,
            exposure_time = expo),
        regexp = "Formula for unexposed"
    )
    
    expect_warning(
        validate_formulas(
            formula_unexposed =  as.formula("~ x1"),
            formula_exposed   = as.formula("~ x1 + I(x2^2) + splines::ns(D_obs, df = 4)"),
            covariates = covars,
            exposure_time = expo),
        regexp = "Formula for unexposed"
    )
})

test_that("Formula for exposed throws warning for missing covariates/time to exposure", {
    covars <- c("x1", "x2")
    expo   <- "D_obs"
    expect_warning(
        validate_formulas(
            formula_unexposed =  as.formula("~ x1 + I(x2^2)"),
            formula_exposed   = as.formula("~ I(x2^2) + splines::ns(D_obs, df = 4)"),
            covariates = covars,
            exposure_time = expo),
        regexp = "Formula for exposed"
    )
    
    expect_warning(
        validate_formulas(
            formula_unexposed =  as.formula("~ x1 + I(x2^2)"),
            formula_exposed   = as.formula("~ x1 + I(x2^2)"),
            covariates = covars,
            exposure_time = expo),
        regexp = "Formula for exposed"
    )
})





