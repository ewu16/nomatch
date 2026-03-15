# Tests for formula helpers -----------------------------------------------
test_that("string_to_formula() passes on clean input", {
    #string without tilde
    f <- string_to_formula("x1 + x2", label = "f")
    expect_setequal(all.vars(f), c("x1", "x2"))
    expect_identical(attr(terms(f), "response"), 0L)
    
    #string with tilde and no outcome
    f <- string_to_formula(" ~ x1 + x2", label = "f")
    expect_setequal(all.vars(f), c("x1", "x2"))
    expect_identical(attr(terms(f), "response"), 0L)
    
    #string with LHS and tilde
    f <- string_to_formula("y ~ x1 + x2", label = "f")
    expect_setequal(all.vars(f), c("y", "x1", "x2"))
    expect_identical(attr(terms(f), "response"), 1L)
    
})

test_that("string_to_formula() requires single string", {
    expect_error(string_to_formula(c("x1", "x2"), label = "f"), "single string")
    expect_error(string_to_formula(1, label = "f"), "single string")
})

test_that("get_rhs_vars() handles one- and two-sided formulas and formulas with special operations", {
    expect_setequal(get_rhs_vars(~ x1 + x2), c("x1", "x2"))
    expect_setequal(get_rhs_vars(y ~ x1 + x2), c("x1", "x2"))
    expect_setequal(get_rhs_vars(y ~ x1 + I(x2*1000)), c("x1", "x2"))
    expect_setequal(get_rhs_vars(y ~ x1*x2), c("x1", "x2"))
    expect_identical(get_rhs_vars(~ 1), character(0))
})

test_that("force_one_sided() passes on one-sided formula and returns original formula", {
    expect_no_warning(f <- force_one_sided(~ x1, "f"))
    expect_identical(f, ~ x1) 
})

test_that("force_one_sided() warns and strips LHS", {
    expect_warning(force_one_sided(y ~ x1, "f"), "not allowed")
    f <- suppressWarnings(force_one_sided(y ~ x1, "f"))
    expect_identical(attr(terms(f), "response"), 0L)
})



# Tests for resolve_hazard_formulas() ---------------------------------------
test_that("resolve_hazard_formulas() autoresolves when both formulas are NULL and covariates are non-NULL", {
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

test_that("resolve_hazard_formulas() produces intercept-only unexposed formula when covariates are NULL", {
    fm <- resolve_hazard_formulas(NULL, NULL, covariates = NULL, exposure_time = "D_obs")
    
    expect_s3_class(fm$formula_unexposed, "formula")
    expect_identical(attr(terms(fm$formula_unexposed), "term.labels"), character(0))
    
    expect_s3_class(fm$formula_exposed, "formula")
    expect_identical(attr(terms(fm$formula_exposed), "term.labels"),"splines::ns(D_obs, df = 4)")
})



test_that("resolve_hazard_formulas() accepts strings without ~", {
    covars <- c("x1", "x2")
    expo   <- "D_obs"
    fm <- resolve_hazard_formulas(
        formula_unexposed =  "x1 + I(x2^2)",
        formula_exposed   = "x1 + I(x2^2) + splines::ns(D_obs, df = 4)",
        covariates = covars,
        exposure_time = expo
    )
    
    expect_s3_class(fm$formula_unexposed, "formula")
    expect_s3_class(fm$formula_exposed, "formula")
})


test_that("resolve_hazard_formulas() accepts strings with ~ and strips LHS", {
    covars <- c("x1", "x2")
    expo   <- "D_obs"
    
    fm <- resolve_hazard_formulas("~ x1 + x2", NULL, covars, expo)
    expect_s3_class(fm$formula_unexposed, "formula")
    expect_identical(attr(terms(fm$formula_unexposed), "response"), 0L)
    
    fm <- suppressWarnings(resolve_hazard_formulas("y ~ x1 + x2", NULL, covars, expo))
    expect_s3_class(fm$formula_unexposed, "formula")
    expect_identical(attr(terms(fm$formula_unexposed), "response"), 0L)
})


test_that("resolve_hazard_formulas() propogates errors from string_to_formula", {
    covars <- c("x1", "x2")
    expo   <- "D_obs"
    
    expect_error(
        resolve_hazard_formulas(NULL, "", covars, expo),
        regexp = "invalid", fixed = FALSE
    )
})

test_that("resolve_hazard_formulas() accepts covariate vectors", {
    covars <- c("x1", "x2")
    expo   <- "D_obs"
    
    fm <- resolve_hazard_formulas(
        formula_unexposed = c("x1", "x2"),
        formula_exposed   = c("x1", "I(x2 == 1)"),
        covariates = covars,
        exposure_time = expo
    )
    expect_s3_class(fm$formula_unexposed, "formula")
    expect_setequal(attr(terms(fm$formula_unexposed), "term.labels"), covars)
    
    expect_s3_class(fm$formula_exposed, "formula")
    expect_setequal(attr(terms(fm$formula_exposed), "term.labels"), c("x1", "I(x2 == 1)"))
})


test_that("resolve_hazard_formulas() accepts formulas", {
    covars <- c("x1", "x2")
    expo   <- "D_obs"
    
    expect_no_error(
        fm <- resolve_hazard_formulas(
            formula_unexposed =  as.formula("~ x1 + I(x2^2)"),
            formula_exposed   = as.formula(" ~ x1 + I(x2^2) + splines::ns(D_obs, df = 4)"),
            covariates = covars,
            exposure_time = expo
        )
    )
})

test_that("resolve_hazard_formulas() returns formulas without LHS", {
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
    # both should be formulas and one-sided (length 2)
    expect_s3_class(fm$formula_unexposed, "formula")
    expect_s3_class(fm$formula_exposed, "formula")
    
    # no response in either
    expect_identical(attr(terms(fm$formula_unexposed), "response"), 0L)
    expect_identical(attr(terms(fm$formula_exposed),   "response"), 0L)
})



# Tests for validate_formulas ---------------------------------------------

test_that("validate_formulas() passes on clean input", {
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

test_that("validate_formulas() warns for missing covariates in unexposed model", {
    covars <- c("x1", "x2")
    expo   <- "D_obs"
    expect_warning(
        validate_formulas(
            formula_unexposed =  as.formula("~  I(x2^2)"),
            formula_exposed   = as.formula("~ x1 + I(x2^2) + splines::ns(D_obs, df = 4)"),
            covariates = covars,
            exposure_time = expo),
        regexp = "does not include"
    )
})

test_that("validate_formulas() warns for missing covariates/time to exposure in exposed model", {
    covars <- c("x1", "x2")
    expo   <- "D_obs"
    expect_warning(
        validate_formulas(
            formula_unexposed =  as.formula("~ x1 + I(x2^2)"),
            formula_exposed   = as.formula("~ I(x2^2) + splines::ns(D_obs, df = 4)"),
            covariates = covars,
            exposure_time = expo),
        regexp = "does not include"
    )
    
    expect_warning(
        validate_formulas(
            formula_unexposed =  as.formula("~ x1 + I(x2^2)"),
            formula_exposed   = as.formula("~ x1 + I(x2^2)"),
            covariates = covars,
            exposure_time = expo),
        regexp = "does not include"
    )
})

test_that("validate_formulas() throws error for extra covariates in unexposed model", {
    expect_error(
        validate_formulas(~ x1 + x_extra, ~ x1 + x2 + D_obs, c("x1","x2"), "D_obs"),
        "`formula_unexposed` contains terms that are not included in `covariates`.*x_extra"
    )
})

test_that("validate_formulas() throws error for extra covariates in exposed model", {
    expect_error(
        validate_formulas(~ x1 + x2, ~ x1 + x2 + D_obs + x_extra, c("x1","x2"), "D_obs"),
        "`formula_exposed` contains terms that are not included in `covariates`.*x_extra"
    )
})

test_that("validate_formulas() works with NULL covariates and intercept-only formulas", {
    expect_no_error(
        validate_formulas(~ 1, ~ D_obs, covariates = NULL, exposure_time = "D_obs")
    )
})




