test_that("check_reserved_vars() passes when there are no conflicts", {
    expect_no_error(check_reserved_vars(vars = c("x1", "x2"), reserved_vars = c(".t0", ".estimate"), vars_label = "covariates"))
})

test_that("check_reserved_vars() errors when vars collide with reserved names", {
    reserved <- c(".t0", ".estimate", "group_id")
    expect_error(
        check_reserved_vars(vars = c("x1", "group_id"), reserved_vars = reserved, vars_label = "covariates"),
        "contain reserved variable names.*group_id", ignore.case = TRUE
    )
})

