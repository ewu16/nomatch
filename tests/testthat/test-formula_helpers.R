# Tests for string_to_formula ---------------------------------------------

test_that("valid RHS returns a one-sided formula", {
    f <- string_to_formula("a + b", label = "0")
    expect_s3_class(f, "formula")
    # one-sided formula has length 2; [[2]] is RHS
    expect_equal(length(f), 2L)
    # RHS variable symbols
    expect_setequal(all.vars(f[[2]]), c("a", "b"))
    
    # terms() should indicate no response
    expect_identical(attr(terms(f), "response"), 0L)
})

test_that("whitespace is tolerated", {
    f <- string_to_formula("   a   +   b   ", label = "f")
    expect_s3_class(f, "formula")
    expect_setequal(all.vars(f[[2]]), c("a", "b"))
})


test_that("numeric-only RHS (e.g., '1') is allowed", {
    f <- string_to_formula("1", label = "f")
    expect_s3_class(f, "formula")
    # There are no variables on RHS, only a constant
    expect_identical(all.vars(f[[2]]), character(0))
    # Still a one-sided formula
    expect_equal(length(f), 2L)
})

test_that("functions and arithmetic inside RHS are handled", {
    f <- string_to_formula("I(b*100) + log(a + 1)", label = "f")
    expect_s3_class(f, "formula")
    expect_setequal(all.vars(f[[2]]), c("a", "b"))
})

test_that("rejects non-character inputs", {
    expect_error(string_to_formula(123, label = "0"),
                 regexp = "`0`?.*must be a single string",
                 fixed = FALSE)
    expect_error(string_to_formula(TRUE, label = "1"),
                 regexp = "`1`?.*must be a single string",
                 fixed = FALSE)
})

test_that("rejects character vectors with length != 1", {
    expect_error(string_to_formula(c("a + b", "c"), label = "0"),
                 regexp = "must be a single string", fixed = FALSE)
})

test_that("rejects strings that include a tilde", {
    expect_error(string_to_formula("y ~ x", label = "0"),
                 regexp = "Do not include '~'.*`0`", fixed = FALSE)
    expect_error(string_to_formula("~ x", label = "1"),
                 regexp = "Do not include '~'.*`1`", fixed = FALSE)
})

test_that("rejects empty or syntactically invalid RHS", {
    # empty
    expect_error(string_to_formula("", label = "0"),
                 regexp = "appears to be invalid", fixed = FALSE)
    # syntactically invalid expression
    expect_error(string_to_formula("a +", label = "0"),
                 regexp = "appears to be invalid", fixed = FALSE)
})

