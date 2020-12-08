context("DGEtools - tests for tidyContrasts.R functions")


test_that('tidyContrasts: tidyContrasts()', {
    suppressWarnings(skip_if(is.null(getType(t_obj1, "topTable"))))

    tidyContrast_one_test <- tidyContrasts(t_obj1)

    expect_true(is.data.frame(tidyContrast_one_test))
    expect_equal(nrow(tidyContrast_one_test), 58545)
    expect_equal(ncol(tidyContrast_one_test), 10)

    tidyContrast_two_test <- tidyContrasts(t_obj1,
                                           rownameColumn = "rownameColumn")

    expect_true(is.data.frame(tidyContrast_two_test))
    expect_equal(nrow(tidyContrast_two_test), 58545)
    expect_equal(ncol(tidyContrast_two_test), 10)
    expect_true("rownameColumn" %in% names(tidyContrast_two_test))

    expect_warning(tidyContrasts(t_obj1,
                                 includeColumns = c("rownames", "logFC", "CI.L", "CI.R")),
                   regexp = "Some requested columns are not present in all dataframes.")

    expect_error(tidyContrasts(contrastObj_test),
                 regexp = "object 'contrastObj_test' not found")
})


test_that('tidyContrasts: incorrect usage', {
    expect_error(tidyContrasts(),
                 regexp = "argument \"contrastObj\" is missing, with no default")
})
