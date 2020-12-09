context("DGEtools - tests for extractCol.R functions")

test_that('extractCol: extractCol()', {
    suppressWarnings(skip_if(is.null(getType(t_obj1, "topTable"))))

    extractCol_contrastList <- getType(t_obj1, "topTable")[1:2]
    extractCol_one_test <- extractCol(extractCol_contrastList, colName = "P.Value")

    expect_true(is.data.frame(extractCol_one_test))
    expect_equal(nrow(extractCol_one_test), 11709)
    expect_equal(ncol(extractCol_one_test), 2)
    expect_equal(names(extractCol_one_test), c("BMTL", "BMTH"))

    extractCol_two_test <- extractCol(extractCol_contrastList, colName = "P.Value", robust = FALSE)

    expect_true(is.data.frame(extractCol_two_test))
    expect_equal(nrow(extractCol_two_test), 11709)
    expect_equal(ncol(extractCol_two_test), 2)
    expect_equal(names(extractCol_two_test), c("BMTL", "BMTH"))

    expect_error(extractCol(extractCol_three_test),
                 regexp = "object 'extractCol_three_test' not found")
})


test_that('extractCol: incorrect usage', {
    expect_error(extractCol(),
                 regexp = "argument \"contrastList\" is missing, with no default")
})
