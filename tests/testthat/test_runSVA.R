context("DGEtools - tests for runSVA.R functions")


test_that("runSVA.R: runSVA()", {
    skip_if(is.null(t_obj1$RG))

    dgeObj_sva <- runSVA(dgeObj = t_obj1, designMatrixName = "RG")

    expect_s3_class(dgeObj_sva, "DGEobj")

    expect_error(runSVA(designMatrixName = "designMatrix"),
                 regexp = "dgeObj must be specified, be of class 'DGEobj', and should have a 'design' attribute.")
    expect_error(runSVA(dgeObj = t_obj1),
                 regexp = "designMatrixName must be specified, should be of class 'character', and must exist as an attribute on the dgeObj.")
    expect_error(runSVA(dgeObj = t_obj1, designMatrixName = "RG", method =  "xyz"),
                 regexp = "method must be one of 'leek' or 'be'.")
})
