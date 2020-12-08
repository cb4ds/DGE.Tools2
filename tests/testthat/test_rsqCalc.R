context("DGEtools - tests for rsqCalc.R functions")


test_that("rsqCalc.R: rsqCalc()", {
    skip_if(is.null(t_obj1$DGEList))

    dgelist <- getItem(t_obj1, "DGEList")
    log2cpm <- cpm(dgelist, log = TRUE)
    rsq <- rsqCalc(log2cpm, t_obj1$RG_fit)
    expect_type(rsq, "double")

    expect_error(rsqCalc(1:10), regexp = "normMatrix must be of class 'data.frame' or 'matrix'.")
    expect_error(rsqCalc(log2cpm, 1:10), regexp = "fit must be of class 'MArrayLM'.")
    expect_error(rsqCalc(as.matrix(LETTERS), t_obj1$RG_fit), regexp = "All of the entries in normMatrix must be numeric.")
})
