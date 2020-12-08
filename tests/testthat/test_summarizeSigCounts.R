context("DGEtools - tests for summarizeSigCounts.R functions")


test_that('summarizeSigCounts.R: summarizeSigCounts()', {
    suppressWarnings(skip_if(is.null(getType(t_obj1, "topTable"))))

    myTopTables <- getType(t_obj1, "topTable")
    summarizedSigCounts <- summarizeSigCounts(myTopTables)

    expect_true(is.matrix(summarizedSigCounts))
    expect_equal(nrow(summarizedSigCounts), 5)
    expect_equal(ncol(summarizedSigCounts), 2)
    expect_equal(rownames(summarizedSigCounts), c("BMTL", "BMTH", "BMTM", "ena", "sham"))
    expect_equal(colnames(summarizedSigCounts), c("P.Value", "adj.P.Val"))

    # specify fcThreshold
    summarizedSigCounts_two <- summarizeSigCounts(myTopTables,
                                                  fcThreshold = 0.1)

    expect_true(is.matrix(summarizedSigCounts_two))
    expect_equal(nrow(summarizedSigCounts_two), 5)
    expect_equal(ncol(summarizedSigCounts_two), 2)
    expect_equal(rownames(summarizedSigCounts_two), c("BMTL", "BMTH", "BMTM", "ena", "sham"))
    expect_equal(colnames(summarizedSigCounts_two), c("P.Value", "adj.P.Val"))

    expect_error(summarizeSigCounts(myTopTable),
                 regexp = "object 'myTopTable' not found")
})


test_that('summarizeSigCounts.R: incorrect usage', {
    expect_error(summarizeSigCounts(),
                 regexp = "argument \"contrastList\" is missing, with no default")
})
