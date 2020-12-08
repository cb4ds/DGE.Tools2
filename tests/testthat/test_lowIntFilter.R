context("DGEtools - tests for lowIntFilter.R functions")


test_that('lowIntFilter: lowIntFilter()', {
    skip_if(is.null(getItem(t_obj1, "geneData")$ExonLength))

    lowIntFilter_one_test <- lowIntFilter(t_obj1, countThreshold = 10)

    expect_s3_class(lowIntFilter_one_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_one_test$counts), 443)

    lowIntFilter_two_test <- lowIntFilter(t_obj1, zfpkmThreshold = -3.0)

    expect_s3_class(lowIntFilter_two_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_two_test$counts), 11699)

    lowIntFilter_three_test <- lowIntFilter(t_obj1, fpkThreshold = 5)

    expect_s3_class(lowIntFilter_three_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_three_test$counts), 11679)

    lowIntFilter_four_test <- lowIntFilter(t_obj1, countThreshold = 10, zfpkmThreshold = -3.0, sampleFraction = 0.75)

    expect_s3_class(lowIntFilter_four_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_four_test$counts), 11542)

    expect_error(lowIntFilter(lowIntFilter_five_test),
                 regexp = "object 'lowIntFilter_five_test' not found")
})


test_that('lowIntFilter: incorrect usage', {
    expect_error(lowIntFilter(),
                 regexp = "argument \"dgeObj\" is missing, with no default")

    expect_error(lowIntFilter(t_obj1, zfpkmThreshold = 3.0, tpmThreshold = 1),
                 regexp = "Must use zfpkmThreshold or tpmThreshold, but not both.")
})
