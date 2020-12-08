context("DGEtools - tests for runIHW.R functions")


test_that('runIHW: runIHW()', {
    suppressWarnings(skip_if(is.null(getType(t_obj1, "topTable"))))

    runIHW_ttList <- getType(t_obj1, "topTable")[1:2]

    runIHW_test_one <- runIHW(runIHW_ttList)

    expect_true(is.list(runIHW_test_one))
    expect_equal(length(runIHW_test_one), 2)
    expect_equal(names(runIHW_test_one), c("contrasts", "ihwObj"))
    expect_s4_class(runIHW_test_one$ihwObj[[1]], "ihwResult")
    expect_s4_class(runIHW_test_one$ihwObj[[2]], "ihwResult")

    expect_error(runIHW(runIWH_test_two),
                 regexp = "object 'runIWH_test_two' not found")
})


test_that('runIHW: incorrect usage', {
    expect_error(runIHW(),
                 regexp = "argument \"contrastList\" is missing, with no default")
})
