context("DGEtools - tests for runVoom.R functions")


test_that('runVoom.R: runVoom()', {
   skip_if(is.null(t_obj1$DGEList))

    dgeObj <- t_obj1
    design <- getItem(dgeObj, "design")
    designMatrix <- model.matrix(~ 0 + ReplicateGroup, design)
    dgeObj <- addItem(dgeObj   = dgeObj,
                      item     = designMatrix,
                      itemName = "designMat",
                      itemType = "designMatrix")

    voom_dgeObj <- runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           runDupCorTwice   = TRUE,
                           qualityWeights   = FALSE,
                           mvPlot           = FALSE,
                           runEBayes        = TRUE,
                           robust           = TRUE,
                           proportion       = 0.01)
    expect_s3_class(voom_dgeObj, "DGEobj")
    expect_true("matrix" %in% class(voom_dgeObj$designMat))
    expect_s4_class(voom_dgeObj$designMat_fit, "MArrayLM")
    expect_s4_class(voom_dgeObj$designMat_Elist, "EList")

    # testing indQW analysis
    voom_dgeObj <- runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           runDupCorTwice   = TRUE,
                           qualityWeights   = TRUE,
                           mvPlot           = FALSE,
                           runEBayes        = FALSE,
                           robust           = TRUE,
                           proportion       = 0.01)
    expect_s3_class(voom_dgeObj, "DGEobj")
    expect_true("matrix" %in% class(voom_dgeObj$designMat))
    expect_s4_class(voom_dgeObj$designMat_fit, "MArrayLM")
    expect_s4_class(voom_dgeObj$designMat_Elist, "EList")

    # testing blockedQW analysis
    voom_dgeObj <- runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           runDupCorTwice   = TRUE,
                           qualityWeights   = TRUE,
                           var.design       = designMatrix,
                           mvPlot           = FALSE,
                           runEBayes        = TRUE,
                           robust           = TRUE,
                           proportion       = 0.01)
    expect_s3_class(voom_dgeObj, "DGEobj")
    expect_true("matrix" %in% class(voom_dgeObj$designMat))
    expect_s4_class(voom_dgeObj$designMat_fit, "MArrayLM")
    expect_s4_class(voom_dgeObj$designMat_Elist, "EList")

    # testing dupcor_indQW analysis
    dupcorBlock <- rep(1:6, 8)
    voom_dgeObj <- runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           runDupCorTwice   = TRUE,
                           dupCorBlock      = dupcorBlock,
                           qualityWeights   = TRUE,
                           mvPlot           = FALSE,
                           runEBayes        = TRUE,
                           robust           = TRUE,
                           proportion       = 0.01)
    expect_s3_class(voom_dgeObj, "DGEobj")
    expect_true("matrix" %in% class(voom_dgeObj$designMat))
    expect_s4_class(voom_dgeObj$designMat_fit, "MArrayLM")
    expect_s4_class(voom_dgeObj$designMat_Elist, "EList")

    # testing dupcor_base analysis
    voom_dgeObj <- runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           dupCorBlock      = dupcorBlock,
                           runDupCorTwice   = TRUE,
                           qualityWeights   = FALSE,
                           mvPlot           = FALSE,
                           runEBayes        = TRUE,
                           robust           = TRUE,
                           proportion       = 0.01)
    expect_s3_class(voom_dgeObj, "DGEobj")
    expect_true("matrix" %in% class(voom_dgeObj$designMat))
    expect_s4_class(voom_dgeObj$designMat_fit, "MArrayLM")
    expect_s4_class(voom_dgeObj$designMat_Elist, "EList")

    # testing dupcor_vdQW analysis
    voom_dgeObj <- runVoom(dgeObj           = dgeObj,
                           designMatrixName = "designMat",
                           runDupCorTwice   = TRUE,
                           dupCorBlock      = dupcorBlock,
                           qualityWeights   = TRUE,
                           var.design       = designMatrix,
                           mvPlot           = FALSE,
                           runEBayes        = TRUE,
                           robust           = TRUE,
                           proportion       = 0.01)
    expect_s3_class(voom_dgeObj, "DGEobj")
    expect_true("matrix" %in% class(voom_dgeObj$designMat))
    expect_s4_class(voom_dgeObj$designMat_fit, "MArrayLM")
    expect_s4_class(voom_dgeObj$designMat_Elist, "EList")

    # testing assert statements
    expect_error(runVoom(dgeObj = NULL),
                 regexp = "dgeObj must be specified and must be of class 'DGEobj'.")
    expect_error(runVoom(dgeObj = dgeObj, designMatrixName = "xyz"),
                 regexp = "designMatrixName must be specified and must be one of the items in dgeObj. Use names(dgeObj) to check for available options.",
                 fixed = TRUE)
})
