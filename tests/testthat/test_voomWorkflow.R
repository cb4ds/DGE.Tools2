context("DGEtools - tests for voomWorkflow.R functions")


test_that('voomWorkflow.R: voomWorkflow()', {
    dgeObj <- DGEobj1
    design <- getItem(dgeObj, "design")
    designMatrix <- model.matrix(~ 0 + ReplicateGroup, design)
    dgeObj <- addItem(dgeObj   = dgeObj,
                      item     = designMatrix,
                      itemName = "designMat",
                      itemType = "designMatrix")
    dgeObj <- rmItem(dgeObj = dgeObj, itemName = "DGEList")

    # testing ellipsisArgs countThreshold and fracThreshold
    voom_dgeObj <- voomWorkflow(dgeObj            = dgeObj,
                                formula           = "~ 0 + ReplicateGroup",
                                designMatrixName  = "DM",
                                proteinCodingOnly = TRUE,
                                countThreshold    = 10,
                                zfpkmThreshold    = -3)
    expect_s3_class(voom_dgeObj, "DGEobj")

    # testing Voom/lmFit
    dupcorBlock <- rep(1:6, 8)
    voom_dgeObj <- voomWorkflow(dgeObj            = dgeObj,
                                formula           = "~ 0 + ReplicateGroup",
                                designMatrixName  = "DM",
                                proteinCodingOnly = TRUE,
                                dupCorBlock       = dupcorBlock,
                                countThreshold    = 10,
                                tpmThreshold      = 5)
    expect_s3_class(voom_dgeObj, "DGEobj")

    # testing assert statements
    expect_error(voomWorkflow(dgeObj = NULL),
                 regexp = "dgeObj must be specified and be of class 'DGEobj'.")
    expect_error(voomWorkflow(dgeObj = dgeObj),
                 regexp = "formula must be specified and be of class 'character'.")
    expect_error(voomWorkflow(dgeObj = dgeObj, formula = "~ 0 + ReplicateGroup"),
                 regexp = "designMatrixName must be specified and be of class 'character'.")
    expect_error(voomWorkflow(dgeObj           = dgeObj,
                              formula          = "~ 0 + ReplicateGroup",
                              designMatrixName = "DM",
                              countThreshold   = 10,
                              tpmThreshold     = 5,
                              zfpkmThreshold   = 3),
                 regexp = "No more than 2 intensity filtering arguments are allowed.")
    expect_error(voomWorkflow(dgeObj           = dgeObj,
                              formula          = "~ 0 + ReplicateGroup",
                              designMatrixName = "DM",
                              tpmThreshold     = 5,
                              zfpkmThreshold   = 3),
                 regexp = "When two intensity filtering arguments are supplied, one must be 'countThreshold'.")
})
