context("DGEtools - tests for plotPValHist.R functions")


test_that("plotPValHist.R: plotPvalHist()", {

    # testing plotPvalHist with savePlot and facet = TRUE
    pvalMatrix <- extractCol(getType(DGEobj1, "topTable"), "P.Value")
    pval_plot <- plotPvalHist(pvalMatrix, facetFontSize = 14)
    expect_s3_class(pval_plot, c("gg","ggplot"))

    # testing plotPvalHist with savePlot and facet = FALSE
    pval_plot <- plotPvalHist(as.matrix(pvalMatrix),
                              facet     = FALSE)
    expect_length(pval_plot, 5)
})
