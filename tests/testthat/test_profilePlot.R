context("DGEtools - tests for profilePlot.R functions")


test_that("profilePlot.R: profilePlot()", {
    skip_if(!("RG_fit" %in% names(t_obj1)))

    contrastDF <- topTable(t_obj1$RG_fit, number = 100)
    profile_plot <- profilePlot(contrastDF, logRatioCol = "adj.P.Val")
    expect_s3_class(profile_plot, c("gg", "ggplot"))

    contrastDF$GeneSym <- rep(c("sym1", "sym2","sym3","sym4"),nrow(contrastDF)/4)
    profile_plot <- profilePlot(contrastDF         = contrastDF,
                                title              = "Plot Title",
                                logRatioCol        = "adj.P.Val",
                                sizeBySignificance = TRUE,
                                rugColor           = "red",
                                geneSymCol         = "GeneSym",
                                geneSymLabels      = c("sym1", "sym2"),
                                xlab               = "XLabel",
                                ylab               = "YLabel",
                                footnote           = "This is footnote",
                                themeStyle         = "bw")
    expect_s3_class(profile_plot, c("gg","ggplot"))

    expect_error(profilePlot(contrastDF, logRatioCol = "xyz"),
                 regexp =  "logRatioCol column not found in contrastDF.")
    expect_error(profilePlot(contrastDF, logRatioCol = "adj.P.Val", logIntCol = "xyz"),
                 regexp =  "logIntCol column not found in contrastDF.")
    expect_error(profilePlot(contrastDF, logRatioCol = "adj.P.Val", pvalCol = "xyz"),
                 regexp = "pvalCol column not found in contrastDF.")
    expect_error(profilePlot(contrastDF, logRatioCol = "adj.P.Val", geneSymCol = "xyz"),
                 regexp = "geneSymCol column not found in contrastDF.")
    expect_error(profilePlot(contrastDF, logRatioCol = "adj.P.Val", symbolSize = 2),
                 regexp = "All specified symbol arguments must be of length 3, including symbolSize, symbolShape, symbolColor, and symbolFill.")
})
