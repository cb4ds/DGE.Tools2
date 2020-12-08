context("DGEtools - tests for volcanoPlot.R functions")


test_that("volcanoPlot.R: volcanoPlot()", {
    skip_if(is.null(t_obj1$RG_fit))

    df <- topTable(t_obj1$RG_fit, number = 100)
    volcano_plot <- volcanoPlot(df, logRatioCol = "adj.P.Val")
    expect_s3_class(volcano_plot, c("gg","ggplot"))

    df$GeneSym <- rep(c("sym1", "sym2","sym3","sym4"),nrow(df)/4)
    volcano_plot <- volcanoPlot(df                 = df,
                                title              = "Plot Title",
                                logRatioCol        = "adj.P.Val",
                                rugColor           = "red",
                                geneSymCol         = "GeneSym",
                                geneSymLabels      = c("sym1", "sym2"),
                                xlab               = "XLabel",
                                ylab               = "YLabel",
                                pthresholdLine     = "blue",
                                footnote           = "This is footnote",
                                themeStyle         = "bw")

    expect_s3_class(volcano_plot, c("gg","ggplot"))

    expect_error(volcanoPlot(df, logRatioCol = "xyz"),
                 regexp =  "logRatioCol column not found in df.")
    expect_error(volcanoPlot(df, logRatioCol = "adj.P.Val", logIntCol = "xyz"),
                 regexp =  "logIntCol column not found in df.")
    expect_error(volcanoPlot(df, logRatioCol = "adj.P.Val", pvalCol = "xyz"),
                 regexp = "pvalCol column not found in df.")
    expect_error(volcanoPlot(df, logRatioCol = "adj.P.Val", symbolSize = 2),
                 regexp = "All specified symbol arguments must be of length 3, including symbolSize, symbolShape, symbolColor, and symbolFill.")
})
