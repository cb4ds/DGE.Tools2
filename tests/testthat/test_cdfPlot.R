context("DGEtools - tests for cdfPlot.R functions")


test_that("cdfPlot.R: cdfPlot()", {
    skip_if(!("RG_fit" %in% names(t_obj1)))

    top_table <- topTable(t_obj1$RG_fit, number = 100)

    # testing plot with default values.
    plot <- cdfPlot(top_table, referenceLine = "blue")
    expect_s3_class(plot$main, c("gg", "ggplot"))
    expect_s3_class(plot$inset, c("gg", "ggplot"))
    expect_s3_class(plot$viewport, "viewport")

    # testing plot with customized aesthetics.
    plot_with_aes <- cdfPlot(top_table,
                             insetTitle    = "Sub plot title",
                             xlab          = "xaxis-title",
                             ylab          = "yaxis-title",
                             title         = "MyPlot",
                             themeStyle    = "bw",
                             referenceLine = "blue",
                             footnote      = "this is footnote of the plot")

    expect_s3_class(plot_with_aes$main, c("gg", "ggplot"))
    expect_s3_class(plot_with_aes$inset, c("gg", "ggplot"))
    expect_s3_class(plot_with_aes$viewport, "viewport")

    expect_setequal(unlist(plot_with_aes$main$labels[c("title", "y", "x")]), c("MyPlot", "yaxis-title", "xaxis-title"))
    expect_setequal(plot_with_aes$inset$labels$title, "Sub plot title")
    expect_equal(plot_with_aes$main$layers[[2]]$geom_params$colour, "blue")
    expect_equal(plot_with_aes$main$theme$strip.background$colour, "grey20")
    expect_equal(plot_with_aes$main$layers[[3]]$geom_params$label, "this is footnote of the plot")

    expect_error({cdfPlot(top_table,pvalCol = "p.value")},
                 regexp = "Specified pvalCol not found in the supplied dataframe (df).",
                 fixed = TRUE)
    expect_error({cdfPlot(top_table, symbolShape = 20)},
                 regexp = "All specified symbol arguments must be of length 2, including symbolSize, symbolShape, symbolColor, and symbolFill.")
})
