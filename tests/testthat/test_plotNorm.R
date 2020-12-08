context("DGEtools - tests for plotNorm.R functions")


test_that("plotNorm.R: plotNorm()", {

    # testing with DGEobj and plotType - box
    norm_plot <- plotNorm(t_obj1$counts, plotType = "box")
    expect_s3_class(norm_plot, c("gg", "ggplot"))

    # testing with count matrix and plotType - density
    norm_plot <- plotNorm(t_obj1, plotType = "density")
    expect_s3_class(norm_plot, c("gg", "ggplot"))

    # testing assert statements
    expect_error(plotNorm(NULL),
                 regexp = "dat must be of either class 'matrix' or 'DGEobj'.")
    expect_error(plotNorm(t_obj1, plotType = "heatmap"),
                 regexp = "plotType must be one of 'box' or 'density'.")
    expect_error(plotNorm(t_obj1, normalize = "xyz"),
                 regexp = "normalize must be one of 'TMM', 'RLE', 'upperquartile', or 'none'.")
})
