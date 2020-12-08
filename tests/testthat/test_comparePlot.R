context("DGEtools - tests for comparePlot.R functions")


test_that("comparePlot.R: comparePlot()", {
    suppressWarnings(skip_if(is.null(getType(t_obj1, "topTable"))))

    ttList <- getType(t_obj1, "topTable")[1:2]
    # Capture the default logFC and P.Value
    compareDat <- comparePrep(ttList)

    # Draw the plot
    cPlot <- comparePlot(compareDat)
    expect_s3_class(cPlot , c("gg", "ggplot"))

    # testing function with df without significance measures supplied
    cPlot <- comparePlot(compareDat[,1:2])
    expect_s3_class(cPlot , c("gg", "ggplot"))

    # testing aesthetics of plots
    cPlot <- comparePlot(compareDat,
                         title    = "MyPlot",
                         xlab     = "xaxis-title",
                         ylab     = "yaxis-title",
                         rugColor = "green",
                         footnote = "This is my footnote")
    expect_setequal(unlist(cPlot$labels[c("title","y", "x")]), c("MyPlot", "yaxis-title", "xaxis-title"))
    expect_setequal(cPlot$layers[[2]]$aes_params$colour, "green")

    # testing assert statement
    expect_error(comparePlot(compareDat[, 1, drop = FALSE]),
                 regexp = "Need at least two numeric columns in df.")
    expect_error(comparePlot(t_obj1$design),
                 regexp = "Need at least two numeric columns in df.")

    # failing function with Symbol argument length < 4
    expect_error(comparePlot(compareDat, symbolSize = 1),
                 regexp = "All specified symbol arguments must be of length 4, including symbolSize, symbolShape, symbolColor, and symbolFill.")
})

test_that("comparePlot.R: comparePrep()", {
    suppressWarnings(skip_if(is.null(getType(t_obj1, "topTable"))))

    ttList <- getType(t_obj1, "topTable")[1:2]
    # Capture the default logFC and P.Value
    compareDat <- comparePrep(ttList)
    expect_s3_class(compareDat,"data.frame")

    expect_error(comparePrep(ttList[[1]]),
                 regexp = "ttList must be a named list of length 2 where both items are of class 'data.frame'.")
    expect_error(comparePrep(ttList, valueCol = "P.val"),
                 regexp = "The valueCol must be included in the colnames of both items of ttList.")
    expect_error(comparePrep(ttList, significanceCol = "P.val"),
                 regexp = "The significanceCol must be included in the colnames of both items of ttList.")
    ttList_uncommon_ids <- list("BMTL" = ttList$BMTL[1:10,], "BMTH" = ttList$BMTH[21:30,])
    expect_error(comparePrep(ttList_uncommon_ids),
                 regexp = "No common gene IDs were found between the two dataframes in ttList.")
})
