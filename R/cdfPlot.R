#' Create deluxe CDF Plots
#'
#' CDF plots are a good complement to p-value histograms as a way to evaluate
#' model performance and examine support for differential expression. Results
#' are ranked by p-value on the x-axis and the p-value plotted on the y-axis.
#' Since p-value distributions should be flat, this type of plot should produce a
#' straight line.  Any observations that fail to meet the null hypothesis will
#' appear as a break in the line at the low end of the curve.
#'
#' This function is designed to take topTable dataframes and display the
#' corresponding CDF plot. Data for the p-values below 0.1 (configurable via
#' pvalMax argument) are show in a full size plot. An inset figure shows the
#' whole p-value scale and highlights the portion shown in the full plot. Points
#' below 0.01 are a different color by default (threshold set by pThreshold
#' argument; shape/color attributes customizable through other arguments).
#'
#' The output can be controlled using printPlot = TRUE, which outputs the compound plot
#' to the console/knitr.
#'
#' \strong{Data Structure for the input dataframe:}
#'
#' The defaults are set for dataframes produced by topTable.  The column
#' "P.Value" is used by default to accommodate the column names used in topTable
#' dataframes.  Any other dataframe can be used with by explicitly defining the
#' p-value column name with the appropriate argument.
#'
#' Sensible defaults are chosen for symbols (Size, Shape, Color, and Fill).
#' There are optional arguments that allow these to be adjusted. A length of 2
#' is required for these arguments which applies the attributes in
#' this order: Significant, Not Significant.
#'
#' @param contrastDF A dataframe with LogRatio and LogIntensity columns and optionally a p-value or FDR column.
#' @param pvalCol Name of the p-value or FDR column (Default = "P.Value")
#' @param pvalMax Limit the range of the main plot (Default = 0.10)
#' @param pThreshold Used to color points (default = 0.01)
#' @param xlab X axis label (default is "LogIntensity column name "Rank")
#' @param ylab Y axis label (default is p-value column name)
#' @param title Plot title (Optional)
#' @param insetTitle Title for the inset plot (Optional)
#' @param symbolSize Size of symbols for up, no change, and down. Default = c(4, 1, 4);
#'        Note: All three cannot be the same size. Decimal values are acceptable to help offset that
#'        (e.g. 4, 4.1, 4.2).
#' @param symbolShape Shape of the symbols for up, no change, and down; Default = c(21, 20, 21)
#'        (20 = filled circle, 21 = fillable open circle); Note: The same symbol shape cannot
#'        be selected for all three symbols.
#'        See \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types}
#' @param symbolColor c(Up, NoChange, Down); default = c("black", "grey25", "grey0")
#'        See \url{http://research.stowers-institute.org/efg/R/Color/Chart}
#'        Note: Colors cannot be duplicated.
#' @param symbolFill Set the fill color for the symbols. Note only symbols 21-25 are fillable.
#'        This will have no effect on other symbols.
#'        Default = c("red3", "grey25", "deepskyblue4")
#'        Note: Colors cannot be duplicated.
#' @param alpha Controls the transparency of the plotted points (0-1; Default =
#'   0.5)
#' @param referenceLine Color for an horizontal line drawn at the p-threshold
#'   (Default = NULL; NULL disables, set to desired color to enable)
#' @param refLineThickness Set thickness of the reference line (Default = 1)
#' @param legendPosition (Default = "se")
#' @param baseFontSize The smallest size font in the figure in points. (Default
#'   = 12)
#' @param viewportX (Default = 0.15)
#' @param viewportY (Default = 0.93)
#' @param viewportWidth (Default = 0.35)
#' @param themeStyle "bw" or "grey" which correspond to theme_bw or theme_grey
#'   respectively. (Default = "grey")
#' @param printPlot Specify printing the combined plot to the console/knitr
#'   (Default = TRUE)
#' @param footnote Optional string placed right justified at bottom of plot.
#' @param footnoteSize Applies to footnote. (Default = 3)
#' @param footnoteColor Applies to footnote. (Default = "black")
#' @param footnoteJust Value 0 - 1. 0 is left justified, 1 is right justified, 0.5 is centered. (Default = 1)
#'
#' @return A list containing main plot, inset plot, and viewport. The plot can be
#'    reconstructed with print(myList$main); print(inset, vp = myList$viewport)
#'
#' @examples
#' \dontrun{
#'    # Plot to console (contrastDF is a topTable dataframe)
#'    cdfPlot(contrastDF, title = "My CDF Plot")
#' }
#' @import ggplot2 magrittr
#' @importFrom grid viewport
#' @importFrom dplyr arrange left_join
#' @importFrom assertthat assert_that
#'
#' @export
cdfPlot <- function(contrastDF,
                    pvalCol = "P.Value",
                    pThreshold = 0.01,
                    xlab = NULL, ylab = NULL,
                    title = NULL, insetTitle = NULL,
                    symbolSize = c(2, 1),
                    symbolShape = c(20, 20),
                    symbolColor = c("red3", "deepskyblue4"),
                    symbolFill = c("red3", "deepskyblue4"),
                    alpha = 1,
                    referenceLine = NULL,
                    refLineThickness = 1,
                    legendPosition = "se",
                    baseFontSize = 12,
                    viewportX = 0.15,
                    viewportY = 0.93,
                    viewportWidth = 0.35,
                    themeStyle = "grey",
                    pvalMax = 0.10,
                    printPlot = TRUE,
                    footnote,
                    footnoteSize = 3,
                    footnoteColor = "black",
                    footnoteJust = 1
) {

    assertthat::assert_that(pvalCol %in% colnames(contrastDF),
                            msg = "Specified pvalCol not found in the supplied dataframe (contrastDF).")

    if (!missing(symbolSize) || !missing(symbolShape) || !missing(symbolColor) || !missing(symbolFill)) {
        assertthat::assert_that(!length(symbolSize) == 2,
                                !length(symbolShape) == 2,
                                !length(symbolColor) == 2,
                                !length(symbolFill) == 2,
                                msg = "All specified symbol arguments must be of length 2, including symbolSize, symbolShape, symbolColor, and symbolFill.")
    }

    names(symbolShape) <- c("Significant", "Not Significant")
    names(symbolSize)  <- c("Significant", "Not Significant")
    names(symbolColor) <- c("Significant", "Not Significant")
    names(symbolFill)  <- c("Significant", "Not Significant")

    ssc <- data.frame(group = c("Significant", "Not Significant"),
                      symbolShape = symbolShape,
                      symbolSize = symbolSize,
                      symbolColor = symbolColor,
                      symbolFill = symbolFill,
                      order = c(1, 2),
                      stringsAsFactors = FALSE) %>%
        dplyr::arrange(order) # Setting the order defines the legend order

    # Columns to plot
    # Capture the labels from the colname
    xlabel <- "Rank"
    ylabel <- pvalCol
    x <- xlabel
    y <- ylabel

    # Combo PLOT: full data inset, most significant data in main plot
    # Rank by p-value
    contrastDF %<>% dplyr::arrange(!!sym(pvalCol))
    contrastDF$Rank <- c(1:nrow(contrastDF))

    # Let's plot the p-value subsets
    Sig <- contrastDF[[pvalCol]] <= pThreshold
    NotSig <- contrastDF[[pvalCol]] > pThreshold

    # Create group factor column in contrastDF
    contrastDF$group <- NA
    contrastDF$group[Sig] <- "Significant"
    contrastDF$group[NotSig] <- "Not Significant"
    contrastDF %<>% dplyr::left_join(ssc)
    contrastDF$group %<>% factor(levels = c("Significant", "Not Significant"))

    # Set an order field to force plotting of NotSig first
    contrastDF$order <- NA
    contrastDF$order[NotSig] <- 1
    contrastDF$order[Sig] <- 2

    # Rows to include in the zoomed in plot
    subsetRows <- sum(contrastDF[[y]] <= pvalMax)
    contrastDFsubset <- contrastDF[1:subsetRows,]

    # Plot subset percent of the data for the main plot
    cdfMain <- ggplot(contrastDFsubset, aes_string(x = x, y = y)) +
        aes(shape = group, size = group, color = group, fill = group, order = order) +
        # Scale lines tell it to use the actual values, not treat them as factors
        scale_shape_manual(name = "Group", guide = "legend", labels = ssc$group, values = ssc$symbolShape) +
        scale_size_manual(name = "Group", guide = "legend", labels = ssc$group, values = ssc$symbolSize) +
        scale_color_manual(name = "Group", guide = "legend", labels = ssc$group, values = ssc$symbolColor) +
        scale_fill_manual(name = "Group", guide = "legend", labels = ssc$group, values = ssc$symbolFill) +
        geom_point(alpha = alpha)

    # Optional Decorations
    if (!is.null(referenceLine)) {
        cdfMain <- cdfMain +
            geom_hline(yintercept = pThreshold, color = referenceLine,
                       size = refLineThickness, alpha = 0.5)
    }

    # Add Labels
    if (is.null(xlab)) { # Use colname unless supplied as argument
        cdfMain <- cdfMain + xlab(xlabel)
    } else {
        cdfMain <- cdfMain + xlab(xlab)
    }
    if (is.null(ylab)) {
        cdfMain <- cdfMain + ylab(ylabel)
    } else {
        cdfMain <- cdfMain + ylab(ylab)
    }
    if (!is.null(title)) {
        cdfMain <- cdfMain +
            ggtitle(title)
    }

    # Set the font size before placing the legend
    if (tolower(themeStyle) == "bw") {
        cdfMain <- cdfMain + theme_bw(baseFontSize)
    } else {
        cdfMain <- cdfMain + theme_grey(baseFontSize)
    }

    cdfMain <- setLegendPosition(cdfMain, legendPosition, themeStyle)

    if (!missing(footnote)) {
        cdfMain <- addFootnote(cdfMain,
                               footnoteText = footnote,
                               footnoteSize = footnoteSize,
                               footnoteColor = "black",
                               footnoteJust = footnoteJust)
    }

    # Set up the inset plot with All Data
    cdfInset <- ggplot(contrastDF, aes_string(x = x, y = y)) +
        aes(shape = group, size = group,
            color = group, fill = group,
            order = order) +
        # Scale lines tell it to use the actual values, not treat them as factors
        scale_shape_manual(name = "Group", guide = "none", labels = ssc$group,
                           values = ssc$symbolShape) +
        scale_size_manual(name = "Group", guide = "none", labels = ssc$group,
                          values = ssc$symbolSize) +
        scale_color_manual(name = "Group", guide = "none", labels = ssc$group,
                           values = ssc$symbolColor) +
        scale_fill_manual(name = "Group", guide = "none", labels = ssc$group,
                          values = ssc$symbolFill) +
        geom_rect(xmin = 0, xmax = subsetRows,
                  ymin = 0, ymax = max(dfsubset[[y]]), color = "lightblue",
                  fill = "lightblue", alpha = 0.2) +
        geom_point(alpha = alpha)

    # Add Labels
    if (is.null(xlab)) { # Use colname unless supplied as argument
        cdfInset <- cdfInset + xlab(xlabel)
    } else {
        cdfInset <- cdfInset + xlab(xlab)
    }
    if (is.null(ylab)) {
        cdfInset <- cdfInset + ylab(ylabel)
    } else {
        cdfInset <- cdfInset + ylab(ylab)
    }
    if (!is.null(insetTitle)) {
        cdfInset <- cdfInset +
            ggtitle(insetTitle)
    }

    # Adjust font size for inset.
    factor <- 4/baseFontSize

    # Set the font size
    if (tolower(themeStyle) == "bw") {
        cdfInset <- cdfInset + theme_bw() + baseTheme(baseFontSize*factor)
    } else {
        cdfInset <- cdfInset + theme_grey() + baseTheme(baseFontSize*factor)
    }

    # Adjust viewport Y if main Title present
    vy <- viewportY
    if (!is.null(title)) {
        adjust <- (baseFontSize/150)
        vy <- viewportY - adjust
    }
    # A viewport taking up a fraction of the plot area (upper left)
    vp <- grid::viewport(width = viewportWidth, height = viewportWidth,
                         x = viewportX, y = vy,
                         just = c("left", "top"))

    if (printPlot == TRUE) {
        print(cdfMain)
        print(cdfInset, vp = vp)
    }

    MyList <- list(main = cdfMain, inset = cdfInset, viewport = vp)
    return(MyList)
}
