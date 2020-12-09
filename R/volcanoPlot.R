#' Create volcano plot
#'
#' A volcano plot shows Log Ratio data on the X axis and Negative Log P-values (NLP) on the
#' Y axis. This function is intended to show the volcano plot from a dataframe
#' created by topTable or topTreat. Properly normalized data will generally be
#' centered around LogRatio = 0.
#'
#' By default, the plot places "logFC" on the X axis and Log10 of the "P.Value" on the Y axis.
#' By default, a reference vertical line is drawn at LogRatio = 0 on the X axis.
#' Optionally, additional reference lines will be drawn at +/- a user supplied Log Ratio threshold.
#' The points are color coded using both the significance and fold-change thresholds supplied by the user.
#' By default, the P.Value field is used with a threshold of 0.01 to color code the points and fold-change
#' threshold of +/- 1.5X.
#'
#' \strong{Data Structure for the input dataframe:}
#'
#' The defaults are set for dataframes produced by topTable and topTreat.  The columns named "logFC"
#' and "P.Value" are used by default to accommodate the column
#' names used in topTable/topTreat dataframes.  Any other dataframe
#' can be used with fold-change, intensity, and significance measures, with appropriate
#' arguments to define the column names to use provided. By default, the
#' column names will be used for the axis labels, but can be overridden with xlab and ylab arguments.
#'
#' A significance measure (which defaults to P.Value <= 0.01) and LogRatio
#' threshold are used to color code genes that are significantly increased or decreased.
#' Use the appropriate arguments to use an FDR measure instead of p-value.
#'
#' Sensible defaults are chosen for symbols (Size, Shape, Color, and Fill), but they can be
#' adjusted through the use of optional arguments. A length of 3 is
#' required for these arguments which applies the attributes in this order:
#' Increased, NoChange, Decreased.
#'
#' @param contrastDF A dataframe with LogRatio and LogIntensity columns and optionally a
#'   p-value or FDR column (typically a topTable dataframe).
#' @param logRatioCol Name of the LogRatio column (Default = "logFC")
#' @param logIntCol Name of the LogIntensity column (Default = "AveExpr")
#' @param pvalCol Name of the p-value or FDR column (Default = "P.Value")
#' @param xlab X axis label (Default is the LogIntensity column name)
#' @param ylab Y axis label (Default is the LogRatio column name)
#' @param title Plot title (optional)
#' @param pthreshold Used to color points (Default = 0.01)
#' @param geneSymLabels A character vector of gene to label (must be the name space of the column
#'   specified by geneSymCol)
#' @param geneSymCol Name of the gene symbol column in contrastDF.  The gene symbol is
#'    not in topTable output by default so the user has to bind this column
#'    to the dataframe in advance.  This column will be used to label
#'    significantly changed points.
#' @param pthresholdLine Color for a horizontal line at the p-threshold (Default
#'   = NULL (disabled))
#' @param symbolSize Size of symbols for Up, no change, and Down. default = c(4,
#'        3.99, 4); Note: All three cannot be the same size. Decimal values are acceptable to help offset that
#'        (e.g. 4, 4.1, 4.2).
#' @param symbolShape Shape of the symbols for Up, no change, and Down; Default =
#'        c(21, 1, 21) (1 = open circle, 21 = fillable open circle); Note: The same symbol shape cannot
#'        be selected for all three symbols. See
#'        \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types}
#' @param symbolColor c(Up, NoChange, Down); default = c("black", "grey25",
#'   "grey0") See \url{http://research.stowers-institute.org/efg/R/Color/Chart}
#'   Note: Colors cannot be duplicated.
#' @param symbolFill Set the fill color for the symbols. Note only symbols 21-25
#'   are fillable. This will have no effect on other symbols. Default =
#'   c("red3", "grey25", "deepskyblue4") Note: Colors cannot be duplicated.
#' @param alpha Controls the transparency of the plotted points (range: 0-1;
#'   default = 0.7)
#' @param sizeByIntensity If TRUE, creates a column to support sizeByIntensity. (Default = TRUE)
#' @param foldChangeLines Position of reference vertical lines for fold change
#'   (Default = log2(1.5); NULL disables)
#' @param legendPosition One of "top", "bottom", "left", "right", "ne", "se",
#'   "nw", "sw", NULL. top/bottom/left/right place the legend outside the
#'   figure.  ne/se/nw/sw place the figure inside the figure. NULL disables the
#'   legend. Default = "right"
#' @param rugColor Specify color for rug density plot along x and y axes. Set to NULL
#'   to disable rug layer. (Default = NULL)
#' @param rugAlpha Sets the transparency for the rug layer.  An alpha <1.0 can take
#'   a long time to draw, so the default = 1.0.
#'   For a final plot try 0.1 or 0.2 for alpha for a more informative rug.
#' @param baseFontSize The smallest size font in the figure in points. Default =
#'   12
#' @param themeStyle "bw" or "grey" which correspond to theme_bw or theme_grey
#'   respectively. Default = bw"
#' @param refLineThickness Set the thickness for all reference lines (Default =
#'   1)
#' @param footnote Optional string placed right justified at bottom of plot.
#' @param footnoteSize Applies to footnote. (Default = 3)
#' @param footnoteColor Applies to footnote. (Default = "black")
#' @param footnoteJust Value 0 - 1. 0 is left justified, 1 is right justified, 0.5 is centered. (Default = 1)
#'
#' @return ggplot object containing volcano plot
#'
#' @examples
#' \dontrun{
#'    # Simple plot with custom title (contrastDF is a topTable dataframe)
#'    myPlot <- volcanoPlot(contrastDF, title = "Plot Title")
#'
#'    # Some options with a custom datafile
#'    myPlot <- volcanoPlot(contrastDF,
#'                          pthreshold = 0.1,
#'                          logRatioCol = "Log2ratio",
#'                          logIntCol = "AverageIntensity",
#'                          pvalCol = "BHFDR",
#'                          xlab = "Log2 Ratio", ylab = "Log10 BHFDR",
#'                          title = "Profile Plot Title",
#'                          referenceLine = "blue",
#'                          legendPosition = "ne")
#' }
#'
#' @import ggplot2 magrittr
#' @importFrom dplyr left_join
#' @importFrom ggrepel geom_text_repel
#' @export
volcanoPlot <- function(contrastDF,
                        logRatioCol = "logFC",
                        logIntCol = "AveExpr",
                        pvalCol = "P.Value",
                        pthreshold = 0.01,
                        geneSymLabels,
                        geneSymCol,
                        xlab = NULL, ylab = NULL, title = NULL,
                        symbolSize = c(4, 3.999, 4),
                        symbolShape = c(21, 1, 21),
                        symbolColor = c("black", "grey25", "grey0"),
                        symbolFill = c("red3", "grey25", "deepskyblue4"),
                        alpha = 0.5,
                        sizeByIntensity = TRUE,
                        pthresholdLine = NULL,
                        foldChangeLines = log2(1.5),
                        refLineThickness = 1,
                        legendPosition = "right",
                        rugColor = NULL,
                        rugAlpha = 1.0,
                        baseFontSize = 12,
                        themeStyle = "grey",
                        footnote,
                        footnoteSize = 3,
                        footnoteColor = "black",
                        footnoteJust = 1) {

    # Make sure specified columns exist
    assertthat::assert_that(logRatioCol %in% colnames(contrastDF),
                            msg = "logRatioCol column not found in contrastDF.")
    assertthat::assert_that(logIntCol %in% colnames(contrastDF),
                            msg = "logIntCol column not found in contrastDF.")
    assertthat::assert_that(pvalCol %in% colnames(contrastDF),
                            msg = "pvalCol column not found in contrastDF.")
    if (!missing(geneSymCol)) {
        assertthat::assert_that(geneSymCol %in% colnames(contrastDF),
                                msg = "geneSymol column not found in contrastDF.")
    }
    if (!missing(symbolSize) || !missing(symbolShape) || !missing(symbolColor) || !missing(symbolFill)) {
        assertthat::assert_that(!length(symbolSize) == 3,
                                !length(symbolShape) == 3,
                                !length(symbolColor) == 3,
                                !length(symbolFill) == 3,
                                msg = "All specified symbol arguments must be of length 3, including symbolSize, symbolShape, symbolColor, and symbolFill.")
    }

    if (sizeByIntensity == TRUE) {
        # Create a column to support sizeByIntensity
        contrastDF$LogInt = contrastDF[[logIntCol]]
        # Set a floor and a ceiling
        contrastDF$LogInt[contrastDF$LogInt < 0] = 0
        contrastDF$LogInt[contrastDF$LogInt > 10] = 10
    }

    names(symbolShape) = c("Increased", "No Change", "Decreased")
    names(symbolSize)  = c("Increased", "No Change", "Decreased")
    names(symbolColor) = c("Increased", "No Change", "Decreased")
    names(symbolFill)  = c("Increased", "No Change", "Decreased")

    ssc = data.frame(group = c("Increased", "No Change", "Decreased"),
                     symbolShape = symbolShape,
                     symbolSize = symbolSize,
                     symbolColor = symbolColor,
                     symbolFill = symbolFill,
                     order = c(1,3,2),
                     stringsAsFactors = FALSE) %>% arrange(order)

    # Capture the labels from the colname
    xlabel = logRatioCol
    ylabel = paste("-log10(", pvalCol, ")", sep = "")
    # Now make the columnames suitable for use with aes_string
    x = make.names(colnames(contrastDF)[colnames(contrastDF) == logRatioCol])
    colnames(contrastDF)[colnames(contrastDF) == logRatioCol] = make.names(colnames(contrastDF)[colnames(contrastDF) == logRatioCol])
    # Make a log10significance column and make that the y column
    contrastDF$NegativeLogP = -log10(contrastDF[,pvalCol])
    y = "NegativeLogP"

    # DELUXE PLOT: plot groups in different colors/shapes
    # Let's plot the subsets
    DEup = contrastDF[[pvalCol]] <= pthreshold & contrastDF[[logRatioCol]] > 0
    DEdn = contrastDF[[pvalCol]] <= pthreshold & contrastDF[[logRatioCol]] < 0
    DEnot = !DEup & !DEdn
    # Create group factor column in contrastDF
    contrastDF$group = NA
    contrastDF$group[DEup] = "Increased"
    contrastDF$group[DEdn] = "Decreased"
    contrastDF$group[DEnot] = "No Change"
    contrastDF %<>% dplyr::left_join(ssc)
    contrastDF$group %<>% factor(levels = c("Increased", "Decreased", "No Change"))

    # Set an order field to force plotting of NoChange first
    contrastDF$order = NA
    contrastDF$order[DEup] = 1
    contrastDF$order[DEdn] = 1
    contrastDF$order[DEnot] = 0

    volcanoPlot <- ggplot(contrastDF, aes_string(x = x, y = y)) +
        aes(shape = group, size = group,
            color = group, fill = group,
            order = order) +
        # Scale lines tell it to use the actual values, not treat them as factors
        scale_shape_manual(name = "Group", guide = "legend", labels = ssc$group,
                           values = ssc$symbolShape) +
        scale_size_manual(name = "Group", guide = "legend", labels = ssc$group,
                          values = ssc$symbolSize) +
        scale_color_manual(name = "Group", guide = "legend", labels = ssc$group,
                           values = ssc$symbolColor) +
        scale_fill_manual(name = "Group", guide = "legend", labels = ssc$group,
                          values = ssc$symbolFill) +
        geom_point(alpha = alpha) +
        # Box around the legend
        theme(legend.background = element_rect(fill = "gray95", size = .5, linetype = "dotted"))

    # Optional Decorations
    if (!is.null(rugColor)) {
        volcanoPlot <- volcanoPlot + geom_rug(data = contrastDF, inherit.aes = FALSE,
                                              color = rugColor,
                                              alpha = rugAlpha,
                                              show.legend = FALSE,
                                              aes_string(x = x, y = y))
    }

    if (sizeByIntensity == TRUE) {
        volcanoPlot <- volcanoPlot + aes(size = LogInt) +
            scale_size_continuous()
    }

    if (!is.null(pthresholdLine)) {
        volcanoPlot <- volcanoPlot +
            geom_hline(yintercept = -log10(pthreshold), color = pthresholdLine,
                       alpha = 0.5, size = refLineThickness)
    }

    if (!is.null(foldChangeLines)) {
        volcanoPlot <- volcanoPlot +
            geom_vline(xintercept = foldChangeLines, color = symbolFill["Increased"],
                       alpha = 0.5, size = refLineThickness) +
            geom_vline(xintercept = -foldChangeLines, color = symbolFill["Decreased"],
                       alpha = 0.5, size = refLineThickness)
    }

    # Add geneSym labels to increased & decreased genes
    if (!missing(geneSymLabels) & !missing(geneSymCol)) {
        # Filter contrastDF to changed genes
        idx <- contrastDF[[geneSymCol]] %in% geneSymLabels
        contrastDFsubset <- contrastDF[idx,]
        volcanoPlot <- volcanoPlot +
            geom_text_repel(data = contrastDFsubset, aes_string(x = x, y = y, label = geneSymCol),
                            show.legend = FALSE)
    }

    # Add Labels
    if (is.null(xlab)) { # Use colname unless supplied as argument
        volcanoPlot <- volcanoPlot + xlab(xlabel)
    } else {
        volcanoPlot <- volcanoPlot + xlab(xlab)
    }
    if (is.null(ylab)) {
        volcanoPlot <- volcanoPlot + ylab(ylabel)
    } else {
        volcanoPlot <- volcanoPlot + ylab(ylab)
    }
    if (!is.null(title)) {
        volcanoPlot <- volcanoPlot + ggtitle(title)
    }

    # Set the font size before placing the legend
    if (tolower(themeStyle) == "bw") {
        volcanoPlot <- volcanoPlot + theme_bw() + baseTheme(baseFontSize)
    } else {
        volcanoPlot <- volcanoPlot + theme_grey() + baseTheme(baseFontSize)
    }

    volcanoPlot <- setLegendPosition(volcanoPlot, legendPosition, themeStyle)

    # Footnote
    if (!missing(footnote)) {
        volcanoPlot <- addFootnote(volcanoPlot,
                                   footnoteText = footnote,
                                   footnoteSize = footnoteSize,
                                   footnoteColor = "black",
                                   footnoteJust = footnoteJust)
    }

    return(volcanoPlot)
}
