#' Plot log intensity versus log ratio
#'
#' A profile plot shows Log Intensity on the X axis and Log Ratio on the Y axis.
#' This function is intended to show the profile plot from a dataframe
#' created by topTable or topTreat.  Properly normalized data will generally be
#' centered around LogRatio = 0.
#'
#' By default, the plot places "logFC" on the Y axis and "AveExpr" on the X
#' axis. By default, a reference horizontal line is shown at LogRatio = 0 on the
#' Y axis. Optionally, additional reference lines will be drawn at +/- a user
#' supplied LogRatio threshold. A Loess line fit is drawn through the actual
#' data. The points are color coded using the significance measure (p-value or
#' FDR threshold) supplied by the user. By default, the P.Value field is used
#' with a threshold of 0.01 to color code the points.
#'
#' The profilePlot has also been implemented
#' with a scalable theme that makes it easy to scale font sizes larger
#' or smaller for PPT or knitr output.  Add baseFont(n),
#' where n equals the base font size (12 works well for knitr, 18 or 24 works well
#' for PPT).  e.g. myPlot + theme_grey(18).
#'
#' \strong{Data Structure for the input dataframe:}
#'
#' The defaults are set for dataframes produced by topTable and topTreat.  The columns named "logFC",
#' "AveExpr", and "P.Value" are used by default to accommodate the column
#' names used in topTable/topTreat dataframes.  Any other dataframe
#' can be used with fold-change, intensity, and significance measures, with appropriate
#' arguments to define the column names to use provided. By default, the
#' column names will be used for the axis labels, but can be overridden with xlab and ylab arguments.
#'
#' A significance measure (which defaults to P.Value <= 0.01) is used to color code genes that
#' are significantly increased or decreased.  Use the appropriate arguments to use an FDR measure instead.
#'
#' Sensible defaults are chosen for symbols (Size, Shape, Color, and Fill), but they can be
#' adjusted through the use of optional arguments. A length of 3 is
#' required for these arguments which applies the attributes in this order:
#' Increased, NoChange, Decreased.
#'
#' @param contrastDF Dataframe with LogIntensity and LogRatio columns and optionally a p-value or FDR column.
#' @param logRatioCol Name of the LogRatio column (Default = "logFC")
#' @param logIntCol Name of the LogIntensity column (Default = "AveExpr")
#' @param pvalCol Name of the p-value or FDR column (Default = "P.Value")
#' @param xlab X axis label (Defaults to the LogIntensity column name)
#' @param ylab Y axis label (Defaults to the LogRatio column name)
#' @param title Plot title (optional)
#' @param pthreshold Used to color points (Default = 0.01)
#' @param geneSymLabels Labels for gene symbols in contrastDF.
#' @param geneSymCol Name of the gene symbol column in contrastDF.  The gene symbol is
#'    not in topTable output by default so the user has to bind this column
#'    to the dataframe in advance.  Then this column will be used to label
#'    significantly changed points.
#' @param rugColor Specify color for rug density plot along x and y axes. Set to NULL
#'    to disable rug layer. (Default = NULL)
#' @param rugAlpha Sets the transparency for the rug layer.  An alpha < 1.0 can take awhile to draw,
#'   thus the default of 1.0, but for a final plot try 0.1 or 0.2 for alpha for more informative rug.
#' @param symbolSize Size of symbols for Up, no change, and Down. default = c(4, 1, 4);
#'        Note: All three cannot be the same size. Decimal values are acceptable to help offset that
#'        (e.g. 4, 4.1, 4.2).
#' @param symbolShape Shape of the symbols for Up, no change, and Down; Default = c(21, 20, 21)
#'        (20 = filled circle, 21 = fillable open circle) Note: The same symbol shape cannot
#'        be selected for all three symbols.
#'        See \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types}
#' @param symbolColor c(Up, NoChange, Down); default = c("black", "grey25", "grey0")
#'        See \url{http://research.stowers-institute.org/efg/R/Color/Chart}
#'        Note: Colors cannot be duplicated.
#' @param symbolFill Set the fill color for the symbols. Note only symbols 21-25 are fillable.
#'        This will have no effect on other symbols.
#'        Default = c("red3", "grey25", "deepskyblue4")
#'        Note: Colors cannot be duplicated.
#' @param alpha Controls the transparency of the plotted points (0-1; Default = 0.5)
#' @param sizeBySignificance Set to TRUE to size points by the negative Log10 of the
#'        Significance measure (Default = FALSE)
#' @param referenceLine Color for an intercept = 0 horizontal reference line
#'        (Default = "blue"; NULL disables)
#' @param foldChangeLines Position of reference horizontal lines for fold change
#'        (Default = log2(1.5); NULL disables)
#' @param refLineThickness Controls size of the horizontal reference line through geom_hline.
#'   (Default = 1)
#' @param lineFitType Enable a line fit through the data (Default = "loess";
#'        "lm" produces a linear fit. NULL disables)
#' @param lineFitColor Color for the fit line (Default = "goldenrod1")
#' @param legendPosition One of "top", "bottom", "left", "right", "ne", "se", "nw", "sw", NULL.
#'        top/bottom/left/right place the legend outside the figure.  ne/se/nw/sw place the figure
#'        inside the figure. NULL disables the legend. Default = "right"
#' @param baseFontSize The smallest size font in the figure in points. Default = 12.
#' @param themeStyle "bw" or "grey" which correspond to theme_bw or theme_grey respectively.
#'        Default = bw"
#' @param footnote Optional string placed right justified at bottom of plot.
#' @param footnoteSize Applies to footnote. (Default = 3)
#' @param footnoteColor Applies to footnote. (Default = "black")
#' @param footnoteJust Value 0-1. 0 is left justified, 1 is right justified, 0.5 is centered. (Default = 1)
#'
#' @return ggplot object
#'
#' @examples
#' \dontrun{
#'    # Simple plot with custom title(contrastDF is a topTable dataframe)
#'    myPlot <- profilePlot(contrastDF, title = "Plot Title")
#'
#'    # Some options with a custom datafile
#'    myPlot <- profilePlot(contrastDF,
#'                          pthreshold = 0.1,
#'                          logRatioCol = "Log2ratio",
#'                          logIntCol = "AverageIntensity",
#'                          pvalCol = "BHFDR",
#'                          xlab = "Log2Intensity", ylab = "Log2Ratio",
#'                          title = "Profile Plot Title",
#'                          referenceLine = "blue",
#'                          legendPosition = "ne")
#' }
#'
#' @import ggplot2 magrittr
#' @importFrom dplyr left_join filter arrange
#' @importFrom assertthat assert_that
#' @importFrom ggrepel geom_text_repel
#'
#' @export
profilePlot <- function(contrastDF,
                        logRatioCol = "logFC",
                        logIntCol = "AveExpr",
                        pvalCol = "P.Value",
                        pthreshold = 0.01,
                        geneSymLabels,
                        geneSymCol,
                        rugColor = NULL,
                        rugAlpha = 1.0,
                        xlab = NULL, ylab = NULL, title = NULL,
                        symbolSize = c(4, 1, 4),
                        symbolShape = c(21, 20, 21),
                        symbolColor = c("black", "grey25", "grey0"),
                        symbolFill = c("red3", "grey25", "deepskyblue4"),
                        alpha = 0.5,
                        sizeBySignificance = FALSE,
                        referenceLine = "grey25",
                        foldChangeLines = log2(1.5),
                        refLineThickness = 1,
                        lineFitType = "loess",
                        lineFitColor = "goldenrod1",
                        legendPosition = "right",
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
                                msg = "geneSymCol column not found in contrastDF.")
    }
    if (!missing(symbolSize) || !missing(symbolShape) || !missing(symbolColor) || !missing(symbolFill)) {
        assertthat::assert_that(!length(symbolSize) == 3,
                                !length(symbolShape) == 3,
                                !length(symbolColor) == 3,
                                !length(symbolFill) == 3,
                                msg = "All specified symbol arguments must be of length 3, including symbolSize, symbolShape, symbolColor, and symbolFill.")
    }

    groupNames <- c("Increased", "No Change", "Decreased")
    names(symbolShape) = groupNames
    names(symbolSize) = groupNames
    names(symbolColor) = groupNames
    names(symbolFill) = groupNames

    ssc = data.frame(group = groupNames,
                     symbolShape = symbolShape,
                     symbolSize = symbolSize,
                     symbolColor = symbolColor,
                     symbolFill = symbolFill,
                     order = c(1,3,2),
                     stringsAsFactors = FALSE) %>% dplyr::arrange(order)

    # Columns to plot
    # Capture the labels from the colname
    xlabel = logIntCol
    ylabel = logRatioCol
    # Now make the columnames suitable for use with aes_string
    x = make.names(colnames(contrastDF)[colnames(contrastDF) == logIntCol])
    y = make.names(colnames(contrastDF)[colnames(contrastDF) == logRatioCol])
    colnames(contrastDF)[colnames(contrastDF) == logIntCol] = make.names(colnames(contrastDF)[colnames(contrastDF) == logIntCol])
    colnames(contrastDF)[colnames(contrastDF) == logRatioCol] = make.names(colnames(contrastDF)[colnames(contrastDF) == logRatioCol])

    # Need a NLP column for sizing
    if (sizeBySignificance == TRUE) {
        contrastDF$negLog10P = -log10(contrastDF[[pvalCol]])
    }

    DEup <- contrastDF[[pvalCol]] <= pthreshold & contrastDF[[logRatioCol]] > 0
    DEdn <- contrastDF[[pvalCol]] <= pthreshold & contrastDF[[logRatioCol]] < 0
    DEnot <- !DEup & !DEdn
    # Create group factor column in contrastDF
    contrastDF$group <- NA
    contrastDF$group[DEup] <- "Increased"
    contrastDF$group[DEdn] <- "Decreased"
    contrastDF$group[DEnot] <- "No Change"
    contrastDF %<>% dplyr::left_join(ssc)
    contrastDF$group %<>% factor(levels = c("Increased", "Decreased", "No Change"))

    # Set an order field to force plotting of NoChange first
    contrastDF$order <- NA
    contrastDF$order[DEup] <- 1
    contrastDF$order[DEdn] <- 1
    contrastDF$order[DEnot] <- 0

    profilePlot <- ggplot(contrastDF, aes_string(x = x, y = y)) +
        aes(shape = group,
            size = group,
            color = group,
            fill = group,
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
        geom_point(alpha = alpha)

    # Optional Decorations
    if (!is.null(rugColor)) {
        profilePlot <- profilePlot + geom_rug(data = contrastDF,
                                              inherit.aes = FALSE,
                                              color = rugColor,
                                              alpha = rugAlpha,
                                              show.legend = FALSE,
                                              aes_string(x = x, y = y))
    }

    if (sizeBySignificance == TRUE) {
        profilePlot <- profilePlot + aes(size = negLog10P) +
            scale_size_continuous()
    }

    if (!is.null(referenceLine)) {
        profilePlot <- profilePlot +
            geom_hline(yintercept = 0,
                       color = referenceLine,
                       size = refLineThickness,
                       alpha = 0.5)
    }

    if (!is.null(foldChangeLines)) {
        profilePlot <- profilePlot +
            geom_hline(yintercept = foldChangeLines,
                       color = symbolFill["Increased"],
                       size = refLineThickness,
                       alpha = 0.5) +
            geom_hline(yintercept = -foldChangeLines,
                       color = symbolFill["Decreased"],
                       size = refLineThickness,
                       alpha = 0.5)
    }

    if (!is.null(lineFitType)) {
        profilePlot <- profilePlot +
            geom_smooth(aes(group = NULL, shape = NULL, size = NULL, color = NULL, fill = NULL),
                        method = tolower(lineFitType),
                        size = refLineThickness,
                        color = lineFitColor,
                        alpha = alpha,
                        se = FALSE,
                        show.legend = FALSE)
    }

    # Add genesym labels to increased, decreased genes
    if (!missing(geneSymLabels) & !missing(geneSymCol)) {
        idx <- contrastDF[[geneSymCol]] %in% geneSymLabels
        contrastDFsubset <- contrastDF[idx,]
        profilePlot <- profilePlot +
            geom_text_repel(data = contrastDFsubset,
                            aes_string(x = x, y = y, label = geneSymCol),
                            show.legend = FALSE)
    }

    # Add axis Labels
    if (is.null(xlab)) {
        profilePlot <- profilePlot + xlab(xlabel)
    } else {
        profilePlot <- profilePlot + xlab(xlab)
    }
    if (is.null(ylab)) {
        profilePlot <- profilePlot + ylab(ylabel)
    } else {
        profilePlot <- profilePlot + ylab(ylab)
    }
    if (!is.null(title)) {
        profilePlot <- profilePlot +
            ggtitle(title)
    }

    # Set the font size before placing the legend
    if (tolower(themeStyle) == "bw") {
        profilePlot <- profilePlot + theme_bw() + baseTheme(baseFontSize)
    } else {
        profilePlot <- profilePlot + theme_grey() + baseTheme(baseFontSize)
    }

    # Footnote
    if (!missing(footnote)) {
        profilePlot <- addFootnote(profilePlot,
                                   footnoteText = footnote,
                                   footnoteSize = footnoteSize,
                                   footnoteColor = "black",
                                   footnoteJust = footnoteJust)
    }

    profilePlot <- setLegendPosition(profilePlot, legendPosition, themeStyle)

    return(profilePlot)
}
