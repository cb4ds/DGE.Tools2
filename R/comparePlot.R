#' Create formatted scatterplot of first two cols of data.frame
#'
#' Creates a nicely formatted scatterplot of the
#' first two columns in a dataframe.  Several formatting options are
#' available to enhance the plots.  If p-values or FDR values and a threshold are supplied,
#' the plot is color coded to show X unique, Y unique, and
#' common differentially expressed (DE) genes in different colors.
#'
#' Other options add an identity line (slope = 1, intercept = 0) and/or a crosshair at (0, 0).
#'
#' comparePlot has also been implemented with a scalable theme that makes it easy to scale font sizes larger
#' or smaller for PPT or knitr output.  Add either theme_grey(n) or theme_bw(n),
#' where n equals the base font size to adjust the text size (12 works well for knitr, 18 or 24 works well
#' for PPT), e.g. MyPlot + theme_grey(18).  However, theme_grey() and theme_bw() have the
#' side effect of resetting the legend position, and should not be used if the custom legend locations
#' need to be preserved. Instead, utilize the baseFontSize
#' argument in this function to set font sizes and preserve a custom legend location.
#'
#' \strong{Data Structure for the input dataframe:}
#'
#' The x and y values should be in the first two columns. By default, their
#' column names will be used for the axis labels. The x and y labels can be changed
#' using the xlab and ylab arguments.
#'
#' Optionally, significance measures in the form of p-values or FDR values can be supplied
#' for X and Y respectively. If provided, these columns \strong{must} be named "xp" and "yp" respectively.
#' Together with a threshold (which defaults to 0.01), these values will
#' be used to color code the plot to show X Unique, Y Unique, and Common DE
#' genes.  Use either p-values or FDR values for the significance columns. Use the
#' pThreshold argument to set a proper threshold for FDR values.
#'
#' Sensible defaults are chosen for symbols (Size, Shape, Color and Fill), but optional
#' arguments allow much about the symbols to be customized. A length of 4 is
#' required for these arguments which applies the attributes in this order:
#' Common, X Unique, Y Unique, and Not Significant.
#'
#' Note: if p-values or FDR values are not used to color the plot, the X Unique color
#' values are used.
#'
#' @param compareDF A dataframe with the first two columns representing the x and y variables.
#'          Optionally add xp and yp columns to hold p-values or FDR values.
#' @param xlab X-axis label (Default to first column name)
#' @param ylab Y-axis label (Default to second column name)
#' @param title Plot title (Optional)
#' @param pThreshold Used to color points (Default = 0.01)
#' @param symbolSize Size of symbols; default = c(4, 4, 4, 2)
#' @param symbolShape Shape of the symbols; Default = c(21, 21, 21, 20)  (20 = filled circle,
#'        21 = fillable open circle) See \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types}
#' @param symbolColor c(Common, xUnique, yUnique, NoChange);
#'        default = c("black", "grey0", "grey1", "grey25")
#' @param symbolFill Set the fill color for the symbols. Note only symbols 21-25 are fillable. This will
#'        have no effect on other symbols. Default = c("darkgoldenrod1", "deepskyblue4", "red3", "grey25")
#' @param alpha Controls the transparency of the plotted points (0-1; Default = 0.5)
#' @param rugColor Specify color for rug density plot along X and Y axes. Set to NULL
#'        to disable rug layer. (Default = NULL)
#' @param rugAlpha Sets the transparency for the rug layer. alpha <1.0 takes
#'        a long time to draw, so for a final plot try 0.1 or 0.2 for alpha for more informative rug.
#' @param crosshair Color for the crosshair; (Default = "grey50", NULL disables)
#'        See \url{http://research.stowers-institute.org/efg/R/Color/Chart}
#' @param referenceLine Color for a slope=1, intercept=0 reference line
#'        (Default = "darkgoldenrod1"; NULL disables)
#' @param refLineThickness Set thickness for crosshair and referenceLine (Default = 1)
#' @param dens2D Overlay a density2D layer (Default = TRUE) Note: only applies
#'        to simple compare plot without p-value coloring.
#' @param legendPosition One of "top", "bottom", "left", "right", "ne", "se", "nw", "sw", NULL.
#'        top/bottom/left/right place the legend outside the figure.  ne/se/nw/sw place the figure
#'        inside the figure. NULL disables the legend. Default = "right"
#' @param baseFontSize The smallest size font in the figure in points. Default = 12
#' @param themeStyle "bw" or "grey" which corresponds to theme_bw or theme_grey respectively.
#'        Default = "bw"
#' @param footnote Optional string placed right justified at bottom of plot.
#' @param footnoteSize Applies to footnote. (Default = 3)
#' @param footnoteColor Applies to footnote. (Default = "black")
#' @param footnoteJust Value 0-1. 0 is left justified, 1 is right justified, 0.5 is centered. (Default = 1)
#'
#' @return ggplot object
#'
#' @examples
#' \dontrun{
#'   # Retrieve the first two contrasts from a DGEobj as a list of dataframes (length = 2; named items)
#'   contrastList <- getType(DGEobj, "topTable")[1:2]
#'
#'   # Capture the default logFC and P.Value
#'   compareDat <- comparePrep(contrastList)
#'
#'   # Switch to an FDR value for the significance measure
#'   compareDat <- comparePrep(contrastList, significanceCol = "adj.P.Val")
#'
#'   # Draw the plot
#'   cPlot <- comparePlot(compareDat, title = "Plot Title")
#'   print(cPlot)
#'
#'   # Deluxe Plot with bells and whistles.
#'   myPlot <- comparePlot(compareDF,
#'                         pThreshold = 0.5,
#'                         xlab = "x Axis Label",
#'                         ylab = "y Axis Label"
#'                         title = "Plot Title",
#'                         crosshair = "red",
#'                         referenceLine = "blue",
#'                         legendPosition = "nw")
#' }
#'
#' @import ggplot2 magrittr
#' @importFrom dplyr left_join
#' @importFrom assertthat assert_that
#'
#' @export
comparePlot <- function(compareDF,
                        pThreshold = 0.01,
                        xlab = NULL, ylab = NULL,
                        title = NULL,
                        symbolSize = c(4, 4, 4, 2),
                        symbolShape = c(21, 21, 21, 20),
                        symbolColor = c("black", "grey0", "grey1", "grey25"),
                        symbolFill = c("darkgoldenrod1", "deepskyblue4", "red3", "grey25"),
                        alpha = 0.5,
                        rugColor = NULL,
                        rugAlpha = 1.0,
                        crosshair = "grey50",
                        referenceLine = "darkgoldenrod1",
                        refLineThickness = 1,
                        dens2D = TRUE,
                        legendPosition = "right",
                        baseFontSize = 12,
                        themeStyle = "bw",
                        footnote,
                        footnoteSize = 3,
                        footnoteColor = "black",
                        footnoteJust = 1
) {

    set.seed(1954)

    assertthat::assert_that(sum(apply(compareDF, 2, FUN = is.numeric)) >= 2,
                            msg = "Need at least two numeric columns in compareDF.")
    if (!missing(symbolSize) || !missing(symbolShape) || !missing(symbolColor) || !missing(symbolFill)) {
        assertthat::assert_that(!length(symbolSize) == 4,
                                !length(symbolShape) == 4,
                                !length(symbolColor) == 4,
                                !length(symbolFill) == 4,
                                msg = "All specified symbol arguments must be of length 4, including symbolSize, symbolShape, symbolColor, and symbolFill.")
    }

    names(symbolShape) = c("Common", "xUnique", "yUnique", "NoChange")
    names(symbolSize)  = c("Common", "xUnique", "yUnique", "NoChange")
    names(symbolColor) = c("Common", "xUnique", "yUnique", "NoChange")
    names(symbolFill)  = c("Common", "xUnique", "yUnique", "NoChange")

    ssc = data.frame(group = c("Common", "X Unique", "Y Unique", "Not Significant"),
                     symbolShape = symbolShape,
                     symbolSize = symbolSize,
                     symbolColor = symbolColor,
                     symbolFill = symbolFill,
                     stringsAsFactors = FALSE)

    # Used to set uniform square scale
    scalemax = compareDF[,1:2] %>% as.matrix %>% abs %>% max %>% multiply_by(1.05)

    # Capture the labels from the colname
    xlabel = colnames(compareDF)[1]
    ylabel = colnames(compareDF)[2]
    # Now make the column names suitable for use with aes_string
    x = make.names(colnames(compareDF)[1])
    y = make.names(colnames(compareDF)[2])
    colnames(compareDF)[1:2] = make.names(colnames(compareDF)[1:2])

    # SIMPLE PLOT: plot all data (no significance measures supplied)
    if (is.null(compareDF[["xp"]]) | is.null(compareDF[["yp"]])) {
        CompPlot = ggplot(compareDF, aes_string(x = x, y = y)) +
            geom_point(shape = 1,
                       size = symbolSize[["xUnique"]],
                       color = symbolFill[["xUnique"]],
                       fill = symbolFill[["xUnique"]],
                       alpha = alpha) +
            coord_equal(xlim = c(-scalemax, scalemax), ylim = c(-scalemax, scalemax))
        if (dens2D == TRUE) {
            CompPlot = CompPlot + geom_density2d(color = symbolFill[["yUnique"]])
        }
    } else
        # DELUXE PLOT: plot groups in different colors/shapes
        if (!is.null(compareDF[["xp"]]) & !is.null(compareDF[["yp"]])) {
            # Plot the subsets
            xindx = compareDF[["xp"]] <= pThreshold
            yindx = compareDF[["yp"]] <= pThreshold

            # Boolean indexes to parse groups
            bothindx = xindx & yindx
            neitherindx = !xindx & !yindx
            xindx = xindx & !bothindx # Unique to X
            yindx = yindx & !bothindx # Unique to y

            # Create group factor column in compareDF
            compareDF$group = NA
            compareDF$group[bothindx] = "Common"
            compareDF$group[xindx] = "X Unique"
            compareDF$group[yindx] = "Y Unique"
            compareDF$group[neitherindx] = "Not Significant"
            compareDF %<>% dplyr::left_join(ssc)
            compareDF$group %<>% factor(levels = c("Common", "X Unique", "Y Unique", "Not Significant"))

            # Set an order field to control order of plotting
            # Plot order is 1) Not Significant plotted first, 2) randomly plot X and Y Unique
            # then plot Common last
            compareDF$order = sample.int(nrow(compareDF)) + # Add random order (seeded for reproducibility)
                nrow(compareDF) * (compareDF$group == "Common") +  # Assign common a high value to sort last
                -nrow(compareDF) * (compareDF$group == "Not Significant") # Assign NotSig group neg values to sort first

            CompPlot <- ggplot(compareDF, aes_string(x = x, y = y)) +
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
                # Make it square with same axis scales
                coord_equal(xlim = c(-scalemax, scalemax), ylim = c(-scalemax, scalemax)) +
                # Box around the legend
                theme(legend.background = element_rect(fill = "gray95", size = .5, linetype = "dotted"))
        }

    # Optional Decorations
    if (!is.null(rugColor)) {
        CompPlot <- CompPlot + geom_rug(data = compareDF,
                                        inherit.aes = FALSE,
                                        color = rugColor,
                                        alpha = rugAlpha,
                                        show.legend = FALSE,
                                        aes_string(x = x, y = y))
    }

    if (!is.null(crosshair)) {
        CompPlot <- CompPlot +
            geom_hline(yintercept = 0,
                       color = crosshair,
                       size = refLineThickness,
                       alpha = 0.5) +
            geom_vline(xintercept = 0,
                       color = crosshair,
                       size = refLineThickness,
                       alpha = 0.5)
    }

    if (!is.null(referenceLine)) {
        CompPlot <- CompPlot +
            geom_abline(slope = 1,
                        intercept = 0,
                        color = referenceLine,
                        size = refLineThickness,
                        alpha = 0.5)
    }

    # Add Labels
    if (is.null(xlab)) { # Use colname unless supplied as argument
        CompPlot <- CompPlot + xlab(xlabel)
    } else {
        CompPlot <- CompPlot + xlab(xlab)
    }
    if (is.null(ylab)) {
        CompPlot <- CompPlot + ylab(ylabel)
    } else {
        CompPlot <- CompPlot + ylab(ylab)
    }
    if (!is.null(title)) {
        CompPlot <- CompPlot + ggtitle(title)
    }

    # Set the font size before placing the legend
    if (tolower(themeStyle) == "bw") {
        CompPlot <- CompPlot + theme_bw() + baseTheme(baseFontSize)
    } else {
        CompPlot <- CompPlot + theme_grey() + baseTheme(baseFontSize)
    }

    CompPlot <- setLegendPosition(CompPlot, legendPosition, themeStyle)

    if (!missing(footnote)) {
        CompPlot <- addFootnote(CompPlot,
                                footnoteText = footnote,
                                footnoteSize = footnoteSize,
                                footnoteColor = "black",
                                footnoteJust = footnoteJust,
                                yoffset = 0.05)
    }

    return(CompPlot)
}


#' Create a data.frame to use with comparePlot from topTable data.frames
#'
#' Takes two topTable dataframes and outputs a dataframe suitable for function
#' comparePlot() (2 columns of LogRatio data and 2 columns of significant
#' measures). Filter the two topTable to contain only the intersecting genes
#' (present in both datasets). The two dataframes must have the same type of
#' gene IDs as rownames.
#'
#' @param contrastList A named list of 2 topTable dataframes (Required). The
#'   names are used as column names for the value columns in the output.
#' @param valueCol Name of column containing values to plot (Default = "logFC")
#' @param significanceCol Name of column to use for significance (Default = "P.Value")
#'
#' @return A data frame with 2 LogRatio measurements and 2 significance columns.  Columns 1 and 3 map
#' to sample 1 and columns 2 and 4 map to sample 2.  The returned dataframe is formatted as expected
#' by the comparePlot function.
#'
#' @examples
#' \dontrun{
#'   # Retrieve the 1st two contrasts from a DGEobj
#'   contrastList <- getType(dgeObj, "topTable")[1:2]
#'
#'   # Capture the default logFC and P.Value
#'   compareDat <- comparePrep(contrastList)
#'
#'   # Switch to an FDR value for the significance measure
#'   compareDat <- comparePrep(contrastList, significanceCol="adj.P.Val")
#'
#'   # Draw the plot
#'   cPlot <- comparePlot(compareDat)
#' }
#'
#' @importFrom assertthat assert_that
#' @import magrittr
#'
#' @export
comparePrep <- function(contrastList,
                        valueCol = "logFC",
                        significanceCol = "P.Value"){

    assertthat::assert_that(length(contrastList) == 2,
                            !is.null(names(contrastList)),
                            "data.frame" %in% class(contrastList[[1]]),
                            "data.frame" %in% class(contrastList[[2]]),
                            msg = "contrastList must be a named list of length 2 where both items are of class 'data.frame'.")
    assertthat::assert_that(valueCol %in% colnames(contrastList[[1]]),
                            valueCol %in% colnames(contrastList[[2]]),
                            msg = "The valueCol must be included in the colnames of both items of contrastList.")
    assertthat::assert_that(significanceCol %in% colnames(contrastList[[1]]),
                            significanceCol %in% colnames(contrastList[[2]]),
                            msg = "The significanceCol must be included in the colnames of both items of contrastList.")

    ttNames <- names(contrastList)
    tt1 <- contrastList[[1]]
    tt2 <- contrastList[[2]]

    commonIDs <- intersect(rownames(tt1), rownames(tt2))
    assertthat::assert_that(!length(commonIDs) < 1,
                            msg = "No common gene IDs were found between the two dataframes in contrastList.")

    # Filter both tables to the same set of genes in the same order
    tt1 %<>% rownames_to_column(var = "geneid") %>%
        dplyr::filter(geneid %in% commonIDs) %>%
        dplyr::arrange(geneid)
    tt2 %<>% rownames_to_column(var = "geneid") %>%
        dplyr::filter(geneid %in% commonIDs) %>%
        dplyr::arrange(geneid)

    assertthat::assert_that(all(tt1$geneid == tt2$geneid),
                            msg = "Gene IDs in the two topTable files in contrastList are not identical.")

    # Assemble the return df
    df <- dplyr::bind_cols(fc1 = tt1[[valueCol]],
                           fc2 = tt2[[valueCol]],
                           xp = tt1[[significanceCol]],
                           yp = tt2[[significanceCol]]) %>%
        as.matrix() %>%
        as.data.frame() %>%
        set_colnames(c(ttNames[1], ttNames[2], "xp", "yp")) %>%
        set_rownames(tt1$geneid)

    return(df)
}
