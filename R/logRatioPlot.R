#' Plot logRatio contrasts
#'
#' Intended to plot a set of contrast results, one plot for each gene of
#' interest. Input is a tidy datafile constructed from topTable output and
#' requires logFC, CI.L, and CI.R columns as well as a gene identifier of choice.
#' Outputs a ggplot2 object faceted by the facetColname or a list of individual
#' ggplots, one for each facetColname value (typically gene).
#'
#' @param contrastsDF A tidy dataframe of data to plot (Required) (see ?tidyContrasts).
#' @param facetColname Define the column name to separate plots (Required) (e.g. GeneID).
#' @param xColname Define the column name to group boxplots by (Required) (e.g. Contrast).
#' @param yColname Define the column name for the output of the boxplots (Default = "logFC")
#' @param CI.R_colname Define name of the CI high value (Default = "CI.R")
#' @param CI.L_colname Define name of the CI low value (Default =  "CI.L")
#' @param xOrder Define the order for the groups in each plot.  Should
#'   contain values in unique(contrastsDF$group) listed in the order that
#'   groups should appear in the plot. (Optional; Default = unique(contrastsDF[xColname]))
#' @param plotType One of "bar" or "point" (Default = "bar")
#' @param refLine Adds a horizontal line at y = 0 (Default = TRUE)
#' @param refLineColor Color for the reference line (Default = "red")
#' @param xlab X axis label (Defaults to xColname)
#' @param ylab Y axis label (Defaults to yColname)
#' @param title Plot title (Optional)
#' @param barColor Color for the bar outline (Default = "dodgerblue4")
#' @param barFill Color for the bar area (Default = "dodgerblue3")
#' @param barSize Set the bar size (thickness of each bar perimeter; Default = 0.1)
#' @param barWidth Set the bar width (Default = 0.8)
#' @param barAlpha Transparency for the bar layer (Default = 1)
#' @param pointColor Color for the point layer (Default = "grey30")
#' @param pointFill Fill color for the point layer (Default = "dodgerblue4")
#' @param pointShape Shape for the point layer (Default = 21; fillable circle)
#' @param pointAlpha Transparency for the box layer (Default = 1)
#' @param pointSize Size of the points (Default = 4)
#' @param lineLayer Add a fitted line layer (Default = FALSE)
#' @param lineColor Color of the line fit (Default = "dodgerblue4")
#' @param lineSize Size of the line fit (Default = 1)
#' @param lineFit Type of fit to use.  One of c("auto", "lm", "glm", "gam",
#'   "loess"). (Default = "loess")
#' @param lineType One of c("solid", "dashed", "dotted", "dotdash", "longdash",
#'   "twodash"). (Default = "solid")
#' @param baseFontSize The smallest size font in the figure in points. (Default =
#'   12)
#' @param themeStyle "bw" or "grey" which correspond to theme_bw or theme_grey
#'   respectively. (Default = "bw")
#' @param facet Specifies whether to facet (TRUE) or print individual plots
#'   (FALSE)  (Default = TRUE)
#' @param facetCol Explicitly set the number of rows for the facet plot. Default
#'   behavior will automatically set the columns. (Default = ceiling(sqrt(length(unique(contrastsDF[facetCol])))))
#' @param xAngle Angle to set the sample labels on the X axis (Default =  45; Range = 0-90)
#' @param scales Specify same scales or independent scales for each subplot (Default = "free_y";
#'   Allowed values: "fixed", "free_x", "free_y", "free")
#'
#' @return ggplot object. If facet = TRUE (default), returns a faceted ggplot object. If
#'   facet = FALSE, returns a list of ggplot objects indexed
#'   by observation (gene) names.
#'
#' @examples
#' \dontrun{
#'   # DGEobj example
#'   # Put contrasts in tidy format keeping logFC, and confidence limits contrastsDF
#'   tidyDat <- tidyContrasts(DGEobj,
#'                            rownameColumn = "EnsgID",
#'                            includeColumns = c("logFC", "CI.R", "CI.L"))
#'
#'   # Add gene symbols from geneData
#'   ens2genesym <- DGEobj$geneData %>%
#'                  rownames_to_column(var = "EnsgID") %>%
#'                  select(EnsgID, GeneSymbol = GeneName)
#'   tidyDat <- left_join(tidyDat, ens2genesym)
#'
#'   # Filter for a small set of genes of interest
#'   idx <- stringr::str_detect(tidyDat$GeneSymbol, "^PPAR")
#'   tidyDat <- tidyDat[idx,]
#'
#'   # Simple barplot
#'   logRatioPlot(tidyDat,
#'                facetColname = "GeneSymbol",
#'                xColname = "Contrast",
#'                facetCol = 2)
#'
#'   # Lineplot with some options
#'   logRatioPlot(tidyDat,
#'                plotType = "point",
#'                facetColname = "GeneSymbol",
#'                xColname = "Contrast",
#'                facetCol = 4,
#'                scales = "fixed",
#'                facet = TRUE,
#'                title = "Test",
#'                pointSize = 4,
#'                lineLayer = TRUE,
#'                lineSize = 0.1,
#'                xAngle = 60)
#' }
#'
#' @import ggplot2 magrittr
#' @importFrom assertthat assert_that
#' @importFrom stringr str_c
#'
#' @export
logRatioPlot <- function(contrastsDF,
                         facetColname,
                         xColname,
                         yColname = "logFC",
                         CI.R_colname = "CI.R",
                         CI.L_colname = "CI.L",
                         xOrder = unique(as.character(contrastsDF[xColname, , drop = TRUE])),
                         plotType = "bar",
                         refLine = TRUE,
                         refLineColor = "red",
                         xlab = xColname, ylab = yColname,
                         title,
                         barColor = "dodgerblue4",
                         barFill = "dodgerblue3",
                         barSize = 0.1,
                         barAlpha = 1,
                         barWidth = 0.9,
                         pointColor = "grey30",
                         pointFill = "dodgerblue4",
                         pointShape = 21,
                         pointAlpha = 1,
                         pointSize = 2,
                         lineLayer = FALSE,
                         lineColor = "dodgerblue4",
                         lineSize = 1,
                         lineType = "solid",
                         lineFit = "loess",
                         baseFontSize = 12,
                         themeStyle = "grey",
                         facet = TRUE,
                         facetCol,
                         xAngle = 45,
                         scales = "free_y") {

    assertthat::assert_that(plotType %in% c("bar", "point"),
                            msg = "plotType must be either 'bar' or 'point'.")

    .addGeoms <- function(myPlot){
        if (tolower(plotType) == "bar") {
            myPlot <- myPlot + geom_bar(stat = "identity",
                                        alpha = barAlpha,
                                        color = barColor,
                                        fill = barFill,
                                        size = barSize,
                                        width = barWidth)
        } else if (tolower(plotType) == "point") {
            myPlot <- myPlot + geom_point(alpha = pointAlpha,
                                          color = pointColor,
                                          fill = pointFill,
                                          size = pointSize,
                                          shape = pointShape)
        }

        # Add error bars if columns present
        if (all(c(CI.L_colname, CI.R_colname) %in% colnames(contrastsDF) )) {
            myPlot <- myPlot + geom_errorbar(aes_string(ymin = CI.L_colname, ymax = CI.R_colname), width = .2)
        } else {
            warning("Confidence limits columns not found.")
        }

        if (lineLayer == TRUE) {
            myPlot <- myPlot + geom_smooth(aes_string(group = facetColname),
                                           method = lineFit,
                                           formula = y ~ x,
                                           color = lineColor,
                                           size = lineSize,
                                           se = FALSE)
        }

        return(myPlot)
    }

    assertthat::assert_that(!missing(contrastsDF),
                            "data.frame" %in% class(contrastsDF),
                            msg = "data must be specified and should be of class 'data.frame'.")
    assertthat::assert_that(facetColname %in% colnames(contrastsDF),
                            msg = "facetColname must be included in the colnames of data.")
    assertthat::assert_that(xColname %in% colnames(contrastsDF),
                            msg = "xColname must be included in the colnames of data.")
    assertthat::assert_that(yColname %in% colnames(contrastsDF),
                            msg = "yColname must be included in the colnames of data.")
    assertthat::assert_that(all(xOrder %in% as.character(contrastsDF[xColname, , drop = TRUE])))

    # Plot code here
    if (facet) {
        # Set facet columns to sqrt of unique observations (rounded up)
        if (missing(facetCol)) {
            numcol <- contrastsDF[facetCol] %>% unique %>% length %>% sqrt %>% ceiling
        } else {
            numcol = facetCol
        }

        myPlot <- ggplot2::ggplot(contrastsDF, aes_string(x = xColname, y = yColname))
        myPlot <- .addGeoms(myPlot)
        facetFormula <- stringr::str_c("~", facetColname, sep = " ")
        myPlot <- myPlot + ggplot2::facet_wrap(facetFormula, ncol = numcol, scales = scales)

        myPlot <- myPlot + ggplot2::xlab(xlab)
        myPlot <- myPlot + ggplot2::ylab(ylab)
        if (!missing(title)) {
            myPlot <- myPlot + ggplot2::ggtitle(title)
        }
        if (tolower(themeStyle) == "bw") {
            myPlot <- myPlot + theme_bw() + baseTheme(baseFontSize)
        } else {
            myPlot <- myPlot + theme_grey() + baseTheme(baseFontSize)
        }

        #rotate xaxis group labels
        if (xAngle > 0) {
            myPlot <- myPlot + theme(axis.text.x = element_text(angle = xAngle, hjust = 1))
        }

        #Add refLine at 0
        if (refLine == TRUE) {
            myPlot <- myPlot + geom_hline(yintercept = 0, color = refLineColor, size = 0.1)
        }

    } else { # Individual plots for each Gene returned in a list

        plotlist <- list()

        for (obs in unique(contrastsDF[[facetColname]])) { # For each gene

            dat <- contrastsDF[contrastsDF[[facetColname]] == obs, ] # Pull data for one gene

            aplot <- ggplot(dat, aes_string(x = xColname, y = yColname)) + # Samples vs Log2CPM
                xlab(xlab) +
                ylab(ylab) +
                ggtitle(obs) +
                theme_grey() + facetTheme(baseFontSize)
            aplot <- .addGeoms(aplot)

            if (!missing(title)) {
                aplot <- aplot + ggplot2::ggtitle(stringr::str_c(title, ": ", obs))
            }

            if (xAngle > 0) {
                aplot <- aplot + theme(axis.text.x = element_text(angle = xAngle, hjust = 1))
            }

            if (refLine == TRUE) {
                aplot <- aplot + geom_hline(yintercept = 0, color = refLineColor, size = 0.1)
            }

            plotlist[[obs]] <- aplot
        }

        myPlot = plotlist
    }

    return(myPlot)
}
