#' Plot edgeR dispersion from a DGEobj
#'
#' Creates an edgeR dispersion plot for RNA-Seq QC purposes.  Takes a counts matrix
#' or DGEList for input.  Dispersion is plotted against AveLogCPM.  Optionally,
#' the plot can instead be Biological Coefficient of Variation (BCV is the square root of
#' dispersion) against AveLogCPM.
#'
#' @param DGEdata Counts matrix or DGEList (Required)
#' @param designMatrix A design matrix created by stats::model.matrix (Required)
#' @param plotType One of "dispersion" or "BCV" (Default = "dispersion")
#' @param symbolSize (Default = 1)
#' @param symbolShape see
#'   http://www.cookbook-r.com/Graphs/Shapes_and_line_types/ (Default = 1)
#' @param symbolColor Default = "darkblue"
#' @param symbolFill  Default = "darkblue"
#' @param symbolAlpha Transparency for the points. Value from 0 to 1. Smaller
#'   indicate more transparency (Default = 0.3)
#' @param linefitSize (Default = 1)
#' @param linefitColor (Default = "yellow")
#' @param lineFit (Default = NULL) Any type supported by geom_smooth (Loess is
#'   recommended).
#' @param lineType Any ggplot2 lineType supported by geom_smooth (Default = 1, solid)
#' @param lineAlpha Alpha desired for geom_smooth line (Default = 1)
#' @param rugColor (Default = NULL)  Set to Null to disable rug plots.
#' @param rugAlpha Transparency for the rug plots (Default = 0.02)
#' @param ... Extra parameters to pass to edgeR::estimateDisp
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#'    myGgplot <- plotDispersion(myDGElist)
#'    myGgplot <- plotDispersion(myDGEobj)
#' }
#'
#' @importFrom assertthat assert_that
#' @importFrom edgeR calcNormFactors estimateDisp DGEList
#'
#' @export
plotDispersion <- function(DGEdata,
                           designMatrix,
                           plotType = "dispersion",
                           symbolSize = 1,
                           symbolShape = 1,
                           symbolColor = "darkblue",
                           symbolFill = "darkblue",
                           symbolAlpha = 0.3,
                           linefitSize = 1,
                           linefitColor = "yellow",
                           lineFit = NULL,
                           lineType = 1,
                           lineAlpha = 1,
                           rugColor = NULL,
                           rugAlpha = 0.02,
                           ...) {

    assertthat::assert_that(!missing(DGEdata),
                            msg = "DGEdata must be specified.")
    if (!missing(designMatrix)) {
        assertthat::assert_that("matrix" %in% class(designMatrix),
                                msg = "designMatrix must be specified and should be of class 'matrix'.")
    }

    if (class(DGEdata)[[1]] == "DGEList") {
        dgelist <- DGEdata %>%
            edgeR::calcNormFactors() %>%
            edgeR::estimateDisp(design = designMatrix, robust = TRUE, ...)
    } else {
        dgelist <- DGEdata %>%  # Process a counts matrix
            as.matrix %>%
            edgeR::DGEList() %>%
            edgeR::calcNormFactors() %>%
            edgeR::estimateDisp(design = designMatrix, robust = TRUE, ...)
    }

    if (tolower(plotType == "dispersion")) {
        plotdata <- data.frame(AveLogCPM = dgelist$AveLogCPM, Dispersion = dgelist$tagwise.dispersion)
    } else {
        plotdata <- data.frame(AveLogCPM = dgelist$AveLogCPM, Dispersion = sqrt(dgelist$tagwise.dispersion))
    }

    MyDispPlot <- ggplot(plotdata, aes(x = AveLogCPM, y = Dispersion)) +
        geom_point(size = symbolSize,
                   shape = symbolShape,
                   fill = symbolFill,
                   color = symbolColor,
                   alpha = symbolAlpha)

    if (!is.null(lineFit)) {
        MyDispPlot <- MyDispPlot +
            geom_smooth(method = lineFit,
                        size = linefitSize,
                        color = linefitColor,
                        linetype = lineType,
                        alpha = lineAlpha)
    }

    if (!is.null(rugColor)) {
        MyDispPlot <- MyDispPlot + geom_rug(data = plotdata,
                                            inherit.aes = FALSE,
                                            color = rugColor,
                                            alpha = rugAlpha,
                                            show.legend = FALSE,
                                            aes(x = AveLogCPM, y = Dispersion))
    }

    MyDispPlot <- MyDispPlot +
        ggtitle("EdgeR Dispersion Plot") +
        expand_limits(x = 0, y = 0) +
        theme_grey()

    if (tolower(plotType) == "bcv") {
        MyDispPlot <- MyDispPlot +
            ylab("BCV") +
            ggtitle("EdgeR BCV Plot")
    }

    return(MyDispPlot)
}
