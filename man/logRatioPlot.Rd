% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logRatioPlot.R
\name{logRatioPlot}
\alias{logRatioPlot}
\title{Function logRatioPlot}
\usage{
logRatioPlot(
  data,
  facetColname,
  xColname,
  yColname = "logFC",
  CI.R_colname = "CI.R",
  CI.L_colname = "CI.L",
  xOrder = unique(as.character(data[xColname, , drop = TRUE])),
  plotType = "bar",
  refLine = TRUE,
  refLineColor = "red",
  xlab = xColname,
  ylab = yColname,
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
  scales = "free_y",
  debug = FALSE
)
}
\arguments{
\item{data}{A tidy dataframe of data to plot (required) (see ?tidyContrasts).}

\item{facetColname}{Define the column name to separate plots (required) (e.g. GeneID).}

\item{xColname}{Define the column name to group boxplots by (required) (e.g. Contrast).}

\item{CI.R_colname}{Define name of the CI high value (Default = "CI.R")}

\item{CI.L_colname}{Define name of the CI low value (Default =  "CI.L")}

\item{xOrder}{Define the order for the groups in each plot.  Should
contain values in unique(data$group) listed in the order that you want the
groups to appear in the plot. (optional; default = unique(data[xColname]))}

\item{plotType}{one of "bar" or "point" (Default = bar")}

\item{refLine}{Adds a horizontal line at y=0 (Default=TRUE)}

\item{refLineColor}{Color for the reference line (Default="red")}

\item{xlab}{X axis label (defaults to xColname)}

\item{ylab}{Y axis label (defaults to yColname)}

\item{title}{Plot title (optional)}

\item{barColor}{Color for the bar outline (default = "dodgerblue4")}

\item{barFill}{Color for the bar area (default = "dodgerblue3")}

\item{barSize}{set the bar size (thickness of each bar perimeter; default = 0.1)}

\item{barAlpha}{Transparency for the bar layer (Default = 1)}

\item{barWidth}{set the bar width (Default = 0.8)}

\item{pointColor}{Color for the point layer (Default = "grey30")}

\item{pointFill}{Fill color for the point layer (Default = "dodgerblue4")}

\item{pointShape}{Shape for the point layer (Default = 21; fillable circle)}

\item{pointAlpha}{Transparency for the box layer (Default = 1)}

\item{pointSize}{Size of the points (Default = 4)}

\item{lineLayer}{Add a fitted line layer (Default = FALSE)}

\item{lineColor}{Color of the line fit (Default = "dodgerblue4")}

\item{lineSize}{Size of the line fit (Default = 1)}

\item{lineType}{One of c("solid", "dashed", "dotted", "dotdash", "longdash",
"twodash"). (Default = "solid")}

\item{lineFit}{Type of fit to use.  One of c("auto", "lm", "glm", "gam",
"loess"). (Default = "loess")}

\item{baseFontSize}{The smallest size font in the figure in points. (Default =
12)}

\item{themeStyle}{"bw" or "grey" which correspond to theme_bw or theme_grey
respectively. (Default = "bw")}

\item{facet}{Specifies whether to facet (TRUE) or print individual plots
(FALSE)  (Default = TRUE)}

\item{facetCol}{Explicitly set the number of Rows for the facet plot. Default
behavior will automatically set the columns. (Default = ceiling(sqrt(length(unique(data[facetCol])))))}

\item{xAngle}{Angle to set the sample labels on the Xaxis (Default =  45; Range = 0-90)}

\item{scales}{Specify same scales or independent scales for each subplot (Default = "free_y";
Allowed values: "fixed", "free_x", "free_y", "free")}

\item{debug}{Turn on debug mode (default = FALSE)}

\item{yColame}{Define the column of values for plotting (Default = "logFC").}

\item{addLine}{Adds a line if point layer chosen (Default = FALSE)}
}
\value{
ggplot If Facet=TRUE (default) returns a facetted ggplot object. If
  facet=FALSE, returns a list of ggplot objects indexed
  by observation (gene) names.
}
\description{
Intended to plot a set of contrast results, one plot for each gene of
interest. Input is a tidy datafile constructed from topTable output and
requires logFC, CI.L and CI.R columns as well as a gene identifier of choice.
Outputs a ggplot object facetted by the facetColname or a list of individual
ggplots, one for each facetColname value (typically gene).
}
\examples{

  #DGEobj example
  #Put contrasts in tidy format keeping logFC, and confidence limits data
  tidyDat <-tidyContrasts(dgeObj, rownameColumn="EnsgID", includeColumns = c("logFC", "CI.R", "CI.L"))

  #add gene symbols from geneData
  ens2genesym <- dgeObj$geneData \%>\%
                 rownames_to_column(var="EnsgID") \%>\%
                 select(EnsgID, GeneSymbol=GeneName)
  tidyDat <-  left_join(tidyDat, ens2genesym)

  #filter for a small set of genes of interest
  idx <- stringr::str_detect(tidyDat$GeneSymbol, "^PPAR")
  tidyDat <- tidyDat[idx,]

  #simple barplot
  logRatioPlot(tidyDat,
              facetColname = "GeneSymbol",
              xColname = "Contrast",
              facetCol = 2)

  #lineplot with some options
  logRatioPlot(tidyDat, plotType="point",
                     facetColname = "GeneSymbol",
                     xColname = "Contrast",
                     facetCol=4,
                     scales="fixed",
                     facet=TRUE,
                     title = "Test",
                     pointSize=4,
                     lineLayer=TRUE,
                     lineSize=0.1,
                     xAngle=60)

}
\author{
John Thompson, \email{john.thompson@bms.com}
}
\keyword{barplot}
\keyword{confidence}
\keyword{contrasts}
\keyword{gplot2}
\keyword{intervals}
\keyword{lineplot}
\keyword{logratio}
