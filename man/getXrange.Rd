% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ThemePack.R
\name{getXrange}
\alias{getXrange}
\title{Function  getXrange}
\usage{
getXrange(p)
}
\arguments{
\item{p}{A ggplot object}
}
\value{
a vector of length 2 with xmin, xmax
}
\description{
Deprecated.  Use xrange instead.
A simple function to spare me from looking up the syntax everytime
I want to get the x range of a plot
}
\examples{
xrange<- getXrange(Myggplot)

}
\author{
John Thompson, \email{john.thompson@bms.com}
}
\keyword{ggplot2}
