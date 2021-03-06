% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractCol.R
\name{extractCol}
\alias{extractCol}
\title{Function  extractCol}
\usage{
extractCol(dflist, colName, robust = "TRUE")
}
\arguments{
\item{dflist}{A list of data.frames which all have the same colnames and same row counts.
The dataframes in the list should have rownames (geneIDs).}

\item{colName}{The name of the data column to extract to a matrix}

\item{robust}{Default = TRUE;  TRUE forces use of a joins to merge columns
which is more reliable, allows you to combine contrasts from different
projects but may nor return items in the same row order as the source
table.  Setting to false invokes a cbind approach that requires all
dataframes to have the same row count and row order but preserves the
original row order}
}
\value{
A dataframe containing the extracted columns
}
\description{
Take a named list of dataframes where each dataframe has the same
column names (e.g. a list of topTable dataframes). Extract
the named column from each dataframe and return a matrix.
}
\details{
The common use case for this is to provide a list of topTable
data frames and extract one column from each file to create
a matrix of LogRatios or Pvalues.

Technically, this should work as long as the requested colName is present
in each dataframe.  The default robust = TRUE should be used unless you
are absolutely certain each dataframe in the input list has the same row count
and row order.
}
\examples{

  MyPvalues  = ExtractCol (TopTableList, colName="P.Value")

}
\author{
John Thompson, \email{john.thompson@bms.com}
}
\keyword{topTable}
