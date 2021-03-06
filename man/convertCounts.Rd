% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convertCounts.R
\name{convertCounts}
\alias{convertCounts}
\title{Function convertCounts}
\usage{
convertCounts(
  counts,
  unit,
  geneLength,
  log = FALSE,
  normalize = "none",
  prior.count = NULL,
  debug = FALSE
)
}
\arguments{
\item{counts}{A numeric matrix or dataframe of N genes x M Samples.  All columns
must be numeric.}

\item{unit}{Required. One of CPM, FPKM, FPK or TPM.}

\item{geneLength}{A vector or matrix of gene lengths. Required for length-normalizes units (TPM, FPKM or FPK).
If you supply a matrix, the rowMeans are calculated and used.}

\item{log}{Default = FALSE.  Set TRUE to return Log2 values.
Employs edgeR functions which use an prior.count of 0.25 scaled by the library size.}

\item{normalize}{Default = "none". Other options: "TMM", "RLE", "upperquartile"
Invokes edgeR::calcNormFactors for normalization. Upperquartile uses the 75th percentile.  Normalize settings are case insensitive.}

\item{prior.count}{Average count to be added to each observation to avoid taking log of zero. Used only if log=TRUE. (Default dependent on method;
0 for TPM, 0.25 for CPM and FPKM)
The prior.count is passed to edgeR cpm and rpkm functions and applies to logTPM, logCPM and logFPKM calculations.}
}
\value{
A matrix in the new unit space
}
\description{
Takes a count matrix as input and converts to other desired units.  Supported
units include CPM, FPKM, FPK and TPM.  Output units can be logged
and/or normalized.  Calculations are performed using edgeR functions except
for the conversion to TPM which is converted from FPKM using the formula provided
by [Harold Pimental](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/).
}
\details{
geneLength is a vector where length(geneLength) == nrow(counts).  If you pass
an RSEM effectiveLength matrix, rowMeans(effectiveLength) is used (because edgeR functions
only accept a vector for effectiveLength).

Note that log2 values for CPM, TPM and FPKM employ edgeR's prior.count handling to avoid divide by zero.
}
\examples{
#TMM normalized Log2FPKM
Log2FPKM = convertCounts(mycounts),
                      unit="fpkm",
                      geneLength=gene.annotation$ExonLength,
                      log=TRUE,
                      normalize="tmm")

#un-normalized CPM (not logged)
RawCPM = convertCounts(MyCounts,
                      unit="CPM",
                      log=FALSE,
                      normalize="none")

}
\author{
John Thompson, \email{john.thompson@bms.com}
}
\keyword{RNA-Seq,}
\keyword{conversion}
\keyword{unit}
