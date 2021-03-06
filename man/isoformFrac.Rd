% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isoformFrac.R
\name{isoformFrac}
\alias{isoformFrac}
\title{Function isoformFrac}
\usage{
isoformFrac(dgeObj, dataType = "fpkm", normalize = "tmm")
}
\arguments{
\item{dgeObj}{An isoform level DGEobj created by function initDGEobj,
Xpress2DGEO or OmicsoftToDgeObj.  Counts and isoformData must be present
in the DGEobj (required).  isoformData$ExonLength must be present or assay =
"effectiveLength" must be present.}

\item{dataType}{One of "fpkm" or "tpm" (default="fpkm")}

\item{normalize}{Default = "TMM" and invokes TMM normalization. Other allowed
values are: "RLE","upperquartile", "none". Invokes edgeR::calcNormFactors for
normalization.  Only invoked when dataType="fpkm".  This is because
applying TPM essentially erases any prior column scaling so TMM and similar
normalizations have no effect.}
}
\value{
A DGEobj with an isoform fraction dataframe added
}
\description{
Takes a dgeObj as input (transcript level data) and returns a matrix
containing isoform fraction data.
}
\details{
Isoform Fraction is calculated using length normalized data (FPKM or TPM), as
length normalized data is required because different isoforms have different
total exon lengths.  If FPKM is specified, you can also specify a
normalization (via edgeR::calcNormFactors). Isoform fraction is calculated
simply as the isoform intensity divided by the summed gene intensity for all
isoforms of a given gene.

TPM or FPKM are calculated directly from counts using all data in the dgeObj.
I recommend performing low intensity filtering at the gene level before
running isoformFrac.
}
\examples{
   MyDgeObj <- isoformFrac(MyDgeObj)

}
\author{
John Thompson, \email{john.thompson@bms.com}
}
\keyword{DGEobj,}
\keyword{RNA-Seq,}
\keyword{fraction}
\keyword{isoform}
