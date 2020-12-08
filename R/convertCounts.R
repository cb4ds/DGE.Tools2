#' Convert count matrix to CPM, FPKM, FPK, or TPM
#'
#' Takes a count matrix as input and converts to other desired units.  Supported
#' units include CPM, FPKM, FPK, and TPM.  Output units can be logged
#' and/or normalized.  Calculations are performed using edgeR functions except
#' for the conversion to TPM which is converted from FPKM using the formula provided
#' by [Harold Pimental](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/).
#'
#' geneLength is a vector where length(geneLength) == nrow(counts). If a RSE effectiveLength
#' matrix is passed as input, rowMeans(effectiveLength) is used (because edgeR functions
#' only accept a vector for effectiveLength).
#'
#' Note that log2 values for CPM, TPM, and FPKM employ edgeR's prior.count handling to avoid divide by zero.
#'
#' @param counts A numeric matrix or dataframe of N genes x M Samples.  All columns
#' must be numeric.
#' @param unit  Required. One of CPM, FPKM, FPK or TPM.
#' @param geneLength A vector or matrix of gene lengths. Required for length-normalized units (TPM, FPKM or FPK).
#'    If geneLength is a matrix, the rowMeans are calculated and used.
#' @param log Default = FALSE.  Set TRUE to return Log2 values.
#'    Employs edgeR functions which use an prior.count of 0.25 scaled by the library size.
#' @param normalize Default = "none". Other options: "TMM", "RLE", "upperquartile"
#'  Invokes edgeR::calcNormFactors() for normalization. Upperquartile uses the 75th percentile.  Normalize settings are case insensitive.
#' @param prior.count Average count to be added to each observation to avoid taking log of zero.
#'  Used only if log = TRUE. (Default dependent on method; 0 for TPM, 0.25 for CPM and FPKM)
#'  The prior.count is passed to edgeR cpm and rpkm functions and applies to logTPM, logCPM, and logFPKM calculations.
#'
#' @return A matrix in the new unit space
#'
#' @examples
#' \dontrun{
#'     # TMM normalized Log2FPKM
#'     Log2FPKM <- convertCounts(mycounts,
#'                               unit = "fpkm",
#'                               geneLength = gene.annotation$ExonLength,
#'                               log = TRUE,
#'                               normalize = "tmm")
#'
#'     # Non-normalized CPM (not logged)
#'     RawCPM <- convertCounts(MyCounts,
#'                             unit = "CPM",
#'                             log = FALSE,
#'                             normalize = "none")
#' }
#'
#' @import magrittr
#' @importFrom edgeR cpm rpkm expandAsMatrix calcNormFactors DGEList
#' @importFrom assertthat assert_that
#'
#' @export
convertCounts <- function(counts,
                          unit,
                          geneLength,
                          log = FALSE,
                          normalize = "none",
                          prior.count = NULL) {

    assertthat::assert_that(!(nrow(counts) == 0),
                            msg = "counts must be specified.")
    assertthat::assert_that(!is.null(unit),
                            msg = "unit must be specified.")

    unit <- toupper(unit)
    if (unit %in% c('FPKM', 'TPM', 'FPK')) {
        # In these cases geneLength is required
        assertthat::assert_that(!missing(geneLength),
                                msg = "geneLength must be specified when unit is 'FPK', 'FPKM', or 'TPM.'")

        if ("matrix" %in% class(geneLength)) {# Flatten to a vector
            geneLength <- rowMeans(geneLength, na.rm = TRUE)
        }
    }
    # Make normalize method case insensitive (calcNormFactors is case sensitive)
    if (toupper(normalize) %in% c("TMM", "RLE")) {
        normalize <- toupper(normalize)
    }
    if (toupper(normalize) %in% c("UPPERQUARTILE", "NONE")) {
        normalize <- tolower(normalize)
    }
    if (toupper(normalize) %in% c("TMMWZP")) {
        normalize <- "TMMwzp"
    }


    # Coerce counts to a matrix
    result <- counts <- as.matrix(counts)
    assertthat::assert_that("matrix" %in% class(counts),
                            msg = "counts must be able to be coerced to a matrix.")

    # Make sure geneLength is correct length
    if (!missing(geneLength)) {
        assertthat::assert_that(length(geneLength) == nrow(counts),
                                msg = "geneLength must be the same length of the number of rows in counts.")
    }

    # Set defaults
    if (missing(log)) {
        log = FALSE
    }
    if (missing(normalize)) {
        normalize = 'none'
    }
    if (is.logical(normalize)) { # Don't encourage logicals; here for backward compatibility
        if (normalize == TRUE) {
            normalize <- 'TMM'
        }
        if (normalize == FALSE) {
            normalize <- 'none'
        }
    }

    if (is.null(prior.count)) {
        if (log == FALSE) {
            prior.count <- 0 # Not used when log = F
        }
        else if (unit == "TPM") {
            prior.count <- 0
        }
        else {
            prior.count <- 0.25
        }
    }

    result <- switch(toupper(unit),
                     "CPM" = calcCPM(counts, log, normalize, prior.count),
                     "FPKM" = calcFPKM(counts, log, normalize, geneLength, prior.count),
                     "FPK" = calcFPK(counts, log, normalize, geneLength, prior.count),
                     "TPM" = calcTPM(counts, log, normalize, geneLength, prior.count)
    )
    return(result)
}


#' Calculate TPM for a subsetted DGEobj
#'
#' Calculates TPM for a heavily subsetted DGEobj. The function will calculate TPM
#' using the original data but returns a DGEobj with the subset.
#'
#' Internally, convertCounts uses edgeR::fpkm() to calculate FPKM and converts to TPM
#' using the formula provided by [Harold Pimental](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/).
#'
#' @param dgeObj A DGEobj data structure
#' @param applyFilter Default = TRUE. If TRUE, reduces to the filtered gene list. FALSE returns
#'   all genes in the raw data.
#'
#' @return A matrix in the new unit space
#'
#' @examples
#' \dontrun{
#'    myTPM <- tpm(DGEobj)
#' }
#' @import DGEobj magrittr
#' @importFrom assertthat assert_that
#' @export
tpm.on.subset <- function(dgeObj, applyFilter = TRUE){

    assertthat::assert_that("DGEobj" %in% class(dgeObj),
                            msg = "dgeObj should be of class 'DGEobj'.")
    assertthat::assert_that(attr(dgeObj, "level") %in% c("isoform", "gene"),
                            msg = "The level of dgeObj should be of type 'isoform' or type 'gene'.")

    # Default to gene level
    level <- "gene"
    if (attr(dgeObj, "level") == "isoform") {
        level <- "isoform"
    }

    if (level == "gene") {
        rowdata <- getItem(dgeObj, "geneData_orig")
    } else {
        rowdata <- getItem(dgeObj, "isoformData_orig")
    }

    # Get geneLength depending on source data
    if (attr(dgeObj, "source") == "Omicsoft") { # Omicsoft data
        geneLength <- rowdata$ExonLength
    } else if ("effectiveLength_orig" %in% names(dgeObj)) { # Use rowMeans(effectiveLength)
        geneLength <- rowMeans(getItem(dgeObj, "effectiveLength_orig"), na.rm = TRUE)
    }

    TPM <- convertCounts(getItem(dgeObj, "counts_orig"),
                         geneLength = geneLength,
                         unit = "tpm",
                         log = FALSE,
                         normalize = FALSE)

    # Remove filtered out genes
    if (applyFilter == TRUE) {
        idx <- rownames(TPM) %in% rownames(getItem(dgeObj, "counts"))
        TPM <- TPM[idx,]
    }
    return(TPM)
}


#' Convert counts and geneLength to TPM units
#'
#' Takes a counts and geneLength as input and converts to TPM units using the equation from
#' [Harold Pimental](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/).
#'
#' The result should be the same as using convertCounts with normalize = 'tpm' and log = FALSE.
#'
#' geneLength can be a vector (length == nrow(counts)) or a matrix (same dim as counts).
#' The geneLength is used as is, or optionally collapsed to a vector by rowMeans.
#'
#' @param counts A numeric matrix of N genes x M samples. All columns must be numeric.
#' @param geneLength Numeric matrix of gene lengths. Often the ExonLength item of a DGEobj.
#' @param collapse Default = FALSE. TRUE or FALSE determines whether to use rowMeans on the geneLength matrix.
#'
#' @return A matrix of TPM values
#'
#' @examples
#' \dontrun{
#'   myTPM <- tpm.direct(myCounts, myGeneLength)
#' }
#'
#' @import magrittr
#' @importFrom edgeR expandAsMatrix
#' @importFrom assertthat assert_that
#'
#' @export
tpm.direct <- function(counts,
                       geneLength,
                       collapse = FALSE) {

    if (!is.matrix(counts)) {
        result <- counts <- as.matrix(counts)
        assertthat::assert_that("matrix" %in% class(result),
                                msg = "counts must be able to be coerced to a matrix.")
    }

    if (is.vector(geneLength)) {
        assertthat::assert_that(length(geneLength) == nrow(counts),
                                msg = "geneLength should be of the same length as the number of rows in counts.")
    } else {
        if (!is.matrix(geneLength)) {
            result <- geneLength <- as.matrix(geneLength)
            assertthat::assert_that("matrix" %in% class(result),
                                    msg = "geneLength must be able to be coerced to a matrix.")
            assertthat::assert_that(all(dim(counts) == dim(geneLength)),
                                    msg = "The dimensions of counts and geneLength should match.")
        }
    }

    if (collapse & is.matrix(geneLength)) {
        geneLength <- rowMeans(geneLength, na.rm = TRUE)
    }

    # Calculation  (fpk / colsum(fpk) ) * 10e6
    fpb <- counts / geneLength
    sumfpb <- colSums(fpb)
    tpm <- fpb / edgeR::expandAsMatrix(sumfpb, byrow = TRUE, dim = dim(fpb)) * 1e6
}

# Helper Functions
calcCPM <- function(counts, log, normalize, prior.count){
    counts %>%
        edgeR::DGEList() %>%
        edgeR::calcNormFactors(method = normalize) %>%
        edgeR::cpm(log = log, prior.count = prior.count)
}

calcFPKM <- function(counts, log, normalize, geneLength, prior.count){
    counts %>%
        edgeR::DGEList() %>%
        edgeR::calcNormFactors(method = normalize) %>%
        edgeR::rpkm(log = log, gene.length = geneLength, prior.count = prior.count)
}

calcTPM <- function(counts, log, normalize, geneLength, prior.count){
    if (normalize != "none") {
        warning(paste('TPM normalization overides', normalize, 'normalization!'))
    }
    if (prior.count != 0 && log == TRUE) {
        warning("Using a prior.count for logTPM calculations is not recommended and may produce unpredictable results!")
    }

    fpkm <- calcFPKM(counts, log = log, normalize = normalize,
                     geneLength = geneLength, prior.count = prior.count)

    # Helper function
    fpkmToTpm <- function(fpkm) {
        colSumMat <- edgeR::expandAsMatrix(colSums(fpkm, na.rm = TRUE), byrow = TRUE, dim = dim(fpkm))
        fpkm / colSumMat * 1e6
    }

    if (log == FALSE) {
        TPM <- fpkmToTpm(fpkm)
    } else {
        TPM <- log2(fpkmToTpm(2^fpkm))
    }
    return(TPM)
}

calcFPK <- function(counts, log, normalize, geneLength, prior.count){
    if (tolower(normalize) == 'none') {
        # Check for zero geneLength just in case
        if (min(geneLength) == 0) {
            geneLength <- geneLength + 1
        }

        FPK <- counts / (geneLength / 1000)

        if (log == TRUE) {
            FPK <- log2(FPK + prior.count)
        }
    }

    return(FPK)
}
