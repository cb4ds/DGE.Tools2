#' Run a typical voom workflow
#'
#' This function runs several steps in the RNA-Seq pipeline in one fell swoop.
#' It supports duplicateCorrelation a blocking vector is provided. It
#' includes low intensity filtering, filtering for protein coding genes, and
#' filters out zero effective length genes (genes shorter than library size).
#' Then it runs TMM normalization, voomWithQualityWeights, lmFit, and eBayes.
#'
#' To incorporate SVA analysis, run SVA first and then add the SVA
#' variables to the design table. After those steps, voomWorkflow can be run.
#'
#' After running this function, define the contrasts and execute runContrast to
#' complete DGE calculations.
#'
#' @param dgeObj A class DGEobj with counts, gene annotation, and sample
#'   annotation.
#' @param formula A text representation of the desired formula (input of class
#'   character, not formula.)
#' @param designMatrixName User defined name for the design matrix.
#' @param dupCorBlock A blocking vector to define which samples belong to the
#'   same subject to be used with the duplicateCorrelation function.
#' @param annotationFile Text file of key/value pairs to populate DGEobj
#'   attributes. (Optional but highly advised)
#' @param proteinCodingOnly Set to TRUE to keep only protein coding genes.
#'   (Default = TRUE)
#' @param sampleFraction Fraction of samples that must meet intensity thresholds
#'   to keep a gene. (Default = 0.5)
#' @param ... Additional named arguments passed to the low intensity filter
#'   function to define the desired intensity filter type (see ?lowIntFilter). Settable
#'   arguments for low intensity filtering are: fracThreshold (Default = 0.5),
#'   countThreshold, fpkThreshold, zfpkmThreshold, tpmThreshold.  countThreshold plus one
#'   other argument can be used.  If no arguments supplied here, the
#'   following defaults apply the following: fracThreshold = 0.5, countThreshold = 10,
#'   fpkThreshold = 5. To disable intensity filtering use sampleFraction = 0.
#'
#' @return A DGEobj with analysis results added.
#'
#' @examples
#' \dontrun{
#'    myDGEobj <- voomWorkflow(myDGEobj,
#'                             formula = "~ 0 + ReplicateGroup",
#'                             designMatrixName = "ReplicateGroup",
#'                             annotationFile = "MyProjectName.txt",
#'                             proteinCodingOnly = TRUE)
#' }
#'
#' @import magrittr
#' @import DGEobj
#' @importFrom assertthat assert_that
#' @importFrom stats as.formula
#'
#' @export
voomWorkflow <- function(dgeObj,
                         formula,
                         designMatrixName,
                         dupCorBlock,
                         annotationFile,
                         proteinCodingOnly = FALSE,
                         sampleFraction = 0.5,
                         ...) {

    ellipsisArgs <- list(...)

    assertthat::assert_that(!missing(dgeObj),
                            "DGEobj" %in% class(dgeObj),
                            msg = "dgeObj must be specified and be of class 'DGEobj'.")
    assertthat::assert_that(!missing(formula),
                            "character" %in% class(formula),
                            msg = "formula must be specified and be of class 'character'.")
    assertthat::assert_that(!missing(designMatrixName),
                            "character" %in% class(designMatrixName),
                            msg = "designMatrixName must be specified and be of class 'character'.")
    assertthat::assert_that(!(length(ellipsisArgs) > 2),
                            msg = "No more than 2 intensity filtering arguments are allowed.")
    if (length(ellipsisArgs) == 2) {
        assertthat::assert_that("countThreshold" %in% names(ellipsisArgs),
                                msg = "When two intensity filtering arguments are supplied, one must be 'countThreshold'.")
    }

    # Add project metadata
    if (!missing(annotationFile)) {
        dgeObj <- DGEobj::annotateDGEobj(dgeObj, annotations = annotationFile)
    }

    # Filter out genes with any zero efflength
    result <- try({el <- DGEobj::getItem(dgeObj, "effectiveLength")}, silent = TRUE)
    if (!("try-error" %in% class(result)) && !is.null(el) ) {
        rowmin <- apply(el, 1, min)
        idx <- rowmin > 0
        dgeObj <- dgeObj[idx,]
    }

    # Intensity filtering args from ellipsis
    thresholds <- list()
    if ("fpkThreshold" %in% names(ellipsisArgs)) {
        thresholds$fpkThreshold <- ellipsisArgs$fpkThreshold
    }

    if ("countThreshold" %in% names(ellipsisArgs)) {
        thresholds$countThreshold <- ellipsisArgs$countThreshold
    }

    if ("zfpkmThreshold" %in% names(ellipsisArgs)) {
        thresholds$zfpkmThreshold <- ellipsisArgs$zfpkmThreshold
    }

    if ("tpmThreshold" %in% names(ellipsisArgs)) {
        thresholds$tpmThreshold <- ellipsisArgs$tpmThreshold
    }

    # Construct filtering command
    cmd <- ("dgeObj <- lowIntFilter(dgeObj, sampleFraction = sampleFraction, ")
    for (i in 1:length(thresholds)) {
        cmd <- stringr::str_c(cmd, names(thresholds)[i], " = ", thresholds[[i]])
        if (i < length(thresholds)) {
            cmd <- stringr::str_c(cmd, ", ")
        }
    }
    cmd <- stringr::str_c(cmd, ")")
    eval(parse(text = cmd))

    # Keep only protein coding genes
    if (proteinCodingOnly == TRUE) {
        idx <- NULL
        if ("Source" %in% colnames(dgeObj$geneData)) {
            idx <- dgeObj$geneData$Source == "protein_coding"
        } else if ("Source" %in% colnames(dgeObj$isoformData)) {
            idx <- dgeObj$isoformData$Source == "protein_coding"
        }

        if (is.null(idx) == FALSE) {
            # Convert those to FALSE
            idx[is.na(idx)] <- FALSE
            dgeObj <- dgeObj[idx,]
        }
    }

    # Set up and save the design matrix
    designMatrix <- model.matrix(as.formula(formula), getItem(dgeObj, "design"))
    # Clean up problem characters in colnames
    colnames(designMatrix) <- make.names(colnames(designMatrix))
    # Capture the formula as an attribute
    designMatrix <- setAttributes(designMatrix, list(formula = formula))
    # Save the modified designMatrix to the DGEobj
    dgeObj <- DGEobj::addItem(dgeObj,
                              item = designMatrix,
                              itemName = designMatrixName,
                              itemType = "designMatrix",
                              parent = "design")
    # Run DGE calculations
    # Normalize
    dgeObj <- runEdgeRNorm(dgeObj,
                           normMethod = "TMM",
                           plotFile = FALSE)

    message("Running Voom... this takes some time...")
    # Voom/lmFit
    if (missing(dupCorBlock) || is.null(dupCorBlock)) {
        dgeObj <- runVoom(dgeObj, designMatrixName,
                          qualityWeights = TRUE,
                          mvPlot = TRUE)
    } else {
        dgeObj <- runVoom(dgeObj, designMatrixName,
                          qualityWeights = TRUE,
                          dupCorBlock = dupCorBlock,
                          mvPlot = TRUE)
    }
}
