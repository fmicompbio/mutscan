#' Get error statistics from constant sequences
#' 
#' From the constant sequence parts, calculate sequencing error statistics by
#' comparison to expected sequences.
#' 
#' @param L SummarizedExperiment, as output by \code{readFastqs}.
#' @param expectedForward,expectedReverse character(1) with the expected sequences
#'     of the forward or reverse constant sequences. If \code{NULL}, the most
#'     frequent sequence will be used.
#' @param verbose logical(1), whether to print out progress messages.
#'
#' @return A list with 6 elements:
#'     \describe{
#'     \item{propErrorsConstantF}{proportion of mismatching bases in forward constant sequence}
#'     \item{propErrorsConstantR}{proportion of mismatching bases in reverse constant sequence}
#'     \item{QcountConstantErrorsF}{numbers of mismatching bases per quality score in forward constant sequence}
#'     \item{QcountConstantErrorsR}{numbers of mismatching bases per quality score in reverse constant sequence}
#'     \item{QcountConstantCorrectF}{numbers of matching bases per quality score in forward constant sequence}
#'     \item{QcountConstantCorrectR}{numbers of matching bases per quality score in reverse constant sequence} }
#'
#' @importFrom BiocGenerics width
#' @importFrom ShortRead tables
#' @importFrom methods is
#'  
#' @author Michael Stadler
#' 
#' @export
calcErrorStats <- function(L, expectedForward = NULL, expectedReverse = NULL,
                           verbose = FALSE) {
  ## --------------------------------------------------------------------------
  ## Pre-flight checks
  ## --------------------------------------------------------------------------
  isValidL(L)
  if (!is.null(expectedForward) && (!is(expectedForward, "character") ||
      length(expectedForward) != 1L)) {
    stop("'expectedForward' must be length-1 character or NULL.")
  }
  if (!is.null(expectedReverse) && (!is(expectedReverse, "character") ||
      length(expectedReverse) != 1L)) {
    stop("'expectedReverse' must be length-1 character or NULL.")
  }
  numberReads <- length(assay(L, "constantSeqForward")$seq)
  constantLengthF <- width(assay(L, "constantSeqForward")$seq)[1]
  constantLengthR <- width(assay(L, "constantSeqReverse")$seq)[1]
  if (any(width(assay(L, "constantSeqForward")$seq) != constantLengthF)) {
    stop("'assay(L, constantSeqForward)' are not all of the same width")
  }
  if (any(width(assay(L, "constantSeqReverse")$seq) != constantLengthR)) {
    stop("'assay(L, constantSeqReverse)' are not all of the same width")
  }
  
  ## --------------------------------------------------------------------------
  ## Expected sequences
  ## --------------------------------------------------------------------------
  if (is.null(expectedForward)) {
    if (verbose) {
      message("identifying expectedForward...", appendLF = FALSE)
    }
    tabF <- ShortRead::tables(assay(L, "constantSeqForward")$seq)
    expectedForward <- names(tabF$top)[1]
    if (verbose) {
      message("done: ", expectedForward, " (", round(100 * tabF$top[1] / numberReads, 1), "%)")
    }
    if (tabF$top[1]/numberReads < 0.5) {
      warning("Auto-detected 'expectedForward' is less than half of the reads - maybe you should set it manually.")
    }
  }
  if (is.null(expectedReverse)) {
    if (verbose) {
      message("identifying expectedReverse...", appendLF = FALSE)
    }
    tabR <- ShortRead::tables(assay(L, "constantSeqReverse")$seq)
    expectedReverse <- names(tabR$top)[1]
    if (verbose) {
      message("done: ", expectedReverse, " (", round(100 * tabR$top[1] / numberReads, 1), "%)")
    }
    if (tabR$top[1]/numberReads < 0.5) {
      warning("Auto-detected 'expectedReverse' is less than half of the reads - maybe you should set it manually.")
    }
  }
  
  ## --------------------------------------------------------------------------
  ## Identify and tabulate mismatches
  ## --------------------------------------------------------------------------
  if (verbose) {
    message("analyzing mismatches in constantSeqForward")
  }
  mmF <- tabulateQualitiesByMatchstate(expectedForward, assay(L, "constantSeqForward")$seq)
  if (verbose) {
    message("analyzing mismatches in constantSeqReverse")
  }
  mmR <- tabulateQualitiesByMatchstate(expectedReverse, assay(L, "constantSeqReverse")$seq)

  list(propErrorsConstantF = sum(mmF$error) / constantLengthF / numberReads,
       propErrorsConstantR = sum(mmR$error) / constantLengthR / numberReads,
       QcountConstantErrorsF = mmF$error,
       QcountConstantErrorsR = mmR$error,
       QcountConstantCorrectF = mmF$correct,
       QcountConstantCorrectR = mmR$correct)
}
