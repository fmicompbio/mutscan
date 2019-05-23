#' Filter FASTQ files from CIS or TRANS experiments
#'
#' @param L SummarizedExperiment, as output by \code{readFastqs}.
#' @param avePhredMin numeric(1) Minimum average Phred score in the variable
#'   region for a read to be retained. If L contains both forward and reverse
#'   variable regions, the minimum average Phred score has to be achieved in
#'   both for a read pair to be retained.
#' @param variableNMax numeric(1) Maximum number of Ns in the variable region
#'   for a read to be retained.
#' @param umiNMax numeric(1) Maximum number of Ns in the UMI for a read to be
#'   retained.
#' @param wildTypeForward,wildTypeReverse character(1), the wild type sequence
#'   for the forward and reverse variable region.
#' @param nbrMutatedCodonsMax numeric(1) Maximum number of mutated codons that
#'   are allowed.
#' @param forbiddenMutatedCodon character(1) A codon (can contain ambiguous
#'   IUPAC characters, see \code{\link[Biostrings]{IUPAC_CODE_MAP}}). If a read
#'   pair contains a mutated codon matching this pattern, it will be filtered
#'   out.
#'
#' @return A filtered SummarizedExperiment object
#'
#' @export
#'
#' @importFrom Biostrings quality vcountPattern vmatchPattern
#' @importFrom ShortRead alphabetScore FastqQuality
#' @importFrom BiocGenerics width
#' @importFrom matrixStats rowMins
#' @importFrom S4Vectors elementNROWS endoapply DataFrame
#' @importFrom IRanges start
#' @importFrom SummarizedExperiment assay assayNames rowData
#'
#' @author Charlotte Soneson
#'   
filterReads <- function(L, avePhredMin = 20, variableNMax = 0, umiNMax = 0,
                        wildTypeForward = NULL, wildTypeReverse = NULL, 
                        nbrMutatedCodonsMax = 1, forbiddenMutatedCodon = NULL) {
  ## --------------------------------------------------------------------------
  ## Pre-flight checks
  ## --------------------------------------------------------------------------
  isValidL(L)
  nbrReads <- length(assay(L, "variableSeqForward")$seq)
  
  ## --------------------------------------------------------------------------
  ## Calculate average Phred score
  ## --------------------------------------------------------------------------
  keepAvePhredVariableForward <- keepAvePhredVariableReverse <- rep(TRUE, nbrReads)
  if (!is.null(avePhredMin)) {
    if ("variableSeqForward" %in% SummarizedExperiment::assayNames(L)) {
      avePhredVariableForward <- ShortRead::alphabetScore(
        Biostrings::quality(assay(L, "variableSeqForward")$seq)
      )/BiocGenerics::width(assay(L, "variableSeqForward")$seq) 
      keepAvePhredVariableForward <- avePhredVariableForward >= avePhredMin
    }
    if ("variableSeqReverse" %in% SummarizedExperiment::assayNames(L)) {
      avePhredVariableReverse <- ShortRead::alphabetScore(
        Biostrings::quality(assay(L, "variableSeqReverse")$seq)
      )/BiocGenerics::width(assay(L, "variableSeqReverse")$seq)
      keepAvePhredVariableReverse <- avePhredVariableReverse >= avePhredMin
    }
  }
  L$nbrReadPairsLowAvePhredVariableForward <- sum(!keepAvePhredVariableForward)
  L$nbrReadPairsLowAvePhredVariableReverse <- sum(!keepAvePhredVariableReverse)
  
  ## --------------------------------------------------------------------------
  ## Count the number of Ns in variable regions
  ## --------------------------------------------------------------------------
  keepVariableNbrNForward <- keepVariableNbrNReverse <- rep(TRUE, nbrReads)
  if (!is.null(variableNMax)) {
    if ("variableSeqForward" %in% SummarizedExperiment::assayNames(L)) {
      variableNbrNForward <- Biostrings::vcountPattern(
        pattern = "N", subject = assay(L, "variableSeqForward")$seq,
        max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, 
        fixed = TRUE, algorithm = "auto")
      keepVariableNbrNForward <- variableNbrNForward <= variableNMax
    }
    if ("variableSeqReverse" %in% SummarizedExperiment::assayNames(L)) {
      variableNbrNReverse <- Biostrings::vcountPattern(
        pattern = "N", subject = assay(L, "variableSeqReverse")$seq,
        max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, 
        fixed = TRUE, algorithm = "auto")
      keepVariableNbrNReverse <- variableNbrNReverse <= variableNMax
    }
  }
  L$nbrReadPairsTooManyNVariableForward <- sum(!keepVariableNbrNForward)
  L$nbrReadPairsTooManyNVariableReverse <- sum(!keepVariableNbrNReverse)
  
  ## --------------------------------------------------------------------------
  ## Count the number of Ns in UMIs
  ## --------------------------------------------------------------------------
  keepUMINbrN <- rep(TRUE, nbrReads)
  if (!is.null(umiNMax)) {
    if ("umis" %in% SummarizedExperiment::assayNames(L)) {
      umiNbrN <- Biostrings::vcountPattern(
        pattern = "N", subject = assay(L, "umis")$seq,
        max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, 
        fixed = TRUE, algorithm = "auto")
      keepUMINbrN <- umiNbrN <= umiNMax
    }
  }
  L$nbrReadPairsTooManyNUMI <- sum(!keepUMINbrN)
  
  ## --------------------------------------------------------------------------
  ## Find positions where variable region has mismatches
  ## --------------------------------------------------------------------------
  keepWildTypeForward <- keepWildTypeReverse <- rep(TRUE, nbrReads)
  qualMutatedForward <- qualMutatedReverse <- NULL
  basesWithMutationsForward <- basesWithMutationsReverse <- NULL
  codonsWithMutationsForward <- codonsWithMutationsReverse <- NULL
  if (!is.null(wildTypeForward)) {
    if ("variableSeqForward" %in% SummarizedExperiment::assayNames(L)) {
      filterForward <- mismatchFilter(wildTypeSeq = wildTypeForward, 
                                      variableSeq = assay(L, "variableSeqForward")$seq, 
                                      nbrMutatedCodonsMax = nbrMutatedCodonsMax,
                                      forbiddenMutatedCodon = forbiddenMutatedCodon)
      SummarizedExperiment::rowData(L)$qualMutatedForward <- filterForward$qualMutated
      SummarizedExperiment::rowData(L)$basesWithMutationsForward <- filterForward$basesWithMutations
      SummarizedExperiment::rowData(L)$codonsWithMutationsForward <- filterForward$codonsWithMutations
      SummarizedExperiment::rowData(L)$encodedMutatedCodonsForward <- filterForward$encodedMutCodons
      keepWildTypeForward <- filterForward$keepMM
    }
  }
  if (!is.null(wildTypeReverse)) {
    if ("variableSeqReverse" %in% SummarizedExperiment::assayNames(L)) {
      filterReverse <- mismatchFilter(wildTypeSeq = wildTypeReverse, 
                                      variableSeq = assay(L, "variableSeqReverse")$seq, 
                                      nbrMutatedCodonsMax = nbrMutatedCodonsMax,
                                      forbiddenMutatedCodon = forbiddenMutatedCodon)
      SummarizedExperiment::rowData(L)$qualMutatedReverse <- filterReverse$qualMutated
      SummarizedExperiment::rowData(L)$basesWithMutationsReverse <- filterReverse$basesWithMutations
      SummarizedExperiment::rowData(L)$codonsWithMutationsReverse <- filterReverse$codonsWithMutations
      SummarizedExperiment::rowData(L)$encodedMutatedCodonsReverse <- filterReverse$encodedMutCodons
      keepWildTypeReverse <- filterReverse$keepMM
    }
  }
  L$nbrReadPairsCodonFilterForward <- sum(!keepWildTypeForward)
  L$nbrReadPairsCodonFilterReverse <- sum(!keepWildTypeReverse)
  
  ## --------------------------------------------------------------------------
  ## Return filtered object
  ## --------------------------------------------------------------------------
  toKeep <- keepAvePhredVariableForward & 
    keepAvePhredVariableReverse & 
    keepVariableNbrNForward & 
    keepVariableNbrNReverse & 
    keepUMINbrN & 
    keepWildTypeForward & 
    keepWildTypeReverse
  
  L$totalNbrReadPairsPassedFilters <- sum(toKeep)
  return(L[toKeep, ])

}

#' Filter reads based on mismatches with the wildtype sequence
#'
#' @param wildTypeSeq character(1), the wild type sequence
#' @param variableSeq QualityScaledDNAStringSet, containing the sequences of the
#'   variable region
#' @param nbrMutatedCodonsMax numeric(1) Maximum number of mutated codons that
#'   are allowed.
#' @param forbiddenMutatedCodon character(1) A codon (can contain ambiguous
#'   IUPAC characters, see \code{\link[Biostrings]{IUPAC_CODE_MAP}}). If a read
#'   pair contains a mutated codon matching this pattern, it will be filtered
#'   out.
#'
#' @author Charlotte Soneson
#' 
#' @keywords internal
#' 
#' @return A list with four elements:
#'   \describe{
#'     \item{qualMutated}{IntegerList with qualities for the mutated bases}
#'     \item{keepMM}{Logical vector, whether to keep each read}
#'     \item{basesWithMutations}{IntegerList, mutated base positions}
#'     \item{codonsWithMutations}{NumericList, mutated codons}
#'  }
#' 
#' @importFrom Biostrings quality vmatchPattern
#' @importFrom ShortRead FastqQuality
#' @importFrom matrixStats rowMins
#' @importFrom S4Vectors elementNROWS endoapply
#' @importFrom IRanges start
#' 
mismatchFilter <- function(wildTypeSeq, variableSeq, nbrMutatedCodonsMax = NULL,
                           forbiddenMutatedCodon = NULL) {
  keepMM <- rep(TRUE, length(variableSeq))
  
  ## --------------------------------------------------------------------------
  ## Find mismatch positions and qualities of mismatched bases
  ## --------------------------------------------------------------------------
  mmp <- findMismatchPositions(pattern = wildTypeSeq, 
                               subject = variableSeq)
  qualMutated <- as(Biostrings::quality(variableSeq)[mmp$nucleotideMismatches], "IntegerList")

  ## --------------------------------------------------------------------------
  ## Count number of mutated codons per read
  ## --------------------------------------------------------------------------
  if (!is.null(nbrMutatedCodonsMax)) {
    nbrMutCodons <- S4Vectors::elementNROWS(unique(mmp$codonMismatches))
    keepMM <- nbrMutCodons <= nbrMutatedCodonsMax
  }
  
  ## --------------------------------------------------------------------------
  ## Derive a characterisation of each sequence, in terms of the codons with 
  ## mismatches
  ## --------------------------------------------------------------------------
  uniqueMutCodons <- unique(mmp$codonMismatches)
  mutCodons <- S4Vectors::endoapply(uniqueMutCodons, 
                                    function(v) rep(3 * v, each = 3) + c(-2, -1, 0))
  seqMutCodons <- variableSeq[mutCodons]
  splitMutCodons <- strsplit(gsub("([^.]{3})", "\\1\\.", as.character(seqMutCodons)), ".", fixed = TRUE)
  encodedMutCodons <- BiocGenerics::mapply(uniqueMutCodons, splitMutCodons, FUN = function(a, b) {
    if (length(a) == 0 || length(b) == 0) {
      return("WT")
    } else {
      out <- paste0("c", a, b)
      return(paste(out, collapse = "_"))
    }
  })
  
  ## --------------------------------------------------------------------------
  ## Find any occurrences of the forbidden mutated codon
  ## --------------------------------------------------------------------------
  if (!is.null(forbiddenMutatedCodon)) {
    matchForbidden <- Biostrings::vmatchPattern(
      pattern = forbiddenMutatedCodon, 
      as(seqMutCodons, "DNAStringSet"), 
      fixed = FALSE)
    keepMM <- keepMM & !vapply(matchForbidden, function(v) {
      any(IRanges::start(v) %% 3 == 1)
    }, FALSE)
  }
  
  list(qualMutated = qualMutated, keepMM = keepMM,
       basesWithMutations = mmp$nucleotideMismatches,
       codonsWithMutations = mmp$codonMismatches,
       encodedMutCodons = encodedMutCodons)
}