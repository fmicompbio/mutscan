#' Filter FASTQ files from CIS or TRANS experiments
#'
#' @param L List, as output by \code{readTransFastqs}.
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
#' @return A filtered list with eight elements: 
#'   \describe{ 
#'     \item{umis}{Merged forward and reverse UMI sequences}
#'     \item{constantSeqForward}{Constant forward sequence}
#'     \item{constantSeqReverse}{Constant reverse sequence}
#'     \item{variableSeqForward}{Variable forward sequence}
#'     \item{variableSeqReverse}{Variable reverse sequence}
#'     \item{minQualMutatedForward}{Minimal base quality of a mutated base in
#'     forward sequence} 
#'     \item{minQualMutatedReverse}{Minimal base quality of a
#'     mutated base in reverse sequence} 
#'     \item{readSummary}{data.frame tabulating the total number of read pairs,
#'     and the number that match any of the filtering criteria.} 
#'  }
#'
#' @export
#'
#' @importFrom Biostrings quality vcountPattern vmatchPattern
#' @importFrom ShortRead alphabetScore FastqQuality
#' @importFrom BiocGenerics width
#' @importFrom matrixStats rowMins
#' @importFrom S4Vectors elementNROWS endoapply
#' @importFrom IRanges start
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
  readSummary <- L$readSummary
  nbrReads <- length(L$variableSeqForward)
  
  ## --------------------------------------------------------------------------
  ## Calculate average Phred score
  ## --------------------------------------------------------------------------
  keepAvePhredVariableForward <- keepAvePhredVariableReverse <- rep(TRUE, nbrReads)
  if (!is.null(avePhredMin)) {
    if (!is.null(L$variableSeqForward)) {
      avePhredVariableForward <- ShortRead::alphabetScore(
        Biostrings::quality(L$variableSeqForward)
      )/BiocGenerics::width(L$variableSeqForward) 
      keepAvePhredVariableForward <- avePhredVariableForward >= avePhredMin
    }
    if (!is.null(L$variableSeqReverse)) {
      avePhredVariableReverse <- ShortRead::alphabetScore(
        Biostrings::quality(L$variableSeqReverse)
      )/BiocGenerics::width(L$variableSeqReverse)
      keepAvePhredVariableReverse <- avePhredVariableReverse >= avePhredMin
    }
  }
  
  ## --------------------------------------------------------------------------
  ## Count the number of Ns in variable regions
  ## --------------------------------------------------------------------------
  keepVariableNbrNForward <- keepVariableNbrNReverse <- rep(TRUE, nbrReads)
  if (!is.null(variableNMax)) {
    if (!is.null(L$variableSeqForward)) {
      variableNbrNForward <- Biostrings::vcountPattern(
        pattern = "N", subject = L$variableSeqForward,
        max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, 
        fixed = TRUE, algorithm = "auto")
      keepVariableNbrNForward <- variableNbrNForward <= variableNMax
    }
    if (!is.null(L$variableSeqReverse)) {
      variableNbrNReverse <- Biostrings::vcountPattern(
        pattern = "N", subject = L$variableSeqReverse,
        max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, 
        fixed = TRUE, algorithm = "auto")
      keepVariableNbrNReverse <- variableNbrNReverse <= variableNMax
    }
  }
  
  ## --------------------------------------------------------------------------
  ## Count the number of Ns in UMIs
  ## --------------------------------------------------------------------------
  keepUMINbrN <- rep(TRUE, nbrReads)
  if (!is.null(umiNMax)) {
    if (!is.null(L$umis)) {
      umiNbrN <- Biostrings::vcountPattern(
        pattern = "N", subject = L$umis,
        max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, 
        fixed = TRUE, algorithm = "auto")
      keepUMINbrN <- umiNbrN <= umiNMax
    }
  }
  
  ## --------------------------------------------------------------------------
  ## Find positions where variable region has mismatches
  ## --------------------------------------------------------------------------
  keepWildTypeForward <- keepWildTypeReverse <- rep(TRUE, nbrReads)
  minQualMutatedForward <- minQualMutatedReverse <- NULL
  if (!is.null(wildTypeForward)) {
    if (!is.null(L$variableSeqForward)) {
      filterForward <- mismatchFilter(wildTypeSeq = wildTypeForward, 
                                      variableSeq = L$variableSeqForward, 
                                      nbrMutatedCodonsMax = nbrMutatedCodonsMax,
                                      forbiddenMutatedCodon = forbiddenMutatedCodon)
      minQualMutatedForward <- filterForward$minQualMutated
      keepWildTypeForward <- filterForward$keepMM
    }
  }
  if (!is.null(wildTypeReverse)) {
    if (!is.null(L$variableSeqReverse)) {
      filterReverse <- mismatchFilter(wildTypeSeq = wildTypeReverse, 
                                      variableSeq = L$variableSeqReverse, 
                                      nbrMutatedCodonsMax = nbrMutatedCodonsMax,
                                      forbiddenMutatedCodon = forbiddenMutatedCodon)
      minQualMutatedReverse <- filterReverse$minQualMutated
      keepWildTypeReverse <- filterReverse$keepMM
    }
  }
  
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
  
  return(list(umis = L$umis[toKeep], 
              constantSeqForward = L$constantSeqForward[toKeep], 
              constantSeqReverse = L$constantSeqReverse[toKeep],
              variableSeqForward = L$variableSeqForward[toKeep],
              variableSeqReverse = L$variableSeqReverse[toKeep],
              minQualMutatedForward = minQualMutatedForward[toKeep],
              minQualMutatedReverse = minQualMutatedReverse[toKeep],
              readSummary = data.frame(readSummary,
                                       nbrReadPairsLowAvePhredVariableForward = sum(!keepAvePhredVariableForward),
                                       nbrReadPairsLowAvePhredVariableReverse = sum(!keepAvePhredVariableReverse),
                                       nbrReadPairsTooManyNVariableForward = sum(!keepVariableNbrNForward),
                                       nbrReadPairsTooManyNVariableReverse = sum(!keepVariableNbrNReverse),
                                       nbrReadPairsTooManyNUMI = sum(!keepUMINbrN),
                                       nbrReadPairsCodonFilterForward = ifelse(is.null(wildTypeForward),
                                                                               NA, sum(!keepWildTypeForward)),
                                       nbrReadPairsCodonFilterReverse = ifelse(is.null(wildTypeReverse),
                                                                               NA, sum(!keepWildTypeReverse)),
                                       totalNbrReadPairsPassedFilters = sum(toKeep))))

}

#' Filter reads based on mismatches with the wildtype sequence
#'
#' @param wildTypeSeq character(1), the wild type sequence
#' @param variableSeq qualityScaledDNAStringSet, containing the sequences of the
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
#' @return A list with two elements:
#'   \describe{
#'     \item{minQualMutated}{Lowest quality for the mutated bases}
#'     \item{keepMM}{Logical vector, whether to keep each read}
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
  mmpForward <- findMismatchPositions(pattern = wildTypeSeq, 
                                      subject = variableSeq)
  minQualMutated <- matrixStats::rowMins(
    as(ShortRead::FastqQuality(Biostrings::quality(variableSeq)[
      mmpForward$nucleotideMismatches]), "matrix"), na.rm = TRUE)
  
  if (!is.null(nbrMutatedCodonsMax)) {
    nbrMutCodons <- S4Vectors::elementNROWS(unique(mmpForward$codonMismatches))
    keepMM <- nbrMutCodons <= nbrMutatedCodonsMax
  }
  if (!is.null(forbiddenMutatedCodon)) {
    mutCodons <- S4Vectors::endoapply(unique(mmpForward$codonMismatches), 
                                      function(v) rep(3 * v, each = 3) + c(-2, -1, 0))
    matchForbidden <- Biostrings::vmatchPattern(
      pattern = forbiddenMutatedCodon, 
      as(variableSeq[mutCodons], "DNAStringSet"), 
      fixed = FALSE)
    keepMM <- keepMM & !vapply(matchForbidden, function(v) {
      any(IRanges::start(v) %% 3 == 1)
    }, FALSE)
  }
  list(minQualMutated = minQualMutated, keepMM = keepMM)
}