#' Filter FASTQ files from TRANS experiment
#' 
#' @param L List, as output by \code{readTransFastqs}.
#' @param avePhredMin Numeric(1) Minimum average Phred score in the 
#'   variable region for a read to be retained.
#' @param variableNMax Numeric(1) Maximum number of Ns in the 
#'   variable region for a read to be retained.
#' @param umiNMax Numeric(1) Maximum number of Ns in the UMI for a 
#'   read to be retained.   
#'
#' @return A filtered list with eight elements:
#' \describe{
#'   \item{umis}{Merged forward and reverse UMI sequences}
#'   \item{constantSeqForward}{Constant forward sequence}
#'   \item{constantSeqReverse}{Constant reverse sequence}
#'   \item{variableSeqForward}{Variable forward sequence}
#'   \item{variableSeqReverse}{Variable reverse sequence}
#'   \item{minQualMutatedForward}{Minimal base quality of a mutated base in forward sequence}
#'   \item{minQualMutatedReverse}{Minimal base quality of a mutated base in reverse sequence}
#'   \item{readSummary}{data.frame tabulating the total number of read pairs, 
#'   and the number that match any of the filtering criteria.} 
#' }
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
filterTrans <- function(L, avePhredMin = 20, variableNMax = 0, umiNMax = 0,
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
    avePhredVariableForward <- ShortRead::alphabetScore(
      Biostrings::quality(L$variableSeqForward)
    )/BiocGenerics::width(L$variableSeqForward) 
    keepAvePhredVariableForward <- avePhredVariableForward >= avePhredMin
    
    avePhredVariableReverse <- ShortRead::alphabetScore(
      Biostrings::quality(L$variableSeqReverse)
    )/BiocGenerics::width(L$variableSeqReverse)
    keepAvePhredVariableReverse <- avePhredVariableReverse >= avePhredMin
  }
  
  ## --------------------------------------------------------------------------
  ## Count the number of Ns in variable regions
  ## --------------------------------------------------------------------------
  keepVariableNbrNForward <- keepVariableNbrNReverse <- rep(TRUE, nbrReads)
  if (!is.null(variableNMax)) {
    variableNbrNForward <- Biostrings::vcountPattern(
      pattern = "N", subject = L$variableSeqForward,
      max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, 
      fixed = TRUE, algorithm = "auto")
    keepVariableNbrNForward <- variableNbrNForward <= variableNMax
    
    variableNbrNReverse <- Biostrings::vcountPattern(
      pattern = "N", subject = L$variableSeqReverse,
      max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, 
      fixed = TRUE, algorithm = "auto")
    keepVariableNbrNReverse <- variableNbrNReverse <= variableNMax
  }
  
  ## --------------------------------------------------------------------------
  ## Count the number of Ns in UMIs
  ## --------------------------------------------------------------------------
  keepUMINbrN <- rep(TRUE, nbrReads)
  if (!is.null(umiNMax)) {
    umiNbrN <- Biostrings::vcountPattern(
      pattern = "N", subject = L$umis,
      max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, 
      fixed = TRUE, algorithm = "auto")
    keepUMINbrN <- umiNbrN <= umiNMax
  }
  
  ## --------------------------------------------------------------------------
  ## Find positions where variable region has mismatches
  ## --------------------------------------------------------------------------
  keepWildTypeForward <- keepWildTypeReverse <- rep(TRUE, nbrReads)
  if (!is.null(wildTypeForward)) {
    mmpForward <- findMismatchPositions(pattern = wildTypeForward, 
                                        subject = L$variableSeqForward)
    minQualMutatedForward <- rowMins(as(FastqQuality(quality(L$variableSeqForward)[
      mmpForward$nucleotideMismatches]), "matrix"), na.rm = TRUE)
    
    if (!is.null(nbrMutatedCodonsMax)) {
      nbrMutCodonsForward <- elementNROWS(unique(mmpForward$codonMismatches))
      keepWildTypeForward <- keepWildTypeForward & nbrMutCodonsForward <= nbrMutatedCodonsMax
    }
    if (!is.null(forbiddenMutatedCodon)) {
      mutCodons <- endoapply(unique(mmpForward$codonMismatches), 
                             function(v) rep(3 * v, each = 3) + c(-2, -1, 0))
      matchForbidden <- vmatchPattern(pattern = forbiddenMutatedCodon, 
                                      as(L$variableSeqForward[mutCodons], "DNAStringSet"), 
                                      fixed = FALSE)
      keepWildTypeForward <- keepWildTypeForward & !vapply(matchForbidden, function(v) {
        any(IRanges::start(v) %% 3 == 1)
      }, FALSE)
    }
  }
  if (!is.null(wildTypeReverse)) {
    mmpReverse <- findMismatchPositions(pattern = wildTypeReverse, 
                                        subject = L$variableSeqReverse)
    minQualMutatedReverse <- rowMins(as(FastqQuality(quality(L$variableSeqReverse)[
      mmpReverse$nucleotideMismatches]), "matrix"), na.rm = TRUE)
    
    if (!is.null(nbrMutatedCodonsMax)) {
      nbrMutCodonsReverse <- elementNROWS(unique(mmpReverse$codonMismatches))
      keepWildTypeReverse <- keepWildTypeReverse & nbrMutCodonsReverse <= nbrMutatedCodonsMax
    }
    if (!is.null(forbiddenMutatedCodon)) {
      mutCodons <- endoapply(unique(mmpReverse$codonMismatches), 
                             function(v) rep(3 * v, each = 3) + c(-2, -1, 0))
      matchForbidden <- vmatchPattern(pattern = forbiddenMutatedCodon, 
                                      as(L$variableSeqReverse[mutCodons], "DNAStringSet"), 
                                      fixed = FALSE)
      keepWildTypeReverse <- keepWildTypeReverse & !vapply(matchForbidden, function(v) {
        any(IRanges::start(v) %% 3 == 1)
      }, FALSE)
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
