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
#' @return A filtered list with five element:
#' \describe{
#'   \item{umis}{Merged forward and reverse UMI sequences}
#'   \item{constantSeqForward}{Constant forward sequence}
#'   \item{constantSeqReverse}{Constant reverse sequence}
#'   \item{variableSeqForward}{Variable forward sequence}
#'   \item{variableSeqReverse}{Variable reverse sequence}
#' }
#' 
#' @export
#' 
#' @importFrom Biostrings quality vcountPattern
#' @importFrom ShortRead alphabetScore
#' @importFrom BiocGenerics width
#' 
#' @author Charlotte Soneson
#'
filterTrans <- function(L, avePhredMin = 20, 
                        variableNMax = 0, umiNMax = 0) {
  ## --------------------------------------------------------------------------
  ## Calculate average Phred score
  ## --------------------------------------------------------------------------
  avePhredVariableForward <- ShortRead::alphabetScore(
    Biostrings::quality(L$variableSeqForward)
  )/BiocGenerics::width(L$variableSeqForward) 
  avePhredVariableReverse <- ShortRead::alphabetScore(
    Biostrings::quality(L$variableSeqReverse)
  )/BiocGenerics::width(L$variableSeqReverse)
  
  ## --------------------------------------------------------------------------
  ## Check if variable region contains an N
  ## --------------------------------------------------------------------------
  variableNbrNForward <- Biostrings::vcountPattern(
    pattern = "N", subject = L$variableSeqForward,
    max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, 
    fixed = TRUE, algorithm = "auto")
  variableNbrNReverse <- Biostrings::vcountPattern(
    pattern = "N", subject = L$variableSeqReverse,
    max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, 
    fixed = TRUE, algorithm = "auto")
  
  ## --------------------------------------------------------------------------
  ## Check if UMI contains an N
  ## --------------------------------------------------------------------------
  umiNbrN <- Biostrings::vcountPattern(
    pattern = "N", subject = L$umis,
    max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, 
    fixed = TRUE, algorithm = "auto")
  
  ## --------------------------------------------------------------------------
  ## Return filtered object
  ## --------------------------------------------------------------------------
  toKeep <- avePhredVariableForward >= avePhredMin & 
    avePhredVariableReverse >= avePhredMin & 
    variableNbrNForward <= variableNMax & 
    variableNbrNReverse <= variableNMax & 
    umiNbrN <= umiNMax
  
  return(list(umis = L$umis[toKeep], 
              constantSeqForward = L$constantSeqForward[toKeep], 
              constantSeqReverse = L$constantSeqReverse[toKeep],
              variableSeqForward = L$variableSeqForward[toKeep],
              variableSeqReverse = L$variableSeqReverse[toKeep]))

}
