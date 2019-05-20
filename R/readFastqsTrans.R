#' Read FASTQ files from TRANS experiment
#'
#' It is assumed that the forward and reverse reads are of the form UMI -
#' constant sequence - variable sequence, and that the length of the UMI
#' sequence, and the length of the constant sequence, is the same across all the
#' reads.
#'
#' @param fastqForward,fastqReverse FASTQ files corresponding to forward and
#'   reverse reads, respectively.
#' @param skipForward,skipReverse Numeric(1), the number of bases to skip in the
#'   start of each read.
#' @param umiLengthForward,umiLengthReverse Numeric(1), the length of the
#'   barcode (UMI) sequence in the forward/reverse reads, respectively, not
#'   including the skipped bases.
#' @param constantLengthForward,constantLengthReverse Numeric(1), the length of
#'   the constant sequence in the forward/reverse reads, respectively.
#' @param adapterForward,adapterReverse Character(1), the adapter sequences, if
#'   found in forward/reverse sequences, the sequence pair will be filtered out.
#'   \code{NULL} does not perform any filtering. The number of filtered read
#'   pairs are reported in the return value.
#'
#' @return A list with five elements: \describe{ \item{umis}{Merged forward and
#'   reverse UMI sequences} \item{constantSeqForward}{Constant forward sequence}
#'   \item{constantSeqReverse}{Constant reverse sequence}
#'   \item{variableSeqForward}{Variable forward sequence}
#'   \item{variableSeqReverse}{Variable reverse sequence}
#'   \item{numberReadPairsFiltered}{because of matches to adaptor sequence(s)} }
#'
#' @export
#'
#' @importFrom Biostrings readQualityScaledDNAStringSet subseq
#'   QualityScaledDNAStringSet xscat quality vcountPattern
#' @importFrom ShortRead readFastq
#'
#' @author Charlotte Soneson
#'   
readFastqsTrans <- function(fastqForward, fastqReverse, skipForward = 1,
                            skipReverse = 1, umiLengthForward = 10, 
                            umiLengthReverse = 8, constantLengthForward = 18,
                            constantLengthReverse = 20,
                            adapterForward = NULL, adapterReverse = NULL) {
  ## --------------------------------------------------------------------------
  ## Read forward and reverse FASTQ files
  ## --------------------------------------------------------------------------
  fqf <- Biostrings::readQualityScaledDNAStringSet(
    fastqForward, quality.scoring = "phred",
    nrec = -1L, skip = 0L, seek.first.rec = FALSE,
    use.names = TRUE)
  
  fqr <- Biostrings::readQualityScaledDNAStringSet(
    fastqReverse, quality.scoring = "phred",
    nrec = -1L, skip = 0L, seek.first.rec = FALSE,
    use.names = TRUE)
  
  ## Check that the files have the same number of reads
  stopifnot(
    length(fqf) == length(fqr)
  )
  
  ## Here we should check that the reads are in the same order, and set the 
  ## names to be the same in both files
  
  ## Search for adapter sequences and filter read pairs
  hasAdapterMatch <- rep(FALSE, length(fqf))
  numberReadPairsFiltered <- 0L
  if (!is.null(adapterForward)) {
    hasAdapterMatch <- vcountPattern(pattern = adapterForward, subject = fqf) > 0
  }
  if (!is.null(adapterReverse)) {
    hasAdapterMatch <- hasAdapterMatch & (vcountPattern(pattern = adapterReverse, subject = fqr) > 0)
  }
  if (any(hasAdapterMatch)) {
    numberReadPairsFiltered <- sum(hasAdapterMatch)
    fqf <- fqf[!hasAdapterMatch]
    fqr <- fqr[!hasAdapterMatch]
  }
  
  ## --------------------------------------------------------------------------
  ## Define start and end positions of the different parts of the reads
  ## --------------------------------------------------------------------------
  umiStartForward <- 1 + skipForward
  umiStartReverse <- 1 + skipReverse
  umiEndForward <- umiStartForward + umiLengthForward - 1
  umiEndReverse <- umiStartReverse + umiLengthReverse - 1
  constantStartForward <- umiEndForward + 1
  constantStartReverse <- umiEndReverse + 1
  constantEndForward <- constantStartForward + constantLengthForward - 1
  constantEndReverse <- constantStartReverse + constantLengthReverse - 1
  variableStartForward <- constantEndForward + 1
  variableStartReverse <- constantEndReverse + 1
  
  ## --------------------------------------------------------------------------
  ## Extract and concatenate UMIs
  ## --------------------------------------------------------------------------
  ## It would be nice if this could be done more cleverly - xscat does not 
  ## return a QualityScaledDNAStringSet if provided with two such object
  umif <- Biostrings::subseq(fqf, start = umiStartForward, 
                             end = umiEndForward)
  umir <- Biostrings::subseq(fqr, start = umiStartReverse, 
                             end = umiEndReverse)
  umis <- Biostrings::QualityScaledDNAStringSet(
    x = Biostrings::xscat(umif, umir),
    quality = as(Biostrings::xscat(Biostrings::quality(umif), 
                                   Biostrings::quality(umir)),
                 "PhredQuality")
  )
  names(umis) <- names(umif)

  ## --------------------------------------------------------------------------
  ## Extract constant sequences
  ## --------------------------------------------------------------------------
  constantSeqForward <- Biostrings::subseq(fqf, start = constantStartForward, 
                                           end = constantEndForward)
  constantSeqReverse <- Biostrings::subseq(fqr, start = constantStartReverse,
                                           end = constantEndReverse)
  
  ## --------------------------------------------------------------------------
  ## Extract variable sequences
  ## --------------------------------------------------------------------------
  variableSeqForward <- Biostrings::subseq(fqf, start = variableStartForward)
  variableSeqReverse <- Biostrings::subseq(fqr, start = variableStartReverse)
  
  ## --------------------------------------------------------------------------
  ## Return values
  ## --------------------------------------------------------------------------
  ## For now we return a list - consider defining a new class if 
  ## that would make things easier. There is also the DNAStringSetList class,
  ## but that would convert the objects to regular DNAStringSet objects
  return(list(
    umis = umis,
    constantSeqForward = constantSeqForward, 
    constantSeqReverse = constantSeqReverse,
    variableSeqForward = variableSeqForward,
    variableSeqReverse = variableSeqReverse,
    numberReadPairsFiltered = numberReadPairsFiltered
  ))
  
}
