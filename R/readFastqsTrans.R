#' Read FASTQ files from TRANS experiment
#'
#' It is assumed that both the forward and reverse reads are of the form UMI -
#' constant sequence - variable sequence, and that the length of the UMI
#' sequence, and the length of the constant sequence, is the same across all the
#' reads in the same file (but it can be different for the forward and reverse 
#' reads).
#'
#' @param fastqForward,fastqReverse Character(1), FASTQ files corresponding to 
#'   forward and reverse reads, respectively.
#' @param skipForward,skipReverse Numeric(1), the number of bases to skip in the
#'   start of each read.
#' @param umiLengthForward,umiLengthReverse Numeric(1), the length of the
#'   barcode (UMI) sequence in the forward/reverse reads, respectively, not
#'   including the skipped bases (defined by \code{skipForward}/\code{skipReverse}.
#' @param constantLengthForward,constantLengthReverse Numeric(1), the lengths of
#'   the constant sequence in the forward/reverse reads, respectively.
#' @param adapterForward,adapterReverse Character(1), the adapter sequences for 
#'   forward/reverse reads, respectively. If a forward/reverse read contains the 
#'   corresponding adapter sequence, the sequence pair will be filtered out.
#'   \code{NULL} does not perform any filtering. The number of filtered read
#'   pairs are reported in the return value.
#' @param verbose Logical(1), whether to print out progress messages.
#'
#' @return A list with five elements: \describe{ \item{umis}{Merged forward and
#'   reverse UMI sequences} \item{constantSeqForward}{Constant forward sequence}
#'   \item{constantSeqReverse}{Constant reverse sequence}
#'   \item{variableSeqForward}{Variable forward sequence}
#'   \item{variableSeqReverse}{Variable reverse sequence}
#'   \item{numberReadPairsFiltered}{Number of read pairs filtered out because of 
#'   matches to adapter sequence(s)} }
#'
#' @export
#'
#' @importFrom Biostrings readQualityScaledDNAStringSet subseq
#'   QualityScaledDNAStringSet xscat quality vcountPattern 
#'   DNA_ALPHABET
#' @importFrom ShortRead readFastq
#'
#' @author Charlotte Soneson
#'   
readFastqsTrans <- function(fastqForward, fastqReverse, skipForward = 1,
                            skipReverse = 1, umiLengthForward = 10, 
                            umiLengthReverse = 8, constantLengthForward = 18,
                            constantLengthReverse = 20, adapterForward = NULL, 
                            adapterReverse = NULL, verbose = FALSE) {
  ## --------------------------------------------------------------------------
  ## Pre-flight checks
  ## --------------------------------------------------------------------------
  if (!is(fastqForward, "character") || length(fastqForward) != 1 || 
      !file.exists(fastqForward)) {
    stop("'fastqForward' must be a length-1 character vector pointing to an 
         existing file.")
  }
  if (!is(fastqReverse, "character") || length(fastqReverse) != 1 || 
      !file.exists(fastqReverse)) {
    stop("'fastqReverse' must be a length-1 character vector pointing to an 
         existing file.")
  }
  
  ## Check that adapter sequences only contain valid letters 
  ## (that would be allowed in a DNAStringSet)
  alph <- Biostrings::DNA_ALPHABET
  alph <- rev(as.character(sort(relevel(as.factor(alph), ref = "-"))))
  rgxdna <- paste0("^[", paste(alph, collapse = ""), "]+$")
  if (!is.null(adapterForward) && !grepl(rgxdna, adapterForward)) {
    stop("'adapterForward can only contain letters from Biostrings::DNA_ALPHABET")
  }
  if (!is.null(adapterReverse) && !grepl(rgxdna, adapterReverse)) {
    stop("'adapterReverse can only contain letters from Biostrings::DNA_ALPHABET")
  }
  
  ## --------------------------------------------------------------------------
  ## Read forward and reverse FASTQ files
  ## --------------------------------------------------------------------------
  if (verbose) {
    message("Reading forward FASTQ file...")
  }
  fqf <-
    as(ShortRead::readFastq(fastqForward, withIds = FALSE),
       "QualityScaledDNAStringSet")
  
  if (verbose) {
    message("Reading reverse FASTQ file...")
  }
  fqr <-
    as(ShortRead::readFastq(fastqReverse, withIds = FALSE),
       "QualityScaledDNAStringSet")
  
  ## Check that the files have the same number of reads
  stopifnot(
    length(fqf) == length(fqr)
  )
  
  ## Here we should check that the reads are in the same order, and set the 
  ## names to be the same in both files
  ## For this, we will need to set readFastq(..., withIds = TRUE)
  
  ## Search for adapter sequences and filter read pairs
  if (verbose && any(!is.null(c(adapterForward, adapterReverse)))) {
    message("Filtering out reads containing adapter sequences...")
  }
  hasAdapterMatch <- rep(FALSE, length(fqf))
  numberReadPairsFiltered <- 0L
  if (!is.null(adapterForward)) {
    hasAdapterMatch <- vcountPattern(pattern = adapterForward, subject = fqf) > 0
  }
  if (!is.null(adapterReverse)) {
    hasAdapterMatch <- hasAdapterMatch | (vcountPattern(pattern = adapterReverse, subject = fqr) > 0)
  }
  if (any(hasAdapterMatch)) {
    numberReadPairsFiltered <- sum(hasAdapterMatch)
    if (verbose) {
      message("Filtering out ", numberReadPairsFiltered, " reads (", 
              round(numberReadPairsFiltered/length(fqf) * 100, 2), "%).")
    }
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
  if (verbose) {
    message("Extracting UMI sequences...")
  }
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
  if (verbose) {
    message("Extracting constant sequences...")
  }
  constantSeqForward <- Biostrings::subseq(fqf, start = constantStartForward, 
                                           end = constantEndForward)
  constantSeqReverse <- Biostrings::subseq(fqr, start = constantStartReverse,
                                           end = constantEndReverse)
  
  ## --------------------------------------------------------------------------
  ## Extract variable sequences
  ## --------------------------------------------------------------------------
  if (verbose) {
    message("Extracting variable sequences")
  }
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
