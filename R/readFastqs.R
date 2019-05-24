#' Read FASTQ files from CIS or TRANS experiment
#'
#' It is assumed that both the forward and reverse reads are of the form UMI -
#' constant sequence - variable sequence, and that the length of the UMI
#' sequence, and the length of the constant sequence, is the same across all the
#' reads in the same file (but it can be different for the forward and reverse
#' reads).
#'
#' @param experimentType character(1), either "cis" or "trans". If this is set
#'   to "cis", the variable sequences from the forward and reverse reads will be
#'   consolidated into one single sequence.
#' @param fastqForward,fastqReverse character(1), paths to FASTQ files
#'   corresponding to forward and reverse reads, respectively.
#' @param skipForward,skipReverse numeric(1), the number of bases to skip in the
#'   start of each forward and reverse read, respectively.
#' @param umiLengthForward,umiLengthReverse numeric(1), the length of the
#'   barcode (UMI) sequence in the forward/reverse reads, respectively, not
#'   including the skipped bases (defined by
#'   \code{skipForward}/\code{skipReverse}).
#' @param constantLengthForward,constantLengthReverse numeric(1), the length of
#'   the constant sequence in the forward/reverse reads, respectively.
#' @param variableLengthForward,variableLengthReverse numeric(1), the length of
#'   the variable sequence in the forward/reverse reads, respectively.
#' @param adapterForward,adapterReverse character(1), the adapter sequence for
#'   forward/reverse reads, respectively. If a forward/reverse read contains the
#'   corresponding adapter sequence, the sequence pair will be filtered out.
#'   If set to \code{NULL}, no adapter filtering is performed. The number of
#'   filtered read pairs are reported in the return value.
#' @param verbose logical(1), whether to print out progress messages.
#'
#' @return A SummarizedExperiment object list with four or five assays:
#'   \describe{ 
#'   \item{umis}{Merged forward and reverse UMI sequences}
#'   \item{constantSeqForward}{Constant forward sequence}
#'   \item{constantSeqReverse}{Constant reverse sequence}
#'   \item{variableSeqForward}{Variable forward sequence}
#'   \item{variableSeqReverse}{Variable reverse sequence, only for TRANS
#'   experiments} 
#'   } 
#'   Each assay is represented as a \code[S4Vectors]{DataFrame} with one column
#'   named \code{seq}, which contains a
#'   \code[Biostrings]{QualityScaledDNAStringSet} object.
#'
#' @export
#'
#' @importFrom Biostrings readQualityScaledDNAStringSet subseq
#'   QualityScaledDNAStringSet xscat quality vcountPattern DNA_ALPHABET
#'   reverseComplement
#' @importFrom ShortRead readFastq
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame
#' @importFrom stats relevel
#'
#' @author Charlotte Soneson
#' 
#' @examples 
#' datadir <- system.file("extdata", package = "mutscan")
#' transInput <- readFastqs(experimentType = "trans", 
#'                          fastqForward = file.path(datadir, "transInput_1.fastq.gz"),
#'                          fastqReverse = file.path(datadir, "transInput_2.fastq.gz"),
#'                          skipForward = 1, skipReverse = 1, umiLengthForward = 10, 
#'                          umiLengthReverse = 8, constantLengthForward = 18,
#'                          constantLengthReverse = 20, adapterForward = "GGAAGAGCACACGTC", 
#'                          adapterReverse = "GGAAGAGCGTCGTGT", verbose = FALSE)
#'
readFastqs <- function(experimentType, fastqForward, fastqReverse, skipForward = 1,
                       skipReverse = 1, umiLengthForward = 10, 
                       umiLengthReverse = 8, constantLengthForward = 18,
                       constantLengthReverse = 20, variableLengthForward = 96,
                       variableLengthReverse = 96, adapterForward = NULL, 
                       adapterReverse = NULL, verbose = FALSE) {
  ## --------------------------------------------------------------------------
  ## Pre-flight checks
  ## --------------------------------------------------------------------------
  readFastqsInputCheck(experimentType = experimentType, fastqForward = fastqForward, 
                       fastqReverse = fastqReverse, skipForward = skipForward,
                       skipReverse = skipReverse, umiLengthForward = umiLengthForward, 
                       umiLengthReverse = umiLengthReverse,
                       constantLengthForward = constantLengthForward,
                       constantLengthReverse = constantLengthReverse,
                       variableLengthForward = variableLengthForward,
                       variableLengthReverse = variableLengthReverse, 
                       adapterForward = adapterForward, 
                       adapterReverse = adapterReverse, verbose = verbose)
  
  ## --------------------------------------------------------------------------
  ## Read forward and reverse FASTQ files
  ## --------------------------------------------------------------------------
  if (verbose) {
    message("Reading forward FASTQ file...")
  }
  fqf <- as(ShortRead::readFastq(fastqForward, withIds = FALSE), 
            "QualityScaledDNAStringSet")
  
  if (verbose) {
    message("Reading reverse FASTQ file...")
  }
  fqr <- as(ShortRead::readFastq(fastqReverse, withIds = FALSE),
            "QualityScaledDNAStringSet")
  
  ## Check that the files have the same number of reads
  stopifnot(
    length(fqf) == length(fqr)
  )
  totalNbrReads <- length(fqf)
  
  ## Here we should check that the reads are in the same order, and set the 
  ## names to be the same in both files
  ## For this, we will need to set readFastq(..., withIds = TRUE)
  
  ## --------------------------------------------------------------------------
  ## Search for adapter sequences and filter read pairs
  ## --------------------------------------------------------------------------
  if (verbose && any(!is.null(c(adapterForward, adapterReverse)))) {
    message("Filtering out reads containing adapter sequences...")
  }
  hasAdapterMatch <- rep(FALSE, totalNbrReads)
  numberReadPairsFiltered <- 0L
  if (!is.null(adapterForward)) {
    hasAdapterMatch <- Biostrings::vcountPattern(pattern = adapterForward, 
                                                 subject = fqf) > 0
  }
  if (!is.null(adapterReverse)) {
    hasAdapterMatch <- hasAdapterMatch | 
      (Biostrings::vcountPattern(pattern = adapterReverse, subject = fqr) > 0)
  }
  if (any(hasAdapterMatch)) {
    numberReadPairsFiltered <- sum(hasAdapterMatch)
    if (verbose) {
      message("   Filtered out ", numberReadPairsFiltered, " reads (", 
              round(numberReadPairsFiltered/totalNbrReads * 100, 2), "%).")
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
  variableEndForward <- variableStartForward + variableLengthForward - 1
  variableEndReverse <- variableStartReverse + variableLengthReverse - 1
  
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
  if (experimentType == "cis") {
    umir <- Biostrings::reverseComplement(umir)
  }
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
  if (experimentType == "cis") {
    constantSeqReverse <- Biostrings::reverseComplement(constantSeqReverse)
  }
  
  ## --------------------------------------------------------------------------
  ## Extract variable sequences
  ## --------------------------------------------------------------------------
  if (verbose) {
    message("Extracting variable sequences...")
  }
  variableSeqForward <- Biostrings::subseq(fqf, start = variableStartForward,
                                           end = variableEndForward)
  variableSeqReverse <- Biostrings::subseq(fqr, start = variableStartReverse,
                                           end = variableEndReverse)
  if (experimentType == "cis") {
    variableSeqForward <- mergeReadPairs(
      readsForward = variableSeqForward,
      readsReverse = Biostrings::reverseComplement(variableSeqReverse)
    )
    variableSeqReverse <- NULL
  }
  
  ## --------------------------------------------------------------------------
  ## Return SummarizedExperiment object
  ## --------------------------------------------------------------------------
  se <- 
    SummarizedExperiment::SummarizedExperiment(
      assays = list(umis = S4Vectors::DataFrame(seq = umis),
                    constantSeqForward = S4Vectors::DataFrame(seq = constantSeqForward),
                    constantSeqReverse = S4Vectors::DataFrame(seq = constantSeqReverse),
                    variableSeqForward = S4Vectors::DataFrame(seq = variableSeqForward)
      ),
      colData = S4Vectors::DataFrame(
        totalNbrReadPairs = totalNbrReads, 
        nbrReadPairsWithAdapter = numberReadPairsFiltered,
        experimentType = experimentType
      )
    )
  if (!is.null(variableSeqReverse)) {
    SummarizedExperiment::assay(se, "variableSeqReverse") <- 
      S4Vectors::DataFrame(seq = variableSeqReverse)
  }
  
  return(se)
}

#' Check inputs for readFastqs
#' 
#' @inheritParams readFastqs
#' 
#' @keywords internal
#' 
readFastqsInputCheck <- function(experimentType, fastqForward, fastqReverse, skipForward,
                                 skipReverse, umiLengthForward, 
                                 umiLengthReverse, constantLengthForward,
                                 constantLengthReverse, variableLengthForward,
                                 variableLengthReverse, adapterForward, 
                                 adapterReverse, verbose) {
  if (!is(experimentType, "character") || length(experimentType) != 1 || 
      !(experimentType %in% c("cis", "trans"))) {
    stop("'experimentType' must be either 'cis' or 'trans'")
  }
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
  if (!all(is(c(skipForward, skipReverse, umiLengthForward, umiLengthReverse,
                constantLengthForward, constantLengthReverse,
                variableLengthForward, variableLengthReverse), "numeric"))) {
    stop("'skipForward', 'skipReverse', 'umiLengthForward', 'umiLengthReverse',
    'constantLengthForward', 'constantLengthReverse', 'variableLengthForward' and 
    'variableLengthReverse' must be numeric.")
  }
  
  ## Check that adapter sequences only contain valid letters 
  ## (that would be allowed in a DNAStringSet)
  alph <- Biostrings::DNA_ALPHABET
  alph <- rev(as.character(sort(stats::relevel(as.factor(alph), ref = "-"))))
  rgxdna <- paste0("^[", paste(alph, collapse = ""), "]+$")
  if (!is.null(adapterForward) && !grepl(rgxdna, adapterForward)) {
    stop("'adapterForward can only contain letters from Biostrings::DNA_ALPHABET")
  }
  if (!is.null(adapterReverse) && !grepl(rgxdna, adapterReverse)) {
    stop("'adapterReverse can only contain letters from Biostrings::DNA_ALPHABET")
  }
}