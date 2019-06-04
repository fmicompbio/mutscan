#' Check numeric input argument
#'
#' Check that input argument is numeric, length 1 and non-negative.
#' 
#' @keywords internal
#' 
checkNumericInput <- function(...) {
  varName <- deparse(substitute(...))
  if (length(...) != 1) {
    stop(varName, " must be of length 1")
  }
  if (!is.numeric(...)) {
    stop(varName, " must be numeric")
  }
  if (... < 0) {
    stop(varName, " must be non-negative")
  }
}
  
#' Read, filter and digest sequences from two fastq files.
#'
#' Read sequences for a pair of fastq files and digest them (extract umis,
#' constant and variable parts, filter, extract mismatch information from
#' constant and count the observed unique variable parts).
#'
#' The processing of a read pair goes as follows:
#' \enumerate{
#'  \item Search for perfect matches to forward/reverse adapter sequences,
#'  filter out the read pair if a match is found in either the forward or
#'  reverse read.
#'  \item For "cis" experiment, collapse forward and reverse variable regions by
#'  retaining, for each position, the base with the highest reported base
#'  quality.
#'  \item Filter out read pair if the average quality in the variable region is
#'  below \code{avePhredMin} (for "trans" experiments, filter out a read pair if
#'  either the forward or reverse variable region falls below the quality
#'  threshold)
#'  \item Filter out read pair if the number of Ns in the variable region (for
#'  "trans" experiments, either the forward or reverse) exceeds
#'  \code{variableNMax}
#'  \item Filter out read pair if the number of Ns in the combined forward and
#'  reverse UMI sequence exceeds \code{umiNMax}
#'  \item If a wild type sequence (for the variable region) is provided, find
#'  the mismatches between the (forward/reverse) variable region and the
#'  provided wild type sequence.
#'  \item Filter out read pair if any mutated base has a quality below
#'  \code{mutatedPhredMin}
#'  \item Filter out read pair if the number of mutated codons exceeds
#'  \code{nbrMutatedCodonsMax}
#'  \item Filter out read pair if any of the mutated codons match any of the
#'  codons encoded by \code{forbiddenMutatedCodons}
#'}
#'
#' Based on the retained reads following this filtering process, count the
#' number of reads, and the number of unique UMIs, for each variable sequence
#' (or pair of variable sequences, for "trans" experiments).
#'
#' @param fastqForward,fastqReverse character(1), paths to gzipped FASTQ files
#'   corresponding to forward and reverse reads, respectively.
#' @param mergeForwardReverse logical(1), whether to fuse the forward and
#'   reverse variable sequences.
#' @param revComplForward,revComplReverse logical(1), whether to reverse
#'   complement the forward/reverse reads, respectively.
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
#' @param wildTypeForward,wildTypeReverse character(1), the wild type sequence
#'   for the forward and reverse variable region.
#' @param constantForward,constantReverse character(1), the expected constant
#'   forward and reverse sequences.
#' @param avePhredMinForward,avePhredMinReverse numeric(1) Minimum average Phred
#'   score in the variable region for a read to be retained. If L contains both
#'   forward and reverse variable regions, the minimum average Phred score has
#'   to be achieved in both for a read pair to be retained.
#' @param variableNMaxForward,variableNMaxReverse numeric(1) Maximum number of
#'   Ns allowed in the variable region for a read to be retained.
#' @param umiNMax numeric(1) Maximum number of Ns allowed in the UMI for a read
#'   to be retained.
#' @param nbrMutatedCodonsMaxForward,nbrMutatedCodonsMaxReverse numeric(1)
#'   Maximum number of mutated codons that are allowed.
#' @param forbiddenMutatedCodonsForward,forbiddenMutatedCodonsReverse character
#'   vector. Codons (can contain ambiguous IUPAC characters, see
#'   \code{\link[Biostrings]{IUPAC_CODE_MAP}}). If a read pair contains a
#'   mutated codon matching this pattern, it will be filtered out.
#' @param mutatedPhredMinForward,mutatedPhredMinReverse numeric(1) Minimum Phred
#'   score of a mutated base for the read to be retained. If any mutated base
#'   has a Phred score lower than \code{mutatedPhredMin}, the read will be
#'   discarded.
#' @param verbose logical(1), whether to print out progress messages.
#'
#' @return A list with four entries:
#' \describe{
#' \item{summaryTable}{A \code{data.frame} that contains, for each observed
#' variable region sequence (or pair of sequences, for trans experiments), a
#' simplified "mutantName", the number of observed such sequences, and the
#' number of unique UMIs obseved for the sequence. The "mutantName" for each
#' mutated codon is of the form "{f/r}xxNNN", where f/r indicates
#' forward/reverse read, xx indicates the mutated codon, and NNN is the observed
#' sequence for the codon. In the case of multiple mutated codons, these are
#' separated by underscores.}
#' \item{filterSummary}{A \code{data.frame} that contains the number of input
#' reads, the number of reads filtered out in the processing, and the number of
#' retained reads. The filters are named according to the convention
#' "fxx_filter", where "xx" indicates the order in which the filters were
#' applied, and "filter" indicates the type of filter. Note that filters are
#' applied successively, and the reads filtered out in one step are not
#' considered for successive filtering steps.}
#' \item{errorStatistics}{A \code{data.frame} that contains, for each Phred
#' quality score between 0 and 99, the number of bases in the extracted constant
#' sequences with that quality score that match/mismatch with the provided
#' reference constant sequence.}
#' \item{parameters}{A \code{list} with all parameter settings that were used in
#' the processing. Also contains the version of the package and the time of
#' processing.}
#' }
#'
#' @export
digestFastqs <- function(fastqForward, fastqReverse,
                         mergeForwardReverse, revComplForward, revComplReverse,
                         skipForward, skipReverse,
                         umiLengthForward, umiLengthReverse,
                         constantLengthForward,
                         constantLengthReverse,
                         variableLengthForward,
                         variableLengthReverse,
                         adapterForward = "", adapterReverse = "",
                         wildTypeForward = "", wildTypeReverse = "", 
                         constantForward = "", constantReverse = "", 
                         avePhredMinForward = 20.0, avePhredMinReverse = 20.0,
                         variableNMaxForward = 0, variableNMaxReverse = 0, 
                         umiNMax = 0,
                         nbrMutatedCodonsMaxForward = 1,
                         nbrMutatedCodonsMaxReverse = 1,
                         forbiddenMutatedCodonsForward = "NNW",
                         forbiddenMutatedCodonsReverse = "NNW",
                         mutatedPhredMinForward = 0.0,
                         mutatedPhredMinReverse = 0.0,
                         verbose = FALSE) {
  ## --------------------------------------------------------------------------
  ## pre-flight checks
  ## --------------------------------------------------------------------------
  ## fastq files exist
  if (length(fastqForward) != 1 || length(fastqReverse) != 1 || 
      !file.exists(fastqForward) || !file.exists(fastqReverse)) {
    stop("'fastqForward' and 'fastqReverse' must point to single, existing files");
  }
  
  ## merging/rev complementing arguments are ok
  if (any(!is.logical(c(mergeForwardReverse, revComplForward, revComplReverse))) || 
      any(c(length(mergeForwardReverse), length(revComplForward), 
            length(revComplReverse)) != 1)) {
    stop("'mergeForwardReverse', 'revComplForward' and 'revComplReverse' must be logical scalars. ")
  }
  
  ## check numeric inputs
  checkNumericInput(skipForward)
  checkNumericInput(skipReverse)
  checkNumericInput(umiLengthForward)
  checkNumericInput(umiLengthReverse)
  checkNumericInput(constantLengthForward)
  checkNumericInput(constantLengthReverse)
  checkNumericInput(variableLengthForward)
  checkNumericInput(variableLengthReverse)
  checkNumericInput(avePhredMinForward)
  checkNumericInput(avePhredMinReverse)
  checkNumericInput(variableNMaxForward)
  checkNumericInput(variableNMaxReverse)
  checkNumericInput(umiNMax)
  checkNumericInput(nbrMutatedCodonsMaxForward)
  checkNumericInput(nbrMutatedCodonsMaxReverse)
  checkNumericInput(mutatedPhredMinForward)
  checkNumericInput(mutatedPhredMinReverse)
  
  ## adapters must be strings, valid DNA characters
  if (!is.character(adapterForward) || length(adapterForward) != 1 ||
      !grepl("^[AaCcGgTt]*$", adapterForward) || !is.character(adapterReverse) ||
      length(adapterReverse) != 1 || !grepl("^[AaCcGgTt]*$", adapterReverse)) {
    stop("adapters must be character strings, only containing valid DNA characters")
  } else {
    adapterForward <- toupper(adapterForward)
    adapterReverse <- toupper(adapterReverse)
  }
  
  ## if wild type sequence is a string, make it into a vector
  if (length(wildTypeForward) == 1) {
    wildTypeForward <- c(f = wildTypeForward)
  }
  if (length(wildTypeReverse) == 1) {
    wildTypeReverse <- c(r = wildTypeReverse)
  }
  
  ## wild type sequences must be given in named vectors
  if (any(names(wildTypeForward) == "") || any(names(wildTypeReverse) == "")) {
    stop('wild type sequences must be given in named vectors')
  }
  
  ## wild type sequences must be strings, valid DNA characters
  if (!all(sapply(wildTypeForward, is.character)) || !all(sapply(wildTypeForward, length) == 1) ||
      !all(sapply(wildTypeForward, function(w) grepl("^[AaCcGgTt]*$", w))) ||
      !all(sapply(wildTypeReverse, is.character)) || !all(sapply(wildTypeReverse, length) == 1) ||
      !all(sapply(wildTypeReverse, function(w) grepl("^[AaCcGgTt]*$", w)))) {
    stop("wild type sequences must be character strings, ", 
         "only containing valid DNA characters")
  } else {
    wildTypeForward <- toupper(wildTypeForward)
    wildTypeReverse <- toupper(wildTypeReverse)
  }
  
  # if (nchar(wildTypeForward) == 0) {
  #   message("'wildTypeForward' is missing, no comparisons to ", 
  #           "wild type sequence will be done.")
  # }
  
  ## wild type sequence lengths must match variable sequence lengths
  if (any(sapply(wildTypeForward, function(w) nchar(w) > 0 && nchar(w) != variableLengthForward))) {
    stop("The lengths of the elements in 'wildTypeForward' (", paste(sapply(wildTypeForward, nchar), collapse = ","), 
         ") do not all correspond to the given 'variableLengthForward' (", 
         variableLengthForward, ")")
  }
  if (any(sapply(wildTypeReverse, function(w) nchar(w) > 0 && nchar(w) != variableLengthReverse))) {
    stop("The lengths of the elements in 'wildTypeReverse' (", paste(sapply(wildTypeReverse, nchar), collapse = ","), 
         ") do not all correspond to the given 'variableLengthReverse' (", 
         variableLengthReverse, ")")
  }
  
  ## cis experiment - should not have wildTypeReverse
  if (mergeForwardReverse && any(sapply(wildTypeReverse, nchar) > 0)) {
    warning("Ignoring 'wildTypeReverse' for CIS experiment")
    wildTypeReverse <- c(r = "")
  }
  
  ## if both constantForward and constantLengthForward are given, check that
  ## the lengths inferred from the two are consistent
  if (!is.character(constantForward) || !is.character(constantReverse) ||
      length(constantForward) != 1 || length(constantReverse) != 1 ||
      !grepl("^[AaCcGgTt]*$", constantForward) || 
      !grepl("^[AaCcGgTt]*$", constantReverse)) {
    stop("constant sequences must be character strings, ", 
         "only containing valid DNA characters.")
  } else {
    constantForward <- toupper(constantForward)
    constantReverse <- toupper(constantReverse)
  }
  
  if (nchar(constantForward) > 0 && nchar(constantForward) != constantLengthForward) {
    stop("'constantLengthForward' (", constantLengthForward, 
         ") does not correspond to the length of the given 'constantForward' (", 
         nchar(constantForward), ")")
  }
  if (nchar(constantReverse) > 0 && nchar(constantReverse) != constantLengthReverse) {
    stop("'constantLengthReverse' (", constantLengthReverse, 
         ") does not correspond to the length of the given 'constantReverse' (", 
         nchar(constantReverse), ")")
  }
  
  if (!all(is.character(forbiddenMutatedCodonsForward)) || 
      !all(grepl("^[ACGTMRWSYKVHDBN]{3}$", toupper(forbiddenMutatedCodonsForward)) | 
           forbiddenMutatedCodonsForward == "") ||
      !all(is.character(forbiddenMutatedCodonsReverse)) || 
      !all(grepl("^[ACGTMRWSYKVHDBN]{3}$", toupper(forbiddenMutatedCodonsReverse)) | 
           forbiddenMutatedCodonsReverse == "")) {
    stop("All elements of 'forbiddenMutatedCodonsForward' and 'forbiddenMutatedCodonsReverse' must be ", 
         "character strings consisting of three valid IUPAC letters.")
  } else {
    forbiddenMutatedCodonsForward <- toupper(forbiddenMutatedCodonsForward)
    forbiddenMutatedCodonsReverse <- toupper(forbiddenMutatedCodonsReverse)
  }
  
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be a logical scalar.")
  }
  
  res <- digestFastqsCpp(fastqForward = fastqForward, 
                         fastqReverse = fastqReverse,
                         mergeForwardReverse = mergeForwardReverse,
                         revComplForward = revComplForward,
                         revComplReverse = revComplReverse,
                         skipForward = skipForward, 
                         skipReverse = skipReverse,
                         umiLengthForward = umiLengthForward, 
                         umiLengthReverse = umiLengthReverse,
                         constantLengthForward = constantLengthForward,
                         constantLengthReverse = constantLengthReverse,
                         variableLengthForward = variableLengthForward,
                         variableLengthReverse = variableLengthReverse,
                         adapterForward = adapterForward, 
                         adapterReverse = adapterReverse,
                         wildTypeForward = wildTypeForward, 
                         wildTypeReverse = wildTypeReverse, 
                         constantForward = constantForward, 
                         constantReverse = constantReverse, 
                         avePhredMinForward = avePhredMinForward,
                         avePhredMinReverse = avePhredMinReverse,
                         variableNMaxForward = variableNMaxForward,
                         variableNMaxReverse = variableNMaxReverse,
                         umiNMax = umiNMax,
                         nbrMutatedCodonsMaxForward = nbrMutatedCodonsMaxForward,
                         nbrMutatedCodonsMaxReverse = nbrMutatedCodonsMaxReverse,
                         forbiddenMutatedCodonsForward = forbiddenMutatedCodonsForward,
                         forbiddenMutatedCodonsReverse = forbiddenMutatedCodonsReverse,
                         mutatedPhredMinForward = mutatedPhredMinForward,
                         mutatedPhredMinReverse = mutatedPhredMinReverse,
                         verbose = verbose)
  
  ## Add package version and processing date
  res$parameters$processingInfo <- paste0(
    "Processed by mutscan v", utils::packageVersion("mutscan"), " on ",
    Sys.time()
  )
  res
}
