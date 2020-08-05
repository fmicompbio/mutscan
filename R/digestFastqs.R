#' Check numeric input argument
#'
#' Check that input argument is numeric, length 1 and non-negative (or -1, if
#' negative values are allowed).
#' 
#' @keywords internal
#' 
checkNumericInput <- function(..., nonnegative) {
  varName <- deparse(substitute(...))
  if (length(...) != 1) {
    stop(varName, " must be of length 1")
  }
  if (!is.numeric(...)) {
    stop(varName, " must be numeric")
  }
  if (nonnegative && ... < 0) {
    stop(varName, " must be non-negative")
  }
  if (!nonnegative && (... < 0 && ... != -1)) {
    stop(varName, " must be non-negative or -1")
  }
}
  
#' Read, filter and digest sequences from fastq file(s).
#'
#' Read sequences for one or a pair of fastq files and digest them (extract
#' umis, constant and variable parts, filter, extract mismatch information from
#' constant and count the observed unique variable parts). Alternatively, primer
#' sequences could be specified, in which case the sequence immediately
#' following the primer will be considered the variable sequence.
#'
#' The processing of a read pair goes as follows:
#' \enumerate{
#'  \item Search for perfect matches to forward/reverse adapter sequences,
#'  filter out the read pair if a match is found in either the forward or
#'  reverse read.
#'  \item If primer sequences are provided, search for perfect matches, and
#'  filter out the read pair if not all provided primer sequences can be found.
#'  \item Extract the variable sequence from forward and reverse reads. If the
#'  lengths of the UMI, constant and variable part of the read are given, they
#'  will be used to extract the variable part. Otherwise, a primer sequence must
#'  be provided, and the variable sequence is assumed to start immediately after
#'  the primer.
#'  \item If requested, collapse forward and reverse variable regions by
#'  retaining, for each position, the base with the highest reported base
#'  quality.
#'  \item Filter out read pair if the average quality in the variable region is
#'  below \code{avePhredMinForward}/\code{avePhredMinReverse}, in either the
#'  forward or reverse read (or the merged read).
#'  \item Filter out read pair if the number of Ns in the variable region
#'  exceeds \code{variableNMaxForward}/\code{variableNMaxReverse}.
#'  \item Filter out read pair if the number of Ns in the combined forward and
#'  reverse UMI sequence exceeds \code{umiNMax}
#'  \item If one or more wild type sequences (for the variable region) are
#'  provided, find the mismatches between the (forward/reverse) variable region
#'  and the provided wild type sequence (if more than one wild type sequence is
#'  provided, first find the one that is closest to the read).
#'  \item Filter out read pair if any mutated base has a quality below
#'  \code{mutatedPhredMinForward}/\code{mutatedPhredMinReverse}.
#'  \item Filter out read pair if the number of mutated codons exceeds
#'  \code{nbrMutatedCodonsMaxForward}/\code{nbrMutatedCodonsMaxReverse}.
#'  \item Filter out read pair if any of the mutated codons match any of the
#'  codons encoded by
#'  \code{forbiddenMutatedCodonsForward}/\code{forbiddenMutatedCodonsReverse}.
#'  \item Assign a 'mutation name' to the read. This name is a combination of
#'  parts of the form XX{.}YY{.}NNN, where XX is the name of the most similar
#'  reference sequence, YY is the mutated codon number, and NNN is the mutated
#'  codon. {.} is a delimiter, specified via \code{mutNameDelimiter}.
#'}
#' 
#' Based on the retained reads following this filtering process, count the
#' number of reads, and the number of unique UMIs, for each variable sequence
#' (or pair of variable sequences).
#'
#' @param fastqForward,fastqReverse character vector, paths to gzipped FASTQ files
#'   corresponding to forward and reverse reads, respectively. If more than one
#'   forward/reverse sequence file is given, they need to be provided in the
#'   same order.
#' @param mergeForwardReverse logical(1), whether to fuse the forward and
#'   reverse variable sequences.
#' @param minOverlap,maxOverlap numeric(1), the minimal and maximal allowed
#'   overlap between the forward and reverse reads when merging. Only used if
#'   \code{mergeForwardReverse} is \code{TRUE}. If set to 0, only overlaps
#'   covering the full length of the shortest of the two reads will be
#'   considered.
#' @param maxFracMismatchOverlap numeric(1), maximal mismatch rate in the
#'   overlap. Only used if \code{mergeForwardReverse} is \code{TRUE}.
#' @param greedyOverlap logical(1). If \code{TRUE}, the first overlap satisfying
#'   \code{minOverlap}, \code{maxOverlap} and \code{maxFracMismatchOverlap} will
#'   be retained. If \code{FALSE}, all valid overlaps will be scored and the one
#'   with the highest score (largest number of matches) will be retained.
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
#' @param primerForward,primerReverse character(1), the primer sequence for
#'   forward/reverse reads, respectively. Only read pairs that contain both the
#'   forward and reverse primers will be retained.
#' @param wildTypeForward,wildTypeReverse character(1) or named character
#'   vector, the wild type sequence for the forward and reverse variable region.
#'   If given as a single string, the reference sequence will be named 'f' (for
#'   forward) or 'r' (for reverse).
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
#' @param mutNameDelimiter character(1) Delimiter used in the naming of mutants.
#'   Generally, mutants will be named as XX{.}YY{.}NNN, where XX is the closest
#'   provided reference sequence, YY is the mutated codon number, and NNN is the
#'   mutated codon. Here, {.} is the provided \code{mutNameDelimiter}. The
#'   delimiter must be a single character (not "_"), and can not appear in any
#'   of the provided reference sequence names.
#' @param variableCollapseMaxDist,umiCollapseMaxDist A \code{numeric} scalar
#'   defining the tolerances for collapsing similar variable (forward + reverse,
#'   if any) or UMI sequences. If the value is in [0, 1), it defines the maximal
#'   Hamming distance in terms of a fraction of sequence length:
#'   (\code{round(variableCollapseMaxDist * nchar(variableSeq))}).
#'   A value greater or equal to 1 is rounded and directly used as the maximum
#'   allowed Hamming distance. Note that variable sequences can only be
#'   collapsed if they are all of the same length and no wild type sequences
#'   (\code{wildTypeForward} or \code{wildTypeReverse}) have been given.
#' @param maxNReads integer(1) Maximal number of reads to process. If set to -1,
#'   all reads will be processed.
#' @param verbose logical(1), whether to print out progress messages.
#'
#' @return A list with four entries:
#' \describe{
#' \item{summaryTable}{A \code{data.frame} that contains, for each observed
#' mutation combination, the corresponding variable region sequences (or pair of
#' sequences), the number of observed such sequences, and the number of unique
#' UMIs obseved for the sequence.}
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
digestFastqs <- function(fastqForward, fastqReverse = NULL,
                         mergeForwardReverse = FALSE, minOverlap = 0, maxOverlap = 0, 
                         maxFracMismatchOverlap = 1, greedyOverlap = TRUE,
                         revComplForward = FALSE, revComplReverse = FALSE,
                         skipForward = 0, skipReverse = -1,
                         umiLengthForward = 0, umiLengthReverse = -1,
                         constantLengthForward = 0,
                         constantLengthReverse = -1,
                         variableLengthForward,
                         variableLengthReverse = -1,
                         adapterForward = "", adapterReverse = "",
                         primerForward = "", primerReverse = "",
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
                         mutNameDelimiter = ".",
                         variableCollapseMaxDist = 0.0,
                         umiCollapseMaxDist = 0.0,
                         maxNReads = -1, verbose = FALSE) {
  ## --------------------------------------------------------------------------
  ## pre-flight checks
  ## --------------------------------------------------------------------------
  ## fastq files exist
  if (length(fastqForward) < 1 || !all(file.exists(fastqForward)) || 
      (!is.null(fastqReverse) && (length(fastqReverse) != length(fastqForward) ||
                                 !all(file.exists(fastqReverse))))) {
    stop("'fastqForward' and 'fastqReverse' must point to one or several matching existing files");
  }
  if (is.null(fastqReverse)) {
    fastqReverse <- rep("", length(fastqForward))
  }
  
  ## merging/rev complementing arguments are ok
  if (any(!is.logical(c(mergeForwardReverse, revComplForward, revComplReverse))) || 
      any(c(length(mergeForwardReverse), length(revComplForward), 
            length(revComplReverse)) != 1)) {
    stop("'mergeForwardReverse', 'revComplForward' and 'revComplReverse' must be logical scalars")
  }
  if (fastqReverse == "" && mergeForwardReverse) {
    stop("Both forward and reverse FASTQ files must be given in order to merge ",
         "forward and reverse reads")
  }
  
  ## check numeric inputs
  checkNumericInput(skipForward, nonnegative = FALSE)
  checkNumericInput(skipReverse, nonnegative = FALSE)
  checkNumericInput(umiLengthForward, nonnegative = FALSE)
  checkNumericInput(umiLengthReverse, nonnegative = FALSE)
  checkNumericInput(constantLengthForward, nonnegative = FALSE)
  checkNumericInput(constantLengthReverse, nonnegative = FALSE)
  checkNumericInput(variableLengthForward, nonnegative = FALSE)
  checkNumericInput(variableLengthReverse, nonnegative = FALSE)
  checkNumericInput(avePhredMinForward, nonnegative = TRUE)
  checkNumericInput(avePhredMinReverse, nonnegative = TRUE)
  checkNumericInput(variableNMaxForward, nonnegative = TRUE)
  checkNumericInput(variableNMaxReverse, nonnegative = TRUE)
  checkNumericInput(umiNMax, nonnegative = TRUE)
  checkNumericInput(nbrMutatedCodonsMaxForward, nonnegative = TRUE)
  checkNumericInput(nbrMutatedCodonsMaxReverse, nonnegative = TRUE)
  checkNumericInput(mutatedPhredMinForward, nonnegative = TRUE)
  checkNumericInput(mutatedPhredMinReverse, nonnegative = TRUE)
  checkNumericInput(variableCollapseMaxDist, nonnegative = TRUE)
  checkNumericInput(umiCollapseMaxDist, nonnegative = TRUE)
  
  ## adapters and primers must be strings, valid DNA characters
  if (!is.character(adapterForward) || length(adapterForward) != 1 ||
      !grepl("^[AaCcGgTt]*$", adapterForward) || !is.character(adapterReverse) ||
      length(adapterReverse) != 1 || !grepl("^[AaCcGgTt]*$", adapterReverse)) {
    stop("Adapters must be character strings, only containing valid DNA characters")
  } else {
    adapterForward <- toupper(adapterForward)
    adapterReverse <- toupper(adapterReverse)
  }
  
  if (!is.character(primerForward) || length(primerForward) != 1 ||
      !grepl("^[AaCcGgTt]*$", primerForward) || !is.character(primerReverse) ||
      length(primerReverse) != 1 || !grepl("^[AaCcGgTt]*$", primerReverse)) {
    stop("Primers must be character strings, only containing valid DNA characters")
  } else {
    primerForward <- toupper(primerForward)
    primerReverse <- toupper(primerReverse)
  }
  
  ## Check that a valid combination of sequence part lengths/primers is provided
  ## If one of skip, umiLength, constantLength is provided (not -1), the others must be given too
  if (any(c(skipForward, umiLengthForward, constantLengthForward) != (-1)) && 
      !all(c(skipForward, umiLengthForward, constantLengthForward) != (-1))) {
    stop("Either none or all of 'skipForward', 'umiLengthForward' and 'constantLengthForward' ",
         "must be specified (> -1)")
  }
  if (any(c(skipReverse, umiLengthReverse, constantLengthReverse) != (-1)) && 
      !all(c(skipReverse, umiLengthReverse, constantLengthReverse) != (-1))) {
    stop("Either none or all of 'skipReverse', 'umiLengthReverse' and 'constantLengthReverse' ",
         "must be specified (> -1)")
  }
  ## Now we know that either all or none of skip, umiLength, constantLength are set
  ## If they are set (!= -1), primers can not be given
  if (skipForward != (-1) && primerForward != "") {
    stop("Both sequence component lengths and primer sequence can not be set (forward)")
  }
  if (skipReverse != (-1) && primerReverse != "") {
    stop("Both sequence component lengths and primer sequence can not be set (reverse)")
  }
  ## If they are not set, primers must be given (unless the fastq file is not given)
  if (skipForward == (-1) && primerForward == "") {
    stop("Either sequence component lengths or primer sequence must be set (forward)")
  }
  if (fastqReverse != "" && skipReverse == (-1) && primerReverse == "") {
    stop("Either sequence component lengths or primer sequence must be set (reverse)")
  }
  
  ## if wild type sequence is a string, make it into a vector
  if (length(wildTypeForward) == 1) {
    wildTypeForward <- c(f = wildTypeForward)
  }
  if (length(wildTypeReverse) == 1) {
    wildTypeReverse <- c(r = wildTypeReverse)
  }
  
  ## wild type sequences must be given in named vectors
  if (is.null(names(wildTypeForward)) || any(names(wildTypeForward) == "") || 
      is.null(names(wildTypeReverse)) || any(names(wildTypeReverse) == "")) {
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
  
  ## wild type sequence lengths must match variable sequence lengths
  if (any(sapply(wildTypeForward, function(w) {
    variableLengthForward != (-1) && nchar(w) > 0 && nchar(w) != variableLengthForward
  }))) {
    stop("The lengths of the elements in 'wildTypeForward' (", paste(sapply(wildTypeForward, nchar), collapse = ","), 
         ") do not all correspond to the given 'variableLengthForward' (", 
         variableLengthForward, ")")
  }
  if (any(sapply(wildTypeReverse, function(w) {
    variableLengthReverse != (-1) && nchar(w) > 0 && nchar(w) != variableLengthReverse
  }))) {
    stop("The lengths of the elements in 'wildTypeReverse' (", paste(sapply(wildTypeReverse, nchar), collapse = ","), 
         ") do not all correspond to the given 'variableLengthReverse' (", 
         variableLengthReverse, ")")
  }
  
  ## cis experiment - should not have wildTypeReverse
  if (mergeForwardReverse && any(sapply(wildTypeReverse, nchar) > 0)) {
    warning("Ignoring 'wildTypeReverse' when forward and reverse reads are merged")
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
  
  ## mutNameDelimiter must be a single character, and can not appear in any of the WT sequence names
  if (!is.character(mutNameDelimiter) || length(mutNameDelimiter) != 1 || 
      nchar(mutNameDelimiter) != 1 || mutNameDelimiter == "_") {
    stop("'mutNameDelimiter' must be a single-letter character scalar, not equal to '_'")
  }
  if (any(grepl(mutNameDelimiter, c(names(wildTypeForward), names(wildTypeReverse)), fixed = TRUE))) {
    stop("'mutNameDelimiter' can not appear in the name of any of the provided wild type sequences.")
  }
  
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be a logical scalar.")
  }
  
  res <- digestFastqsCpp(fastqForwardVect = fastqForward, 
                         fastqReverseVect = fastqReverse,
                         mergeForwardReverse = mergeForwardReverse,
                         minOverlap = minOverlap, 
                         maxOverlap = maxOverlap, 
                         maxFracMismatchOverlap = maxFracMismatchOverlap,
                         greedyOverlap = greedyOverlap,
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
                         primerForward = primerForward,
                         primerReverse = primerReverse,
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
                         mutNameDelimiter = mutNameDelimiter,
                         variableCollapseMaxDist = variableCollapseMaxDist,
                         umiCollapseMaxDist = umiCollapseMaxDist,
                         maxNReads = maxNReads,
                         verbose = verbose)
  
  ## Add package version and processing date
  res$parameters$processingInfo <- paste0(
    "Processed by mutscan v", utils::packageVersion("mutscan"), " on ",
    Sys.time()
  )
  res
}
