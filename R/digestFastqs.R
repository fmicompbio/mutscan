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
#'  \item Extract the UMI, constant and variable sequence from forward and
#'  reverse reads, based on the definition of the respective read composition.
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
#'  codon. {.} is a delimiter, specified via \code{mutNameDelimiter}. If no
#'  wildtype sequences are provided, the read sequence will be used as the
#'  'mutation name'.
#'}
#'
#' Based on the retained reads following this filtering process, count the
#' number of reads, and the number of unique UMIs, for each variable sequence
#' (or pair of variable sequences).
#'
#' @param fastqForward,fastqReverse character vector, paths to gzipped FASTQ files
#'   corresponding to forward and reverse reads, respectively. If more than one
#'   forward/reverse sequence file is given, they need to be provided in the
#'   same order. Note that if multiple fastq files are provided, they are all
#'   assumed to correspond to the same sample, and will effectively be concatenated.
#' @param mergeForwardReverse logical(1), whether to fuse the forward and
#'   reverse variable sequences.
#' @param minOverlap,maxOverlap numeric(1), the minimal and maximal allowed
#'   overlap between the forward and reverse reads when merging. Only used if
#'   \code{mergeForwardReverse} is \code{TRUE}. If set to 0, only overlaps
#'   covering the full length of the shortest of the two reads will be
#'   considered.
#' @param minMergedLength,maxMergedLength numeric(1), the minimal and maximal
#'   allowed total length of the merged product (if \code{mergeForwardReverse} is
#'   \code{TRUE}). If set to 0, any length is allowed.
#' @param maxFracMismatchOverlap numeric(1), maximal mismatch rate in the
#'   overlap. Only used if \code{mergeForwardReverse} is \code{TRUE}.
#' @param greedyOverlap logical(1). If \code{TRUE}, the first overlap satisfying
#'   \code{minOverlap}, \code{maxOverlap}, \code{minMergedLength},
#'   \code{maxMergedLength} and \code{maxFracMismatchOverlap} will
#'   be retained. If \code{FALSE}, all valid overlaps will be scored and the one
#'   with the highest score (largest number of matches) will be retained.
#' @param revComplForward,revComplReverse logical(1), whether to reverse
#'   complement the forward/reverse variable and constant sequences, respectively.
#' @param elementsForward,elementsReverse Character strings representing the
#'   composition of the forward and reverse reads, respectively. The strings should
#'   consist only of the letters S (skip), C (constant), U (umi), P (primer),
#'   V (variable), and cover the full extent of the read. Most combinations
#'   are allowed (and a given letter can appear multiple times), but there
#'   can be at most one occurrence of P. If a given letter is included
#'   multiple times, the corresponding sequences will be concatenated in the output.
#' @param elementLengthsForward,elementLengthsReverse Numeric vectors containing
#'   the lengths of each read component from \code{elementsForward}/\code{elementsReverse},
#'   respectively. If the length of one element is set to -1, it will be inferred
#'   from the other lengths (as the remainder of the read). At most one number
#'   (or one number on each side of the primer P) can be set to -1. The indicated
#'   length of the primer is not used (instead it's inferred from the provided primer
#'   sequence) and can also be set to -1.
#' @param adapterForward,adapterReverse character(1), the adapter sequence for
#'   forward/reverse reads, respectively. If a forward/reverse read contains the
#'   corresponding adapter sequence, the sequence pair will be filtered out.
#'   If set to \code{NULL}, no adapter filtering is performed. The number of
#'   filtered read pairs are reported in the return value.
#' @param primerForward,primerReverse Character vectors, representing the primer
#'   sequence(s) for forward/reverse reads, respectively. Only read pairs that
#'   contain perfect matches to both the forward and reverse primers (if given)
#'   will be retained. Multiple primers can be specified - they will be
#'   considered in order and the first match will be used.
#' @param wildTypeForward,wildTypeReverse character(1) or named character
#'   vector, the wild type sequence for the forward and reverse variable region.
#'   If given as a single string, the reference sequence will be named 'f' (for
#'   forward) or 'r' (for reverse).
#' @param constantForward,constantReverse character(1) or character vector,
#'   the expected constant forward and reverse sequences. If a vector, all
#'   elements must have the same length.
#' @param avePhredMinForward,avePhredMinReverse numeric(1) Minimum average Phred
#'   score in the variable region for a read to be retained. If L contains both
#'   forward and reverse variable regions, the minimum average Phred score has
#'   to be achieved in both for a read pair to be retained.
#' @param variableNMaxForward,variableNMaxReverse numeric(1) Maximum number of
#'   Ns allowed in the variable region for a read to be retained.
#' @param umiNMax numeric(1) Maximum number of Ns allowed in the UMI for a read
#'   to be retained.
#' @param nbrMutatedCodonsMaxForward,nbrMutatedCodonsMaxReverse numeric(1)
#'   Maximum number of mutated codons that are allowed. Note that for the
#'   forward and reverse sequence, respectively, exactly one of
#'   \code{nbrMutatedCodonsMax} and \code{nbrMutatedBasesMax} must be -1,
#'   and the other must be a non-negative number. The one that is not -1
#'   will be used to filter and name the identified mutants.
#' @param nbrMutatedBasesMaxForward,nbrMutatedBasesMaxReverse numeric(1)
#'   Maximum number of mutated bases that are allowed. Note that for the
#'   forward and reverse sequence, respectively, exactly one of
#'   \code{nbrMutatedCodonsMax} and \code{nbrMutatedBasesMax} must be -1,
#'   and the other must be a non-negative number. The one that is not -1
#'   will be used to filter and name the identified mutants.
#' @param forbiddenMutatedCodonsForward,forbiddenMutatedCodonsReverse character
#'   vector of codons (can contain ambiguous IUPAC characters, see
#'   \code{\link[Biostrings]{IUPAC_CODE_MAP}}). If a read pair contains a
#'   mutated codon matching this pattern, it will be filtered out.
#' @param useTreeWTmatch logical(1). Should a tree-based matching
#'   to wild type sequences be used if possible? If the number of allowed
#'   mismatches is small, and the number of wild type sequences is large,
#'   this is typically faster.
#' @param collapseToWTForward,collapseToWTReverse Logical scalar, indicating
#'   whether to just represent the observed variable sequence by the
#'   closest wildtype sequence rather than retaining the information about
#'   the mutations.
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
#' @param constantMaxDistForward,constantMaxDistReverse numeric(1) maximum
#'   allowed Hamming distance between the extracted and expected constant sequence.
#'   If multiple constant sequences are provided, the most similar one is used.
#'   Reads with a larger distance to the expected constant sequence are
#'   discarded. If set to -1, no filtering is done.
#' @param variableCollapseMaxDist,umiCollapseMaxDist A \code{numeric} scalar
#'   defining the tolerances for collapsing similar variable (forward + reverse,
#'   if any) or UMI sequences. If the value is in [0, 1), it defines the maximal
#'   Hamming distance in terms of a fraction of sequence length:
#'   (\code{round(variableCollapseMaxDist * nchar(variableSeq))}).
#'   A value greater or equal to 1 is rounded and directly used as the maximum
#'   allowed Hamming distance. Note that variable sequences can only be
#'   collapsed if they are all of the same length and no wild type sequences
#'   (\code{wildTypeForward} or \code{wildTypeReverse}) have been given.
#' @param variableCollapseMinReads A \code{numeric} scalar. Indicates the minimum
#'   number of reads for the read to be considered for collapsing with similar
#'   sequences.
#' @param variableCollapseMinRatio A \code{numeric} scalar. During collapsing of
#'   similar variable sequences, a low-frequency sequence will be collapsed with
#'   a higher-frequency sequence only if the ratio between the high-frequency
#'   and the low-frequency sequences is at least this high. The default value
#'   of 0 indicates that no such check is performed.
#' @param filteredReadsFastqForward,filteredReadsFastqReverse character(1), the
#'   name of a (pair of) FASTQ file(s) where filtered-out reads will be written.
#'   The name(s) should end in .gz (the output will always be compressed). If
#'   empty, filtered reads will not be written to a file.
#' @param maxNReads integer(1) Maximal number of reads to process. If set to -1,
#'   all reads will be processed.
#' @param verbose logical(1), whether to print out progress messages.
#' @param nThreads numeric(1), number of threads to use for parallel processing.
#' @param chunkSize numeric(1), number of read (pairs) to keep in memory for
#'   parallel processing. Reduce from the default value if you run out of memory.
#'
#' @return A list with four entries:
#' \describe{
#' \item{summaryTable}{A \code{data.frame} that contains, for each observed
#' mutation combination, the corresponding variable region sequences (or pair of
#' sequences), the number of observed such sequences, and the number of unique
#' UMIs observed for the sequence. It also has a column named 'maxNbrReads',
#' which contains the number of reads for the most frequent observed sequence
#' represented by the feature (only relevant if similar variable regions are
#' collapsed).}
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
#' @import zlibbioc
digestFastqs <- function(fastqForward, fastqReverse = NULL,
                         mergeForwardReverse = FALSE, minOverlap = 0, maxOverlap = 0,
                         minMergedLength = 0, maxMergedLength = 0,
                         maxFracMismatchOverlap = 1, greedyOverlap = TRUE,
                         revComplForward = FALSE, revComplReverse = FALSE,
                         adapterForward = "", adapterReverse = "",
                         elementsForward = "",
                         elementLengthsForward = numeric(0),
                         elementsReverse = "",
                         elementLengthsReverse = numeric(0),
                         primerForward = c(""), primerReverse = c(""),
                         wildTypeForward = "", wildTypeReverse = "",
                         constantForward = c(""), constantReverse = c(""),
                         avePhredMinForward = 20.0, avePhredMinReverse = 20.0,
                         variableNMaxForward = 0, variableNMaxReverse = 0,
                         umiNMax = 0,
                         nbrMutatedCodonsMaxForward = 1,
                         nbrMutatedCodonsMaxReverse = 1,
                         nbrMutatedBasesMaxForward = -1,
                         nbrMutatedBasesMaxReverse = -1,
                         forbiddenMutatedCodonsForward = "",
                         forbiddenMutatedCodonsReverse = "",
                         useTreeWTmatch = FALSE,
                         collapseToWTForward = FALSE,
                         collapseToWTReverse = FALSE,
                         mutatedPhredMinForward = 0.0,
                         mutatedPhredMinReverse = 0.0,
                         mutNameDelimiter = ".",
                         constantMaxDistForward = -1,
                         constantMaxDistReverse = -1,
                         variableCollapseMaxDist = 0.0,
                         variableCollapseMinReads = 0,
                         variableCollapseMinRatio = 0,
                         umiCollapseMaxDist = 0.0,
                         filteredReadsFastqForward = "",
                         filteredReadsFastqReverse = "",
                         maxNReads = -1, verbose = FALSE,
                         nThreads = 1, chunkSize = 100000) {
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
  .assertScalar(x = mergeForwardReverse, type = "logical")
  .assertScalar(x = revComplForward, type = "logical")
  .assertScalar(x = revComplReverse, type = "logical")
  if (any(fastqReverse == "") && mergeForwardReverse) {
    stop("Both forward and reverse FASTQ files must be given in order to merge ",
         "forward and reverse reads")
  }

  ## check numeric inputs
  .assertVector(x = elementLengthsForward, type = "numeric",
                rngIncl = c(0, Inf), validValues = -1)
  .assertVector(x = elementLengthsReverse, type = "numeric",
                rngIncl = c(0, Inf), validValues = -1)
  .assertScalar(x = minOverlap, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = maxOverlap, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = minMergedLength, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = maxMergedLength, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = avePhredMinForward, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = avePhredMinReverse, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = variableNMaxForward, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = variableNMaxReverse, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = umiNMax, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = nbrMutatedCodonsMaxForward, type = "numeric",
                rngIncl = c(0, Inf), validValues = -1)
  .assertScalar(x = nbrMutatedCodonsMaxReverse, type = "numeric",
                rngIncl = c(0, Inf), validValues = -1)
  .assertScalar(x = nbrMutatedBasesMaxForward, type = "numeric",
                rngIncl = c(0, Inf), validValues = -1)
  .assertScalar(x = nbrMutatedBasesMaxReverse, type = "numeric",
                rngIncl = c(0, Inf), validValues = -1)
  .assertScalar(x = mutatedPhredMinForward, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = mutatedPhredMinReverse, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = constantMaxDistForward, type = "numeric",
                rngIncl = c(0, Inf), validValues = -1)
  .assertScalar(x = constantMaxDistReverse, type = "numeric",
                rngIncl = c(0, Inf), validValues = -1)
  .assertScalar(x = variableCollapseMaxDist, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = variableCollapseMinReads, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = variableCollapseMinRatio, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = umiCollapseMaxDist, type = "numeric", rngIncl = c(0, Inf))
  .assertScalar(x = maxNReads, type = "numeric", rngIncl = c(0, Inf),
                validValues = -1)
  .assertScalar(x = nThreads, type = "numeric", rngExcl = c(0, Inf))
  .assertScalar(x = chunkSize, type = "numeric", rngExcl = c(0, Inf))

  ## If a wildtype sequence is provided, it must be unambiguous how to identify and name mutants
  if (any(wildTypeForward != "")) {
    if ((nbrMutatedCodonsMaxForward == (-1) && nbrMutatedBasesMaxForward == (-1)) ||
        (nbrMutatedCodonsMaxForward != (-1) && nbrMutatedBasesMaxForward != (-1))) {
      stop("Exactly one of 'nbrMutatedCodonsMaxForward' and 'nbrMutatedBasesMaxForward' must be -1")
    }
  }
  if (any(wildTypeReverse != "")) {
    if ((nbrMutatedCodonsMaxReverse == (-1) && nbrMutatedBasesMaxReverse == (-1)) ||
        (nbrMutatedCodonsMaxReverse != (-1) && nbrMutatedBasesMaxReverse != (-1))) {
      stop("Exactly one of 'nbrMutatedCodonsMaxReverse' and 'nbrMutatedBasesMaxReverse' must be -1")
    }
  }

  ## adapters must be strings, valid DNA characters
  if (length(adapterForward) != 1 || !is.character(adapterForward) ||
      !grepl("^[AaCcGgTt]*$", adapterForward) || length(adapterReverse) != 1 ||
      !is.character(adapterReverse) || !grepl("^[AaCcGgTt]*$", adapterReverse)) {
    stop("Adapters must be character strings, only containing valid DNA characters")
  } else {
    adapterForward <- toupper(adapterForward)
    adapterReverse <- toupper(adapterReverse)
  }

  ## Check provided read composition
  if (!(length(elementsForward) == 1 && is.character(elementsForward))) {
    stop("'elementsForward' must be a character scalar")
  }
  if (!(length(elementsReverse) == 1 && is.character(elementsReverse))) {
    stop("'elementsReverse' must be a character scalar")
  }

  if (nchar(elementsForward) == 0) {
    stop("'elementsForward' must be a non-empty character scalar")
  }
  if (any(fastqReverse != "") && nchar(elementsReverse) == 0) {
    stop("'elementsReverse' must be a non-empty character scalar")
  }

  if (!all(strsplit(elementsForward, "")[[1]] %in% c("C", "U", "S", "V", "P"))) {
    stop("'elementsForward' can only contain letters 'CUSVP'")
  }
  if (!all(strsplit(elementsReverse, "")[[1]] %in% c("C", "U", "S", "V", "P"))) {
    stop("'elementsReverse' can only contain letters 'CUSVP'")
  }

  if (nchar(elementsForward) != length(elementLengthsForward)) {
    stop("'elementsForward' and 'elementsLengthsForward' must have the same length")
  }
  if (nchar(elementsReverse) != length(elementLengthsReverse)) {
    stop("'elementsReverse' and 'elementsLengthsReverse' must have the same length")
  }

  ## Max one 'P'
  PposFwd <- gregexpr(pattern = "P", elementsForward, fixed = TRUE)[[1]]
  PposRev <- gregexpr(pattern = "P", elementsReverse, fixed = TRUE)[[1]]

  if (length(PposFwd) > 1) {
    stop("'elementsForward' can contain max one 'P'")
  }
  if (length(PposRev) > 1) {
    stop("'elementsReverse' can contain max one 'P'")
  }

  ## If no 'P', max one -1 length
  if (length(PposFwd) == 1 && PposFwd == -1 &&
      sum(elementLengthsForward == -1) > 1) {
    stop("Max one element length (forward) can be -1")
  }
  if (length(PposRev) == 1 && PposRev == -1 &&
      sum(elementLengthsReverse == -1) > 1) {
    stop("Max one element length (reverse) can be -1")
  }

  ## If a 'P', max one -1 length on each side
  if (length(PposFwd) == 1 && PposFwd != -1) {
    if ((PposFwd != 1 && sum(elementLengthsForward[1:(PposFwd - 1)] == -1) > 1) ||
        (PposFwd != nchar(elementsForward) &&
         sum(elementLengthsForward[(PposFwd + 1):nchar(elementsForward)] == -1) > 1)) {
      stop("Max one element length (forward) on each side of the primer can be -1")
    }
  }
  if (length(PposRev) == 1 && PposRev != -1) {
    if ((PposRev != 1 && sum(elementLengthsReverse[1:(PposRev - 1)] == -1) > 1) ||
        (PposRev != nchar(elementsReverse) &&
         sum(elementLengthsReverse[(PposRev + 1):nchar(elementsReverse)] == -1) > 1)) {
      stop("Max one element length (reverse) on each side of the primer can be -1")
    }
  }

  if (!is.character(primerForward) || length(primerForward) < 1 ||
      !all(grepl("^[AaCcGgTt]*$", primerForward)) || !is.character(primerReverse) ||
      length(primerReverse) < 1 || !all(grepl("^[AaCcGgTt]*$", primerReverse))) {
    stop("Primers must be character vectors, only containing valid DNA characters")
  } else {
    primerForward <- vapply(primerForward, toupper, "")
    primerReverse <- vapply(primerReverse, toupper, "")
  }

  ## if wild type sequence is a string, make it into a vector
  if (length(wildTypeForward) == 1 && is.null(names(wildTypeForward))) {
    wildTypeForward <- c(f = wildTypeForward)
  }
  if (length(wildTypeReverse) == 1 && is.null(names(wildTypeReverse))) {
    wildTypeReverse <- c(r = wildTypeReverse)
  }

  ## wild type sequences must be given in named vectors
  if (any(is.null(names(wildTypeForward))) || any(names(wildTypeForward) == "") ||
      any(is.null(names(wildTypeReverse))) || any(names(wildTypeReverse) == "")) {
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

  ## all wild type sequences should be unique
  if (any(duplicated(wildTypeForward)) ||
      any(duplicated(wildTypeReverse))) {
    stop("Duplicated wild type sequences are not allowed")
  }

  ## cis experiment - should not have wildTypeReverse
  if (mergeForwardReverse && any(sapply(wildTypeReverse, nchar) > 0)) {
    warning("Ignoring 'wildTypeReverse' when forward and reverse reads are merged")
    wildTypeReverse <- c(r = "")
  }

  ## if both constantForward and constantLengthForward are given, check that
  ## the lengths inferred from the two are consistent
  if (!is.character(constantForward) || !is.character(constantReverse) ||
      length(constantForward) < 1 || length(constantReverse) < 1 ||
      !all(grepl("^[AaCcGgTt]*$", constantForward)) ||
      !all(grepl("^[AaCcGgTt]*$", constantReverse))) {
    stop("constant sequences must be character strings, ",
         "only containing valid DNA characters.")
  } else {
    constantForward <- toupper(constantForward)
    constantReverse <- toupper(constantReverse)
  }

  if (!all(vapply(constantForward, nchar, 0L) == nchar(constantForward[1])) ||
      !all(vapply(constantReverse, nchar, 0L) == nchar(constantReverse[1]))) {
    stop("All constant sequences must be of the same length")
  }

  CposFwd <- gregexpr(pattern = "C", elementsForward)[[1]]
  if (nchar(constantForward[1]) > 0 &&
      sum(elementLengthsForward[CposFwd]) != nchar(constantForward[1])) {
    stop("The sum of the constant sequence lengths in elementsForward (",
         sum(elementLengthsForward[CposFwd]),
         ") does not correspond to the length of the given 'constantForward' (",
         nchar(constantForward[1]), ")")
  }
  CposRev <- gregexpr(pattern = "C", elementsReverse)[[1]]
  if (nchar(constantReverse[1]) > 0 &&
      sum(elementLengthsReverse[CposRev]) != nchar(constantReverse[1])) {
    stop("The sum of the constant sequence lengths in elementsReverse (",
         sum(elementLengthsReverse[CposRev]),
         ") does not correspond to the length of the given 'constantReverse' (",
         nchar(constantReverse[1]), ")")
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

  if (!is.logical(useTreeWTmatch) || length(useTreeWTmatch) != 1) {
    stop("'useTreeWTmatch' must be a logical scalar.")
  }

  if (!is.logical(collapseToWTForward) || length(collapseToWTForward) != 1) {
    stop("'collapseToWTForward' must be a logical scalar.")
  }
  if (!is.logical(collapseToWTReverse) || length(collapseToWTReverse) != 1) {
    stop("'collapseToWTReverse' must be a logical scalar.")
  }

  ## mutNameDelimiter must be a single character, and can not appear in any of the WT sequence names
  if (!is.character(mutNameDelimiter) || length(mutNameDelimiter) != 1 ||
      nchar(mutNameDelimiter) != 1 || mutNameDelimiter == "_") {
    stop("'mutNameDelimiter' must be a single-letter character scalar, not equal to '_'")
  }
  if (any(grepl(mutNameDelimiter, c(names(wildTypeForward), names(wildTypeReverse)), fixed = TRUE))) {
    stop("'mutNameDelimiter' can not appear in the name of any of the provided wild type sequences.")
  }

  if (!is.character(filteredReadsFastqForward) || !is.character(filteredReadsFastqReverse) ||
      length(filteredReadsFastqForward) != 1 || length(filteredReadsFastqReverse) != 1 ||
      (filteredReadsFastqForward != "" && !grepl("\\.gz$", filteredReadsFastqForward)) ||
      (filteredReadsFastqReverse != "" && !grepl("\\.gz$", filteredReadsFastqReverse))) {
    stop("'filteredReadsFastqForward' and 'filteredReadsFastqReverse' must be character ",
         "scalars ending with .gz.")
  }

  if ((any(fastqReverse == "") && filteredReadsFastqReverse != "") ||
      (all(fastqReverse != "") && filteredReadsFastqForward != "" && filteredReadsFastqReverse == "") ||
      (all(fastqForward != "") && filteredReadsFastqForward == "" && filteredReadsFastqReverse != "")) {
    stop("The pairing of the output FASTQ files must be compatible with that of the input files.")
  }

  .assertScalar(x = verbose, type = "logical")

  ## Represent the wildtype sequences as pairs of vectors (one with sequences,
  ## one with names), to make things faster on the C++ side
  wildTypeForwardNames <- names(wildTypeForward)
  wildTypeForward <- unname(wildTypeForward)
  wildTypeReverseNames <- names(wildTypeReverse)
  wildTypeReverse <- unname(wildTypeReverse)

  res <- digestFastqsCpp(fastqForwardVect = fastqForward,
                         fastqReverseVect = fastqReverse,
                         mergeForwardReverse = mergeForwardReverse,
                         minOverlap = minOverlap,
                         maxOverlap = maxOverlap,
                         minMergedLength = minMergedLength,
                         maxMergedLength = maxMergedLength,
                         maxFracMismatchOverlap = maxFracMismatchOverlap,
                         greedyOverlap = greedyOverlap,
                         revComplForward = revComplForward,
                         revComplReverse = revComplReverse,
                         elementsForward = elementsForward,
                         elementLengthsForward = as.numeric(elementLengthsForward),
                         elementsReverse = elementsReverse,
                         elementLengthsReverse = as.numeric(elementLengthsReverse),
                         adapterForward = adapterForward,
                         adapterReverse = adapterReverse,
                         primerForward = primerForward,
                         primerReverse = primerReverse,
                         wildTypeForward = wildTypeForward,
                         wildTypeForwardNames = wildTypeForwardNames,
                         wildTypeReverse = wildTypeReverse,
                         wildTypeReverseNames = wildTypeReverseNames,
                         constantForward = constantForward,
                         constantReverse = constantReverse,
                         avePhredMinForward = avePhredMinForward,
                         avePhredMinReverse = avePhredMinReverse,
                         variableNMaxForward = variableNMaxForward,
                         variableNMaxReverse = variableNMaxReverse,
                         umiNMax = umiNMax,
                         nbrMutatedCodonsMaxForward = nbrMutatedCodonsMaxForward,
                         nbrMutatedCodonsMaxReverse = nbrMutatedCodonsMaxReverse,
                         nbrMutatedBasesMaxForward = nbrMutatedBasesMaxForward,
                         nbrMutatedBasesMaxReverse = nbrMutatedBasesMaxReverse,
                         forbiddenMutatedCodonsForward = forbiddenMutatedCodonsForward,
                         forbiddenMutatedCodonsReverse = forbiddenMutatedCodonsReverse,
                         useTreeWTmatch = useTreeWTmatch,
                         collapseToWTForward = collapseToWTForward,
                         collapseToWTReverse = collapseToWTReverse,
                         mutatedPhredMinForward = mutatedPhredMinForward,
                         mutatedPhredMinReverse = mutatedPhredMinReverse,
                         mutNameDelimiter = mutNameDelimiter,
                         constantMaxDistForward = constantMaxDistForward,
                         constantMaxDistReverse = constantMaxDistReverse,
                         variableCollapseMaxDist = variableCollapseMaxDist,
                         variableCollapseMinReads = as.integer(ceiling(variableCollapseMinReads)),
                         variableCollapseMinRatio = variableCollapseMinRatio,
                         umiCollapseMaxDist = umiCollapseMaxDist,
                         filteredReadsFastqForward = filteredReadsFastqForward,
                         filteredReadsFastqReverse = filteredReadsFastqReverse,
                         maxNReads = maxNReads,
                         verbose = verbose,
                         nThreads = as.integer(nThreads),
                         chunkSize = as.integer(chunkSize))

  ## Add package version and processing date
  res$parameters$processingInfo <- paste0(
    "Processed by mutscan v", utils::packageVersion("mutscan"), " on ",
    Sys.time()
  )
  res
}
