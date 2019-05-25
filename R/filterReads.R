#' Filter FASTQ files from CIS or TRANS experiments
#'
#' @param L SummarizedExperiment, as output by \code{\link{readFastqs}}.
#' @param avePhredMin numeric(1) Minimum average Phred score in the variable
#'   region for a read to be retained. If L contains both forward and reverse
#'   variable regions, the minimum average Phred score has to be achieved in
#'   both for a read pair to be retained.
#' @param variableNMax numeric(1) Maximum number of Ns allowed in the variable
#'   region for a read to be retained.
#' @param umiNMax numeric(1) Maximum number of Ns allowed in the UMI for a read
#'   to be retained.
#' @param wildTypeForward,wildTypeReverse character(1), the wild type sequence
#'   for the forward and reverse variable region.
#' @param nbrMutatedCodonsMax numeric(1) Maximum number of mutated codons that
#'   are allowed.
#' @param forbiddenMutatedCodon character(1) A codon (can contain ambiguous
#'   IUPAC characters, see \code{\link[Biostrings]{IUPAC_CODE_MAP}}). If a read
#'   pair contains a mutated codon matching this pattern, it will be filtered
#'   out.
#' @param maxChunkSize numeric(1), largest allowed chunk size for steps where
#'   the reads are processed in chunks.
#' @param verbose logical(1), whether to print out progress messages.
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
#' @examples
#' datadir <- system.file("extdata", package = "mutscan")
#' transInput <- readFastqs(experimentType = "trans",
#'                          fastqForward = file.path(datadir, "transInput_1.fastq.gz"),
#'                          fastqReverse = file.path(datadir, "transInput_2.fastq.gz"),
#'                          skipForward = 1, skipReverse = 1, umiLengthForward = 10,
#'                          umiLengthReverse = 8, constantLengthForward = 18,
#'                          constantLengthReverse = 20, adapterForward = "GGAAGAGCACACGTC",
#'                          adapterReverse = "GGAAGAGCGTCGTGT", verbose = TRUE)
#' transInputFiltered <- filterReads(transInput, avePhredMin = 20,
#'                                   variableNMax = 0, umiNMax = 0,
#'                                   wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
#'                                   wildTypeReverse = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT",
#'                                   nbrMutatedCodonsMax = 1, forbiddenMutatedCodon = "NNW")
#'                                   
filterReads <- function(L, avePhredMin = 20, variableNMax = 0, umiNMax = 0,
                        wildTypeForward = NULL, wildTypeReverse = NULL, 
                        nbrMutatedCodonsMax = 1, forbiddenMutatedCodon = NULL,
                        maxChunkSize = 1e5, verbose = FALSE) {
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
    if (verbose) {
      message("Calculating average phred score in variable regions...")
    }
    if ("variableSeqForward" %in% SummarizedExperiment::assayNames(L)) {
      if (verbose) {
        message("...forward...")
      }
      avePhredVariableForward <- ShortRead::alphabetScore(
        Biostrings::quality(assay(L, "variableSeqForward")$seq)
      )/BiocGenerics::width(assay(L, "variableSeqForward")$seq) 
      keepAvePhredVariableForward <- avePhredVariableForward >= avePhredMin
    }
    if ("variableSeqReverse" %in% SummarizedExperiment::assayNames(L)) {
      if (verbose) {
        message("...reverse...")
      }
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
    if (verbose) {
      message("Counting the number of Ns in variable regions...")
    }
    if ("variableSeqForward" %in% SummarizedExperiment::assayNames(L)) {
      if (verbose) {
        message("...forward...")
      }
      variableNbrNForward <- Biostrings::vcountPattern(
        pattern = "N", subject = assay(L, "variableSeqForward")$seq,
        max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, 
        fixed = TRUE, algorithm = "auto")
      keepVariableNbrNForward <- variableNbrNForward <= variableNMax
    }
    if ("variableSeqReverse" %in% SummarizedExperiment::assayNames(L)) {
      message("...reverse...")
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
    if (verbose) {
      message("Counting the number of Ns in UMIs...")
    }
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
  if (!is.null(wildTypeForward)) {
    if ("variableSeqForward" %in% SummarizedExperiment::assayNames(L)) {
      if (verbose) {
        message("Finding mismatch positions, forward...")
      }
      filterForward <- mismatchFilter(wildTypeSeq = wildTypeForward, 
                                      variableSeq = assay(L, "variableSeqForward")$seq, 
                                      nbrMutatedCodonsMax = nbrMutatedCodonsMax,
                                      forbiddenMutatedCodon = forbiddenMutatedCodon,
                                      maxChunkSize = maxChunkSize, verbose = verbose)
      SummarizedExperiment::rowData(L)$qualMutatedForward <- filterForward$qualMutated
      SummarizedExperiment::rowData(L)$basesWithMutationsForward <- filterForward$basesWithMutations
      SummarizedExperiment::rowData(L)$codonsWithMutationsForward <- filterForward$codonsWithMutations
      SummarizedExperiment::rowData(L)$encodedMutatedCodonsForward <- filterForward$encodedMutCodons
      keepWildTypeForward <- filterForward$keepMM
    }
  }
  if (!is.null(wildTypeReverse)) {
    if ("variableSeqReverse" %in% SummarizedExperiment::assayNames(L)) {
      if (verbose) {
        message("Finding mismatch positions, reverse...")
      }
      filterReverse <- mismatchFilter(wildTypeSeq = wildTypeReverse, 
                                      variableSeq = assay(L, "variableSeqReverse")$seq, 
                                      nbrMutatedCodonsMax = nbrMutatedCodonsMax,
                                      forbiddenMutatedCodon = forbiddenMutatedCodon,
                                      maxChunkSize = maxChunkSize, verbose = verbose)
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
  if (verbose) {
    message("Filtering...")
  }
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
#' @param maxChunkSize numeric(1), largest allowed chunk size for steps where
#'   the reads are processed in chunks.
#' @param verbose logical(1), whether to print out progress messages.
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
#' @importFrom S4Vectors elementNROWS endoapply unstrsplit
#' @importFrom IRanges start
#' 
mismatchFilter <- function(wildTypeSeq, variableSeq, nbrMutatedCodonsMax = NULL,
                           forbiddenMutatedCodon = NULL, maxChunkSize = 1e5, 
                           verbose = FALSE) {
  keepMM <- rep(TRUE, length(variableSeq))
  
  ## --------------------------------------------------------------------------
  ## Find mismatch positions and qualities of mismatched bases
  ## --------------------------------------------------------------------------
  if (verbose) {
    message("Finding mismatch positions...")
  }
  mmp <- findMismatchPositions(pattern = wildTypeSeq, 
                               subject = variableSeq,
                               maxChunkSize = maxChunkSize)
  
  if (verbose) {
    message("Extracting qualities of mutated bases...")
  }
  m <- factor(seq_along(variableSeq) %/% maxChunkSize)
  qualMutated <- do.call(c, lapply(levels(m), function(i) {
    as(Biostrings::quality(variableSeq[m == i])[mmp$nucleotideMismatches[which(m == i)]], "IntegerList")
  }))

  ## --------------------------------------------------------------------------
  ## Count number of mutated codons per read
  ## --------------------------------------------------------------------------
  if (!is.null(nbrMutatedCodonsMax)) {
    if (verbose) {
      message("Counting number of mutated codons per read...")
    }
    nbrMutCodons <- S4Vectors::elementNROWS(unique(mmp$codonMismatches))
    keepMM <- nbrMutCodons <= nbrMutatedCodonsMax
  }
  
  ## --------------------------------------------------------------------------
  ## Derive a characterisation of each sequence, in terms of the codons with 
  ## mismatches
  ## --------------------------------------------------------------------------
  if (verbose) {
    message("Encoding sequence variants...")
  }
  uniqueMutCodons <- unique(mmp$codonMismatches)
  mutCodons <- S4Vectors::endoapply(uniqueMutCodons, 
                                    function(v) rep(3 * v, each = 3) + c(-2, -1, 0))
  
  variableSeq <- as(variableSeq, "DNAStringSet")
  m <- factor(seq_along(variableSeq) %/% maxChunkSize)
  seqMutCodons <- do.call(c, lapply(levels(m), function(i) {
    variableSeq[m == i][mutCodons[which(m == i)]]
  }))
    
  splitMutCodons <- as(strsplit(gsub("([^.]{3})", "\\1\\.", as.character(seqMutCodons)), ".", fixed = TRUE), "CharacterList")
  encodedMutCodons <- S4Vectors::unstrsplit(BiocGenerics::relist(paste0("c", unlist(uniqueMutCodons), unlist(splitMutCodons)), 
                                                                 uniqueMutCodons),
                                            sep = "_")

  ## --------------------------------------------------------------------------
  ## Find any occurrences of the forbidden mutated codon
  ## --------------------------------------------------------------------------
  if (!is.null(forbiddenMutatedCodon)) {
    if (verbose) {
      message("Finding occurrences of forbidden codons...")
    }
    matchForbidden <- Biostrings::vmatchPattern(
      pattern = forbiddenMutatedCodon, 
      as(seqMutCodons, "DNAStringSet"), 
      fixed = FALSE)
    keepMM <- keepMM & !any(relist(start(unlist(matchForbidden)) %% 3==1, 
                                   matchForbidden))
  }
  
  list(qualMutated = qualMutated, keepMM = keepMM,
       basesWithMutations = mmp$nucleotideMismatches,
       codonsWithMutations = mmp$codonMismatches,
       encodedMutCodons = encodedMutCodons)
}