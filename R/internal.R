#' Check if "sequences" argument is valid
#' 
#' Currently: check if it is a list with expected elements.
#' 
#' @param L A list as returned by \code{readFastqsTrans}
#' 
isValidL <- function(L) {
  if (!is(L, "list") ||
      !identical(names(L), c("umis","constantSeqForward","constantSeqReverse",
                             "variableSeqForward","variableSeqReverse","readSummary")) ||
      !is(L$umis, "QualityScaledDNAStringSet") ||
      !is(L$constantSeqForward, "QualityScaledDNAStringSet") ||
      !is(L$constantSeqReverse, "QualityScaledDNAStringSet") ||
      !is(L$variableSeqForward, "QualityScaledDNAStringSet") ||
      !is(L$variableSeqReverse, "QualityScaledDNAStringSet") ||
      !is(L$readSummary, "data.frame")) {
    stop("'L' must be a 6-element list as returned by 'readFastqsTrans'.")
  }
  return(invisible(TRUE))
}


#' Find positions of mismatching bases in DNAStringSet
#' 
#' Compare sequences in a DNAStringSet to expected sequence
#' and return a list with positions of mismatches.
#' 
#' @param pattern A nucleotide sequence corresponding to the 
#'   wildtype sequence
#' @param subject A DNAStringSet with observed variable sequences
#' 
#' @return A two-element list with base/codon positions of mismatches. The
#'   list elements have classes InterList and NumericList, respectively.
#'   
#' @importFrom Biostrings DNAStringSet quality
#' @importFrom BiocGenerics unlist width relist
#' @importFrom IRanges PartitioningByEnd
#' 
findMismatchPositions <- function(pattern, subject) {
  ## Code from https://support.bioconductor.org/p/63460/
  if (!is(pattern, "character") || length(pattern) != 1L) {
    stop("'pattern' must be a length-1 character vector")
  }
  if (!is(subject, "DNAStringSet")) {
    stop("'subject' must be a DNAStringSet object")
  }
  
  pattern <- DNAStringSet(rep(pattern, length(subject)))
  
  pattern_width <- width(pattern)
  subject_width <- width(subject)
  
  unlisted_ans <- which(as.raw(unlist(pattern)) != as.raw(unlist(subject)))
  breakpoints <- cumsum(pattern_width)
  ans_eltlens <- tabulate(findInterval(unlisted_ans - 1L,
                                       breakpoints) + 1L,
                          nbins = length(pattern))
  skeleton <- IRanges::PartitioningByEnd(cumsum(ans_eltlens))
  nucleotideMismatches <- relist(unlisted_ans, skeleton)
  offsets <- c(0L, breakpoints[-length(breakpoints)])
  nucleotideMismatches <- nucleotideMismatches - offsets
  
  codonMismatches <- (nucleotideMismatches - 1) %/% 3 + 1
  
  list(nucleotideMismatches = nucleotideMismatches,
       codonMismatches = codonMismatches)
}


#' Count numbers of matching/mismatching bases per quality score
#' 
#' Based on the code from findMismatchpositions, tabulate base qualities
#' for matching/mismatching bases.
#' 
#' @param pattern A nucleotide sequence corresponding to the 
#'   expected sequence
#' @param subject A DNAStringSet with observed variable sequences
#' @param maxQ Maximum quality to tabulate (\code{nbins} argument to \code{tabulate})
#' 
#' @return a two-element list with the numbers of matching/mismatching bases of a given quality
#'  
#' @importFrom Biostrings DNAStringSet quality
#' @importFrom BiocGenerics unlist
#' 
tabulateQualitiesByMatchstate <- function(pattern, subject, maxQ = 99L) {
  if (!is(pattern, "character") || length(pattern) != 1L) {
    stop("'pattern' must be a length-1 character vector")
  }
  if (!is(subject, "DNAStringSet")) {
    stop("'subject' must be a DNAStringSet object")
  }
  
  pattern <- DNAStringSet(rep(pattern, length(subject)))
  mismatch_pos <- which(as.raw(unlist(pattern)) != as.raw(unlist(subject)))
  unlisted_quals <- unlist(as(Biostrings::quality(subject), "IntegerList"))
  
  tabErr <- tabulate(unlisted_quals[mismatch_pos], nbins = maxQ)
  tabCor <- tabulate(unlisted_quals[-mismatch_pos], nbins = maxQ)
  names(tabErr) <- names(tabCor) <- as.character(seq.int(maxQ))
  
  list(error = tabErr, correct = tabCor)
}
