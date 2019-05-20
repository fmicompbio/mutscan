#' @param pattern A nucleotide sequence corresponding to the 
#'   wildtype sequence
#' @param subject A DNAStringSet with observed variable sequences
#' 
findMismatchPositions <- function(pattern, subject) {
  ## Code from https://support.bioconductor.org/p/63460/
  pattern <- DNAStringSet(rep(pattern, length(subject)))
  
  if (!is(subject, "DNAStringSet")) {
    stop("'subject' must be a DNAStringSet object")
  }
  
  pattern_width <- width(pattern)
  subject_width <- width(subject)
  
  unlisted_ans <- which(as.raw(unlist(pattern)) != as.raw(unlist(subject)))
  breakpoints <- cumsum(pattern_width)
  ans_eltlens <- tabulate(findInterval(unlisted_ans - 1L,
                                       breakpoints) + 1L,
                          nbins = length(pattern))
  skeleton <- PartitioningByEnd(cumsum(ans_eltlens))
  nucleotideMismatches <- relist(unlisted_ans, skeleton)
  offsets <- c(0L, breakpoints[-length(breakpoints)])
  nucleotideMismatches <- nucleotideMismatches - offsets
  
  codonMismatches <- (nucleotideMismatches - 1) %/% 3
  
  list(nucleotideMismatches = nucleotideMismatches,
       codonMismatches = codonMismatches)
}