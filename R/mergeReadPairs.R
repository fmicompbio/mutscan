#' Merge forward and reverse read pairs
#'
#' Generate a consensus sequence for each pair of forward and reverse reads, for
#' each position choosing the base (forward/reverse) with the highest base
#' quality. If both reads have equal base quality, the forward read base is
#' selected. The function assumes that all reads (forward/reverse, and within
#' each object) have the same length.
#'
#' @param readsForward,readsReverse
#'   \code{\link[Biostrings]{QualityScaledDNAStringSet}} objects containing
#'   matched forward and reverse reads.
#'
#' @author Charlotte Soneson
#'
#' @keywords internal 
#' 
#' @return A \code{\link[Biostrings]{QualityScaledDNAStringSet}} with consensus
#'   sequences/qualities
#'
#' @importFrom ShortRead FastqQuality
#' @importFrom Biostrings quality QualityScaledDNAStringSet replaceLetterAt
#'   PhredQuality
#' @importFrom seqTools ascii2char
#'   
mergeReadPairs <- function(readsForward, readsReverse) {
  ## --------------------------------------------------------------------------
  ## Find positions where reverse read has higher quality than forward read
  ## --------------------------------------------------------------------------
  replaceAt <- as(ShortRead::FastqQuality(Biostrings::quality(readsReverse)), "matrix") > 
    as(ShortRead::FastqQuality(Biostrings::quality(readsForward)), "matrix")
  replaceAtList <- as(lapply(seq_len(nrow(replaceAt)),
                             function(x) which(replaceAt[x, ])), "IntegerList")
  
  ## --------------------------------------------------------------------------
  ## Replace bases and qualities in forward read with those from reverse read
  ## --------------------------------------------------------------------------
  consensusSeq <- Biostrings::replaceLetterAt(x = readsForward, at = replaceAt, 
                                              letter = readsReverse[replaceAtList])
  qualsForward <- as(ShortRead::FastqQuality(Biostrings::quality(readsForward)), "matrix")
  qualsReverse <- as(ShortRead::FastqQuality(Biostrings::quality(readsReverse)), "matrix")
  qualsForward[replaceAt] <- qualsReverse[replaceAt]
  phredQuals <- apply(qualsForward, 1, function(v) seqTools::ascii2char(v + 33))
  consensusQuals <- Biostrings::PhredQuality(phredQuals)
  
  QualityScaledDNAStringSet(consensusSeq, consensusQuals)
  
}

