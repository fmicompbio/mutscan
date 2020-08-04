# collapsing is done directly from C++, so no need to provide access from R
#
# #' @title Collapse similar sequences
# #' 
# #' @description Group sequences of equal length into sets of similar sequences,
# #'   based on Hamming distance. An example application is the collapsing of
# #'   unique molecular identifier variants that arise due to sequencing errors.
# #' 
# #' @param seqs A \code{character} vector with sequences to be grouped. The
# #'   sequences have to be of equal length.
# #' @param counts [optional] A \code{numeric} vector with observed counts for
# #'   each sequence in \code{seqs}. If provided, \code{seqs} will be sorted
# #'   decreasingly by \code{counts} (start collapsoing from the most frequent
# #'   sequence). If \code{NULL}, the sequences are collapsed in the order given.
# #' @param method A \code{character} scalar selecting the collapsing method.
# #'   One of: greedy_clustering.
# #' @param tol A \code{numeric} scalar defining the tolerance for similar
# #'   sequences. If \code{tol} is in [0, 1), it defines the maximal Hamming
# #'   distance in terms of a fraction of sequence length:
# #'   (\code{round(tol * nchar(seqs[1]))}). A value greater or equal to 1 is
# #'   rounded and directly used as the maximum allowed Hamming distance.
# #' @param verbose A \code{logical} scalar. If \code{TRUE}, report on progress on
# #'   the R console.
# #' 
# #' @details To speed-up the search for similar sequences, the function
# #'   stores the sequences in a BK tree (https://en.wikipedia.org/wiki/BK-tree),
# #'   which reduced the number of string distance calculations. For a general
# #'   discussion of this problem with further details on BK trees see
# #'   Liu, D. (2019).
# #'   
# #'   The following algorithms can be selected using the \code{method} argument:
# #'   \itemize{
# #'     \item greedy_clustering: staring from the first (most frequent)
# #'       sequence, group it with all similar sequences that are within
# #'       Hamming distance \code{tol}. Remove these sequences from \code{seqs}
# #'       and continue with the next frequent sequence, until \code{seqs} is
# #'       empty.
# #'   }
# #'
# #' @return A \code{list} of \code{character} vectors. List elements are the
# #'   grouped sequences, list names are the query sequences that defined the
# #'   group.
# #'
# #' @examples 
# #' seqs <- c("AAAA","AAAC", "TTTT", "TTTC")
# #' collapseSeqs(seqs, tol = 1)
# #' 
# #' @references Daniel Liu. "Algorithms for efficiently collapsing reads with
# #'   Unique Molecular Identifiers". PeerJ. 2019; 7: e8275.
# #'   doi: 10.7717/peerj.8275
# #'   PMID: 31871845
# #'   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6921982/
# #'
# #' @author Michael Stadler
# #' 
# #' @export
# collapseSeqs <- function(seqs,
#                          counts = NULL,
#                          method = c("greedy_clustering"),
#                          tol = 0.05,
#                          verbose = TRUE) {
#   ## pre-flight checks
#   method = match.arg(method, several.ok = FALSE)
#   stopifnot(exprs = {
#     is.character(seqs)
#     length(seqs) > 0L
#     all(nchar(seqs) == nchar(seqs[1]))
#     is.null(counts) || (is.numeric(counts) && length(counts) == length(seqs))
#     is.numeric(tol)
#     length(tol) == 1L
#     tol >= 0
#     is.logical(verbose)
#     length(verbose) == 1L
#   })
#   
#   ## calculate allowed number of mismatches
#   if (tol < 1.0) {
#     tol <- round(tol * nchar(seqs[1]))
#   } else {
#     tol <- round(tol)
#   }
#   if (verbose)
#     message("analyzing ", nchar(seqs[1]),
#             "-mers, will tolerate a Hamming distance up to ", tol)
#   
#   ## sort sequences
#   if (!is.null(counts)) {
#     if (is.unsorted(-counts)) {
#       if (verbose)
#         message("sorting 'seqs' according to 'counts'...", appendLF = FALSE)
#       o <- order(counts, decreasing = TRUE)
#       seqs <- seqs[o]
#       counts <- counts[o]
#       if (verbose)
#         message("done")
#     }
#   }
#   
#   ## store sequences in BK tree
#   if (verbose)
#     message("storing sequences in BK tree...", appendLF = FALSE)
#   n <- mutscan:::bk_new(seqs) # BKtree::items keeps order of seqs
#   if (verbose)
#     message("done (", n, " unique sequences added)")
#   
#   ## create collapse groups
#   groups <- list()
#   if (method == "greedy_clustering") {
#     if (verbose)
#       message("start grouping (greedy_clustering)...", appendLF = FALSE)
#     groups <- mutscan:::greedy_clustering(tol, verbose);
#   }
#   if (verbose)
#     message("done (", n, " sequences grouped into ", length(groups), " sets)")
# 
#   ## return results
#   return(groups)
# }
