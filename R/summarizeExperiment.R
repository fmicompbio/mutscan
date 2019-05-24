#' Summarize and collapse multiple mutational scanning experiments
#' 
#' Combine multiple sequence lists (as returned by \code{\link{readFastqs}} or
#' \code{\link{filterReads}}) into a \code{\link[SummarizedExperiment]{SummarizedExperiment}},
#' with observed variable sequences (sequence pairs) in rows and samples in
#' columns.
#' 
#' @param x A named list of SummarizedExperiment objects as returned by \code{\link{readFastqs}}
#'   or \code{\link{filterReads}}). Names are used to link the objects to the
#'   metadata provided in \code{meta}.
#' @param meta A \code{data.frame} with at least two columns ("Name" and "Condition").
#'   "Name" will be used to link to objects in \code{x}, and a potentially subset
#'   and reordered version of \code{meta} is stored in \code{colData} or the
#'   returned \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} with
#'   \describe{
#'     \item{assays(x)$counts}{containing the observed number of sequences or sequence pairs.}
#'     \item{rowData(x)}{containing the unique sequences or sequence pairs.}
#'     \item{colData(x)}{containing the metadata provided by \code{meta}.}
#'   }
#'
#' @importFrom SummarizedExperiment SummarizedExperiment colData
#' @importFrom BiocGenerics paste
#' @importFrom S4Vectors DataFrame
#' @importFrom methods is
#' 
#' @author Michael Stadler
#' 
#' @export
summarizeExperiment <- function(x, meta) {
  ## --------------------------------------------------------------------------
  ## Pre-flight checks
  ## --------------------------------------------------------------------------
  if (!is(x, "list") || is.null(names(x)) ||
      !all(sapply(x, isValidL)) ||
      length(unique(sapply(x, function(x) SummarizedExperiment::colData(x)[1, "experimentType"]))) != 1L) {
    stop("'x' must be a named list with either only 'cis' or only 'trans' objects")
  }
  if (any(duplicated(names(x)))) {
    stop("duplicated names in 'x' (e.g. technical replicated to be merged) is not supported yet")
  }
  if (!is(meta, "data.frame") || !all(c("Name","Condition") %in% colnames(meta))) {
    stop("'meta' must be a data.frame with at least two columns with names ",
         "'Name' and 'Condition'.")
  }
  experimentType <- SummarizedExperiment::colData(x[[1]])[1, "experimentType"]
  
  ## --------------------------------------------------------------------------
  ## Link elements in x with meta
  ## --------------------------------------------------------------------------
  nms <- intersect(names(x), meta$Name)
  if (length(nms) == 0) {
    stop("names in 'x' do not match 'meta$Name'")
  } else if (length(nms) < length(x)) {
    nmsNotUsed <- setdiff(names(x), nms)
    warning("skipping ", length(nmsNotUsed), " elements in 'x' (",
            paste(nmsNotUsed, collapse = ", "), ") because they did not ",
            "match any 'meta$Name'.")
  }
  x <- x[nms]
  meta <- meta[ match(nms, meta$Name), ]
  
  ## --------------------------------------------------------------------------
  ## Count observed reads (read pairs)
  ## --------------------------------------------------------------------------
  ## ... collapse raw counts, umi counts and names per unique sequence
  xCounts <- lapply(x, function(x1) {
    if (experimentType == "cis") {
      seq <- as.character(assay(x1, "variableSeqForward")$seq)
      umiL <- split(as.character(assay(x1, "umis")$seq), seq)
      cnms <- as.character(rowData(x1)$encodedMutatedCodonsForward)[match(names(umiL), seq)]
    } else {
      seq <- paste(assay(x1, "variableSeqForward")$seq, 
                   assay(x1, "variableSeqReverse")$seq, 
                   sep = "__")
      umiL <- split(as.character(assay(x1, "umis")$seq), seq)
      cnms <- paste(rowData(x1)$encodedMutatedCodonsForward,
                    rowData(x1)$encodedMutatedCodonsReverse,
                    sep = ".")[match(names(umiL), seq)]
    }
    list(umi = unlist(lapply(umiL, function(u) length(unique(u))), use.names = FALSE),
         cnt = lengths(umiL, use.names = FALSE),
         seq = names(umiL),
         cnms = cnms)
  })
  uniseqs <- unique(do.call(c, lapply(xCounts, "[[", "seq")))
  cntall <- umiall <- matrix(0L, nrow = length(uniseqs), ncol = length(nms), dimnames = list(NULL, nms))
  cnms <- rep(NA, length(uniseqs))
  for (nm in nms) {
    i <- match(xCounts[[nm]]$seq, uniseqs)
    cntall[i, nm] <- as.integer(xCounts[[nm]]$cnt)
    umiall[i, nm] <- as.integer(xCounts[[nm]]$umi)
    cnms[i] <- xCounts[[nm]]$cnms
  }
  rownames(cntall) <- rownames(umiall) <- cnms

  ## REMARK: this does not consider x[[nm]]minQualMutatedForward or x[[nm]]minQualMutatedReverse yet
  ##         it could be
  ##         1. propagated to rowData(se), combining values for multiple identical
  ##            reads with min, max, mean, median, ...
  ##         2. used directly to further filter reads (would need an additional argument)
  ##
  ## REMARK: this currently only "rescues" encodedMutatedCodonsForward/Reverse from rowData(x) to
  ##         rowData(se) - more could be added
  
  ## --------------------------------------------------------------------------
  ## Create SummarizedExperiment object
  ## --------------------------------------------------------------------------
  if (experimentType == "cis") {
    rd <- DataFrame(variableSeqForward = uniseqs,
                    variableSeqReverse = NULL)
  } else {
    rd <- DataFrame(variableSeqForward = unlist(lapply(strsplit(uniseqs, "__"), "[", 1), use.names = FALSE),
                    variableSeqReverse = unlist(lapply(strsplit(uniseqs, "__"), "[", 2), use.names = FALSE))
  }
  se <- SummarizedExperiment(assays = list(counts = cntall, umis = umiall),
                             rowData = rd,
                             colData = as(meta, "DataFrame"))
  return(se)
}