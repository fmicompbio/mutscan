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
      length(unique(sapply(x, function(w) w$experimentType))) != 1L) {
    stop("'x' must be a named list with either only 'cis' or only 'trans' objects")
  }
  if (any(duplicated(names(x)))) {
    stop("duplicated names in 'x' (e.g. technical replicated to be merged) is not supported yet")
  }
  if (!is(meta, "data.frame") || !all(c("Name", "Condition") %in% colnames(meta))) {
    stop("'meta' must be a data.frame with at least two columns with names ",
         "'Name' and 'Condition'.")
  }
  experimentType <- x[[1]]$experimentType
  
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
  
  allSequences <- do.call(
    dplyr::bind_rows, 
    lapply(x, function(w) w$summaryTable[, c("sequence", "mutantName")])) %>%
    dplyr::distinct()
  allSamples <- names(x)
  names(allSamples) <- allSamples
  
  ## Create a sparse matrix
  tmp <- do.call(rbind, lapply(allSamples, function(s) {
    st <- x[[s]]$summaryTable
    data.frame(i = match(st$mutantName, allSequences$mutantName),
               j = match(s, allSamples),
               x = as.numeric(st$nbrUmis))
  }))
  umiCounts <- new("dgTMatrix", i = tmp$i - 1L, j = tmp$j - 1L, 
                   x = tmp$x, Dim = c(nrow(allSequences), length(allSamples)))
  
  addMeta <- do.call(rbind, lapply(allSamples, function(s) {
    x[[s]]$filterSummary %>% dplyr::mutate(Name = s)
  }))
  meta <- dplyr::left_join(meta, addMeta, by = "Name")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = umiCounts),
    colData = meta[match(allSamples, meta$Name), ],
    rowData = DataFrame(sequence = allSequences$sequence),
    metadata = lapply(allSamples, function(w) x[[w]]$parameters)
  )
  
  rownames(se) <- allSequences$mutantName
  colnames(se) <- allSamples
  
  return(se)
}