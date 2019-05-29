#' Summarize and collapse multiple mutational scanning experiments
#' 
#' Combine multiple sequence lists (as returned by \code{\link{digestFastqs}}
#' into a \code{\link[SummarizedExperiment]{SummarizedExperiment}}, with
#' observed variable sequences (sequence pairs) in rows and samples in columns.
#' 
#' @param x A named list of SummarizedExperiment objects as returned by
#'   \code{\link{digestFastqs}}. Names are used to link the objects to the
#'   metadata provided in \code{coldata}.
#' @param meta A \code{data.frame} with at least one column ("Name"). "Name"
#'   will be used to link to objects in \code{x}, and a potentially subset and
#'   reordered version of \code{coldata} is stored in \code{colData} or the
#'   returned \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#' @param countType Either "reads" or "umis". If "reads", the "count" assay of
#'   the returned object will contain the observed number of reads for each
#'   sequence (pair). If "umis", the "count" assay will contain the number of
#'   unique UMIs observed for each sequence (pair).
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} with
#'   \describe{
#'     \item{assays(x)$counts}{containing the observed number of sequences or
#'     sequence pairs (if \code{countType} = "reads"), or the observed number of
#'     unique UMIs for each sequence or sequence pair (if \code{countType} =
#'     "umis").}
#'     \item{rowData(x)}{containing the unique sequences or sequence pairs.}
#'     \item{colData(x)}{containing the metadata provided by \code{coldata}.}
#'   }
#'
#' @importFrom SummarizedExperiment SummarizedExperiment colData
#' @importFrom BiocGenerics paste
#' @importFrom S4Vectors DataFrame
#' @importFrom methods is new
#' @importFrom dplyr bind_rows distinct left_join
#' 
#' @author Michael Stadler
#' 
#' @export
summarizeExperiment <- function(x, coldata, countType = "umis") {
  ## --------------------------------------------------------------------------
  ## Pre-flight checks
  ## --------------------------------------------------------------------------
  if (!is(x, "list") || is.null(names(x)) ||
      length(unique(vapply(x, function(w) w$experimentType, ""))) != 1L) {
    stop("'x' must be a named list with either only 'cis' or only 'trans' objects")
  }
  if (any(duplicated(names(x)))) {
    stop("duplicated names in 'x' (e.g. technical replicated to be merged) is not supported yet")
  }
  if (!is(coldata, "data.frame") || !("Name" %in% colnames(coldata))) {
    stop("'coldata' must be a data.frame with at least one column named ",
         "'Name'.")
  }
  if (!is.character(countType) || length(countType) != 1 || 
      !(countType %in% c("umis", "reads"))) {
    stop("'countType' must be either 'umis' or 'reads'")
  }
  coldata$Name <- as.character(coldata$Name)

  ## --------------------------------------------------------------------------
  ## Link elements in x with coldata
  ## --------------------------------------------------------------------------
  nms <- intersect(names(x), coldata$Name)
  if (length(nms) == 0) {
    stop("names in 'x' do not match 'coldata$Name'")
  } else if (length(nms) < length(x)) {
    nmsNotUsed <- setdiff(names(x), nms)
    warning("skipping ", length(nmsNotUsed), " elements in 'x' (",
            paste(nmsNotUsed, collapse = ", "), ") because they did not ",
            "match any 'coldata$Name'.")
  }
  x <- x[nms]
  coldata <- coldata[match(nms, coldata$Name), ]
  
  ## --------------------------------------------------------------------------
  ## Get all sequences, and all sample names
  ## --------------------------------------------------------------------------
  allSequences <- do.call(
    dplyr::bind_rows, 
    lapply(x, function(w) w$summaryTable[, c("sequence", "mutantName")])) %>%
    dplyr::distinct()
  allSamples <- names(x)
  names(allSamples) <- allSamples
  
  ## --------------------------------------------------------------------------
  ## Create a sparse matrix
  ## --------------------------------------------------------------------------
  countCol <- ifelse(countType == "umis", "nbrUmis", "nbrReads")
  tmp <- do.call(dplyr::bind_rows, lapply(allSamples, function(s) {
    st <- x[[s]]$summaryTable
    data.frame(i = match(st$mutantName, allSequences$mutantName),
               j = match(s, allSamples),
               x = as.numeric(st[, countCol]))
  }))
  countMat <- methods::new("dgTMatrix", i = tmp$i - 1L, j = tmp$j - 1L, 
                           x = tmp$x, Dim = c(nrow(allSequences), length(allSamples)))
  
  ## --------------------------------------------------------------------------
  ## Create the colData
  ## --------------------------------------------------------------------------
  addMeta <- do.call(dplyr::bind_rows, lapply(allSamples, function(s) {
    x[[s]]$filterSummary %>% dplyr::mutate(Name = s)
  }))
  coldata <- dplyr::left_join(coldata, addMeta, by = "Name")
  
  ## --------------------------------------------------------------------------
  ## Create SummarizedExperiment object
  ## --------------------------------------------------------------------------
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = countMat),
    colData = coldata[match(allSamples, coldata$Name), ],
    rowData = S4Vectors::DataFrame(sequence = allSequences$sequence),
    metadata = lapply(allSamples, function(w) x[[w]]$parameters)
  )
  
  if (!any(allSequences$mutantName == "")) {
    rownames(se) <- allSequences$mutantName
  } 
  colnames(se) <- allSamples
  
  return(se)
}