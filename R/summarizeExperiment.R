.hasReadComponent <- function(composition, component) {
    tmp <- gregexpr(pattern = component, text = composition, fixed = TRUE)[[1]]
    if (length(tmp) == 1 && tmp == -1) {
        FALSE
    } else {
        TRUE
    }
}

#' Summarize and collapse multiple mutational scanning experiments
#'
#' Combine multiple sequence lists (as returned by \code{\link{digestFastqs}}
#' into a \code{\link[SummarizedExperiment]{SummarizedExperiment}}, with
#' observed variable sequences (sequence pairs) in rows and samples in columns.
#'
#' @param x A named list of objects returned by \code{\link{digestFastqs}}.
#'     Names are used to link the objects to the metadata provided in
#'     \code{coldata}.
#' @param coldata A \code{data.frame} with at least one column "Name", which
#'     will be used to link to objects in \code{x}. A potentially subset and
#'     reordered version of \code{coldata} is stored in the \code{colData} of the
#'     returned \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#' @param countType Either "reads" or "umis". If "reads", the "count" assay of
#'     the returned object will contain the observed number of reads for each
#'     sequence (pair). If "umis", the "count" assay will contain the number of
#'     unique UMIs observed for each sequence (pair).
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} \code{x}
#'     with
#'     \describe{
#'         \item{assays(x)$counts}{containing the observed number of sequences or
#'         sequence pairs (if \code{countType} = "reads"), or the observed number of
#'         unique UMIs for each sequence or sequence pair (if \code{countType} =
#'         "umis").}
#'         \item{rowData(x)}{containing the unique sequences or sequence pairs.}
#'         \item{colData(x)}{containing the metadata provided by \code{coldata}.}
#'     }
#'
#' @author Michael Stadler, Charlotte Soneson
#'
#' @export
#'
#' @importFrom SummarizedExperiment SummarizedExperiment colData
#' @importFrom BiocGenerics paste
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges IntegerList
#' @importFrom methods is new as
#' @importFrom dplyr bind_rows distinct left_join %>% mutate
#' @importFrom tidyr unite
#' @importFrom rlang .data
#'
summarizeExperiment <- function(x, coldata, countType = "umis") {
    ## --------------------------------------------------------------------------
    ## Pre-flight checks
    ## --------------------------------------------------------------------------
    if (!is(x, "list") || is.null(names(x))) {
        stop("'x' must be a named list")
    }
    if (any(duplicated(names(x)))) {
        stop("duplicated names in 'x' (e.g. technical replicated to be merged) is not supported yet")
    }
    if (!is(coldata, "data.frame") || !("Name" %in% colnames(coldata))) {
        stop("'coldata' must be a data.frame with at least one column, named ",
             "'Name'.")
    }
    .assertScalar(x = countType, type = "character",
                  validValues = c("umis", "reads"))
    ## If no UMI sequences were given, then countType = "umis" should not be allowed
    if (countType == "umis" &&
        !(all(sapply(x, function(w) .hasReadComponent(w$parameters$elementsForward, "U") ||
                     .hasReadComponent(w$parameters$elementsReverse, "U"))))) {
        stop("'countType' is set to 'umis', but no UMI sequences ",
             "were provided when quantifying. ",
             "Set 'countType' to 'reads' instead.")
    }
    ## Get the mutNameDelimiters, and make sure that they are the same
    mutnamedel <- unique(sapply(x, function(w) w$parameters$mutNameDelimiter))
    if (length(mutnamedel) > 1) {
        stop("All samples must have the same 'mutNameDelimiter'")
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
    coldata <- coldata[match(nms, coldata$Name), , drop = FALSE]

    ## --------------------------------------------------------------------------
    ## Get all sequences, and all sample names
    ## --------------------------------------------------------------------------
    allSequences <- Reduce(function(...) dplyr::full_join(..., by = "mutantName"),
                           lapply(x, function(w) w$summaryTable[, c("mutantName", "sequence")])) %>%
        replace(., is.na(.), "") %>%
        tidyr::unite(sequence, -.data$mutantName, sep = ",") %>%
        dplyr::mutate(sequence = gsub("[,]+", ",", sequence)) %>%
        dplyr::mutate(sequence = sub(",$", "", sequence)) %>%
        dplyr::mutate(sequence = sub("^,", "", sequence)) %>%
        dplyr::mutate(sequence = vapply(sequence, function(x) {
            paste(unique(strsplit(x, ",")[[1]]), collapse = ",")
        }, ""))
    allSequences <- S4Vectors::DataFrame(allSequences)
    allSamples <- names(x)
    names(allSamples) <- allSamples

    ## Add info about nbr mutated bases/codons/AAs
    for (v in c("nbrMutBases", "nbrMutCodons", "nbrMutAAs")) {
        tmp <- do.call(dplyr::bind_rows,
                       lapply(x, function(w)
                           w$summaryTable[, c("mutantName", v)])) %>%
            dplyr::distinct()
        tmp <- methods::as(split(tmp[[v]], f = tmp$mutantName),
                           "IntegerList")
        allSequences[[v]] <- tmp[allSequences$mutantName]
        allSequences[[paste0("min", sub("^n", "N", v))]] <- min(allSequences[[v]])
        allSequences[[paste0("max", sub("^n", "N", v))]] <- max(allSequences[[v]])
    }

    ## varLengths, mutantNameAA, mutationTypes
    for (v in c("varLengths", "mutantNameAA", "mutationTypes")) {
        tmp <- do.call(dplyr::bind_rows, 
                       lapply(x, function(w) 
                           w$summaryTable[, c("mutantName", v)])) %>%
            dplyr::distinct()
        tmp <- methods::as(split(tmp[[v]], f = tmp$mutantName),
                           "CharacterList")
        allSequences[[v]] <- tmp[allSequences$mutantName]
    }

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
        colData = coldata[match(allSamples, coldata$Name), , drop = FALSE],
        rowData = allSequences,
        metadata = list(parameters = lapply(allSamples, function(w) x[[w]]$parameters),
                        countType = countType,
                        mutNameDelimiter = mutnamedel)
    )

    if (!any(allSequences$mutantName == "")) {
        rownames(se) <- allSequences$mutantName
    }
    colnames(se) <- allSamples

    return(se)
}