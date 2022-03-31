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
#' @importFrom dplyr bind_rows distinct left_join %>% mutate filter group_by
#'     summarize
#' @importFrom tidyr unite separate_rows
#' @importFrom rlang .data
#' @importFrom stats setNames
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

    ## All sample names
    allSamples <- names(x)
    names(allSamples) <- allSamples
    
    ## ------------------------------------------------------------------------
    ## Get all observed sequences for each mutant name
    ## For a given sample, the same mutant name can correspond to multiple 
    ## sequences, separated by ,
    ## ------------------------------------------------------------------------
    tmpdf <- do.call(dplyr::bind_rows, lapply(x, function(w) w$summaryTable))
    
    allSequences <- mergeValues(tmpdf$mutantName, tmpdf$sequence) %>%
        stats::setNames(c("mutantName", "sequence"))
    allSequences <- S4Vectors::DataFrame(allSequences)
    
    ## ------------------------------------------------------------------------
    ## Add info about nbr mutated bases/codons/AAs
    ## Also here, each column can contain multiple values separated with ,
    ## (e.g. if variable sequences were collapsed to WT in digestFastqs)
    ## ------------------------------------------------------------------------
    for (v in c("nbrMutBases", "nbrMutCodons", "nbrMutAAs")) {
        tmp <- mergeValues(tmpdf$mutantName, tmpdf[[v]]) %>%
            stats::setNames(c("mutantName", v))
        tmp[[v]] <- methods::as(
            lapply(strsplit(tmp[[v]], ","), function(w) sort(as.integer(w))),
            "IntegerList")
        allSequences[[v]] <- tmp[[v]][match(allSequences$mutantName,
                                            tmp$mutantName)]
        allSequences[[paste0("min", sub("^n", "N", v))]] <- min(allSequences[[v]])
        allSequences[[paste0("max", sub("^n", "N", v))]] <- max(allSequences[[v]])
    }

    ## ------------------------------------------------------------------------
    ## Add info about sequenceAA, mutantNameAA, mutationTypes, varLengths
    ## ------------------------------------------------------------------------
    for (v in c("sequenceAA", "mutantNameAA", "mutationTypes")) {
        tmp <- mergeValues(tmpdf$mutantName, tmpdf[[v]]) %>%
            stats::setNames(c("mutantName", v))
        allSequences[[v]] <- tmp[[v]][match(allSequences$mutantName,
                                            tmp$mutantName)]
    }
    
    ## There is only one varLengths value per sample, by construction 
    ## This is only useful when there's no wildtype sequence, so we use the 
    ## representative sequence only
    tmp <- mergeValues(tmpdf$mutantName, tmpdf$varLengths) %>%
        stats::setNames(c("mutantName", "varLengths"))
    tmp$varLengths <- methods::as(
        lapply(strsplit(tmp$varLengths, ","), function(w) unique(w)),
        "CharacterList")
    if (any(lengths(tmp$varLengths) != 1)) {
        warning("There are sequences with multiple different values for ",
                "varLengths. The uniqueVarLengths column will contain a ",
                "randomly selected one.")
    }
    allSequences$varLengths <- tmp$varLengths[match(allSequences$mutantName,
                                                    tmp$mutantName)]
    allSequences$uniqueVarLengths <- vapply(allSequences$varLengths,
                                            "[", 1, FUN.VALUE = "")


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