#' Process an experiment with multiple variable sequences
#'
#' This function enables the processing of data sets with multiple
#' variable sequences, which should potentially be handled in different
#' ways. For example, a barcode association experiment
#' with two variable sequences (the barcode and the biological variant)
#' that need to be processed differently, e.g. in terms of matching to
#' wildtype sequences or collapsing of similar sequences.
#' In contrast, while \code{digestFastqs} allow the specification
#' of multiple variable sequences (within each of the forward and reverse
#' reads), they will be concatenated and processed as a single unit.
#'
#' linkMultipleVariants will process the input in the following way:
#' \itemize{
#' \item First, run \code{digestFastqs} with the parameters provided
#'     in \code{combinedDigestParams}. Typically, this will be a
#'     "naive" counting run, where the frequencies of all observed
#'     variants are tabulated. The variable sequences
#'     within the forward and reverse reads, respectively, will be
#'     processed as a single sequence. 
#' \item Next, run \code{digestFastqs} with each of the additional
#'     parameter sets provided (\code{...}). Each of these should
#'     correspond to a single variable sequence from the combined
#'     run (i.e., if there are two Vs in the element specifications
#'     in the combined run, there should be two additional
#'     parameter sets provided, each corresponding to the
#'     processing of one variable sequence part). It is assumed
#'     that the order of the additional arguments correspond to the
#'     order of the variable sequences in the combined run, in such a way
#'     that if the variable sequences extracted in each of the separate
#'     runs are concatenated in the order that the parameter sets are
#'     provided to \code{linkMultipleVariants}, they will form the variable
#'     sequence extracted in the combined run.
#' \item The result of each of the separate runs is a 'conversion table',
#'     containing the final set of identified sequence variants as well
#'     as all individual sequences corresponding to each of them. This
#'     is then combined with the count table from the combined, "naive"
#'     run in order to create an aggregated count table. More precisely,
#'     each sequence in the combined run is split into the constituent
#'     variable sequences, and
#'     each variable sequence is then matched to the output from the right
#'     separate run, from which the final feature ID (mutant name, or
#'     collapsed sequence) will be extracted and used to replace the original
#'     sequence in the combined count table. Once all the matches are done,
#'     rows with NAs (where no match could be found in the separate run)
#'     are removed and the counts are aggregated across all identical
#'     combinations of variable sequences.
#' }
#'
#' In order to define the \code{elementsForward} and \code{elementsReverse}
#' arguments for the separate runs, a strategy that often works is to simply
#' copy the arguments from the combined run, and successively replace all
#' but one of the 'V's by 'S'. This will effectively process one variable
#' sequence at the time, while keeping all other elements of the reads
#' consistent (since this can affect e.g. filtering criteria). Note that
#' to process individual variable sequences in the reverse read, you also
#' need to swap the 'forward' and 'reverse' specifications (since
#' \code{digestFastqs} requires a forward read).
#'
#' @param combinedDigestParams A named list of arguments to
#'     \code{digestFastqs} for the combined ("naive") run.
#' @param ... Additional arguments providing arguments to \code{digestFastqs}
#'     for the separate runs (processing each variable sequence in turn).
#'     Each argument must be a named list of arguments to \code{digestFastqs}.
#' @param collapseToAA Either \code{TRUE}, \code{FALSE}, or a character vector
#'     indicating for which of the separate runs the sequences should be 
#'     collapsed to the amino acid mutant name rather than the codon- or 
#'     nucleotide-level name. \code{TRUE} is equivalent to providing the 
#'     names of all lists provided in \code{...}. 
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @return A list with the following elements: 
#' \itemize{
#' \item countAggregated - a \code{tibble} with columns corresponding to 
#'     each of the variable sequences, and a column with the total observed 
#'     read count for the combination.
#' \item convSeparate - a list of conversion tables from the respective 
#'     separate runs.
#' \item outCombined - the \code{digestFastqs} output for the combined run.
#' \item filtSeparate - a list of filtering tables for the separate runs.
#' }
#' 
#' @importFrom dplyr select rename group_by summarize across matches mutate
#' @importFrom tidyr separate separate_rows
#' @importFrom rlang .data
#'
linkMultipleVariants <- function(combinedDigestParams = list(), ..., 
                                 collapseToAA = FALSE) {

    ## Process additional arguments
    paramsSeparate <- list(...)

    ## --------------------------------------------------------------------- ##
    ## Initial checks
    ## --------------------------------------------------------------------- ##
    ## combinedDigestParams and each element of paramsSeparate must be named
    ## lists
    .assertVector(x = combinedDigestParams, type = "list")
    defaults <- formals(digestFastqs)
    .assertVector(x = names(combinedDigestParams), type = "character",
                  allowNULL = FALSE, validValues = names(defaults))
    for (parms in paramsSeparate) {
        .assertVector(x = parms, type = "list")
        .assertVector(x = names(parms), type = "character", allowNULL = FALSE,
                      validValues = names(defaults))
    }

    ## Fill the parameter lists with default values for all arguments
    ## that are not provided, to avoid having to check for possible
    ## NULL values downstream
    paramsSeparate <- lapply(paramsSeparate, function(parm) {
        c(parm, as.list(defaults[!(names(defaults) %in% names(parm))]))
    })
    combinedDigestParams <- c(combinedDigestParams,
                              as.list(defaults[!(names(defaults) %in%
                                                     names(combinedDigestParams))]))

    ## Make sure collapseToAA is a vector containing the names of the 
    ## runs where sequences should be collapsed to AAs
    if (length(collapseToAA) == 1) {
        if (is.logical(collapseToAA)) {
            if (collapseToAA) {
                collapseToAA <- names(paramsSeparate)
            } else {
                collapseToAA <- character(0)
            }
        } else {
            .assertVector(x = collapseToAA, type = "character",
                          validValues = names(paramsSeparate))
        }
    } else {
        .assertVector(x = collapseToAA, type = "character",
                      validValues = names(paramsSeparate))
    }
    
    ## --------------------------------------------------------------------- ##
    ## Checks
    ## --------------------------------------------------------------------- ##
    ## We can't have a fwd and a rev V in the same run unless they are merged
    ## since this would disagree with the order in the concatenated V
    ## Currently, each separate run can only have one V in the fwd and/or one 
    ## in the rev (both allowed if they are matched); otherwise the lengths
    ## extracted from the combined run don't match with those from the 
    ## separate runs
    for (parm in paramsSeparate) {
        nbrFwdV <- sum(strsplit(parm$elementsForward, "")[[1]] == "V")
        nbrRevV <- sum(strsplit(parm$elementsReverse, "")[[1]] == "V")
        if (nbrFwdV == 0 && nbrRevV == 0) {
            stop("Each separate run must contain at least one variable ",
                 "segment.")
        }
        if (nbrFwdV > 0 && nbrRevV > 0 && 
            !parm$mergeForwardReverse) {
            stop("A separate run can not have variable sequences in both ",
                 "the forward and reverse reads unless 'mergeForwardReverse' ",
                 "is TRUE.")
        }
        if (nbrFwdV > 1 || nbrRevV > 1) {
            stop("A separate run can have at most one variable segment ",
                 "in each of the forward and reverse reads, respectively.")
        }
    }
    
    ## Summing UMI counts after the processing may not be accurate - will
    ## always use read counts
    if (length(grep("U", paste0(combinedDigestParams$elementsForward,
                                combinedDigestParams$elementsReverse))) > 0) {
        warning("Aggregating UMI counts may not be accurate - will use ",
                "read counts instead")
    }

    ## --------------------------------------------------------------------- ##
    ## Combined quantification - just tabulate combined variable sequences
    ## --------------------------------------------------------------------- ##
    outCombined <- do.call(digestFastqs, combinedDigestParams)

    ## Get count matrix with "raw" (uncorrected) sequences
    countCombined <- outCombined$summaryTable %>%
        dplyr::select(.data$sequence, .data$nbrReads, .data$varLengths) %>%
        dplyr::mutate(idx = paste0("I", seq_along(.data$sequence)))

    ## If applicable, separate into forward and reverse sequences
    if (any(grepl("_", countCombined$sequence))) {
        countCombined <- countCombined %>%
            tidyr::separate(sequence, into = c("sequenceForward", 
                                               "sequenceReverse"), 
                            sep = "_") %>%
            tidyr::separate(.data$varLengths, into = c("varLengthsForward", 
                                                       "varLengthsReverse"),
                            sep = "_") %>%
            dplyr::mutate(
                nCompForward = vapply(strsplit(.data$varLengthsForward, ","), 
                                      length, 0),
                nCompReverse = vapply(strsplit(.data$varLengthsReverse, ","), 
                                      length, 0))
    } else {
        countCombined <- countCombined %>%
            dplyr::rename(sequenceForward = .data$sequence,
                          varLengthsForward = .data$varLengths) %>%
            dplyr::mutate(
                nCompForward = vapply(strsplit(.data$varLengthsForward, ","), 
                                      length, 0))
    }
    
    offsetForward <- 0
    ## Split forward and reverse sequences, respectively
    if ("sequenceForward" %in% colnames(countCombined)) {
        stopifnot("varLengthsForward" %in% colnames(countCombined))
        ## Check that all sequences have the same number of components
        stopifnot(length(unique(countCombined$nCompForward)) == 1)
        offsetForward <- countCombined$nCompForward[1]
        tmp <- split(countCombined, countCombined$varLengthsForward)
        countCombined <- unsplit(lapply(tmp, function(df) {
            w <- as.numeric(strsplit(df$varLengthsForward[1], ",")[[1]])
            df <- df %>%
                tidyr::separate(
                    .data$sequenceForward, 
                    into = names(paramsSeparate)[seq_along(w)],
                    # into = paste0("V", seq_along(w)),
                    sep = cumsum(w[seq_len(length(w) - 1)]))
        }), countCombined$varLengthsForward)
    }
    if ("sequenceReverse" %in% colnames(countCombined)) {
        stopifnot("varLengthsReverse" %in% colnames(countCombined))
        ## Check that all sequences have the same number of components
        stopifnot(length(unique(countCombined$nCompReverse)) == 1)
        tmp <- split(countCombined, countCombined$varLengthsReverse)
        countCombined <- unsplit(lapply(tmp, function(df) {
            w <- as.numeric(strsplit(df$varLengthsReverse[1], ",")[[1]])
            df <- df %>%
                tidyr::separate(
                    .data$sequenceReverse, 
                    into = names(paramsSeparate)[offsetForward + seq_along(w)],
                    # into = paste0("V", offsetForward + seq_along(w)),
                    sep = cumsum(w[seq_len(length(w) - 1)]))
        }), countCombined$varLengthsReverse)
    }
    
    ## --------------------------------------------------------------------- ##
    ## Quantify each variable sequence separately
    ## --------------------------------------------------------------------- ##
    outSeparate <- lapply(paramsSeparate, function(ps) {
        do.call(digestFastqs, ps)
    })

    ## Filter tables
    filtSeparate <- lapply(outSeparate, function(out) out$filterSummary)
    
    ## Conversion tables
    convSeparate <- lapply(outSeparate, function(out) {
        out$summaryTable %>%
            dplyr::select(.data$mutantName, .data$mutantNameAA,
                          .data$sequence) %>%
            tidyr::separate_rows(.data$sequence, sep = ",")
    })

    ## --------------------------------------------------------------------- ##
    ## Replace naive sequences with corrected ones
    ## --------------------------------------------------------------------- ##
    for (i in names(paramsSeparate)) {
        if (i %in% collapseToAA) {
            countCombined[[i]] <- 
                convSeparate[[i]]$mutantNameAA[match(countCombined[[i]],
                                                     convSeparate[[i]]$sequence)]
        } else {
            countCombined[[i]] <- 
                convSeparate[[i]]$mutantName[match(countCombined[[i]],
                                                   convSeparate[[i]]$sequence)]
        }
    }

    ## Filter out rows with NAs (the sequence was not retained)
    countCombined <- countCombined[rowSums(is.na(countCombined)) == 0, ,
                                   drop = FALSE]

    ## Aggregate counts
    countAggregated <- countCombined %>%
        dplyr::group_by(dplyr::across(names(paramsSeparate))) %>%
        dplyr::summarize(nbrReads = sum(.data$nbrReads), .groups = "drop")

    list(countAggregated = countAggregated, convSeparate = convSeparate,
         outCombined = outCombined, filtSeparate = filtSeparate)
}
