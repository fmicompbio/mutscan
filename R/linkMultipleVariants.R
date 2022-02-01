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
#' reads), they will be merged and processed as a single unit.
#'
#' linkMultipleVariants will process the input in the following way:
#' \itemize{
#' \item First, run \code{digestFastqs} with the parameters provided
#'     in \code{combinedDigestParams}. Typically, this will be a
#'     "naive" counting run, where the frequencies of all observed
#'     variants are tabulated. Moreover, the variable sequences
#'     within the forward and reverse reads, respectively, will be
#'     concatenated into a single sequence.
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
#'     variable sequences (as specified by the element lengths), and
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
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @return A \code{tibble} with columns corresponding to each of the variable
#' sequences, and a column with the total observed read count for the
#' combination.
#'
#' @importFrom dplyr select rename group_by summarize across matches
#' @importFrom tidyr separate separate_rows
#' @importFrom rlang .data
#'
linkMultipleVariants <- function(combinedDigestParams = list(), ...) {

    ## Process additional arguments
    paramsSeparate <- list(...)

    ## --------------------------------------------------------------------- ##
    ## Initial checks
    ## --------------------------------------------------------------------- ##
    ## combinedDigestParams and each element of paramsSeparate must be named
    ## lists
    .assertVector(x = combinedDigestParams, type = "list")
    .assertVector(x = names(combinedDigestParams), type = "character",
                  allowNULL = FALSE)
    for (parms in paramsSeparate) {
        .assertVector(x = parms, type = "list")
        .assertVector(x = names(parms),
                      type = "character", allowNULL = FALSE)
    }

    ## Fill the parameter lists with default values for all arguments
    ## that are not provided, to avoid having to check for possible
    ## NULL values downstream
    defaults <- formals(digestFastqs)
    paramsSeparate <- lapply(paramsSeparate, function(parm) {
        c(parm, as.list(defaults[!(names(defaults) %in% names(parm))]))
    })
    combinedDigestParams <- c(combinedDigestParams,
                              as.list(defaults[!(names(defaults) %in%
                                                     names(combinedDigestParams))]))

    ## --------------------------------------------------------------------- ##
    ## Checks
    ## --------------------------------------------------------------------- ##
    ## Arguments that must be the same in the combined run and all the
    ## separate runs
    argsEqual <- c("mergeForwardReverse")
    for (ae in argsEqual) {
        for (parm in paramsSeparate) {
            if (parm[[ae]] != combinedDigestParams[[ae]]) {
                stop("The argument ", ae, " must be the same in the ",
                     "combined and all the separate runs")
            }
        }
    }

    ## Additional arguments that are recommended to be the same (e.g. to get
    ## the same filtering)
    argsRecEqual <- c("adapterForward", "adapterReverse", "primerForward",
                      "primerReverse", "constantForward", "constantForward",
                      "avePhredMinForward", "avePhredMinReverse",
                      "constantMaxDistForward", "constantMaxDistReverse")
    for (are in argsRecEqual) {
        for (parm in paramsSeparate) {
            if (parm[[are]] != combinedDigestParams[[are]]) {
                message("We recommend that the argument ", are,
                        " is set to the same value in the ",
                        "combined and all the separate runs")
            }
        }
    }

    ## Adding UMI counts after the processing may not be accurate - will
    ## always use read counts
    if (length(grep("U", paste0(combinedDigestParams$elementLengthsForward,
                                combinedDigestParams$elementLengthsReverse))) > 0) {
        warning("Aggregating UMI counts may not be accurate - will use ",
                "read counts instead")
    }

    ## Get number and lengths of forward and reverse variable sequences
    nbrFwdV <- sum(strsplit(combinedDigestParams$elementsForward, "")[[1]] == "V")
    nbrRevV <- sum(strsplit(combinedDigestParams$elementsReverse, "")[[1]] == "V")
    ## Check each of the separate runs
    lengthsV <- unlist(lapply(paramsSeparate, function(parm) {
        if (!parm$mergeForwardReverse) {
            ## Don't merge - there can be only one V in total, and the length
            ## must be given
            els <- c(strsplit(parm$elementsForward, "")[[1]],
                     strsplit(parm$elementsReverse, "")[[1]])
            ell <- c(parm$elementLengthsForward, parm$elementLengthsReverse)
            if (sum(els == "V") != 1) {
                stop("Must have exactly one V in each separate run")
            }
            if (any(ell[which(els == "V")] == -1)) {
                stop("The length of the variable sequence must be explicitly ",
                     "specified in the separate runs")
            }
            ell[which(els == "V")]
        } else {
            ## Merge - there must be a V in both fwd and reverse, with the
            ## same length
            elsfwd <- strsplit(parm$elementsForward, "")[[1]]
            elsrev <- strsplit(parm$elementsReverse, "")[[1]]
            ellfwd <- parm$elementLengthsForward
            ellrev <- parm$elementLengthsReverse
            if (sum(elsfwd == "V") != 1 || sum(elsrev == "V") != 1) {
                stop("Must have exactly one V in the forward and one V in ",
                     "the reverse read in each separate run")
            }
            if ((ellfwd[which(elsfwd == "V")] != ellrev[which(elsrev == "V")]) ||
                any(c(ellfwd[which(elsfwd == "V")],
                      ellrev[which(elsrev == "V")]) == -1)) {
                stop("The length of the variable sequence must be explicitly ",
                     "specified and identical between the forward and reverse ",
                     "reads in the separate runs")
            }
        }
    }))

    ## Check that the number of Vs in the combined run agrees with the number of
    ## separate runs
    if (!combinedDigestParams$mergeForwardReverse) {
        if (nbrFwdV + nbrRevV != length(lengthsV)) {
            stop("The number of variable sequences in the combined run does not ",
                 "correspond to the number of separate runs")
        }
    } else {
        if ((nbrFwdV != nbrRevV) || (nbrFwdV != length(lengthsV))) {
            stop("The number of variable sequences in the combined run does not ",
                 "correspond to the number of separate runs")
        }
    }

    ## --------------------------------------------------------------------- ##
    ## Combined quantification - just tabulate combined variable sequences
    ## --------------------------------------------------------------------- ##
    outCombined <- do.call(digestFastqs, combinedDigestParams)

    ## Get count matrix with "raw" (uncorrected) sequences
    countCombined <- outCombined$summaryTable %>%
        dplyr::select(.data$sequence, .data$nbrReads)

    ## If applicable, separate into forward and reverse sequences
    ## (separated by _). Also check that all sequences are the same (and
    ## correct) length
    if (nbrFwdV > 0 && nbrRevV > 0 && !combinedDigestParams$mergeForwardReverse) {
        countCombined <- countCombined %>%
            tidyr::separate(sequence, into = c("forward", "reverse"), sep = "_")
        if (!all(nchar(countCombined$forward) ==
                 sum(lengthsV[seq_len(nbrFwdV)]))) {
            stop("Forward sequences have the wrong length")
        }
        if (!all(nchar(countCombined$reverse) ==
                 sum(lengthsV[nbrFwdV + seq_len(nbrRevV)]))) {
            stop("Reverse sequences have the wrong length")
        }
    } else if (nbrFwdV > 0) {
        countCombined <- countCombined %>%
            dplyr::rename(forward = sequence)
        if (!all(nchar(countCombined$forward) ==
                 sum(lengthsV[seq_len(nbrFwdV)]))) {
            stop("Forward sequences have the wrong length")
        }
    } else if (nbrRevV > 0) {
        countCombined <- countCombined %>%
            dplyr::rename(reverse = sequence)
        if (!all(nchar(countCombined$reverse) ==
                 sum(lengthsV[nbrFwdV + seq_len(nbrRevV)]))) {
            stop("Reverse sequences have the wrong length")
        }
    }

    ## Split forward and reverse sequences, respectively, into the parts
    ## corresponding to the separate runs
    if ("forward" %in% colnames(countCombined)) {
        countCombined <- countCombined %>%
            tidyr::separate(forward, into = paste0("V", seq_len(nbrFwdV)),
                            sep = cumsum(lengthsV[seq_len(nbrFwdV - 1)]))
    } else {
        colnames(countCombined)[colnames(countCombined) == "forward"] <- "V1"
    }

    if ("reverse" %in% colnames(countCombined)) {
        countCombined <- countCombined %>%
            tidyr::separate(reverse, into = paste0("V", nbrFwdV + seq_len(nbrRevV)),
                            sep = cumsum(lengthsV[nbrFwdV + seq_len(nbrRevV - 1)]))
    } else {
        colnames(countCombined)[colnames(countCombined) == "reverse"] <-
            paste0("V", nbrFwdV + 1)
    }

    ## --------------------------------------------------------------------- ##
    ## Quantify each variable sequence separately
    ## --------------------------------------------------------------------- ##
    outSeparate <- lapply(paramsSeparate, function(ps) {
        do.call(digestFastqs, ps)
    })

    convSeparate <- lapply(outSeparate, function(out) {
        out$summaryTable %>%
            dplyr::select(.data$mutantName, .data$sequence) %>%
            tidyr::separate_rows(.data$sequence, sep = ",")
    })

    ## --------------------------------------------------------------------- ##
    ## Replace naive sequences with corrected ones
    ## --------------------------------------------------------------------- ##
    for (i in seq_along(convSeparate)) {
        countCombined[[paste0("V", i)]] <-
            convSeparate[[i]]$mutantName[match(countCombined[[paste0("V", i)]],
                                               convSeparate[[i]]$sequence)]
    }

    ## Filter out rows with NAs (the sequence was not retained)
    countCombined <- countCombined[rowSums(is.na(countCombined)) == 0, ,
                                   drop = FALSE]

    ## Aggregate counts
    countAggregated <- countCombined %>%
        dplyr::group_by(dplyr::across(dplyr::matches("^V[0-9]+"))) %>%
        dplyr::summarize(nbrReads = sum(.data$nbrReads), .groups = "drop")

    countAggregated
}
