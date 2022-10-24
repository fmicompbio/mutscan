#' Visualize the filtering procedure 
#' 
#' Display the number (or fraction) of reads remaining after each step 
#' of the internal \code{mutscan} filtering. 
#' 
#' The function assumes that the number of reads filtered out in each step 
#' are provided as columns of \code{colData(se)}, with column names 
#' of the form \code{f[0-9]_filteringreason}, and that all filtering columns 
#' occur between the columns named \code{nbrTotal} and \code{nbrRetained}. 
#' 
#' @author Charlotte Soneson
#' 
#' @export
#' 
#' @return A \code{ggplot} object. 
#' 
#' @param se A \code{SummarizedExperiment} object, e.g. from 
#'     \code{summarizeExperiment}.
#' @param valueType Either "reads" or "fractions", indicating whether to 
#'     plot the number of reads, or the fraction of the total number of reads,
#'     that are retained after/filtered out in each filtering step.
#' @param onlyActiveFilters Logical scalar, whether to only include the 
#'     active filters (i.e., where any read was filtered out in any of the samples). 
#'     Defaults to \code{TRUE}. 
#' @param displayNumbers Logical scalar, indicating whether to display the 
#'     number (or fraction) of reads retained at every filtering step. 
#' @param numberSize Numeric scalar, indicating the size of the displayed 
#'     numbers (if \code{displayNumbers} is \code{TRUE}). 
#' @param plotType Character scalar, indicating what to show in the plot. 
#'     Either \code{"remaining"} or \code{"filtered"}.
#' @param facetBy Character scalar, indicating the variable by which the plots 
#'     should be facetted. Either \code{"sample"} or \code{"step"}.
#' 
#' @importFrom tidyselect matches
#' @importFrom dplyr select %>% mutate group_by summarize pull ungroup filter
#' @importFrom tibble rownames_to_column 
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap theme theme_bw labs
#'     geom_text element_text
#' @importFrom rlang .data
#' 
#' @examples 
#' se <- readRDS(system.file("extdata", "GSE102901_cis_se.rds", 
#'                           package = "mutscan"))[1:200, ]
#' plotFiltering(se)
#' 
plotFiltering <- function(se, valueType = "reads", onlyActiveFilters = TRUE,
                          displayNumbers = TRUE, numberSize = 4,
                          plotType = "remaining", facetBy = "sample") {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = valueType, type = "character", 
                  validValues = c("reads", "fractions"))
    .assertScalar(x = onlyActiveFilters, type = "logical")
    .assertScalar(x = displayNumbers, type = "logical")
    .assertScalar(x = numberSize, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = plotType, type = "character",
                  validValues = c("remaining", "filtered"))
    .assertScalar(x = facetBy, type = "character", 
                  validValues = c("sample", "step"))

    ## ----------------------------------------------------------------------- ##
    ## Extract relevant columns from colData(se)
    ## ----------------------------------------------------------------------- ##
    cd <- as.data.frame(colData(se)) %>%
        dplyr::select("nbrTotal":"nbrRetained") %>%
        dplyr::select(c("nbrTotal", tidyselect::matches("^f.*_"), "nbrRetained"))
    
    ## Remove inactive filters if desired
    if (onlyActiveFilters) {
        cd <- cd[, colSums(cd) > 0, drop = FALSE]
    }
    
    ## Check that filtering columns are in the right order
    nbrs <- gsub("^f", "", sapply(strsplit(colnames(cd), "_"), .subset, 1))
    nbrs <- nbrs[-c(1, length(nbrs))]
    nbrs <- gsub("a$|b$", "", nbrs)
    stopifnot(all(diff(as.numeric(nbrs)) >= 0))
    
    ## Check that all columns are numeric
    stopifnot(all(apply(cd, 2, is.numeric)))
    
    ## ----------------------------------------------------------------------- ##
    ## Reshape into long format
    ## ----------------------------------------------------------------------- ##
    cd <- cd %>% 
        tibble::rownames_to_column("sample") %>%
        tidyr::gather(key = "step", value = "nbrReads", -"sample") %>%
        dplyr::mutate(step = factor(.data$step, levels = colnames(cd)))
    
    ## Check that numbers add up (total - all filters = retained)
    stopifnot(all(cd %>% 
                      dplyr::group_by(sample) %>%
                      dplyr::summarize(remdiff = .data$nbrReads[.data$step == "nbrTotal"] - 
                                           sum(.data$nbrReads[grepl("^f", .data$step)]),
                                       remlist = .data$nbrReads[.data$step == "nbrRetained"]) %>%
                      dplyr::mutate(obsdiff = .data$remdiff - .data$remlist) %>%
                      dplyr::pull(.data$obsdiff) == 0)) 
    
    ## ----------------------------------------------------------------------- ##
    ## Calculate number of remaining reads at each step
    ## ----------------------------------------------------------------------- ##
    cd <- cd %>%
        dplyr::group_by(sample) %>%
        dplyr::mutate(fracReads = signif(.data$nbrReads / 
                                             .data$nbrReads[.data$step == "nbrTotal"],
                                        digits = 3)) %>%
        dplyr::mutate(cumulsum = sapply(seq_along(.data$step), function(i) {
            sum(.data$nbrReads[as.numeric(.data$step) <= as.numeric(.data$step[i]) & 
                                   .data$step != "nbrTotal"])
        })) %>%
        dplyr::mutate(nbrRemaining = .data$nbrReads[.data$step == "nbrTotal"] - 
                          .data$cumulsum) %>%
        dplyr::mutate(fracRemaining = signif(.data$nbrRemaining /
                                                 .data$nbrReads[.data$step == "nbrTotal"],
                                            digits = 3)) %>%
        dplyr::filter(.data$step != "nbrRetained") %>%
        dplyr::ungroup()
    
    ## ----------------------------------------------------------------------- ##
    ## Create plot
    ## ----------------------------------------------------------------------- ##
    
    yvar <- switch(
        paste0(valueType, "_", plotType),
        reads_remaining = "nbrRemaining",
        fractions_remaining = "fracRemaining",
        reads_filtered = "nbrReads",
        fractions_filtered = "fracReads"
    )
    ylab1 <- switch(
        valueType,
        reads = "Number of reads",
        fractions = "Fraction of reads"
    )
    ylab2 <- switch(
        plotType,
        remaining = " remaining after",
        filtered = " filtered out in"
    )
    if (plotType == "filtered") {
        cd <- cd %>% 
            dplyr::filter(.data$step != "nbrTotal")
    }
    gg <- ggplot2::ggplot(cd, ggplot2::aes(x = .data[[setdiff(c("step", "sample"), facetBy)]], 
                                           y = .data[[yvar]], 
                                           label = .data[[yvar]])) + 
        ggplot2::geom_bar(stat = "identity") + 
        ggplot2::facet_wrap(~ .data[[facetBy]], ncol = 1, scales = "free_y") + 
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, 
                                                           vjust = 0.5, size = 12),
                       axis.text.y = ggplot2::element_text(size = 12),
                       axis.title = ggplot2::element_text(size = 14)) + 
        ggplot2::labs(x = "", y = ylab1, 
                      title = paste0(ylab1, ylab2, " each filtering step"))
    
    if (displayNumbers) {
        gg <- gg + ggplot2::geom_text(vjust = 1.5, color = "white", size = numberSize)
    }
    gg
}