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
#'   \code{summarizeExperiment}.
#' @param valueType Either "reads" or "fractions", indicating whether to 
#'   plot the number of reads, or the fraction of the total number of reads,
#'   that are retained after each filtering step.
#' @param onlyActiveFilters Logical scalar, whether to only include the 
#'   active filters (i.e., where any read was filtered out in any of the samples). 
#'   Defaults to \code{TRUE}. 
#' @param displayNumbers Logical scalar, indicating whether to display the 
#'   number (or fraction) of reads retained at every filtering step. 
#' @param numberSize Numeric scalar, indicating the size of the displayed 
#'   numbers (if \code{displayNumbers} is \code{TRUE}). 
#' 
#' @importFrom tidyselect matches
#' @importFrom dplyr select %>% mutate group_by summarize pull ungroup filter
#' @importFrom tibble rownames_to_column 
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap theme theme_bw labs
#'   geom_text element_text
#' @importFrom rlang .data
#' 
plotFiltering <- function(se, valueType = "reads", onlyActiveFilters = TRUE,
                          displayNumbers = TRUE, numberSize = 4) {
  stopifnot(is(se, "SummarizedExperiment"))
  stopifnot(is.character(valueType) && length(valueType) == 1 && 
              valueType %in% c("reads", "fractions"))
  stopifnot(is.logical(onlyActiveFilters) && length(onlyActiveFilters) == 1)
  stopifnot(is.logical(displayNumbers) && length(displayNumbers) == 1)
  stopifnot(is.numeric(numberSize) && length(numberSize) == 1 && 
              numberSize > 0)
  
  ## ----------------------------------------------------------------------- ##
  ## Extract relevant columns from colData(se)
  ## ----------------------------------------------------------------------- ##
  cd <- as.data.frame(colData(se)) %>%
    dplyr::select(.data$nbrTotal:.data$nbrRetained) %>%
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
    tidyr::gather(key = "step", value = "nbrReads", -.data$sample) %>%
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
    dplyr::mutate(cumulsum = sapply(seq_along(.data$step), function(i) {
      sum(.data$nbrReads[as.numeric(.data$step) <= as.numeric(.data$step[i]) & 
                           .data$step != "nbrTotal"])
    })) %>%
    dplyr::mutate(nbrRemaining = .data$nbrReads[.data$step == "nbrTotal"] - 
                    .data$cumulsum) %>%
    dplyr::mutate(fracRemaining = round(.data$nbrRemaining /
                                          .data$nbrReads[.data$step == "nbrTotal"],
                                        digits = 3)) %>%
    dplyr::filter(.data$step != "nbrRetained") %>%
    dplyr::ungroup()
  
  ## ----------------------------------------------------------------------- ##
  ## Create plot
  ## ----------------------------------------------------------------------- ##
  yvar <- switch(
    valueType,
    reads = "nbrRemaining",
    fractions = "fracRemaining"
  )
  ylab <- switch(
    valueType,
    reads = "Number of reads",
    fractions = "Fraction of reads"
  )
  gg <- ggplot2::ggplot(cd, ggplot2::aes(x = .data$step, y = .data[[yvar]], 
                                         label = .data[[yvar]])) + 
    ggplot2::geom_bar(stat = "identity") + 
    ggplot2::facet_wrap(~ sample, ncol = 1, scales = "free_y") + 
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, 
                                                       vjust = 0.5)) + 
    ggplot2::labs(x = "", y = ylab, 
                  title = paste0(ylab, " remaining after each filtering step"))
  
  if (displayNumbers) {
    gg <- gg + ggplot2::geom_text(vjust = 1.5, color = "white", size = numberSize)
  }
  gg
}