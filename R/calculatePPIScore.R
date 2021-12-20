#' Calculate protein-protein interaction scores (PPI).
#' 
#' Using sequence counts before and after selection, calculate protein-protein
#' interaction scores as described in Diss et al. 2018.
#' 
#' @param se SummarizedExperiment object as returned by \code{\link{digestFastqs}}.
#' @param pairingCol Name of column in \code{colData(se)} with replicate/pairing
#'   information. Samples with the same value in this column will be paired.
#' @param ODCols Name(s) of column(s) in \code{colData(se)} with OD values (numeric),
#'   used to normalize for different numbers of cells.
#' @param comparison 3-element character vector of the form \code{(column, numerator,
#'   denominator)}. \code{column} is the name of the column in \code{colData(se)}
#'   with experimental conditions. \code{numerator} and \code{denominator} define
#'   the comparison, e.g. \code{c("cond", "output", "input")} will look in the
#'   \code{"cond"} column and calculate PPI for the ratio of \code{"output"} over
#'   \code{"input"} counts.
#' @param WTrows Vector of row names that will be used as the reference when
#'   calculating PPI scores. If more than one value is provided, the average of
#'   the corresponding PPI scores is used as a reference. If NULL, no division
#'   by WT scores will be done.
#'   
#' @return A numeric vector with PPI scores.
#' 
#' @author Michael Stadler and Charlotte Soneson
#' 
#' @references "The genetic landscape of a physical interaction."
#'   Diss G and Lehner B. Elife. 2018;7:e32472. doi: 10.7554/eLife.32472.
#' 
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData assay
#' @importFrom Matrix colSums
#' 
#' @export
calculatePPIScore <- function(se, pairingCol, ODCols, comparison, WTrows) {
    ## ------------------------------------------------------------------------
    ## pre-flight checks
    ## ------------------------------------------------------------------------
    .assertVector(x = se, type = "SummarizedExperiment")
    
    ## pairingCol is in colData(se)
    .assertScalar(x = pairingCol, type = "character", 
                  validValues = colnames(SummarizedExperiment::colData(se)))

    ## ODCols are all in colData(se) and contain numeric values
    .assertVector(x = ODCols, type = "character", rngLen = c(1, Inf),
                  validValues = colnames(SummarizedExperiment::colData(se)))
    for (odc in ODCols) {
        .assertVector(x = SummarizedExperiment::colData(se)[[odc]], 
                      type = "numeric")
    }

    ## comparison is length(3)-character with column and values in colData(se)
    .assertVector(x = comparison, type = "character", len = 3)
    .assertScalar(x = comparison[1], type = "character", 
                  validValues = colnames(SummarizedExperiment::colData(se)))
    .assertVector(x = comparison[2:3], type = "character", 
                  validValues = SummarizedExperiment::colData(se)[[comparison[1]]])
    
    ## there is exactly one observation per pairing and condition
    if (any(table(colData(se)[colData(se)[, comparison[1]] %in% 
                              comparison[2:3], pairingCol],
                  colData(se)[colData(se)[, comparison[1]] %in% 
                              comparison[2:3], comparison[1]]) != 1)) {
        stop("There must be exactly one sample for every ", 
             "replicate-condition combination ",
             "defined by 'pairingCol' and 'comparison'")
    }
    
    ## WTrows exist in the SE
    .assertVector(x = WTrows, type = "character", validValues = rownames(se), 
                  allowNULL = TRUE)

    ## ------------------------------------------------------------------------
    ## subset se and reorder samples by shared replicates
    ## ------------------------------------------------------------------------
    se_numerator <- se[, colData(se)[, comparison[1]] == comparison[2]]
    se_denominator <- se[, colData(se)[, comparison[1]] == comparison[3]]
    
    shared_repl <- intersect(colData(se_numerator)[, pairingCol], 
                             colData(se_denominator)[, pairingCol])
    se_numerator <- se_numerator[, match(shared_repl, colData(se_numerator)[, pairingCol])]
    se_denominator <- se_denominator[, match(shared_repl, 
                                             colData(se_denominator)[, pairingCol])]
    
    ## ------------------------------------------------------------------------
    ## calculate normalized counts (n_i)
    ## ------------------------------------------------------------------------
    norm_counts_numerator <- sweep(
        as.matrix(assay(se_numerator, "counts")), MARGIN = 2, 
        STATS = apply(colData(se_numerator)[, ODCols, drop = FALSE], MARGIN = 1, prod)/
            Matrix::colSums(assay(se_numerator, "counts")), 
        FUN = "*")
    norm_counts_denominator <- sweep(
        as.matrix(assay(se_denominator, "counts")), MARGIN = 2, 
        STATS = apply(colData(se_denominator)[, ODCols, drop = FALSE], MARGIN = 1, prod)/
            Matrix::colSums(assay(se_denominator, "counts")), 
        FUN = "*")
    n <- log2(norm_counts_numerator/norm_counts_denominator)
    n[!is.finite(n)] <- NA
    colnames(n) <- paste0(comparison[2], "_vs_", comparison[3], "_repl", shared_repl)
    
    
    ## ------------------------------------------------------------------------
    ## calculate PPI = n_i / n_WT
    ## ------------------------------------------------------------------------
    if (is.null(WTrows)) {
        nWT <- rep(1, ncol(n))
    } else {
        ## geometric mean
        # tmp0 <- n[WTrows, , drop = FALSE]
        # tmp0 <- tmp0[apply(tmp0, 1, min) > 0, , drop = FALSE]
        # nWT <- apply(tmp0, 2, function(s) exp(mean(log(s))))
        ## arithmetic mean
        nWT <- colMeans(n[WTrows, , drop = FALSE])
    }
    PPI <- sweep(n, MARGIN = 2, STATS = nWT, FUN = "/")
    return(PPI)
}