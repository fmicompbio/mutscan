#' Calculate logFCs relative to WT using edgeR
#'
#' Calculate logFCs and associated p-values for a given comparison, using the
#' Negative Binomial quasi-likelihood framework provided by edgeR. The observed
#' counts for the WT variant will be used as offsets in the model.
#'
#' @param se SummarizedExperiment object.
#' @param design Design matrix. The rows of the design matrix must be in the
#'   same order as the columns in \code{se}.
#' @param coef Coefficient(s) to test with edgeR. 
#' @param selAssay Assay to select from \code{se} for the analysis.
#' @param pseudocount Pseudocount to add when calculating log-fold changes. 
#' 
#' @author Charlotte Soneson, Michael Stadler
#' 
#' @export
#' 
#' @importFrom edgeR DGEList scaleOffset estimateDisp glmQLFit glmQLFTest topTags predFC
#' @importFrom SummarizedExperiment colData assay assayNames
#' 
calculateRelativeFC <- function(se, design, coef, selAssay = "counts", pseudocount = 1) {
  if (!is(se, "SummarizedExperiment")) {
    stop("'se' must be a SummarizedExperiment object.")
  }
  
  if (!is.character(selAssay) || length(selAssay) != 1) {
    stop("'selAssay' must be a character scalar.")
  }

  if (!(selAssay %in% SummarizedExperiment::assayNames(se))) {
    if (is.null(SummarizedExperiment::assayNames(se)) && 
        length(SummarizedExperiment::assays(se)) == 1) {
      warning("No assayNames provided in 'se', but only one assay present - using that.")
    } else {
      stop("Provided 'selAssay' not present in 'se'.")
    }
  }
  
  if (!("WT" %in% rownames(se))) {
    stop("'se' must have a row named 'WT'.")
  }
  
  if (nrow(design) != ncol(se)) {
    stop("The number of rows in 'design' (", nrow(design), ") is not equal to the number", 
         " of columns in 'se' (", ncol(se), ")")
  }
  
  if (!is.numeric(pseudocount) || length(pseudocount) != 1 || pseudocount < 0) {
    stop("'pseudocount' must be a non-negative scalar value.")
  }
  
  dge <- edgeR::DGEList(counts = as.matrix(SummarizedExperiment::assay(se, selAssay)),
                        samples = SummarizedExperiment::colData(se))
  offsets <- dge$counts["WT", ]
  dge <- edgeR::scaleOffset(dge, log(offsets))
  dge <- edgeR::estimateDisp(dge, design = design)
  fit <- edgeR::glmQLFit(dge, design = design)
  qlf <- edgeR::glmQLFTest(fit, coef = coef)
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  if (length(coef) == 1) {
    predfc <- edgeR::predFC(dge, design = design, prior.count = pseudocount)
    tt$logFC_shrunk <- predfc[, coef]
  }
  tt
}
