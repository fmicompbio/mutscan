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
#' @param WTrows Vector of row names that will be used as the reference when
#'   calculating logFCs and statistics. If more than one value is provided, the
#'   sum of the corresponding counts is used to generate offsets. If NULL,
#'   offsets will be defined as the effective library sizes (using TMM
#'   normalization factors).
#' @param selAssay Assay to select from \code{se} for the analysis.
#' @param pseudocount Pseudocount to add when calculating log-fold changes. 
#' @param method Either 'edgeR' or 'limma'. If set to 'limma', voom is used to
#'   transform the counts and estimate observation weights before applying
#'   limma. In this case, the results also contain the standard errors of the
#'   logFCs.
#' 
#' @author Charlotte Soneson, Michael Stadler
#' 
#' @export
#' 
#' @importFrom edgeR DGEList scaleOffset estimateDisp glmQLFit glmQLFTest
#'   topTags predFC topTags calcNormFactors
#' @importFrom SummarizedExperiment colData assay assayNames
#' @importFrom limma voom eBayes topTable lmFit
#' 
calculateRelativeFC <- function(se, design, coef, WTrows, selAssay = "counts", 
                                pseudocount = 1, method = "edgeR") {
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
      stop("The provided 'selAssay' not present in 'se'.")
    }
  }
  
  if (!all(WTrows %in% rownames(se))) {
    stop("'se' must have rows named '", paste(WTrows, collapse = ", "), "'.")
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
  if (is.null(WTrows)) {
    offsets <- colSums(dge$counts) * edgeR::calcNormFactors(dge$counts)
  } else {
    ## geometric mean
    tmp0 <- dge$counts[WTrows, , drop = FALSE]
    tmp0 <- tmp0[apply(tmp0, 1, min) > 0, , drop = FALSE]
    offsets <- apply(tmp0, 2, function(s) exp(mean(log(s))))
    ## sum
    # offsets <- colSums(dge$counts[WTrows, , drop = FALSE])
  }
  dge <- edgeR::scaleOffset(dge, log(offsets))
  
  if (method == "edgeR") {
    dge <- edgeR::estimateDisp(dge, design = design)
    fit <- edgeR::glmQLFit(dge, design = design)
    qlf <- edgeR::glmQLFTest(fit, coef = coef)
    tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
    if (length(coef) == 1) {
      predfc <- edgeR::predFC(dge, design = design, prior.count = pseudocount)
      tt$logFC_shrunk <- predfc[, coef]
    }
  } else if (method == "limma") {
    vm <- limma::voom(dge, design = design, lib.size = exp(dge$offset))
    fit <- limma::lmFit(vm, design = design)
    fit <- limma::eBayes(fit)
    tt <- limma::topTable(fit, coef = coef, confint = TRUE, n = Inf, sort.by = "none")
    if (length(coef) == 1) {
      tt$se.logFC <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
    }
  } else {
    stop("'method' must be either 'edgeR' or 'limma'")
  }
  tt
}
