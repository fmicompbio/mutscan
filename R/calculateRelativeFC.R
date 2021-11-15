#' Calculate logFCs relative to WT using edgeR
#'
#' Calculate logFCs and associated p-values for a given comparison, using the
#' Negative Binomial quasi-likelihood framework provided by edgeR. The observed
#' counts for the WT variant will be used as offsets in the model.
#'
#' @param se SummarizedExperiment object.
#' @param design Design matrix. The rows of the design matrix must be in the
#'     same order as the columns in \code{se}.
#' @param coef Coefficient(s) to test with edgeR or limma. 
#' @param contrast Numeric contrast to test with edgeR or limma. 
#' @param WTrows Vector of row names that will be used as the reference when
#'     calculating logFCs and statistics. If more than one value is provided, the
#'     sum of the corresponding counts is used to generate offsets. If NULL,
#'     offsets will be defined as the effective library sizes (using TMM
#'     normalization factors).
#' @param selAssay Assay to select from \code{se} for the analysis.
#' @param pseudocount Pseudocount to add when calculating log-fold changes. 
#' @param method Either 'edgeR' or 'limma'. If set to 'limma', voom is used to
#'     transform the counts and estimate observation weights before applying
#'     limma. In this case, the results also contain the standard errors of the
#'     logFCs.
#' @param csawnorm Logical scalar. If \code{TRUE}, uses 
#'     \code{csaw::normOffsets()} to estimate offsets for use with \code{edgeR}, 
#'     rather than the default TMM-based size factors.
#' 
#' @author Charlotte Soneson, Michael Stadler
#' 
#' @export
#' 
#' @importFrom edgeR DGEList scaleOffset estimateDisp glmQLFit glmQLFTest
#'     topTags predFC topTags calcNormFactors
#' @importFrom SummarizedExperiment colData assay assayNames
#' @importFrom limma voom eBayes topTable lmFit
#' @importFrom csaw normOffsets
#' 
calculateRelativeFC <- function(se, design, coef = NULL, contrast = NULL, 
                                WTrows, selAssay = "counts", 
                                pseudocount = 1, method = "edgeR", 
                                csawnorm = FALSE) {
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
             " of columns in 'se' (", ncol(se), ").")
    }
    
    if (!is.numeric(pseudocount) || length(pseudocount) != 1 || pseudocount < 0) {
        stop("'pseudocount' must be a non-negative scalar value.")
    }
    
    if (csawnorm && !is.null(WTrows)) {
        stop("'csawnorm = TRUE' can only be used with WTrows = NULL.")
    }
    
    if (csawnorm && method == "limma") {
        stop("'csawnorm = TRUE' can only be used with method = 'edgeR'.")
    }
    
    if (is.null(coef) && is.null(contrast)) {
        stop("Both 'coef' and 'contrast' can not be NULL.")
    }
    
    ## Create DGEList from SummarizedExperiment
    dge <- edgeR::DGEList(counts = as.matrix(SummarizedExperiment::assay(se, selAssay)),
                          samples = SummarizedExperiment::colData(se))
    if (csawnorm) {
        ## csaw normalization - also calculate normalization factors since 
        ## aveLogCPM does not use provided offsets
        ## In this case, we know that WTrows is NULL, so all features 
        ## will be used for the normalization
        dge <- edgeR::calcNormFactors(dge)
        dge <- csaw::normOffsets(dge)
    } else if (is.null(WTrows)) {
        ## TMM normalization, with all features
        dge <- edgeR::calcNormFactors(dge)
    } else {
        ## Use size factors (offsets) derived from the geometric mean 
        ## of the WT rows
        tmp0 <- dge$counts[WTrows, , drop = FALSE]
        tmp0 <- tmp0[apply(tmp0, 1, min) > 0, , drop = FALSE]
        offsets <- apply(tmp0, 2, function(s) exp(mean(log(s))))
        offsets <- matrix(offsets, nrow = nrow(dge), ncol = length(offsets), byrow = TRUE)
        dge <- edgeR::scaleOffset(dge, log(offsets))
    }
    
    ## Fit model and perform test
    if (method == "edgeR") {
        dge <- edgeR::estimateDisp(dge, design = design)
        fit <- edgeR::glmQLFit(dge, design = design)
        qlf <- edgeR::glmQLFTest(fit, coef = coef, contrast = contrast)
        tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
        ## Calculate shrunken fold changes. Only when testing 
        ## a single coefficient or contrast
        predfc <- edgeR::predFC(dge, design = design, prior.count = pseudocount)
        if (length(coef) == 1 && is.null(contrast)) {
            tt$logFC_shrunk <- predfc[, coef]
        } else if (!is.null(contrast)) {
            tt$logFC_shrunk <- predfc %*% cbind(contrast)
        }
    } else if (method == "limma") {
        if (!is.null(dge$offset)) {
            vm <- limma::voom(dge, design = design, lib.size = exp(dge$offset))
        } else {
            vm <- limma::voom(dge, design = design, lib.size = effectiveLibSizes(dge))
        }
        fit <- limma::lmFit(vm, design = design)
        if (!is.null(contrast)) {
            fit <- contrasts.fit(fit, contrasts = contrast)
            coef <- 1
        }
        fit <- limma::eBayes(fit)
        tt <- limma::topTable(fit, coef = coef,  
                              confint = TRUE, number = Inf, sort.by = "none")
        if (length(coef) == 1 && is.null(contrast)) {
            tt$se.logFC <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
        }
    } else {
        stop("'method' must be either 'edgeR' or 'limma'")
    }
    tt
}
