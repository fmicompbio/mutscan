#' Calculate logFCs relative to WT using edgeR
#'
#' Calculate logFCs and associated p-values for a given comparison, using 
#' either limma or the Negative Binomial quasi-likelihood framework of edgeR. 
#' The observed counts for the WT variants can be used as offsets in the model.
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
#' @param normMethod Character scalar indicating which normalization method 
#'     should be used to calculate size factors. Should be either \code{"TMM"} or 
#'     \code{"csaw"} when \code{WTrows} is \code{NULL}, and \code{"geomean"} or 
#'     \code{"sum"} when \code{WTrows} is provided.
#' 
#' @author Charlotte Soneson, Michael Stadler
#' 
#' @export
#' 
#' @importFrom edgeR DGEList scaleOffset estimateDisp glmQLFit glmQLFTest
#'     topTags predFC topTags calcNormFactors effectiveLibSizes
#' @importFrom SummarizedExperiment colData assay assayNames
#' @importFrom limma voom eBayes topTable lmFit contrasts.fit
#' @importFrom csaw normOffsets
#' 
#' @examples 
#' se <- readRDS(system.file("extdata", "GSE102901_cis_se.rds", 
#'                           package = "mutscan"))[1:200, ]
#' design <- model.matrix(~ Replicate + Condition, 
#'                        data = SummarizedExperiment::colData(se))
#' res <- calculateRelativeFC(se, design, coef = "Conditioncis_output")
#' 
calculateRelativeFC <- function(se, design, coef = NULL, contrast = NULL, 
                                WTrows = NULL, selAssay = "counts", 
                                pseudocount = 1, method = "edgeR", 
                                normMethod = ifelse(is.null(WTrows), 
                                                    "TMM", "sum")) {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = selAssay, type = "character")
    if (!(selAssay %in% SummarizedExperiment::assayNames(se))) {
        if (is.null(SummarizedExperiment::assayNames(se)) && 
            length(SummarizedExperiment::assays(se)) == 1) {
            warning("No assayNames provided in 'se', but only one ", 
                    "assay present - using that.")
            selAssay <- 1L
        } else {
            stop("The provided 'selAssay' not present in 'se'.")
        }
    }
    
    if (!is.null(WTrows)) {
        .assertVector(x = WTrows, type = "character", validValues = rownames(se))
    }
    
    if (nrow(design) != ncol(se)) {
        stop("The number of rows in 'design' (", nrow(design), 
             ") is not equal to the number", 
             " of columns in 'se' (", ncol(se), ").")
    }
    
    .assertScalar(x = pseudocount, type = "numeric", rngIncl = c(0, Inf))

    if (normMethod %in% c("csaw", "TMM") && !is.null(WTrows)) {
        stop("normMethod = '", normMethod, 
             "' can only be used when WTrows is NULL.")
    }

    if (normMethod %in% c("sum", "geomean") && is.null(WTrows)) {
        stop("normMethod = '", normMethod,
             "' can only be used when WTrows is not NULL.")
    }
    
    if (normMethod == "csaw" && method == "limma") {
        stop("normMethod = 'csaw' can only be used with method = 'edgeR'.")
    }
    
    .assertScalar(x = normMethod, type = "character", 
                  validValues = c("csaw", "TMM", "geomean", "sum"))
    
    if (is.null(coef) && is.null(contrast)) {
        stop("'coef' and 'contrast' can not both be NULL.")
    }
    
    if (!is.null(contrast) && !is.null(dim(contrast))) {
        stop("'contrast' must be a vector.")
    }
    
    ## Create DGEList from SummarizedExperiment
    dge <- edgeR::DGEList(counts = as.matrix(SummarizedExperiment::assay(se, selAssay)),
                          samples = SummarizedExperiment::colData(se))
    if (normMethod == "csaw") {
        ## csaw normalization - also calculate normalization factors since 
        ## aveLogCPM does not use provided offsets
        ## In this case, we know that WTrows is NULL, so all features 
        ## will be used for the normalization
        dge <- edgeR::calcNormFactors(dge)
        dge <- csaw::normOffsets(dge)
    } else if (normMethod == "TMM") {
        ## TMM normalization, with all features
        dge <- edgeR::calcNormFactors(dge)
    } else if (normMethod == "geomean") {
        ## Use size factors (offsets) derived from the geometric mean 
        ## of the WT rows
        tmp0 <- dge$counts[WTrows, , drop = FALSE]
        tmp0 <- tmp0[apply(tmp0, 1, min) > 0, , drop = FALSE]
        logoffsets <- apply(tmp0, 2, function(s) mean(log(s)))
        dge <- edgeR::scaleOffset(dge, logoffsets)
    } else if (normMethod == "sum") {
        ## Use size factors (offsets) derived from the sum of the 
        ## WT rows
        tmp0 <- dge$counts[WTrows, , drop = FALSE]
        dge <- edgeR::scaleOffset(dge, log(colSums(tmp0)))
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
            tt$logFC_shrunk <- c(predfc %*% cbind(contrast))
        }
    } else if (method == "limma") {
        if (!is.null(dge$offset)) {
            vm <- limma::voom(dge, design = design, lib.size = exp(dge$offset))
        } else {
            vm <- limma::voom(dge, design = design, 
                              lib.size = edgeR::effectiveLibSizes(dge))
        }
        fit <- limma::lmFit(vm, design = design)
        if (!is.null(contrast)) {
            fit <- limma::contrasts.fit(fit, contrasts = contrast)
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
