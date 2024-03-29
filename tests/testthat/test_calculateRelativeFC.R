se <- readRDS(system.file("extdata/GSE102901_cis_se.rds", package = "mutscan"))
se <- se[1:200, ]
design <- stats::model.matrix(~ Replicate + Condition, 
                              data = SummarizedExperiment::colData(se))

test_that("calculateRelativeFC fails with incorrect arguments", {
    expect_error(calculateRelativeFC(se = NULL, design = design,
                                     coef = "Conditioncis_output"),
                 "'se' must not be NULL")
    expect_error(calculateRelativeFC(se = as.matrix(assay(se)), design = design,
                                     coef = "Conditioncis_output"),
                 "'se' must be of class 'SummarizedExperiment'")
    expect_error(calculateRelativeFC(se = se, design = NULL,
                                     coef = "Conditioncis_output"),
                 "argument is of length zero")
    expect_error(calculateRelativeFC(se = se, design = design[1:3, ],
                                     coef = "Conditioncis_output"),
                 "The number of rows in")
    expect_error(calculateRelativeFC(se = se, design = design,
                                     coef = "missing"),
                 "One or more named coef arguments")
    expect_error(calculateRelativeFC(se = se, design = design,
                                     coef = NULL, contrast = NULL),
                 "can not both be NULL")
    expect_error(calculateRelativeFC(se = se, design = design,
                                     coef = NULL,
                                     contrast = matrix(sample(c(0, 1),
                                                              ncol(se) * 2,
                                                              replace = TRUE),
                                                       ncol = 2)),
                 "'contrast' must be a vector")
    expect_error(calculateRelativeFC(se = se, design = design,
                                     coef = "Conditioncis_output",
                                     WTrows = "missing"),
                 "All values in 'WTrows' must be")
    expect_error(calculateRelativeFC(se = se, design = design,
                                     coef = "Conditioncis_output",
                                     selAssay = "missing"),
                 "The provided 'selAssay' not present")
    expect_error(calculateRelativeFC(se = se, design = design,
                                     coef = "Conditioncis_output",
                                     selAssay = 1),
                 "'selAssay' must be of class 'character'")
    expect_error(calculateRelativeFC(se = se, design = design,
                                     coef = "Conditioncis_output",
                                     pseudocount = -1),
                 "must be within")
    expect_error(calculateRelativeFC(se = se, design = design,
                                     coef = "Conditioncis_output",
                                     method = "missing"),
                 "'method' must be one of: edgeR, limma")
    expect_error(calculateRelativeFC(se = se, design = design,
                                     coef = "Conditioncis_output",
                                     WTrows = NULL, normMethod = "sum"),
                 "normMethod = 'sum' can only be used when WTrows is not NULL")
    expect_error(calculateRelativeFC(se = se, design = design,
                                     coef = "Conditioncis_output",
                                     WTrows = NULL, normMethod = "geomean"),
                 "normMethod = 'geomean' can only be used when WTrows is not NULL")
    expect_error(calculateRelativeFC(se = se, design = design,
                                     coef = "Conditioncis_output",
                                     WTrows = rownames(se)[1], normMethod = "TMM"),
                 "normMethod = 'TMM' can only be used when WTrows is NULL")
    expect_error(calculateRelativeFC(se = se, design = design,
                                     coef = "Conditioncis_output",
                                     WTrows = rownames(se)[1], normMethod = "csaw"),
                 "normMethod = 'csaw' can only be used when WTrows is NULL")
    expect_error(calculateRelativeFC(se = se, design = design,
                                     coef = "Conditioncis_output",
                                     WTrows = NULL, normMethod = "csaw",
                                     method = "limma"),
                 "normMethod = 'csaw' can only be used with method = 'edgeR'")
})

test_that("calculateRelativeFC works", {
    ## With a single WT row, normMethod = "sum" or "geomean" should be identical

    res1a <- calculateRelativeFC(se, design, coef = "Conditioncis_output",
                                 contrast = NULL, WTrows = "f.0.WT",
                                 selAssay = "counts", pseudocount = 1,
                                 method = "edgeR", normMethod = "sum")
    res1b <- calculateRelativeFC(se, design, coef = "Conditioncis_output",
                                 contrast = NULL, WTrows = "f.0.WT",
                                 selAssay = "counts", pseudocount = 1,
                                 method = "edgeR", normMethod = "geomean")
    ## If there's only one assay, it doesn't need to be named
    senoname <- SummarizedExperiment(
        assays = list(assay(se, "counts")),
        rowData = rowData(se),
        colData = colData(se),
        metadata = metadata(se)
    )
    expect_warning(
        res1c <- calculateRelativeFC(senoname, design, coef = "Conditioncis_output",
                             contrast = NULL, WTrows = "f.0.WT", selAssay = "counts",
                             pseudocount = 1, method = "edgeR",
                             normMethod = "sum"),
        "No assayNames provided in 'se'")
    expect_equal(res1a$logFC, res1b$logFC)
    expect_equal(res1a$logFC, res1c$logFC)
    expect_equal(res1a$logFC_shrunk, res1b$logFC_shrunk)
    expect_equal(res1a$logFC_shrunk, res1c$logFC_shrunk)
    expect_equal(res1a$logCPM, res1b$logCPM)
    expect_equal(res1a$logCPM, res1c$logCPM)
    expect_equal(res1a$F, res1b$F)
    expect_equal(res1a$F, res1c$F)
    expect_equal(res1a$PValue, res1b$PValue)
    expect_equal(res1a$PValue, res1c$PValue)
    expect_equal(res1a$FDR, res1b$FDR)
    expect_equal(res1a$FDR, res1c$FDR)

    ## Specifying the coef or contrast should be equivalent
    res2b <- calculateRelativeFC(se, design, coef = NULL,
                                 contrast = as.numeric(colnames(design) == "Conditioncis_output"),
                                 WTrows = "f.0.WT", selAssay = "counts",
                                 pseudocount = 1, method = "edgeR",
                                 normMethod = "sum")
    expect_equal(res1a$logFC, res2b$logFC)
    expect_equal(res1a$logFC_shrunk, res2b$logFC_shrunk)
    expect_equal(res1a$logCPM, res2b$logCPM)
    expect_equal(res1a$F, res2b$F)
    expect_equal(res1a$PValue, res2b$PValue)
    expect_equal(res1a$FDR, res2b$FDR)

    ## With a single WT row, that one should not show any change
    expect_equal(res1a["f.0.WT", "logFC"], 0, tolerance = 1e-4)
    expect_equal(res1a["f.0.WT", "logFC_shrunk"], 0, tolerance = 1e-4)
    expect_equal(res1a["f.0.WT", "PValue"], 1, tolerance = 1e-4)
    expect_equal(res1a["f.0.WT", "FDR"], 1, tolerance = 1e-4)
    expect_named(res1a, c("logFC", "logCPM", "F", "PValue", "FDR",
                          "logFC_shrunk", "df.total", "df.prior", "df.test"))

    ## Should also work with limma
    res3a <- calculateRelativeFC(se, design, coef = "Conditioncis_output",
                                 contrast = NULL, WTrows = "f.0.WT", selAssay = "counts",
                                 pseudocount = 1, method = "limma",
                                 normMethod = "sum")
    expect_equal(res3a["f.0.WT", "logFC"], 0, tolerance = 1e-4)
    expect_equal(res3a["f.0.WT", "P.Value"], 1, tolerance = 1e-4)
    expect_equal(res3a["f.0.WT", "adj.P.Val"], 1, tolerance = 1e-4)
    expect_named(res3a, c("logFC", "CI.L", "CI.R", "AveExpr", "t",
                          "P.Value", "adj.P.Val", "B", "se.logFC", 
                          "df.total", "df.prior"))
    
    ## Also with contrast
    res3b <- calculateRelativeFC(se, design, coef = NULL,
                                 contrast = as.numeric(colnames(design) == "Conditioncis_output"),
                                 WTrows = "f.0.WT", selAssay = "counts",
                                 pseudocount = 1, method = "limma",
                                 normMethod = "sum")
    expect_equal(res3b["f.0.WT", "logFC"], 0, tolerance = 1e-4)
    expect_equal(res3b["f.0.WT", "P.Value"], 1, tolerance = 1e-4)
    expect_equal(res3b["f.0.WT", "adj.P.Val"], 1, tolerance = 1e-4)
    expect_named(res3b, c("logFC", "CI.L", "CI.R", "AveExpr", "t",
                          "P.Value", "adj.P.Val", "B", "se.logFC",
                          "df.total", "df.prior"))
    
    ## Coef and contrast should be equivalent also for limma
    expect_equal(res3a$logFC, res3b$logFC)
    expect_equal(res3a$CI.L, res3b$CI.L)
    expect_equal(res3a$CI.R, res3b$CI.R)
    expect_equal(res3a$AveExpr, res3b$AveExpr)
    expect_equal(res3a$t, res3b$t)
    expect_equal(res3a$P.Value, res3b$P.Value)
    expect_equal(res3a$adj.P.Val, res3b$adj.P.Val)
    expect_equal(res3a$B, res3b$B)
    expect_true(all(res3a$se.logFC >= 0))
    expect_equal(res3a$se.logFC, res3b$se.logFC)
    
    ## Correlation between edgeR and limma logFCs should be high
    expect_gt(cor(res1a$logFC, res3a$logFC), 0.98)

    ## Run without WTrows
    res4a <- calculateRelativeFC(se, design, coef = "Conditioncis_output",
                                 contrast = NULL, WTrows = NULL, selAssay = "counts",
                                 pseudocount = 1, method = "edgeR",
                                 normMethod = "TMM")
    res4b <- calculateRelativeFC(se, design, coef = "Conditioncis_output",
                                 contrast = NULL, WTrows = NULL, selAssay = "counts",
                                 pseudocount = 1, method = "limma",
                                 normMethod = "TMM")
    res4c <- calculateRelativeFC(se, design, coef = "Conditioncis_output",
                                 contrast = NULL, WTrows = NULL, selAssay = "counts",
                                 pseudocount = 1, method = "edgeR",
                                 normMethod = "csaw")
    expect_gt(cor(res4a$logFC, res4b$logFC), 0.98)
})

