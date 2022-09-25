se <- readRDS(system.file("extdata/GSE102901_cis_se.rds", package = "mutscan"))
se <- se[1:200, ]
design <- stats::model.matrix(~ Replicate + Condition, 
                              data = SummarizedExperiment::colData(se))
res0 <- calculateRelativeFC(se, design, coef = "Conditioncis_output",
                            contrast = NULL, WTrows = "f.0.WT", selAssay = "counts",
                            pseudocount = 1, method = "edgeR",
                            normMethod = "sum")

test_that("plotResults functions fail with incorrect input", {

    expect_error(plotMeanDiff(as.matrix(res0)))
    expect_error(plotVolcano(as.matrix(res0)))

    expect_error(plotMeanDiff(res0, meanCol = 1))
    expect_error(plotMeanDiff(res0, meanCol = "missing"))
    expect_error(plotMeanDiff(res0, meanCol = c("AveExpr", "logCPM")))
    expect_error(plotMeanDiff(res0, logFCCol = 1))
    expect_error(plotMeanDiff(res0, logFCCol = "missing"))
    expect_error(plotMeanDiff(res0, logFCCol = c("AveExpr", "logFC")))
    expect_error(plotMeanDiff(res0, pvalCol = 1))
    expect_error(plotMeanDiff(res0, pvalCol = "missing"))
    expect_error(plotMeanDiff(res0, pvalCol = c("PValue", "P.Value")))
    expect_error(plotMeanDiff(res0, padjCol = 1))
    expect_error(plotMeanDiff(res0, padjCol = "missing"))
    expect_error(plotMeanDiff(res0, padjCol = c("FDR", "adj.P.Val")))

    expect_error(plotVolcano(res0, logFCCol = 1))
    expect_error(plotVolcano(res0, logFCCol = "missing"))
    expect_error(plotVolcano(res0, logFCCol = c("AveExpr", "logFC")))
    expect_error(plotVolcano(res0, pvalCol = 1))
    expect_error(plotVolcano(res0, pvalCol = "missing"))
    expect_error(plotVolcano(res0, pvalCol = c("PValue", "P.Value")))
    expect_error(plotVolcano(res0, padjCol = 1))
    expect_error(plotVolcano(res0, padjCol = "missing"))
    expect_error(plotVolcano(res0, padjCol = c("FDR", "adj.P.Val")))

    expect_error(plotMeanDiff(res0, padjThreshold = "0.05"))
    expect_error(plotMeanDiff(res0, padjThreshold = c(0.01, 0.05)))
    expect_error(plotVolcano(res0, padjThreshold = "0.05"))
    expect_error(plotVolcano(res0, padjThreshold = c(0.01, 0.05)))

    expect_error(plotMeanDiff(res0, pointSize = 1))
    expect_error(plotMeanDiff(res0, pointSize = c("small", "large")))
    expect_error(plotMeanDiff(res0, pointSize = "wrong"))
    expect_error(plotVolcano(res0, pointSize = 1))
    expect_error(plotVolcano(res0, pointSize = c("small", "large")))
    expect_error(plotVolcano(res0, pointSize = "wrong"))

    expect_error(plotMeanDiff(res0, interactivePlot = 1))
    expect_error(plotMeanDiff(res0, interactivePlot = c(TRUE, FALSE)))
    expect_error(plotMeanDiff(res0, interactivePlot = "TRUE"))
    expect_error(plotVolcano(res0, interactivePlot = 1))
    expect_error(plotVolcano(res0, interactivePlot = c(TRUE, FALSE)))
    expect_error(plotVolcano(res0, interactivePlot = "TRUE"))

    expect_error(plotMeanDiff(res0, nTopToLabel = "1"))
    expect_error(plotMeanDiff(res0, nTopToLabel = c(1, 2)))
    expect_error(plotMeanDiff(res0, nTopToLabel = -1))
    expect_error(plotVolcano(res0, nTopToLabel = "1"))
    expect_error(plotVolcano(res0, nTopToLabel = c(1, 2)))
    expect_error(plotVolcano(res0, nTopToLabel = -1))

    res1 <- res0
    colnames(res1)[colnames(res1) == "logCPM"] <- "wrong"
    expect_error(plotMeanDiff(res1),
                 "No suitable column found for x-axis")

    res1 <- res0
    colnames(res1)[colnames(res1) == "logFC"] <- "wrong"
    expect_error(plotMeanDiff(res1),
                 "No suitable column found for y-axis")
    expect_error(plotVolcano(res1),
                 "No suitable column found for x-axis")

    res1 <- res0
    colnames(res1)[colnames(res1) == "FDR"] <- "wrong"
    expect_error(plotMeanDiff(res1),
                 "No suitable column found for highlighting")
    expect_error(plotVolcano(res1),
                 "No suitable column found for highlighting")

    res1 <- res0
    colnames(res1)[colnames(res1) == "PValue"] <- "wrong"
    expect_error(plotVolcano(res1),
                 "No suitable column found for y-axis")

    res1 <- res0
    res1$AveExpr <- res1$logCPM
    expect_warning(plotMeanDiff(res1),
                   "Multiple suitable columns found for x-axis")

    res1 <- res0
    res1$P.Value <- res1$PValue
    expect_warning(plotVolcano(res1),
                   "Multiple suitable columns found for y-axis")

    res1 <- res0
    res1$adj.P.Val <- res1$FDR
    expect_warning(plotMeanDiff(res1),
                   "Multiple suitable columns found for highlighting")
    expect_warning(plotVolcano(res1),
                   "Multiple suitable columns found for highlighting")
})

test_that("plotResults functions work as expected", {
    res1 <- calculateRelativeFC(se, design, coef = "Conditioncis_output",
                                contrast = NULL, WTrows = "f.0.WT", selAssay = "counts",
                                pseudocount = 1, method = "edgeR",
                                normMethod = "sum")
    res2 <- calculateRelativeFC(se, design, coef = "Conditioncis_output",
                                contrast = NULL, WTrows = "f.0.WT", selAssay = "counts",
                                pseudocount = 1, method = "limma",
                                normMethod = "sum")

    expect_s3_class(plotMeanDiff(res1), "ggplot")
    expect_s3_class(plotMeanDiff(res2), "ggplot")
    expect_s3_class(plotVolcano(res1), "ggplot")
    expect_s3_class(plotVolcano(res2), "ggplot")

    expect_s3_class(plotMeanDiff(res1, meanCol = "logCPM", logFCCol = "logFC",
                                 padjCol = "FDR", padjThreshold = 0.1), "ggplot")
    expect_s3_class(plotMeanDiff(res2, meanCol = "AveExpr", logFCCol = "logFC",
                                 padjCol = "adj.P.Val", padjThreshold = 0.1), "ggplot")
    expect_s3_class(plotVolcano(res1, logFCCol = "logFC", pvalCol = "PValue",
                                padjCol = "FDR", padjThreshold = 0.1), "ggplot")
    expect_s3_class(plotVolcano(res2, logFCCol = "logFC", pvalCol = "P.Value",
                                padjCol = "adj.P.Val", padjThreshold = 0.1), "ggplot")

    expect_s3_class(plotMeanDiff(res1, nTopToLabel = 5), "ggplot")
    expect_s3_class(plotMeanDiff(res2, nTopToLabel = 5), "ggplot")
    expect_s3_class(plotVolcano(res1, nTopToLabel = 5), "ggplot")
    expect_s3_class(plotVolcano(res2, nTopToLabel = 5), "ggplot")

    expect_s3_class(plotMeanDiff(res1, pointSize = "large"), "ggplot")
    expect_s3_class(plotMeanDiff(res2, pointSize = "large"), "ggplot")
    expect_s3_class(plotVolcano(res1, pointSize = "large"), "ggplot")
    expect_s3_class(plotVolcano(res2, pointSize = "large"), "ggplot")

    ## These tests won't run if the X11 display connection can not be opened
    # skip_if_not_installed("plotly")
    # expect_s3_class(plotMeanDiff(res1, interactivePlot = TRUE), "plotly")
    # expect_s3_class(plotMeanDiff(res2, interactivePlot = TRUE), "plotly")
    # expect_s3_class(plotVolcano(res1, interactivePlot = TRUE), "plotly")
    # expect_s3_class(plotVolcano(res2, interactivePlot = TRUE), "plotly")
})
