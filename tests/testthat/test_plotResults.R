context("plotResults")

se <- readRDS(system.file("extdata/GSE102901_cis_se.rds", package = "mutscan"))
se <- se[1:100, ]
design <- model.matrix(~ Replicate + Condition, data = colData(se))
res0 <- calculateRelativeFC(se, design, coef = "Conditioncis_output", 
                            contrast = NULL, WTrows = "f.0.WT", selAssay = "counts", 
                            pseudocount = 1, method = "edgeR", 
                            normMethod = "sum")

test_that("plotResults functions fail with incorrect input", {
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
    
    expect_is(plotMeanDiff(res1), "ggplot")
    expect_is(plotMeanDiff(res2), "ggplot")
    expect_is(plotVolcano(res1), "ggplot")
    expect_is(plotVolcano(res2), "ggplot")

    expect_is(plotMeanDiff(res1, pointSize = "large"), "ggplot")
    expect_is(plotMeanDiff(res2, pointSize = "large"), "ggplot")
    expect_is(plotVolcano(res1, pointSize = "large"), "ggplot")
    expect_is(plotVolcano(res2, pointSize = "large"), "ggplot")
    
    ## These tests won't run if the X11 display connection can not be opened
    # skip_if_not_installed("plotly")
    # expect_is(plotMeanDiff(res1, interactivePlot = TRUE), "plotly")
    # expect_is(plotMeanDiff(res2, interactivePlot = TRUE), "plotly")
    # expect_is(plotVolcano(res1, interactivePlot = TRUE), "plotly")
    # expect_is(plotVolcano(res2, interactivePlot = TRUE), "plotly")
})
