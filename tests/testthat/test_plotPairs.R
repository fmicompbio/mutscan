se <- readRDS(system.file("extdata/GSE102901_cis_se.rds", package = "mutscan"))
se <- se[1:1000, 1:3]

test_that("plotPairs fails with incorrect arguments", {
    expect_error(plotPairs(se = 1))
    expect_error(plotPairs(se = "x"))
    expect_error(plotPairs(se = matrix(1:6), 2, 3))
    
    expect_error(plotPairs(se = se, selAssay = 1))
    expect_error(plotPairs(se = se, selAssay = c("counts", "logcounts")))
    expect_error(plotPairs(se = se, selAssay = "logcounts"))
    
    expect_error(plotPairs(se = se, selAssay = "counts", doLog = 1))
    expect_error(plotPairs(se = se, selAssay = "counts", doLog = "TRUE"))
    expect_error(plotPairs(se = se, selAssay = "counts", doLog = c(TRUE, FALSE)))
    
    expect_error(plotPairs(se = se, selAssay = "counts", pseudocount = "1"))
    expect_error(plotPairs(se = se, selAssay = "counts", pseudocount = c(1, 2)))
    expect_error(plotPairs(se = se, selAssay = "counts", pseudocount = -1))
    
    expect_error(plotPairs(se = se, selAssay = "counts", corMethod = 1))
    expect_error(plotPairs(se = se, selAssay = "counts", 
                           corMethod = c("spearman", "pearson")))
    expect_error(plotPairs(se = se, selAssay = "counts", corMethod = "pearman"))
    
    expect_error(plotPairs(se = se, selAssay = "counts", histBreaks = "1"))
    expect_error(plotPairs(se = se, selAssay = "counts", histBreaks = c(1, 2)))
    expect_error(plotPairs(se = se, selAssay = "counts", histBreaks = 0))
    
    expect_error(plotPairs(se = se, selAssay = "counts", pointsType = 1))
    expect_error(plotPairs(se = se, selAssay = "counts", 
                           pointsType = c("smoothscatter", "points")))
    expect_error(plotPairs(se = se, selAssay = "counts", pointsType = "oints"))
    
    expect_error(plotPairs(se = se, selAssay = "counts", corSizeMult = "1"))
    expect_error(plotPairs(se = se, selAssay = "counts", corSizeMult = c(1, 2)))
    expect_error(plotPairs(se = se, selAssay = "counts", corSizeMult = 0))
    
    expect_error(plotPairs(se = se, selAssay = "counts", corSizeAdd = "1"))
    expect_error(plotPairs(se = se, selAssay = "counts", corSizeAdd = c(1, 2)))
    expect_error(plotPairs(se = se, selAssay = "counts", corSizeAdd = -1))
    
    expect_error(plotPairs(se = se, selAssay = "counts", pointSize = "1"))
    expect_error(plotPairs(se = se, selAssay = "counts", pointSize = c(1, 2)))
    expect_error(plotPairs(se = se, selAssay = "counts", pointSize = 0))
    
    expect_error(plotPairs(se = se, selAssay = "counts", pointAlpha = "1"))
    expect_error(plotPairs(se = se, selAssay = "counts", pointAlpha = c(1, 2)))
    expect_error(plotPairs(se = se, selAssay = "counts", pointAlpha = 0))
    
    expect_error(plotPairs(se = se, selAssay = "counts", colorByCorrelation = 1))
    expect_error(plotPairs(se = se, selAssay = "counts", colorByCorrelation = "TRUE"))
    expect_error(plotPairs(se = se, selAssay = "counts", colorByCorrelation = c(TRUE, FALSE)))
})

test_that("plotPairs works with correct arguments", {
    ## Defaults
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson", 
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE), "ggmatrix")
    
    ## Not log-transformed
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = FALSE,
                              pseudocount = 1, corMethod = "pearson", 
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE), "ggmatrix")
    
    ## No pseudocount
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 0, corMethod = "pearson", 
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE), "ggmatrix")
    
    ## Spearman correlation
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "spearman", 
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE), "ggmatrix")
    
    ## Fewer breaks
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson", 
                              histBreaks = 5, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE), "ggmatrix")
    
    ## smoothScatter
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson", 
                              histBreaks = 40, pointsType = "smoothscatter", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE), "ggmatrix")
    
    ## Change font size to correlation relation
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson", 
                              histBreaks = 40, pointsType = "points", corSizeMult = 2,
                              corSizeAdd = 0, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE), "ggmatrix")
    
    ## Point size, alpha
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson", 
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 1, pointAlpha = 3,
                              colorByCorrelation = TRUE), "ggmatrix")
    
    ## Don't color by correlation
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson", 
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = FALSE), "ggmatrix")
    
    ## Make one correlation negative
    seneg <- se
    SummarizedExperiment::assay(seneg, 1)[, 2] <- 
        1/(SummarizedExperiment::assay(seneg, 1)[, 2] + 1)
    expect_s3_class(plotPairs(se = seneg, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson", 
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE), "ggmatrix")
})