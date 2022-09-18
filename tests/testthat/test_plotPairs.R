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
    
    expect_error(plotPairs(se = se, selAssay = "counts", corrColorRange = 1))
    expect_error(plotPairs(se = se, selAssay = "counts", corrColorRange = c(-1, 0.5)))
    expect_error(plotPairs(se = se, selAssay = "counts", corrColorRange = c(1, 0.5)))
    expect_error(plotPairs(se = se, selAssay = "counts", corrColorRange = c("0.5", "1")))
    expect_error(plotPairs(se = se, selAssay = "counts", corrColorRange = c(0, 2)))
    
    expect_error(plotPairs(se = se, selAssay = "counts", addIdentityLine = 1))
    expect_error(plotPairs(se = se, selAssay = "counts", addIdentityLine = "TRUE"))
    expect_error(plotPairs(se = se, selAssay = "counts", addIdentityLine = c(TRUE, FALSE)))
})

test_that("plotPairs works with correct arguments", {
    ## Defaults
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson",
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE, corrColorRange = NULL,
                              addIdentityLine = FALSE), "ggmatrix")

    ## Not log-transformed
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = FALSE,
                              pseudocount = 1, corMethod = "pearson",
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE, corrColorRange = NULL, 
                              addIdentityLine = FALSE), "ggmatrix")

    ## No pseudocount
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 0, corMethod = "pearson",
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE, corrColorRange = NULL, 
                              addIdentityLine = FALSE), "ggmatrix")

    ## Spearman correlation
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "spearman",
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE, corrColorRange = NULL, 
                              addIdentityLine = FALSE), "ggmatrix")

    ## Fewer breaks
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson",
                              histBreaks = 5, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE, corrColorRange = NULL, 
                              addIdentityLine = FALSE), "ggmatrix")

    ## smoothScatter
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson",
                              histBreaks = 40, pointsType = "smoothscatter", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE, corrColorRange = NULL, 
                              addIdentityLine = FALSE), "ggmatrix")

    ## Change font size to correlation relation
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson",
                              histBreaks = 40, pointsType = "points", corSizeMult = 2,
                              corSizeAdd = 0, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE, corrColorRange = NULL, 
                              addIdentityLine = FALSE), "ggmatrix")

    ## Point size, alpha
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson",
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 1, pointAlpha = 3,
                              colorByCorrelation = TRUE, corrColorRange = NULL, 
                              addIdentityLine = FALSE), "ggmatrix")

    ## Don't color by correlation
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson",
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = FALSE, corrColorRange = NULL, 
                              addIdentityLine = FALSE), "ggmatrix")
    
    ## Change color range
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson",
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE, corrColorRange = c(0.7, 1), 
                              addIdentityLine = FALSE), "ggmatrix")
    
    ## Add identity line
    expect_s3_class(plotPairs(se = se, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson",
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE, corrColorRange = NULL, 
                              addIdentityLine = TRUE), "ggmatrix")

    ## Make one correlation negative
    seneg <- se
    SummarizedExperiment::assay(seneg, 1)[, 2] <-
        1/(SummarizedExperiment::assay(seneg, 1)[, 2] + 1)
    expect_s3_class(plotPairs(se = seneg, selAssay = "counts", doLog = TRUE,
                              pseudocount = 1, corMethod = "pearson",
                              histBreaks = 40, pointsType = "points", corSizeMult = 5,
                              corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                              colorByCorrelation = TRUE, corrColorRange = NULL,
                              addIdentityLine = FALSE), "ggmatrix")
})