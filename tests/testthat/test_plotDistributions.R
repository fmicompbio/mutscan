se <- readRDS(system.file("extdata/GSE102901_cis_se.rds", package = "mutscan"))
se <- se[1:1000, 1:3]

test_that("plotDistributions fails with incorrect arguments", {
    expect_error(plotDistributions(se = 1))
    expect_error(plotDistributions(se = "x"))
    expect_error(plotDistributions(se = matrix(1:6), 2, 3))
    
    expect_error(plotDistributions(se = se, selAssay = 1))
    expect_error(plotDistributions(se = se, selAssay = c("counts", "counts")))
    expect_error(plotDistributions(se = se, selAssay = TRUE))
    expect_error(plotDistributions(se = se, selAssay = "wrong"))
    
    expect_error(plotDistributions(se = se, groupBy = 1))
    expect_error(plotDistributions(se = se, groupBy = "wrong"))
    expect_error(plotDistributions(se = se, groupBy = c("Name", "Condition")))

    expect_error(plotDistributions(se = se, plotType = 1))
    expect_error(plotDistributions(se = se, plotType = TRUE))
    expect_error(plotDistributions(se = se, plotType = c("density", "histogram")))
    expect_error(plotDistributions(se = se, plotType = "wrong"))
    
    expect_error(plotDistributions(se = se, facet = 1))
    expect_error(plotDistributions(se = se, facet = c(TRUE, FALSE)))
    expect_error(plotDistributions(se = se, facet = "wrong"))
    
    expect_error(plotDistributions(se = se, pseudocount = TRUE))
    expect_error(plotDistributions(se = se, pseudocount = "1"))
    expect_error(plotDistributions(se = se, pseudocount = c(0, 1)))
    expect_error(plotDistributions(se = se, pseudocount = -1))
})

test_that("plotDistributions works as expected", {
    ## Defaults
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = NULL, plotType = "density",
                                      facet = FALSE, pseudocount = 0), "ggplot")

    ## Change plot type
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = NULL, plotType = "histogram",
                                      facet = FALSE, pseudocount = 0), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = NULL, plotType = "knee",
                                      facet = FALSE, pseudocount = 0), "ggplot")
    
    ## groupBy Name
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Name", plotType = "density",
                                      facet = FALSE, pseudocount = 0), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Name", plotType = "histogram",
                                      facet = FALSE, pseudocount = 0), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Name", plotType = "knee",
                                      facet = FALSE, pseudocount = 0), "ggplot")
    
    ## groupBy Condition
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Condition", plotType = "density",
                                      facet = FALSE, pseudocount = 0), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Condition", plotType = "histogram",
                                      facet = FALSE, pseudocount = 0), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Condition", plotType = "knee",
                                      facet = FALSE, pseudocount = 0), "ggplot")
    
    ## Facet
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = NULL, plotType = "density",
                                      facet = TRUE, pseudocount = 0), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = NULL, plotType = "histogram",
                                      facet = TRUE, pseudocount = 0), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = NULL, plotType = "knee",
                                      facet = TRUE, pseudocount = 0), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Name", plotType = "density",
                                      facet = TRUE, pseudocount = 0), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Name", plotType = "histogram",
                                      facet = TRUE, pseudocount = 0), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Name", plotType = "knee",
                                      facet = TRUE, pseudocount = 0), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Condition", plotType = "density",
                                      facet = TRUE, pseudocount = 0), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Condition", plotType = "histogram",
                                      facet = TRUE, pseudocount = 0), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Condition", plotType = "knee",
                                      facet = TRUE, pseudocount = 0), "ggplot")
    
    ## Increase pseudocount
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Condition", plotType = "density",
                                      facet = FALSE, pseudocount = 2), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Condition", plotType = "histogram",
                                      facet = FALSE, pseudocount = 3), "ggplot")
    expect_s3_class(plotDistributions(se = se, selAssay = "counts",
                                      groupBy = "Condition", plotType = "knee",
                                      facet = FALSE, pseudocount = 4), "ggplot")
    
})