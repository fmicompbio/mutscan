se <- readRDS(system.file("extdata/GSE102901_cis_se.rds", package = "mutscan"))
se <- se[1:1000, 1:3]
SummarizedExperiment::rowData(se)$categ <- sample(LETTERS[1:3], nrow(se), replace = TRUE)

test_that("plotTotals fails with incorrect arguments", {
    expect_error(plotTotals(se = 1))
    expect_error(plotTotals(se = "x"))
    expect_error(plotTotals(se = matrix(1:6), 2, 3))
    
    expect_error(plotTotals(se = se, selAssay = 1))
    expect_error(plotTotals(se = se, selAssay = c("counts", "counts")))
    expect_error(plotTotals(se = se, selAssay = TRUE))
    expect_error(plotTotals(se = se, selAssay = "wrong"))
    
    expect_error(plotTotals(se = se, groupBy = 1))
    expect_error(plotTotals(se = se, groupBy = "wrong"))
    expect_error(plotTotals(se = se, groupBy = c("Name", "Condition")))
})

test_that("plotTotals works as expected", {
    ## Defaults
    expect_s3_class(plotTotals(se = se, selAssay = "counts",
                               groupBy = NULL), "ggplot")

    ## Group by column
    expect_s3_class(plotTotals(se = se, selAssay = "counts",
                               groupBy = "categ"), "ggplot")
    
})