context("plotFiltering")

se <- readRDS(system.file("extdata/GSE102901_cis_se.rds", package = "mutscan"))
se <- se[1:1000, 1:3]

test_that("plotFiltering fails with incorrect arguments", {
    expect_error(plotFiltering(se = 1))
    expect_error(plotFiltering(se = "x"))
    expect_error(plotFiltering(se = matrix(1:6), 2, 3))
    
    expect_error(plotFiltering(se = se, valueType = 1))
    expect_error(plotFiltering(se = se, valueType = c("reads", "fractions")))
    expect_error(plotFiltering(se = se, valueType = "wrong"))
    
    expect_error(plotFiltering(se = se, onlyActiveFilters = 1))
    expect_error(plotFiltering(se = se, onlyActiveFilters = "TRUE"))
    expect_error(plotFiltering(se = se, onlyActiveFilters = c(TRUE, FALSE)))
    
    expect_error(plotFiltering(se = se, displayNumbers = 1))
    expect_error(plotFiltering(se = se, displayNumbers = "TRUE"))
    expect_error(plotFiltering(se = se, displayNumbers = c(TRUE, FALSE)))
    
    expect_error(plotFiltering(se = se, numberSize = "1"))
    expect_error(plotFiltering(se = se, numberSize = TRUE))
    expect_error(plotFiltering(se = se, numberSize = c(1, 2)))
    expect_error(plotFiltering(se = se, numberSize = 0))
    
    expect_error(plotFiltering(se = se, plotType = 1))
    expect_error(plotFiltering(se = se, plotType = TRUE))
    expect_error(plotFiltering(se = se, plotType = c("remaining", "filtered")))
    expect_error(plotFiltering(se = se, plotType = "wrong"))
    
    expect_error(plotFiltering(se = se, facetBy = 1))
    expect_error(plotFiltering(se = se, facetBy = TRUE))
    expect_error(plotFiltering(se = se, facetBy = c("sample", "step")))
    expect_error(plotFiltering(se = se, facetBy = "wrong"))
})

test_that("plotFiltering works with correct arguments", {
    ## Defaults
    expect_is(plotFiltering(se = se, valueType = "reads", 
                            onlyActiveFilters = FALSE, 
                            displayNumbers = TRUE, numberSize = 4,
                            plotType = "remaining", facetBy = "sample"), "ggplot")
    
    ## Fractions
    expect_is(plotFiltering(se = se, valueType = "fractions", 
                            onlyActiveFilters = FALSE, 
                            displayNumbers = TRUE, numberSize = 4,
                            plotType = "remaining", facetBy = "sample"), "ggplot")
    
    ## Only active filters
    expect_is(plotFiltering(se = se, valueType = "reads", 
                            onlyActiveFilters = TRUE, 
                            displayNumbers = TRUE, numberSize = 4,
                            plotType = "remaining", facetBy = "sample"), "ggplot")
    
    ## Don't display numbers
    expect_is(plotFiltering(se = se, valueType = "reads", 
                            onlyActiveFilters = FALSE, 
                            displayNumbers = FALSE, numberSize = 4,
                            plotType = "remaining", facetBy = "sample"), "ggplot")
    
    ## Change size of displayed numbers
    expect_is(plotFiltering(se = se, valueType = "reads", 
                            onlyActiveFilters = FALSE, 
                            displayNumbers = FALSE, numberSize = 2,
                            plotType = "remaining", facetBy = "sample"), "ggplot")
    
    ## Reads + Filtered + Sample
    expect_is(plotFiltering(se = se, valueType = "reads", 
                            onlyActiveFilters = FALSE, 
                            displayNumbers = TRUE, numberSize = 4,
                            plotType = "filtered", facetBy = "sample"), "ggplot")

    ## Reads + Filtered + Step
    expect_is(plotFiltering(se = se, valueType = "reads", 
                            onlyActiveFilters = FALSE, 
                            displayNumbers = TRUE, numberSize = 4,
                            plotType = "filtered", facetBy = "step"), "ggplot")

    ## Reads + Remaining + Step
    expect_is(plotFiltering(se = se, valueType = "reads", 
                            onlyActiveFilters = FALSE, 
                            displayNumbers = TRUE, numberSize = 4,
                            plotType = "remaining", facetBy = "step"), "ggplot")

    ## Fractions + Filtered + Sample
    expect_is(plotFiltering(se = se, valueType = "fractions", 
                            onlyActiveFilters = FALSE, 
                            displayNumbers = TRUE, numberSize = 4,
                            plotType = "filtered", facetBy = "sample"), "ggplot")
    
    ## Fractions + Filtered + Step
    expect_is(plotFiltering(se = se, valueType = "fractions", 
                            onlyActiveFilters = FALSE, 
                            displayNumbers = TRUE, numberSize = 4,
                            plotType = "filtered", facetBy = "step"), "ggplot")
    
    ## Fractions + Remaining + Step
    expect_is(plotFiltering(se = se, valueType = "fractions", 
                            onlyActiveFilters = FALSE, 
                            displayNumbers = TRUE, numberSize = 4,
                            plotType = "remaining", facetBy = "step"), "ggplot")
    
})
