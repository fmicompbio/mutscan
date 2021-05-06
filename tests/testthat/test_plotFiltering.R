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
})

test_that("plotFiltering works with correct arguments", {
  ## Defaults
  expect_is(plotFiltering(se = se, valueType = "reads", 
                          onlyActiveFilters = FALSE, 
                          displayNumbers = TRUE, numberSize = 4), "ggplot")
  
  ## Fractions
  expect_is(plotFiltering(se = se, valueType = "fractions", 
                          onlyActiveFilters = FALSE, 
                          displayNumbers = TRUE, numberSize = 4), "ggplot")
  
  ## Only active filters
  expect_is(plotFiltering(se = se, valueType = "reads", 
                          onlyActiveFilters = TRUE, 
                          displayNumbers = TRUE, numberSize = 4), "ggplot")
  
  ## Don't display numbers
  expect_is(plotFiltering(se = se, valueType = "reads", 
                          onlyActiveFilters = FALSE, 
                          displayNumbers = FALSE, numberSize = 4), "ggplot")

  ## Change size of displayed numbers
  expect_is(plotFiltering(se = se, valueType = "reads", 
                          onlyActiveFilters = FALSE, 
                          displayNumbers = FALSE, numberSize = 2), "ggplot")
  
})
