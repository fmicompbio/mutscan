context("generateQCReport")

se <- readRDS(system.file("extdata/GSE102901_cis_se.rds", package = "mutscan"))
se <- se[1:100, ]
outFile <- tempfile(fileext = ".html")

test_that("generateQCReport fails with incorrect input", {
    expect_error(generateQCReport(se = NULL, outFile = outFile))
    expect_error(generateQCReport(se = "file", outFile = outFile))
    
    expect_error(generateQCReport(se = se, outFile = 1))
    expect_error(generateQCReport(se = se, outFile = "out.pdf"))
    expect_error(generateQCReport(se = se, outFile = c("out1.html", "out2.html")))

    expect_error(generateQCReport(se = se, outFile = outFile, forceOverwrite = 1))
    expect_error(generateQCReport(se = se, outFile = outFile, forceOverwrite = "wrong"))
    expect_error(generateQCReport(se = se, outFile = outFile, forceOverwrite = c(TRUE, FALSE)))
})

test_that("generateQCReport works as expected", {
    expect_equal(generateQCReport(se = se, outFile = outFile, forceOverwrite = FALSE), 
                 outFile)
    expect_true(file.exists(outFile))
    expect_error(generateQCReport(se = se, outFile = outFile, forceOverwrite = FALSE))
    expect_equal(generateQCReport(se = se, outFile = outFile, forceOverwrite = TRUE), 
                 outFile)
})