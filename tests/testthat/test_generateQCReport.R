se <- readRDS(system.file("extdata/GSE102901_cis_se.rds", package = "mutscan"))
se <- se[1:100, ]
outFile <- tempfile(fileext = ".html")
outFile2 <- file.path(tempdir(), "newfolder", "newreport.html")

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
    expect_warning(expect_equal(generateQCReport(se = se, outFile = outFile,
                                                 forceOverwrite = FALSE,
                                                 quiet = TRUE,
                                                 run_pandoc = FALSE),
                                outFile),
                   "Will ignore arguments")
    expect_true(file.exists(outFile))
    expect_error(generateQCReport(se = se, outFile = outFile, forceOverwrite = FALSE))
    expect_message(expect_equal(generateQCReport(se = se, outFile = outFile,
                                                 forceOverwrite = TRUE,
                                                 quiet = TRUE),
                                outFile),
                   "overwriting")
    expect_equal(generateQCReport(se = se, outFile = outFile2,
                                  forceOverwrite = FALSE, quiet = TRUE),
                 outFile2)
})