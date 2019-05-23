context("calcErrorStats")

test_that("calcErrorStats works as expected on trans experiments", {
  fnameR1 <- system.file("extdata", "transInput_1.fastq.gz", package = "mutscan")
  fnameR2 <- system.file("extdata", "transInput_2.fastq.gz", package = "mutscan")
  
  transInput <- readFastqs(experimentType = "trans", fastqForward = fnameR1, fastqReverse = fnameR2,
                           skipForward = 1, skipReverse = 1, umiLengthForward = 10, 
                           umiLengthReverse = 8, constantLengthForward = 18,
                           constantLengthReverse = 20, adapterForward = "GGAAGAGCACACGTC", 
                           adapterReverse = "GGAAGAGCGTCGTGT", verbose = FALSE)
  expect_message(errStats <- calcErrorStats(transInput, verbose = TRUE))
  expect_is(errStats, "list")
  expect_length(errStats, 6L)
  expect_identical(sum(errStats$QcountConstantErrorsF + errStats$QcountConstantCorrectF),
                   length(assay(transInput, "constantSeqForward")$seq) * Biostrings::width(assay(transInput, "constantSeqForward")$seq)[1])
  expect_equal(errStats$propErrorsConstantF, 0.03336573)
  expect_equal(errStats$propErrorsConstantR, 0.04635569)
})

