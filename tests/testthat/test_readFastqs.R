context("readFastqs")

test_that("readFastqs works as expected for trans experiments", {
  fnameR1 <- system.file("extdata", "transInput_1.fastq.gz", package = "mutscan")
  fnameR2 <- system.file("extdata", "transInput_2.fastq.gz", package = "mutscan")
  
  transInput <- readFastqs(experimentType = "trans", fastqForward = fnameR1, fastqReverse = fnameR2,
                           skipForward = 1, skipReverse = 1, umiLengthForward = 10, 
                           umiLengthReverse = 8, constantLengthForward = 18,
                           constantLengthReverse = 20, variableLengthForward = 96,
                           variableLengthReverse = 96, adapterForward = "GGAAGAGCACACGTC", 
                           adapterReverse = "GGAAGAGCGTCGTGT", verbose = FALSE)
  s1 <- ShortRead::readFastq(fnameR1)
  
  expect_true(isValidL(transInput))
  expect_identical(length(s1), colData(transInput)[1, "totalNbrReadPairs"])
  expect_true(length(s1) == length(assay(transInput, "umis")$seq) + colData(transInput)[1, "nbrReadPairsWithAdapter"])
  expect_true(all(c("umis", "constantSeqForward", "constantSeqReverse", "variableSeqForward", 
                    "variableSeqReverse") %in% SummarizedExperiment::assayNames(transInput)))
})

test_that("readFastqs works as expected for cis experiments", {
  fnameR1 <- system.file("extdata", "cisInput_1.fastq.gz", package = "mutscan")
  fnameR2 <- system.file("extdata", "cisInput_2.fastq.gz", package = "mutscan")
  
  cisInput <- readFastqs(experimentType = "cis", fastqForward = fnameR1, fastqReverse = fnameR2,
                         skipForward = 1, skipReverse = 1, umiLengthForward = 10, 
                         umiLengthReverse = 7, constantLengthForward = 18,
                         constantLengthReverse = 17, variableLengthForward = 96,
                         variableLengthReverse = 96, adapterForward = "GGAAGAGCACACGTC", 
                         adapterReverse = "GGAAGAGCGTCGTGT", verbose = FALSE)
  s1 <- ShortRead::readFastq(fnameR1)
  
  expect_true(isValidL(cisInput))
  expect_identical(length(s1), colData(cisInput)[1, "totalNbrReadPairs"])
  expect_true(length(s1) == length(assay(cisInput, "umis")$seq) + colData(cisInput)[1, "nbrReadPairsWithAdapter"])
  expect_true(all(c("umis", "constantSeqForward", "constantSeqReverse", "variableSeqForward") %in% 
                    SummarizedExperiment::assayNames(cisInput)))
  expect_false("variableSeqReverse" %in% assaynames(cisInput))
})

test_that("readFastqs fails when the wrong inputs are given", {
  fnameR1 <- system.file("extdata", "transInput_1.fastq.gz", package = "mutscan")
  fnameR2 <- system.file("extdata", "transInput_2.fastq.gz", package = "mutscan")
  expect_error(readFastqs(experimentType = "unknown", fastqForward = fnameR1, fastqReverse = fnameR2))
  expect_error(readFastqs(experimentType = "trans", fastqForward = paste0(fnameR1, "_unknown"),
                          fastqReverse = fnameR2))
  expect_error(readFastqs(experimentType = "trans", fastqForward = fnameR1, fastqReverse = fnameR2,
                          skipForward = "1"))
  expect_error(readFastqs(experimentType = "trans", fastqForward = fnameR1, fastqReverse = fnameR2,
                          adapterForward = "EE "))
})