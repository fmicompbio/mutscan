context("readFastqsTrans")

test_that("readFastqsTrans works as expected", {
  fnameR1 <- system.file("extdata", "transInput_1.fastq.gz", package = "mutscan")
  fnameR2 <- system.file("extdata", "transInput_2.fastq.gz", package = "mutscan")
  
  transInput <- readFastqsTrans(fastqForward = fnameR1, fastqReverse = fnameR2,
                                skipForward = 1, skipReverse = 1, umiLengthForward = 10, 
                                umiLengthReverse = 8, constantLengthForward = 18,
                                constantLengthReverse = 20, adapterForward = "GGAAGAGCACACGTC", 
                                adapterReverse = "GGAAGAGCGTCGTGT", verbose = FALSE)
  s1 <- ShortRead::readFastq(fnameR1)
  
  expect_true(isValidL(transInput))
  expect_identical(length(s1), transInput$readSummary[1,"totalNbrReadPairs"])
  expect_true(length(s1) == length(transInput$umis) + transInput$readSummary[1,"nbrReadPairsWithAdapter"])
})

