context("filterReads")

test_that("filterReads works as expected for trans experiments", {
  fnameR1 <- system.file("extdata", "transInput_1.fastq.gz", package = "mutscan")
  fnameR2 <- system.file("extdata", "transInput_2.fastq.gz", package = "mutscan")
  
  transInput <- readFastqs(experimentType = "trans", fastqForward = fnameR1, fastqReverse = fnameR2,
                           skipForward = 1, skipReverse = 1, umiLengthForward = 10, 
                           umiLengthReverse = 8, constantLengthForward = 18,
                           constantLengthReverse = 20, adapterForward = "GGAAGAGCACACGTC", 
                           adapterReverse = "GGAAGAGCGTCGTGT", verbose = FALSE)
  transInputFiltered1 <- filterReads(transInput, avePhredMin = 20, 
                                     variableNMax = 0, umiNMax = 0, 
                                     wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
                                     wildTypeReverse = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT", 
                                     nbrMutatedCodonsMax = 1, forbiddenMutatedCodon = "NNW", 
                                     maxChunkSize = 100)
  transInputFiltered2 <- filterReads(transInput, avePhredMin = 20, 
                                     variableNMax = 0, umiNMax = 0, 
                                     wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
                                     wildTypeReverse = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT", 
                                     nbrMutatedCodonsMax = 2, maxChunkSize = 100)
  expect_true(isValidL(transInputFiltered1))
  expect_true(isValidL(transInputFiltered2))
  expect_identical(colData(transInputFiltered1)[1, "totalNbrReadPairsPassedFilters"], 279L)
  expect_identical(colData(transInputFiltered2)[1, "totalNbrReadPairsPassedFilters"], 489L)
})

