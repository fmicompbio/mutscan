context("digestFastqs - helpers")
test_that("compareCodonPositions works", {
  expect_true(compareCodonPositions("f1AAA_", "f2ACT_"))
  expect_true(compareCodonPositions("f2ACT_", "f10TCA_"))
})

context("digestFastqs - inputs")
test_that("digestFastqs fails with incorrect arguments", {
  ## example "trans" fastq files 
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    experimentType = "trans", fastqForward = fqt1, 
    fastqReverse = fqt2, skipForward = 1, skipReverse = 1, 
    umiLengthForward = 10, umiLengthReverse = 8, 
    constantLengthForward = 18, constantLengthReverse = 20, 
    variableLengthForward = 96, variableLengthReverse = 96,
    adapterForward = "GGAAGAGCACACGTC", 
    adapterReverse = "GGAAGAGCGTCGTGT",
    wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
    wildTypeReverse = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT", 
    constantForward = "AACCGGAGGAGGGAGCTG", 
    constantReverse = "GAAAAAGGAAGCTGGAGAGA", 
    avePhredMin = 20.0, variableNMax = 0, umiNMax = 0,
    nbrMutatedCodonsMax = 1,
    forbiddenMutatedCodons = "NNW",
    mutatedPhredMin = 0.0,
    verbose = FALSE
  )
  
  ## Incorrect experimentType
  L <- Ldef; L$experimentType <- "unknown"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$experimentType <- c("cis", "trans")
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$experimentType <- 3
  expect_error(do.call(digestFastqs, L))

  ## Nonexistent fastqForward/fastqReverse file
  L <- Ldef; L$fastqForward <- "nonexistent.fastq.gz"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqReverse <- "nonexistent.fastq.gz"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqForward <- c(fqt1, fqt2)
  expect_error(do.call(digestFastqs, L))
  
  ## Wrong type of numeric argument
  for (var in c("skipForward", "skipReverse", "umiLengthForward", 
                "umiLengthRevese", "constantLengthForward", "constantLengthReverse",
                "variableLengthForward", "variableLengthReverse", 
                "avePhredMin", "variableNMax", "umiNMax", 
                "nbrMutatedCodonsMax", "mutatedPhredMin")) {
    L <- Ldef; L[[var]] <- "str"
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- c(1, 2)
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- -1
    expect_error(do.call(digestFastqs, L))
  }
  
  ## Invalid sequences
  for (var in c("adapterForward", "adapterReverse", "wildTypeForward",
                "wildTypeReverse", "constantForward", "constantReverse")) {
    L <- Ldef; L[[var]] <- c("ACGT", "ACGT")
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- 1
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- "EF"
    expect_error(do.call(digestFastqs, L))
  }
  
  ## Wild type sequence of wrong length
  for (var in c("wildTypeForward", "wildTypeReverse")) {
    L <- Ldef; L[[var]] <- substr(L[[var]], 1, 10)
    expect_error(do.call(digestFastqs, L))
  }
  for (var in c("variableLengthForward", "variableLengthReverse")) {
    L <- Ldef; L[[var]] <- 10
    expect_error(do.call(digestFastqs, L))
  }
  L <- Ldef; L[["experimentType"]] <- "cis"
  expect_warning(do.call(digestFastqs, L))
  
  ## Constant sequence of wrong length
  for (var in c("constantForward", "constantReverse")) {
    L <- Ldef; L[[var]] <- substr(L[[var]], 1, 10)
    expect_error(do.call(digestFastqs, L))
  }
  for (var in c("constantLengthForward", "constantLengthReverse")) {
    L <- Ldef; L[[var]] <- 10
    expect_error(do.call(digestFastqs, L))
  }
  
  ## Invalid forbidden codons
  L <- Ldef; L[["forbiddenMutatedCodons"]] <- c("EFI")
  expect_error(do.call(digestFastqs, L))
  
  ## Invalid value of verbose
  L <- Ldef; L[["verbose"]] <- 2
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["verbose"]] <- "TRUE"
  expect_error(do.call(digestFastqs, L))
})

test_that("readFastqs works as expected for trans experiments", {
  fnameR1 <- system.file("extdata", "transInput_1.fastq.gz", package = "mutscan")
  fnameR2 <- system.file("extdata", "transInput_2.fastq.gz", package = "mutscan")
  
  transInput <- readFastqs(experimentType = "trans", fastqForward = fnameR1, fastqReverse = fnameR2,
                           skipForward = 1, skipReverse = 1, umiLengthForward = 10, 
                           umiLengthReverse = 8, constantLengthForward = 18,
                           constantLengthReverse = 20, variableLengthForward = 96,
                           variableLengthReverse = 96, adapterForward = "GGAAGAGCACACGTC", 
                           adapterReverse = "GGAAGAGCGTCGTGT", maxChunkSize = 100, verbose = FALSE)
  s1 <- ShortRead::readFastq(fnameR1)
  
  expect_true(isValidL(transInput))
  expect_identical(length(s1), SummarizedExperiment::colData(transInput)[1, "totalNbrReadPairs"])
  expect_true(length(s1) == length(assay(transInput, "umis")$seq) + 
                SummarizedExperiment::colData(transInput)[1, "nbrReadPairsWithAdapter"])
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
  expect_identical(length(s1), SummarizedExperiment::colData(cisInput)[1, "totalNbrReadPairs"])
  expect_true(length(s1) == length(assay(cisInput, "umis")$seq) + 
                SummarizedExperiment::colData(cisInput)[1, "nbrReadPairsWithAdapter"])
  expect_true(all(c("umis", "constantSeqForward", "constantSeqReverse", "variableSeqForward") %in% 
                    SummarizedExperiment::assayNames(cisInput)))
  expect_false("variableSeqReverse" %in% assayNames(cisInput))
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