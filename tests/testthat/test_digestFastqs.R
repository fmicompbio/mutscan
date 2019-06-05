context("digestFastqs - helpers")
test_that("compareCodonPositions works", {
  expect_true(compareCodonPositions("f.1.AAA_", "f.2.ACT_", "."))
  expect_true(compareCodonPositions("f.2.ACT_", "f.10.TCA_", "."))
})

context("digestFastqs - inputs")
test_that("digestFastqs fails with incorrect arguments", {
  ## example "trans" fastq files 
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2, 
    mergeForwardReverse = FALSE, revComplForward = FALSE, revComplReverse = FALSE,
    skipForward = 1, skipReverse = 1, 
    umiLengthForward = 10, umiLengthReverse = 8, 
    constantLengthForward = 18, constantLengthReverse = 20, 
    variableLengthForward = 96, variableLengthReverse = 96,
    adapterForward = "GGAAGAGCACACGTC", 
    adapterReverse = "GGAAGAGCGTCGTGT",
    primerForward = "",
    primerReverse = "",
    wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
    wildTypeReverse = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT", 
    constantForward = "AACCGGAGGAGGGAGCTG", 
    constantReverse = "GAAAAAGGAAGCTGGAGAGA", 
    avePhredMinForward = 20.0, avePhredMinReverse = 20.0, 
    variableNMaxForward = 0, variableNMaxReverse = 0, 
    umiNMax = 0,
    nbrMutatedCodonsMaxForward = 1,
    nbrMutatedCodonsMaxReverse = 1,
    forbiddenMutatedCodonsForward = "NNW",
    forbiddenMutatedCodonsReverse = "NNW",
    mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".",
    verbose = FALSE
  )
  
  ## Incorrect specification of merging/rev complementing arguments
  for (var in c("mergeForwardReverse", "revComplForward", "revComplReverse")) {
    L <- Ldef; L[[var]] <- "str"
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- c(TRUE, TRUE)
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- c(TRUE, FALSE)
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- 1
    expect_error(do.call(digestFastqs, L))
  }
  
  ## Nonexistent fastqForward/fastqReverse file
  L <- Ldef; L$fastqForward <- "nonexistent.fastq.gz"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqReverse <- "nonexistent.fastq.gz"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqForward <- c(fqt1, fqt2)
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqReverse <- c(fqt1, fqt2)
  expect_error(do.call(digestFastqs, L))
  
  ## Wrong type of numeric argument
  for (var in c("skipForward", "skipReverse", "umiLengthForward", 
                "umiLengthReverse", "constantLengthForward", "constantLengthReverse",
                "variableLengthForward", "variableLengthReverse", 
                "avePhredMinForward", "avePhredMinReverse", "variableNMaxForward",
                "variableNMaxReverse", "umiNMax", 
                "nbrMutatedCodonsMaxForward", "nbrMutatedCodonsMaxReverse", 
                "mutatedPhredMinForward", "mutatedPhredMinReverse")) {
    L <- Ldef; L[[var]] <- "str"
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- c(1, 2)
    expect_error(do.call(digestFastqs, L))
    if (!(var %in% c("umiLengthForward", "umiLengthReverse", "constantLengthForward",
                     "constantLengthReverse", "variableLengthForward", "variableLengthReverse"))) {
      L <- Ldef; L[[var]] <- -1
      expect_error(do.call(digestFastqs, L))
    }
    L <- Ldef; L[[var]] <- TRUE
    expect_error(do.call(digestFastqs, L))
  }
  
  ## Either all or none of the sequence lengths must be specified
  for (var in c("skipForward", "skipReverse", "umiLengthForward",
                "umiLengthReverse", "constantLengthForward", "constantLengthReverse")) {
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
    L <- Ldef; L[[var]] <- paste0(L[[var]], " ")
    expect_error(do.call(digestFastqs, L))
  }
  
  ## Wild type sequence not in named vector (or unnamed string)
  L <- Ldef
  L$wildTypeForward <- c("ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA", 
                         "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT")
  expect_error(do.call(digestFastqs, L))
  
  ## Wild type sequence of wrong length
  for (var in c("wildTypeForward", "wildTypeReverse")) {
    L <- Ldef; L[[var]] <- substr(L[[var]], 1, 10)
    expect_error(do.call(digestFastqs, L))
  }
  for (var in c("variableLengthForward", "variableLengthReverse")) {
    L <- Ldef; L[[var]] <- 10
    expect_error(do.call(digestFastqs, L))
  }
  L <- Ldef; L[["mergeForwardReverse"]] <- TRUE
  expect_warning(do.call(digestFastqs, L))
  
  ## Constant sequence of wrong length
  for (var in c("constantForward", "constantReverse")) {
    L <- Ldef; L[[var]] <- substr(L[[var]], 1, 10)
    expect_error(do.call(digestFastqs, L))
  }
  for (var in c("constantLengthForward", "constantLengthReverse")) {
    L <- Ldef; L[[var]] <- 10
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- -1
    expect_error(do.call(digestFastqs, L))
  }
  
  ## mergeForwardReverse=TRUE requires equal-length fwd/rev variable sequences
  L <- Ldef; L[["mergeForwardReverse"]] <- TRUE; L[["variableLengthForward"]] <- -1; L[["wildTypeReverse"]] <- ""
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["mergeForwardReverse"]] <- TRUE; L[["variableLengthReverse"]] <- -1; L[["wildTypeReverse"]] <- ""
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["mergeForwardReverse"]] <- TRUE; L[["variableLengthForward"]] <- 10; L[["wildTypeReverse"]] <- ""
  L[["wildTypeForward"]] <- substr(L[["wildTypeForward"]], 1, 10)
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["mergeForwardReverse"]] <- TRUE; L[["variableLengthReverse"]] <- 10; L[["wildTypeReverse"]] <- ""
  L[["wildTypeReverse"]] <- substr(L[["wildTypeReverse"]], 1, 10)
  expect_error(do.call(digestFastqs, L))
  
  ## Invalid forbidden codons
  L <- Ldef; L[["forbiddenMutatedCodonsForward"]] <- c("EFI")
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["forbiddenMutatedCodonsReverse"]] <- c("EFI")
  expect_error(do.call(digestFastqs, L))
  
  ## Invalid mutNameDelimiter
  L <- Ldef; L[["mutNameDelimiter"]] <- ""
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["mutNameDelimiter"]] <- "__"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["mutNameDelimiter"]] <- "f"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["mutNameDelimiter"]] <- "r"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["mutNameDelimiter"]] <- 1
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["mutNameDelimiter"]] <- c("a", "B")
  expect_error(do.call(digestFastqs, L))
  
  ## Invalid value of verbose
  L <- Ldef; L[["verbose"]] <- 2
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["verbose"]] <- "TRUE"
  expect_error(do.call(digestFastqs, L))
})

context("digestFastqs - trans")
test_that("digestFastqs works as expected for trans experiments", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2, 
    mergeForwardReverse = FALSE, revComplForward = FALSE, revComplReverse = FALSE,
    skipForward = 1, skipReverse = 1, 
    umiLengthForward = 10, umiLengthReverse = 8, 
    constantLengthForward = 18, constantLengthReverse = 20, 
    variableLengthForward = 96, variableLengthReverse = 96,
    adapterForward = "GGAAGAGCACACGTC", 
    adapterReverse = "GGAAGAGCGTCGTGT",
    primerForward = "",
    primerReverse = "",
    wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
    wildTypeReverse = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT", 
    constantForward = "AACCGGAGGAGGGAGCTG", 
    constantReverse = "GAAAAAGGAAGCTGGAGAGA", 
    avePhredMinForward = 20.0, avePhredMinReverse = 20.0,
    variableNMaxForward = 0, variableNMaxReverse = 0, 
    umiNMax = 0,
    nbrMutatedCodonsMaxForward = 1,
    nbrMutatedCodonsMaxReverse = 1,
    forbiddenMutatedCodonsForward = "NNW",
    forbiddenMutatedCodonsReverse = "NNW",
    mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".",
    verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f3_nTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f4_nTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f5_nMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f6_nTooManyMutCodons, 287L + 105L)
  expect_equal(res$filterSummary$f7_nForbiddenCodons, 6L + 2L)
  expect_equal(res$filterSummary$nbrRetained, 279L)
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(sum(res$summaryTable$nbrReads == 2), 2L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 2]),
               sort(c("f.13.GAG", "f.26.TAG_r.8.AGG")))
  expect_true(all(res$summaryTable$nbrReads == res$summaryTable$nbrUmis))
  
  ## Check that mutant naming worked (compare to manual matching)
  example_seq <- paste0("ACTGATACAACCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTG", 
                        "CAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA_ATCGCCCGGCT", 
                        "GGAGGAAAAAGTGGGCACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGC", 
                        "CAACATGCTCAGGGAACAGGTGGCACAGCTT")
  expect_equal(res$summaryTable$mutantName[res$summaryTable$sequence == example_seq], 
               "f.4.ACC_r.9.GGC")
  
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 14], 160L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 22], 52L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 27], 302L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 33], 472L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 37], 4020L)
  
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 14], 11L)
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 27], 4L)
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 37], 1L)

  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 14], 204L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 22], 14L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 27], 474L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 33], 462L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 37], 4405L)
  
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 14], 17L)
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 27], 3L)
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 33], 1L)
  
})

context("digestFastqs - trans - do not specify variable sequence length")
test_that("digestFastqs works as expected for trans experiments, when variable sequence length is not given", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2, 
    mergeForwardReverse = FALSE, revComplForward = FALSE, revComplReverse = FALSE,
    skipForward = 1, skipReverse = 1, 
    umiLengthForward = 10, umiLengthReverse = 8, 
    constantLengthForward = 18, constantLengthReverse = 20, 
    variableLengthForward = -1, variableLengthReverse = -1,
    adapterForward = "GGAAGAGCACACGTC", 
    adapterReverse = "GGAAGAGCGTCGTGT",
    primerForward = "",
    primerReverse = "",
    wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
    wildTypeReverse = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT", 
    constantForward = "AACCGGAGGAGGGAGCTG", 
    constantReverse = "GAAAAAGGAAGCTGGAGAGA", 
    avePhredMinForward = 20.0, avePhredMinReverse = 20.0,
    variableNMaxForward = 0, variableNMaxReverse = 0, 
    umiNMax = 0,
    nbrMutatedCodonsMaxForward = 1,
    nbrMutatedCodonsMaxReverse = 1,
    forbiddenMutatedCodonsForward = "NNW",
    forbiddenMutatedCodonsReverse = "NNW",
    mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".",
    verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f3_nTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f4_nTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f5_nMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f6_nTooManyMutCodons, 287L + 105L)
  expect_equal(res$filterSummary$f7_nForbiddenCodons, 6L + 2L)
  expect_equal(res$filterSummary$nbrRetained, 279L)
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(sum(res$summaryTable$nbrReads == 2), 2L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 2]),
               sort(c("f.13.GAG", "f.26.TAG_r.8.AGG")))
  expect_true(all(res$summaryTable$nbrReads == res$summaryTable$nbrUmis))
  
  ## Check that mutant naming worked (compare to manual matching)
  example_seq <- paste0("ACTGATACAACCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTG", 
                        "CAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA_ATCGCCCGGCT", 
                        "GGAGGAAAAAGTGGGCACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGC", 
                        "CAACATGCTCAGGGAACAGGTGGCACAGCTT")
  expect_equal(res$summaryTable$mutantName[res$summaryTable$sequence == example_seq], 
               "f.4.ACC_r.9.GGC")
  
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 14], 160L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 22], 52L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 27], 302L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 33], 472L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 37], 4020L)
  
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 14], 11L)
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 27], 4L)
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 37], 1L)
  
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 14], 204L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 22], 14L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 27], 474L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 33], 462L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 37], 4405L)
  
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 14], 17L)
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 27], 3L)
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 33], 1L)
  
})

context("digestFastqs - trans - multiple reference sequences")
test_that("digestFastqs works as expected for trans experiments when multiple reference sequences are provided", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2, 
    mergeForwardReverse = FALSE, revComplForward = FALSE, revComplReverse = FALSE,
    skipForward = 1, skipReverse = 1, 
    umiLengthForward = 10, umiLengthReverse = 8, 
    constantLengthForward = 18, constantLengthReverse = 20, 
    variableLengthForward = 96, variableLengthReverse = 96,
    adapterForward = "GGAAGAGCACACGTC", 
    adapterReverse = "GGAAGAGCGTCGTGT",
    primerForward = "",
    primerReverse = "",
    wildTypeForward = c(forward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",  ## this is the right one
                        wrong = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT"),
    wildTypeReverse = c(wrong = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
                        reverse1 = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT"), ## this is the right one
    constantForward = "AACCGGAGGAGGGAGCTG", 
    constantReverse = "GAAAAAGGAAGCTGGAGAGA", 
    avePhredMinForward = 20.0, avePhredMinReverse = 30.0,
    variableNMaxForward = 0, variableNMaxReverse = 0, 
    umiNMax = 0,
    nbrMutatedCodonsMaxForward = 1,
    nbrMutatedCodonsMaxReverse = 1,
    forbiddenMutatedCodonsForward = "NNA",
    forbiddenMutatedCodonsReverse = "NNW",
    mutatedPhredMinForward = 25.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = "=",
    verbose = TRUE
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nAvgVarQualTooLow, 88L)
  expect_equal(res$filterSummary$f3_nTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f4_nTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f5_nMutQualTooLow, 333L)
  expect_equal(res$filterSummary$f6_nTooManyMutCodons, 10L + 59L)
  expect_equal(res$filterSummary$f7_nForbiddenCodons, 1L + 2L)
  expect_equal(res$filterSummary$nbrRetained, 193L)
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(sum(res$summaryTable$nbrReads == 2), 1L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 2]),
               sort(c("forward=26=TAG_reverse1=8=AGG")))
  expect_true(all(res$summaryTable$nbrReads == res$summaryTable$nbrUmis))
  
  ## Check that mutant naming worked (compare to manual matching)
  example_seq <- paste0("ACTGATACAACCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTG", 
                        "CAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA_ATCGCCCGGCT", 
                        "GGAGGAAAAAGTGGGCACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGC", 
                        "CAACATGCTCAGGGAACAGGTGGCACAGCTT")
  expect_equal(res$summaryTable$mutantName[res$summaryTable$sequence == example_seq], 
               "forward=4=ACC_reverse1=9=GGC")
  
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 14], 96L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 22], 26L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 27], 179L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 33], 302L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 37], 2865L)
  
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 14], 2L)
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 27], 3L)
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 37], 1L)
  
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 14], 116L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 22], 5L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 27], 308L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 33], 322L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 37], 3101L)
  
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 14], 6L)
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 27], 1L)
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 33], 1L)
})

context("digestFastqs - cis")
test_that("digestFastqs works as expected for cis experiments", {
  fqc1 <- system.file("extdata/cisInput_1.fastq.gz", package = "mutscan")
  fqc2 <- system.file("extdata/cisInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqc1, fastqReverse = fqc2, 
    mergeForwardReverse = TRUE, revComplForward = FALSE, revComplReverse = TRUE,
    skipForward = 1, skipReverse = 1, 
    umiLengthForward = 10, umiLengthReverse = 7, 
    constantLengthForward = 18, constantLengthReverse = 17, 
    variableLengthForward = 96, variableLengthReverse = 96,
    adapterForward = "GGAAGAGCACACGTC", 
    adapterReverse = "GGAAGAGCGTCGTGT",
    primerForward = "",
    primerReverse = "",
    wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
    wildTypeReverse = "", 
    constantForward = "AACCGGAGGAGGGAGCTG", 
    constantReverse = "GAGTTCATCCTGGCAGC", 
    avePhredMinForward = 20.0, avePhredMinReverse = 20.0,
    variableNMaxForward = 0, variableNMaxReverse = 0, 
    umiNMax = 0,
    nbrMutatedCodonsMaxForward = 1,
    nbrMutatedCodonsMaxReverse = 1,
    forbiddenMutatedCodonsForward = "NNW",
    forbiddenMutatedCodonsReverse = "NNW",
    mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".",
    verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 126L)
  expect_equal(res$filterSummary$f2_nAvgVarQualTooLow, 0L)
  expect_equal(res$filterSummary$f3_nTooManyNinVar, 44L)
  expect_equal(res$filterSummary$f4_nTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f5_nMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f6_nTooManyMutCodons, 581L)
  expect_equal(res$filterSummary$f7_nForbiddenCodons, 82L)
  expect_equal(res$filterSummary$nbrRetained, 167)
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(sum(res$summaryTable$nbrReads == 2), 11L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 2]),
               sort(c("f.14.AAG", "f.15.ATG", "f.19.CAC", "f.1.ACC", "f.20.AAC", "f.21.GTG",
                      "f.24.AGC", "f.4.CGC", "f.7.GTG", "f.9.GGC", "f.9.GTC")))
  expect_true(all(res$summaryTable$nbrReads == res$summaryTable$nbrUmis))
  
  ## Check that mutant naming worked (compare to manual matching)
  example_seq <- paste0("ACTGATACACGCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCT", 
                        "GCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA")
  expect_equal(res$summaryTable$mutantName[res$summaryTable$sequence == example_seq], 
               "f.4.CGC")
  
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 14], 51L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 22], 19L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 27], 79L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 33], 163L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 37], 2690L)
  
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 14], 1L)
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 27], 1L)
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 37], 2L)
  
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 14], 41L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 27], 56L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 33], 96L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 37], 2620L)
  
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 2], 4L)
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 14], 7L)
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 33], 1L)
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 37], 14L)
  
})

