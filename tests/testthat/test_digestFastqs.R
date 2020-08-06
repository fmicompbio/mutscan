context("digestFastqs - helpers")
test_that("compareCodonPositions works", {
  expect_true(compareCodonPositions("f.1.AAA_", "f.2.ACT_", "."))
  expect_true(compareCodonPositions("f.2.ACT_", "f.10.TCA_", "."))
})

test_that("findClosestRefSeq works", {
  expect_equal(findClosestRefSeq(varSeq = "ACGT", wtSeq = c("ATTT", "ACTT")), 1L)
  expect_equal(findClosestRefSeq(varSeq = "ACGT", wtSeq = c("ACGT", "ACTT")), 0L)
  expect_equal(findClosestRefSeq(varSeq = "ACGT", wtSeq = c("ACG", "ACGT")), 1L)
  expect_equal(findClosestRefSeq(varSeq = "ACGT", wtSeq = c("AACGT", "ACCTA")), 1L)
})

test_that("mergeReadPairsPartial works", {
  ## don't count N as mismatch, pick base with higher quality
  sF1 <- "AAAANA"; qF1 <- rep(40L, nchar(sF1))
  sR1 <- "AACCCC"; qR1 <- rep(42L, nchar(sR1))
  res1a <- mutscan:::test_mergeReadPairPartial(sF1, qF1, sR1, qR1, 1, 6, 1/4, TRUE)
  res1b <- mutscan:::test_mergeReadPairPartial(sF1, qF1, sR1, qR1, 1, 6, 1/4, FALSE)
  res1c <- mutscan:::test_mergeReadPairPartial(sF1, qR1, sR1, qF1, 2, 2, 1/4, TRUE)
  res1d <- mutscan:::test_mergeReadPairPartial(sF1, qF1, sR1, qR1, 2, 2, 0, TRUE)
  expect_is(res1a, "list")
  expect_is(res1b, "list")
  expect_is(res1c, "list")
  expect_is(res1d, "list")
  expect_identical(res1a$mergedSeq, "AAAACCCC")
  expect_identical(res1a$mergedQual, rep(c(40L, 42L), c(2, 6)))
  expect_identical(res1a, res1b)
  expect_identical(res1c$mergedSeq, "AAAANACCCC")
  expect_identical(res1d$mergedSeq, "AAAAAACCCC")
  
  ## no valid overlap, multiple possible overlaps with single valid overlap
  sF2 <- "TTACACG"; qF2 <- rep(10L, nchar(sF2))
  sR2 <- "ACACACA"; qR2 <- rep(40L, nchar(sR2))
  res2a <- mutscan:::test_mergeReadPairPartial(sF2, qF2, sR2, qR2, maxFracMismatchOverlap = 3/7, TRUE)
  res2b <- mutscan:::test_mergeReadPairPartial(sF2, qF2, sR2, qR2, 1, 7, 1/5, TRUE)
  expect_is(res2a, "list")
  expect_is(res2b, "list")
  expect_identical(res2a$mergedSeq, sR2)
  expect_identical(res2a$mergedQual, qR2)
  expect_identical(res2b$mergedSeq, "TTACACACA")
  expect_identical(res2b$mergedQual, rep(c(10L,40L), c(2,7)))
  
  ## padded reads
  for (i in 1:10) {
    sF <- paste(rep(c("C","A"), c(i, 6)), collapse = "")
    qF <- rep(30L, nchar(sF))
    sR <- paste(rep(c("A","C"), c(6, i)), collapse = "")
    qR <- rep(30L, nchar(sR))
    res <- mutscan:::test_mergeReadPairPartial(sF, qF, sR, qR, 6, 6, 0)
    expect_identical(res$mergedSeq, paste(rep(c("C","A","C"), c(i, 6, i)), collapse = ""))
  }

  ## reads of unequal length
  sR <- paste(rep("A", 6), collapse = "")
  qR <- rep(30L, nchar(sR))
  for (i in 1:10) {
    sF <- paste(rep(c("C","A"), c(i, 6)), collapse = "")
    qF <- rep(30L, nchar(sF))
    res <- mutscan:::test_mergeReadPairPartial(sF, qF, sR, qR, 6, 6, 0)
    expect_identical(res$mergedSeq, sF)
    res <- mutscan:::test_mergeReadPairPartial(sR, qR, sF, qF, 6, 6, 0)
    expect_identical(res$mergedSeq, sR)
  }
})


context("digestFastqs - inputs")
test_that("digestFastqs fails with incorrect arguments", {
  ## example "trans" fastq files 
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2, 
    mergeForwardReverse = FALSE, 
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
    revComplForward = FALSE, revComplReverse = FALSE,
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
    maxNReads = -1, verbose = FALSE
  )
  
  ## Incorrect specification of merging/rev complementing arguments
  for (var in c("mergeForwardReverse", "revComplForward", "revComplReverse")) {
    L <- Ldef; L[[var]] <- "str"
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- ""
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
  L <- Ldef; L$fastqForward <- ""
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqForward <- NULL
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqReverse <- "nonexistent.fastq.gz"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqReverse <- ""
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqForward <- c(fqt1, fqt2)
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqReverse <- c(fqt1, fqt2)
  expect_error(do.call(digestFastqs, L))
  
  ## No reverse reads if mergeForwardReverse is TRUE
  L <- Ldef; L[["fastqReverse"]] <- NULL; L[["mergeForwardReverse"]] <- TRUE
  expect_error(do.call(digestFastqs, L))
  
  ## Wrong type of numeric argument
  for (var in c("skipForward", "skipReverse", "umiLengthForward", 
                "umiLengthReverse", "constantLengthForward", "constantLengthReverse",
                "variableLengthForward", "variableLengthReverse", 
                "avePhredMinForward", "avePhredMinReverse", "variableNMaxForward",
                "variableNMaxReverse", "umiNMax", 
                "nbrMutatedCodonsMaxForward", "nbrMutatedCodonsMaxReverse", 
                "mutatedPhredMinForward", "mutatedPhredMinReverse",
                "variableCollapseMaxDist", "umiCollapseMaxDist", "variableCollapseMinReads")) {
    L <- Ldef; L[[var]] <- "str"
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- c(1, 2)
    expect_error(do.call(digestFastqs, L))
    if (!(var %in% c("skipForward", "skipReverse", "umiLengthForward",
                     "umiLengthReverse", "constantLengthForward",
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
  
  ## If sequence part lengths are given, primers can not be (and opposite)
  for (var in c("skip", "umiLength", "constantLength")) {
    for (di in c("Forward", "Reverse")) {
      L <- Ldef; L[[paste0("primer", di)]] <- "ACGT"
      expect_error(do.call(digestFastqs, L))
      L <- Ldef; L[[paste0(var, di)]] <- -1; L[[paste0("primer", di)]] <- ""
      expect_error(do.call(digestFastqs, L))
    }
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
    mergeForwardReverse = FALSE, 
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
    revComplForward = FALSE, revComplReverse = FALSE,
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
    maxNReads = -1, verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)

  L <- Ldef; L$fastqForward <- c(fqt1, fqt1); L$fastqReverse <- c(fqt2, fqt2)
  res2 <- do.call(digestFastqs, L)
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadTooShort, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f9_nbrTooManyMutCodons, 287L + 105L)
  expect_equal(res$filterSummary$f10_nbrForbiddenCodons, 6L + 2L)
  expect_equal(res$filterSummary$nbrRetained, 279L)
  
  expect_equal(res$filterSummary * 2, res2$filterSummary)
  expect_equal(res$summaryTable$nbrReads * 2, res2$summaryTable$nbrReads)
  expect_equal(res$summaryTable$nbrUmis, res2$summaryTable$nbrUmis)
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(res$summaryTable$nbrReads, res$summaryTable$maxNbrReads)
  expect_equal(sum(res$summaryTable$nbrReads == 2), 2L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 2]),
               sort(c("f.13.GAG_r.0.WT", "f.26.TAG_r.8.AGG")))
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
    mergeForwardReverse = FALSE, 
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
    revComplForward = FALSE, revComplReverse = FALSE,
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
    maxNReads = -1, verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadTooShort, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f9_nbrTooManyMutCodons, 287L + 105L)
  expect_equal(res$filterSummary$f10_nbrForbiddenCodons, 6L + 2L)
  expect_equal(res$filterSummary$nbrRetained, 279L)
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(res$summaryTable$nbrReads, res$summaryTable$maxNbrReads)
  expect_equal(sum(res$summaryTable$nbrReads == 2), 2L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 2]),
               sort(c("f.13.GAG_r.0.WT", "f.26.TAG_r.8.AGG")))
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
    mergeForwardReverse = FALSE, 
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
    revComplForward = FALSE, revComplReverse = FALSE,
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
    maxNReads = -1, verbose = TRUE
  )
  
  expect_output(res <- do.call(digestFastqs, Ldef))
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadTooShort, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 88L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrMutQualTooLow, 333L)
  expect_equal(res$filterSummary$f9_nbrTooManyMutCodons, 10L + 59L)
  expect_equal(res$filterSummary$f10_nbrForbiddenCodons, 1L + 2L)
  expect_equal(res$filterSummary$nbrRetained, 193L)
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(res$summaryTable$nbrReads, res$summaryTable$maxNbrReads)
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

context("digestFastqs - only forward read given")
test_that("digestFastqs works as expected for experiments with only forward read", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = NULL, 
    mergeForwardReverse = FALSE, 
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
    revComplForward = FALSE, revComplReverse = FALSE,
    skipForward = 1, skipReverse = 1, 
    umiLengthForward = 10, umiLengthReverse = 8, 
    constantLengthForward = 18, constantLengthReverse = 20, 
    variableLengthForward = 96, variableLengthReverse = 96,
    adapterForward = "GGAAGAGCACACGTC", 
    adapterReverse = "GGAAGAGCGTCGTGT",
    primerForward = "",
    primerReverse = "",
    wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
    wildTypeReverse = "", 
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
    mutNameDelimiter = ".",
    maxNReads = -1, verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)

  L <- Ldef; L$fastqForward <- c(fqt1, fqt1); L$fastqReverse <- NULL
  res2 <- do.call(digestFastqs, L)
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 297L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadTooShort, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 0L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrMutQualTooLow, 418L)
  expect_equal(res$filterSummary$f9_nbrTooManyMutCodons, 12L)
  expect_equal(res$filterSummary$f10_nbrForbiddenCodons, 1L)
  expect_equal(res$filterSummary$nbrRetained, 272L)
  
  expect_equal(res$filterSummary * 2, res2$filterSummary)
  expect_equal(res$summaryTable$nbrReads * 2, res2$summaryTable$nbrReads)
  expect_equal(res$summaryTable$nbrUmis, res2$summaryTable$nbrUmis)
  
  for (nm in setdiff(names(Ldef), c("fastqReverse", "forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
  }
  expect_equivalent(res$parameters[["fastqReverse"]], "")
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(res$summaryTable$nbrReads, res$summaryTable$maxNbrReads)
  expect_equal(sum(res$summaryTable$nbrReads == 2), 29L)

  ## Check that mutant naming worked (compare to manual matching)
  example_seq <- paste0("ACTGATACAACCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTG", 
                        "CAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA")
  expect_equal(res$summaryTable$mutantName[res$summaryTable$sequence == example_seq], 
               "f.4.ACC")
  
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 14], 187L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 22], 47L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 27], 314L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 33], 469L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 37], 3868L)
  
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 14], 6L)
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 27], 4L)
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 37], 1L)
  
})

context("digestFastqs - trans - collapsing")
test_that("digestFastqs works as expected for trans experiments, when similar sequences are collapsed", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2, 
    mergeForwardReverse = FALSE, 
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
    revComplForward = FALSE, revComplReverse = FALSE,
    skipForward = 1, skipReverse = 1, 
    umiLengthForward = 10, umiLengthReverse = 8, 
    constantLengthForward = 18, constantLengthReverse = 20, 
    variableLengthForward = 96, variableLengthReverse = 96,
    adapterForward = "GGAAGAGCACACGTC", 
    adapterReverse = "GGAAGAGCGTCGTGT",
    primerForward = "",
    primerReverse = "",
    wildTypeForward = "",
    wildTypeReverse = "", 
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
    variableCollapseMaxDist = 6, umiCollapseMaxDist = 4, variableCollapseMinReads = 0,
    maxNReads = -1, verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  L <- Ldef; L$variableCollapseMaxDist = 0; L$umiCollapseMaxDist = 0; res0 <- do.call(digestFastqs, L)
  expect_equal(res$summaryTable$maxNbrReads, res0$summaryTable$nbrReads[match(res$summaryTable$mutantName, 
                                                                              res0$summaryTable$mutantName)])
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadTooShort, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f9_nbrTooManyMutCodons, 0L)
  expect_equal(res$filterSummary$f10_nbrForbiddenCodons, 0L)
  expect_equal(res$filterSummary$nbrRetained, 679L)
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(nrow(res$summaryTable), 294L)
  expect_equal(sum(res$summaryTable$nbrUmis), 677)
})

test_that("digestFastqs works as expected for trans experiments, when similar sequences are collapsed - extreme case", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2, 
    mergeForwardReverse = FALSE, 
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
    revComplForward = FALSE, revComplReverse = FALSE,
    skipForward = 1, skipReverse = 1, 
    umiLengthForward = 10, umiLengthReverse = 8, 
    constantLengthForward = 18, constantLengthReverse = 20, 
    variableLengthForward = 96, variableLengthReverse = 96,
    adapterForward = "GGAAGAGCACACGTC", 
    adapterReverse = "GGAAGAGCGTCGTGT",
    primerForward = "",
    primerReverse = "",
    wildTypeForward = "",
    wildTypeReverse = "", 
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
    variableCollapseMaxDist = 500, umiCollapseMaxDist = 100, variableCollapseMinReads = 0,
    maxNReads = -1, verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadTooShort, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f9_nbrTooManyMutCodons, 0L)
  expect_equal(res$filterSummary$f10_nbrForbiddenCodons, 0L)
  expect_equal(res$filterSummary$nbrRetained, 679L)
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(nrow(res$summaryTable), 1L)
  expect_equal(sum(res$summaryTable$nbrUmis), 1)
})

test_that("digestFastqs works as expected for trans experiments, when similar sequences are collapsed - only forward read", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = NULL, 
    mergeForwardReverse = FALSE, 
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
    revComplForward = FALSE, revComplReverse = FALSE,
    skipForward = 1, skipReverse = 1, 
    umiLengthForward = 10, umiLengthReverse = 8, 
    constantLengthForward = 18, constantLengthReverse = 20, 
    variableLengthForward = 96, variableLengthReverse = 96,
    adapterForward = "GGAAGAGCACACGTC", 
    adapterReverse = "GGAAGAGCGTCGTGT",
    primerForward = "",
    primerReverse = "",
    wildTypeForward = "",
    wildTypeReverse = "", 
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
    variableCollapseMaxDist = 10, umiCollapseMaxDist = 5, variableCollapseMinReads = 0,
    maxNReads = -1, verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 297L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadTooShort, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 0L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f9_nbrTooManyMutCodons, 0L)
  expect_equal(res$filterSummary$f10_nbrForbiddenCodons, 0L)
  expect_equal(res$filterSummary$nbrRetained, 703L)
  
  for (nm in setdiff(names(Ldef), c("fastqReverse", "forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(nrow(res$summaryTable), 52L)
  expect_equal(sum(res$summaryTable$nbrUmis), 97)
})

test_that("digestFastqs works as expected for trans experiments, when similar sequences are collapsed (only abundant ones)", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2, 
    mergeForwardReverse = FALSE, 
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
    revComplForward = FALSE, revComplReverse = FALSE,
    skipForward = 1, skipReverse = 1, 
    umiLengthForward = 10, umiLengthReverse = 8, 
    constantLengthForward = 18, constantLengthReverse = 20, 
    variableLengthForward = 96, variableLengthReverse = 96,
    adapterForward = "GGAAGAGCACACGTC", 
    adapterReverse = "GGAAGAGCGTCGTGT",
    primerForward = "",
    primerReverse = "",
    wildTypeForward = "",
    wildTypeReverse = "", 
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
    variableCollapseMaxDist = 2, umiCollapseMaxDist = 0,
    variableCollapseMinReads = 1.5,
    maxNReads = -1, verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  L <- Ldef; L$variableCollapseMaxDist = 0; L$umiCollapseMaxDist = 0; res0 <- do.call(digestFastqs, L)
  expect_equal(res$summaryTable$maxNbrReads, res0$summaryTable$nbrReads[match(res$summaryTable$mutantName, 
                                                                              res0$summaryTable$mutantName)])
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadTooShort, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f9_nbrTooManyMutCodons, 0L)
  expect_equal(res$filterSummary$f10_nbrForbiddenCodons, 0L)
  expect_equal(res$filterSummary$nbrRetained, 679L)
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    if (nm == "variableCollapseMinReads") {
      expect_equivalent(res$parameters[[nm]], as.integer(ceiling(Ldef[[nm]])))
    } else {
      expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
    }
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(nrow(res$summaryTable), 673L)
  expect_equal(sum(res$summaryTable$nbrUmis), 679)
})

context("digestFastqs - cis")
test_that("digestFastqs works as expected for cis experiments", {
  fqc1 <- system.file("extdata/cisInput_1.fastq.gz", package = "mutscan")
  fqc2 <- system.file("extdata/cisInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqc1, fastqReverse = fqc2, 
    mergeForwardReverse = TRUE, 
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 1, greedyOverlap = TRUE, 
    revComplForward = FALSE, revComplReverse = TRUE,
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
    maxNReads = -1, verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 126L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadTooShort, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 0L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 44L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f9_nbrTooManyMutCodons, 581L)
  expect_equal(res$filterSummary$f10_nbrForbiddenCodons, 82L)
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

context("digestFastqs - primers")
test_that("digestFastqs works as expected for primer experiments", {
  fqt1 <- system.file("extdata/leujunt0_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/leujunt0_2.fastq.gz", package = "mutscan")
  leu <- c(ATF2 = "GATCCTGATGAAAAAAGGAGAAAGTTTTTAGAGCGAAATAGAGCAGCAGCTTCAAGATGCCGACAAAAAAGGAAAGTCTGGGTTCAGTCTTTAGAGAAGAAAGCTGAAGACTTGAGTTCATTAAATGGTCAGCTGCAGAGTGAAGTCACCCTGCTGAGAAATGAAGTGGCACAGCTGAAACAGCTTCTTCTGGCT",
           ATF7 = "GATCCAGATGAGCGACGGCAGCGCTTTCTGGAGCGCAACCGGGCTGCAGCCTCCCGCTGCCGCCAAAAGCGAAAGCTGTGGGTGTCCTCCCTAGAGAAGAAGGCCGAAGAACTCACTTCTCAGAACATTCAGCTGAGTAATGAAGTCACATTACTACGCAATGAGGTGGCCCAGTTGAAACAGCTACTGTTAGCT",
           CREB5 = "GATCCGGACGAGAGGCGGCGGAAATTTCTGGAACGGAACCGGGCAGCTGCCACCCGCTGCAGACAGAAGAGGAAGGTCTGGGTGATGTCATTGGAAAAGAAAGCAGAAGAACTCACCCAGACAAACATGCAGCTTCAGAATGAAGTGTCTATGTTGAAAAATGAGGTGGCCCAGCTGAAACAGTTGTTGTTAACA",
           ATF3 = "GAAGAAGATGAAAGGAAAAAGAGGCGACGAGAAAGAAATAAGATTGCAGCTGCAAAGTGCCGAAACAAGAAGAAGGAGAAGACGGAGTGCCTGCAGAAAGAGTCGGAGAAGCTGGAAAGTGTGAATGCTGAACTGAAGGCTCAGATTGAGGAGCTCAAGAACGAGAAGCAGCATTTGATATACATGCTCAACCTT",
           JDP2 = "GAGGAAGAGGAGCGAAGGAAAAGGCGCCGGGAGAAGAACAAAGTCGCAGCAGCCCGATGCCGGAACAAGAAGAAGGAGCGCACGGAGTTTCTGCAGCGGGAATCCGAGCGGCTGGAACTCATGAACGCAGAGCTGAAGACCCAGATTGAGGAGCTGAAGCAGGAGCGGCAGCAGCTCATCCTGATGCTGAACCGA",
           ATF4 = "GAGAAACTGGATAAGAAGCTGAAAAAAATGGAGCAAAACAAGACAGCAGCCACTAGGTACCGCCAGAAGAAGAGGGCGGAGCAGGAGGCTCTTACTGGTGAGTGCAAAGAGCTGGAAAAGAAGAACGAGGCTCTAAAAGAGAGGGCGGATTCCCTGGCCAAGGAGATCCAGTACCTGAAAGATTTGATAGAAGAG",
           ATF5 = "ACCCGAGGGGACCGCAAGCAAAAGAAGAGAGACCAGAACAAGTCGGCGGCTCTGAGGTACCGCCAGCGGAAGCGGGCAGAGGGTGAGGCCCTGGAGGGCGAGTGCCAGGGGCTGGAGGCACGGAATCGCGAGCTGAAGGAACGGGCAGAGTCCGTGGAGCGCGAGATCCAGTACGTCAAGGACCTGCTCATCGAG",
           CREBZF = "AGTCCCCGGAAGGCGGCGGCGGCCGCTGCCCGCCTTAATCGACTGAAGAAGAAGGAGTACGTGATGGGGCTGGAGAGTCGAGTCCGGGGTCTGGCAGCCGAGAACCAGGAGCTGCGGGCCGAGAATCGGGAGCTGGGCAAACGCGTACAGGCACTGCAGGAGGAGAGTCGCTACCTACGGGCAGTCTTAGCCAAC",
           BATF2 = "CCCAAGGAGCAACAAAGGCAGCTGAAGAAGCAGAAGAACCGGGCAGCCGCCCAGCGAAGCCGGCAGAAGCACACAGACAAGGCAGACGCCCTGCACCAGCAGCACGAGTCTCTGGAAAAAGACAACCTCGCCCTGCGGAAGGAGATCCAGTCCCTGCAGGCCGAGCTGGCGTGGTGGAGCCGGACCCTGCACGTG",
           BATF3 = "GAGGATGATGACAGGAAGGTCCGAAGGAGAGAAAAAAACCGAGTTGCTGCTCAGAGAAGTCGGAAGAAGCAGACCCAGAAGGCTGACAAGCTCCATGAGGAATATGAGAGCCTGGAGCAAGAAAACACCATGCTGCGGAGAGAGATCGGGAAGCTGACAGAGGAGCTGAAGCACCTGACAGAGGCACTGAAGGAG",
           CEBPE = "AAAGATAGCCTTGAGTACCGGCTGAGGCGGGAGCGCAACAACATCGCCGTGCGCAAGAGCCGAGACAAGGCCAAGAGGCGCATTCTGGAGACGCAGCAGAAGGTGCTGGAGTACATGGCAGAGAACGAGCGCCTCCGCAGCCGCGTGGAGCAGCTCACCCAGGAGCTAGACACCCTCCGCAACCTCTTCCGCCAG",
           BACH1 = "CTGGATTGTATCCATGATATTCGAAGAAGAAGTAAAAACAGAATTGCTGCACAGCGCTGTCGCAAGAGAAAACTTGACTGTATACAGAATCTTGAATCAGAAATTGAGAAGCTGCAAAGTGAAAAGGAGAGCTTGTTGAAGGAAAGAGATCACATTTTGTCAACTCTGGGTGAGACAAAGCAGAACCTAACTGGA",
           BACH2 = "TTAGAGTTTATTCATGATGTCCGACGGCGCAGCAAGAACCGCATCGCGGCCCAGCGCTGCCGCAAAAGGAAACTGGACTGTATTCAGAATTTAGAATGTGAAATCCGCAAATTGGTGTGTGAGAAAGAGAAACTGTTGTCAGAGAGGAATCAACTGAAAGCATGCATGGGGGAACTGTTGGACAACTTCTCCTGC",
           NFE2L1 = "CTGAGCCTCATCCGAGACATCCGGCGCCGGGGCAAGAACAAGATGGCGGCGCAGAACTGCCGCAAGCGCAAGCTGGACACCATCCTGAATCTGGAGCGTGATGTGGAGGACCTGCAGCGTGACAAAGCCCGGCTGCTGCGGGAGAAAGTGGAGTTCCTGCGCTCCCTGCGACAGATGAAGCAGAAGGTCCAGAGC",
           NFE2 = "CTAGCGCTAGTCCGGGACATCCGACGACGGGGCAAAAACAAGGTGGCAGCCCAGAACTGCCGCAAGAGGAAGCTGGAAACCATTGTGCAGCTGGAGCGGGAGCTGGAGCGGCTGACCAATGAACGGGAGCGGCTTCTCAGGGCCCGCGGGGAGGCAGACCGGACCCTGGAGGTCATGCGCCAACAGCTGACAGAG",
           NFIL3 = "AAGAAAGATGCTATGTATTGGGAAAAAAGGCGGAAAAATAATGAAGCTGCCAAAAGATCTCGTGAGAAGCGTCGACTGAATGACCTGGTTTTAGAGAACAAACTAATTGCACTGGGAGAAGAAAACGCCACTTTAAAAGCTGAGCTGCTTTCACTAAAATTAAAGTTTGGTTTAATTAGCTCCACAGCATATGCT",
           FOS = "GAAGAAGAAGAGAAAAGGAGAATCCGAAGGGAAAGGAATAAGATGGCTGCAGCCAAATGCCGCAACCGGAGGAGGGAGCTGACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTAGAGTTCATCCTGGCAGCT",
           FOSB = "GAGGAAGAGGAGAAGCGAAGGGTGCGCCGGGAACGAAATAAACTAGCAGCAGCTAAATGCAGGAACCGGCGGAGGGAGCTGACCGACCGACTCCAGGCGGAGACAGATCAGTTGGAGGAAGAAAAAGCAGAGCTGGAGTCGGAGATCGCCGAGCTCCAAAAGGAGAAGGAACGTCTGGAGTTTGTGCTGGTGGCC",
           FOSL1 = "GAGGAAGAGGAGCGCCGCCGAGTAAGGCGCGAGCGGAACAAGCTGGCTGCGGCCAAGTGCAGGAACCGGAGGAAGGAACTGACCGACTTCCTGCAGGCGGAGACTGACAAACTGGAAGATGAGAAATCTGGGCTGCAGCGAGAGATTGAGGAGCTGCAGAAGCAGAAGGAGCGCCTAGAGCTGGTGCTGGAAGCC",
           FOSL2 = "GAAGAGGAGGAGAAGCGTCGCATCCGGCGGGAGAGGAACAAGCTGGCTGCAGCCAAGTGCCGGAACCGACGCCGGGAGCTGACAGAGAAGCTGCAGGCGGAGACAGAGGAGCTGGAGGAGGAGAAGTCAGGCCTGCAGAAGGAGATTGCTGAGCTGCAGAAGGAGAAGGAGAAGCTGGAGTTCATGTTGGTGGCT",
           MAFB = "GTGATCCGCCTGAAGCAGAAGCGGCGGACCCTGAAGAACCGGGGCTACGCCCAGTCTTGCAGGTATAAACGCGTCCAGCAGAAGCACCACCTGGAGAATGAGAAGACGCAGCTCATTCAGCAGGTGGAGCAGCTTAAGCAGGAGGTGTCCCGGCTGGCCCGCGAGAGAGACGCCTACAAGGTCAAGTGCGAGAAA",
           JUN = "CAGGAGCGGATCAAGGCGGAGAGGAAGCGCATGAGGAACCGCATCGCTGCCTCCAAGTGCCGAAAAAGGAAGCTGGAGAGAATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTTAAACAGAAAGTCATGAAC",
           JUNB = "CAAGAGCGCATCAAAGTGGAGCGCAAGCGGCTGCGGAACCGGCTGGCGGCCACCAAGTGCCGGAAGCGGAAGCTGGAGCGCATCGCGCGCCTGGAGGACAAGGTGAAGACGCTCAAGGCCGAGAACGCGGGGCTGTCGAGTACCGCCGGCCTCCTCCGGGAGCAGGTGGCCCAGCTCAAACAGAAGGTCATGACC",
           JUND = "CAGGAGCGCATCAAGGCGGAGCGCAAGCGGCTGCGCAACCGCATCGCCGCCTCCAAGTGCCGCAAGCGCAAGCTGGAGCGCATCTCGCGCCTGGAAGAGAAAGTGAAGACCCTCAAGAGTCAGAACACGGAGCTGGCGTCCACGGCGAGCCTGCTGCGCGAGCAGGTGGCGCAGCTCAAGCAGAAAGTCCTCAGC",
           CREB3 = "GAACAAATTCTGAAACGTGTGCGGAGGAAGATTCGAAATAAAAGATCTGCTCAAGAGAGCCGCAGGAAAAAGAAGGTGTATGTTGGGGGTTTAGAGAGCAGGGTCTTGAAATACACAGCCCAGAATATGGAGCTTCAGAACAAAGTACAGCTTCTGGAGGAACAGAATTTGTCCCTTCTAGATCAACTGAGGAAA",
           HLF = "CTGAAGGATGACAAGTACTGGGCAAGGCGCAGAAAGAACAACATGGCAGCCAAGCGCTCCCGCGACGCCCGGAGGCTGAAAGAGAACCAGATCGCCATCCGGGCCTCGTTCCTGGAGAAGGAGAACTCGGCCCTCCGCCAGGAGGTGGCTGACTTGAGGAAGGAGCTGGGCAAATGCAAGAACATACTTGCCAAG",
           MAFG = "ATCGTCCAGCTGAAGCAGCGCCGGCGCACGCTCAAGAACCGCGGCTACGCTGCCAGCTGCCGCGTGAAGCGGGTGACGCAGAAGGAGGAGCTGGAGAAGCAGAAGGCGGAGCTGCAGCAGGAGGTGGAGAAGCTGGCCTCAGAGAACGCCAGCATGAAGCTGGAGCTCGACGCGCTGCGCTCCAAGTACGAGGCG",
           MAFK = "GTGACCCGCCTGAAGCAGCGTCGGCGCACACTCAAGAACCGCGGCTACGCGGCCAGCTGCCGCATCAAGCGGGTGACGCAGAAGGAGGAGCTGGAGCGGCAGCGCGTGGAGCTGCAGCAGGAGGTGGAGAAGCTGGCGCGTGAGAACAGCAGCATGCGGCTGGAGCTGGACGCCCTGCGCTCCAAGTACGAGGCG",
           XBP1 = "AGCCCCGAGGAGAAGGCGCTGAGGAGGAAACTGAAAAACAGAGTAGCAGCTCAGACTGCCAGAGATCGAAAGAAGGCTCGAATGAGTGAGCTGGAACAGCAAGTGGTAGATTTAGAAGAAGAGAACCAAAAACTTTTGCTAGAAAATCAGCTTTTACGAGAGAAAACTCATGGCCTTGTAGTTGAGAACCAGGAG",
           ATF6 = "ATTGCTGTGCTAAGGAGACAGCAACGTATGATAAAAAATCGAGAATCCGCTTGTCAGTCTCGCAAGAAGAAGAAAGAATATATGCTAGGGTTAGAGGCGAGATTAAAGGCTGCCCTCTCAGAAAACGAGCAACTGAAGAAAGAAAATGGAACACTGAAGCGGCAGCTGGATGAAGTTGTGTCAGAGAACCAGAGG",
           ATF6B = "GCAAAGCTGCTGAAGCGGCAGCAGCGAATGATCAAGAACCGGGAGTCAGCCTGCCAGTCCCGGAGAAAGAAGAAAGAGTATCTGCAGGGACTGGAGGCTCGGCTGCAAGCAGTACTGGCTGACAACCAGCAGCTCCGCCGAGAGAATGCTGCCCTCCGGCGGCGGCTGGAGGCCCTGCTGGCTGAAAACAGCGAG",
           CEBPA = "AAGAACAGCAACGAGTACCGGGTGCGGCGCGAGCGCAACAACATCGCGGTGCGCAAGAGCCGCGACAAGGCCAAGCAGCGCAACGTGGAGACGCAGCAGAAGGTGCTGGAGCTGACCAGTGACAATGACCGCCTGCGCAAGCGGGTGGAACAGCTGAGCCGCGAACTGGACACGCTGCGGGGCATCTTCCGCCAG",
           CEBPB = "AAGCACAGCGACGAGTACAAGATCCGGCGCGAGCGCAACAACATCGCCGTGCGCAAGAGCCGCGACAAGGCCAAGATGCGCAACCTGGAGACGCAGCACAAGGTCCTGGAGCTCACGGCCGAGAACGAGCGGCTGCAGAAGAAGGTGGAGCAGCTGTCGCGCGAGCTCAGCACCCTGCGGAACTTGTTCAAGCAG",
           CEBPD = "CGCGGCAGCCCCGAGTACCGGCAGCGGCGCGAGCGCAACAACATCGCCGTGCGCAAGAGCCGCGACAAGGCCAAGCGGCGCAACCAGGAGATGCAGCAGAAGTTGGTGGAGCTGTCGGCTGAGAACGAGAAGCTGCACCAGCGCGTGGAGCAGCTCACGCGGGACCTGGCCGGCCTCCGGCAGTTCTTCAAGCAG",
           CEBPG = "CGAAACAGTGACGAGTATCGGCAACGCCGAGAGAGGAACAACATGGCTGTGAAAAAGAGCCGGTTGAAAAGCAAGCAGAAAGCACAAGACACACTGCAGAGAGTCAATCAGCTCAAAGAAGAGAATGAACGGTTGGAAGCAAAAATCAAATTGCTGACCAAGGAATTAAGTGTACTCAAAGATTTGTTTCTTGAG",
           CREB1 = "GAAGCAGCACGAAAGAGAGAGGTCCGTCTAATGAAGAACAGGGAAGCAGCTCGAGAGTGTCGTAGAAAGAAGAAAGAATATGTGAAATGTTTAGAAAACAGAGTGGCAGTGCTTGAAAATCAAAACAAGACATTGATTGAGGAGCTAAAAGCACTTAAGGACCTTTACTGCCACAAATCAGAT",
           CREB3L1 = "GAGAAGGCCTTGAAGAGAGTCCGGAGGAAAATCAAGAACAAGATCTCAGCCCAGGAGAGCCGTCGTAAGAAGAAGGAGTATGTGGAGTGTCTAGAAAAGAAGGTGGAGACATTTACATCTGAGAACAATGAACTGTGGAAGAAGGTGGAGACCCTGGAGAATGCCAACAGGACCCTGCTCCAGCAGCTGCAGAAA",
           CREB3L2 = "GAGAAGGCCCTGAAGAAAATTCGGAGGAAGATCAAGAATAAGATTTCTGCTCAGGAAAGTAGGAGAAAGAAGAAAGAATACATGGACAGCCTGGAGAAAAAAGTGGAGTCTTGTTCAACTGAGAACTTGGAGCTTCGGAAGAAGGTAGAGGTTCTAGAGAACACTAATAGGACTCTCCTTCAGCAACTCCAGAAG",
           CREB3L3 = "GAGCGAGTGCTGAAAAAAATCCGCCGGAAAATCCGGAACAAGCAGTCGGCGCAAGAAAGCAGGAAGAAGAAGAAGGAATATATCGATGGCCTGGAGACTCGGATGTCAGCTTGCACTGCTCAGAATCAGGAGTTACAGAGGAAAGTCTTGCATCTCGAGAAGCAAAACCTGTCCCTCTTGGAGCAACTGAAGAAA",
           CREB3L4 = "GAGAGGGTCCTCAAGAAGGTCAGGAGGAAAATCCGTAACAAGCAGTCAGCTCAGGACAGTCGGCGGCGGAAGAAGGAGTACATTGATGGGCTGGAGAGCAGGGTGGCAGCCTGTTCTGCACAGAACCAAGAATTACAGAAAAAAGTCCAGGAGCTGGAGAGGCACAACATCTCCTTGGTAGCTCAGCTCCGCCAG",
           CREBL2 = "CCAGCCAAAATTGACTTGAAAGCAAAACTTGAGAGGAGCCGGCAGAGTGCAAGAGAATGCCGAGCCCGAAAAAAGCTGAGATATCAGTATTTGGAAGAGTTGGTATCCAGTCGAGAAAGAGCTATATGTGCCCTCAGAGAGGAACTGGAAATGTACAAGCAGTGGTGCATGGCAATGGACCAAGGAAAAATCCCT",
           CREBRF = "CCCTTAACAGCCCGACCAAGGTCAAGGAAGGAAAAAAATAAGCTGGCTTCCAGAGCTTGTCGGTTAAAGAAGAAAGCCCAGTATGAAGCTAATAAAGTGAAATTATGGGGCCTCAACACAGAATATGATAATTTATTGTTTGTAATCAACTCCATCAAGCAAGAGATTGTAAACCGGGTACAGAATCCAAGAGAT",
           DBP = "CAGAAGGATGAGAAATACTGGAGCCGGCGGTACAAGAACAACGAGGCAGCCAAGCGGTCCCGTGACGCCCGGCGGCTCAAGGAGAACCAGATATCGGTGCGGGCGGCCTTCCTGGAGAAGGAGAACGCCCTGCTGCGGCAGGAAGTTGTGGCCGTGCGCCAGGAGCTGTCCCACTACCGCGCCGTGCTGTCCCGA",
           NFE2L2 = "CTTGCATTAATTCGGGATATACGTAGGAGGGGTAAGAATAAAGTGGCTGCTCAGAATTGCAGAAAAAGAAAACTGGAAAATATAGTAGAACTAGAGCAAGATTTAGATCATTTGAAAGATGAAAAAGAAAAATTGCTCAAAGAAAAAGGAGAAAATGACAAAAGCCTTCACCTACTGAAAAAACAACTCAGCACC",
           NFE2L3 = "GTCTCACTTATCCGTGACATCAGACGAAGAGGGAAAAATAAAGTTGCTGCGCAGAACTGTCGTAAACGCAAATTGGACATAATTTTGAATTTAGAAGATGATGTATGTAACTTGCAAGCAAAGAAGGAAACTCTTAAGAGAGAGCAAGCACAATGTAACAAAGCTATTAACATAATGAAACAGAAACTGCATGAC",
           TEF = "CAGAAGGATGAAAAGTACTGGACAAGACGCAAGAAGAACAACGTGGCAGCTAAACGGTCACGGGATGCCCGGCGCCTGAAAGAGAATCAGATCACCATCCGGGCAGCCTTCCTGGAGAAGGAGAACACAGCCCTGCGGACGGAGGTGGCCGAGCTACGCAAGGAGGTGGGCAAGTGCAAGACCATCGTGTCCAAG"
  )
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2, 
    mergeForwardReverse = FALSE, 
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
    revComplForward = FALSE, revComplReverse = FALSE,
    skipForward = -1, skipReverse = -1, 
    umiLengthForward = -1, umiLengthReverse = -1, 
    constantLengthForward = -1, constantLengthReverse = -1, 
    variableLengthForward = -1, variableLengthReverse = 96,
    adapterForward = "", 
    adapterReverse = "",
    primerForward = "GTCAGGTGGAGGCGGATCC",
    primerReverse = "GAAAAAGGAAGCTGGAGAGA",
    wildTypeForward = leu,
    wildTypeReverse = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT", 
    constantForward = "", 
    constantReverse = "", 
    avePhredMinForward = 20.0, avePhredMinReverse = 20.0,
    variableNMaxForward = 0, variableNMaxReverse = 0, 
    umiNMax = 0,
    nbrMutatedCodonsMaxForward = 0,
    nbrMutatedCodonsMaxReverse = 1,
    forbiddenMutatedCodonsForward = "",
    forbiddenMutatedCodonsReverse = "NNW",
    mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".", maxNReads = -1, verbose = FALSE
  )

  res <- do.call(digestFastqs, Ldef)
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 0L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 126L)
  expect_equal(res$filterSummary$f3_nbrReadTooShort, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 0L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 76L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f9_nbrTooManyMutCodons, 58L + 137L)
  expect_equal(res$filterSummary$f10_nbrForbiddenCodons, 3L)
  expect_equal(res$filterSummary$nbrRetained, 600L)
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(sum(res$summaryTable$nbrReads == 2), 11L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 2]),
               sort(c("BACH2.0.WT_r.6.CCC", "BATF2.0.WT_r.21.TGG", "BATF2.0.WT_r.7.GGG", "CEBPB.0.WT_r.12.CTC",
                      "CEBPD.0.WT_r.13.CCC", "CEBPD.0.WT_r.23.GAG", "CREB3L3.0.WT_r.22.GGG",
                      "FOSL1.0.WT_r.22.CCC", "FOSL2.0.WT_r.13.GGG", "JUNB.0.WT_r.0.WT", "XBP1.0.WT_r.0.WT")))
  expect_equal(res$summaryTable$nbrUmis, rep(1L, nrow(res$summaryTable))) ## no UMIs in this experiment
})
