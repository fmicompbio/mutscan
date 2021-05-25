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
    elementsForward = "SUCV", elementsReverse = "SUCVS",
    elementLengthsForward = c(1, 10, 18, 96),
    elementLengthsReverse = c(1, 7, 17, 96, -1),
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
    nbrMutatedBasesMaxForward = -1,
    nbrMutatedBasesMaxReverse = -1,
    forbiddenMutatedCodonsForward = "NNW",
    forbiddenMutatedCodonsReverse = "NNW",
    useTreeWTmatch = FALSE,
    mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".",
    constantMaxDistForward = -1,
    constantMaxDistReverse = -1,
    variableCollapseMaxDist = 0,
    variableCollapseMinReads = 0,
    variableCollapseMinRatio = 0,
    umiCollapseMaxDist = 0,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = FALSE,
    nThreads = 1, chunkSize = 1000
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 126L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 0L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 44L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 581L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 82L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
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
  
  ## Check the number of reads with a given number of mismatches
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 0]), 77L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 1]), 88L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 2]), 2L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutCodons == 0]), 77L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutCodons == 1]), 90L)
  
  ## Similarly but for the number of unique sequences
  expect_equal(sum(res$summaryTable$nbrMutBases == 0), 1L)
  expect_equal(sum(res$summaryTable$nbrMutBases == 1), 64L)
  expect_equal(sum(res$summaryTable$nbrMutBases == 2), 2L)
  expect_equal(sum(res$summaryTable$nbrMutCodons == 0), 1L)
  expect_equal(sum(res$summaryTable$nbrMutCodons == 1), 66L)
  

  ## Check that mutant naming worked (compare to manual matching)
  example_seq <- paste0("ACTGATACACGCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCT", 
                        "GCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA")
  expect_equal(res$summaryTable$mutantName[res$summaryTable$sequence == example_seq], 
               "f.4.CGC")
  
  expect_equal(sum(res$errorStatistics$nbrMatchForward + res$errorStatistics$nbrMismatchForward),
               nchar(Ldef$constantForward[1]) * res$filterSummary$nbrRetained)
  expect_equal(sum(res$errorStatistics$nbrMatchReverse + res$errorStatistics$nbrMismatchReverse),
               nchar(Ldef$constantReverse[1]) * res$filterSummary$nbrRetained)
  
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

test_that("digestFastqs works as expected when specifying max nbr of mutated bases rather than codons", {
  fqc1 <- system.file("extdata/cisInput_1.fastq.gz", package = "mutscan")
  fqc2 <- system.file("extdata/cisInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqc1, fastqReverse = fqc2, 
    mergeForwardReverse = TRUE, 
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 1, greedyOverlap = TRUE, 
    revComplForward = FALSE, revComplReverse = TRUE,
    elementsForward = "SUCV", elementsReverse = "SUCVS",
    elementLengthsForward = c(1, 10, 18, 96),
    elementLengthsReverse = c(1, 7, 17, 96, -1),
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
    nbrMutatedCodonsMaxForward = -1,
    nbrMutatedCodonsMaxReverse = -1,
    nbrMutatedBasesMaxForward = 2,
    nbrMutatedBasesMaxReverse = 2,
    forbiddenMutatedCodonsForward = "NNW",
    forbiddenMutatedCodonsReverse = "NNW",
    useTreeWTmatch = FALSE,
    mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".",
    constantMaxDistForward = -1,
    constantMaxDistReverse = -1,
    variableCollapseMaxDist = 0,
    variableCollapseMinReads = 0,
    variableCollapseMinRatio = 0,
    umiCollapseMaxDist = 0,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = FALSE,
    nThreads = 1, chunkSize = 1000
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 126L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 0L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 44L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 0L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 416L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 212L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
  expect_equal(res$filterSummary$nbrRetained, 202L)
  
  expect_named(res$summaryTable, c("mutantName", "sequence", "nbrMutBases", "nbrMutCodons",
                                   "nbrReads", "maxNbrReads", "nbrUmis"))
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(sum(res$summaryTable$nbrReads == 3), 5L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 3]),
               sort(c("f.12.G", "f.17.T", "f.60.G", "f.85.T", "f.96.G")))
  expect_true(all(res$summaryTable$nbrReads == res$summaryTable$nbrUmis))
  
  ## Check the number of reads with a given number of mismatches
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 0]), 77L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 1]), 88L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 2]), 37L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutCodons == 0]), 77L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutCodons == 1]), 90L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutCodons == 2]), 35L)

  ## Similarly but for the number of unique sequences
  expect_equal(sum(res$summaryTable$nbrMutBases == 0), 1L)
  expect_equal(sum(res$summaryTable$nbrMutBases == 1), 64L)
  expect_equal(sum(res$summaryTable$nbrMutBases == 2), 36L)
  expect_equal(sum(res$summaryTable$nbrMutCodons == 0), 1L)
  expect_equal(sum(res$summaryTable$nbrMutCodons == 1), 66L)
  expect_equal(sum(res$summaryTable$nbrMutCodons == 2), 34L)
  
  ## Check that mutant naming worked (compare to manual matching)
  example_seq <- paste0("ACTGATACACGCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCT",
                        "GCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA")
  expect_equal(res$summaryTable$mutantName[res$summaryTable$sequence == example_seq],
               "f.11.G")

  expect_equal(sum(res$errorStatistics$nbrMatchForward + res$errorStatistics$nbrMismatchForward),
               nchar(Ldef$constantForward[1]) * res$filterSummary$nbrRetained)
  expect_equal(sum(res$errorStatistics$nbrMatchReverse + res$errorStatistics$nbrMismatchReverse),
               nchar(Ldef$constantReverse[1]) * res$filterSummary$nbrRetained)
  
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 14], 61L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 22], 22L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 27], 97L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 33], 203L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 37], 3249L)

  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 14], 1L)
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 27], 1L)
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 37], 2L)

  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 14], 42L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 27], 63L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 33], 117L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 37], 3184L)

  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 2], 5L)
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 14], 7L)
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 33], 2L)
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 37], 14L)
  
})

test_that("digestFastqs gives the same results regardless of how WT matching is done", {
  fqc1 <- system.file("extdata/cisInput_1.fastq.gz", package = "mutscan")
  fqc2 <- system.file("extdata/cisInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqc1, fastqReverse = fqc2, 
    mergeForwardReverse = TRUE, 
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 1, greedyOverlap = TRUE, 
    revComplForward = FALSE, revComplReverse = TRUE,
    elementsForward = "SUCV", elementsReverse = "SUCVS",
    elementLengthsForward = c(1, 10, 18, 96),
    elementLengthsReverse = c(1, 7, 17, 96, -1),
    adapterForward = "GGAAGAGCACACGTC", 
    adapterReverse = "GGAAGAGCGTCGTGT",
    primerForward = "",
    primerReverse = "",
    wildTypeForward = c(wt1 = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
                        wt2 = "CCTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
                        wt3 = "ACTGATACACTCCAAGCGGAGACTGACCAACTAGAAGATGAGAAGTCTTCTTTGCAGACCGAGATTGCCAACGTGCTGAAGGAGAAGGAAAAACTA"),
    wildTypeReverse = "", 
    constantForward = "AACCGGAGGAGGGAGCTG", 
    constantReverse = "GAGTTCATCCTGGCAGC", 
    avePhredMinForward = 20.0, avePhredMinReverse = 20.0,
    variableNMaxForward = 0, variableNMaxReverse = 0, 
    umiNMax = 0,
    nbrMutatedCodonsMaxForward = -1,
    nbrMutatedCodonsMaxReverse = -1,
    nbrMutatedBasesMaxForward = 2,
    nbrMutatedBasesMaxReverse = 2,
    forbiddenMutatedCodonsForward = "NNW",
    forbiddenMutatedCodonsReverse = "NNW",
    useTreeWTmatch = FALSE,
    mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".",
    constantMaxDistForward = -1,
    constantMaxDistReverse = -1,
    variableCollapseMaxDist = 0,
    variableCollapseMinReads = 0,
    variableCollapseMinRatio = 0,
    umiCollapseMaxDist = 0,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = FALSE,
    nThreads = 1, chunkSize = 1000
  )
  
  res1 <- do.call(digestFastqs, Ldef)
  Ldef$useTreeWTmatch <- TRUE
  res2 <- do.call(digestFastqs, Ldef)
  
  expect_equal(res1$filterSummary$nbrTotal, res2$filterSummary$nbrTotal)
  expect_equal(res1$filterSummary$f1_nbrAdapter, res2$filterSummary$f1_nbrAdapter)
  expect_equal(res1$filterSummary$f2_nbrNoPrimer, res2$filterSummary$f2_nbrNoPrimer)
  expect_equal(res1$filterSummary$f3_nbrReadWrongLength, res2$filterSummary$f3_nbrReadWrongLength)
  expect_equal(res1$filterSummary$f4_nbrNoValidOverlap, res2$filterSummary$f4_nbrNoValidOverlap)
  expect_equal(res1$filterSummary$f5_nbrAvgVarQualTooLow, res2$filterSummary$f5_nbrAvgVarQualTooLow)
  expect_equal(res1$filterSummary$f6_nbrTooManyNinVar, res2$filterSummary$f6_nbrTooManyNinVar)
  expect_equal(res1$filterSummary$f7_nbrTooManyNinUMI, res2$filterSummary$f7_nbrTooManyNinUMI)
  expect_equal(res1$filterSummary$f8_nbrTooManyBestWTHits, res2$filterSummary$f8_nbrTooManyBestWTHits)
  expect_equal(res1$filterSummary$f9_nbrMutQualTooLow, res2$filterSummary$f9_nbrMutQualTooLow)
  expect_equal(res1$filterSummary$f10a_nbrTooManyMutCodons, res2$filterSummary$f10a_nbrTooManyMutCodons)
  expect_equal(res1$filterSummary$f10b_nbrTooManyMutBases, res2$filterSummary$f10b_nbrTooManyMutBases)
  expect_equal(res1$filterSummary$f11_nbrForbiddenCodons, res2$filterSummary$f11_nbrForbiddenCodons)
  expect_equal(res1$filterSummary$f12_nbrTooManyMutConstant, res2$filterSummary$f12_nbrTooManyMutConstant)
  expect_equal(res1$filterSummary$f13_nbrTooManyBestConstantHits, res2$filterSummary$f13_nbrTooManyBestConstantHits)
  expect_equal(res1$filterSummary$nbrRetained, res2$filterSummary$nbrRetained)
  
  expect_equal(res1$summaryTable$nbrReads, res2$summaryTable$nbrReads)
})

