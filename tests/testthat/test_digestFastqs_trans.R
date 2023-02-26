test_that("digestFastqs works as expected for trans experiments", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2,
    mergeForwardReverse = FALSE,
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE,
    revComplForward = FALSE, revComplReverse = FALSE,
    elementsForward = "SUCVV", elementsReverse = "SUCVV",
    elementLengthsForward = c(1, 10, 18, 80, 16),
    elementLengthsReverse = c(1, 8, 20, 16, 80),
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

  L <- Ldef; L$fastqForward <- c(fqt1, fqt1); L$fastqReverse <- c(fqt2, fqt2)
  res2 <- do.call(digestFastqs, L)

  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 192L + 95L + 68L + 37L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 6L + 2L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
  expect_equal(res$filterSummary$nbrRetained, 279L)

  expect_equal(res$filterSummary * 2, res2$filterSummary)
  expect_equal(res$summaryTable$nbrReads * 2, res2$summaryTable$nbrReads)
  expect_equal(res$summaryTable$nbrUmis, res2$summaryTable$nbrUmis)

  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose", "fastqForward", "fastqReverse"))) {
    expect_equal(res$parameters[[nm]], Ldef[[nm]], ignore_attr = TRUE)
  }
  for (nm in c("fastqForward", "fastqReverse")) {
    expect_equal(res$parameters[[nm]], normalizePath(Ldef[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(res$summaryTable$nbrReads, res$summaryTable$maxNbrReads)
  expect_equal(sum(res$summaryTable$nbrReads == 2), 2L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 2]),
               sort(c("f.13.GAG_r.0.WT", "f.26.TAG_r.8.AGG")))
  expect_equal(sort(res$summaryTable$mutantNameCodon[res$summaryTable$nbrReads == 2]),
               sort(c("f.13.GAG_r.0.WT", "f.26.TAG_r.8.AGG")))
  expect_equal(sort(res$summaryTable$mutantNameBase[res$summaryTable$nbrReads == 2]),
               sort(c("f.39.G_r.0.WT", "f.76.T_f.77.A_r.22.A_r.23.G")))
  expect_equal(sort(res$summaryTable$mutantNameBaseHGVS[res$summaryTable$nbrReads == 2]),
               sort(c("f:c.39T>G_r:c", "f:c.76_77delinsTA_r:c.22_23delinsAG")))
  expect_equal(sort(res$summaryTable$mutantNameAAHGVS[res$summaryTable$nbrReads == 2]),
               sort(c("f:p.(Asp13Glu)_r:p", "f:p.(Leu26*)_r:p.(Val8Arg)")))
  expect_true(all(res$summaryTable$nbrReads == res$summaryTable$nbrUmis))

  ## Check the number of reads with a given number of mismatches
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 1]), 4L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 2]), 14L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 3]), 52L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 4]), 83L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 5]), 87L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 6]), 39L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutCodons == 1]), 15L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutCodons == 2]), 264L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutAAs == 0]), 1L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutAAs == 1]), 37L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutAAs == 2]), 241L)

  ## Similarly but for the number of unique sequences
  expect_equal(sum(res$summaryTable$nbrMutBases == 1), 3L)
  expect_equal(sum(res$summaryTable$nbrMutBases == 2), 14L)
  expect_equal(sum(res$summaryTable$nbrMutBases == 3), 52L)
  expect_equal(sum(res$summaryTable$nbrMutBases == 4), 82L)
  expect_equal(sum(res$summaryTable$nbrMutBases == 5), 87L)
  expect_equal(sum(res$summaryTable$nbrMutBases == 6), 39L)
  expect_equal(sum(res$summaryTable$nbrMutCodons == 1), 14L)
  expect_equal(sum(res$summaryTable$nbrMutCodons == 2), 263L)
  expect_equal(sum(res$summaryTable$nbrMutAAs == 0), 1L)
  expect_equal(sum(res$summaryTable$nbrMutAAs == 1), 36L)
  expect_equal(sum(res$summaryTable$nbrMutAAs == 2), 240L)

  ## check that variable segment lengths are reported correctly
  expect_true(all(res$summaryTable$varLengths == "80,16_16,80"))

  ## Check that mutant naming worked (compare to manual matching)
  example_seq <- paste0("ACTGATACAACCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTG",
                        "CAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA_ATCGCCCGGCT",
                        "GGAGGAAAAAGTGGGCACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGC",
                        "CAACATGCTCAGGGAACAGGTGGCACAGCTT")
  expect_equal(res$summaryTable$mutantName[res$summaryTable$sequence == example_seq],
               "f.4.ACC_r.9.GGC")
  expect_equal(res$summaryTable$mutantNameAA[res$summaryTable$sequence == example_seq],
               "f.4.T_r.9.G")
  expect_equal(res$summaryTable$mutantNameCodon[res$summaryTable$sequence == example_seq],
               "f.4.ACC_r.9.GGC")
  expect_equal(res$summaryTable$mutantNameBase[res$summaryTable$sequence == example_seq],
               "f.10.A_f.11.C_r.25.G_r.26.G_r.27.C")
  expect_equal(res$summaryTable$mutantNameBaseHGVS[res$summaryTable$sequence == example_seq],
               "f:c.10_11delinsAC_r:c.25_27delinsGGC")
  expect_equal(res$summaryTable$mutantNameAAHGVS[res$summaryTable$sequence == example_seq],
               "f:p.(Leu4Thr)_r:p.(Lys9Gly)")
  
  expect_equal(sum(grepl("WT", res$summaryTable$mutantNameAA) & res$summaryTable$nbrMutBases > 0), 37L)
  expect_true(all(grepl("silent", res$summaryTable$mutationTypes[grepl("WT", res$summaryTable$mutantNameAA) & grepl("f\\.[1-9]", res$summaryTable$mutantName) & grepl("r\\.[1-9]", res$summaryTable$mutantName)])))
  expect_true(all(grepl("nonsynonymous|stop", res$summaryTable$mutationTypes[grepl("\\.[1-9][0-9]*\\.", res$summaryTable$mutantNameAA)])))
  expect_true(all(grepl("stop", res$summaryTable$mutationTypes[grepl("TAG|TAA|TGA", res$summaryTable$mutantName)])))

  expect_equal(sum(res$errorStatistics$nbrMatchForward + res$errorStatistics$nbrMismatchForward),
               nchar(Ldef$constantForward[1]) * res$filterSummary$nbrRetained)
  expect_equal(sum(res$errorStatistics$nbrMatchReverse + res$errorStatistics$nbrMismatchReverse),
               nchar(Ldef$constantReverse[1]) * res$filterSummary$nbrRetained)

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

## ----------------------------------------------------------------------------
## Constant seq mismatch filtering
## ----------------------------------------------------------------------------
test_that("digestFastqs works as expected for trans experiments, filter based on constant seq mismatches", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2,
    mergeForwardReverse = FALSE,
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE,
    revComplForward = FALSE, revComplReverse = FALSE,
    elementsForward = "SUCV", elementsReverse = "SUCV",
    elementLengthsForward = c(1, 10, 18, 96),
    elementLengthsReverse = c(1, 8, 20, 96),
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
    nbrMutatedBasesMaxForward = -1,
    nbrMutatedBasesMaxReverse = -1,
    forbiddenMutatedCodonsForward = "NNW",
    forbiddenMutatedCodonsReverse = "NNW",
    useTreeWTmatch = FALSE,
    mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".",
    constantMaxDistForward = 0,
    constantMaxDistReverse = 0,
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
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 287L + 105L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 6L + 2L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 28L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
  expect_equal(res$filterSummary$nbrRetained, 251L)

  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose", "fastqForward", "fastqReverse"))) {
    expect_equal(res$parameters[[nm]], Ldef[[nm]], ignore_attr = TRUE)
  }
  for (nm in c("fastqForward", "fastqReverse")) {
    expect_equal(res$parameters[[nm]], normalizePath(Ldef[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }
  

  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(res$summaryTable$nbrReads, res$summaryTable$maxNbrReads)
  expect_true(all(res$summaryTable$nbrReads == res$summaryTable$nbrUmis))
  expect_equal(sum(res$summaryTable$nbrReads == 2), 1L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 2]),
               sort(c("f.26.TAG_r.8.AGG")))

  ## Check the number of reads with a given number of mismatches
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 1]), 3L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 2]), 14L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 3]), 46L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 4]), 76L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 5]), 75L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 6]), 37L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutCodons == 1]), 14L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutCodons == 2]), 237L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutAAs == 0]), 1L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutAAs == 1]), 35L)
  expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutAAs == 2]), 215L)
  
  ## Similarly but for the number of unique sequences
  expect_equal(sum(res$summaryTable$nbrMutBases == 1), 3L)
  expect_equal(sum(res$summaryTable$nbrMutBases == 2), 14L)
  expect_equal(sum(res$summaryTable$nbrMutBases == 3), 46L)
  expect_equal(sum(res$summaryTable$nbrMutBases == 4), 75L)
  expect_equal(sum(res$summaryTable$nbrMutBases == 5), 75L)
  expect_equal(sum(res$summaryTable$nbrMutBases == 6), 37L)
  expect_equal(sum(res$summaryTable$nbrMutCodons == 1), 14L)
  expect_equal(sum(res$summaryTable$nbrMutCodons == 2), 236L)
  expect_equal(sum(res$summaryTable$nbrMutAAs == 0), 1L)
  expect_equal(sum(res$summaryTable$nbrMutAAs == 1), 35L)
  expect_equal(sum(res$summaryTable$nbrMutAAs == 2), 214L)
  
  
  expect_equal(sum(res$errorStatistics$nbrMatchForward + res$errorStatistics$nbrMismatchForward),
               nchar(Ldef$constantForward[1]) * res$filterSummary$nbrRetained)
  expect_equal(sum(res$errorStatistics$nbrMatchReverse + res$errorStatistics$nbrMismatchReverse),
               nchar(Ldef$constantReverse[1]) * res$filterSummary$nbrRetained)

  expect_equal(res$errorStatistics$nbrMismatchForward, rep(0L, nrow(res$errorStatistics)))
  expect_equal(res$errorStatistics$nbrMismatchReverse, rep(0L, nrow(res$errorStatistics)))

  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 14], 126L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 22], 43L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 27], 248L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 33], 398L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 37], 3703L)

  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 14], 137L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 22], 11L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 27], 380L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 33], 407L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 37], 4085L)

  ## --------------------------------------------------------------------------
  ## Allow 1 mismatch
  ## --------------------------------------------------------------------------
  L <- Ldef; L$constantMaxDistForward <- 1L; L$constantMaxDistReverse <- 1L
  res <- do.call(digestFastqs, L)

  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 287L + 105L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 6L + 2L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 4L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
  expect_equal(res$filterSummary$nbrRetained, 275L)

  for (nm in setdiff(names(L), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose", "fastqForward", "fastqReverse"))) {
    expect_equal(res$parameters[[nm]], L[[nm]], ignore_attr = TRUE)
  }
  for (nm in c("fastqForward", "fastqReverse")) {
    expect_equal(res$parameters[[nm]], normalizePath(L[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(res$summaryTable$nbrReads, res$summaryTable$maxNbrReads)
  expect_true(all(res$summaryTable$nbrReads == res$summaryTable$nbrUmis))
  expect_equal(sum(res$summaryTable$nbrReads == 2), 2L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 2]),
               sort(c("f.13.GAG_r.0.WT", "f.26.TAG_r.8.AGG")))

  expect_equal(sum(res$errorStatistics$nbrMatchForward + res$errorStatistics$nbrMismatchForward),
               nchar(L$constantForward[1]) * res$filterSummary$nbrRetained)
  expect_equal(sum(res$errorStatistics$nbrMatchReverse + res$errorStatistics$nbrMismatchReverse),
               nchar(L$constantReverse[1]) * res$filterSummary$nbrRetained)

  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 14], 147L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 22], 49L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 27], 288L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 33], 459L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 37], 3994L)

  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 14], 9L)
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 27], 3L)
  expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 37], 1L)

  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 14], 178L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 22], 12L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 27], 450L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 33], 454L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 37], 4394L)

  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 14], 10L)
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 27], 1L)
  expect_equal(res$errorStatistics$nbrMismatchReverse[res$errorStatistics$PhredQuality == 33], 1L)

  ## --------------------------------------------------------------------------
  ## Provide 2 constant sequences
  ## --------------------------------------------------------------------------
  L <- Ldef; L$constantForward <- c("AACCGGAGGAGGGAGCTG", "AACCGGCGGAGGGAGCTG")
  res <- do.call(digestFastqs, L)

  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 287L + 105L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 6L + 2L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 26L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
  expect_equal(res$filterSummary$nbrRetained, 253L)

  for (nm in setdiff(names(L), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose", "fastqForward", "fastqReverse"))) {
    expect_equal(res$parameters[[nm]], L[[nm]], ignore_attr = TRUE)
  }
  for (nm in c("fastqForward", "fastqReverse")) {
    expect_equal(res$parameters[[nm]], normalizePath(L[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }
  

  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(res$summaryTable$nbrReads, res$summaryTable$maxNbrReads)
  expect_true(all(res$summaryTable$nbrReads == res$summaryTable$nbrUmis))
  expect_equal(sum(res$summaryTable$nbrReads == 2), 1L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 2]),
               sort(c("f.26.TAG_r.8.AGG")))

  expect_equal(sum(res$errorStatistics$nbrMatchForward + res$errorStatistics$nbrMismatchForward),
               nchar(L$constantForward[1]) * res$filterSummary$nbrRetained)
  expect_equal(sum(res$errorStatistics$nbrMatchReverse + res$errorStatistics$nbrMismatchReverse),
               nchar(L$constantReverse[1]) * res$filterSummary$nbrRetained)

  ## --------------------------------------------------------------------------
  ## Split constant sequence in two and let the function concatenate them
  ## --------------------------------------------------------------------------
  L <- Ldef; L$elementsForward <- "SUCCV"; L$elementsReverse <- "SUCCV"
  L$elementLengthsForward <- c(1, 10, 10, 8, 96); L$elementLengthsReverse <- c(1, 8, 9, 11, 96)

  res <- do.call(digestFastqs, L)

  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 287L + 105L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 6L + 2L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 28L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
  expect_equal(res$filterSummary$nbrRetained, 251L)

  for (nm in setdiff(names(L), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose", "fastqForward", "fastqReverse"))) {
    expect_equal(res$parameters[[nm]], L[[nm]], ignore_attr = TRUE)
  }
  for (nm in c("fastqForward", "fastqReverse")) {
    expect_equal(res$parameters[[nm]], normalizePath(L[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }
  
  

  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(res$summaryTable$nbrReads, res$summaryTable$maxNbrReads)
  expect_true(all(res$summaryTable$nbrReads == res$summaryTable$nbrUmis))
  expect_equal(sum(res$summaryTable$nbrReads == 2), 1L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 2]),
               sort(c("f.26.TAG_r.8.AGG")))

  expect_equal(sum(res$errorStatistics$nbrMatchForward + res$errorStatistics$nbrMismatchForward),
               nchar(L$constantForward[1]) * res$filterSummary$nbrRetained)
  expect_equal(sum(res$errorStatistics$nbrMatchReverse + res$errorStatistics$nbrMismatchReverse),
               nchar(L$constantReverse[1]) * res$filterSummary$nbrRetained)

  expect_equal(res$errorStatistics$nbrMismatchForward, rep(0L, nrow(res$errorStatistics)))
  expect_equal(res$errorStatistics$nbrMismatchReverse, rep(0L, nrow(res$errorStatistics)))

  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 14], 126L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 22], 43L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 27], 248L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 33], 398L)
  expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 37], 3703L)

  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 14], 137L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 22], 11L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 27], 380L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 33], 407L)
  expect_equal(res$errorStatistics$nbrMatchReverse[res$errorStatistics$PhredQuality == 37], 4085L)

})

## ----------------------------------------------------------------------------
## Read composition specification
## ----------------------------------------------------------------------------
test_that("digestFastqs works as expected for trans experiments - no UMI specified", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2,
    mergeForwardReverse = FALSE,
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE,
    revComplForward = FALSE, revComplReverse = FALSE,
    elementsForward = "SSCV", elementsReverse = "SSCV",
    elementLengthsForward = c(1, 10, 18, 96),
    elementLengthsReverse = c(1, 8, 20, 96),
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
    nThreads = 1, chunkSize = 500
  )

  res <- do.call(digestFastqs, Ldef)

  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 287L + 105L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 6L + 2L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
  expect_equal(res$filterSummary$nbrRetained, 279L)

  expect_equal(res$summaryTable$nbrUmis, rep(0, nrow(res$summaryTable)))
})

test_that("digestFastqs works as expected for trans experiments, when variable sequence length is not given", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2,
    mergeForwardReverse = FALSE,
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE,
    revComplForward = FALSE, revComplReverse = FALSE,
    elementsForward = "SUCV", elementsReverse = "SUCV",
    elementLengthsForward = c(1, 10, 18, -1),
    elementLengthsReverse = c(1, 8, 20, -1),
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
    nThreads = 1, chunkSize = 700
  )

  res <- do.call(digestFastqs, Ldef)

  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 287L + 105L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 6L + 2L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
  expect_equal(res$filterSummary$nbrRetained, 279L)

  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose", "fastqForward", "fastqReverse"))) {
    expect_equal(res$parameters[[nm]], Ldef[[nm]], ignore_attr = TRUE)
  }
  for (nm in c("fastqForward", "fastqReverse")) {
    expect_equal(res$parameters[[nm]], normalizePath(Ldef[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }
  

  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(res$summaryTable$nbrReads, res$summaryTable$maxNbrReads)
  expect_equal(sum(res$summaryTable$nbrReads == 2), 2L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 2]),
               sort(c("f.13.GAG_r.0.WT", "f.26.TAG_r.8.AGG")))
  expect_equal(sort(res$summaryTable$mutantNameCodon[res$summaryTable$nbrReads == 2]),
               sort(c("f.13.GAG_r.0.WT", "f.26.TAG_r.8.AGG")))
  expect_equal(sort(res$summaryTable$mutantNameBase[res$summaryTable$nbrReads == 2]),
               sort(c("f.39.G_r.0.WT", "f.76.T_f.77.A_r.22.A_r.23.G")))
  expect_equal(sort(res$summaryTable$mutantNameBaseHGVS[res$summaryTable$nbrReads == 2]),
               sort(c("f:c.39T>G_r:c", "f:c.76_77delinsTA_r:c.22_23delinsAG")))
  expect_equal(sort(res$summaryTable$mutantNameAAHGVS[res$summaryTable$nbrReads == 2]),
               sort(c("f:p.(Asp13Glu)_r:p", "f:p.(Leu26*)_r:p.(Val8Arg)")))
  expect_true(all(res$summaryTable$nbrReads == res$summaryTable$nbrUmis))
  
  expect_true(all(res$summaryTable$nbrReads == res$summaryTable$nbrUmis))

  ## Check that mutant naming worked (compare to manual matching)
  example_seq <- paste0("ACTGATACAACCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTG",
                        "CAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA_ATCGCCCGGCT",
                        "GGAGGAAAAAGTGGGCACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGC",
                        "CAACATGCTCAGGGAACAGGTGGCACAGCTT")
  expect_equal(res$summaryTable$mutantName[res$summaryTable$sequence == example_seq],
               "f.4.ACC_r.9.GGC")
  expect_equal(res$summaryTable$mutantNameAA[res$summaryTable$sequence == example_seq],
               "f.4.T_r.9.G")
  expect_equal(res$summaryTable$mutantNameCodon[res$summaryTable$sequence == example_seq],
               "f.4.ACC_r.9.GGC")
  expect_equal(res$summaryTable$mutantNameBase[res$summaryTable$sequence == example_seq],
               "f.10.A_f.11.C_r.25.G_r.26.G_r.27.C")
  expect_equal(res$summaryTable$mutantNameBaseHGVS[res$summaryTable$sequence == example_seq],
               "f:c.10_11delinsAC_r:c.25_27delinsGGC")
  expect_equal(res$summaryTable$mutantNameAAHGVS[res$summaryTable$sequence == example_seq],
               "f:p.(Leu4Thr)_r:p.(Lys9Gly)")

  expect_equal(sum(res$errorStatistics$nbrMatchForward + res$errorStatistics$nbrMismatchForward),
               nchar(Ldef$constantForward[1]) * res$filterSummary$nbrRetained)
  expect_equal(sum(res$errorStatistics$nbrMatchReverse + res$errorStatistics$nbrMismatchReverse),
               nchar(Ldef$constantReverse[1]) * res$filterSummary$nbrRetained)

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

## ----------------------------------------------------------------------------
## Multiple wildtype sequences
## ----------------------------------------------------------------------------
test_that("digestFastqs works as expected for trans experiments when multiple reference sequences are provided", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2,
    mergeForwardReverse = FALSE,
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE,
    revComplForward = FALSE, revComplReverse = FALSE,
    elementsForward = "SUCV", elementsReverse = "SUCV",
    elementLengthsForward = c(1, 10, 18, 96),
    elementLengthsReverse = c(1, 8, 20, 96),
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
    nbrMutatedBasesMaxForward = -1,
    nbrMutatedBasesMaxReverse = -1,
    forbiddenMutatedCodonsForward = "NNA",
    forbiddenMutatedCodonsReverse = "NNW",
    useTreeWTmatch = FALSE,
    mutatedPhredMinForward = 25.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = "=",
    constantMaxDistForward = -1,
    constantMaxDistReverse = -1,
    variableCollapseMaxDist = 0,
    variableCollapseMinReads = 0,
    variableCollapseMinRatio = 0,
    umiCollapseMaxDist = 0,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = TRUE,
    nThreads = 1, chunkSize = 1000
  )

  expect_output(res <- do.call(digestFastqs, Ldef))

  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 88L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 189L)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 149L + 5L + 34L + 25L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 1L + 2L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
  expect_equal(res$filterSummary$nbrRetained, 193L)

  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "nThreads",
                                    "forbiddenMutatedCodonsReverse", "verbose",
                                    "fastqForward", "fastqReverse"))) {
    expect_equal(res$parameters[[nm]], Ldef[[nm]], ignore_attr = TRUE)
  }
  for (nm in c("fastqForward", "fastqReverse")) {
    expect_equal(res$parameters[[nm]], normalizePath(Ldef[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }
  

  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(res$summaryTable$nbrReads, res$summaryTable$maxNbrReads)
  expect_equal(sum(res$summaryTable$nbrReads == 2), 1L)
  expect_equal(sort(res$summaryTable$mutantName[res$summaryTable$nbrReads == 2]),
               sort(c("forward=26=TAG_reverse1=8=AGG")))
  expect_equal(sort(res$summaryTable$mutantNameCodon[res$summaryTable$nbrReads == 2]),
               sort(c("forward=26=TAG_reverse1=8=AGG")))
  expect_equal(sort(res$summaryTable$mutantNameBase[res$summaryTable$nbrReads == 2]),
               sort(c("forward=76=T_forward=77=A_reverse1=22=A_reverse1=23=G")))
  expect_equal(sort(res$summaryTable$mutantNameBaseHGVS[res$summaryTable$nbrReads == 2]),
               sort(c("forward:c.76_77delinsTA_reverse1:c.22_23delinsAG")))
  expect_equal(sort(res$summaryTable$mutantNameAAHGVS[res$summaryTable$nbrReads == 2]),
               sort(c("forward:p.(Leu26*)_reverse1:p.(Val8Arg)")))
  expect_true(all(res$summaryTable$nbrReads == res$summaryTable$nbrUmis))

  ## Check that mutant naming worked (compare to manual matching)
  example_seq <- paste0("ACTGATACAACCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTG",
                        "CAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA_ATCGCCCGGCT",
                        "GGAGGAAAAAGTGGGCACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGC",
                        "CAACATGCTCAGGGAACAGGTGGCACAGCTT")
  expect_equal(res$summaryTable$mutantName[res$summaryTable$sequence == example_seq],
               "forward=4=ACC_reverse1=9=GGC")
  expect_equal(res$summaryTable$mutantNameAA[res$summaryTable$sequence == example_seq],
               "forward=4=T_reverse1=9=G")
  expect_equal(res$summaryTable$mutantNameCodon[res$summaryTable$sequence == example_seq],
               "forward=4=ACC_reverse1=9=GGC")
  expect_equal(res$summaryTable$mutantNameBase[res$summaryTable$sequence == example_seq],
               "forward=10=A_forward=11=C_reverse1=25=G_reverse1=26=G_reverse1=27=C")
  expect_equal(res$summaryTable$mutantNameBaseHGVS[res$summaryTable$sequence == example_seq],
               "forward:c.10_11delinsAC_reverse1:c.25_27delinsGGC")
  expect_equal(res$summaryTable$mutantNameAAHGVS[res$summaryTable$sequence == example_seq],
               "forward:p.(Leu4Thr)_reverse1:p.(Lys9Gly)")

  expect_equal(sum(res$errorStatistics$nbrMatchForward + res$errorStatistics$nbrMismatchForward),
               nchar(Ldef$constantForward[1]) * res$filterSummary$nbrRetained)
  expect_equal(sum(res$errorStatistics$nbrMatchReverse + res$errorStatistics$nbrMismatchReverse),
               nchar(Ldef$constantReverse[1]) * res$filterSummary$nbrRetained)

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

  ## Test that we get the same results with useTreeWTmatch = TRUE
  Ldef$useTreeWTmatch <- TRUE
  expect_output(res2 <- do.call(digestFastqs, Ldef))

  expect_equal(res$filterSummary$nbrTotal, res2$filterSummary$nbrTotal)
  expect_equal(res$filterSummary$f1_nbrAdapter, res2$filterSummary$f1_nbrAdapter)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, res2$filterSummary$f2_nbrNoPrimer)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, res2$filterSummary$f3_nbrReadWrongLength)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, res2$filterSummary$f4_nbrNoValidOverlap)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, res2$filterSummary$f5_nbrAvgVarQualTooLow)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, res2$filterSummary$f6_nbrTooManyNinVar)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, res2$filterSummary$f7_nbrTooManyNinUMI)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, res2$filterSummary$f8_nbrTooManyBestWTHits)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, res2$filterSummary$f9_nbrMutQualTooLow)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, res2$filterSummary$f10a_nbrTooManyMutCodons)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, res2$filterSummary$f10b_nbrTooManyMutBases)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, res2$filterSummary$f11_nbrForbiddenCodons)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, res2$filterSummary$f12_nbrTooManyMutConstant)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, res2$filterSummary$f13_nbrTooManyBestConstantHits)
  expect_equal(res$filterSummary$nbrRetained, res2$filterSummary$nbrRetained)

  expect_equal(res$summaryTable$nbrReads, res2$summaryTable$nbrReads)
  expect_equal(res$summaryTable$mutantName, res2$summaryTable$mutantName)
  expect_equal(res$summaryTable$sequence, res2$summaryTable$sequence)

})

## ----------------------------------------------------------------------------
## Only forward read
## ----------------------------------------------------------------------------
test_that("digestFastqs works as expected for experiments with only forward read", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = NULL,
    mergeForwardReverse = FALSE,
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE,
    revComplForward = FALSE, revComplReverse = FALSE,
    elementsForward = "SUCV", elementsReverse = "",
    elementLengthsForward = c(1, 10, 18, 96),
    elementLengthsReverse = numeric(0),
    adapterForward = "GGAAGAGCACACGTC",
    adapterReverse = "",
    primerForward = "",
    primerReverse = "",
    wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
    wildTypeReverse = "",
    constantForward = "AACCGGAGGAGGGAGCTG",
    constantReverse = "",
    avePhredMinForward = 20.0, avePhredMinReverse = 30.0,
    variableNMaxForward = 0, variableNMaxReverse = 0,
    umiNMax = 0,
    nbrMutatedCodonsMaxForward = 1,
    nbrMutatedCodonsMaxReverse = 1,
    nbrMutatedBasesMaxForward = -1,
    nbrMutatedBasesMaxReverse = -1,
    forbiddenMutatedCodonsForward = "NNA",
    forbiddenMutatedCodonsReverse = "NNW",
    useTreeWTmatch = FALSE,
    mutatedPhredMinForward = 25.0, mutatedPhredMinReverse = 0.0,
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
    nThreads = 1, chunkSize = 300
  )

  res <- do.call(digestFastqs, Ldef)

  L <- Ldef; L$fastqForward <- c(fqt1, fqt1); L$fastqReverse <- NULL
  res2 <- do.call(digestFastqs, L)

  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 297L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 0L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 211L)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 212L + 7L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 1L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
    expect_equal(res$filterSummary$nbrRetained, 272L)
    
    expect_equal(res$filterSummary * 2, res2$filterSummary)
    expect_equal(res$summaryTable$nbrReads * 2, res2$summaryTable$nbrReads)
    expect_equal(res$summaryTable$nbrUmis, res2$summaryTable$nbrUmis)
    
    for (nm in setdiff(names(Ldef), c("fastqReverse", "forbiddenMutatedCodonsForward",
                                      "forbiddenMutatedCodonsReverse", "verbose",
                                      "elementLengthsReverse", "fastqForward"))) {
        expect_equal(res$parameters[[nm]], Ldef[[nm]], ignore_attr = TRUE)
    }
    expect_equal(res$parameters[["fastqReverse"]], "", ignore_attr = TRUE)
    expect_equal(res$parameters[["fastqForward"]], 
                 normalizePath(Ldef[["fastqForward"]], mustWork = FALSE), 
                 ignore_attr = TRUE)

    
    expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
    expect_equal(res$summaryTable$nbrReads, res$summaryTable$maxNbrReads)
    expect_equal(sum(res$summaryTable$nbrReads == 2), 29L)
    
    ## Check the number of reads with a given number of mismatches
    expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 0]), 10L)
    expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 1]), 47L)
    expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 2]), 123L)
    expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == 3]), 92L)
    expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutCodons == 0]), 10L)
    expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutCodons == 1]), 262L)
    expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutAAs == 0]), 21L)
    expect_equal(sum(res$summaryTable$nbrReads[res$summaryTable$nbrMutAAs == 1]), 251L)
    
    ## Similarly but for the number of unique sequences
    expect_equal(sum(res$summaryTable$nbrMutBases == 0), 1L)
    expect_equal(sum(res$summaryTable$nbrMutBases == 1), 43L)
    expect_equal(sum(res$summaryTable$nbrMutBases == 2), 102L)
    expect_equal(sum(res$summaryTable$nbrMutBases == 3), 78L)
    expect_equal(sum(res$summaryTable$nbrMutCodons == 0), 1L)
    expect_equal(sum(res$summaryTable$nbrMutCodons == 1), 223L)
    expect_equal(sum(res$summaryTable$nbrMutAAs == 0), 10L)
    expect_equal(sum(res$summaryTable$nbrMutAAs == 1), 214L)
    
    
    ## Check that mutant naming worked (compare to manual matching)
    example_seq <- paste0("ACTGATACAACCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTG",
                          "CAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA")
    expect_equal(res$summaryTable$mutantName[res$summaryTable$sequence == example_seq],
                 "f.4.ACC")
    
    expect_equal(sum(res$errorStatistics$nbrMatchForward + res$errorStatistics$nbrMismatchForward),
                 nchar(Ldef$constantForward[1]) * res$filterSummary$nbrRetained)
    expect_equal(sum(res$errorStatistics$nbrMatchReverse + res$errorStatistics$nbrMismatchReverse),
                 nchar(Ldef$constantReverse[1]) * res$filterSummary$nbrRetained)
    
    expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 14], 187L)
    expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 22], 47L)
    expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 27], 314L)
    expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 33], 469L)
    expect_equal(res$errorStatistics$nbrMatchForward[res$errorStatistics$PhredQuality == 37], 3868L)
    
    expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 14], 6L)
    expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 27], 4L)
    expect_equal(res$errorStatistics$nbrMismatchForward[res$errorStatistics$PhredQuality == 37], 1L)
    
})

## ----------------------------------------------------------------------------
## Only forward read, provide a too short WT sequence
## ----------------------------------------------------------------------------
test_that("digestFastqs works as expected when the WT sequence is too short", {
    fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
    fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
    ## default arguments
    Ldef <- list(
        fastqForward = fqt1, fastqReverse = NULL,
        mergeForwardReverse = FALSE,
        minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE,
        revComplForward = FALSE, revComplReverse = FALSE,
        elementsForward = "SUCV", elementsReverse = "",
        elementLengthsForward = c(1, 10, 18, 96),
        elementLengthsReverse = numeric(0),
        adapterForward = "GGAAGAGCACACGTC",
        adapterReverse = "",
        primerForward = "",
        primerReverse = "",
        wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACT",
        wildTypeReverse = "",
        constantForward = "AACCGGAGGAGGGAGCTG",
        constantReverse = "",
        avePhredMinForward = 20.0, avePhredMinReverse = 30.0,
        variableNMaxForward = 0, variableNMaxReverse = 0,
        umiNMax = 0,
        nbrMutatedCodonsMaxForward = 1,
        nbrMutatedCodonsMaxReverse = 1,
        nbrMutatedBasesMaxForward = -1,
        nbrMutatedBasesMaxReverse = -1,
        forbiddenMutatedCodonsForward = "NNA",
        forbiddenMutatedCodonsReverse = "NNW",
        useTreeWTmatch = FALSE,
        mutatedPhredMinForward = 25.0, mutatedPhredMinReverse = 0.0,
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
        nThreads = 1, chunkSize = 300
    )
    
    res <- do.call(digestFastqs, Ldef)
    
    expect_equal(res$filterSummary$nbrTotal, 1000L)
    expect_equal(res$filterSummary$f1_nbrAdapter, 297L)
    expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
    expect_equal(res$filterSummary$f3_nbrReadWrongLength, 284L) ## generated as 1000-297-419 (all var seqs are too long)
    expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
    expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 0L)
    expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
    expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
    expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
    expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 0L)
    expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 419L) ## this value was generated and checked manually
    expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
    expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 0L)
    expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
    expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
    expect_equal(res$filterSummary$nbrRetained, 0L)
    
    expect_equal(nrow(res$summaryTable), 0L)
    
    for (nm in setdiff(names(Ldef), c("fastqReverse", "forbiddenMutatedCodonsForward",
                                      "forbiddenMutatedCodonsReverse", "verbose",
                                      "elementLengthsReverse", "fastqForward"))) {
        expect_equal(res$parameters[[nm]], Ldef[[nm]], ignore_attr = TRUE)
    }
    expect_equal(res$parameters[["fastqReverse"]], "", ignore_attr = TRUE)
    expect_equal(res$parameters[["fastqForward"]], 
                 normalizePath(Ldef[["fastqForward"]], mustWork = FALSE), 
                 ignore_attr = TRUE)
    
})

## ----------------------------------------------------------------------------
## Save filtered reads to fastq files
## ----------------------------------------------------------------------------
test_that("writing filtered reads to file works", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")

  outfile1 <- file.path(tempdir(), "mutscan_test_1.fastq.gz")
  outfile2 <- gsub("_1.fastq", "_2.fastq", outfile1)

  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2,
    mergeForwardReverse = FALSE,
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE,
    revComplForward = FALSE, revComplReverse = FALSE,
    elementsForward = "SUCV", elementsReverse = "SUCV",
    elementLengthsForward = c(1, 10, 18, 96),
    elementLengthsReverse = c(1, 8, 20, 96),
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
    filteredReadsFastqForward = outfile1,
    filteredReadsFastqReverse = outfile2,
    maxNReads = -1, verbose = FALSE,
    nThreads = 1, chunkSize = 1000
  )

  res <- do.call(digestFastqs, Ldef)

  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 287L + 105L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 6L + 2L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
  expect_equal(res$filterSummary$nbrRetained, 279L)

  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose", "fastqForward", "fastqReverse",
                                    "filteredReadsFastqForward", 
                                    "filteredReadsFastqReverse"))) {
    expect_equal(res$parameters[[nm]], Ldef[[nm]], ignore_attr = TRUE)
  }
  for (nm in c("fastqForward", "fastqReverse", "filteredReadsFastqForward",
               "filteredReadsFastqReverse")) {
    expect_equal(normalizePath(res$parameters[[nm]], mustWork = FALSE), 
                 normalizePath(Ldef[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }
  

  # Biostrings 2.63.1 on arm64 produces the following warning:
  # Warning message:
  # In XStringSet("DNA", x, start = start, end = end, width = width,  :
  #   metadata columns on input DNAStringSet object were dropped
  suppressWarnings({
      out1 <- Biostrings::readQualityScaledDNAStringSet(outfile1)
      out2 <- Biostrings::readQualityScaledDNAStringSet(outfile2)
  })

  expect_equal(length(out1), length(out2))
  expect_equal(names(out1), gsub(" 2", " 1", names(out2)))
  expect_equal(length(out1), res$filterSummary$nbrTotal - res$filterSummary$nbrRetained)
  expect_equal(length(grep("adapter", names(out1))), res$filterSummary$f1_nbrAdapter)
  expect_equal(length(grep("noPrimer", names(out1))), res$filterSummary$f2_nbrNoPrimer)
  expect_equal(length(grep("readWrongLength", names(out1))), res$filterSummary$f3_nbrReadWrongLength)
  expect_equal(length(grep("noValidOverlap", names(out1))), res$filterSummary$f4_nbrNoValidOverlap)
  expect_equal(length(grep("avgVarQualTooLow", names(out1))), res$filterSummary$f5_nbrAvgVarQualTooLow)
  expect_equal(length(grep("tooManyNinVar", names(out1))), res$filterSummary$f6_nbrTooManyNinVar)
  expect_equal(length(grep("tooManyNinUMI", names(out1))), res$filterSummary$f7_nbrTooManyNinUMI)
  expect_equal(length(grep("tooManyBestWTHits", names(out1))), res$filterSummary$f8_nbrTooManyBestWTHits)
  expect_equal(length(grep("mutQualTooLow", names(out1))),
               res$filterSummary$f9_nbrMutQualTooLow + res$filterSummary$f10a_nbrTooManyMutCodons +
                 res$filterSummary$f11_nbrForbiddenCodons)
  expect_equal(length(grep("tooManyMutConstant", names(out1))), res$filterSummary$f12_nbrTooManyMutConstant)
  expect_equal(length(grep("tooManyBestConstantHits", names(out1))), res$filterSummary$f13_nbrTooManyBestConstantHits)

  ## Test that we get the same results with useTreeWTmatch = TRUE
  Ldef$useTreeWTmatch <- TRUE
  res2 <- do.call(digestFastqs, Ldef)

  expect_equal(res$filterSummary$nbrTotal, res2$filterSummary$nbrTotal)
  expect_equal(res$filterSummary$f1_nbrAdapter, res2$filterSummary$f1_nbrAdapter)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, res2$filterSummary$f2_nbrNoPrimer)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, res2$filterSummary$f3_nbrReadWrongLength)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, res2$filterSummary$f4_nbrNoValidOverlap)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, res2$filterSummary$f5_nbrAvgVarQualTooLow)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, res2$filterSummary$f6_nbrTooManyNinVar)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, res2$filterSummary$f7_nbrTooManyNinUMI)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, res2$filterSummary$f8_nbrTooManyBestWTHits)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, res2$filterSummary$f9_nbrMutQualTooLow)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, res2$filterSummary$f10a_nbrTooManyMutCodons)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, res2$filterSummary$f10b_nbrTooManyMutBases)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, res2$filterSummary$f11_nbrForbiddenCodons)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, res2$filterSummary$f12_nbrTooManyMutConstant)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, res2$filterSummary$f13_nbrTooManyBestConstantHits)
  expect_equal(res$filterSummary$nbrRetained, res2$filterSummary$nbrRetained)

  expect_equal(res$summaryTable$nbrReads, res2$summaryTable$nbrReads)
  expect_equal(res$summaryTable$mutantName, res2$summaryTable$mutantName)
  expect_equal(res$summaryTable$sequence, res2$summaryTable$sequence)

})

test_that("digestFastqs works as expected for trans experiments, if collapsing to WT", {
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2,
    mergeForwardReverse = FALSE,
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE,
    revComplForward = FALSE, revComplReverse = FALSE,
    elementsForward = "SUCVV", elementsReverse = "SUCVV",
    elementLengthsForward = c(1, 10, 18, 80, 16),
    elementLengthsReverse = c(1, 8, 20, 16, 80),
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
    nbrMutatedBasesMaxForward = -1,
    nbrMutatedBasesMaxReverse = -1,
    forbiddenMutatedCodonsForward = "NNW",
    forbiddenMutatedCodonsReverse = "NNW",
    useTreeWTmatch = FALSE,
    collapseToWTForward = TRUE,
    collapseToWTReverse = TRUE, 
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
  expect_equal(res$filterSummary$f1_nbrAdapter, 314L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 192L + 95L + 68L + 37L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 6L + 2L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
  expect_equal(res$filterSummary$nbrRetained, 279L)
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose", "fastqForward", "fastqReverse"))) {
    expect_equal(res$parameters[[nm]], Ldef[[nm]], ignore_attr = TRUE)
  }
  for (nm in c("fastqForward", "fastqReverse")) {
    expect_equal(res$parameters[[nm]], normalizePath(Ldef[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }
  
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(res$summaryTable$nbrReads, res$summaryTable$maxNbrReads)
  expect_equal(nrow(res$summaryTable), 1L)
  expect_true(all(res$summaryTable$nbrReads == res$summaryTable$nbrUmis))
  
  ## Check the number of reads with a given number of mismatches
  expect_equal(res$summaryTable$nbrReads[res$summaryTable$nbrMutBases == "1,2,3,4,5,6"], 279L)
  expect_equal(res$summaryTable$nbrReads[res$summaryTable$nbrMutCodons == "1,2"], 279L)
  expect_equal(res$summaryTable$nbrReads[res$summaryTable$nbrMutAAs == "0,1,2"], 279L)
  expect_equal(res$summaryTable$mutationTypes, "nonsynonymous,silent,stop")
  expect_equal(res$summaryTable$mutantNameAA, "f_r")
  expect_equal(res$summaryTable$mutantNameBase, "f_r")
  expect_equal(res$summaryTable$mutantNameCodon, "f_r")
  expect_equal(res$summaryTable$mutantNameBaseHGVS, "f:c_r:c")
  expect_equal(res$summaryTable$mutantNameAAHGVS, "f:p_r:p")
  
  ## check that variable segment lengths are reported correctly
  expect_equal(res$summaryTable$varLengths, "80,16_16,80")
})

