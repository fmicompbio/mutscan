test_that("parallel processing works (cis)", {
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
    nThreads = 1, chunkSize = 1
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  ## Increase number of threads
  Ldef1 <- Ldef; Ldef1$nThreads <- 3
  ## Capture warnings from missing OpenMP
  a <- capture_warnings({
      res1 <- do.call(digestFastqs, Ldef1)
  })
  expect_equal(sum(grepl("OpenMP", a)), length(a))

  expect_equal(res$filterSummary$nbrTotal, res1$filterSummary$nbrTotal)
  expect_equal(res$filterSummary$f1_nbrAdapter, res1$filterSummary$f1_nbrAdapter)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, res1$filterSummary$f2_nbrNoPrimer)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, res1$filterSummary$f3_nbrReadWrongLength)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, res1$filterSummary$f4_nbrNoValidOverlap)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, res1$filterSummary$f5_nbrAvgVarQualTooLow)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, res1$filterSummary$f6_nbrTooManyNinVar)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, res1$filterSummary$f7_nbrTooManyNinUMI)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, res1$filterSummary$f8_nbrTooManyBestWTHits)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, res1$filterSummary$f9_nbrMutQualTooLow)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, res1$filterSummary$f10a_nbrTooManyMutCodons)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, res1$filterSummary$f10b_nbrTooManyMutBases)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, res1$filterSummary$f11_nbrForbiddenCodons)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, res1$filterSummary$f12_nbrTooManyMutConstant)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, res1$filterSummary$f13_nbrTooManyBestConstantHits)
  expect_equal(res$filterSummary$nbrRetained, res1$filterSummary$nbrRetained)
  
  expect_equal(res$summaryTable$nbrReads, res1$summaryTable$nbrReads)
  expect_equal(res$summaryTable$mutantName, res1$summaryTable$mutantName)
  expect_equal(res$summaryTable$sequence, res1$summaryTable$sequence)
  
  ## Increase chunk size
  Ldef1 <- Ldef; Ldef1$nThreads <- 3; Ldef1$chunkSize <- 100
  ## Capture warnings from missing OpenMP
  a <- capture_warnings({
      res1 <- do.call(digestFastqs, Ldef1)
  })
  expect_equal(sum(grepl("OpenMP", a)), length(a))
  
  expect_equal(res$filterSummary$nbrTotal, res1$filterSummary$nbrTotal)
  expect_equal(res$filterSummary$f1_nbrAdapter, res1$filterSummary$f1_nbrAdapter)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, res1$filterSummary$f2_nbrNoPrimer)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, res1$filterSummary$f3_nbrReadWrongLength)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, res1$filterSummary$f4_nbrNoValidOverlap)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, res1$filterSummary$f5_nbrAvgVarQualTooLow)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, res1$filterSummary$f6_nbrTooManyNinVar)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, res1$filterSummary$f7_nbrTooManyNinUMI)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, res1$filterSummary$f8_nbrTooManyBestWTHits)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, res1$filterSummary$f9_nbrMutQualTooLow)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, res1$filterSummary$f10a_nbrTooManyMutCodons)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, res1$filterSummary$f10b_nbrTooManyMutBases)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, res1$filterSummary$f11_nbrForbiddenCodons)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, res1$filterSummary$f12_nbrTooManyMutConstant)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, res1$filterSummary$f13_nbrTooManyBestConstantHits)
  expect_equal(res$filterSummary$nbrRetained, res1$filterSummary$nbrRetained)
  
  expect_equal(res$summaryTable$nbrReads, res1$summaryTable$nbrReads)
  expect_equal(res$summaryTable$mutantName, res1$summaryTable$mutantName)
  expect_equal(res$summaryTable$sequence, res1$summaryTable$sequence)
  
  ## With tree matching
  Ldef1 <- Ldef; Ldef1$useTreeWTmatch <- TRUE; Ldef1$nThreads <- 3; Ldef1$chunkSize <- 100
  ## Capture warnings from missing OpenMP
  a <- capture_warnings({
      res1 <- do.call(digestFastqs, Ldef1)
  })
  expect_equal(sum(grepl("OpenMP", a)), length(a))
  
  expect_equal(res$filterSummary$nbrTotal, res1$filterSummary$nbrTotal)
  expect_equal(res$filterSummary$f1_nbrAdapter, res1$filterSummary$f1_nbrAdapter)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, res1$filterSummary$f2_nbrNoPrimer)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, res1$filterSummary$f3_nbrReadWrongLength)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, res1$filterSummary$f4_nbrNoValidOverlap)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, res1$filterSummary$f5_nbrAvgVarQualTooLow)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, res1$filterSummary$f6_nbrTooManyNinVar)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, res1$filterSummary$f7_nbrTooManyNinUMI)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, res1$filterSummary$f8_nbrTooManyBestWTHits)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, res1$filterSummary$f9_nbrMutQualTooLow)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, res1$filterSummary$f10a_nbrTooManyMutCodons)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, res1$filterSummary$f10b_nbrTooManyMutBases)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, res1$filterSummary$f11_nbrForbiddenCodons)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, res1$filterSummary$f12_nbrTooManyMutConstant)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, res1$filterSummary$f13_nbrTooManyBestConstantHits)
  expect_equal(res$filterSummary$nbrRetained, res1$filterSummary$nbrRetained)
  
  expect_equal(res$summaryTable$nbrReads, res1$summaryTable$nbrReads)
  expect_equal(res$summaryTable$mutantName, res1$summaryTable$mutantName)
  expect_equal(res$summaryTable$sequence, res1$summaryTable$sequence)
})

test_that("parallel processing works (trans)", {
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
    nThreads = 1, chunkSize = 1
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  ## Increase number of threads
  Ldef1 <- Ldef; Ldef1$nThreads <- 3
  ## Capture warnings from missing OpenMP
  a <- capture_warnings({
      res1 <- do.call(digestFastqs, Ldef1)
  })
  expect_equal(sum(grepl("OpenMP", a)), length(a))
  
  expect_equal(res$filterSummary$nbrTotal, res1$filterSummary$nbrTotal)
  expect_equal(res$filterSummary$f1_nbrAdapter, res1$filterSummary$f1_nbrAdapter)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, res1$filterSummary$f2_nbrNoPrimer)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, res1$filterSummary$f3_nbrReadWrongLength)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, res1$filterSummary$f4_nbrNoValidOverlap)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, res1$filterSummary$f5_nbrAvgVarQualTooLow)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, res1$filterSummary$f6_nbrTooManyNinVar)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, res1$filterSummary$f7_nbrTooManyNinUMI)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, res1$filterSummary$f8_nbrTooManyBestWTHits)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, res1$filterSummary$f9_nbrMutQualTooLow)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, res1$filterSummary$f10a_nbrTooManyMutCodons)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, res1$filterSummary$f10b_nbrTooManyMutBases)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, res1$filterSummary$f11_nbrForbiddenCodons)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, res1$filterSummary$f12_nbrTooManyMutConstant)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, res1$filterSummary$f13_nbrTooManyBestConstantHits)
  expect_equal(res$filterSummary$nbrRetained, res1$filterSummary$nbrRetained)
  
  expect_equal(res$summaryTable$nbrReads, res1$summaryTable$nbrReads)
  expect_equal(res$summaryTable$mutantName, res1$summaryTable$mutantName)
  expect_equal(res$summaryTable$sequence, res1$summaryTable$sequence)
  
  ## Increase chunk size
  Ldef1 <- Ldef; Ldef1$nThreads <- 3; Ldef1$chunkSize <- 100
  ## Capture warnings from missing OpenMP
  a <- capture_warnings({
      res1 <- do.call(digestFastqs, Ldef1)
  })
  expect_equal(sum(grepl("OpenMP", a)), length(a))
  
  expect_equal(res$filterSummary$nbrTotal, res1$filterSummary$nbrTotal)
  expect_equal(res$filterSummary$f1_nbrAdapter, res1$filterSummary$f1_nbrAdapter)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, res1$filterSummary$f2_nbrNoPrimer)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, res1$filterSummary$f3_nbrReadWrongLength)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, res1$filterSummary$f4_nbrNoValidOverlap)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, res1$filterSummary$f5_nbrAvgVarQualTooLow)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, res1$filterSummary$f6_nbrTooManyNinVar)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, res1$filterSummary$f7_nbrTooManyNinUMI)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, res1$filterSummary$f8_nbrTooManyBestWTHits)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, res1$filterSummary$f9_nbrMutQualTooLow)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, res1$filterSummary$f10a_nbrTooManyMutCodons)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, res1$filterSummary$f10b_nbrTooManyMutBases)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, res1$filterSummary$f11_nbrForbiddenCodons)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, res1$filterSummary$f12_nbrTooManyMutConstant)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, res1$filterSummary$f13_nbrTooManyBestConstantHits)
  expect_equal(res$filterSummary$nbrRetained, res1$filterSummary$nbrRetained)
  
  expect_equal(res$summaryTable$nbrReads, res1$summaryTable$nbrReads)
  expect_equal(res$summaryTable$mutantName, res1$summaryTable$mutantName)
  expect_equal(res$summaryTable$sequence, res1$summaryTable$sequence)
  
  ## With tree matching
  Ldef1 <- Ldef; Ldef1$useTreeWTmatch <- TRUE; Ldef1$nThreads <- 3; Ldef1$chunkSize <- 100
  ## Capture warnings from missing OpenMP
  a <- capture_warnings({
      res1 <- do.call(digestFastqs, Ldef1)
  })
  expect_equal(sum(grepl("OpenMP", a)), length(a))
  
  expect_equal(res$filterSummary$nbrTotal, res1$filterSummary$nbrTotal)
  expect_equal(res$filterSummary$f1_nbrAdapter, res1$filterSummary$f1_nbrAdapter)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, res1$filterSummary$f2_nbrNoPrimer)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, res1$filterSummary$f3_nbrReadWrongLength)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, res1$filterSummary$f4_nbrNoValidOverlap)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, res1$filterSummary$f5_nbrAvgVarQualTooLow)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, res1$filterSummary$f6_nbrTooManyNinVar)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, res1$filterSummary$f7_nbrTooManyNinUMI)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, res1$filterSummary$f8_nbrTooManyBestWTHits)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, res1$filterSummary$f9_nbrMutQualTooLow)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, res1$filterSummary$f10a_nbrTooManyMutCodons)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, res1$filterSummary$f10b_nbrTooManyMutBases)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, res1$filterSummary$f11_nbrForbiddenCodons)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, res1$filterSummary$f12_nbrTooManyMutConstant)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, res1$filterSummary$f13_nbrTooManyBestConstantHits)
  expect_equal(res$filterSummary$nbrRetained, res1$filterSummary$nbrRetained)
  
  expect_equal(res$summaryTable$nbrReads, res1$summaryTable$nbrReads)
  expect_equal(res$summaryTable$mutantName, res1$summaryTable$mutantName)
  expect_equal(res$summaryTable$sequence, res1$summaryTable$sequence)
  
})
