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
    elementsForward = "SUCV", elementsReverse = "SUCV",
    elementLengthsForward = c(1, 10, 18, 96),
    elementLengthsReverse = c(1, 8, 20, 96),
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
    nbrMutatedBasesMaxForward = -1,
    nbrMutatedBasesMaxReverse = -1,
    forbiddenMutatedCodonsForward = "NNW",
    forbiddenMutatedCodonsReverse = "NNW",
    useTreeWTmatch = FALSE,
    mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".",
    constantMaxDistForward = -1,
    constantMaxDistReverse = -1,
    variableCollapseMaxDist = 6, 
    umiCollapseMaxDist = 4, 
    variableCollapseMinReads = 0,
    variableCollapseMinRatio = 0,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  L <- Ldef; L$variableCollapseMaxDist = 0; L$umiCollapseMaxDist = 0; res0 <- do.call(digestFastqs, L)
  expect_equal(res$summaryTable$maxNbrReads, res0$summaryTable$nbrReads[match(res$summaryTable$mutantName, 
                                                                              res0$summaryTable$mutantName)])
  
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
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 0L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 0L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
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
    elementsForward = "SUCV", elementsReverse = "SUCV",
    elementLengthsForward = c(1, 10, 18, 96),
    elementLengthsReverse = c(1, 8, 20, 96),
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
    nbrMutatedBasesMaxForward = -1,
    nbrMutatedBasesMaxReverse = -1,
    forbiddenMutatedCodonsForward = "NNW",
    forbiddenMutatedCodonsReverse = "NNW",
    useTreeWTmatch = FALSE, 
    mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".",
    constantMaxDistForward = -1,
    constantMaxDistReverse = -1,
    variableCollapseMaxDist = 500, 
    umiCollapseMaxDist = 100, 
    variableCollapseMinReads = 0,
    variableCollapseMinRatio = 0,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = FALSE
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
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 0L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 0L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
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
    elementsForward = "SUCV", elementsReverse = "SUCV",
    elementLengthsForward = c(1, 10, 18, 96),
    elementLengthsReverse = c(1, 8, 20, 96),
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
    nbrMutatedBasesMaxForward = -1,
    nbrMutatedBasesMaxReverse = -1,
    forbiddenMutatedCodonsForward = "NNW",
    forbiddenMutatedCodonsReverse = "NNW",
    useTreeWTmatch = FALSE, 
    mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".",
    constantMaxDistForward = -1,
    constantMaxDistReverse = -1,
    variableCollapseMaxDist = 10, 
    umiCollapseMaxDist = 5, 
    variableCollapseMinReads = 0,
    variableCollapseMinRatio = 0, 
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  expect_equal(res$filterSummary$nbrTotal, 1000L)
  expect_equal(res$filterSummary$f1_nbrAdapter, 297L)
  expect_equal(res$filterSummary$f2_nbrNoPrimer, 0L)
  expect_equal(res$filterSummary$f3_nbrReadWrongLength, 0L)
  expect_equal(res$filterSummary$f4_nbrNoValidOverlap, 0L)
  expect_equal(res$filterSummary$f5_nbrAvgVarQualTooLow, 0L)
  expect_equal(res$filterSummary$f6_nbrTooManyNinVar, 0L)
  expect_equal(res$filterSummary$f7_nbrTooManyNinUMI, 0L)
  expect_equal(res$filterSummary$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(res$filterSummary$f9_nbrMutQualTooLow, 0L)
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 0L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 0L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
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
    elementsForward = "SUCV", elementsReverse = "SUCV",
    elementLengthsForward = c(1, 10, 18, 96),
    elementLengthsReverse = c(1, 8, 20, 96),
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
    nbrMutatedBasesMaxForward = -1,
    nbrMutatedBasesMaxReverse = -1,
    forbiddenMutatedCodonsForward = "NNW",
    forbiddenMutatedCodonsReverse = "NNW",
    useTreeWTmatch = FALSE, 
    mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".",
    constantMaxDistForward = -1,
    constantMaxDistReverse = -1,
    variableCollapseMaxDist = 2, 
    umiCollapseMaxDist = 0,
    variableCollapseMinReads = 1.5,
    variableCollapseMinRatio = 0,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  L <- Ldef; L$variableCollapseMaxDist = 0; L$umiCollapseMaxDist = 0; res0 <- do.call(digestFastqs, L)
  expect_equal(res$summaryTable$maxNbrReads, res0$summaryTable$nbrReads[match(res$summaryTable$mutantName, 
                                                                              res0$summaryTable$mutantName)])
  
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
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 0L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 0L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
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

test_that("digestFastqs works as expected for trans experiments, when similar sequences are collapsed (only high ratio)", {
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
    wildTypeForward = "",
    wildTypeReverse = "", 
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
    variableCollapseMaxDist = 3, 
    umiCollapseMaxDist = 0, 
    variableCollapseMinReads = 1, 
    variableCollapseMinRatio = 1.5,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = FALSE
  )
  
  res <- do.call(digestFastqs, Ldef)
  
  L <- Ldef; L$variableCollapseMaxDist = 0; L$umiCollapseMaxDist = 0
  L$variableCollapseMinRatio <- 0; res0 <- do.call(digestFastqs, L)
  expect_equal(res$summaryTable$maxNbrReads, res0$summaryTable$nbrReads[match(res$summaryTable$mutantName, 
                                                                              res0$summaryTable$mutantName)])
  
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
  expect_equal(res$filterSummary$f10a_nbrTooManyMutCodons, 0L)
  expect_equal(res$filterSummary$f10b_nbrTooManyMutBases, 0L)
  expect_equal(res$filterSummary$f11_nbrForbiddenCodons, 0L)
  expect_equal(res$filterSummary$f12_nbrTooManyMutConstant, 0L)
  expect_equal(res$filterSummary$f13_nbrTooManyBestConstantHits, 0L)
  expect_equal(res$filterSummary$nbrRetained, 679L)
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose"))) {
    if (nm == "variableCollapseMinReads") {
      expect_equivalent(res$parameters[[nm]], as.integer(ceiling(Ldef[[nm]])))
    } else {
      expect_equivalent(res$parameters[[nm]], Ldef[[nm]])
    }
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(nrow(res$summaryTable), 656L)
  expect_equal(sum(res$summaryTable$nbrUmis), 679)
})
