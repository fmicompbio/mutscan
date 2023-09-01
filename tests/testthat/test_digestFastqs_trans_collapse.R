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
    umiCollapseMaxDist = 4,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = TRUE,
    nThreads = 1, chunkSize = 1000, 
    maxReadLength = 1024
  )

  res <- do.call(digestFastqs, Ldef)
  
  ## First check that we get deprecation warnings if deprecated arguments are 
  ## specified. Results should be the same as if these arguments are not used
  lifecycle::expect_deprecated({
      resdep1 <- do.call(digestFastqs, c(Ldef, list(variableCollapseMaxDist = 2)))})
  lifecycle::expect_deprecated({
      resdep2 <- do.call(digestFastqs, c(Ldef, list(variableCollapseMinReads = 2)))})
  lifecycle::expect_deprecated({
      resdep3 <- do.call(digestFastqs, c(Ldef, list(variableCollapseMinRatio = 2)))})
  
  expect_equal(res$filterSummary, resdep1$filterSummary)
  expect_equal(res$summaryTable, resdep1$summaryTable)
  expect_equal(res$errorStatistics, resdep1$errorStatistics)
  expect_equal(res$filterSummary, resdep2$filterSummary)
  expect_equal(res$summaryTable, resdep2$summaryTable)
  expect_equal(res$errorStatistics, resdep2$errorStatistics)
  expect_equal(res$filterSummary, resdep3$filterSummary)
  expect_equal(res$summaryTable, resdep3$summaryTable)
  expect_equal(res$errorStatistics, resdep3$errorStatistics)
  
  ## Summarize single sample and collapse
  se <- summarizeExperiment(list(s1 = res), coldata = data.frame(Name = "s1"), 
                            countType = "reads")
  secoll <- collapseMutantsBySimilarity(se, assayName = "counts", 
                                        collapseMaxDist = 6, 
                                        collapseMinScore = 0, collapseMinRatio = 0, 
                                        verbose = FALSE)
  seumi <- summarizeExperiment(list(s1 = res), coldata = data.frame(Name = "s1"), 
                               countType = "umis")
  secollumi <- collapseMutantsBySimilarity(seumi, assayName = "counts", 
                                           scoreMethod = "rowMean",
                                           collapseMaxDist = 6, 
                                           collapseMinScore = 0, collapseMinRatio = 0, 
                                           verbose = TRUE)

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
  
  expect_equal(SummarizedExperiment::colData(secoll)$nbrTotal, 1000L)
  expect_equal(SummarizedExperiment::colData(secoll)$f1_nbrAdapter, 314L)
  expect_equal(SummarizedExperiment::colData(secoll)$f2_nbrNoPrimer, 0L)
  expect_equal(SummarizedExperiment::colData(secoll)$f3_nbrReadWrongLength, 0L)
  expect_equal(SummarizedExperiment::colData(secoll)$f4_nbrNoValidOverlap, 0L)
  expect_equal(SummarizedExperiment::colData(secoll)$f5_nbrAvgVarQualTooLow, 7L)
  expect_equal(SummarizedExperiment::colData(secoll)$f6_nbrTooManyNinVar, 0L)
  expect_equal(SummarizedExperiment::colData(secoll)$f7_nbrTooManyNinUMI, 0L)
  expect_equal(SummarizedExperiment::colData(secoll)$f8_nbrTooManyBestWTHits, 0L)
  expect_equal(SummarizedExperiment::colData(secoll)$f9_nbrMutQualTooLow, 0L)
  expect_equal(SummarizedExperiment::colData(secoll)$f10a_nbrTooManyMutCodons, 0L)
  expect_equal(SummarizedExperiment::colData(secoll)$f10b_nbrTooManyMutBases, 0L)
  expect_equal(SummarizedExperiment::colData(secoll)$f11_nbrForbiddenCodons, 0L)
  expect_equal(SummarizedExperiment::colData(secoll)$f12_nbrTooManyMutConstant, 0L)
  expect_equal(SummarizedExperiment::colData(secoll)$f13_nbrTooManyBestConstantHits, 0L)
  expect_equal(SummarizedExperiment::colData(secoll)$nbrRetained, 679L)

  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose", "fastqForward", "fastqReverse"))) {
    expect_equal(res$parameters[[nm]], Ldef[[nm]], ignore_attr = TRUE)
  }
  for (nm in c("fastqForward", "fastqReverse")) {
    expect_equal(res$parameters[[nm]], normalizePath(Ldef[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(sum(SummarizedExperiment::assay(secoll, "counts")), 
               res$filterSummary$nbrRetained)
  expect_equal(nrow(res$summaryTable), 677L)
  expect_equal(sum(res$summaryTable$nbrUmis), 679L)
  expect_equal(nrow(secoll), 294L)
  expect_equal(nrow(secollumi), 294L)
  expect_true(all(res$summaryTable$varLengths == "96_96"))
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
    umiCollapseMaxDist = 100,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = FALSE,
    nThreads = 1, chunkSize = 1000, 
    maxReadLength = 125
  )

  res <- do.call(digestFastqs, Ldef)
  se <- summarizeExperiment(list(s1 = res), coldata = data.frame(Name = "s1"), 
                            countType = "reads")
  secoll <- collapseMutantsBySimilarity(se, assayName = "counts", 
                                        collapseMaxDist = 500, 
                                        collapseMinScore = 0, collapseMinRatio = 0, 
                                        verbose = FALSE)
  seumi <- summarizeExperiment(list(s1 = res), coldata = data.frame(Name = "s1"), 
                               countType = "umis")
  secollumi <- collapseMutantsBySimilarity(seumi, assayName = "counts", 
                                           collapseMaxDist = 500, 
                                           collapseMinScore = 0, collapseMinRatio = 0, 
                                           verbose = FALSE)

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

  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose", "fastqForward", "fastqReverse"))) {
    expect_equal(res$parameters[[nm]], Ldef[[nm]], ignore_attr = TRUE)
  }
  for (nm in c("fastqForward", "fastqReverse")) {
    expect_equal(res$parameters[[nm]], normalizePath(Ldef[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(nrow(res$summaryTable), 677L)
  expect_equal(sum(res$summaryTable$nbrUmis), 677L)
  expect_equal(nrow(secoll), 1L)
  expect_equal(nrow(secollumi), 1L)
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
    umiCollapseMaxDist = 5,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = FALSE,
    nThreads = 1, chunkSize = 1000, 
    maxReadLength = 1024
  )

  res <- do.call(digestFastqs, Ldef)
  se <- summarizeExperiment(list(s1 = res), coldata = data.frame(Name = "s1"), 
                            countType = "reads")
  secoll <- collapseMutantsBySimilarity(se, assayName = "counts", 
                                        collapseMaxDist = 10, 
                                        collapseMinScore = 0, collapseMinRatio = 0, 
                                        verbose = FALSE)
  seumi <- summarizeExperiment(list(s1 = res), coldata = data.frame(Name = "s1"), 
                               countType = "umis")
  secollumi <- collapseMutantsBySimilarity(seumi, assayName = "counts", 
                                           collapseMaxDist = 10, 
                                           collapseMinScore = 0, collapseMinRatio = 0, 
                                           verbose = FALSE)
  
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

  for (nm in setdiff(names(Ldef), c("fastqReverse", "forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose", "fastqForward"))) {
    expect_equal(res$parameters[[nm]], Ldef[[nm]], ignore_attr = TRUE)
  }
  for (nm in c("fastqForward")) {
    expect_equal(res$parameters[[nm]], normalizePath(Ldef[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(nrow(res$summaryTable), 613L)
  expect_equal(sum(res$summaryTable$nbrUmis), 687L)
  expect_equal(nrow(se), 613L)
  expect_equal(nrow(secoll), 52L)
  expect_equal(nrow(seumi), 613L)
  expect_equal(nrow(secollumi), 52L)
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
    umiCollapseMaxDist = 0,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = FALSE,
    nThreads = 1, chunkSize = 500, 
    maxReadLength = 1024
  )

  res <- do.call(digestFastqs, Ldef)
  se <- summarizeExperiment(list(s1 = res), coldata = data.frame(Name = "s1"), 
                            countType = "reads")
  secoll <- collapseMutantsBySimilarity(se, assayName = "counts", 
                                        collapseMaxDist = 2, 
                                        collapseMinScore = 1.5, collapseMinRatio = 0, 
                                        verbose = FALSE)
  seumi <- summarizeExperiment(list(s1 = res), coldata = data.frame(Name = "s1"), 
                               countType = "umis")
  secollumi <- collapseMutantsBySimilarity(seumi, assayName = "counts", 
                                           collapseMaxDist = 2, 
                                           collapseMinScore = 1.5, collapseMinRatio = 0, 
                                           verbose = FALSE)

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

  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "nThreads",
                                    "forbiddenMutatedCodonsReverse", "verbose",
                                    "fastqForward", "fastqReverse"))) {
    if (nm == "variableCollapseMinReads") {
      expect_equal(res$parameters[[nm]], as.integer(ceiling(Ldef[[nm]])),
                   ignore_attr = TRUE)
    } else {
      expect_equal(res$parameters[[nm]], Ldef[[nm]], ignore_attr = TRUE)
    }
  }
  for (nm in c("fastqForward", "fastqReverse")) {
    expect_equal(res$parameters[[nm]], normalizePath(Ldef[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }

  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(nrow(res$summaryTable), 677L)
  expect_equal(sum(res$summaryTable$nbrUmis), 679L)
  expect_equal(nrow(secoll), 673L)
  expect_equal(nrow(secollumi), 673L)
  
  ## Don't consider mutations here since we're collapsing (i.e., we have no wildtype)
  expect_equal(sum(res$summaryTable$nbrMutBases == 0), nrow(res$summaryTable))
  expect_equal(sum(res$summaryTable$nbrMutCodons == 0), nrow(res$summaryTable))
  expect_equal(sum(res$summaryTable$nbrMutAAs == 0), nrow(res$summaryTable))
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
    umiCollapseMaxDist = 0,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = FALSE,
    nThreads = 1, chunkSize = 1000, 
    maxReadLength = 1024
  )

  res <- do.call(digestFastqs, Ldef)
  se <- summarizeExperiment(list(s1 = res), coldata = data.frame(Name = "s1"), 
                            countType = "reads")
  secoll <- collapseMutantsBySimilarity(se, assayName = "counts", 
                                        collapseMaxDist = 3, 
                                        collapseMinScore = 1, collapseMinRatio = 1.5, 
                                        verbose = FALSE)
  seumi <- summarizeExperiment(list(s1 = res), coldata = data.frame(Name = "s1"), 
                               countType = "umis")
  secollumi <- collapseMutantsBySimilarity(seumi, assayName = "counts", 
                                           collapseMaxDist = 3, 
                                           collapseMinScore = 1, collapseMinRatio = 1.5, 
                                           verbose = FALSE)
  
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

  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose", "fastqForward", "fastqReverse"))) {
    if (nm == "variableCollapseMinReads") {
      expect_equal(res$parameters[[nm]], as.integer(ceiling(Ldef[[nm]])),
                   ignore_attr = TRUE)
    } else {
      expect_equal(res$parameters[[nm]], Ldef[[nm]], ignore_attr = TRUE)
    }
  }
  for (nm in c("fastqForward", "fastqReverse")) {
    expect_equal(res$parameters[[nm]], normalizePath(Ldef[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(nrow(res$summaryTable), 677L)
  expect_equal(sum(res$summaryTable$nbrUmis), 679L)
  expect_equal(nrow(secoll), 656L)
  expect_equal(nrow(secollumi), 656L)
})

test_that("digestFastqs works as expected for trans experiments, when similar sequences are collapsed (only high ratio), specify distance threshold as ratio", {
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
    umiCollapseMaxDist = 0,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = FALSE,
    nThreads = 1, chunkSize = 1000, 
    maxReadLength = 1024
  )
  
  res <- do.call(digestFastqs, Ldef)
  se <- summarizeExperiment(list(s1 = res), coldata = data.frame(Name = "s1"), 
                            countType = "reads")
  secoll <- collapseMutantsBySimilarity(se, assayName = "counts", 
                                        collapseMaxDist = 0.016, 
                                        collapseMinScore = 1, collapseMinRatio = 1.5, 
                                        verbose = FALSE)
  seumi <- summarizeExperiment(list(s1 = res), coldata = data.frame(Name = "s1"), 
                               countType = "umis")
  secollumi <- collapseMutantsBySimilarity(seumi, assayName = "counts", 
                                           collapseMaxDist = 0.016, 
                                           collapseMinScore = 1, collapseMinRatio = 1.5, 
                                           verbose = FALSE)
  
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
  
  for (nm in setdiff(names(Ldef), c("forbiddenMutatedCodonsForward", "forbiddenMutatedCodonsReverse", "verbose", "fastqForward", "fastqReverse"))) {
    if (nm == "variableCollapseMinReads") {
      expect_equal(res$parameters[[nm]], as.integer(ceiling(Ldef[[nm]])),
                   ignore_attr = TRUE)
    } else {
      expect_equal(res$parameters[[nm]], Ldef[[nm]], ignore_attr = TRUE)
    }
  }
  for (nm in c("fastqForward", "fastqReverse")) {
    expect_equal(res$parameters[[nm]], normalizePath(Ldef[[nm]], mustWork = FALSE), 
                 ignore_attr = TRUE)
  }
  
  expect_equal(sum(res$summaryTable$nbrReads), res$filterSummary$nbrRetained)
  expect_equal(nrow(res$summaryTable), 677L)
  expect_equal(sum(res$summaryTable$nbrUmis), 679L)
  expect_equal(nrow(secoll), 656L)
  expect_equal(nrow(secollumi), 656L)
})
