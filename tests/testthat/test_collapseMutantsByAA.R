## groupSimilarSequences
seqs <- c("AACGTAGCA", "ACCGTAGCA", "AACGGAGCA", "ATCGGAGCA", "TGAGGCATA")
scores <- c(5, 1, 3, 1, 8)
gr <- groupSimilarSequences(seqs = seqs, scores = scores, collapseMaxDist = 1, 
                            collapseMinScore = 0, collapseMinRatio = 0, verbose = FALSE)
expect_equal(gr$sequence, seqs)
expect_equal(gr$representative, seqs[c(1, 1, 1, 4, 5)])

gr <- groupSimilarSequences(seqs = seqs, scores = scores, collapseMaxDist = 0, 
                            collapseMinScore = 0, collapseMinRatio = 0, verbose = FALSE)
expect_equal(gr$sequence, seqs)
expect_equal(gr$representative, seqs[c(1, 2, 3, 4, 5)])

gr <- groupSimilarSequences(seqs = seqs, scores = scores, collapseMaxDist = 1, 
                            collapseMinScore = 0, collapseMinRatio = 2, verbose = FALSE)
expect_equal(gr$sequence, seqs)
expect_equal(gr$representative, seqs[c(1, 1, 3, 3, 5)])

gr <- groupSimilarSequences(seqs = seqs, scores = scores, collapseMaxDist = 1, 
                            collapseMinScore = 4, collapseMinRatio = 2, verbose = FALSE)
expect_equal(gr$sequence, seqs)
expect_equal(gr$representative, seqs[c(1, 1, 3, 4, 5)])

gr <- groupSimilarSequences(seqs = seqs, scores = scores, collapseMaxDist = 9, 
                            collapseMinScore = 4, collapseMinRatio = 2, verbose = FALSE)
expect_equal(gr$sequence, seqs)
expect_equal(gr$representative, seqs[c(1, 5, 5, 5, 5)])


## collapseMutants
Ldef <- list(
  mergeForwardReverse = TRUE, 
  minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
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
  umiCollapseMaxDist = 0,
  filteredReadsFastqForward = "",
  filteredReadsFastqReverse = "",
  maxNReads = -1, 
  verbose = FALSE,
  nThreads = 1, chunkSize = 1000, 
  maxReadLength = 1024
)
Ldef1 <- c(
  list(fastqForward = system.file("extdata/cisInput_1.fastq.gz", package = "mutscan"),
       fastqReverse = system.file("extdata/cisInput_2.fastq.gz", package = "mutscan")
  ), Ldef)
Ldef2 <- c(
  list(fastqForward = system.file("extdata/cisOutput_1.fastq.gz", package = "mutscan"),
       fastqReverse = system.file("extdata/cisOutput_2.fastq.gz", package = "mutscan")
  ), Ldef)

out1 <- do.call(digestFastqs, Ldef1)
out2 <- do.call(digestFastqs, Ldef2)

coldata <- data.frame(Name = c("sample1", "sample2"), 
                      Condition = c("input", "output"),
                      Replicate = c(1, 1), OD = c(0.05, 1.5),
                      stringsAsFactors = FALSE)

se <- summarizeExperiment(x = list(sample1 = out1, sample2 = out2),
                          coldata = coldata, countType = "umis")

test_that("collapseMutantsByAA fails with incorrect arguments", {
    expect_error(collapseMutantsByAA(se = 1))
    expect_error(collapseMutantsByAA(se = coldata))
    expect_error(collapseMutantsByAA(se = out1))
    expect_error(collapseMutantsByAA(se = list(sample1 = out1, sample2 = out2)))
    expect_error(collapseMutantsByAA(se = se, nameCol = 1))
    expect_error(collapseMutantsByAA(se = se, nameCol = c("sequence", "mutantNameAA")))
    expect_error(collapseMutantsByAA(se = se, nameCol = "missing"))
})

test_that("collapseMutantsBySimilarity fails with incorrect arguments", {
    expect_error(collapseMutantsBySimilarity(se = 1))
    expect_error(collapseMutantsBySimilarity(se = coldata))
    expect_error(collapseMutantsBySimilarity(se = out1))
    expect_error(collapseMutantsBySimilarity(se = list(sample1 = out1, sample2 = out2)))
    expect_error(collapseMutantsBySimilarity(se = se, assayName = 1),
                 "'assayName' must be of class 'character'")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "random"),
                 "All values in 'assayName'")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = c("counts", "counts")),
                 "'assayName' must have length 1")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "counts", 
                                             scoreMethod = 1),
                 "'scoreMethod' must be of class 'character'")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "counts", 
                                             scoreMethod = c("rowSum", "rowMean")),
                 "'scoreMethod' must have length 1")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "counts", 
                                             scoreMethod = "random"),
                 "All values in 'scoreMethod' must be one of")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "counts", 
                                             scoreMethod = "rowSum", sequenceCol = 1),
                 "'sequenceCol' must be of class 'character'")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "counts", 
                                             scoreMethod = "rowSum", 
                                             sequenceCol = c("sequence", "sequence")),
                 "'sequenceCol' must have length 1")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "counts", 
                                             scoreMethod = "rowSum", sequenceCol = "abc"),
                 "All values in 'sequenceCol' must be one of")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "counts", 
                                             scoreMethod = "rowSum", sequenceCol = "sequence",
                                             collapseMaxDist = "1"),
                 "'collapseMaxDist' must be of class 'numeric'")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "counts", 
                                             scoreMethod = "rowSum", sequenceCol = "sequence",
                                             collapseMaxDist = c(1, 2)),
                 "'collapseMaxDist' must have length 1")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "counts", 
                                             scoreMethod = "rowSum", sequenceCol = "sequence",
                                             collapseMaxDist = 1, collapseMinScore = "1"),
                 "'collapseMinScore' must be of class 'numeric'")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "counts", 
                                             scoreMethod = "rowSum", sequenceCol = "sequence",
                                             collapseMaxDist = 1, collapseMinScore = c(1, 2)),
                 "'collapseMinScore' must have length 1")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "counts", 
                                             scoreMethod = "rowSum", sequenceCol = "sequence",
                                             collapseMaxDist = 1, collapseMinScore = 1,
                                             collapseMinRatio = "1"),
                 "'collapseMinRatio' must be of class 'numeric'")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "counts", 
                                             scoreMethod = "rowSum", sequenceCol = "sequence",
                                             collapseMaxDist = 1, collapseMinScore = 1,
                                             collapseMinRatio = c(1, 2)),
                 "'collapseMinRatio' must have length 1")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "counts", 
                                             scoreMethod = "rowSum", sequenceCol = "sequence",
                                             collapseMaxDist = 1, collapseMinScore = 1,
                                             collapseMinRatio = 0, verbose = "TRUE"),
                 "'verbose' must be of class 'logical'")
    expect_error(collapseMutantsBySimilarity(se = se, assayName = "counts", 
                                             scoreMethod = "rowSum", sequenceCol = "sequence",
                                             collapseMaxDist = 1, collapseMinScore = 1,
                                             collapseMinRatio = 0, verbose = c(TRUE, FALSE)),
                 "'verbose' must have length 1")
    
    se0 <- se
    SummarizedExperiment::rowData(se0)$sequence[1] <- "ACGT"
    expect_error(collapseMutantsBySimilarity(se = se0, assayName = "counts", 
                                             scoreMethod = "rowSum", sequenceCol = "sequence",
                                             collapseMaxDist = 1, collapseMinScore = 1,
                                             collapseMinRatio = 0, verbose = TRUE),
                 "All entries in the sequence column must have the same length")
    se0 <- se
    SummarizedExperiment::rowData(se0)$sequence[1] <- gsub("A", "B", rowData(se0)$sequence[1])
    expect_error(collapseMutantsBySimilarity(se = se0, assayName = "counts", 
                                             scoreMethod = "rowSum", sequenceCol = "sequence",
                                             collapseMaxDist = 1, collapseMinScore = 1,
                                             collapseMinRatio = 0, verbose = TRUE),
                 "All entries in the sequence column must consist of")
})


test_that("collapseMutantsByAA works as expected", {
  secoll <- collapseMutantsByAA(se)
  expect_s4_class(secoll, "SummarizedExperiment")
  expect_equal(SummarizedExperiment::colData(se), 
               SummarizedExperiment::colData(secoll))
  expect_equal(S4Vectors::metadata(se), 
               S4Vectors::metadata(secoll))
  expect_equal(colnames(se), colnames(secoll))
  
  aapos <- unique(unlist(rowData(se)$mutantNameAA))
  expect_equal(nrow(secoll), length(aapos))
  expect_equal(sort(rownames(secoll)), sort(unique(aapos)))
  
  tmp <- table(rep(unlist(SummarizedExperiment::rowData(se)$mutantNameAA), 
                   SummarizedExperiment::assay(se, "counts")[, 1]))
  tmp <- tmp[rownames(secoll)]
  tmp[is.na(tmp)] <- 0
  expect_equal(SummarizedExperiment::assay(secoll, "counts")[, 1],
               as.numeric(tmp), ignore_attr = TRUE)
  
  tmp <- table(rep(unlist(SummarizedExperiment::rowData(se)$mutantNameAA), 
                   SummarizedExperiment::assay(se, "counts")[, 2]))
  tmp <- tmp[rownames(secoll)]
  tmp[is.na(tmp)] <- 0
  expect_equal(SummarizedExperiment::assay(secoll, "counts")[, 2],
               as.numeric(tmp), ignore_attr = TRUE)
  
  for (mn in rownames(secoll)) {
    expect_equal(SummarizedExperiment::assay(secoll, "counts")[mn, ],
                 colSums(SummarizedExperiment::assay(se, "counts")[SummarizedExperiment::rowData(se)$mutantNameAA == mn, , drop = FALSE]))
  }
  expect_equal(min(methods::as(
      lapply(strsplit(SummarizedExperiment::rowData(se)$nbrMutBases, ","), function(w) sort(as.integer(w))),
      "IntegerList")), 
      SummarizedExperiment::rowData(se)$minNbrMutBases)
  expect_equal(min(methods::as(
      lapply(strsplit(SummarizedExperiment::rowData(se)$nbrMutCodons, ","), function(w) sort(as.integer(w))),
      "IntegerList")), 
      SummarizedExperiment::rowData(se)$minNbrMutCodons)
  expect_equal(min(methods::as(
      lapply(strsplit(SummarizedExperiment::rowData(se)$nbrMutAAs, ","), function(w) sort(as.integer(w))),
      "IntegerList")), 
      SummarizedExperiment::rowData(se)$minNbrMutAAs)
  expect_equal(max(methods::as(
      lapply(strsplit(SummarizedExperiment::rowData(se)$nbrMutBases, ","), function(w) sort(as.integer(w))),
      "IntegerList")), 
      SummarizedExperiment::rowData(se)$maxNbrMutBases)
  expect_equal(max(methods::as(
      lapply(strsplit(SummarizedExperiment::rowData(se)$nbrMutCodons, ","), function(w) sort(as.integer(w))),
      "IntegerList")), 
      SummarizedExperiment::rowData(se)$maxNbrMutCodons)
  expect_equal(max(methods::as(
      lapply(strsplit(SummarizedExperiment::rowData(se)$nbrMutAAs, ","), function(w) sort(as.integer(w))),
      "IntegerList")), 
      SummarizedExperiment::rowData(se)$maxNbrMutAAs)
})

test_that("collapseMutantsByAA works as expected - collapseToWT", {
  Ldef1$collapseToWTForward <- TRUE
  Ldef2$collapseToWTForward <- TRUE
  out1 <- do.call(digestFastqs, Ldef1)
  out2 <- do.call(digestFastqs, Ldef2)
  
  se <- summarizeExperiment(x = list(sample1 = out1, sample2 = out2),
                            coldata = coldata, countType = "umis")
  
  secoll <- collapseMutantsByAA(se)
  expect_equal(SummarizedExperiment::colData(se), 
               SummarizedExperiment::colData(secoll))
  expect_equal(S4Vectors::metadata(se), 
               S4Vectors::metadata(secoll))
  expect_equal(colnames(se), colnames(secoll))
  
  aapos <- unique(unlist(rowData(se)$mutantNameAA))
  expect_equal(nrow(secoll), length(aapos))
  expect_equal(sort(rownames(secoll)), sort(unique(aapos)))
  
  tmp <- table(rep(unlist(SummarizedExperiment::rowData(se)$mutantNameAA), 
                   SummarizedExperiment::assay(se, "counts")[, 1]))
  tmp <- tmp[rownames(secoll)]
  tmp[is.na(tmp)] <- 0
  expect_equal(SummarizedExperiment::assay(secoll, "counts")[, 1],
               as.numeric(tmp), ignore_attr = TRUE)
  
  tmp <- table(rep(unlist(SummarizedExperiment::rowData(se)$mutantNameAA), 
                   SummarizedExperiment::assay(se, "counts")[, 2]))
  tmp <- tmp[rownames(secoll)]
  tmp[is.na(tmp)] <- 0
  expect_equal(SummarizedExperiment::assay(secoll, "counts")[, 2],
               as.numeric(tmp), ignore_attr = TRUE)
  
  for (mn in rownames(secoll)) {
    expect_equal(SummarizedExperiment::assay(secoll, "counts")[mn, ],
                 colSums(SummarizedExperiment::assay(se, "counts")[SummarizedExperiment::rowData(se)$mutantNameAA == mn, , drop = FALSE]))
  }
  expect_equal(min(methods::as(
      lapply(strsplit(SummarizedExperiment::rowData(se)$nbrMutBases, ","), function(w) sort(as.integer(w))),
      "IntegerList")), 
      SummarizedExperiment::rowData(se)$minNbrMutBases)
  expect_equal(min(methods::as(
      lapply(strsplit(SummarizedExperiment::rowData(se)$nbrMutCodons, ","), function(w) sort(as.integer(w))),
      "IntegerList")), 
      SummarizedExperiment::rowData(se)$minNbrMutCodons)
  expect_equal(min(methods::as(
      lapply(strsplit(SummarizedExperiment::rowData(se)$nbrMutAAs, ","), function(w) sort(as.integer(w))),
      "IntegerList")), 
      SummarizedExperiment::rowData(se)$minNbrMutAAs)
  expect_equal(max(methods::as(
      lapply(strsplit(SummarizedExperiment::rowData(se)$nbrMutBases, ","), function(w) sort(as.integer(w))),
      "IntegerList")), 
      SummarizedExperiment::rowData(se)$maxNbrMutBases)
  expect_equal(max(methods::as(
      lapply(strsplit(SummarizedExperiment::rowData(se)$nbrMutCodons, ","), function(w) sort(as.integer(w))),
      "IntegerList")), 
      SummarizedExperiment::rowData(se)$maxNbrMutCodons)
  expect_equal(max(methods::as(
      lapply(strsplit(SummarizedExperiment::rowData(se)$nbrMutAAs, ","), function(w) sort(as.integer(w))),
      "IntegerList")), 
      SummarizedExperiment::rowData(se)$maxNbrMutAAs)
  
  expect_equal(nrow(secoll), 1L)
  expect_equal(ncol(secoll), 2L)
  expect_equal(SummarizedExperiment::rowData(secoll)$nbrMutBases[[1]], "0,1,2")
  expect_equal(SummarizedExperiment::rowData(secoll)$nbrMutCodons[[1]], "0,1")
  expect_equal(SummarizedExperiment::rowData(secoll)$nbrMutAAs[[1]], "0,1")
  expect_equal(SummarizedExperiment::rowData(secoll)$mutationTypes, 
               "nonsynonymous,silent,stop")
})



