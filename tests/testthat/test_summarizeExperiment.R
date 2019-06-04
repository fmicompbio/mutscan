context("summarizeExperiment")

Ldef <- list(
  mergeForwardReverse = FALSE, revComplForward = FALSE, revComplReverse = FALSE,
  skipForward = 1, skipReverse = 1, 
  umiLengthForward = 10, umiLengthReverse = 8, 
  constantLengthForward = 18, constantLengthReverse = 20, 
  variableLengthForward = 96, variableLengthReverse = 96,
  adapterForward = "GGAAGAGCACACGTC", 
  adapterReverse = "GGAAGAGCGTCGTGT",
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
  verbose = FALSE
)
Ldef1 <- c(
  list(fastqForward = system.file("extdata/transInput_1.fastq.gz", package = "mutscan"),
       fastqReverse = system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  ), Ldef)
Ldef2 <- c(
  list(fastqForward = system.file("extdata/transOutput_1.fastq.gz", package = "mutscan"),
       fastqReverse = system.file("extdata/transOutput_2.fastq.gz", package = "mutscan")
  ), Ldef)

out1 <- do.call(digestFastqs, Ldef1)
out2 <- do.call(digestFastqs, Ldef2)

coldata <- data.frame(Name = c("sample1", "sample2"), 
                      Condition = c("input", "output"),
                      Replicate = c(1, 1), OD = c(0.05, 1.5),
                      stringsAsFactors = FALSE)

test_that("summarizeExperiment fails with incorrect arguments", {
  ## x must be a named list, with names matching coldata$Name
  expect_error(summarizeExperiment(x = 1, coldata = coldata))
  expect_error(summarizeExperiment(x = list(out1), coldata = coldata))
  expect_error(summarizeExperiment(x = out1, coldata = coldata))
  expect_error(summarizeExperiment(x = list(s1 = out1, s1 = out2), coldata = coldata))
  
  ## coldata must be a data.frame with a column Name
  expect_error(summarizeExperiment(x = list(sample1 = out1, sample2 = out2), 
                                   coldata = 1))
  expect_error(summarizeExperiment(x = list(sample1 = out1, sample2 = out2), 
                                   coldata = coldata[, c("Condition", "OD")]))
  expect_error(summarizeExperiment(x = list(sample1 = out1, sample2 = out2), 
                                   coldata = as.list(coldata)))
  
  ## countType must be reads or umis
  expect_error(summarizeExperiment(x = list(sample1 = out1, sample2 = out2),
                                   coldata = coldata, countType = 1))
  expect_error(summarizeExperiment(x = list(sample1 = out1, sample2 = out2),
                                   coldata = coldata, countType = "umi"))
  expect_error(summarizeExperiment(x = list(sample1 = out1, sample2 = out2),
                                   coldata = coldata, countType = c("reads", "umis")))
  
  ## names must match
  expect_error(summarizeExperiment(x = list(s1 = out1, s2 = out2),
                                   coldata = coldata, countType = "umis"))
  expect_warning(summarizeExperiment(x = list(sample1 = out1, sample2 = out2,
                                              sample3 = out1),
                                     coldata = coldata, countType = "umis"))
})

test_that("summarizeExperiment works as expected with reads output", {
  se <- summarizeExperiment(x = list(sample1 = out1, sample2 = out2),
                            coldata = coldata, countType = "reads")
  expect_equal(nrow(se), length(union(out1$summaryTable$mutantName,
                                      out2$summaryTable$mutantName)))
  expect_equal(ncol(se), 2)
  expect_equal(sort(rownames(se)), 
               sort(union(out1$summaryTable$mutantName, out2$summaryTable$mutantName)))
  expect_equivalent(Matrix::colSums(SummarizedExperiment::assay(se, "counts")),
                    c(out1$filterSummary$nbrRetained, out2$filterSummary$nbrRetained))
  expect_equivalent(Matrix::colSums(SummarizedExperiment::assay(se, "counts")),
                    c(sum(out1$summaryTable$nbrReads), sum(out2$summaryTable$nbrReads)))

  for (cn in colnames(out1$filterSummary)) {
    expect_equal(SummarizedExperiment::colData(se)["sample1", cn],
                 out1$filterSummary[, cn])
    expect_equal(SummarizedExperiment::colData(se)["sample2", cn],
                 out2$filterSummary[, cn])
  }
  
  for (cn in colnames(coldata)) {
    expect_equal(SummarizedExperiment::colData(se)["sample1", cn],
                 coldata[1, cn])
    expect_equal(SummarizedExperiment::colData(se)["sample2", cn],
                 coldata[2, cn])
  }
  
  expect_equal(S4Vectors::metadata(se)[["sample1"]], out1$parameters)
  expect_equal(S4Vectors::metadata(se)[["sample2"]], out2$parameters)
})

test_that("summarizeExperiment works as expected with umis output", {
  se <- summarizeExperiment(x = list(sample1 = out1, sample2 = out2),
                            coldata = coldata, countType = "umis")
  expect_equal(nrow(se), length(union(out1$summaryTable$mutantName,
                                      out2$summaryTable$mutantName)))
  expect_equal(ncol(se), 2)
  expect_equal(sort(rownames(se)), 
               sort(union(out1$summaryTable$mutantName, out2$summaryTable$mutantName)))
  expect_equivalent(Matrix::colSums(SummarizedExperiment::assay(se, "counts")),
                    c(sum(out1$summaryTable$nbrUmis), sum(out2$summaryTable$nbrUmis)))
  
  for (cn in colnames(out1$filterSummary)) {
    expect_equal(SummarizedExperiment::colData(se)["sample1", cn],
                 out1$filterSummary[, cn])
    expect_equal(SummarizedExperiment::colData(se)["sample2", cn],
                 out2$filterSummary[, cn])
  }
  
  for (cn in colnames(coldata)) {
    expect_equal(SummarizedExperiment::colData(se)["sample1", cn],
                 coldata[1, cn])
    expect_equal(SummarizedExperiment::colData(se)["sample2", cn],
                 coldata[2, cn])
  }
  
  expect_equal(S4Vectors::metadata(se)[["sample1"]], out1$parameters)
  expect_equal(S4Vectors::metadata(se)[["sample2"]], out2$parameters)
})

test_that("summarizeExperiment orders samples equally in count matrix/colData", {
  se1 <- summarizeExperiment(x = list(sample1 = out1, sample2 = out2),
                             coldata = coldata, countType = "umis")
  se2 <- summarizeExperiment(x = list(sample2 = out2, sample1 = out1),
                             coldata = coldata, countType = "umis")
  m1 <- se1[, match(coldata$Name, colnames(se1))]
  m2 <- se2[rownames(m1), match(coldata$Name, colnames(se2))]
  expect_equal(SummarizedExperiment::colData(m1), SummarizedExperiment::colData(m2))
  expect_equal(SummarizedExperiment::assay(m1, "counts"),
               SummarizedExperiment::assay(m2, "counts"))
})









