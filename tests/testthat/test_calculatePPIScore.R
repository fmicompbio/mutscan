context("PPI score calculations")

Ldef <- list(
  mergeForwardReverse = TRUE, 
  minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
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
  verbose = FALSE
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

coldata <- data.frame(Name = c("sample1", "sample2", "sample3", "sample4"), 
                      Condition = c("input", "output", "input", "output"),
                      Replicate = c(1, 1, 2, 2), OD = c(0.05, 1.5, 0.05, 1.5),
                      stringsAsFactors = FALSE)

se <- summarizeExperiment(x = list(sample1 = out1, sample2 = out2,
                                   sample3 = out2, sample4 = out1),
                          coldata = coldata, countType = "umis")
secoll <- collapseMutantsByAA(se)

test_that("calculatePPIScore fails with incorrect arguments", {
  expect_error(calculatePPIScore(se = 1, pairingCol = "Replicate", ODCols = "OD",
                                 comparison = c("Condition", "output", "input")))
  expect_error(calculatePPIScore(se = out1, pairingCol = "Replicate", ODCols = "OD",
                                 comparison = c("Condition", "output", "input")))
  expect_error(calculatePPIScore(se = list(sample1 = out1, sample2 = out2), 
                                 pairingCol = "Replicate", ODCols = "OD",
                                 comparison = c("Condition", "output", "input")))
  expect_error(calculatePPIScore(se = secoll, pairingCol = "Unknown", ODCols = "OD",
                                 comparison = c("Condition", "output", "input")))
  expect_error(calculatePPIScore(se = secoll, pairingCol = 1, ODCols = "OD",
                                 comparison = c("Condition", "output", "input")))
  expect_error(calculatePPIScore(se = secoll, pairingCol = "OD", ODCols = "OD",
                                 comparison = c("Condition", "output", "input")))
  expect_error(calculatePPIScore(se = secoll, pairingCol = "Replicate", ODCols = "Unknown",
                                 comparison = c("Condition", "output", "input")))
  expect_error(calculatePPIScore(se = secoll, pairingCol = "Replicate", ODCols = "Condition",
                                 comparison = c("Condition", "output", "input")))
  expect_error(calculatePPIScore(se = secoll, pairingCol = "Replicate", ODCols = 1,
                                 comparison = c("Condition", "output", "input")))
  expect_error(calculatePPIScore(se = secoll, pairingCol = "Replicate", ODCols = "OD",
                                 comparison = c("Unknown", "output", "input")))
  expect_error(calculatePPIScore(se = secoll, pairingCol = "Replicate", ODCols = "OD",
                                 comparison = c("Condition", "unknown", "input")))
  expect_error(calculatePPIScore(se = secoll, pairingCol = "Replicate", ODCols = "OD",
                                 comparison = c("Condition", "output", "unknown")))
})

test_that("calculatePPIScore works as expected", {
  ppi <- calculatePPIScore(se = secoll, pairingCol = "Replicate", ODCols = "OD",
                           comparison = c("Condition", "output", "input"), WTrows = "f.0.NA")
  
  expect_equal(nrow(ppi), nrow(secoll))
  expect_equal(rownames(ppi), rownames(secoll))
  expect_equal(ncol(ppi), 2)
  expect_equal(colnames(ppi), c("output_vs_input_repl1", "output_vs_input_repl2"))
  expect_equivalent(ppi["f.0.NA", ], c(1, 1))
  
  ## Test "replicate 1"
  w <- SummarizedExperiment::colData(secoll)$Condition == "output" & 
    SummarizedExperiment::colData(secoll)$Replicate == 1
  ncout <- as.matrix(SummarizedExperiment::assay(secoll, "counts"))[, w]
  ncout <- ncout * SummarizedExperiment::colData(secoll)$OD[w]/sum(ncout)
  
  w <- SummarizedExperiment::colData(secoll)$Condition == "input" & 
    SummarizedExperiment::colData(secoll)$Replicate == 1
  ncin <- as.matrix(SummarizedExperiment::assay(secoll, "counts"))[, w]
  ncin <- ncin * SummarizedExperiment::colData(secoll)$OD[w]/sum(ncin)
  
  ratios <- log2(ncout/ncin)
  ratios[!is.finite(ratios)] <- NA
  ratios <- ratios/ratios["f.0.NA"]
  
  expect_equal(ppi[, "output_vs_input_repl1"], ratios)
  
  ## Test "replicate 2"
  w <- SummarizedExperiment::colData(secoll)$Condition == "output" & 
    SummarizedExperiment::colData(secoll)$Replicate == 2
  ncout <- as.matrix(SummarizedExperiment::assay(secoll, "counts"))[, w]
  ncout <- ncout * SummarizedExperiment::colData(secoll)$OD[w]/sum(ncout)
  
  w <- SummarizedExperiment::colData(secoll)$Condition == "input" & 
    SummarizedExperiment::colData(secoll)$Replicate == 2
  ncin <- as.matrix(SummarizedExperiment::assay(secoll, "counts"))[, w]
  ncin <- ncin * SummarizedExperiment::colData(secoll)$OD[w]/sum(ncin)
  
  ratios <- log2(ncout/ncin)
  ratios[!is.finite(ratios)] <- NA
  ratios <- ratios/ratios["f.0.NA"]
  
  expect_equal(ppi[, "output_vs_input_repl2"], ratios)
})




