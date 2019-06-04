context("collapse mutants by amino acid")

Ldef <- list(
  mergeForwardReverse = TRUE, revComplForward = FALSE, revComplReverse = TRUE,
  skipForward = 1, skipReverse = 1, 
  umiLengthForward = 10, umiLengthReverse = 7, 
  constantLengthForward = 18, constantLengthReverse = 17, 
  variableLengthForward = 96, variableLengthReverse = 96,
  adapterForward = "GGAAGAGCACACGTC", 
  adapterReverse = "GGAAGAGCGTCGTGT",
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
})

test_that("collapseMutantsByAA works as expected", {
  secoll <- collapseMutantsByAA(se)
  expect_equal(SummarizedExperiment::colData(se), 
               SummarizedExperiment::colData(secoll))
  expect_equal(S4Vectors::metadata(se), 
               S4Vectors::metadata(secoll))
  expect_equal(colnames(se), colnames(secoll))
  
  ## Number of unique amino acids/position combinations
  aapos <- paste0(substr(rownames(se), start = 1, stop = nchar(rownames(se)) - 3), 
                  GENETIC_CODE[substr(rownames(se), start = nchar(rownames(se)) - 2, 
                                      stop = nchar(rownames(se)))])
  aapos[aapos == "NA"] <- "WT"
  expect_equal(nrow(secoll), length(unique(aapos)))
  expect_equal(sort(rownames(secoll)), sort(unique(aapos)))
  
  tmp <- table(rep(aapos, SummarizedExperiment::assay(se, "counts")[, 1]))
  tmp <- tmp[rownames(secoll)]
  tmp[is.na(tmp)] <- 0
  expect_equivalent(SummarizedExperiment::assay(secoll, "counts")[, 1],
                    as.numeric(tmp))
  
  tmp <- table(rep(aapos, SummarizedExperiment::assay(se, "counts")[, 2]))
  tmp <- tmp[rownames(secoll)]
  tmp[is.na(tmp)] <- 0
  expect_equivalent(SummarizedExperiment::assay(secoll, "counts")[, 2],
                    as.numeric(tmp))
})


