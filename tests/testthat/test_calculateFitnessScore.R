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

coldata <- data.frame(Name = c("sample1", "sample2", "sample3", "sample4"), 
                      Condition = c("input", "output", "input", "output"),
                      Replicate = c(1, 1, 2, 2), OD = c(0.05, 1.5, 0.05, 1.5),
                      stringsAsFactors = FALSE)

se <- summarizeExperiment(x = list(sample1 = out1, sample2 = out2,
                                   sample3 = out2, sample4 = out1),
                          coldata = coldata, countType = "umis")
secoll <- collapseMutantsByAA(se)

test_that("calculateFitnessScore fails with incorrect arguments", {
    expect_error(calculateFitnessScore(se = 1, pairingCol = "Replicate", ODCols = "OD",
                                       comparison = c("Condition", "output", "input"),
                                       WTrows = "f.0.NA"))
    expect_error(calculateFitnessScore(se = out1, pairingCol = "Replicate", ODCols = "OD",
                                       comparison = c("Condition", "output", "input"),
                                       WTrows = "f.0.NA"))
    expect_error(calculateFitnessScore(se = list(sample1 = out1, sample2 = out2), 
                                       pairingCol = "Replicate", ODCols = "OD",
                                       comparison = c("Condition", "output", "input"),
                                       WTrows = "f.0.NA"))
    expect_error(calculateFitnessScore(se = secoll, pairingCol = "Unknown", ODCols = "OD",
                                       comparison = c("Condition", "output", "input"),
                                       WTrows = "f.0.NA"))
    expect_error(calculateFitnessScore(se = secoll, pairingCol = 1, ODCols = "OD",
                                       comparison = c("Condition", "output", "input"),
                                       WTrows = "f.0.NA"))
    expect_error(calculateFitnessScore(se = secoll, pairingCol = "OD", ODCols = "OD",
                                       comparison = c("Condition", "output", "input"),
                                       WTrows = "f.0.NA"))
    expect_error(calculateFitnessScore(se = secoll, pairingCol = "Replicate", ODCols = "Unknown",
                                       comparison = c("Condition", "output", "input"),
                                       WTrows = "f.0.NA"))
    expect_error(calculateFitnessScore(se = secoll, pairingCol = "Replicate", ODCols = "Condition",
                                       comparison = c("Condition", "output", "input"),
                                       WTrows = "f.0.NA"))
    expect_error(calculateFitnessScore(se = secoll, pairingCol = "Replicate", ODCols = 1,
                                       comparison = c("Condition", "output", "input"),
                                       WTrows = "f.0.NA"))
    expect_error(calculateFitnessScore(se = secoll, pairingCol = "Replicate", ODCols = "OD",
                                       comparison = c("Unknown", "output", "input"),
                                       WTrows = "f.0.NA"))
    expect_error(calculateFitnessScore(se = secoll, pairingCol = "Replicate", ODCols = "OD",
                                       comparison = c("Condition", "unknown", "input"),
                                       WTrows = "f.0.NA"))
    expect_error(calculateFitnessScore(se = secoll, pairingCol = "Replicate", ODCols = "OD",
                                       comparison = c("Condition", "output", "unknown"),
                                       WTrows = "f.0.NA"))
    expect_error(calculateFitnessScore(se = secoll, pairingCol = "Replicate", ODCols = "OD",
                                       comparison = c("Condition", "output", "input"),
                                       WTrows = "f.0.NA", selAssay = "missing"))
})

test_that("calculateFitnessScore works as expected", {
    ppi <- calculateFitnessScore(se = secoll, pairingCol = "Replicate", ODCols = "OD",
                                 comparison = c("Condition", "output", "input"), 
                                 WTrows = "f.0.WT")
    
    expect_equal(nrow(ppi), nrow(secoll))
    expect_equal(rownames(ppi), rownames(secoll))
    expect_equal(ncol(ppi), 2)
    expect_equal(colnames(ppi), c("output_vs_input_repl1", "output_vs_input_repl2"))
    expect_equal(ppi["f.0.WT", ], c(1, 1), ignore_attr = TRUE)
    
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
    ratios <- ratios/ratios["f.0.WT"]
    
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
    
    ## Test performance if WTrows = NULL
    ppinull <- calculateFitnessScore(se = secoll, pairingCol = "Replicate",
                                     ODCols = "OD",
                                     comparison = c("Condition", "output", "input"),
                                     WTrows = NULL)
    expect_equal(ppinull[, "output_vs_input_repl2"] , ratios)
    
    ## Use the WTrows
    ratios <- ratios/ratios["f.0.WT"]
    expect_equal(ppi[, "output_vs_input_repl2"], ratios)
    
    ## Test that PPI scores of input/output gives the same as output/input
    ppirev <- calculateFitnessScore(se = secoll, pairingCol = "Replicate", 
                                    ODCols = "OD",
                                    comparison = c("Condition", "input", "output"), 
                                    WTrows = "f.0.WT")
    expect_equal(ppirev[, "input_vs_output_repl1"], ppi[, "output_vs_input_repl1"])
    
})

