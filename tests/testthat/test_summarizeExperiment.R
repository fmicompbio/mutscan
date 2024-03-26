Ldef <- list(
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
    collapseToWTForward = FALSE,
    collapseToWTReverse = FALSE, 
    mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".",
    constantMaxDistForward = -1,
    constantMaxDistReverse = -1,
    umiCollapseMaxDist = 0,
    maxNReads = -1, verbose = FALSE
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

    ## samples must have the same mutNameDelimiter
    tmpout2 <- out2
    tmpout2$parameters$mutNameDelimiter <- ":"
    expect_error(summarizeExperiment(x = list(sample1 = out1, sample2 = tmpout2),
                                     coldata = coldata, countType = "reads"))
})

test_that("summarizeExperiment works as expected with reads output", {
    se <- summarizeExperiment(x = list(sample1 = out1, sample2 = out2),
                              coldata = coldata, countType = "reads")
    expect_equal(nrow(se), length(union(out1$summaryTable$mutantName,
                                        out2$summaryTable$mutantName)))
    expect_equal(ncol(se), 2)
    expect_equal(sort(rownames(se)),
                 sort(union(out1$summaryTable$mutantName, out2$summaryTable$mutantName)))
    expect_equal(Matrix::colSums(SummarizedExperiment::assay(se, "counts")),
                 c(out1$filterSummary$nbrRetained, out2$filterSummary$nbrRetained),
                 ignore_attr = TRUE)
    expect_equal(Matrix::colSums(SummarizedExperiment::assay(se, "counts")),
                 c(sum(out1$summaryTable$nbrReads), sum(out2$summaryTable$nbrReads)),
                 ignore_attr = TRUE)

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

    expect_equal(S4Vectors::metadata(se)$parameters[["sample1"]], out1$parameters)
    expect_equal(S4Vectors::metadata(se)$parameters[["sample2"]], out2$parameters)
    expect_equal(S4Vectors::metadata(se)$errorStatistics[["sample1"]], out1$errorStatistics)
    expect_equal(S4Vectors::metadata(se)$errorStatistics[["sample2"]], out2$errorStatistics)

    expect_true(all(c("mutantName", "sequence", "nbrMutBases",
                      "minNbrMutBases", "maxNbrMutBases",
                      "nbrMutCodons", "minNbrMutCodons", "maxNbrMutCodons",
                      "nbrMutAAs", "minNbrMutAAs", "maxNbrMutAAs",
                      "sequenceAA", "mutantNameAA", "mutationTypes",
                      "varLengths") %in% 
                        colnames(SummarizedExperiment::rowData(se))))
    
    ## Check that the number of mutated codons are correct (=equal to the number of
    ## entries in the mutant name that don't contain WT)
    expect_equal(sapply(strsplit(SummarizedExperiment::rowData(se)$mutantName, "_"),
                        function(w) length(w[!grepl("WT", w)])),
                 SummarizedExperiment::rowData(se)$minNbrMutCodons, ignore_attr = TRUE)
    expect_equal(sapply(strsplit(SummarizedExperiment::rowData(se)$mutantName, "_"),
                        function(w) length(w[!grepl("WT", w)])),
                 SummarizedExperiment::rowData(se)$maxNbrMutCodons, ignore_attr = TRUE)
    
    ## Check that the number of mutated bases is not larger than 3x the number of mutated codons
    expect_true(all(SummarizedExperiment::rowData(se)$maxNbrMutBases <=  
                        3 * SummarizedExperiment::rowData(se)$maxNbrMutCodons))
    expect_true(all(SummarizedExperiment::rowData(se)$maxNbrMutAAs <= 
                        SummarizedExperiment::rowData(se)$maxNbrMutCodons))
    
    ## variable lengths
    expect_equal(SummarizedExperiment::rowData(se)$varLengths, 
                 rep("96_96", nrow(se)), ignore_attr = TRUE)
    
    ## All variants with no mutated AAs must have a WT in the name
    expect_true(all(grepl("WT", SummarizedExperiment::rowData(se)$mutantNameAA[SummarizedExperiment::rowData(se)$maxNbrMutAAs == 0])))
    ## No mutation for the complete WT
    expect_equal(SummarizedExperiment::rowData(se)["f.0.WT_r.0.WT", ]$mutationTypes, "")
    ## Mutation types for variants with no mutated AAs, but mutated bases, should be silent
    expect_true(all(SummarizedExperiment::rowData(se)$mutationTypes[SummarizedExperiment::rowData(se)$maxNbrMutAAs == 0 & SummarizedExperiment::rowData(se)$maxNbrMutBases > 0] == "silent"))
    ## check translation
    expect_equal(SummarizedExperiment::rowData(se)$sequenceAA[3], 
                 mutscan:::translateString(SummarizedExperiment::rowData(se)$sequence[3]))
    
    ## Spot checks
    expect_equal(SummarizedExperiment::rowData(se)$minNbrMutBases[SummarizedExperiment::rowData(se)$mutantName == "f.0.WT_r.13.CTC"], 3)  ## WT: GCT
    expect_equal(SummarizedExperiment::rowData(se)$minNbrMutCodons[SummarizedExperiment::rowData(se)$mutantName == "f.0.WT_r.13.CTC"], 1)  ## WT: GCT
    expect_equal(SummarizedExperiment::rowData(se)$minNbrMutAAs[SummarizedExperiment::rowData(se)$mutantName == "f.0.WT_r.13.CTC"], 1)  ## WT: A
    expect_equal(SummarizedExperiment::rowData(se)$mutantNameAA[SummarizedExperiment::rowData(se)$mutantName == "f.0.WT_r.13.CTC"], "f.0.WT_r.13.L")  ## WT: A
    expect_equal(SummarizedExperiment::rowData(se)$mutationTypes[SummarizedExperiment::rowData(se)$mutantName == "f.0.WT_r.13.CTC"], "nonsynonymous")  ## WT: A
    
    expect_equal(SummarizedExperiment::rowData(se)$minNbrMutBases[SummarizedExperiment::rowData(se)$mutantName == "f.0.WT_r.13.GCG"], 1)  ## WT: GCT
    expect_equal(SummarizedExperiment::rowData(se)$minNbrMutCodons[SummarizedExperiment::rowData(se)$mutantName == "f.0.WT_r.13.GCG"], 1)  ## WT: GCT
    expect_equal(SummarizedExperiment::rowData(se)$minNbrMutAAs[SummarizedExperiment::rowData(se)$mutantName == "f.0.WT_r.13.GCG"], 0)  ## WT: A
    expect_equal(SummarizedExperiment::rowData(se)$mutantNameAA[SummarizedExperiment::rowData(se)$mutantName == "f.0.WT_r.13.GCG"], "f.0.WT_r.0.WT")  ## WT: A
    expect_equal(SummarizedExperiment::rowData(se)$mutationTypes[SummarizedExperiment::rowData(se)$mutantName == "f.0.WT_r.13.GCG"], "silent")  ## WT: A
    
})

test_that("summarizeExperiment works as expected with umis output", {
    se <- summarizeExperiment(x = list(sample1 = out1, sample2 = out2),
                              coldata = coldata, countType = "umis")
    expect_equal(nrow(se), length(union(out1$summaryTable$mutantName,
                                        out2$summaryTable$mutantName)))
    expect_equal(ncol(se), 2)
    expect_equal(sort(rownames(se)),
                 sort(union(out1$summaryTable$mutantName, out2$summaryTable$mutantName)))
    expect_equal(Matrix::colSums(SummarizedExperiment::assay(se, "counts")),
                 c(sum(out1$summaryTable$nbrUmis), sum(out2$summaryTable$nbrUmis)),
                 ignore_attr = TRUE)

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

    expect_equal(S4Vectors::metadata(se)$parameters[["sample1"]], out1$parameters)
    expect_equal(S4Vectors::metadata(se)$parameters[["sample2"]], out2$parameters)
    expect_equal(S4Vectors::metadata(se)$errorStatistics[["sample1"]], out1$errorStatistics)
    expect_equal(S4Vectors::metadata(se)$errorStatistics[["sample2"]], out2$errorStatistics)
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

test_that("summarizeExperiment recognizes the presence of UMI counts correctly", {
    L1 <- Ldef1; L1$elementsForward <- "SSCV"; L1$elementsReverse <- "SSCV"
    L2 <- Ldef2; L2$elementsForward <- "SSCV"; L2$elementsReverse <- "SSCV"
    outl1 <- do.call(digestFastqs, L1)
    outl2 <- do.call(digestFastqs, L2)
    expect_error(summarizeExperiment(x = list(sample1 = outl1, sample2 = outl2),
                                     coldata = coldata, countType = "umis"))

})

test_that("summarizeExperiment works as expected when collapsing to WT", {
    ## Also test that errorStatistics is empty if there is no constant sequence
    Ldef1$collapseToWTForward <- TRUE
    Ldef2$collapseToWTForward <- TRUE
    Ldef1$elementsForward = "SUSV"
    Ldef1$elementsReverse = "SUSV"
    Ldef2$elementsForward = "SUSV"
    Ldef2$elementsReverse = "SUSV"
    
    out1 <- do.call(digestFastqs, Ldef1)
    out2 <- do.call(digestFastqs, Ldef2)
    
    se <- summarizeExperiment(x = list(sample1 = out1, sample2 = out2),
                              coldata = coldata, countType = "reads")
    expect_equal(nrow(se), length(union(out1$summaryTable$mutantName,
                                        out2$summaryTable$mutantName)))
    expect_equal(ncol(se), 2)
    expect_equal(sort(rownames(se)),
                 sort(union(out1$summaryTable$mutantName, out2$summaryTable$mutantName)))
    expect_equal(Matrix::colSums(SummarizedExperiment::assay(se, "counts")),
                 c(out1$filterSummary$nbrRetained, out2$filterSummary$nbrRetained),
                 ignore_attr = TRUE)
    expect_equal(Matrix::colSums(SummarizedExperiment::assay(se, "counts")),
                 c(sum(out1$summaryTable$nbrReads), sum(out2$summaryTable$nbrReads)),
                 ignore_attr = TRUE)
    
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
    
    expect_equal(S4Vectors::metadata(se)$parameters[["sample1"]], out1$parameters)
    expect_equal(S4Vectors::metadata(se)$parameters[["sample2"]], out2$parameters)
    expect_equal(S4Vectors::metadata(se)$errorStatistics[["sample1"]], out1$errorStatistics)
    expect_equal(S4Vectors::metadata(se)$errorStatistics[["sample2"]], out2$errorStatistics)
    expect_equal(sum(S4Vectors::metadata(se)$errorStatistics[["sample1"]][, c("nbrMatchForward", "nbrMismatchForward", "nbrMatchReverse", "nbrMismatchReverse")]), 0)
    expect_equal(sum(S4Vectors::metadata(se)$errorStatistics[["sample2"]][, c("nbrMatchForward", "nbrMismatchForward", "nbrMatchReverse", "nbrMismatchReverse")]), 0)
    
    expect_true(all(c("mutantName", "sequence", "nbrMutBases",
                      "minNbrMutBases", "maxNbrMutBases",
                      "nbrMutCodons", "minNbrMutCodons", "maxNbrMutCodons",
                      "nbrMutAAs", "minNbrMutAAs", "maxNbrMutAAs",
                      "sequenceAA", "mutantNameAA", "mutationTypes",
                      "varLengths") %in% 
                        colnames(SummarizedExperiment::rowData(se))))
    
    ## variable lengths
    expect_equal(SummarizedExperiment::rowData(se)$varLengths, 
                 rep("96_96", nrow(se)), ignore_attr = TRUE)
    
    expect_false(any(grepl("^,", SummarizedExperiment::rowData(se)$mutationTypes)))
    expect_false(any(grepl(",$", SummarizedExperiment::rowData(se)$mutationTypes)))
    expect_equal(SummarizedExperiment::rowData(se)$sequenceAA[3], 
                 mutscan:::translateString(SummarizedExperiment::rowData(se)$sequence[3]))
    expect_type(SummarizedExperiment::rowData(se)$nbrMutBases, "character")
    expect_type(SummarizedExperiment::rowData(se)$nbrMutCodons, "character")
    expect_type(SummarizedExperiment::rowData(se)$nbrMutAAs, "character")
    expect_equal(SummarizedExperiment::rowData(se)$nbrMutBases[[which(SummarizedExperiment::rowData(se)$mutantName == "f_r.0.WT")]], "0,1,2,3")
    expect_equal(SummarizedExperiment::rowData(se)$nbrMutCodons[[which(SummarizedExperiment::rowData(se)$mutantName == "f_r.0.WT")]], "0,1")
    expect_equal(SummarizedExperiment::rowData(se)$nbrMutAAs[[which(SummarizedExperiment::rowData(se)$mutantName == "f_r.0.WT")]], "0,1")
    
    expect_true(all(grepl("stop", SummarizedExperiment::rowData(se)$mutationTypes[grep("\\*", SummarizedExperiment::rowData(se)$mutantNameAA)])))
})

test_that("mergeValues works", {
    res <- mergeValues(c("A", "B", "C", "A", "D", "B"),
                       c("a,b", "b,c", "c", "b,c", "b,a", "d"))
    expect_s3_class(res, "data.frame")
    expect_named(res, c("mutantNameColl", "valueColl"))
    expect_equal(res$mutantNameColl, c("A", "B", "C", "D"))
    expect_equal(res$value, c("a,b,c", "b,c,d", "c", "a,b"))
    
    res <- mergeValues(c("B", "A", "C", "D"),
                       c("a,b", "c", "c", "d"))
    expect_s3_class(res, "data.frame")
    expect_named(res, c("mutantNameColl", "valueColl"))
    expect_equal(res$mutantNameColl, c("A", "B", "C", "D"))
    expect_equal(res$value, c("c", "a,b", "c", "d"))
    
    res <- mergeValues(c("A", "B", "C", "A", "D", "B"),
                       c("a:b", "b:c", "c", "b:c", "b:a", "d"),
                       delimiter = ":")
    expect_s3_class(res, "data.frame")
    expect_named(res, c("mutantNameColl", "valueColl"))
    expect_equal(res$mutantNameColl, c("A", "B", "C", "D"))
    expect_equal(res$value, c("a:b:c", "b:c:d", "c", "a:b"))
    
    res <- mergeValues(c("A", "B", "C", "A", "D", "B"),
                       c("a,b", "b:c", "c", "b:c", "b:a", "d"),
                       delimiter = ":")
    expect_s3_class(res, "data.frame")
    expect_named(res, c("mutantNameColl", "valueColl"))
    expect_equal(res$mutantNameColl, c("A", "B", "C", "D"))
    expect_equal(res$value, c("a,b:b:c", "b:c:d", "c", "a:b"))
})




