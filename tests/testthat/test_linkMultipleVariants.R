test_that("linkMultipleVariants works", {
    ## Fails with the wrong input
    ## ------------------------------------------------------------------------
    expect_error(linkMultipleVariants(combinedDigestParams = 1),
                 "'combinedDigestParams' must be of class 'list'")
    expect_error(linkMultipleVariants(combinedDigestParams = list(1)),
                 "'namescombinedDigestParams' must not be NULL")
    expect_error(linkMultipleVariants(combinedDigestParams = list(unknown = 1)),
                 "All values in 'namescombinedDigestParams' must be one of")
    expect_error(linkMultipleVariants(combinedDigestParams = list(fastqForward = 1),
                                      secondArg = 1),
                 "'parms' must be of class 'list'")
    expect_error(linkMultipleVariants(combinedDigestParams = list(fastqForward = 1),
                                      secondArg = list(1)),
                 "'namesparms' must not be NULL")
    expect_error(linkMultipleVariants(combinedDigestParams = list(fastqForward = 1),
                                      secondArg = list(unknown = 1)),
                 "All values in 'namesparms' must be one of")
    
    expect_error(linkMultipleVariants(combinedDigestParams = list(fastqForward = 1),
                                      secondArg = list(elementsForward = "CC")),
                 "Each separate run must contain at least one variable")
    expect_error(linkMultipleVariants(combinedDigestParams = list(fastqForward = 1),
                                      secondArg = list(elementsForward = "V",
                                                       elementsReverse = "V")),
                 "A separate run can not have variable sequences in both")
    expect_error(linkMultipleVariants(combinedDigestParams = list(fastqForward = 1),
                                      secondArg = list(elementsForward = "VV")),
                 "A separate run can have at most one variable segment")
    
    ## Works with correct input
    ## ------------------------------------------------------------------------
    ## Read truth
    truth <- readRDS(system.file("extdata", "multipleVariableRegions_truth.rds",
                                 package = "mutscan"))
    verbose <- FALSE
    
    ## ------------------------------------------------------------------------
    ## Don't filter by mutations in constant sequence
    res <- linkMultipleVariants(
        combinedDigestParams = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            fastqReverse = system.file("extdata", "multipleVariableRegions_R2.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CVCVCV", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            elementsReverse = "CVCV", elementLengthsReverse = c(6, 40, 10, -1),
            mergeForwardReverse = TRUE, minOverlap = 20, 
            maxOverlap = 30, maxFracMismatchOverlap = 0,
            revComplForward = FALSE, revComplReverse = TRUE,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            verbose = verbose),
        barcode = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CVCSCS", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            avePhredMinForward = 20, variableCollapseMaxDist = 1, 
            variableCollapseMinReads = 1, variableCollapseMinRatio = 1,
            verbose = verbose),
        V2 = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            fastqReverse = system.file("extdata", "multipleVariableRegions_R2.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CSCVCS", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            elementsReverse = "CSCV", elementLengthsReverse = c(6, 40, 10, -1),
            mergeForwardReverse = TRUE, minOverlap = 10, 
            maxOverlap = 20, maxFracMismatchOverlap = 0,
            revComplForward = FALSE, revComplReverse = TRUE,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            wildTypeForward = truth$WTV2,
            nbrMutatedCodonsMaxForward = -1, nbrMutatedBasesMaxForward = 1, 
            forbiddenMutatedCodonsForward = "", collapseToWTForward = FALSE,
            verbose = verbose, useTreeWTmatch = TRUE),
        V3 = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            fastqReverse = system.file("extdata", "multipleVariableRegions_R2.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CSCSCV", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            elementsReverse = "CVCS", elementLengthsReverse = c(6, 40, 10, -1),
            mergeForwardReverse = TRUE, minOverlap = 5, 
            maxOverlap = 15, maxFracMismatchOverlap = 0,
            revComplForward = FALSE, revComplReverse = TRUE,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            wildTypeForward = truth$WTV3,
            nbrMutatedCodonsMaxForward = -1, nbrMutatedBasesMaxForward = 1, 
            forbiddenMutatedCodonsForward = "", collapseToWTForward = FALSE,
            verbose = verbose, useTreeWTmatch = TRUE)
    )
    
    ## Filter out only reads with N/deletions in the variable sequence
    keep <- !grepl("N|del", truth$truth$status)
    
    ## Truth
    correct <- truth$truth[keep, ] %>%
        dplyr::group_by(trueBarcode, trueV2, trueV3) %>% dplyr::tally() %>%
        dplyr::arrange(dplyr::desc(n), trueBarcode, trueV2, trueV3)
    
    ## Obs
    obs <- res$countAggregated %>% dplyr::arrange(desc(nbrReads), barcode, V2, V3)
    
    expect_equal(sum(obs$nbrReads), sum(keep))
    expect_equal(correct$trueBarcode, obs$barcode)
    expect_equal(correct$trueV2, obs$V2)
    expect_equal(correct$trueV3, obs$V3)
    expect_equal(correct$n, obs$nbrReads)
    expect_equal(res$outCombined$filterSummary[, "f4_nbrNoValidOverlap"], 
                 sum(grepl("del", truth$truth$status)))  ## del
    expect_equal(res$outCombined$filterSummary[, "f6_nbrTooManyNinVar"], 
                 sum(grepl("N", truth$truth$status)))   ## N
    expect_equal(res$outCombined$filterSummary[, "nbrRetained"], 
                 sum(!grepl("N|del", truth$truth$status)))
    expect_equal(res$filtSeparate$barcode[, "f6_nbrTooManyNinVar"], 
                 sum(grepl("V1N", truth$truth$status)))
    expect_equal(res$filtSeparate$barcode[, "nbrRetained"], 
                 sum(!grepl("V1N", truth$truth$status)))
    expect_equal(res$filtSeparate$V2[, "f4_nbrNoValidOverlap"], 
                 sum(grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V2[, "nbrRetained"], 
                 sum(!grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V3[, "f4_nbrNoValidOverlap"], 
                 sum(grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V3[, "nbrRetained"], 
                 sum(!grepl("del", truth$truth$status)))
    
    ## ------------------------------------------------------------------------
    ## Filter out mutations in constant sequences
    res <- linkMultipleVariants(
        combinedDigestParams = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            fastqReverse = system.file("extdata", "multipleVariableRegions_R2.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CVCVCV", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            elementsReverse = "CVCV", elementLengthsReverse = c(6, 40, 10, -1),
            mergeForwardReverse = TRUE, minOverlap = 20, 
            maxOverlap = 30, maxFracMismatchOverlap = 0,
            revComplForward = FALSE, revComplReverse = TRUE,
            constantForward = truth$constFwd,
            constantReverse = truth$constRev,
            constantMaxDistForward = 0,
            constantMaxDistReverse = 0,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            verbose = verbose),
        barcode = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CVCSCS", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            constantForward = truth$constFwd,
            constantMaxDistForward = 0,
            avePhredMinForward = 20, variableCollapseMaxDist = 1, 
            variableCollapseMinReads = 1, variableCollapseMinRatio = 1,
            verbose = verbose),
        V2 = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            fastqReverse = system.file("extdata", "multipleVariableRegions_R2.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CSCVCS", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            elementsReverse = "CSCV", elementLengthsReverse = c(6, 40, 10, -1),
            mergeForwardReverse = TRUE, minOverlap = 10, 
            maxOverlap = 20, maxFracMismatchOverlap = 0,
            revComplForward = FALSE, revComplReverse = TRUE,
            constantForward = truth$constFwd,
            constantReverse = truth$constRev,
            constantMaxDistForward = 0,
            constantMaxDistReverse = 0,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            wildTypeForward = truth$WTV2,
            nbrMutatedCodonsMaxForward = -1, nbrMutatedBasesMaxForward = 1, 
            forbiddenMutatedCodonsForward = "", collapseToWTForward = FALSE,
            verbose = verbose, useTreeWTmatch = TRUE),
        V3 = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            fastqReverse = system.file("extdata", "multipleVariableRegions_R2.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CSCSCV", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            elementsReverse = "CVCS", elementLengthsReverse = c(6, 40, 10, -1),
            mergeForwardReverse = TRUE, minOverlap = 5, 
            maxOverlap = 15, maxFracMismatchOverlap = 0,
            revComplForward = FALSE, revComplReverse = TRUE,
            constantForward = truth$constFwd,
            constantReverse = truth$constRev,
            constantMaxDistForward = 0,
            constantMaxDistReverse = 0,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            wildTypeForward = truth$WTV3,
            nbrMutatedCodonsMaxForward = -1, nbrMutatedBasesMaxForward = 1, 
            forbiddenMutatedCodonsForward = "", collapseToWTForward = FALSE,
            verbose = verbose, useTreeWTmatch = TRUE)
    )
    
    ## Filter out reads with N/deletions in the variable sequence, or mutations in constant
    keep <- !grepl("C.*mut|N|del", truth$truth$status)
    
    ## Truth
    correct <- truth$truth[keep, ] %>%
        dplyr::group_by(trueBarcode, trueV2, trueV3) %>% dplyr::tally() %>%
        dplyr::arrange(dplyr::desc(n), trueBarcode, trueV2, trueV3)
    
    ## Obs
    obs <- res$countAggregated %>% dplyr::arrange(desc(nbrReads), barcode, V2, V3)
    
    expect_equal(sum(obs$nbrReads), sum(keep))
    expect_equal(correct$trueBarcode, obs$barcode)
    expect_equal(correct$trueV2, obs$V2)
    expect_equal(correct$trueV3, obs$V3)
    expect_equal(correct$n, obs$nbrReads)
    
    ## first filter out reads with deletion -> no valid overlap
    ## then reads with N in barcode -> too many N
    ## then reads with mutation in constant
    expect_equal(res$outCombined$filterSummary[, "f4_nbrNoValidOverlap"], 
                 sum(grepl("del", truth$truth$status)))
    expect_equal(res$outCombined$filterSummary[, "f6_nbrTooManyNinVar"], 
                 sum(grepl("N", truth$truth$status)))
    expect_equal(res$outCombined$filterSummary[, "f12_nbrTooManyMutConstant"], 
                 sum(grepl("C.*mut", truth$truth$status) & 
                         !grepl("N", truth$truth$status) & 
                         !grepl("del", truth$truth$status)))
    expect_equal(res$outCombined$filterSummary[, "nbrRetained"], 
                 sum(!grepl("C.*mut|N|del", truth$truth$status)))
    expect_equal(res$filtSeparate$barcode[, "f6_nbrTooManyNinVar"], 
                 sum(grepl("V1N", truth$truth$status)))
    expect_equal(res$filtSeparate$barcode[, "f12_nbrTooManyMutConstant"], 
                 sum(grepl("C[1-3]mut|del", truth$truth$status) & 
                         !grepl("V1N", truth$truth$status)))
    expect_equal(res$filtSeparate$barcode[, "nbrRetained"], 
                 sum(!grepl("V1N|C[1-3]mut|del", truth$truth$status)))
    expect_equal(res$filtSeparate$V2[, "f4_nbrNoValidOverlap"], 
                 sum(grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V2[, "f12_nbrTooManyMutConstant"], 
                 sum(grepl("C.*mut", truth$truth$status) & 
                         !grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V2[, "nbrRetained"], 
                 sum(!grepl("C.*mut|del", truth$truth$status)))
    expect_equal(res$filtSeparate$V3[, "f4_nbrNoValidOverlap"], 
                 sum(grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V3[, "f12_nbrTooManyMutConstant"], 
                 sum(grepl("C.*mut", truth$truth$status) & 
                         !grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V3[, "nbrRetained"], 
                 sum(!grepl("C.*mut|del", truth$truth$status)))

    ## ------------------------------------------------------------------------
    ## Collapse to WT
    res <- linkMultipleVariants(
        combinedDigestParams = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            fastqReverse = system.file("extdata", "multipleVariableRegions_R2.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CVCVCV", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            elementsReverse = "CVCV", elementLengthsReverse = c(6, 40, 10, -1),
            mergeForwardReverse = TRUE, minOverlap = 20, 
            maxOverlap = 30, maxFracMismatchOverlap = 0,
            revComplForward = FALSE, revComplReverse = TRUE,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            verbose = verbose),
        barcode = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CVCSCS", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            avePhredMinForward = 20, variableCollapseMaxDist = 1, 
            variableCollapseMinReads = 1, variableCollapseMinRatio = 1,
            verbose = verbose),
        V2 = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            fastqReverse = system.file("extdata", "multipleVariableRegions_R2.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CSCVCS", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            elementsReverse = "CSCV", elementLengthsReverse = c(6, 40, 10, -1),
            mergeForwardReverse = TRUE, minOverlap = 10, 
            maxOverlap = 20, maxFracMismatchOverlap = 0,
            revComplForward = FALSE, revComplReverse = TRUE,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            wildTypeForward = truth$WTV2,
            nbrMutatedCodonsMaxForward = -1, nbrMutatedBasesMaxForward = 1, 
            forbiddenMutatedCodonsForward = "", collapseToWTForward = TRUE,
            verbose = verbose, useTreeWTmatch = TRUE),
        V3 = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            fastqReverse = system.file("extdata", "multipleVariableRegions_R2.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CSCSCV", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            elementsReverse = "CVCS", elementLengthsReverse = c(6, 40, 10, -1),
            mergeForwardReverse = TRUE, minOverlap = 5, 
            maxOverlap = 15, maxFracMismatchOverlap = 0,
            revComplForward = FALSE, revComplReverse = TRUE,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            wildTypeForward = truth$WTV3,
            nbrMutatedCodonsMaxForward = -1, nbrMutatedBasesMaxForward = 1, 
            forbiddenMutatedCodonsForward = "", collapseToWTForward = TRUE,
            verbose = verbose, useTreeWTmatch = TRUE)
    )
    
    ## Filter out only reads with N/deletions in the variable sequence
    keep <- !grepl("N|del", truth$truth$status)
    
    ## Truth
    correct <- truth$truth[keep, ] %>%
        dplyr::mutate(trueV2 = gsub("\\..*", "", trueV2),
                      trueV3 = gsub("\\..*", "", trueV3)) %>%
        dplyr::group_by(trueBarcode, trueV2, trueV3) %>% dplyr::tally() %>%
        dplyr::arrange(dplyr::desc(n), trueBarcode, trueV2, trueV3)
    
    ## Obs
    obs <- res$countAggregated %>% dplyr::arrange(desc(nbrReads), barcode, V2, V3)
    
    expect_equal(sum(obs$nbrReads), sum(keep))
    expect_equal(correct$trueBarcode, obs$barcode)
    expect_equal(correct$trueV2, obs$V2)
    expect_equal(correct$trueV3, obs$V3)
    expect_equal(correct$n, obs$nbrReads)
    expect_equal(res$outCombined$filterSummary[, "f4_nbrNoValidOverlap"], 
                 sum(grepl("del", truth$truth$status)))  ## del
    expect_equal(res$outCombined$filterSummary[, "f6_nbrTooManyNinVar"], 
                 sum(grepl("N", truth$truth$status)))   ## N
    expect_equal(res$outCombined$filterSummary[, "nbrRetained"], 
                 sum(!grepl("N|del", truth$truth$status)))
    expect_equal(res$filtSeparate$barcode[, "f6_nbrTooManyNinVar"], 
                 sum(grepl("V1N", truth$truth$status)))
    expect_equal(res$filtSeparate$barcode[, "nbrRetained"], 
                 sum(!grepl("V1N", truth$truth$status)))
    expect_equal(res$filtSeparate$V2[, "f4_nbrNoValidOverlap"], 
                 sum(grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V2[, "nbrRetained"], 
                 sum(!grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V3[, "f4_nbrNoValidOverlap"], 
                 sum(grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V3[, "nbrRetained"], 
                 sum(!grepl("del", truth$truth$status)))
    
    ## ------------------------------------------------------------------------
    ## Don't merge fwd and rev
    res <- linkMultipleVariants(
        combinedDigestParams = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            fastqReverse = system.file("extdata", "multipleVariableRegions_R2.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CVCVCS", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            elementsReverse = "CVCS", elementLengthsReverse = c(6, 40, 10, -1),
            mergeForwardReverse = FALSE,
            revComplForward = FALSE, revComplReverse = TRUE,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            verbose = verbose),
        barcode = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CVCSCS", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            avePhredMinForward = 20, variableCollapseMaxDist = 1, 
            variableCollapseMinReads = 1, variableCollapseMinRatio = 1,
            verbose = verbose),
        V2 = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CSCVCS", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            mergeForwardReverse = FALSE,
            revComplForward = FALSE, revComplReverse = TRUE,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            wildTypeForward = truth$WTV2,
            nbrMutatedCodonsMaxForward = -1, nbrMutatedBasesMaxForward = 1, 
            forbiddenMutatedCodonsForward = "", collapseToWTForward = FALSE,
            verbose = verbose, useTreeWTmatch = TRUE),
        V3 = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R2.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CVCS", elementLengthsForward = c(6, 40, 10, -1),
            mergeForwardReverse = FALSE,
            revComplForward = TRUE,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            wildTypeForward = truth$WTV3,
            nbrMutatedCodonsMaxForward = -1, nbrMutatedBasesMaxForward = 1, 
            forbiddenMutatedCodonsForward = "", collapseToWTForward = FALSE,
            verbose = verbose, useTreeWTmatch = TRUE)
    )
    
    ## Filter out only reads with N/deletions in the variable sequence
    keep <- !grepl("N|del", truth$truth$status)
    
    ## Truth
    correct <- truth$truth[keep, ] %>%
        dplyr::group_by(trueBarcode, trueV2, trueV3) %>% dplyr::tally() %>%
        dplyr::arrange(dplyr::desc(n), trueBarcode, trueV2, trueV3)
    
    ## Obs
    obs <- res$countAggregated %>% dplyr::arrange(desc(nbrReads), barcode, V2, V3)
    
    expect_equal(sum(obs$nbrReads), sum(keep))
    expect_equal(correct$trueBarcode, obs$barcode)
    expect_equal(correct$trueV2, obs$V2)
    expect_equal(correct$trueV3, obs$V3)
    expect_equal(correct$n, obs$nbrReads)
    expect_equal(res$outCombined$filterSummary[, "f6_nbrTooManyNinVar"], 
                 sum(grepl("N", truth$truth$status)))
    expect_equal(res$outCombined$filterSummary[, "nbrRetained"], 
                 sum(!grepl("N", truth$truth$status)))
    expect_equal(res$filtSeparate$barcode[, "f6_nbrTooManyNinVar"], 
                 sum(grepl("V1N", truth$truth$status)))
    expect_equal(res$filtSeparate$barcode[, "nbrRetained"], 
                 sum(!grepl("V1N", truth$truth$status)))
    expect_equal(res$filtSeparate$V2[, "f10b_nbrTooManyMutBases"], 
                 sum(grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V2[, "nbrRetained"], 
                 sum(!grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V3[, "nbrRetained"], 100)
    
    ## ------------------------------------------------------------------------
    ## Don't allow any mutations
    res <- linkMultipleVariants(
        combinedDigestParams = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            fastqReverse = system.file("extdata", "multipleVariableRegions_R2.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CVCVCV", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            elementsReverse = "CVCV", elementLengthsReverse = c(6, 40, 10, -1),
            mergeForwardReverse = TRUE, minOverlap = 20, 
            maxOverlap = 30, maxFracMismatchOverlap = 0,
            revComplForward = FALSE, revComplReverse = TRUE,
            constantForward = truth$constFwd,
            constantReverse = truth$constRev,
            constantMaxDistForward = 0,
            constantMaxDistReverse = 0,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            verbose = verbose),
        barcode = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CVCSCS", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            constantForward = truth$constFwd,
            constantMaxDistForward = 0,
            avePhredMinForward = 20, variableCollapseMaxDist = 0, 
            variableCollapseMinReads = 1, variableCollapseMinRatio = 1,
            verbose = verbose),
        V2 = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            fastqReverse = system.file("extdata", "multipleVariableRegions_R2.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CSCVCS", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            elementsReverse = "CSCV", elementLengthsReverse = c(6, 40, 10, -1),
            mergeForwardReverse = TRUE, minOverlap = 10, 
            maxOverlap = 20, maxFracMismatchOverlap = 0,
            revComplForward = FALSE, revComplReverse = TRUE,
            constantForward = truth$constFwd,
            constantReverse = truth$constRev,
            constantMaxDistForward = 0,
            constantMaxDistReverse = 0,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            wildTypeForward = truth$WTV2,
            nbrMutatedCodonsMaxForward = -1, nbrMutatedBasesMaxForward = 0, 
            forbiddenMutatedCodonsForward = "", collapseToWTForward = FALSE,
            verbose = verbose, useTreeWTmatch = TRUE),
        V3 = list(
            fastqForward = system.file("extdata", "multipleVariableRegions_R1.fastq.gz",
                                       package = "mutscan"),
            fastqReverse = system.file("extdata", "multipleVariableRegions_R2.fastq.gz",
                                       package = "mutscan"),
            elementsForward = "CSCSCV", elementLengthsForward = c(6, 24, 10, 30, 10, -1),
            elementsReverse = "CVCS", elementLengthsReverse = c(6, 40, 10, -1),
            mergeForwardReverse = TRUE, minOverlap = 5, 
            maxOverlap = 15, maxFracMismatchOverlap = 0,
            revComplForward = FALSE, revComplReverse = TRUE,
            constantForward = truth$constFwd,
            constantReverse = truth$constRev,
            constantMaxDistForward = 0,
            constantMaxDistReverse = 0,
            avePhredMinForward = 20, avePhredMinReverse = 20,
            wildTypeForward = truth$WTV3,
            nbrMutatedCodonsMaxForward = -1, nbrMutatedBasesMaxForward = 0, 
            forbiddenMutatedCodonsForward = "", collapseToWTForward = FALSE,
            verbose = verbose, useTreeWTmatch = TRUE)
    )
    
    ## Filter out reads with N/deletions in the variable sequence, or mutations in constant
    keep <- !grepl("C.*mut|N|del|V2|V3", truth$truth$status)
    
    ## Truth
    correct <- truth$truth[keep, ] %>%
        dplyr::group_by(obsBarcode, trueV2, trueV3) %>% dplyr::tally() %>%
        dplyr::arrange(dplyr::desc(n), obsBarcode, trueV2, trueV3)
    
    ## Obs
    obs <- res$countAggregated %>% dplyr::arrange(desc(nbrReads), barcode, V2, V3)
    
    expect_equal(sum(obs$nbrReads), sum(keep))
    expect_equal(correct$obsBarcode, obs$barcode)
    expect_equal(correct$trueV2, obs$V2)
    expect_equal(correct$trueV3, obs$V3)
    expect_equal(correct$n, obs$nbrReads)
    
    expect_equal(res$outCombined$filterSummary[, "f4_nbrNoValidOverlap"], 
                 sum(grepl("del", truth$truth$status)))
    expect_equal(res$outCombined$filterSummary[, "f6_nbrTooManyNinVar"], 
                 sum(grepl("N", truth$truth$status)))
    expect_equal(res$outCombined$filterSummary[, "f12_nbrTooManyMutConstant"], 
                 sum(grepl("C.*mut", truth$truth$status) & 
                         !grepl("N", truth$truth$status) & 
                         !grepl("del", truth$truth$status)))
    expect_equal(res$outCombined$filterSummary[, "nbrRetained"], 
                 sum(!grepl("C.*mut|N|del", truth$truth$status)))
    expect_equal(res$filtSeparate$barcode[, "f6_nbrTooManyNinVar"], 
                 sum(grepl("V1N", truth$truth$status)))
    expect_equal(res$filtSeparate$barcode[, "f12_nbrTooManyMutConstant"], 
                 sum(grepl("C[1-3]mut|del", truth$truth$status) & 
                         !grepl("V1N", truth$truth$status)))
    expect_equal(res$filtSeparate$barcode[, "nbrRetained"], 
                 sum(!grepl("V1N|C[1-3]mut|del", truth$truth$status)))
    expect_equal(res$filtSeparate$V2[, "f4_nbrNoValidOverlap"], 
                 sum(grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V2[, "f10b_nbrTooManyMutBases"],
                 sum(grepl("V2v[0-9]\\.", truth$truth$status) & 
                         !grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V2[, "f12_nbrTooManyMutConstant"], 
                 sum(grepl("C.*mut", truth$truth$status) & 
                         !grepl("del", truth$truth$status) & 
                         !grepl("V2v[0-9]\\.", truth$truth$status)))
    expect_equal(res$filtSeparate$V2[, "nbrRetained"], 
                 sum(!grepl("C.*mut|del|V2v[0-9]\\.", truth$truth$status)))
    expect_equal(res$filtSeparate$V3[, "f4_nbrNoValidOverlap"], 
                 sum(grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V3[, "f10b_nbrTooManyMutBases"],
                 sum(grepl("V3\\.", truth$truth$status) & 
                         !grepl("del", truth$truth$status)))
    expect_equal(res$filtSeparate$V3[, "f12_nbrTooManyMutConstant"], 
                 sum(grepl("C.*mut", truth$truth$status) & 
                         !grepl("del", truth$truth$status) & 
                         !grepl("V3\\.", truth$truth$status)))
    expect_equal(res$filtSeparate$V3[, "nbrRetained"], 
                 sum(!grepl("C.*mut|del|V3\\.", truth$truth$status)))
    
})


