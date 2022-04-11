test_that("associateBarcodes works", {
    ## Read truth (for WT sequences)
    truth <- readRDS(system.file("extdata", "multipleVariableRegions_truth.rds",
                                 package = "mutscan"))
    verbose <- FALSE
    
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
            verbose = verbose, useTreeWTmatch = TRUE),
        collapseToAA = FALSE
    )
    
    ## Fails with the wrong input
    ## ------------------------------------------------------------------------
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = 1,
                                   variantCol = "V2", countCol = "nbrReads", minCount = 2, 
                                   minFracOfBarcodeCount = 0.25),
                 "'barcodeCol' must be of class 'character'")
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = c("barcode", "V2"),
                                   variantCol = "V2", countCol = "nbrReads", minCount = 2, 
                                   minFracOfBarcodeCount = 0.25),
                 "'barcodeCol' must have length 1")
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "missing",
                                   variantCol = "V2", countCol = "nbrReads", minCount = 2, 
                                   minFracOfBarcodeCount = 0.25),
                 "All values in 'barcodeCol' must be one of")
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "nbrReads",
                                   variantCol = "V2", countCol = "nbrReads", minCount = 2, 
                                   minFracOfBarcodeCount = 0.25),
                 "'[[countTablebarcodeCol' must be of class 'character'", fixed = TRUE)
    
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "barcode",
                                   variantCol = 1, countCol = "nbrReads", minCount = 2, 
                                   minFracOfBarcodeCount = 0.25),
                 "'variantCol' must be of class 'character'")
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "barcode",
                                   variantCol = c("V2", "barcode"), 
                                   countCol = "nbrReads", minCount = 2, 
                                   minFracOfBarcodeCount = 0.25),
                 "'variantCol' must have length 1")
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "barcode",
                                   variantCol = "missing", countCol = "nbrReads", minCount = 2, 
                                   minFracOfBarcodeCount = 0.25),
                 "All values in 'variantCol' must be one of")
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "barcode",
                                   variantCol = "nbrReads", countCol = "nbrReads", minCount = 2, 
                                   minFracOfBarcodeCount = 0.25),
                 "'[[countTablevariantCol' must be of class 'character'", fixed = TRUE)
    
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "barcode",
                                   variantCol = "V2", countCol = 1, minCount = 2, 
                                   minFracOfBarcodeCount = 0.25),
                 "'countCol' must be of class 'character'")
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "barcode",
                                   variantCol = "V2", countCol = c("barcode", "V2"), minCount = 2, 
                                   minFracOfBarcodeCount = 0.25),
                 "'countCol' must have length 1")
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "barcode",
                                   variantCol = "V2", countCol = "missing", minCount = 2, 
                                   minFracOfBarcodeCount = 0.25),
                 "All values in 'countCol' must be one of")
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "barcode",
                                   variantCol = "V2", countCol = "barcode", minCount = 2, 
                                   minFracOfBarcodeCount = 0.25),
                 "'[[countTablecountCol' must be of class 'numeric'", fixed = TRUE)
    
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "barcode",
                                   variantCol = "V2", countCol = "nbrReads", minCount = "2", 
                                   minFracOfBarcodeCount = 0.25),
                 "'minCount' must be of class 'numeric'")
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "barcode",
                                   variantCol = "V2", countCol = "nbrReads", minCount = c(1, 2), 
                                   minFracOfBarcodeCount = 0.25),
                 "'minCount' must have length 1")
    
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "barcode",
                                   variantCol = "V2", countCol = "nbrReads", minCount = 2, 
                                   minFracOfBarcodeCount = "0.25"),
                 "'minFracOfBarcodeCount' must be of class 'numeric'")
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "barcode",
                                   variantCol = "V2", countCol = "nbrReads", minCount = 2, 
                                   minFracOfBarcodeCount = c(0.1, 0.25)),
                 "'minFracOfBarcodeCount' must have length 1")
    expect_error(associateBarcodes(res$countAggregated, barcodeCol = "barcode",
                                   variantCol = "V2", countCol = "nbrReads", minCount = 2, 
                                   minFracOfBarcodeCount = 2.5),
                 "'minFracOfBarcodeCount' must be within [0,1]", fixed = TRUE)

    
    ## Works with the correct input
    ## ------------------------------------------------------------------------
    ## Run barcode association
    out <- associateBarcodes(res$countAggregated, barcodeCol = "barcode",
                             variantCol = "V2", minCount = 2, 
                             minFracOfBarcodeCount = 0.25)
    
    ## Calculate expected values manually
    correct <- truth$truth[!grepl("N|del", truth$truth$status), ] %>%
        ## get frequency for each barcode/variant pair
        dplyr::group_by(trueBarcode, trueV2) %>%
        dplyr::tally() %>% 
        ## filter rare pairs
        dplyr::filter(n >= 2) %>%
        ## filter pairs with too low fraction of total bc counts
        dplyr::group_by(trueBarcode) %>% 
        dplyr::mutate(relFrac = n/sum(n)) %>%
        dplyr::filter(relFrac >= 0.25) %>%
        ## keep top variant for each barcode %>%
        dplyr::arrange(dplyr::desc(n)) %>%
        dplyr::filter(!duplicated(trueBarcode))
        
    expect_equal(sum(out$counts$keep), nrow(correct))
    expect_true(all(paste0(out$counts$barcode[out$counts$keep], 
                           out$counts$V2[out$counts$keep]) %in% 
                        paste0(correct$trueBarcode, correct$trueV2)))
    tmp <- out$counts[match(paste0(correct$trueBarcode, correct$trueV2),
                            paste0(out$counts$barcode, out$counts$V2)), ]
    expect_equal(tmp$nbrReads, correct$n)
    expect_equal(tmp$fracOfTotalBarcodeCountFilt, correct$relFrac)
})
