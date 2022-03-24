test_that("relabeling mutants works", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(exprs = matrix(rnorm(14), nrow = 7))
    )
    rownames(se) <- c("f.0.WT", "f.1.CGG_r.2.AGT", "f.1.AGT_f.3.GTA",
                      "f.2.G", "f.3.G_r.2.A", "f.3.G_r.0.WT",
                      "f.2.GGT_f.3.AGT_r.1.GAG")
    colnames(se) <- c("sample1", "sample2")
    
    convTable <- data.frame(seqname = c("f", "f", "f", "f", "r", "r", "r"),
                            position = c(0, 1, 2, 3, 0, 1, 2),
                            name = LETTERS[seq_len(7)])
    
    ## Fails with the wrong arguments
    expect_error(relabelMutPositions(se = 1, conversionTable = convTable,
                                     mutNameDelimiter = "."),
                 "'se' must be of class 'SummarizedExperiment'")
    expect_error(relabelMutPositions(se = se, conversionTable = 1,
                                     mutNameDelimiter = "."),
                 "'colnamesconversionTable' must not be NULL")
    expect_error(relabelMutPositions(se = se, conversionTable = convTable[, -1],
                                     mutNameDelimiter = "."),
                 "%in% colnames(conversionTable)) is not TRUE", fixed = TRUE)
    expect_error(relabelMutPositions(se = se, conversionTable = convTable[, -2],
                                     mutNameDelimiter = "."),
                 "%in% colnames(conversionTable)) is not TRUE", fixed = TRUE)
    expect_error(relabelMutPositions(se = se, conversionTable = convTable[, -3],
                                     mutNameDelimiter = "."),
                 "%in% colnames(conversionTable)) is not TRUE", fixed = TRUE)
    expect_error(relabelMutPositions(se = se, conversionTable = convTable,
                                     mutNameDelimiter = 1),
                 "'mutNameDelimiter' must be of class 'character'")
    expect_error(relabelMutPositions(se = se, conversionTable = convTable,
                                     mutNameDelimiter = c(".", "_")),
                 "'mutNameDelimiter' must have length 1")
    
    ## Works with the right arguments
    out <- relabelMutPositions(se, convTable, ".")
    expect_equal(rownames(out),
                 c("f.A.WT", "f.B.CGG_r.G.AGT", "f.B.AGT_f.D.GTA",
                   "f.C.G", "f.D.G_r.G.A", "f.D.G_r.E.WT", 
                   "f.C.GGT_f.D.AGT_r.F.GAG"))
    expect_equal(colnames(out), colnames(se))
    
    ## Missing entries in conversion table
    convTable <- convTable[c(2, 4, 6, 7), ]
    
    out <- relabelMutPositions(se, convTable, ".")
    expect_equal(rownames(out),
                 c("f.0.WT", "f.B.CGG_r.G.AGT", "f.B.AGT_f.D.GTA",
                   "f.2.G", "f.D.G_r.G.A", "f.D.G_r.0.WT", 
                   "f.2.GGT_f.D.AGT_r.F.GAG"))
    expect_equal(colnames(out), colnames(se))
})
