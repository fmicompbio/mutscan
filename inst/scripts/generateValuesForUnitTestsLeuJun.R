suppressPackageStartupMessages({
  library(Biostrings)
  library(ShortRead)
})

processReadsLeuJun <- function(Ldef) {
  ## Read data
  fq1 <- Biostrings::readQualityScaledDNAStringSet(Ldef$fastqForward)
  fq2 <- Biostrings::readQualityScaledDNAStringSet(Ldef$fastqReverse)
  
  ## Number of reads
  message("Number of reads: ", length(fq1))  ## 1000
  
  ## Primer sequences
  prim <- Biostrings::vcountPattern(Ldef$primerForward, fq1) == 0 |
    Biostrings::vcountPattern(Ldef$primerReverse, fq2) == 0
  message("Number of reads without primers: ", sum(prim))  ## 126
  fq1 <- fq1[!prim]
  fq2 <- fq2[!prim]
  
  ## Average quality in variable sequences too low
  endprimFwd <- unlist(end(Biostrings::vmatchPattern(Ldef$primerForward, fq1)))
  varForward <- Biostrings::subseq(fq1, start = endprimFwd + 1, width = width(fq1) - endprimFwd)
  endprimRev <- unlist(end(Biostrings::vmatchPattern(Ldef$primerReverse, fq2)))
  varReverse <- Biostrings::subseq(fq2, start = endprimRev + 1, width = Ldef$variableLengthReverse)
  avgQualForward <- ShortRead::alphabetScore(Biostrings::quality(varForward))/
    BiocGenerics::width(varForward)
  avgQualReverse <- ShortRead::alphabetScore(Biostrings::quality(varReverse))/
    BiocGenerics::width(varReverse)
  lowqual <- avgQualForward < Ldef$avePhredMinForward | avgQualReverse < Ldef$avePhredMinReverse
  message("Number of reads with low average quality: ", sum(lowqual))  # 0
  fq1 <- fq1[!lowqual]
  fq2 <- fq2[!lowqual]
  varForward <- varForward[!lowqual]
  varReverse <- varReverse[!lowqual]
  
  ## Too many Ns in variable sequences
  nbrNForward <- Biostrings::vcountPattern("N", varForward, fixed = TRUE)
  nbrNReverse <- Biostrings::vcountPattern("N", varReverse, fixed = TRUE)
  toomanynvar <- nbrNForward > Ldef$variableNMaxForward | nbrNReverse > Ldef$variableNMaxReverse
  message("Number of reads with too many N in variable sequence: ", sum(toomanynvar))  # 76
  fq1 <- fq1[!toomanynvar]
  fq2 <- fq2[!toomanynvar]
  varForward <- varForward[!toomanynvar]
  varReverse <- varReverse[!toomanynvar]
  
  ## Compare to wildtype, forward
  bestmatch <- DNAStringSet()
  for (i in seq_along(varForward)) {
    s <- DNAStringSet(Ldef$wildTypeForward)
    s <- Biostrings::subseq(s, start = 1, width = length(varForward[[i]]))
    p <- DNAStringSet(rep(as.character(varForward[[i]]), length(s)))
    pattern_width <- width(p)
    subject_width <- width(s)
    unlisted_ans <- which(as.raw(unlist(p)) == as.raw(unlist(s)))
    breakpoints <- cumsum(pattern_width)
    ans_eltlens <- tabulate(findInterval(unlisted_ans - 1L,
                                         breakpoints) + 1L,
                            nbins = length(p))
    skeleton <- IRanges::PartitioningByEnd(cumsum(ans_eltlens))
    nucleotideMatches <- BiocGenerics::relist(unlisted_ans, skeleton)
    offsets <- c(0L, breakpoints[-length(breakpoints)])
    matchPositions <- nucleotideMatches - offsets
    bestmatch <- c(bestmatch, s[which.max(sapply(matchPositions, length))])
  }
  p <- bestmatch
  pattern_width <- width(p)
  subject_width <- width(varForward)
  unlisted_ans <- which(as.raw(unlist(p)) != as.raw(unlist(varForward)))
  breakpoints <- cumsum(pattern_width)
  ans_eltlens <- tabulate(findInterval(unlisted_ans - 1L,
                                       breakpoints) + 1L,
                          nbins = length(p))
  skeleton <- IRanges::PartitioningByEnd(cumsum(ans_eltlens))
  nucleotideMismatches <- BiocGenerics::relist(unlisted_ans, skeleton)
  offsets <- c(0L, breakpoints[-length(breakpoints)])
  mismatchPositions <- nucleotideMismatches - offsets
  mismatchQualities <- Biostrings::quality(varForward)[mismatchPositions]
  lowqualmm <- min(as(mismatchQualities, "IntegerList")) < Ldef$mutatedPhredMinForward
  message("Number of reads with low-quality mismatch bases, forward: ", sum(lowqualmm))  ## 0
  fq1 <- fq1[!lowqualmm]
  fq2 <- fq2[!lowqualmm]
  varForward <- varForward[!lowqualmm]
  varReverse <- varReverse[!lowqualmm]
  mismatchPositions <- mismatchPositions[!lowqualmm]
  p <- p[!lowqualmm]
  ## Number of mutated codons
  codonMismatches <- (mismatchPositions - 1) %/% 3 + 1
  nbrMutCodons <- S4Vectors::elementNROWS(unique(codonMismatches))
  toomanycodons <- nbrMutCodons > Ldef$nbrMutatedCodonsMaxForward
  message("Number of reads with too many mutated codons, forward: ", sum(toomanycodons))  ## 58
  fq1 <- fq1[!toomanycodons]
  fq2 <- fq2[!toomanycodons]
  codonMismatches <- codonMismatches[!toomanycodons]
  varForward <- varForward[!toomanycodons]
  varReverse <- varReverse[!toomanycodons]
  p <- p[!toomanycodons]
  ## Forbidden codons
  # uniqueMutCodons <- unique(codonMismatches)
  # mutCodons <- S4Vectors::endoapply(uniqueMutCodons, 
  #                                   function(v) rep(3 * v, each = 3) + c(-2, -1, 0))
  # mutCodons <- as(varForward, "DNAStringSet")[mutCodons]
  # matchForbidden <- Biostrings::vmatchPattern(
  #   pattern = Ldef$forbiddenMutatedCodonsForward, 
  #   as(mutCodons, "DNAStringSet"), 
  #   fixed = FALSE)
  # forbiddencodon <- sapply(BiocGenerics::relist(start(unlist(matchForbidden)) %% 3==1, 
  #                                               matchForbidden), function(i) any(i))
  # message("Number of reads with forbidden codons, forward: ", sum(forbiddencodon))  ## 6
  # fq1 <- fq1[!forbiddencodon]
  # fq2 <- fq2[!forbiddencodon]
  # mutCodons <- mutCodons[!forbiddencodon]
  # uniqueMutCodons <- uniqueMutCodons[!forbiddencodon]
  # splitMutCodons <- as(strsplit(gsub("([^.]{3})", "\\1\\.", 
  #                                    as.character(mutCodons)), ".", fixed = TRUE), 
  #                      "CharacterList")
  # encodedMutCodonsForward <- S4Vectors::unstrsplit(BiocGenerics::relist(
  #   paste0("f", unlist(uniqueMutCodons), unlist(splitMutCodons)), 
  #   uniqueMutCodons),
  #   sep = "_")
  encodedMutCodonsForward <- names(p)
  
  
  ## Compare to wildtype, reverse
  p <- DNAStringSet(rep(Ldef$wildTypeReverse, length(varReverse)))
  pattern_width <- width(p)
  subject_width <- width(varReverse)
  unlisted_ans <- which(as.raw(unlist(p)) != as.raw(unlist(varReverse)))
  breakpoints <- cumsum(pattern_width)
  ans_eltlens <- tabulate(findInterval(unlisted_ans - 1L,
                                       breakpoints) + 1L,
                          nbins = length(p))
  skeleton <- IRanges::PartitioningByEnd(cumsum(ans_eltlens))
  nucleotideMismatches <- BiocGenerics::relist(unlisted_ans, skeleton)
  offsets <- c(0L, breakpoints[-length(breakpoints)])
  mismatchPositions <- nucleotideMismatches - offsets
  mismatchQualities <- Biostrings::quality(varReverse)[mismatchPositions]
  lowqualmm <- min(as(mismatchQualities, "IntegerList")) < Ldef$mutatedPhredMinReverse
  message("Number of reads with low-quality mismatch bases, reverse: ", sum(lowqualmm))  ## 0
  fq1 <- fq1[!lowqualmm]
  fq2 <- fq2[!lowqualmm]
  mismatchPositions <- mismatchPositions[!lowqualmm]
  varForward <- varForward[!lowqualmm]
  varReverse <- varReverse[!lowqualmm]
  encodedMutCodonsForward <- encodedMutCodonsForward[!lowqualmm]
  ## Number of mutated codons
  codonMismatches <- (mismatchPositions - 1) %/% 3 + 1
  nbrMutCodons <- S4Vectors::elementNROWS(unique(codonMismatches))
  toomanycodons <- nbrMutCodons > Ldef$nbrMutatedCodonsMaxReverse
  message("Number of reads with too many mutated codons, reverse: ", sum(toomanycodons))  ## 137
  fq1 <- fq1[!toomanycodons]
  fq2 <- fq2[!toomanycodons]
  codonMismatches <- codonMismatches[!toomanycodons]
  varForward <- varForward[!toomanycodons]
  varReverse <- varReverse[!toomanycodons]
  encodedMutCodonsForward <- encodedMutCodonsForward[!toomanycodons]
  ## Forbidden codons
  uniqueMutCodons <- unique(codonMismatches)
  mutCodons <- S4Vectors::endoapply(uniqueMutCodons, 
                                    function(v) rep(3 * v, each = 3) + c(-2, -1, 0))
  mutCodons <- as(varReverse, "DNAStringSet")[mutCodons]
  matchForbidden <- Biostrings::vmatchPattern(
    pattern = Ldef$forbiddenMutatedCodonsReverse, 
    as(mutCodons, "DNAStringSet"), 
    fixed = FALSE)
  forbiddencodon <- sapply(BiocGenerics::relist(start(unlist(matchForbidden)) %% 3==1, 
                                                matchForbidden), function(i) any(i))
  message("Number of reads with forbidden codons, reverse: ", sum(forbiddencodon))  ## 3
  fq1 <- fq1[!forbiddencodon]
  fq2 <- fq2[!forbiddencodon]
  encodedMutCodonsForward <- encodedMutCodonsForward[!forbiddencodon]
  mutCodons <- mutCodons[!forbiddencodon]
  uniqueMutCodons <- uniqueMutCodons[!forbiddencodon]
  splitMutCodons <- as(strsplit(gsub("([^.]{3})", "\\1\\.", 
                                     as.character(mutCodons)), ".", fixed = TRUE), 
                       "CharacterList")
  encodedMutCodonsReverse <- S4Vectors::unstrsplit(BiocGenerics::relist(
    paste0("r", unlist(uniqueMutCodons), unlist(splitMutCodons)), 
    uniqueMutCodons),
    sep = "_")
  mutNames <- paste(encodedMutCodonsForward, encodedMutCodonsReverse, sep = "_")
  mutNames <- gsub("^_", "", gsub("_$", "", mutNames))
  mutNames[mutNames <- ""] <- "WT"
  message("Number of variants seen twice: ", sum(table(mutNames) == 2))  ## 5
  message("Variants seen twice: ")
  print(which(table(mutNames) == 2))  ## ATF6B_r13CCC, ATF7_r26CTC, CEBPA, FOSB, FOSL1_r3GGC
  
  ## Retained reads
  message("Number of retained reads: ", length(fq1))  ## 600
}

## ----------------------------------------------------------------------------
leu <- c(ATF2 = "GATCCTGATGAAAAAAGGAGAAAGTTTTTAGAGCGAAATAGAGCAGCAGCTTCAAGATGCCGACAAAAAAGGAAAGTCTGGGTTCAGTCTTTAGAGAAGAAAGCTGAAGACTTGAGTTCATTAAATGGTCAGCTGCAGAGTGAAGTCACCCTGCTGAGAAATGAAGTGGCACAGCTGAAACAGCTTCTTCTGGCT",
         ATF7 = "GATCCAGATGAGCGACGGCAGCGCTTTCTGGAGCGCAACCGGGCTGCAGCCTCCCGCTGCCGCCAAAAGCGAAAGCTGTGGGTGTCCTCCCTAGAGAAGAAGGCCGAAGAACTCACTTCTCAGAACATTCAGCTGAGTAATGAAGTCACATTACTACGCAATGAGGTGGCCCAGTTGAAACAGCTACTGTTAGCT",
         CREB5 = "GATCCGGACGAGAGGCGGCGGAAATTTCTGGAACGGAACCGGGCAGCTGCCACCCGCTGCAGACAGAAGAGGAAGGTCTGGGTGATGTCATTGGAAAAGAAAGCAGAAGAACTCACCCAGACAAACATGCAGCTTCAGAATGAAGTGTCTATGTTGAAAAATGAGGTGGCCCAGCTGAAACAGTTGTTGTTAACA",
         ATF3 = "GAAGAAGATGAAAGGAAAAAGAGGCGACGAGAAAGAAATAAGATTGCAGCTGCAAAGTGCCGAAACAAGAAGAAGGAGAAGACGGAGTGCCTGCAGAAAGAGTCGGAGAAGCTGGAAAGTGTGAATGCTGAACTGAAGGCTCAGATTGAGGAGCTCAAGAACGAGAAGCAGCATTTGATATACATGCTCAACCTT",
         JDP2 = "GAGGAAGAGGAGCGAAGGAAAAGGCGCCGGGAGAAGAACAAAGTCGCAGCAGCCCGATGCCGGAACAAGAAGAAGGAGCGCACGGAGTTTCTGCAGCGGGAATCCGAGCGGCTGGAACTCATGAACGCAGAGCTGAAGACCCAGATTGAGGAGCTGAAGCAGGAGCGGCAGCAGCTCATCCTGATGCTGAACCGA",
         ATF4 = "GAGAAACTGGATAAGAAGCTGAAAAAAATGGAGCAAAACAAGACAGCAGCCACTAGGTACCGCCAGAAGAAGAGGGCGGAGCAGGAGGCTCTTACTGGTGAGTGCAAAGAGCTGGAAAAGAAGAACGAGGCTCTAAAAGAGAGGGCGGATTCCCTGGCCAAGGAGATCCAGTACCTGAAAGATTTGATAGAAGAG",
         ATF5 = "ACCCGAGGGGACCGCAAGCAAAAGAAGAGAGACCAGAACAAGTCGGCGGCTCTGAGGTACCGCCAGCGGAAGCGGGCAGAGGGTGAGGCCCTGGAGGGCGAGTGCCAGGGGCTGGAGGCACGGAATCGCGAGCTGAAGGAACGGGCAGAGTCCGTGGAGCGCGAGATCCAGTACGTCAAGGACCTGCTCATCGAG",
         CREBZF = "AGTCCCCGGAAGGCGGCGGCGGCCGCTGCCCGCCTTAATCGACTGAAGAAGAAGGAGTACGTGATGGGGCTGGAGAGTCGAGTCCGGGGTCTGGCAGCCGAGAACCAGGAGCTGCGGGCCGAGAATCGGGAGCTGGGCAAACGCGTACAGGCACTGCAGGAGGAGAGTCGCTACCTACGGGCAGTCTTAGCCAAC",
         BATF2 = "CCCAAGGAGCAACAAAGGCAGCTGAAGAAGCAGAAGAACCGGGCAGCCGCCCAGCGAAGCCGGCAGAAGCACACAGACAAGGCAGACGCCCTGCACCAGCAGCACGAGTCTCTGGAAAAAGACAACCTCGCCCTGCGGAAGGAGATCCAGTCCCTGCAGGCCGAGCTGGCGTGGTGGAGCCGGACCCTGCACGTG",
         BATF3 = "GAGGATGATGACAGGAAGGTCCGAAGGAGAGAAAAAAACCGAGTTGCTGCTCAGAGAAGTCGGAAGAAGCAGACCCAGAAGGCTGACAAGCTCCATGAGGAATATGAGAGCCTGGAGCAAGAAAACACCATGCTGCGGAGAGAGATCGGGAAGCTGACAGAGGAGCTGAAGCACCTGACAGAGGCACTGAAGGAG",
         CEBPE = "AAAGATAGCCTTGAGTACCGGCTGAGGCGGGAGCGCAACAACATCGCCGTGCGCAAGAGCCGAGACAAGGCCAAGAGGCGCATTCTGGAGACGCAGCAGAAGGTGCTGGAGTACATGGCAGAGAACGAGCGCCTCCGCAGCCGCGTGGAGCAGCTCACCCAGGAGCTAGACACCCTCCGCAACCTCTTCCGCCAG",
         BACH1 = "CTGGATTGTATCCATGATATTCGAAGAAGAAGTAAAAACAGAATTGCTGCACAGCGCTGTCGCAAGAGAAAACTTGACTGTATACAGAATCTTGAATCAGAAATTGAGAAGCTGCAAAGTGAAAAGGAGAGCTTGTTGAAGGAAAGAGATCACATTTTGTCAACTCTGGGTGAGACAAAGCAGAACCTAACTGGA",
         BACH2 = "TTAGAGTTTATTCATGATGTCCGACGGCGCAGCAAGAACCGCATCGCGGCCCAGCGCTGCCGCAAAAGGAAACTGGACTGTATTCAGAATTTAGAATGTGAAATCCGCAAATTGGTGTGTGAGAAAGAGAAACTGTTGTCAGAGAGGAATCAACTGAAAGCATGCATGGGGGAACTGTTGGACAACTTCTCCTGC",
         NFE2L1 = "CTGAGCCTCATCCGAGACATCCGGCGCCGGGGCAAGAACAAGATGGCGGCGCAGAACTGCCGCAAGCGCAAGCTGGACACCATCCTGAATCTGGAGCGTGATGTGGAGGACCTGCAGCGTGACAAAGCCCGGCTGCTGCGGGAGAAAGTGGAGTTCCTGCGCTCCCTGCGACAGATGAAGCAGAAGGTCCAGAGC",
         NFE2 = "CTAGCGCTAGTCCGGGACATCCGACGACGGGGCAAAAACAAGGTGGCAGCCCAGAACTGCCGCAAGAGGAAGCTGGAAACCATTGTGCAGCTGGAGCGGGAGCTGGAGCGGCTGACCAATGAACGGGAGCGGCTTCTCAGGGCCCGCGGGGAGGCAGACCGGACCCTGGAGGTCATGCGCCAACAGCTGACAGAG",
         NFIL3 = "AAGAAAGATGCTATGTATTGGGAAAAAAGGCGGAAAAATAATGAAGCTGCCAAAAGATCTCGTGAGAAGCGTCGACTGAATGACCTGGTTTTAGAGAACAAACTAATTGCACTGGGAGAAGAAAACGCCACTTTAAAAGCTGAGCTGCTTTCACTAAAATTAAAGTTTGGTTTAATTAGCTCCACAGCATATGCT",
         FOS = "GAAGAAGAAGAGAAAAGGAGAATCCGAAGGGAAAGGAATAAGATGGCTGCAGCCAAATGCCGCAACCGGAGGAGGGAGCTGACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTAGAGTTCATCCTGGCAGCT",
         FOSB = "GAGGAAGAGGAGAAGCGAAGGGTGCGCCGGGAACGAAATAAACTAGCAGCAGCTAAATGCAGGAACCGGCGGAGGGAGCTGACCGACCGACTCCAGGCGGAGACAGATCAGTTGGAGGAAGAAAAAGCAGAGCTGGAGTCGGAGATCGCCGAGCTCCAAAAGGAGAAGGAACGTCTGGAGTTTGTGCTGGTGGCC",
         FOSL1 = "GAGGAAGAGGAGCGCCGCCGAGTAAGGCGCGAGCGGAACAAGCTGGCTGCGGCCAAGTGCAGGAACCGGAGGAAGGAACTGACCGACTTCCTGCAGGCGGAGACTGACAAACTGGAAGATGAGAAATCTGGGCTGCAGCGAGAGATTGAGGAGCTGCAGAAGCAGAAGGAGCGCCTAGAGCTGGTGCTGGAAGCC",
         FOSL2 = "GAAGAGGAGGAGAAGCGTCGCATCCGGCGGGAGAGGAACAAGCTGGCTGCAGCCAAGTGCCGGAACCGACGCCGGGAGCTGACAGAGAAGCTGCAGGCGGAGACAGAGGAGCTGGAGGAGGAGAAGTCAGGCCTGCAGAAGGAGATTGCTGAGCTGCAGAAGGAGAAGGAGAAGCTGGAGTTCATGTTGGTGGCT",
         MAFB = "GTGATCCGCCTGAAGCAGAAGCGGCGGACCCTGAAGAACCGGGGCTACGCCCAGTCTTGCAGGTATAAACGCGTCCAGCAGAAGCACCACCTGGAGAATGAGAAGACGCAGCTCATTCAGCAGGTGGAGCAGCTTAAGCAGGAGGTGTCCCGGCTGGCCCGCGAGAGAGACGCCTACAAGGTCAAGTGCGAGAAA",
         JUN = "CAGGAGCGGATCAAGGCGGAGAGGAAGCGCATGAGGAACCGCATCGCTGCCTCCAAGTGCCGAAAAAGGAAGCTGGAGAGAATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTTAAACAGAAAGTCATGAAC",
         JUNB = "CAAGAGCGCATCAAAGTGGAGCGCAAGCGGCTGCGGAACCGGCTGGCGGCCACCAAGTGCCGGAAGCGGAAGCTGGAGCGCATCGCGCGCCTGGAGGACAAGGTGAAGACGCTCAAGGCCGAGAACGCGGGGCTGTCGAGTACCGCCGGCCTCCTCCGGGAGCAGGTGGCCCAGCTCAAACAGAAGGTCATGACC",
         JUND = "CAGGAGCGCATCAAGGCGGAGCGCAAGCGGCTGCGCAACCGCATCGCCGCCTCCAAGTGCCGCAAGCGCAAGCTGGAGCGCATCTCGCGCCTGGAAGAGAAAGTGAAGACCCTCAAGAGTCAGAACACGGAGCTGGCGTCCACGGCGAGCCTGCTGCGCGAGCAGGTGGCGCAGCTCAAGCAGAAAGTCCTCAGC",
         CREB3 = "GAACAAATTCTGAAACGTGTGCGGAGGAAGATTCGAAATAAAAGATCTGCTCAAGAGAGCCGCAGGAAAAAGAAGGTGTATGTTGGGGGTTTAGAGAGCAGGGTCTTGAAATACACAGCCCAGAATATGGAGCTTCAGAACAAAGTACAGCTTCTGGAGGAACAGAATTTGTCCCTTCTAGATCAACTGAGGAAA",
         HLF = "CTGAAGGATGACAAGTACTGGGCAAGGCGCAGAAAGAACAACATGGCAGCCAAGCGCTCCCGCGACGCCCGGAGGCTGAAAGAGAACCAGATCGCCATCCGGGCCTCGTTCCTGGAGAAGGAGAACTCGGCCCTCCGCCAGGAGGTGGCTGACTTGAGGAAGGAGCTGGGCAAATGCAAGAACATACTTGCCAAG",
         MAFG = "ATCGTCCAGCTGAAGCAGCGCCGGCGCACGCTCAAGAACCGCGGCTACGCTGCCAGCTGCCGCGTGAAGCGGGTGACGCAGAAGGAGGAGCTGGAGAAGCAGAAGGCGGAGCTGCAGCAGGAGGTGGAGAAGCTGGCCTCAGAGAACGCCAGCATGAAGCTGGAGCTCGACGCGCTGCGCTCCAAGTACGAGGCG",
         MAFK = "GTGACCCGCCTGAAGCAGCGTCGGCGCACACTCAAGAACCGCGGCTACGCGGCCAGCTGCCGCATCAAGCGGGTGACGCAGAAGGAGGAGCTGGAGCGGCAGCGCGTGGAGCTGCAGCAGGAGGTGGAGAAGCTGGCGCGTGAGAACAGCAGCATGCGGCTGGAGCTGGACGCCCTGCGCTCCAAGTACGAGGCG",
         XBP1 = "AGCCCCGAGGAGAAGGCGCTGAGGAGGAAACTGAAAAACAGAGTAGCAGCTCAGACTGCCAGAGATCGAAAGAAGGCTCGAATGAGTGAGCTGGAACAGCAAGTGGTAGATTTAGAAGAAGAGAACCAAAAACTTTTGCTAGAAAATCAGCTTTTACGAGAGAAAACTCATGGCCTTGTAGTTGAGAACCAGGAG",
         ATF6 = "ATTGCTGTGCTAAGGAGACAGCAACGTATGATAAAAAATCGAGAATCCGCTTGTCAGTCTCGCAAGAAGAAGAAAGAATATATGCTAGGGTTAGAGGCGAGATTAAAGGCTGCCCTCTCAGAAAACGAGCAACTGAAGAAAGAAAATGGAACACTGAAGCGGCAGCTGGATGAAGTTGTGTCAGAGAACCAGAGG",
         ATF6B = "GCAAAGCTGCTGAAGCGGCAGCAGCGAATGATCAAGAACCGGGAGTCAGCCTGCCAGTCCCGGAGAAAGAAGAAAGAGTATCTGCAGGGACTGGAGGCTCGGCTGCAAGCAGTACTGGCTGACAACCAGCAGCTCCGCCGAGAGAATGCTGCCCTCCGGCGGCGGCTGGAGGCCCTGCTGGCTGAAAACAGCGAG",
         CEBPA = "AAGAACAGCAACGAGTACCGGGTGCGGCGCGAGCGCAACAACATCGCGGTGCGCAAGAGCCGCGACAAGGCCAAGCAGCGCAACGTGGAGACGCAGCAGAAGGTGCTGGAGCTGACCAGTGACAATGACCGCCTGCGCAAGCGGGTGGAACAGCTGAGCCGCGAACTGGACACGCTGCGGGGCATCTTCCGCCAG",
         CEBPB = "AAGCACAGCGACGAGTACAAGATCCGGCGCGAGCGCAACAACATCGCCGTGCGCAAGAGCCGCGACAAGGCCAAGATGCGCAACCTGGAGACGCAGCACAAGGTCCTGGAGCTCACGGCCGAGAACGAGCGGCTGCAGAAGAAGGTGGAGCAGCTGTCGCGCGAGCTCAGCACCCTGCGGAACTTGTTCAAGCAG",
         CEBPD = "CGCGGCAGCCCCGAGTACCGGCAGCGGCGCGAGCGCAACAACATCGCCGTGCGCAAGAGCCGCGACAAGGCCAAGCGGCGCAACCAGGAGATGCAGCAGAAGTTGGTGGAGCTGTCGGCTGAGAACGAGAAGCTGCACCAGCGCGTGGAGCAGCTCACGCGGGACCTGGCCGGCCTCCGGCAGTTCTTCAAGCAG",
         CEBPG = "CGAAACAGTGACGAGTATCGGCAACGCCGAGAGAGGAACAACATGGCTGTGAAAAAGAGCCGGTTGAAAAGCAAGCAGAAAGCACAAGACACACTGCAGAGAGTCAATCAGCTCAAAGAAGAGAATGAACGGTTGGAAGCAAAAATCAAATTGCTGACCAAGGAATTAAGTGTACTCAAAGATTTGTTTCTTGAG",
         CREB1 = "GAAGCAGCACGAAAGAGAGAGGTCCGTCTAATGAAGAACAGGGAAGCAGCTCGAGAGTGTCGTAGAAAGAAGAAAGAATATGTGAAATGTTTAGAAAACAGAGTGGCAGTGCTTGAAAATCAAAACAAGACATTGATTGAGGAGCTAAAAGCACTTAAGGACCTTTACTGCCACAAATCAGAT",
         CREB3L1 = "GAGAAGGCCTTGAAGAGAGTCCGGAGGAAAATCAAGAACAAGATCTCAGCCCAGGAGAGCCGTCGTAAGAAGAAGGAGTATGTGGAGTGTCTAGAAAAGAAGGTGGAGACATTTACATCTGAGAACAATGAACTGTGGAAGAAGGTGGAGACCCTGGAGAATGCCAACAGGACCCTGCTCCAGCAGCTGCAGAAA",
         CREB3L2 = "GAGAAGGCCCTGAAGAAAATTCGGAGGAAGATCAAGAATAAGATTTCTGCTCAGGAAAGTAGGAGAAAGAAGAAAGAATACATGGACAGCCTGGAGAAAAAAGTGGAGTCTTGTTCAACTGAGAACTTGGAGCTTCGGAAGAAGGTAGAGGTTCTAGAGAACACTAATAGGACTCTCCTTCAGCAACTCCAGAAG",
         CREB3L3 = "GAGCGAGTGCTGAAAAAAATCCGCCGGAAAATCCGGAACAAGCAGTCGGCGCAAGAAAGCAGGAAGAAGAAGAAGGAATATATCGATGGCCTGGAGACTCGGATGTCAGCTTGCACTGCTCAGAATCAGGAGTTACAGAGGAAAGTCTTGCATCTCGAGAAGCAAAACCTGTCCCTCTTGGAGCAACTGAAGAAA",
         CREB3L4 = "GAGAGGGTCCTCAAGAAGGTCAGGAGGAAAATCCGTAACAAGCAGTCAGCTCAGGACAGTCGGCGGCGGAAGAAGGAGTACATTGATGGGCTGGAGAGCAGGGTGGCAGCCTGTTCTGCACAGAACCAAGAATTACAGAAAAAAGTCCAGGAGCTGGAGAGGCACAACATCTCCTTGGTAGCTCAGCTCCGCCAG",
         CREBL2 = "CCAGCCAAAATTGACTTGAAAGCAAAACTTGAGAGGAGCCGGCAGAGTGCAAGAGAATGCCGAGCCCGAAAAAAGCTGAGATATCAGTATTTGGAAGAGTTGGTATCCAGTCGAGAAAGAGCTATATGTGCCCTCAGAGAGGAACTGGAAATGTACAAGCAGTGGTGCATGGCAATGGACCAAGGAAAAATCCCT",
         CREBRF = "CCCTTAACAGCCCGACCAAGGTCAAGGAAGGAAAAAAATAAGCTGGCTTCCAGAGCTTGTCGGTTAAAGAAGAAAGCCCAGTATGAAGCTAATAAAGTGAAATTATGGGGCCTCAACACAGAATATGATAATTTATTGTTTGTAATCAACTCCATCAAGCAAGAGATTGTAAACCGGGTACAGAATCCAAGAGAT",
         DBP = "CAGAAGGATGAGAAATACTGGAGCCGGCGGTACAAGAACAACGAGGCAGCCAAGCGGTCCCGTGACGCCCGGCGGCTCAAGGAGAACCAGATATCGGTGCGGGCGGCCTTCCTGGAGAAGGAGAACGCCCTGCTGCGGCAGGAAGTTGTGGCCGTGCGCCAGGAGCTGTCCCACTACCGCGCCGTGCTGTCCCGA",
         NFE2L2 = "CTTGCATTAATTCGGGATATACGTAGGAGGGGTAAGAATAAAGTGGCTGCTCAGAATTGCAGAAAAAGAAAACTGGAAAATATAGTAGAACTAGAGCAAGATTTAGATCATTTGAAAGATGAAAAAGAAAAATTGCTCAAAGAAAAAGGAGAAAATGACAAAAGCCTTCACCTACTGAAAAAACAACTCAGCACC",
         NFE2L3 = "GTCTCACTTATCCGTGACATCAGACGAAGAGGGAAAAATAAAGTTGCTGCGCAGAACTGTCGTAAACGCAAATTGGACATAATTTTGAATTTAGAAGATGATGTATGTAACTTGCAAGCAAAGAAGGAAACTCTTAAGAGAGAGCAAGCACAATGTAACAAAGCTATTAACATAATGAAACAGAAACTGCATGAC",
         TEF = "CAGAAGGATGAAAAGTACTGGACAAGACGCAAGAAGAACAACGTGGCAGCTAAACGGTCACGGGATGCCCGGCGCCTGAAAGAGAATCAGATCACCATCCGGGCAGCCTTCCTGGAGAAGGAGAACACAGCCCTGCGGACGGAGGTGGCCGAGCTACGCAAGGAGGTGGGCAAGTGCAAGACCATCGTGTCCAAG"
)
## Generate values for comparison in unit tests - trans experiment
fqt1 <- system.file("extdata/leujunt0_1.fastq.gz", package = "mutscan")
fqt2 <- system.file("extdata/leujunt0_2.fastq.gz", package = "mutscan")
## default arguments
Ldef <- list(
  fastqForward = fqt1, fastqReverse = fqt2, 
  mergeForwardReverse = FALSE, revComplForward = FALSE, revComplReverse = FALSE,
  skipForward = -1, skipReverse = -1, 
  umiLengthForward = -1, umiLengthReverse = -1, 
  constantLengthForward = -1, constantLengthReverse = -1, 
  variableLengthForward = -1, variableLengthReverse = 96,
  adapterForward = "", 
  adapterReverse = "",
  primerForward = "GTCAGGTGGAGGCGGATCC",
  primerReverse = "GAAAAAGGAAGCTGGAGAGA",
  wildTypeForward = leu,
  wildTypeReverse = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT", 
  constantForward = "", 
  constantReverse = "", 
  avePhredMinForward = 20.0, avePhredMinReverse = 20.0,
  variableNMaxForward = 0, variableNMaxReverse = 0, 
  umiNMax = 0,
  nbrMutatedCodonsMaxForward = 0,
  nbrMutatedCodonsMaxReverse = 1,
  forbiddenMutatedCodonsForward = "",
  forbiddenMutatedCodonsReverse = "NNW",
  mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
  mutNameDelimiter = ".", verbose = FALSE
)
processReadsLeuJun(Ldef)

