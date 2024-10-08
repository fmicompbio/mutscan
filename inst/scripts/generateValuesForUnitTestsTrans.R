suppressPackageStartupMessages({
  library(Biostrings)
  library(ShortRead)
  library(dplyr)
})

collapse_seqs <- function(seqs, umis, reads, maxdist, umimaxdist, minreads, minratio) {
  tokeep <- list()
  while (length(seqs) > 0) {
    if (reads[1] >= minreads) {
      d <- stringdist::stringdist(seqs[1], seqs, method = "hamming")
      w <- union(1, which(d <= maxdist & reads[1] >= minratio * reads))
    } else {
      w <- 1
    }
    tokeep[[seqs[1]]] <- paste(umis[w], collapse = ",")
    seqs <- seqs[-w]
    umis <- umis[-w]
    reads <- reads[-w]
  }
  tokeep <- lapply(tokeep, function(tk) {
    tk <- sort(unique(strsplit(tk, ",")[[1]]))
    if (length(tk) == 1) {
      tk
    } else {
      k <- c()
      while (length(tk) > 0) {
        d <- stringdist::stringdist(tk[1], tk, method = "hamming")
        w <- which(d <= umimaxdist)
        k <- c(k, tk[1])
        tk <- tk[-w]
      }
      k
    }
  })
  tokeep
}


processReadsTrans <- function(Ldef) {
  ## Read data
  fq1 <- Biostrings::readQualityScaledDNAStringSet(Ldef$fastqForward)
  if (!is.null(Ldef$fastqReverse)) {
    fq2 <- Biostrings::readQualityScaledDNAStringSet(Ldef$fastqReverse)
  }
  
  ## Number of reads
  message("Number of reads: ", length(fq1))  ## 1000
  
  ## Adapter sequences
  adapt <- Biostrings::vcountPattern(Ldef$adapterForward, fq1) > 0
  if (!is.null(Ldef$fastqReverse)) {
    adapt2 <- Biostrings::vcountPattern(Ldef$adapterReverse, fq2) > 0
    adapt <- adapt | adapt2
  }
  message("Number of reads with adapters: ", sum(adapt))  ## 314
  fq1 <- fq1[!adapt]
  if (!is.null(Ldef$fastqReverse)) {
    fq2 <- fq2[!adapt]
  }
  
  ## Average quality in variable sequences too low
  varForward <-  Biostrings::subseq(fq1, start = Ldef$skipForward + Ldef$umiLengthForward + 
                                      Ldef$constantLengthForward + 1, width = Ldef$variableLengthForward)
  avgQualForward <- ShortRead::alphabetScore(Biostrings::quality(varForward))/
    BiocGenerics::width(varForward)
  lowqual <- avgQualForward < Ldef$avePhredMinForward
  if (!is.null(Ldef$fastqReverse)) {
    varReverse <-  Biostrings::subseq(fq2, start = Ldef$skipReverse + Ldef$umiLengthReverse + 
                                        Ldef$constantLengthReverse + 1, width = Ldef$variableLengthReverse)
    avgQualReverse <- ShortRead::alphabetScore(Biostrings::quality(varReverse))/
      BiocGenerics::width(varReverse)
    lowqual <- lowqual | avgQualReverse < Ldef$avePhredMinReverse
  }
  message("Number of reads with low average quality: ", sum(lowqual))  # 7
  fq1 <- fq1[!lowqual]
  if (!is.null(Ldef$fastqReverse)) {
    fq2 <- fq2[!lowqual]
  }
  
  ## Too many Ns in variable sequences
  varForward <- Biostrings::subseq(fq1, start = Ldef$skipForward + Ldef$umiLengthForward + 
                                     Ldef$constantLengthForward + 1, width = Ldef$variableLengthForward)
  nbrNForward <- Biostrings::vcountPattern("N", varForward, fixed = TRUE)
  toomanynvar <- nbrNForward > Ldef$variableNMaxForward
  if (!is.null(Ldef$fastqReverse)) {
    varReverse <- Biostrings::subseq(fq2, start = Ldef$skipReverse + Ldef$umiLengthReverse + 
                                       Ldef$constantLengthReverse + 1, width = Ldef$variableLengthReverse)
    nbrNReverse <- Biostrings::vcountPattern("N", varReverse, fixed = TRUE)
    toomanynvar <- toomanynvar | nbrNReverse > Ldef$variableNMaxReverse
  }
  message("Number of reads with too many N in variable sequence: ", sum(toomanynvar))  # 0
  fq1 <- fq1[!toomanynvar]
  if (!is.null(Ldef$fastqReverse)) {
    fq2 <- fq2[!toomanynvar]
  }
  
  ## Too many Ns in UMIs
  umiForward <- Biostrings::subseq(fq1, start = Ldef$skipForward + 1, width = Ldef$umiLengthForward)
  if (!is.null(Ldef$fastqReverse)) {
    umiReverse <- Biostrings::subseq(fq2, start = Ldef$skipReverse + 1, width = Ldef$umiLengthReverse)
    umis <- Biostrings::xscat(umiForward, umiReverse)
  } else {
    umis <- umiForward
  }
  nbrNumi <- Biostrings::vcountPattern("N", umis, fixed = TRUE)
  toomanynumi <- nbrNumi > Ldef$umiNMax
  message("Number of reads with too many N in UMI: ", sum(toomanynumi))  # 0
  fq1 <- fq1[!toomanynumi]
  if (!is.null(Ldef$fastqReverse)) {
    fq2 <- fq2[!toomanynumi]
  }
  
  ## Compare to wildtype, forward
  mutCodons <- NULL
  uniqueMutCodons <- NULL
  mutBases <- NULL
  uniqueMutBases <- NULL
  nbrMutBasesTot <- NULL
  if (Ldef$wildTypeForward != "") {
    varForward <- Biostrings::subseq(fq1, start = Ldef$skipForward + Ldef$umiLengthForward + 
                                       Ldef$constantLengthForward + 1, width = Ldef$variableLengthForward)
    p <- DNAStringSet(rep(Ldef$wildTypeForward, length(varForward)))
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
    ## First check, excluding reads with more than 3x(max nbr mut codons) mismatches
    baseMismatches <- mismatchPositions
    nbrMutBases <- S4Vectors::elementNROWS(unique(baseMismatches))
    toomanybases <- nbrMutBases > 3 * Ldef$nbrMutatedCodonsMaxForward
    message("Number of reads with too many mutated codons, forward (1): ", sum(toomanybases))
    fq1 <- fq1[!toomanybases]
    if (!is.null(Ldef$fastqReverse)) {
      fq2 <- fq2[!toomanybases]
    }
    mismatchPositions <- mismatchPositions[!toomanybases]
    varForward <- varForward[!toomanybases]
    
    mismatchQualities <- Biostrings::quality(varForward)[mismatchPositions]
    lowqualmm <- min(as(mismatchQualities, "IntegerList")) < Ldef$mutatedPhredMinForward
    message("Number of reads with low-quality mismatch bases, forward: ", sum(lowqualmm))  ## 0
    fq1 <- fq1[!lowqualmm]
    if (!is.null(Ldef$fastqReverse)) {
      fq2 <- fq2[!lowqualmm]
    }
    mismatchPositions <- mismatchPositions[!lowqualmm]
    varForward <- varForward[!lowqualmm]
    ## Number of mutated codons
    codonMismatches <- (mismatchPositions - 1) %/% 3 + 1
    nbrMutCodons <- S4Vectors::elementNROWS(unique(codonMismatches))
    toomanycodons <- nbrMutCodons > Ldef$nbrMutatedCodonsMaxForward
    message("Number of reads with too many mutated codons, forward: ", sum(toomanycodons))  ## 287
    fq1 <- fq1[!toomanycodons]
    if (!is.null(Ldef$fastqReverse)) {
      fq2 <- fq2[!toomanycodons]
    }
    codonMismatches <- codonMismatches[!toomanycodons]
    varForward <- varForward[!toomanycodons]
    nbrMutCodons <- nbrMutCodons[!toomanycodons]
    nbrMutBases <- S4Vectors::elementNROWS(unique(mismatchPositions[!toomanycodons]))
    ## Forbidden codons
    uniqueMutCodons <- unique(codonMismatches)
    mutCodons <- S4Vectors::endoapply(uniqueMutCodons, 
                                      function(v) rep(3 * v, each = 3) + c(-2, -1, 0))
    mutCodons <- as(varForward, "DNAStringSet")[mutCodons]
    matchForbidden <- Biostrings::vmatchPattern(
      pattern = Ldef$forbiddenMutatedCodonsForward, 
      as(mutCodons, "DNAStringSet"), 
      fixed = FALSE)
    forbiddencodon <- sapply(BiocGenerics::relist(start(unlist(matchForbidden)) %% 3==1, 
                                                  matchForbidden), function(i) any(i))
    message("Number of reads with forbidden codons, forward: ", sum(forbiddencodon))  ## 6
    fq1 <- fq1[!forbiddencodon]
    if (!is.null(Ldef$fastqReverse)) {
      fq2 <- fq2[!forbiddencodon]
    }
    mutCodons <- mutCodons[!forbiddencodon]
    nbrMutCodonsTot <- nbrMutCodons[!forbiddencodon]
    nbrMutBasesTot <- nbrMutBases[!forbiddencodon]
    uniqueMutCodons <- uniqueMutCodons[!forbiddencodon]
    splitMutCodons <- as(strsplit(gsub("([^.]{3})", "\\1\\.", 
                                       as.character(mutCodons)), ".", fixed = TRUE), 
                         "CharacterList")
    encodedMutCodonsForward <- S4Vectors::unstrsplit(BiocGenerics::relist(
      paste0("f", unlist(uniqueMutCodons), unlist(splitMutCodons)), 
      uniqueMutCodons),
      sep = "_")
  } else {
    encodedMutCodonsForward <- ""
  }
  
  ## Compare to wildtype, reverse
  if (!is.null(Ldef$fastqReverse) && Ldef$wildTypeReverse != "") {
    varReverse <- Biostrings::subseq(fq2, start = Ldef$skipReverse + Ldef$umiLengthReverse + 
                                       Ldef$constantLengthReverse + 1, width = Ldef$variableLengthReverse)
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
    
    ## First check, excluding reads with more than 3x(max nbr mut codons) mismatches
    baseMismatches <- mismatchPositions
    nbrMutBases <- S4Vectors::elementNROWS(unique(baseMismatches))
    toomanybases <- nbrMutBases > 3 * Ldef$nbrMutatedCodonsMaxReverse
    message("Number of reads with too many mutated codons, reverse (1): ", sum(toomanybases))
    fq1 <- fq1[!toomanybases]
    if (!is.null(Ldef$fastqReverse)) {
      fq2 <- fq2[!toomanybases]
    }
    mismatchPositions <- mismatchPositions[!toomanybases]
    varReverse <- varReverse[!toomanybases]
    encodedMutCodonsForward <- encodedMutCodonsForward[!toomanybases]
    nbrMutBasesTot <- nbrMutBasesTot[!toomanybases]
    nbrMutCodonsTot <- nbrMutCodonsTot[!toomanybases]
    
    mismatchQualities <- Biostrings::quality(varReverse)[mismatchPositions]
    lowqualmm <- min(as(mismatchQualities, "IntegerList")) < Ldef$mutatedPhredMinReverse
    message("Number of reads with low-quality mismatch bases, reverse: ", sum(lowqualmm))  ## 0
    fq1 <- fq1[!lowqualmm]
    fq2 <- fq2[!lowqualmm]
    mismatchPositions <- mismatchPositions[!lowqualmm]
    varReverse <- varReverse[!lowqualmm]
    encodedMutCodonsForward <- encodedMutCodonsForward[!lowqualmm]
    nbrMutBasesTot <- nbrMutBasesTot[!lowqualmm]
    nbrMutCodonsTot <- nbrMutCodonsTot[!lowqualmm]
    ## Number of mutated codons
    codonMismatches <- (mismatchPositions - 1) %/% 3 + 1
    nbrMutCodons <- S4Vectors::elementNROWS(unique(codonMismatches))
    toomanycodons <- nbrMutCodons > Ldef$nbrMutatedCodonsMaxReverse
    message("Number of reads with too many mutated codons, reverse: ", sum(toomanycodons))  ## 105
    fq1 <- fq1[!toomanycodons]
    fq2 <- fq2[!toomanycodons]
    codonMismatches <- codonMismatches[!toomanycodons]
    varReverse <- varReverse[!toomanycodons]
    encodedMutCodonsForward <- encodedMutCodonsForward[!toomanycodons]
    nbrMutBasesTot <- nbrMutBasesTot[!toomanycodons]
    nbrMutCodonsTot <- nbrMutCodonsTot[!toomanycodons]
    nbrMutCodons <- nbrMutCodons[!toomanycodons]
    nbrMutBases <- S4Vectors::elementNROWS(unique(mismatchPositions[!toomanycodons]))
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
    message("Number of reads with forbidden codons, reverse: ", sum(forbiddencodon))  ## 2
    fq1 <- fq1[!forbiddencodon]
    fq2 <- fq2[!forbiddencodon]
    encodedMutCodonsForward <- encodedMutCodonsForward[!forbiddencodon]
    mutCodons <- mutCodons[!forbiddencodon]
    nbrMutCodonsTot <- nbrMutCodonsTot[!forbiddencodon]
    nbrMutBasesTot <- nbrMutBasesTot[!forbiddencodon]
    nbrMutCodons <- nbrMutCodons[!forbiddencodon]
    nbrMutBases <- nbrMutBases[!forbiddencodon]
    uniqueMutCodons <- uniqueMutCodons[!forbiddencodon]
    splitMutCodons <- as(strsplit(gsub("([^.]{3})", "\\1\\.", 
                                       as.character(mutCodons)), ".", fixed = TRUE), 
                         "CharacterList")
    encodedMutCodonsReverse <- S4Vectors::unstrsplit(BiocGenerics::relist(
      paste0("r", unlist(uniqueMutCodons), unlist(splitMutCodons)), 
      uniqueMutCodons),
      sep = "_")
    mutNames <- paste(encodedMutCodonsForward, encodedMutCodonsReverse, sep = "_")
    
    nbrMutCodonsTot <- nbrMutCodonsTot + nbrMutCodons
    nbrMutBasesTot <- nbrMutBasesTot + nbrMutBases
    
  } else {
    mutNames <- encodedMutCodonsForward
  }
  
  ## Constant sequence filtering
  keep <- rep(TRUE, length(fq1))
  if (!all(Ldef$constantForward == "") && Ldef$constantMaxDistForward > (-1)) {
    constForward <- Biostrings::subseq(fq1, start = Ldef$skipForward + Ldef$umiLengthForward + 1, 
                                       width = Ldef$constantLengthForward)
    dists <- as.matrix(pwalign::stringDist(
      c(DNAStringSet(structure(Ldef$constantForward, 
                               names = paste0("C", length(Ldef$constantForward)))),
        DNAStringSet(constForward)
      )
    ))
    dists <- dists[-(1:length(Ldef$constantForward)), 1:length(Ldef$constantForward), drop = FALSE]
    keep <- keep & (apply(dists, 1, min) <= Ldef$constantMaxDistForward)
  }
  if (!all(Ldef$constantReverse == "") && Ldef$constantMaxDistReverse > (-1)) {
    constReverse <- Biostrings::subseq(fq2, start = Ldef$skipReverse + Ldef$umiLengthReverse + 1, 
                                       width = Ldef$constantLengthReverse)
    dists <- as.matrix(pwalign::stringDist(
      c(DNAStringSet(structure(Ldef$constantReverse, 
                               names = paste0("C", length(Ldef$constantReverse)))),
        DNAStringSet(constReverse)
      )
    ))
    dists <- dists[-(1:length(Ldef$constantReverse)), 1:length(Ldef$constantReverse), drop = FALSE]
    keep <- keep & (apply(dists, 1, min) <= Ldef$constantMaxDistReverse)
  }
  message("Number of reads with too many mismatches in constant seq: ", sum(!keep))
  fq1 <- fq1[keep]
  if (!is.null(Ldef$fastqReverse)) {
    fq2 <- fq2[keep]
  }
  mutCodons <- mutCodons[keep]
  uniqueMutCodons <- uniqueMutCodons[keep]
  mutNames <- mutNames[keep]

  if (!is.null(nbrMutBasesTot)) {
    nbrMutBasesTot <- nbrMutBasesTot[keep]
    nbrMutCodonsTot <- nbrMutCodonsTot[keep]
    ## Number of mutated bases and codons
    varForwardFinal <- Biostrings::subseq(
      fq1, start = Ldef$skipForward + Ldef$umiLengthForward + 
        Ldef$constantLengthForward + 1, width = Ldef$variableLengthForward)
    if (!is.null(Ldef$fastqReverse)) {
      varReverseFinal <- Biostrings::subseq(
        fq2, start = Ldef$skipReverse + Ldef$umiLengthReverse + 
          Ldef$constantLengthReverse + 1, width = Ldef$variableLengthReverse)
    } else {
      varReverseFinal <- ""
    }
    message("Number of mutated codons per retained read: ")
    print(table(nbrMutCodonsTot))
    message("Number of different reads with each number of mismatched codons: ")
    print(table(nbrMutCodonsTot[!duplicated(paste0(varForwardFinal, "_", varReverseFinal))]))
    message("Number of mutated bases per retained read: ")
    print(table(nbrMutBasesTot))
    message("Number of different reads with each number of mismatched bases: ")
    print(table(nbrMutBasesTot[!duplicated(paste0(varForwardFinal, "_", varReverseFinal))]))
    
    ## Number of mutated amino acids
    ## Forward
    p_aa <- Biostrings::translate(DNAStringSet(rep(Ldef$wildTypeForward, length(varForwardFinal))))
    pattern_width_aa <- width(p_aa)
    varForwardFinal_aa <- Biostrings::translate(varForwardFinal)
    subject_width_aa <- width(varForwardFinal_aa)
    unlisted_ans_aa <- which(as.raw(unlist(p_aa)) != as.raw(unlist(varForwardFinal_aa)))
    breakpoints_aa <- cumsum(pattern_width_aa)
    ans_eltlens_aa <- tabulate(findInterval(unlisted_ans_aa - 1L,
                                            breakpoints_aa) + 1L,
                               nbins = length(p_aa))
    skeleton_aa <- IRanges::PartitioningByEnd(cumsum(ans_eltlens_aa))
    aaMismatchesForward <- BiocGenerics::relist(unlisted_ans_aa, skeleton_aa)
    ## Reverse
    if (!is.null(Ldef$fastqReverse)) {
      p_aa <- Biostrings::translate(DNAStringSet(rep(Ldef$wildTypeReverse, length(varReverseFinal))))
      pattern_width_aa <- width(p_aa)
      varReverseFinal_aa <- Biostrings::translate(varReverseFinal)
      subject_width_aa <- width(varReverseFinal_aa)
      unlisted_ans_aa <- which(as.raw(unlist(p_aa)) != as.raw(unlist(varReverseFinal_aa)))
      breakpoints_aa <- cumsum(pattern_width_aa)
      ans_eltlens_aa <- tabulate(findInterval(unlisted_ans_aa - 1L,
                                              breakpoints_aa) + 1L,
                                 nbins = length(p_aa))
      skeleton_aa <- IRanges::PartitioningByEnd(cumsum(ans_eltlens_aa))
      aaMismatchesReverse <- BiocGenerics::relist(unlisted_ans_aa, skeleton_aa)
     
      message("Number of mutated aas per retained read: ")
      print(table(lengths(aaMismatchesForward) + lengths(aaMismatchesReverse)))
      message("Number of different reads with each number of mismatched aas: ")
      print(table((lengths(aaMismatchesForward) + lengths(aaMismatchesReverse))[!duplicated(paste0(varForwardFinal, "_", varReverseFinal))]))
    } else {
      message("Number of mutated aas per retained read: ")
      print(table(lengths(aaMismatchesForward)))
      message("Number of different reads with each number of mismatched aas: ")
      print(table(lengths(aaMismatchesForward)[!duplicated(varForwardFinal)])) 
      
    }
  }

  mutNames <- gsub("^_", "", gsub("_$", "", mutNames))
  mutNames[mutNames == ""] <- "WT"
  message("Number of variants seen twice: ", sum(table(mutNames) == 2))  ## 2
  message("Variants seen twice: ")
  print(which(table(mutNames) == 2))
  
  ## Retained reads
  message("Number of retained reads: ", length(fq1))  ## 279
  
  ## Error statistics
  ## Forward
  if (length(Ldef$constantForward) == 1) {
    constForward <- Biostrings::subseq(fq1, start = Ldef$skipForward + Ldef$umiLengthForward + 1, 
                                       width = Ldef$constantLengthForward)
    p <- DNAStringSet(rep(Ldef$constantForward, length(constForward)))
    pattern_width <- width(p)
    subject_width <- width(constForward)
    
    unlisted_ans <- which(as.raw(unlist(p)) != as.raw(unlist(constForward)))
    breakpoints <- cumsum(pattern_width)
    ans_eltlens <- tabulate(findInterval(unlisted_ans - 1L,
                                         breakpoints) + 1L,
                            nbins = length(p))
    skeleton <- IRanges::PartitioningByEnd(cumsum(ans_eltlens))
    nucleotideMismatches <- BiocGenerics::relist(unlisted_ans, skeleton)
    offsets <- c(0L, breakpoints[-length(breakpoints)])
    mismatchPositions <- nucleotideMismatches - offsets
    mismatchQualities <- unlist(as(quality(constForward)[mismatchPositions], "IntegerList"))
    
    unlisted_ans <- which(as.raw(unlist(p)) == as.raw(unlist(constForward)))
    breakpoints <- cumsum(pattern_width)
    ans_eltlens <- tabulate(findInterval(unlisted_ans - 1L,
                                         breakpoints) + 1L,
                            nbins = length(p))
    skeleton <- IRanges::PartitioningByEnd(cumsum(ans_eltlens))
    nucleotideMatches <- BiocGenerics::relist(unlisted_ans, skeleton)
    offsets <- c(0L, breakpoints[-length(breakpoints)])
    matchPositions <- nucleotideMatches - offsets
    matchQualities <- unlist(as(quality(constForward)[matchPositions], "IntegerList"))
    
    message("Match qualities, forward:")
    print(table(matchQualities))
    message("Mismatch qualities, forward:")
    print(table(mismatchQualities))
  }
  
  ## Reverse
  if (!is.null(Ldef$fastqReverse) && length(Ldef$constantReverse) == 1) {
    constReverse <- Biostrings::subseq(fq2, start = Ldef$skipReverse + Ldef$umiLengthReverse + 1,
                                       width = Ldef$constantLengthReverse)
    #constReverse <- Biostrings::reverseComplement(constReverse)
    p <- DNAStringSet(rep(Ldef$constantReverse, length(constReverse)))
    pattern_width <- width(p)
    subject_width <- width(constReverse)
    
    unlisted_ans <- which(as.raw(unlist(p)) != as.raw(unlist(constReverse)))
    breakpoints <- cumsum(pattern_width)
    ans_eltlens <- tabulate(findInterval(unlisted_ans - 1L,
                                         breakpoints) + 1L,
                            nbins = length(p))
    skeleton <- IRanges::PartitioningByEnd(cumsum(ans_eltlens))
    nucleotideMismatches <- BiocGenerics::relist(unlisted_ans, skeleton)
    offsets <- c(0L, breakpoints[-length(breakpoints)])
    mismatchPositions <- nucleotideMismatches - offsets
    mismatchQualities <- unlist(as(quality(constReverse)[mismatchPositions], "IntegerList"))
    
    unlisted_ans <- which(as.raw(unlist(p)) == as.raw(unlist(constReverse)))
    breakpoints <- cumsum(pattern_width)
    ans_eltlens <- tabulate(findInterval(unlisted_ans - 1L,
                                         breakpoints) + 1L,
                            nbins = length(p))
    skeleton <- IRanges::PartitioningByEnd(cumsum(ans_eltlens))
    nucleotideMatches <- BiocGenerics::relist(unlisted_ans, skeleton)
    offsets <- c(0L, breakpoints[-length(breakpoints)])
    matchPositions <- nucleotideMatches - offsets
    matchQualities <- unlist(as(quality(constReverse)[matchPositions], "IntegerList"))
    
    message("Match qualities, reverse:")
    print(table(matchQualities))
    message("Mismatch qualities, reverse:")
    print(table(mismatchQualities))
  }
  
  ## Collapsing
  if (Ldef$wildTypeForward == "" && Ldef$wildTypeReverse == "") {
    umicoll <- list()
    if (!is.null(Ldef$fastqReverse)) {
      fqseq <- paste0(as.character(varForward), "_", as.character(varReverse))
    } else {
      fqseq <- as.character(varForward)
    }
    tbl <- as.data.frame(table(fqseq)) %>% dplyr::arrange(desc(Freq), fqseq) %>%
      dplyr::mutate(fqseq = as.character(fqseq))
    tbl$umis <- ""
    for (i in seq_len(nrow(tbl))) {
      tbl$umis[i] <- paste(umis[fqseq == tbl$fqseq[i]], collapse = ",")
    }
    collseqs <- collapse_seqs(tbl$fqseq, tbl$umis, tbl$Freq, 
                              maxdist = Ldef$variableCollapseMaxDist,
                              umimaxdist = Ldef$umiCollapseMaxDist,
                              minreads = Ldef$variableCollapseMinReads, 
                              minratio = Ldef$variableCollapseMinRatio)
    message("Number of unique sequences: ", nrow(tbl))
    message("Number of collapsed sequences: ", length(collseqs))
    message("Total UMI count: ", length(unlist(collseqs)))
    return(tbl)
  }
  return(NULL)
}

## ----------------------------------------------------------------------------
## Generate values for comparison in unit tests - trans experiment
fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
## default arguments
Ldef <- list(
  fastqForward = fqt1, fastqReverse = fqt2, 
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
  variableCollapseMaxDist = 0, umiCollapseMaxDist = 0, 
  variableCollapseMinReads = 0, variableCollapseMinRatio = 0,
  constantMaxDistForward = -1, constantMaxDistReverse = -1,
  verbose = FALSE
)
processReadsTrans(Ldef)

## ----------------------------------------------------------------------------
## Change some parameter values
fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
Ldef <- list(
  fastqForward = fqt1, fastqReverse = fqt2, 
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
  avePhredMinForward = 20.0, avePhredMinReverse = 30.0,
  variableNMaxForward = 0, variableNMaxReverse = 0, 
  umiNMax = 0,
  nbrMutatedCodonsMaxForward = 1,
  nbrMutatedCodonsMaxReverse = 1,
  forbiddenMutatedCodonsForward = "NNA",
  forbiddenMutatedCodonsReverse = "NNW",
  mutatedPhredMinForward = 25.0, mutatedPhredMinReverse = 0.0,
  variableCollapseMaxDist = 0, umiCollapseMaxDist = 0, 
  variableCollapseMinReads = 0, variableCollapseMinRatio = 0,
  constantMaxDistForward = -1, constantMaxDistReverse = -1,
  verbose = FALSE
)
processReadsTrans(Ldef)

## ----------------------------------------------------------------------------
## Only one input file
fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
Ldef <- list(
  fastqForward = fqt1, fastqReverse = NULL, 
  mergeForwardReverse = FALSE, revComplForward = FALSE, revComplReverse = FALSE,
  skipForward = 1, skipReverse = 1, 
  umiLengthForward = 10, umiLengthReverse = 8, 
  constantLengthForward = 18, constantLengthReverse = 20, 
  variableLengthForward = 96, variableLengthReverse = 96,
  adapterForward = "GGAAGAGCACACGTC", 
  adapterReverse = "GGAAGAGCGTCGTGT",
  wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
  wildTypeReverse = "", 
  constantForward = "AACCGGAGGAGGGAGCTG", 
  constantReverse = "GAAAAAGGAAGCTGGAGAGA", 
  avePhredMinForward = 20.0, avePhredMinReverse = 30.0,
  variableNMaxForward = 0, variableNMaxReverse = 0, 
  umiNMax = 0,
  nbrMutatedCodonsMaxForward = 1,
  nbrMutatedCodonsMaxReverse = 1,
  forbiddenMutatedCodonsForward = "NNA",
  forbiddenMutatedCodonsReverse = "NNW",
  mutatedPhredMinForward = 25.0, mutatedPhredMinReverse = 0.0,
  variableCollapseMaxDist = 0, umiCollapseMaxDist = 0, 
  variableCollapseMinReads = 0, variableCollapseMinRatio = 0,
  constantMaxDistForward = -1, constantMaxDistReverse = -1,
  verbose = FALSE
)
processReadsTrans(Ldef)

## Test collapsing
fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
## default arguments
Ldef <- list(
  fastqForward = fqt1, fastqReverse = fqt2, 
  mergeForwardReverse = FALSE, 
  minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
  revComplForward = FALSE, revComplReverse = FALSE,
  skipForward = 1, skipReverse = 1, 
  umiLengthForward = 10, umiLengthReverse = 8, 
  constantLengthForward = 18, constantLengthReverse = 20, 
  variableLengthForward = 96, variableLengthReverse = 96,
  adapterForward = "GGAAGAGCACACGTC", 
  adapterReverse = "GGAAGAGCGTCGTGT",
  primerForward = "",
  primerReverse = "",
  wildTypeForward = "",
  wildTypeReverse = "", 
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
  mutNameDelimiter = ".",
  maxNReads = -1, variableCollapseMaxDist = 6, umiCollapseMaxDist = 4, 
  variableCollapseMinReads = 0, variableCollapseMinRatio = 0,
  constantMaxDistForward = -1, constantMaxDistReverse = -1,
  verbose = FALSE
)
processReadsTrans(Ldef)

fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
## default arguments
Ldef <- list(
  fastqForward = fqt1, fastqReverse = NULL, 
  mergeForwardReverse = FALSE, 
  minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
  revComplForward = FALSE, revComplReverse = FALSE,
  skipForward = 1, skipReverse = 1, 
  umiLengthForward = 10, umiLengthReverse = 8, 
  constantLengthForward = 18, constantLengthReverse = 20, 
  variableLengthForward = 96, variableLengthReverse = 96,
  adapterForward = "GGAAGAGCACACGTC", 
  adapterReverse = "GGAAGAGCGTCGTGT",
  primerForward = "",
  primerReverse = "",
  wildTypeForward = "",
  wildTypeReverse = "", 
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
  mutNameDelimiter = ".",
  variableCollapseMaxDist = 10, umiCollapseMaxDist = 5, 
  variableCollapseMinReads = 0, variableCollapseMinRatio = 0,
  constantMaxDistForward = -1, constantMaxDistReverse = -1,
  maxNReads = -1, verbose = FALSE
)
processReadsTrans(Ldef)

## don't add anything between the above and below runs
## same as above but with new collapsing approach
Ldef$variableCollapseMaxDist <- 0
Ldef$umiCollapseMaxDist <- 5 
Ldef$variableCollapseMinReads <- 0
Ldef$variableCollapseMinRatio <- 0
v <- processReadsTrans(Ldef)
cseq <- collapse_seqs(v$fqseq, v$umis, v$Freq, maxdist = 10, 
                      umimaxdist = 0, minreads = 0, minratio = 0)
length(cseq)

## Collapse only highly abundant features
fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
## default arguments
Ldef <- list(
  fastqForward = fqt1, fastqReverse = fqt2, 
  mergeForwardReverse = FALSE, 
  minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
  revComplForward = FALSE, revComplReverse = FALSE,
  skipForward = 1, skipReverse = 1, 
  umiLengthForward = 10, umiLengthReverse = 8, 
  constantLengthForward = 18, constantLengthReverse = 20, 
  variableLengthForward = 96, variableLengthReverse = 96,
  adapterForward = "GGAAGAGCACACGTC", 
  adapterReverse = "GGAAGAGCGTCGTGT",
  primerForward = "",
  primerReverse = "",
  wildTypeForward = "",
  wildTypeReverse = "", 
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
  mutNameDelimiter = ".",
  maxNReads = -1, variableCollapseMaxDist = 2, umiCollapseMaxDist = 0, 
  variableCollapseMinReads = 2, variableCollapseMinRatio = 0,
  constantMaxDistForward = -1, constantMaxDistReverse = -1,
  verbose = FALSE
)
v <- processReadsTrans(Ldef)

## don't add anything between the above and below runs
## same as above but with new collapsing approach
Ldef$variableCollapseMaxDist <- 0
Ldef$umiCollapseMaxDist <- 0 
Ldef$variableCollapseMinReads <- 0
Ldef$variableCollapseMinRatio <- 0
v <- processReadsTrans(Ldef)
cseq <- collapse_seqs(v$fqseq, v$umis, v$Freq, maxdist = 2, 
                      umimaxdist = 0, minreads = 1.5, minratio = 0)
length(cseq)

## Collapse only features with high enough count ratio
fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
## default arguments
Ldef <- list(
  fastqForward = fqt1, fastqReverse = fqt2, 
  mergeForwardReverse = FALSE, 
  minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
  revComplForward = FALSE, revComplReverse = FALSE,
  skipForward = 1, skipReverse = 1, 
  umiLengthForward = 10, umiLengthReverse = 8, 
  constantLengthForward = 18, constantLengthReverse = 20, 
  variableLengthForward = 96, variableLengthReverse = 96,
  adapterForward = "GGAAGAGCACACGTC", 
  adapterReverse = "GGAAGAGCGTCGTGT",
  primerForward = "",
  primerReverse = "",
  wildTypeForward = "",
  wildTypeReverse = "", 
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
  mutNameDelimiter = ".",
  maxNReads = -1, variableCollapseMaxDist = 3, umiCollapseMaxDist = 0, 
  variableCollapseMinReads = 1, variableCollapseMinRatio = 1.5,
  constantMaxDistForward = -1, constantMaxDistReverse = -1,
  verbose = FALSE
)
v <- processReadsTrans(Ldef)

## don't add anything between the above and below runs
## same as above but with new collapsing approach
Ldef$variableCollapseMaxDist <- 0
Ldef$umiCollapseMaxDist <- 0 
Ldef$variableCollapseMinReads <- 0
Ldef$variableCollapseMinRatio <- 0
v <- processReadsTrans(Ldef)
cseq <- collapse_seqs(v$fqseq, v$umis, v$Freq, maxdist = 3, 
                      umimaxdist = 0, minreads = 1, minratio = 1.5)
length(cseq)

## filter based on constant sequences
Ldef <- list(
  fastqForward = fqt1, fastqReverse = fqt2, 
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
  variableCollapseMaxDist = 0, umiCollapseMaxDist = 0, 
  variableCollapseMinReads = 0, variableCollapseMinRatio = 0,
  constantMaxDistForward = 0, constantMaxDistReverse = 0,
  verbose = FALSE
)
processReadsTrans(Ldef)

## 1 mismatch
Ldef <- list(
  fastqForward = fqt1, fastqReverse = fqt2, 
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
  variableCollapseMaxDist = 0, umiCollapseMaxDist = 0, 
  variableCollapseMinReads = 0, variableCollapseMinRatio = 0,
  constantMaxDistForward = 1, constantMaxDistReverse = 1,
  verbose = FALSE
)
processReadsTrans(Ldef)

## Multiple constant sequences
Ldef <- list(
  fastqForward = fqt1, fastqReverse = fqt2, 
  mergeForwardReverse = FALSE, revComplForward = FALSE, revComplReverse = FALSE,
  skipForward = 1, skipReverse = 1, 
  umiLengthForward = 10, umiLengthReverse = 8, 
  constantLengthForward = 18, constantLengthReverse = 20, 
  variableLengthForward = 96, variableLengthReverse = 96,
  adapterForward = "GGAAGAGCACACGTC", 
  adapterReverse = "GGAAGAGCGTCGTGT",
  wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
  wildTypeReverse = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT", 
  constantForward = c("AACCGGAGGAGGGAGCTG", "AACCGGCGGAGGGAGCTG"), 
  constantReverse = "GAAAAAGGAAGCTGGAGAGA", 
  avePhredMinForward = 20.0, avePhredMinReverse = 20.0,
  variableNMaxForward = 0, variableNMaxReverse = 0, 
  umiNMax = 0,
  nbrMutatedCodonsMaxForward = 1,
  nbrMutatedCodonsMaxReverse = 1,
  forbiddenMutatedCodonsForward = "NNW",
  forbiddenMutatedCodonsReverse = "NNW",
  mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
  variableCollapseMaxDist = 0, umiCollapseMaxDist = 0, 
  variableCollapseMinReads = 0, variableCollapseMinRatio = 0,
  constantMaxDistForward = 0, constantMaxDistReverse = 0,
  verbose = FALSE
)
processReadsTrans(Ldef)
