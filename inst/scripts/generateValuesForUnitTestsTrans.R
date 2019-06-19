suppressPackageStartupMessages({
  library(Biostrings)
  library(ShortRead)
})

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
  uniqueMutCodons <- uniqueMutCodons[!forbiddencodon]
  splitMutCodons <- as(strsplit(gsub("([^.]{3})", "\\1\\.", 
                                     as.character(mutCodons)), ".", fixed = TRUE), 
                       "CharacterList")
  encodedMutCodonsForward <- S4Vectors::unstrsplit(BiocGenerics::relist(
    paste0("f", unlist(uniqueMutCodons), unlist(splitMutCodons)), 
    uniqueMutCodons),
    sep = "_")
  
  
  ## Compare to wildtype, reverse
  if (!is.null(Ldef$fastqReverse)) {
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
    mismatchQualities <- Biostrings::quality(varReverse)[mismatchPositions]
    lowqualmm <- min(as(mismatchQualities, "IntegerList")) < Ldef$mutatedPhredMinReverse
    message("Number of reads with low-quality mismatch bases, reverse: ", sum(lowqualmm))  ## 0
    fq1 <- fq1[!lowqualmm]
    fq2 <- fq2[!lowqualmm]
    mismatchPositions <- mismatchPositions[!lowqualmm]
    varReverse <- varReverse[!lowqualmm]
    encodedMutCodonsForward <- encodedMutCodonsForward[!lowqualmm]
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
    uniqueMutCodons <- uniqueMutCodons[!forbiddencodon]
    splitMutCodons <- as(strsplit(gsub("([^.]{3})", "\\1\\.", 
                                       as.character(mutCodons)), ".", fixed = TRUE), 
                         "CharacterList")
    encodedMutCodonsReverse <- S4Vectors::unstrsplit(BiocGenerics::relist(
      paste0("r", unlist(uniqueMutCodons), unlist(splitMutCodons)), 
      uniqueMutCodons),
      sep = "_")
    mutNames <- paste(encodedMutCodonsForward, encodedMutCodonsReverse, sep = "_")
  } else {
    mutNames <- encodedMutCodonsForward
  }
  mutNames <- gsub("^_", "", gsub("_$", "", mutNames))
  mutNames[mutNames <- ""] <- "WT"
  message("Number of variants seen twice: ", sum(table(mutNames) == 2))  ## 2
  message("Variants seen twice: ")
  print(which(table(mutNames) == 2))
  
  ## Retained reads
  message("Number of retained reads: ", length(fq1))  ## 279
  
  ## Error statistics
  ## Forward
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
  
  ## Reverse
  if (!is.null(Ldef$fastqReverse)) {
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
  verbose = FALSE
)
processReadsTrans(Ldef)
