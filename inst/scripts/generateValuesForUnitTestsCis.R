suppressPackageStartupMessages({
  library(Biostrings)
  library(ShortRead)
})

processReadsCis <- function(Ldef) {
  ## Read data
  fq1 <- Biostrings::readQualityScaledDNAStringSet(Ldef$fastqForward)
  fq2 <- Biostrings::readQualityScaledDNAStringSet(Ldef$fastqReverse)
  
  ## Number of reads
  message("Number of reads: ", length(fq1))  ## 1000
  
  ## Adapter sequences
  adapt <- Biostrings::vcountPattern(Ldef$adapterForward, fq1) > 0 |
    Biostrings::vcountPattern(Ldef$adapterReverse, fq2) > 0
  message("Number of reads with adapters: ", sum(adapt))  ## 126
  fq1 <- fq1[!adapt]
  fq2 <- fq2[!adapt]
  
  ## Merge sequences
  varForward <-  Biostrings::subseq(fq1, start = Ldef$skipForward + Ldef$umiLengthForward + 
                                      Ldef$constantLengthForward + 1, width = Ldef$variableLengthForward)
  varReverse <-  Biostrings::subseq(fq2, start = Ldef$skipReverse + Ldef$umiLengthReverse + 
                                      Ldef$constantLengthReverse + 1, width = Ldef$variableLengthReverse)
  varReverse <- Biostrings::reverseComplement(varReverse)
  replaceAt <- as(ShortRead::FastqQuality(Biostrings::quality(varReverse)), "matrix") > 
    as(ShortRead::FastqQuality(Biostrings::quality(varForward)), "matrix")
  replaceAtList <- as(lapply(seq_len(nrow(replaceAt)),
                             function(x) which(replaceAt[x, ])), "IntegerList")
  consensusSeq <- Biostrings::replaceLetterAt(x = varForward, at = replaceAt, 
                                              letter = varReverse[replaceAtList])
  qualsForward <- as(ShortRead::FastqQuality(Biostrings::quality(varForward)), "matrix")
  qualsReverse <- as(ShortRead::FastqQuality(Biostrings::quality(varReverse)), "matrix")
  qualsForward[replaceAt] <- qualsReverse[replaceAt]
  phredQuals <- apply(qualsForward, 1, function(v) seqTools::ascii2char(v + 33))
  consensusQuals <- Biostrings::PhredQuality(phredQuals)
  varForward <- Biostrings::QualityScaledDNAStringSet(consensusSeq, consensusQuals)
  
  ## Average quality in variable sequences too low
  avgQualForward <- ShortRead::alphabetScore(Biostrings::quality(varForward))/
    BiocGenerics::width(varForward)
  lowqual <- avgQualForward < Ldef$avePhredMinForward
  message("Number of reads with low average quality: ", sum(lowqual))  # 0
  fq1 <- fq1[!lowqual]
  fq2 <- fq2[!lowqual]
  varForward <- varForward[!lowqual]
  
  ## Too many Ns in variable sequences
  nbrNForward <- Biostrings::vcountPattern("N", varForward, fixed = TRUE)
  toomanynvar <- nbrNForward > Ldef$variableNMaxForward
  message("Number of reads with too many N in variable sequence: ", sum(toomanynvar))  # 44
  fq1 <- fq1[!toomanynvar]
  fq2 <- fq2[!toomanynvar]
  varForward <- varForward[!toomanynvar]
  
  ## Too many Ns in UMIs
  umiForward <- Biostrings::subseq(fq1, start = Ldef$skipForward + 1, width = Ldef$umiLengthForward)
  umiReverse <- Biostrings::subseq(fq2, start = Ldef$skipReverse + 1, width = Ldef$umiLengthReverse)
  umis <- Biostrings::xscat(umiForward, umiReverse)
  nbrNumi <- Biostrings::vcountPattern("N", umis, fixed = TRUE)
  toomanynumi <- nbrNumi > Ldef$umiNMax
  message("Number of reads with too many N in UMI: ", sum(toomanynumi))  # 0
  fq1 <- fq1[!toomanynumi]
  fq2 <- fq2[!toomanynumi]
  varForward <- varForward[!toomanynumi]
  
  ## Compare to wildtype, forward
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
  message("Number of reads with low-quality mismatch bases: ", sum(lowqualmm))  ## 0
  fq1 <- fq1[!lowqualmm]
  fq2 <- fq2[!lowqualmm]
  mismatchPositions <- mismatchPositions[!lowqualmm]
  varForward <- varForward[!lowqualmm]
  ## Number of mutated codons/bases
  codonMismatches <- (mismatchPositions - 1) %/% 3 + 1
  nbrMutCodons <- S4Vectors::elementNROWS(unique(codonMismatches))
  baseMismatches <- mismatchPositions
  nbrMutBases <- S4Vectors::elementNROWS(unique(baseMismatches))
  if (Ldef$nbrMutatedCodonsMaxForward != (-1)) {
    toomanycodons <- nbrMutCodons > Ldef$nbrMutatedCodonsMaxForward
  } else {
    toomanycodons <- nbrMutBases > Ldef$nbrMutatedBasesMaxForward
  }
  message("Number of reads with too many mutated codons: ", sum(toomanycodons))  ## 581
  fq1 <- fq1[!toomanycodons]
  fq2 <- fq2[!toomanycodons]
  codonMismatches <- codonMismatches[!toomanycodons]
  baseMismatches <- baseMismatches[!toomanycodons]
  varForward <- varForward[!toomanycodons]
  ## Forbidden codons
  # if (Ldef$nbrMutatedCodonsMaxForward != (-1)) {
  uniqueMutCodons <- unique(codonMismatches)
  uniqueMutBases <- unique(baseMismatches)
  mutBases <- S4Vectors::endoapply(uniqueMutBases,
                                    function(v) v)
  mutBases <- as(varForward, "DNAStringSet")[mutBases]
  mutCodons <- S4Vectors::endoapply(uniqueMutCodons, 
                                    function(v) rep(3 * v, each = 3) + c(-2, -1, 0))
  mutCodons <- as(varForward, "DNAStringSet")[mutCodons]
  matchForbidden <- Biostrings::vmatchPattern(
    pattern = Ldef$forbiddenMutatedCodonsForward, 
    as(mutCodons, "DNAStringSet"), 
    fixed = FALSE)
  forbiddencodon <- sapply(BiocGenerics::relist(start(unlist(matchForbidden)) %% 3==1, 
                                                matchForbidden), function(i) any(i))
  message("Number of reads with forbidden codons: ", sum(forbiddencodon))  ## 82
  fq1 <- fq1[!forbiddencodon]
  fq2 <- fq2[!forbiddencodon]
  varForward <- varForward[!forbiddencodon]
  mutCodons <- mutCodons[!forbiddencodon]
  mutBases <- mutBases[!forbiddencodon]
  uniqueMutCodons <- uniqueMutCodons[!forbiddencodon]
  uniqueMutBases <- uniqueMutBases[!forbiddencodon]
  if (Ldef$nbrMutatedCodonsMaxForward != (-1)) {
    splitMutCodons <- as(strsplit(gsub("([^.]{3})", "\\1\\.", 
                                       as.character(mutCodons)), ".", fixed = TRUE), 
                         "CharacterList")
  } else {
    splitMutCodons <- as(strsplit(gsub("([^.]{1})", "\\1\\.",
                                       as.character(mutBases)), ".", fixed = TRUE),
                         "CharacterList")
    uniqueMutCodons <- uniqueMutBases
    mutCodons <- mutBases
  }
  
  encodedMutCodonsForward <- S4Vectors::unstrsplit(BiocGenerics::relist(
    paste0("f", unlist(uniqueMutCodons), unlist(splitMutCodons)), 
    uniqueMutCodons),
    sep = "_")
  
  mutNames <- encodedMutCodonsForward
  mutNames <- gsub("^_", "", gsub("_$", "", mutNames))
  mutNames[mutNames <- ""] <- "WT"
  message("Number of variants seen thrice: ", sum(table(mutNames) == 3))  ## 11
  message("Variants seen thrice: ")
  print(which(table(mutNames) == 3))  ## f12G f17T f60G f85T f91G f96G
  
  ## Retained reads
  message("Number of retained reads: ", length(fq1))  ## 167
  
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
  constReverse <- Biostrings::subseq(fq2, start = Ldef$skipReverse + Ldef$umiLengthReverse + 1,
                                     width = Ldef$constantLengthReverse)
  constReverse <- Biostrings::reverseComplement(constReverse)
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

## ----------------------------------------------------------------------------
## Generate values for comparison in unit tests - cis experiment
fqc1 <- system.file("extdata/cisInput_1.fastq.gz", package = "mutscan")
fqc2 <- system.file("extdata/cisInput_2.fastq.gz", package = "mutscan")
## default arguments
Ldef <- list(
  fastqForward = fqc1, fastqReverse = fqc2, 
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
processReadsCis(Ldef)

## Limit number of mismatching bases rather than codons
Ldef <- list(
  fastqForward = fqc1, fastqReverse = fqc2, 
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
  nbrMutatedCodonsMaxForward = -1,
  nbrMutatedCodonsMaxReverse = -1,
  nbrMutatedBasesMaxForward = 2,
  nbrMutatedBasesMaxReverse = 2,
  forbiddenMutatedCodonsForward = "NNW",
  forbiddenMutatedCodonsReverse = "NNW",
  mutatedPhredMinForward = 0.0, mutatedPhredMinReverse = 0.0,
  verbose = FALSE
)
processReadsCis(Ldef)
