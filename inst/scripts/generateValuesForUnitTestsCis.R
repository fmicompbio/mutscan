## Generate values for comparison in unit tests - cis experiment
fqc1 <- system.file("extdata/cisInput_1.fastq.gz", package = "mutscan")
fqc2 <- system.file("extdata/cisInput_2.fastq.gz", package = "mutscan")
## default arguments
Ldef <- list(
  experimentType = "cis", fastqForward = fqc1, 
  fastqReverse = fqc2, skipForward = 1, skipReverse = 1, 
  umiLengthForward = 10, umiLengthReverse = 7, 
  constantLengthForward = 18, constantLengthReverse = 17, 
  variableLengthForward = 96, variableLengthReverse = 96,
  adapterForward = "GGAAGAGCACACGTC", 
  adapterReverse = "GGAAGAGCGTCGTGT",
  wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
  wildTypeReverse = "", 
  constantForward = "AACCGGAGGAGGGAGCTG", 
  constantReverse = "GAGTTCATCCTGGCAGC", 
  avePhredMin = 20.0, variableNMax = 0, umiNMax = 0,
  nbrMutatedCodonsMax = 1,
  forbiddenMutatedCodons = "NNW",
  mutatedPhredMin = 0.0,
  verbose = FALSE
)

## Read data
fq1 <- Biostrings::readQualityScaledDNAStringSet(fqc1)
fq2 <- Biostrings::readQualityScaledDNAStringSet(fqc2)

## Number of reads
length(fq1)  ## 1000

## Adapter sequences
adapt <- Biostrings::vcountPattern(Ldef$adapterForward, fq1) > 0 |
  Biostrings::vcountPattern(Ldef$adapterReverse, fq2) > 0
sum(adapt)  ## 126
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
varForward <- QualityScaledDNAStringSet(consensusSeq, consensusQuals)

## Average quality in variable sequences too low
avgQualForward <- ShortRead::alphabetScore(Biostrings::quality(varForward))/
  BiocGenerics::width(varForward)
lowqual <- avgQualForward < Ldef$avePhredMin
sum(lowqual)  # 0
fq1 <- fq1[!lowqual]
fq2 <- fq2[!lowqual]
varForward <- varForward[!lowqual]

## Too many Ns in variable sequences
nbrNForward <- Biostrings::vcountPattern("N", varForward, fixed = TRUE)
toomanynvar <- nbrNForward > Ldef$variableNMax
sum(toomanynvar)  # 44
fq1 <- fq1[!toomanynvar]
fq2 <- fq2[!toomanynvar]
varForward <- varForward[!toomanynvar]

## Too many Ns in UMIs
umiForward <- Biostrings::subseq(fq1, start = Ldef$skipForward + 1, width = Ldef$umiLengthForward)
umiReverse <- Biostrings::subseq(fq2, start = Ldef$skipReverse + 1, width = Ldef$umiLengthReverse)
umis <- Biostrings::xscat(umiForward, umiReverse)
nbrNumi <- Biostrings::vcountPattern("N", umis, fixed = TRUE)
toomanynumi <- nbrNumi > Ldef$umiNMax
sum(toomanynumi)  # 0
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
nucleotideMismatches <- relist(unlisted_ans, skeleton)
offsets <- c(0L, breakpoints[-length(breakpoints)])
mismatchPositions <- nucleotideMismatches - offsets
mismatchQualities <- Biostrings::quality(varForward)[mismatchPositions]
lowqualmm <- min(as(mismatchQualities, "IntegerList")) < Ldef$mutatedPhredMin
sum(lowqualmm)  ## 0
fq1 <- fq1[!lowqualmm]
fq2 <- fq2[!lowqualmm]
mismatchPositions <- mismatchPositions[!lowqualmm]
varForward <- varForward[!lowqualmm]
## Number of mutated codons
codonMismatches <- (mismatchPositions - 1) %/% 3 + 1
nbrMutCodons <- S4Vectors::elementNROWS(unique(codonMismatches))
toomanycodons <- nbrMutCodons > Ldef$nbrMutatedCodonsMax
sum(toomanycodons)  ## 581
fq1 <- fq1[!toomanycodons]
fq2 <- fq2[!toomanycodons]
codonMismatches <- codonMismatches[!toomanycodons]
varForward <- varForward[!toomanycodons]
## Forbidden codons
uniqueMutCodons <- unique(codonMismatches)
mutCodons <- S4Vectors::endoapply(uniqueMutCodons, 
                                  function(v) rep(3 * v, each = 3) + c(-2, -1, 0))
mutCodons <- as(varForward, "DNAStringSet")[mutCodons]
matchForbidden <- Biostrings::vmatchPattern(
  pattern = Ldef$forbiddenMutatedCodon, 
  as(mutCodons, "DNAStringSet"), 
  fixed = FALSE)
forbiddencodon <- sapply(relist(start(unlist(matchForbidden)) %% 3==1, 
                                matchForbidden), function(i) any(i))
sum(forbiddencodon)  ## 82
fq1 <- fq1[!forbiddencodon]
fq2 <- fq2[!forbiddencodon]
varForward <- varForward[!forbiddencodon]
mutCodons <- mutCodons[!forbiddencodon]
uniqueMutCodons <- uniqueMutCodons[!forbiddencodon]
splitMutCodons <- as(strsplit(gsub("([^.]{3})", "\\1\\.", 
                                   as.character(mutCodons)), ".", fixed = TRUE), 
                     "CharacterList")
encodedMutCodonsForward <- S4Vectors::unstrsplit(BiocGenerics::relist(
  paste0("f", unlist(uniqueMutCodons), unlist(splitMutCodons)), 
  uniqueMutCodons),
  sep = "_")

mutNames <- encodedMutCodonsForward
mutNames <- gsub("^_", "", gsub("_$", "", mutNames))
mutNames[mutNames <- ""] <- "WT"
sum(table(mutNames) == 2)  ## 11
which(table(mutNames) == 2)  ## f14AAG, f15ATG, f19CAC, f1ACC, f20AAC, f21GTG, f24AGC, 
                             ## f4CGC f7GTG, f9GGC, f9GTC

## Retained reads
length(fq1)  ## 167

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
nucleotideMismatches <- relist(unlisted_ans, skeleton)
offsets <- c(0L, breakpoints[-length(breakpoints)])
mismatchPositions <- nucleotideMismatches - offsets
mismatchQualities <- unlist(as(quality(constForward)[mismatchPositions], "IntegerList"))

unlisted_ans <- which(as.raw(unlist(p)) == as.raw(unlist(constForward)))
breakpoints <- cumsum(pattern_width)
ans_eltlens <- tabulate(findInterval(unlisted_ans - 1L,
                                     breakpoints) + 1L,
                        nbins = length(p))
skeleton <- IRanges::PartitioningByEnd(cumsum(ans_eltlens))
nucleotideMatches <- relist(unlisted_ans, skeleton)
offsets <- c(0L, breakpoints[-length(breakpoints)])
matchPositions <- nucleotideMatches - offsets
matchQualities <- unlist(as(quality(constForward)[matchPositions], "IntegerList"))

table(matchQualities)
table(mismatchQualities)

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
nucleotideMismatches <- relist(unlisted_ans, skeleton)
offsets <- c(0L, breakpoints[-length(breakpoints)])
mismatchPositions <- nucleotideMismatches - offsets
mismatchQualities <- unlist(as(quality(constReverse)[mismatchPositions], "IntegerList"))

unlisted_ans <- which(as.raw(unlist(p)) == as.raw(unlist(constReverse)))
breakpoints <- cumsum(pattern_width)
ans_eltlens <- tabulate(findInterval(unlisted_ans - 1L,
                                     breakpoints) + 1L,
                        nbins = length(p))
skeleton <- IRanges::PartitioningByEnd(cumsum(ans_eltlens))
nucleotideMatches <- relist(unlisted_ans, skeleton)
offsets <- c(0L, breakpoints[-length(breakpoints)])
matchPositions <- nucleotideMatches - offsets
matchQualities <- unlist(as(quality(constReverse)[matchPositions], "IntegerList"))

table(matchQualities)
table(mismatchQualities)




