suppressPackageStartupMessages({
  library(ShortRead)
})

## Trans, Input
r1 <- ShortRead::readFastq("../../../data/GSE102901/FASTQ/SRR5952429_1.fastq.gz", withIds = TRUE)
r2 <- ShortRead::readFastq("../../../data/GSE102901/FASTQ/SRR5952429_2.fastq.gz", withIds = TRUE)
ShortRead::writeFastq(r1[1:1000], file = "transInput_1.fastq.gz", mode = "w", full = FALSE, compress = TRUE)
ShortRead::writeFastq(r2[1:1000], file = "transInput_2.fastq.gz", mode = "w", full = FALSE, compress = TRUE)

## Trans, Output
r1 <- ShortRead::readFastq("../../../data/GSE102901/FASTQ/SRR5952432_1.fastq.gz", withIds = TRUE)
r2 <- ShortRead::readFastq("../../../data/GSE102901/FASTQ/SRR5952432_2.fastq.gz", withIds = TRUE)
ShortRead::writeFastq(r1[1:1000], file = "transOutput_1.fastq.gz", mode = "w", full = FALSE, compress = TRUE)
ShortRead::writeFastq(r2[1:1000], file = "transOutput_2.fastq.gz", mode = "w", full = FALSE, compress = TRUE)

## Cis, Input
r1 <- ShortRead::readFastq("../../../data/GSE102901/FASTQ/SRR5952435_1.fastq.gz", withIds = TRUE)
r2 <- ShortRead::readFastq("../../../data/GSE102901/FASTQ/SRR5952435_2.fastq.gz", withIds = TRUE)
ShortRead::writeFastq(r1[1:1000], file = "cisInput_1.fastq.gz", mode = "w", full = FALSE, compress = TRUE)
ShortRead::writeFastq(r2[1:1000], file = "cisInput_2.fastq.gz", mode = "w", full = FALSE, compress = TRUE)

## Cis, Output
r1 <- ShortRead::readFastq("../../../data/GSE102901/FASTQ/SRR5952438_1.fastq.gz", withIds = TRUE)
r2 <- ShortRead::readFastq("../../../data/GSE102901/FASTQ/SRR5952438_2.fastq.gz", withIds = TRUE)
ShortRead::writeFastq(r1[1:1000], file = "cisOutput_1.fastq.gz", mode = "w", full = FALSE, compress = TRUE)
ShortRead::writeFastq(r2[1:1000], file = "cisOutput_2.fastq.gz", mode = "w", full = FALSE, compress = TRUE)

## Count matrix
meta <- read.delim("../../../data/GSE102901/meta/GSE102901_meta.tsv", header = TRUE, as.is = TRUE)
sampleIds <- meta$Run
names(sampleIds) <- sampleIds
allSamples <- lapply(sampleIds, function(id) {
  print(id)
  expType <- strsplit(meta$title[match(id, meta$Run)], " ")[[1]][1]
  if (expType == "trans") {
    tmp <- mutscan::readFastqs(experimentType = expType, 
                               fastqForward = paste0("../../../data/GSE102901/FASTQ/", id, "_1.fastq.gz"),
                               fastqReverse = paste0("../../../data/GSE102901/FASTQ/", id, "_2.fastq.gz"),
                               skipForward = 1, skipReverse = 1, umiLengthForward = 10, 
                               umiLengthReverse = 8, constantLengthForward = 18,
                               constantLengthReverse = 20, adapterForward = "GGAAGAGCACACGTC", 
                               adapterReverse = "GGAAGAGCGTCGTGT", verbose = TRUE)
    tmp <- filterReads(tmp, avePhredMin = 20, variableNMax = 0, umiNMax = 0,
                       wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA", 
                       wildTypeReverse = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT", 
                       nbrMutatedCodonsMax = 1, forbiddenMutatedCodon = "NNW",
                       verbose = TRUE)
  } else {
    tmp <- mutscan::readFastqs(experimentType = expType, 
                               fastqForward = paste0("../../../data/GSE102901/FASTQ/", id, "_1.fastq.gz"),
                               fastqReverse = paste0("../../../data/GSE102901/FASTQ/", id, "_2.fastq.gz"),
                               skipForward = 1, skipReverse = 1, umiLengthForward = 10, 
                               umiLengthReverse = 7, constantLengthForward = 18,
                               constantLengthReverse = 17, adapterForward = "GGAAGAGCACACGTC", 
                               adapterReverse = "GGAAGAGCGTCGTGT", verbose = TRUE)
    tmp <- filterReads(tmp, avePhredMin = 20, variableNMax = 0, umiNMax = 0,
                       wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA", 
                       wildTypeReverse = NULL, 
                       nbrMutatedCodonsMax = 1, forbiddenMutatedCodon = "NNW", 
                       verbose = TRUE)
  }
  tmp
})
  