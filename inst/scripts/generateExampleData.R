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

## Count matrix, trans
# meta <- read.delim("../../../data/GSE102901/meta/GSE102901_trans_coldata.tsv", header = TRUE, as.is = TRUE)
# sampleIds <- meta$Name
# names(sampleIds) <- sampleIds
# allSamples <- mclapply(sampleIds, function(id) {
#   print(id)
#   expType <- "trans"
#   mutscan::digestFastqs(experimentType = expType, 
#                         fastqForward = paste0("../../../data/GSE102901/FASTQ/", id, "_1.fastq.gz"),
#                         fastqReverse = paste0("../../../data/GSE102901/FASTQ/", id, "_2.fastq.gz"),
#                         skipForward = 1, skipReverse = 1, umiLengthForward = 10, 
#                         umiLengthReverse = 8, constantLengthForward = 18,
#                         constantLengthReverse = 20, variableLengthForward = 96,
#                         variableLengthReverse = 96, constantForward = "AACCGGAGGAGGGAGCTG",
#                         constantReverse = "GAAAAAGGAAGCTGGAGAGA", adapterForward = "GGAAGAGCACACGTC", 
#                         adapterReverse = "GGAAGAGCGTCGTGT", 
#                         avePhredMin = 20, variableNMax = 0, umiNMax = 0,
#                         wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA", 
#                         wildTypeReverse = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT", 
#                         nbrMutatedCodonsMax = 1, forbiddenMutatedCodon = "NNW",verbose = TRUE)
# }, mc.cores = 6)
# se <- mutscan::summarizeExperiment(x = allSamples, coldata = meta)
# seaa <- mutscan::collapseMutantsByAA(se)
# saveRDS(list(se = se, se_collapsed = seaa), file = "../../../data/GSE102901/processed_data/GSE102901_trans_se.rds")

## Count matrix, cis
meta <- read.delim("../../../data/GSE102901/meta/GSE102901_cis_coldata.tsv", header = TRUE, as.is = TRUE)
sampleIds <- meta$Name
names(sampleIds) <- sampleIds
allSamples <- mclapply(sampleIds, function(id) {
  print(id)
  mutscan::digestFastqs(fastqForward = paste0("../../../data/GSE102901/FASTQ/", id, "_1.fastq.gz"),
                        fastqReverse = paste0("../../../data/GSE102901/FASTQ/", id, "_2.fastq.gz"),
                        mergeForwardReverse = TRUE, revComplForward = FALSE, revComplReverse = TRUE,
                        skipForward = 1, skipReverse = 1, umiLengthForward = 10, 
                        umiLengthReverse = 7, constantLengthForward = 18,
                        constantLengthReverse = 17, variableLengthForward = 96,
                        variableLengthReverse = 96, constantForward = "AACCGGAGGAGGGAGCTG",
                        constantReverse = "GAGTTCATCCTGGCAGC", adapterForward = "GGAAGAGCACACGTC", 
                        adapterReverse = "GGAAGAGCGTCGTGT", 
                        avePhredMinForward = 20, avePhredMinReverse = 20,
                        variableNMaxForward = 0, variableNMaxReverse = 0, umiNMax = 0,
                        wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA", 
                        wildTypeReverse = "", 
                        nbrMutatedCodonsMaxForward = 1, nbrMutatedCodonsMaxReverse = 1,
                        forbiddenMutatedCodonsForward = "NNW", forbiddenMutatedCodonsReverse = "NNW", 
                        verbose = TRUE)
}, mc.cores = 6)
se <- mutscan::summarizeExperiment(x = allSamples, coldata = meta, countType = "umis")
seaa <- mutscan::collapseMutantsByAA(se)
saveRDS(se, file = "GSE102901_cis_se.rds")
