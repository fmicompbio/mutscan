suppressPackageStartupMessages({
  library(ShortRead)
  library(BiocParallel)
})

## Trans, Input
r1 <- ShortRead::readFastq("/tungstenfs/groups/gbioinfo/sonechar/FMIgroups/GroupDiss/Diss_deep_mutational_scanning/data/GSE102901/FASTQ/SRR5952429_1.fastq.gz", withIds = TRUE)
r2 <- ShortRead::readFastq("/tungstenfs/groups/gbioinfo/sonechar/FMIgroups/GroupDiss/Diss_deep_mutational_scanning/data/GSE102901/FASTQ/SRR5952429_2.fastq.gz", withIds = TRUE)
ShortRead::writeFastq(r1[1:1000], file = "transInput_1.fastq.gz", mode = "w", full = FALSE, compress = TRUE)
ShortRead::writeFastq(r2[1:1000], file = "transInput_2.fastq.gz", mode = "w", full = FALSE, compress = TRUE)

## Trans, Output
r1 <- ShortRead::readFastq("/tungstenfs/groups/gbioinfo/sonechar/FMIgroups/GroupDiss/Diss_deep_mutational_scanning/data/GSE102901/FASTQ/SRR5952432_1.fastq.gz", withIds = TRUE)
r2 <- ShortRead::readFastq("/tungstenfs/groups/gbioinfo/sonechar/FMIgroups/GroupDiss/Diss_deep_mutational_scanning/data/GSE102901/FASTQ/SRR5952432_2.fastq.gz", withIds = TRUE)
ShortRead::writeFastq(r1[1:1000], file = "transOutput_1.fastq.gz", mode = "w", full = FALSE, compress = TRUE)
ShortRead::writeFastq(r2[1:1000], file = "transOutput_2.fastq.gz", mode = "w", full = FALSE, compress = TRUE)

## Cis, Input
r1 <- ShortRead::readFastq("/tungstenfs/groups/gbioinfo/sonechar/FMIgroups/GroupDiss/Diss_deep_mutational_scanning/data/GSE102901/FASTQ/SRR5952435_1.fastq.gz", withIds = TRUE)
r2 <- ShortRead::readFastq("/tungstenfs/groups/gbioinfo/sonechar/FMIgroups/GroupDiss/Diss_deep_mutational_scanning/data/GSE102901/FASTQ/SRR5952435_2.fastq.gz", withIds = TRUE)
ShortRead::writeFastq(r1[1:1000], file = "cisInput_1.fastq.gz", mode = "w", full = FALSE, compress = TRUE)
ShortRead::writeFastq(r2[1:1000], file = "cisInput_2.fastq.gz", mode = "w", full = FALSE, compress = TRUE)

## Cis, Output
r1 <- ShortRead::readFastq("/tungstenfs/groups/gbioinfo/sonechar/FMIgroups/GroupDiss/Diss_deep_mutational_scanning/data/GSE102901/FASTQ/SRR5952438_1.fastq.gz", withIds = TRUE)
r2 <- ShortRead::readFastq("/tungstenfs/groups/gbioinfo/sonechar/FMIgroups/GroupDiss/Diss_deep_mutational_scanning/data/GSE102901/FASTQ/SRR5952438_2.fastq.gz", withIds = TRUE)
ShortRead::writeFastq(r1[1:1000], file = "cisOutput_1.fastq.gz", mode = "w", full = FALSE, compress = TRUE)
ShortRead::writeFastq(r2[1:1000], file = "cisOutput_2.fastq.gz", mode = "w", full = FALSE, compress = TRUE)

## Count matrix, cis
meta <- read.delim("/tungstenfs/groups/gbioinfo/sonechar/FMIgroups/GroupDiss/Diss_deep_mutational_scanning/data/GSE102901/meta/GSE102901_cis_coldata.tsv", header = TRUE, as.is = TRUE)
sampleIds <- meta$Name
names(sampleIds) <- sampleIds
allSamples <- lapply(sampleIds, function(id) {
  print(id)
  mutscan::digestFastqs(fastqForward = paste0("/tungstenfs/groups/gbioinfo/sonechar/FMIgroups/GroupDiss/Diss_deep_mutational_scanning/data/GSE102901/FASTQ/", id, "_1.fastq.gz"),
                        fastqReverse = paste0("/tungstenfs/groups/gbioinfo/sonechar/FMIgroups/GroupDiss/Diss_deep_mutational_scanning/data/GSE102901/FASTQ/", id, "_2.fastq.gz"),
                        mergeForwardReverse = TRUE, minOverlap = 0, maxOverlap = 0,
                        maxFracMismatchOverlap = 1, greedyOverlap = TRUE, 
                        revComplForward = FALSE, revComplReverse = TRUE,
                        elementsForward = "SUCV", elementsReverse = "SUCVS",
                        elementLengthsForward = c(1, 10, 18, 96),
                        elementLengthsReverse = c(1, 7, 17, 96, -1),
                        constantForward = "AACCGGAGGAGGGAGCTG",
                        constantReverse = "GAGTTCATCCTGGCAGC", adapterForward = "GGAAGAGCACACGTC", 
                        adapterReverse = "GGAAGAGCGTCGTGT", primerForward = "",
                        primerReverse = "", 
                        avePhredMinForward = 20, avePhredMinReverse = 20,
                        variableNMaxForward = 0, variableNMaxReverse = 0, umiNMax = 0,
                        wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA", 
                        wildTypeReverse = "", 
                        nbrMutatedCodonsMaxForward = 1, nbrMutatedCodonsMaxReverse = 1,
                        forbiddenMutatedCodonsForward = "NNW", forbiddenMutatedCodonsReverse = "NNW", 
                        mutNameDelimiter = ".", maxNReads = -1, verbose = TRUE, nThreads = 6)
})
se <- mutscan::summarizeExperiment(x = allSamples, coldata = meta, countType = "umis")
seaa <- mutscan::collapseMutantsByAA(se)
saveRDS(se, file = "GSE102901_cis_se.rds")

## Leu+Jun
r1 <- ShortRead::readFastq("/tungstenfs/groups/gbioinfo/sonechar/FMIgroups/GroupDiss/Diss_deep_mutational_scanning/data/LEU_JUN_6rep_3tp/FASTQ/t0/GD09_26674_AACGGT_read1.fastq.gz", withIds = TRUE)
r2 <- ShortRead::readFastq("/tungstenfs/groups/gbioinfo/sonechar/FMIgroups/GroupDiss/Diss_deep_mutational_scanning/data/LEU_JUN_6rep_3tp/FASTQ/t0/GD09_26674_AACGGT_read2.fastq.gz", withIds = TRUE)
ShortRead::writeFastq(r1[22501:23500], file = "leujunt0_1.fastq.gz", mode = "w", full = FALSE, compress = TRUE)
ShortRead::writeFastq(r2[22501:23500], file = "leujunt0_2.fastq.gz", mode = "w", full = FALSE, compress = TRUE)
