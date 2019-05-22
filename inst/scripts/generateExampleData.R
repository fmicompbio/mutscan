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
