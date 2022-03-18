## Generate test data for linkMultipleVariants
## Construct: C(6) - V(24) - C(10) - V(30) - C(10) - V(40) - C(6)
## Forward read length 90
## Reverse read length 70

nReads <- 100
status <- rep("OK", nReads)

## First constant region
c1 <- rep("ACGTGA", nReads)
set.seed(1)
i1 <- sample(seq_len(nReads), 4)  ## reads with mutation
c1[i1] <- "ACTTGA"
status[i1] <- gsub("OK,", "", paste0(status[i1], ",C1mut"))

## Second constant region
c2 <- rep("TTGATCCCGG", nReads)
set.seed(2)
i2 <- sample(seq_len(nReads), 5)  ## reads with mutation
c2[i2] <- "TTGATCGCGG"
status[i2] <- gsub("OK,", "", paste0(status[i2], ",C2mut"))

## Third constant region
c3 <- rep("TAAATACCGG", nReads)
set.seed(3)
i3 <- sample(seq_len(nReads), 5)  ## reads with mutation
c3[i3] <- "TAGATACCGG"
status[i3] <- gsub("OK,", "", paste0(status[i3], ",C3mut"))

## Fourth constant region
c4 <- rep("ACCTGA", nReads)
set.seed(4)
i4 <- sample(seq_len(nReads), 4)  ## reads with mutation
c4[i4] <- "ACCCGA"
status[i4] <- gsub("OK,", "", paste0(status[i4], ",C4mut"))

## First variable region (barcodes)
tmpstatus <- rep("", nReads)
set.seed(1)
v1true <- rep(c(paste(sample(c("A", "C", "G", "T"), 24, replace = TRUE), collapse = ""),
                paste(sample(c("A", "C", "G", "T"), 24, replace = TRUE), collapse = ""),
                paste(sample(c("A", "C", "G", "T"), 24, replace = TRUE), collapse = "")),
              c(floor(nReads/3), floor(nReads/3), nReads - 2 * floor(nReads/3)))
v1 <- v1true
i1 <- sample(seq_len(nReads), 10)  ## reads with mutation
for (i in i1) {
    tmp <- strsplit(v1[i], "")[[1]]
    i0 <- sample(seq_along(tmp), 1)
    tmp[i0] <- ifelse(tmp[i0] == "A", "T", "A")
    v1[i] <- paste(tmp, collapse = "")
    tmpstatus[i] <- ",V1mut"
}
i2 <- sample(seq_len(nReads), 4)  ## reads with N
for (i in i2) {
    tmp <- strsplit(v1[i], "")[[1]]
    tmp[sample(seq_along(tmp), 1)] <- "N"
    v1[i] <- paste(tmp, collapse = "")
    tmpstatus[i] <- paste0(tmpstatus[i], ",V1N")
}
ordr <- sample(seq_along(v1true), length(v1true))
v1true <- v1true[ordr]
v1 <- v1[ordr]
tmpstatus <- tmpstatus[ordr]
status <- gsub("OK,", "", paste0(status, tmpstatus))

## Second variable region (multiple WT sequences)
tmpstatus <- rep("", nReads)
set.seed(2)
V2v1 <- paste(sample(c("A", "C", "G", "T"), 30, replace = TRUE), collapse = "")
V2v2 <- paste(sample(c("A", "C", "G", "T"), 30, replace = TRUE), collapse = "")
V2v3 <- paste(sample(c("A", "C", "G", "T"), 30, replace = TRUE), collapse = "")
v2true <- rep(c(V2v1, V2v2, V2v3),
              c(floor(nReads/3), floor(nReads/3), nReads - 2 * floor(nReads/3)))
n2true <- rep(c("V2v1", "V2v2", "V2v3"),
              c(floor(nReads/3), floor(nReads/3), nReads - 2 * floor(nReads/3)))
v2 <- v2true
i1 <- sample(seq_len(nReads), 20)   ## reads with mutation
for (i in i1) {
    tmp <- strsplit(v2[i], "")[[1]]
    i0 <- sample(seq_along(tmp), 1)
    repl <- ifelse(tmp[i0] == "A", "T", "A")
    tmp[i0] <- repl
    v2[i] <- paste(tmp, collapse = "")
    tmpstatus[i] <- paste0(",", n2true[i], ".", i0, ".", repl)
}
i3 <- sample(seq_len(nReads), 3)   ## reads with deletion
## don't modify the reads here, since that won't generate an invalid overlap 
## (the modification will be in both the fwd and rev reads). Just record
## which reads to modify
for (i in i3) {
    tmpstatus[i] <- paste0(tmpstatus[i], ",V2del")
}
ordr <- sample(seq_along(v2true), length(v2true))
v2true <- v2true[ordr]
n2true <- n2true[ordr]
v2 <- v2[ordr]
tmpstatus <- tmpstatus[ordr]
status <- gsub("OK,", "", paste0(status, tmpstatus))

## Third variable region (one WT sequence)
tmpstatus <- rep("", nReads)
set.seed(3)
V3 <- paste(sample(c("A", "C", "G", "T"), 40, replace = TRUE), collapse = "")
v3true <- rep(V3, nReads)
v3 <- v3true
i1 <- sample(seq_len(nReads), 30)   ## reads with mutation
for (i in i1) {
    tmp <- strsplit(v3[i], "")[[1]]
    i0 <- sample(seq_along(tmp), 1)
    repl <- ifelse(tmp[i0] == "A", "T", "A")
    tmp[i0] <- repl
    v3[i] <- paste(tmp, collapse = "")
    tmpstatus[i] <- paste0(",V3.", i0, ".", repl)
}
ordr <- sample(seq_along(v3true), length(v3true))
v3true <- v3true[ordr]
v3 <- v3[ordr]
tmpstatus <- tmpstatus[ordr]
status <- gsub("OK,", "", paste0(status, tmpstatus))

## Final reads
reads <- paste0(c1, v1, c2, v2, c3, v3, c4)
qualities <- rep(paste(rep("I", 126), collapse = ""), nReads)

readsR1 <- substr(reads, start = 1, stop = 90)
readsR2 <- substr(reads, start = 57, stop = 126)
qualitiesR1 <- substr(qualities, start = 1, stop = 90)
qualitiesR2 <- substr(qualities, start = 57, stop = 126)

## Add the deletion to the forward reads
idx <- grep("V2del", status)
for (i in idx) {
    tmp <- strsplit(readsR1[i], "")[[1]]
    tmp[41] <- ""
    tmp <- c(tmp, substr(reads[i], start = 91, stop = 91))
    tmp <- paste(tmp, collapse = "")
    readsR1[i] <- tmp
}


fastqR1 <- Biostrings::QualityScaledDNAStringSet(
    x = Biostrings::DNAStringSet(
        x = structure(readsR1, names = paste0("R", seq_len(nReads)))
    ), 
    quality = Biostrings::PhredQuality(
        Biostrings::BStringSet(
            x = structure(qualitiesR1, names = paste0("R", seq_len(nReads)))
        )
    )
)
fastqR2 <- Biostrings::QualityScaledDNAStringSet(
    x = Biostrings::reverseComplement(
        Biostrings::DNAStringSet(
            x = structure(readsR2, names = paste0("R", seq_len(nReads)))
        )
    ), 
    quality = Biostrings::PhredQuality(
        Biostrings::BStringSet(
            x = structure(qualitiesR2, names = paste0("R", seq_len(nReads)))
        )
    )
)

Biostrings::writeQualityScaledXStringSet(fastqR1, filepath = "multipleVariableRegions_R1.fastq.gz", compress = TRUE)
Biostrings::writeQualityScaledXStringSet(fastqR2, filepath = "multipleVariableRegions_R2.fastq.gz", compress = TRUE)

## True sequences
barcode <- v1true
V2var <- vapply(strsplit(status, ","), function(x) {
    i <- grep("V2.*\\.", x)
    if (length(i) > 0) {
        x[i]
    } else {
        NA_character_
    }}, NA_character_)
V2var[is.na(V2var)] <- paste0(n2true[is.na(V2var)], ".0.WT")
V3var <- vapply(strsplit(status, ","), function(x) {
    i <- grep("V3.*\\.", x)
    if (length(i) > 0) {
        x[i]
    } else {
        NA_character_
    }}, NA_character_)
V3var[is.na(V3var)] <- "V3.0.WT"
WTV2 <- c(V2v1 = V2v1, V2v2 = V2v2, V2v3 = V2v3)
WTV3 <- c(V3 = V3)
constFwd <- "ACGTGATTGATCCCGGTAAATACCGG"
constRev <- "TAAATACCGGACCTGA"

truth <- data.frame(read = paste0("R", seq_len(nReads)),
                    trueBarcode = barcode,
                    trueV2 = V2var,
                    trueV3 = V3var,
                    obsBarcode = v1,
                    status = status)

saveRDS(list(truth = truth, WTV2 = WTV2, WTV3 = WTV3,
             constFwd = constFwd, constRev = constRev),
        "multipleVariableRegions_truth.rds")

