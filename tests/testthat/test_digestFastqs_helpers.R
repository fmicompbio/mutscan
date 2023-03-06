## ----------------------------------------------------------------------------
## compareCodonPositions
## ----------------------------------------------------------------------------
test_that("compareCodonPositions works", {
  expect_true(compareCodonPositions("f.1.AAA_", "f.2.ACT_", "."))
  expect_true(compareCodonPositions("f.2.ACT_", "f.10.TCA_", "."))
})

## ----------------------------------------------------------------------------
## findClosestRefSeq
## ----------------------------------------------------------------------------
test_that("findClosestRefSeq works", {
  expect_equal(findClosestRefSeq(varSeq = "ACGT", wtSeq = c("ATTT", "ACTT"),
                                 upperBoundMismatch = 4L, sim = 0L), 1L)
  expect_equal(findClosestRefSeq(varSeq = "ACGT", wtSeq = c("ACGT", "ACTT"),
                                 upperBoundMismatch = 4L, sim = 0L), 0L)
  expect_equal(findClosestRefSeq(varSeq = "ACGT", wtSeq = c("ACG", "ACGT"),
                                 upperBoundMismatch = 4L, sim = 0L), 1L)
  expect_equal(findClosestRefSeq(varSeq = "ACGT", wtSeq = c("AACGT", "ACCTA"),
                                 upperBoundMismatch = 4L, sim = 0L), 1L)
  expect_equal(findClosestRefSeq(varSeq = "ACGT", wtSeq = c("ATTT", "ACTT"),
                                 upperBoundMismatch = 0L, sim = 0L), -1L)
  expect_equal(findClosestRefSeq(varSeq = "ACGT", wtSeq = c("ATGT", "ACTT"),
                                 upperBoundMismatch = 1L, sim = 0L), -2L)

  expect_equal(findClosestRefSeqEarlyStop(varSeq = "ACGT", wtSeq = c("ATTT", "ACTT"),
                                          upperBoundMismatch = 4L, sim = 0L), 1L)
  expect_equal(findClosestRefSeqEarlyStop(varSeq = "ACGT", wtSeq = c("ACGT", "ACTT"),
                                          upperBoundMismatch = 4L, sim = 0L), 0L)
  expect_equal(findClosestRefSeqEarlyStop(varSeq = "ACGT", wtSeq = c("ACG", "ACGT"),
                                          upperBoundMismatch = 4L, sim = 0L), 1L)
  expect_equal(findClosestRefSeqEarlyStop(varSeq = "ACGT", wtSeq = c("AACGT", "ACCTA"),
                                          upperBoundMismatch = 4L, sim = 0L), 1L)
  expect_equal(findClosestRefSeqEarlyStop(varSeq = "ACGT", wtSeq = c("ATTT", "ACTT"),
                                          upperBoundMismatch = 0L, sim = 0L), -1L)
  expect_equal(findClosestRefSeqEarlyStop(varSeq = "ACGT", wtSeq = c("ATGT", "ACTT"),
                                          upperBoundMismatch = 1L, sim = 0L), -2L)
})

## ----------------------------------------------------------------------------
## decomposeRead
## ----------------------------------------------------------------------------
test_that("decomposeRead works", {
  ## All element lengths specified, no primer
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "USCVCUSV",
    elementLengths = c(2, 1, 3, 3, 1, 2, 1, 1), primerSeqs = "",
    varSeq = "", varQual = "", varLengths = integer(0L), umiSeq = "",
    constSeq = "", constQual = "", nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "ACCT")
  expect_equal(res$varSeq, "AAGG")
  expect_equal(res$varQual, "gABG")
  expect_equal(res$varLengths, c(3L, 1L))
  expect_equal(res$constSeq, "TTAC")
  expect_equal(res$constQual, "defC")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)

  ## All element lengths specified, no primer, 0 lengths included
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTG", squal = "abcdefgABCDEF", elements = "SUSCVUCUS",
    elementLengths = c(0, 2, 1, 3, 3, 0, 1, 2, 1), primerSeqs = "",
    varSeq = "", varQual = "", varLengths = integer(0L), umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "ACCT")
  expect_equal(res$varSeq, "AAG")
  expect_equal(res$varQual, "gAB")
  expect_equal(res$varLengths, 3L)
  expect_equal(res$constSeq, "TTAC")
  expect_equal(res$constQual, "defC")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)

  ## Derive variable length from other lengths
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTG", squal = "abcdefgABCDEF", elements = "USCVCUS",
    elementLengths = c(2, 1, 3, -1, 1, 2, 1), primerSeqs = "",
    varSeq = "", varQual = "", varLengths = integer(0L), umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "ACCT")
  expect_equal(res$varSeq, "AAG")
  expect_equal(res$varQual, "gAB")
  expect_equal(res$varLengths, 3L)
  expect_equal(res$constSeq, "TTAC")
  expect_equal(res$constQual, "defC")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)

  ## Primer
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "USCVPS",
    elementLengths = c(2, 1, 3, 3, 3, -1), primerSeqs = c("CCT"),
    varSeq = "", varQual = "", varLengths = integer(0L), umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "AC")
  expect_equal(res$varSeq, "AAG")
  expect_equal(res$varQual, "gAB")
  expect_equal(res$varLengths, 3L)
  expect_equal(res$constSeq, "TTA")
  expect_equal(res$constQual, "def")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)

  ## Only variable sequence
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "V",
    elementLengths = -1, primerSeqs = c("CCT"),
    varSeq = "", varQual = "", varLengths = integer(0L), umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "ACGTTAAAGCCTGG")
  expect_equal(res$varQual, "abcdefgABCDEFG")
  expect_equal(res$varLengths, 14L)
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)

  ## Primer + variable sequence
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "PV",
    elementLengths = c(-1, -1), primerSeqs = c("ACG"),
    varSeq = "", varQual = "", varLengths = integer(0L), umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "TTAAAGCCTGG")
  expect_equal(res$varQual, "defgABCDEFG")
  expect_equal(res$varLengths, 11L)
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)

  ## Variable sequence + primer
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "VP",
    elementLengths = c(-1, -1), primerSeqs = c("ACG"),
    varSeq = "A", varQual = "E", varLengths = 1L, umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "A")
  expect_equal(res$varQual, "E")
  expect_equal(res$varLengths, c(1L, 0L))
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)

  ## Primer, not the right length of sequence
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "VPU",
    elementLengths = c(10, 3, 2), primerSeqs = c("CCT"),
    varSeq = "", varQual = "", varLengths = integer(0L), umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "")
  expect_equal(res$varQual, "")
  expect_equal(res$varLengths, integer(0))
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 1)

  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "VPU",
    elementLengths = c(8, 3, 2), primerSeqs = c("CCT"),
    varSeq = "", varQual = "", varLengths = integer(0), umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "")
  expect_equal(res$varQual, "")
  expect_equal(res$varLengths, integer(0))
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 1)

  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "CVPU",
    elementLengths = c(10, -1, 3, 2), primerSeqs = c("CCT"),
    varSeq = "", varQual = "", varLengths = integer(0), umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "")
  expect_equal(res$varQual, "")
  expect_equal(res$varLengths, integer(0))
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 1)

  ## Primer, but don't explicitly specify that the part after the primer should be skipped
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "USCVP",
    elementLengths = c(2, 1, 3, 3, 3), primerSeqs = c("CCT"),
    varSeq = "", varQual = "", varLengths = integer(0), umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "AC")
  expect_equal(res$varSeq, "AAG")
  expect_equal(res$varQual, "gAB")
  expect_equal(res$varLengths, 3L)
  expect_equal(res$constSeq, "TTA")
  expect_equal(res$constQual, "def")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)

  ## Multiple primers, only one present
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "USCVP",
    elementLengths = c(2, 1, 3, 3, 3), primerSeqs = c("TTT", "CCT"),
    varSeq = "", varQual = "", varLengths = integer(0), umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "AC")
  expect_equal(res$varSeq, "AAG")
  expect_equal(res$varQual, "gAB")
  expect_equal(res$varLengths, 3L)
  expect_equal(res$constSeq, "TTA")
  expect_equal(res$constQual, "def")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)

  ## Multiple primers, none present
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "USCVP",
    elementLengths = c(2, 1, 3, 3, 3), primerSeqs = c("TTT", "GGG"),
    varSeq = "", varQual = "", varLengths = integer(0), umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "")
  expect_equal(res$varQual, "")
  expect_equal(res$varLengths, integer(0))
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 1)
  expect_equal(res$nReadWrongLength, 0)

  ## Multiple primers, all present (use the first one)
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "USCVP",
    elementLengths = c(2, 1, 3, 3, 3), primerSeqs = c("CCT", "TTA"),
    varSeq = "", varQual = "", varLengths = integer(0), umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "AC")
  expect_equal(res$varSeq, "AAG")
  expect_equal(res$varQual, "gAB")
  expect_equal(res$varLengths, 3L)
  expect_equal(res$constSeq, "TTA")
  expect_equal(res$constQual, "def")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)

  ## Primer, infer length of sequences
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "UPV",
    elementLengths = c(-1, 3, -1), primerSeqs = c("CCT"),
    varSeq = "", varQual = "", varLengths = integer(0), umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "ACGTTAAAG")
  expect_equal(res$varSeq, "GG")
  expect_equal(res$varQual, "FG")
  expect_equal(res$varLengths, 2L)
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)

  ## Primer, don't specify primer length
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "UPV",
    elementLengths = c(-1, -1, -1), primerSeqs = c("CCT"),
    varSeq = "", varQual = "", varLengths = integer(0), umiSeq = "",
    constSeq = "", constQual = "",
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "ACGTTAAAG")
  expect_equal(res$varSeq, "GG")
  expect_equal(res$varQual, "FG")
  expect_equal(res$varLengths, 2L)
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)

  ## Multiple variable segments with different lengths
  res <- mutscan:::test_decomposeRead(
      sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "VPV",
      elementLengths = c(-1, -1, -1), primerSeqs = c("CCT"),
      varSeq = "", varQual = "", varLengths = integer(0), umiSeq = "",
      constSeq = "", constQual = "",
      nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "ACGTTAAAGGG")
  expect_equal(res$varQual, "abcdefgABFG")
  expect_equal(res$varLengths, c(9L, 2L))
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)
})

## ----------------------------------------------------------------------------
## mergeReadPairsPartial
## ----------------------------------------------------------------------------
test_that("mergeReadPairsPartial works", {
  ## don't count N as mismatch, pick base with higher quality
  sF1 <- "AAAANA"; qF1 <- rep(40L, nchar(sF1)); lF1 <- 6L
  sR1 <- "AACCCC"; qR1 <- rep(42L, nchar(sR1)); lR1 <- 6L
  res1a <- mutscan:::test_mergeReadPairPartial(sF1, qF1, sR1, qR1, lF1, lR1, 1, 6, 0, 0, 1/4, TRUE)
  res1b <- mutscan:::test_mergeReadPairPartial(sF1, qF1, sR1, qR1, lF1, lR1, 1, 6, 0, 0, 1/4, FALSE)
  res1c <- mutscan:::test_mergeReadPairPartial(sF1, qR1, sR1, qF1, lF1, lR1, 2, 2, 0, 0, 1/4, TRUE)
  res1d <- mutscan:::test_mergeReadPairPartial(sF1, qF1, sR1, qR1, lF1, lR1, 2, 2, 0, 0, 0, TRUE)
  expect_type(res1a, "list")
  expect_type(res1b, "list")
  expect_type(res1c, "list")
  expect_type(res1d, "list")
  expect_identical(res1a$mergedSeq, "AAAACCCC")
  expect_identical(res1a$mergedQual, rep(c(40L, 42L), c(2, 6)))
  expect_identical(res1a$mergedLengths, 8L)
  expect_identical(res1a, res1b)
  expect_identical(res1c$mergedSeq, "AAAANACCCC")
  expect_identical(res1d$mergedSeq, "AAAAAACCCC")
  expect_identical(res1c$mergedLengths, 10L)
  expect_identical(res1d$mergedLengths, 10L)

  ## no valid overlap, multiple possible overlaps with single mismatch
  sF2 <- "TTACACG"; qF2 <- rep(10L, nchar(sF2)); lF2 <- 7L
  sR2 <- "ACACACA"; qR2 <- rep(40L, nchar(sR2)); lR2 <- 7L
  res2a <- mutscan:::test_mergeReadPairPartial(sF2, qF2, sR2, qR2, lF2, lR2, 0, 0, 0, 0, maxFracMismatchOverlap = 3/7, TRUE)
  res2b <- mutscan:::test_mergeReadPairPartial(sF2, qF2, sR2, qR2, lF2, lR2, 1, 7, 0, 0, 1/5, TRUE)
  expect_type(res2a, "list")
  expect_type(res2b, "list")
  expect_identical(res2a$mergedSeq, sR2)
  expect_identical(res2a$mergedQual, qR2)
  expect_identical(res2a$mergedLengths, 7L)
  expect_identical(res2b$mergedSeq, "TTACACACA")
  expect_identical(res2b$mergedQual, rep(c(10L,40L), c(2,7)))
  expect_identical(res2b$mergedLengths, 9L)

  ## specify min/max merged length instead
  res2c <- mutscan:::test_mergeReadPairPartial(sF2, qF2, sR2, qR2, lF2, lR2, 1, 7, 9, 9, 1, FALSE)
  res2d <- mutscan:::test_mergeReadPairPartial(sF2, qF2, sR2, qR2, lF2, lR2, 1, 7, 11, 11, 1, FALSE)
  res2e <- mutscan:::test_mergeReadPairPartial(sF2, qF2, sR2, qR2, lF2, lR2, 1, 7, 12, 12, 1, FALSE)
  expect_type(res2c, "list")
  expect_type(res2d, "list")
  expect_type(res2e, "list")
  expect_identical(res2c$mergedSeq, "TTACACACA")
  expect_identical(res2c$mergedQual, rep(c(10L, 40L), c(2, 7)))
  expect_identical(res2c$mergedLengths, 9L)
  expect_identical(res2d$mergedSeq, "TTACACACACA")
  expect_identical(res2d$mergedQual, rep(c(10L, 40L), c(4, 7)))
  expect_identical(res2d$mergedLengths, 11L)
  expect_identical(res2e$mergedSeq, "TTACAACACACA")
  expect_identical(res2e$mergedQual, rep(c(10L, 40L), c(5, 7)))
  expect_identical(res2e$mergedLengths, 12L)
  
  ## invalid overlaps specified/overlaps modified internally
  ## minOverlap > lenF -> no valid overlap
  res0a <- mutscan:::test_mergeReadPairPartial(sF2, qF2, sR2, qR2, lF2, lR2, 8, 7, 9, 9, 1, FALSE)
  expect_type(res0a, "list")
  expect_identical(res0a$mergedSeq, sF2)
  expect_identical(res0a$mergedQual, qF2)
  expect_identical(res0a$mergedLengths, 7L)
  expect_true(res0a$return)
  ## minOverlap = 0, lenF > lenR -> minOverlap = lenR
  sF0 <- "TTACACG"; qF0 <- rep(10L, nchar(sF0)); lF0 <- nchar(sF0)
  sR0 <- "ACACG"; qR0 <- rep(40L, nchar(sR0)); lR0 <- nchar(sR0)
  res0b <- mutscan:::test_mergeReadPairPartial(sF0, qF0, sR0, qR0, lF0, lR0, 0, 7, 5, 9, 1, FALSE)
  expect_type(res0b, "list")
  expect_identical(res0b$mergedSeq, sF0)
  expect_identical(res0b$mergedQual, rep(c(10L, 40L), c(2, 5)))
  expect_identical(res0b$mergedLengths, 7L)
  expect_false(res0b$return)
  ## maxOverlap = 0, lenF > lenR -> maxOverlap = lenR
  res0c <- mutscan:::test_mergeReadPairPartial(sF0, qF0, sR0, qR0, lF0, lR0, 1, 0, 5, 9, 1, FALSE)
  expect_type(res0c, "list")
  expect_identical(res0c$mergedSeq, sF0)
  expect_identical(res0c$mergedQual, rep(c(10L, 40L), c(2, 5)))
  expect_identical(res0c$mergedLengths, 7L)
  expect_false(res0c$return)
  ## minOverlap > maxOverlap -> no valid overlap
  res0d <- mutscan:::test_mergeReadPairPartial(sF2, qF2, sR2, qR2, lF2, lR2, 4, 3, 9, 9, 1, FALSE)
  expect_type(res0d, "list")
  expect_identical(res0d$mergedSeq, sF2)
  expect_identical(res0d$mergedQual, qF2)
  expect_identical(res0d$mergedLengths, 7L)
  expect_true(res0d$return)

  ## padded reads
  for (i in 1:10) {
    sF <- paste(rep(c("C","A"), c(i, 6)), collapse = "")
    qF <- rep(30L, nchar(sF))
    sR <- paste(rep(c("A","C"), c(6, i)), collapse = "")
    qR <- rep(30L, nchar(sR))
    lF <- i + 6
    lR <- 6 + i
    res <- mutscan:::test_mergeReadPairPartial(sF, qF, sR, qR, lF, lR, 6, 6, 0, 0, 0)
    expect_identical(res$mergedSeq, paste(rep(c("C","A","C"), c(i, 6, i)), collapse = ""))
    expect_equal(res$mergedLengths, i + 6 + i)
  }

  ## reads of unequal length
  sR <- paste(rep("A", 6), collapse = "")
  qR <- rep(30L, nchar(sR))
  lR <- 6L
  for (i in 1:10) {
    sF <- paste(rep(c("C","A"), c(i, 6)), collapse = "")
    qF <- rep(30L, nchar(sF))
    lF <- i + 6
    res <- mutscan:::test_mergeReadPairPartial(sF, qF, sR, qR, lF, lR, 6, 6, 0, 0, 0)
    expect_identical(res$mergedSeq, sF)
    expect_equal(res$mergedLengths, lF)
    res <- mutscan:::test_mergeReadPairPartial(sR, qR, sF, qF, lR, lF, 6, 6, 0, 0, 0)
    expect_identical(res$mergedSeq, sR)
    expect_equal(res$mergedLengths, lR)
  }

  ## test variable segment lengths merging
  ## ... that must work
  sF3 <- "AAAAAAAAACGTCCCAAC";   qF3 <- rep(40L, nchar(sF3)); lF3 <- c(13L,  5L)
  sR3 <- "CGTCCCAACCCGGGGGGGGG"; qR3 <- rep(42L, nchar(sR3)); lR3 <- c( 4L, 16L)
  res3a <- mutscan:::test_mergeReadPairPartial(sF3, qF3, sR3, qR3, lF3, lR3, 1, 0, 1, 0, 0.1, FALSE)
  res3b <- mutscan:::test_mergeReadPairPartial(sF3, qF3, sR3, qR3, c(5L, 13L), 20L, 1, 0, 1, 0, 0.1, FALSE)
  res3c <- mutscan:::test_mergeReadPairPartial(sF3, qF3, sR3, qR3, c(5L, 13L), c(14L, 6L), 1, 0, 1, 0, 0.1, FALSE)
  res3d <- mutscan:::test_mergeReadPairPartial(sF3, qF3, sR3, qR3, c(5L, 8L, 5L), c(4L, 13L, 3L), 1, 0, 1, 0, 0.1, FALSE)
  res3e <- mutscan:::test_mergeReadPairPartial(sF3, qF3, sR3, qR3, c(3L, 2L, 8L, 2L, 3L), c(4L, 2L, 11L, 3L), 1, 0, 1, 0, 0.1, FALSE)
  res3f <- mutscan:::test_mergeReadPairPartial(sF3, qF3, sR3, qR3, c(5L, 13L), c(9L, 5L, 6L), 1, 0, 1, 0, 0.1, FALSE)
  expect_false(res3a$return)
  expect_identical(res3a$mergedSeq, "AAAAAAAAACGTCCCAACCCGGGGGGGGG")
  expect_equal(res3a$mergedLengths, c(13L, 16L))
  expect_identical(res3b$mergedSeq, "AAAAAAAAACGTCCCAACCCGGGGGGGGG")
  expect_equal(res3b$mergedLengths, c(5L, 24L))
  expect_identical(res3c$mergedSeq, "AAAAAAAAACGTCCCAACCCGGGGGGGGG")
  expect_equal(res3c$mergedLengths, c(5L, 18L, 6L))
  expect_identical(res3d$mergedSeq, "AAAAAAAAACGTCCCAACCCGGGGGGGGG")
  expect_equal(res3d$mergedLengths, c(5L, 8L, 13L, 3L))
  expect_identical(res3e$mergedSeq, "AAAAAAAAACGTCCCAACCCGGGGGGGGG")
  expect_equal(res3e$mergedLengths, c(3L, 2L, 8L, 2L, 11L, 3L))
  expect_identical(res3f$mergedSeq, "AAAAAAAAACGTCCCAACCCGGGGGGGGG")
  expect_equal(res3f$mergedLengths, c(5L, 13L, 5L, 6L))
  ## ... that must fail
  res3g <- mutscan:::test_mergeReadPairPartial(sF3, qF3, sR3, qR3, c(12L, 6L), c( 4L, 16L), 1, 0, 1, 0, 0.1, FALSE)
  res3h <- mutscan:::test_mergeReadPairPartial(sF3, qF3, sR3, qR3, c(5L, 13L), c(2L, 18L), 1, 0, 1, 0, 0.1, FALSE)
  res3i <- mutscan:::test_mergeReadPairPartial(sF3, qF3, sR3, qR3, c(12L, 6L), c(3L, 4L, 2L, 11L), 1, 0, 1, 0, 0.1, FALSE)
  expect_identical(res3g$mergedSeq, sF3)
  expect_true(res3g$return)
  expect_identical(res3h$mergedSeq, sF3)
  expect_true(res3h$return)
  expect_identical(res3i$mergedSeq, sF3)
  expect_true(res3i$return)
  ## ... that completely overlaps and must work
  res3j <- mutscan:::test_mergeReadPairPartial(sF3, qF3, sF3, qF3, 18L, 18L, 1, 0, 1, 0, 0.1, FALSE)
  res3k <- mutscan:::test_mergeReadPairPartial(sF3, qF3, sF3, qF3, c(5L, 13L), c(5L, 13L), 1, 0, 1, 0, 0.1, FALSE)
  expect_identical(res3j$mergedSeq, sF3)
  expect_false(res3j$return)
  expect_equal(res3j$mergedLengths, 18L)
  expect_identical(res3k$mergedSeq, sF3)
  expect_false(res3k$return)
  expect_equal(res3k$mergedLengths, c(5L, 13L))
})

## ----------------------------------------------------------------------------
## makeBaseHGVS/makeAAHGVS
## ----------------------------------------------------------------------------
test_that("makeBaseHGVS works", {
    expect_equal(makeBaseHGVS(c("FOS.1.A", "FOS.2.T"), ".", "TAG", "ATG"), "1_2delinsAT_")
    expect_equal(makeBaseHGVS(c("f.1.A", "r.4.A"), ".", "TAGT", "AAGA"), "[1T>A;4T>A]_")
    expect_equal(makeBaseHGVS(c("f.1.A", "r.4.A", "f.6.C"), ".", "TAGTGTAGTCCGT", "AAGAGCAGTCCGT"),
                 "[1T>A;4_6delinsAGC]_")
    expect_equal(makeBaseHGVS(c("r.4.A"), ".", "TAGTGTAGTCCGT", "TAGAGTAGTCCGT"), "4T>A_")
})

test_that("makeAAHGVS works", {
    expect_equal(test_makeAAHGVS("FOS.1.H", ".", "TDTLQAETDQLEDEKSALQTEIANLLKEKEKL"),
                 "(Thr1His)_")
    expect_equal(test_makeAAHGVS(c("f.1.L", "r.4.M"), ".", "TDTLQAETDQLEDEKSALQTEIANLLKEKEKL"),
                 "[(Thr1Leu);(Leu4Met)]_")
    expect_equal(test_makeAAHGVS(c("f.1.L", "f.2.M"), ".", "TDTLQAETDQLEDEKSALQTEIANLLKEKEKL"),
                 "[(Thr1Leu);(Asp2Met)]_")
    expect_equal(test_makeAAHGVS(c("f.1.L", "f.2.M", "f.7.M"), ".", "TDTLQAETDQLEDEKSALQTEIANLLKEKEKL"),
                 "[(Thr1Leu);(Asp2Met);(Glu7Met)]_")
})
