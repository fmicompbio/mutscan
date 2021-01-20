context("digestFastqs - helpers")

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
})

## ----------------------------------------------------------------------------
## decomposeRead
## ----------------------------------------------------------------------------
test_that("decomposeRead works", {
  ## All element lengths specified, no primer
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "USCVCUSV", 
    elementLengths = c(2, 1, 3, 3, 1, 2, 1, 1), primerSeqs = "", 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "ACCT")
  expect_equal(res$varSeq, "AAGG")
  expect_equal(res$varQual, "gABG")
  expect_equal(res$constSeq, "TTAC")
  expect_equal(res$constQual, "defC")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)
  
  ## All element lengths specified, no primer, 0 lengths included
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTG", squal = "abcdefgABCDEF", elements = "SUSCVUCUS", 
    elementLengths = c(0, 2, 1, 3, 3, 0, 1, 2, 1), primerSeqs = "", 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "ACCT")
  expect_equal(res$varSeq, "AAG")
  expect_equal(res$varQual, "gAB")
  expect_equal(res$constSeq, "TTAC")
  expect_equal(res$constQual, "defC")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)
  
  ## Derive variable length from other lengths
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTG", squal = "abcdefgABCDEF", elements = "USCVCUS", 
    elementLengths = c(2, 1, 3, -1, 1, 2, 1), primerSeqs = "", 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "ACCT")
  expect_equal(res$varSeq, "AAG")
  expect_equal(res$varQual, "gAB")
  expect_equal(res$constSeq, "TTAC")
  expect_equal(res$constQual, "defC")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)
  
  ## Primer
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "USCVPS", 
    elementLengths = c(2, 1, 3, 3, 3, -1), primerSeqs = c("CCT"), 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "AC")
  expect_equal(res$varSeq, "AAG")
  expect_equal(res$varQual, "gAB")
  expect_equal(res$constSeq, "TTA")
  expect_equal(res$constQual, "def")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)
  
  ## Only variable sequence
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "V", 
    elementLengths = -1, primerSeqs = c("CCT"), 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "ACGTTAAAGCCTGG")
  expect_equal(res$varQual, "abcdefgABCDEFG")
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)
  
  ## Primer + variable sequence
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "PV", 
    elementLengths = c(-1, -1), primerSeqs = c("ACG"), 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "TTAAAGCCTGG")
  expect_equal(res$varQual, "defgABCDEFG")
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)
  
  ## Variable sequence + primer
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "VP", 
    elementLengths = c(-1, -1), primerSeqs = c("ACG"), 
    varSeq = "A", varQual = "E", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "A")
  expect_equal(res$varQual, "E")
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)
  
  ## Primer, not the right length of sequence
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "VPU", 
    elementLengths = c(10, 3, 2), primerSeqs = c("CCT"), 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "")
  expect_equal(res$varQual, "")
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 1)
  
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "VPU", 
    elementLengths = c(8, 3, 2), primerSeqs = c("CCT"), 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "")
  expect_equal(res$varQual, "")
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 1)
  
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "CVPU", 
    elementLengths = c(10, -1, 3, 2), primerSeqs = c("CCT"), 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "")
  expect_equal(res$varQual, "")
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 1)
  
  ## Primer, but don't explicitly specify that the part after the primer should be skipped
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "USCVP", 
    elementLengths = c(2, 1, 3, 3, 3), primerSeqs = c("CCT"), 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "AC")
  expect_equal(res$varSeq, "AAG")
  expect_equal(res$varQual, "gAB")
  expect_equal(res$constSeq, "TTA")
  expect_equal(res$constQual, "def")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)
  
  ## Multiple primers, only one present
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "USCVP", 
    elementLengths = c(2, 1, 3, 3, 3), primerSeqs = c("TTT", "CCT"), 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "AC")
  expect_equal(res$varSeq, "AAG")
  expect_equal(res$varQual, "gAB")
  expect_equal(res$constSeq, "TTA")
  expect_equal(res$constQual, "def")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)
  
  ## Multiple primers, none present
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "USCVP", 
    elementLengths = c(2, 1, 3, 3, 3), primerSeqs = c("TTT", "GGG"), 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "")
  expect_equal(res$varSeq, "")
  expect_equal(res$varQual, "")
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 1)
  expect_equal(res$nReadWrongLength, 0)
  
  ## Multiple primers, all present (use the first one)
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "USCVP", 
    elementLengths = c(2, 1, 3, 3, 3), primerSeqs = c("CCT", "TTA"), 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "AC")
  expect_equal(res$varSeq, "AAG")
  expect_equal(res$varQual, "gAB")
  expect_equal(res$constSeq, "TTA")
  expect_equal(res$constQual, "def")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)
  
  ## Primer, infer length of sequences
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "UPV", 
    elementLengths = c(-1, 3, -1), primerSeqs = c("CCT"), 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "ACGTTAAAG")
  expect_equal(res$varSeq, "GG")
  expect_equal(res$varQual, "FG")
  expect_equal(res$constSeq, "")
  expect_equal(res$constQual, "")
  expect_equal(res$nNoPrimer, 0)
  expect_equal(res$nReadWrongLength, 0)
  
  ## Primer, don't specify primer length
  res <- mutscan:::test_decomposeRead(
    sseq = "ACGTTAAAGCCTGG", squal = "abcdefgABCDEFG", elements = "UPV", 
    elementLengths = c(-1, -1, -1), primerSeqs = c("CCT"), 
    varSeq = "", varQual = "", umiSeq = "", constSeq = "", constQual = "", 
    nNoPrimer = 0, nReadWrongLength = 0)
  expect_equal(res$umiSeq, "ACGTTAAAG")
  expect_equal(res$varSeq, "GG")
  expect_equal(res$varQual, "FG")
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
  sF1 <- "AAAANA"; qF1 <- rep(40L, nchar(sF1))
  sR1 <- "AACCCC"; qR1 <- rep(42L, nchar(sR1))
  res1a <- mutscan:::test_mergeReadPairPartial(sF1, qF1, sR1, qR1, 1, 6, 1/4, TRUE)
  res1b <- mutscan:::test_mergeReadPairPartial(sF1, qF1, sR1, qR1, 1, 6, 1/4, FALSE)
  res1c <- mutscan:::test_mergeReadPairPartial(sF1, qR1, sR1, qF1, 2, 2, 1/4, TRUE)
  res1d <- mutscan:::test_mergeReadPairPartial(sF1, qF1, sR1, qR1, 2, 2, 0, TRUE)
  expect_is(res1a, "list")
  expect_is(res1b, "list")
  expect_is(res1c, "list")
  expect_is(res1d, "list")
  expect_identical(res1a$mergedSeq, "AAAACCCC")
  expect_identical(res1a$mergedQual, rep(c(40L, 42L), c(2, 6)))
  expect_identical(res1a, res1b)
  expect_identical(res1c$mergedSeq, "AAAANACCCC")
  expect_identical(res1d$mergedSeq, "AAAAAACCCC")
  
  ## no valid overlap, multiple possible overlaps with single valid overlap
  sF2 <- "TTACACG"; qF2 <- rep(10L, nchar(sF2))
  sR2 <- "ACACACA"; qR2 <- rep(40L, nchar(sR2))
  res2a <- mutscan:::test_mergeReadPairPartial(sF2, qF2, sR2, qR2, maxFracMismatchOverlap = 3/7, TRUE)
  res2b <- mutscan:::test_mergeReadPairPartial(sF2, qF2, sR2, qR2, 1, 7, 1/5, TRUE)
  expect_is(res2a, "list")
  expect_is(res2b, "list")
  expect_identical(res2a$mergedSeq, sR2)
  expect_identical(res2a$mergedQual, qR2)
  expect_identical(res2b$mergedSeq, "TTACACACA")
  expect_identical(res2b$mergedQual, rep(c(10L,40L), c(2,7)))
  
  ## padded reads
  for (i in 1:10) {
    sF <- paste(rep(c("C","A"), c(i, 6)), collapse = "")
    qF <- rep(30L, nchar(sF))
    sR <- paste(rep(c("A","C"), c(6, i)), collapse = "")
    qR <- rep(30L, nchar(sR))
    res <- mutscan:::test_mergeReadPairPartial(sF, qF, sR, qR, 6, 6, 0)
    expect_identical(res$mergedSeq, paste(rep(c("C","A","C"), c(i, 6, i)), collapse = ""))
  }
  
  ## reads of unequal length
  sR <- paste(rep("A", 6), collapse = "")
  qR <- rep(30L, nchar(sR))
  for (i in 1:10) {
    sF <- paste(rep(c("C","A"), c(i, 6)), collapse = "")
    qF <- rep(30L, nchar(sF))
    res <- mutscan:::test_mergeReadPairPartial(sF, qF, sR, qR, 6, 6, 0)
    expect_identical(res$mergedSeq, sF)
    res <- mutscan:::test_mergeReadPairPartial(sR, qR, sF, qF, 6, 6, 0)
    expect_identical(res$mergedSeq, sR)
  }
})
