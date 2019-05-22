context("isValidL")

test_that("isValidL works as expected", {
  s <- QualityScaledDNAStringSet(c("AAAAAAAAAA", "ACAAAAAAAC"),
                                 Biostrings::PhredQuality(c("FFFFFFFFFF", "FBFFFFFFF<")))
  L <- list(umis = s, constantSeqForward = s, constantSeqReverse = s,
            variableSeqForward = s, variableSeqReverse = s, readSummary = data.frame())
  expect_true(isValidL(L))
  expect_error(isValidL(L[-1]))
  expect_error(isValidL(unname(L)))
})


context("findMismatchPositions")

test_that("findMismatchPositions works as expected", {
  patt <- "AAAAAAAAAA"
  subj <- DNAStringSet(c("AAAAAAAAAA", "ACAAAAAAAC", "AACAAAACCA", "AAACAAAAAA", "AAAACAAAAA"))

  expect_error(findMismatchPositions(NULL, subj))
  expect_error(findMismatchPositions(patt, NULL))

  res <- findMismatchPositions(patt, subj)
  expect_is(res, "list")
  expect_length(res, 2L)
  expect_identical(names(res), c("nucleotideMismatches", "codonMismatches"))
  expect_is(res$nucleotideMismatches, "IntegerList")
  expect_length(unlist(res$nucleotideMismatches), 7L)
  expect_identical(sum(unlist(res$nucleotideMismatches)), 41L)
})


context("tabulateQualitiesByMatchstate")

test_that("tabulateQualitiesByMatchstate works as expected", {
  patt <- "AAAAAAAAAA"
  subj <- QualityScaledDNAStringSet(c("AAAAAAAAAA", "ACAAAAAAAC", "AACAAAACCA", "AAACAAAAAA", "AAAACAAAAA"),
                                    Biostrings::PhredQuality(c("FFFFFFFFFF", "FBFFFFFFF<", "FF<FFFF//F", "FFF7FFFFFF", "FFFFFFFFFF")))
  
  expect_error(tabulateQualitiesByMatchstate(NULL, subj))
  expect_error(tabulateQualitiesByMatchstate(patt, NULL))
  
  res <- tabulateQualitiesByMatchstate(patt, subj)
  expect_is(res, "list")
  expect_length(res, 2L)
  expect_identical(names(res), c("error", "correct"))
  expect_length(res$error, 99L)
  expect_identical(sum(res$error), 7L)
  expect_equal(unname(which(res$error > 0)), c(14, 22, 27, 33, 37))
  
  res2 <- findMismatchPositions(patt, subj)
  expect_identical(sum(lengths(res2$nucleotideMismatches)), sum(res$error))
})
