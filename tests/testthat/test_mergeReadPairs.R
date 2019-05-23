context("mergeReadPairs")

test_that("mergeReadPairs works as expected", {
  read1 <- Biostrings::QualityScaledDNAStringSet(
    x = c("AAGCTNGGC", "NTGCAGGCT"), 
    quality = Biostrings::PhredQuality(
      x = Biostrings::BStringSet(x = c(seqTools::ascii2char(c(17, 24, 37, 37, 32, 15, 37, 26, 12) + 33), 
                                       seqTools::ascii2char(c(28, 22, 12, 37, 37, 17, 32, 19, 12) + 33)))))
  read2 <- Biostrings::QualityScaledDNAStringSet(
    x = c("ATCCTTGGT", "NTNCGGACG"), 
    quality = Biostrings::PhredQuality(
      x = Biostrings::BStringSet(x = c(seqTools::ascii2char(c(28, 22, 12, 37, 37, 17, 32, 19, 12) + 33), 
                                       seqTools::ascii2char(c(17, 24, 37, 37, 32, 15, 37, 26, 12) + 33)))))
  readmerged <- mergeReadPairs(readsForward = read1, readsReverse = read2)
  expect_equal(readmerged, Biostrings::QualityScaledDNAStringSet(
    x = c("AAGCTTGGC", "NTNCAGACT"),
    quality = Biostrings::PhredQuality(
      x = Biostrings::BStringSet(x = c(seqTools::ascii2char(c(28, 24, 37, 37, 37, 17, 37, 26, 12) + 33),
                                       seqTools::ascii2char(c(28, 24, 37, 37, 37, 17, 37, 26, 12) + 33))))
  ))
})
