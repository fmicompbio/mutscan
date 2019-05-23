context("mergeReadPairs")

test_that("mergeReadPairs works as expected", {
  read1 <- QualityScaledDNAStringSet(
    x = c("AAGCTNGGC", "NTGCAGGCT"), 
    quality = PhredQuality(x = BStringSet(x = c(ascii2char(c(17, 24, 37, 37, 32, 15, 37, 26, 12) + 33), 
                                                ascii2char(c(28, 22, 12, 37, 37, 17, 32, 19, 12) + 33)))))
  read2 <- QualityScaledDNAStringSet(
    x = c("ATCCTTGGT", "NTNCGGACG"), 
    quality = PhredQuality(x = BStringSet(x = c(ascii2char(c(28, 22, 12, 37, 37, 17, 32, 19, 12) + 33), 
                                                ascii2char(c(17, 24, 37, 37, 32, 15, 37, 26, 12) + 33)))))
  readmerged <- mergeReadPairs(readsForward = read1, readsReverse = read2)
  expect_equal(readmerged, QualityScaledDNAStringSet(
    x = c("AAGCTTGGC", "NTNCAGACT"),
    quality = PhredQuality(x = BStringSet(x = c(ascii2char(c(28, 24, 37, 37, 37, 17, 37, 26, 12) + 33),
                                                ascii2char(c(28, 24, 37, 37, 37, 17, 37, 26, 12) + 33))))
  ))
})
