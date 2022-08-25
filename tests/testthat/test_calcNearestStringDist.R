test_that("calcNearestStringDist fails with incorrect arguments", {
    expect_error(calcNearestStringDist(), "argument \"x\" is missing")
    expect_error(calcNearestStringDist(x = 1:3), "Expecting a string vector")
    expect_error(calcNearestStringDist(x = letters, metric = "error"), "unknown distance metric 'error'")
    expect_error(calcNearestStringDist(x = letters, nThreads = c(1L, 2L)), "Expecting a single value")
})

test_that("calcNearestStringDist works as expected", {
    # strs1 <- c("lazy", "hazy", "craz")  # same lengths
    # strs2 <- c("lazy", "hazy", "crazy") # unequal lengths

    set.seed(42)
    strs1 <- sample(x = Biostrings::mkAllStrings(alphabet = Biostrings::DNA_BASES,
                                                 width = 7),
                    size = 200)
    strs2 <- c(sample(x = Biostrings::mkAllStrings(alphabet = Biostrings::DNA_BASES,
                                                   width = 7),
                      size = 200),
               sample(x = Biostrings::mkAllStrings(alphabet = Biostrings::DNA_BASES,
                                                   width = 8),
                      size = 200))

    d0 <- Biostrings::stringDist(x = strs1, method = "hamming")
    e0 <- Biostrings::stringDist(x = strs2, method = "levenshtein")

    dm <- as.matrix(d0)
    em <- as.matrix(e0)

    diag(dm) <- diag(em) <- NA

    d <- as.integer(unname(apply(dm, 1, min, na.rm = TRUE)))
    e <- as.integer(unname(apply(em, 1, min, na.rm = TRUE)))

    expect_type(res1 <- calcNearestStringDist(x = strs1, metric = "hamming", nThreads = 1L), "integer")
    expect_type(res2 <- calcNearestStringDist(x = strs1, metric = "hamming", nThreads = 4L), "integer")
    expect_type(res3 <- calcNearestStringDist(x = strs2, metric = "levenshtein", nThreads = 1L), "integer")
    expect_type(res4 <- calcNearestStringDist(x = strs2, metric = "levenshtein", nThreads = 4L), "integer")
    expect_type(res5 <- calcNearestStringDist(x = strs1, metric = "hamming_shift", nThreads = 1L), "integer")
    expect_type(res6 <- calcNearestStringDist(x = strs1, metric = "hamming_shift", nThreads = 4L), "integer")

    expect_identical(d, res1)
    expect_identical(e, res3)
    expect_true(all(res5 <= res1))
    expect_identical(res1, res2)
    expect_identical(res3, res4)
    expect_identical(res5, res6)
})




