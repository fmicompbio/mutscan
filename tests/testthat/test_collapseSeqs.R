context("collapse sequences by similarity")

test_that("low-level BK tree functions work as expected", {
  # create random k-mers
  set.seed(42)
  k <- 30L
  n <- 20L
  alph <- c("A","C","G","T","N")
  seqs <- unlist(lapply(seq_len(n - length(alph) + 1),
                        function(i) paste(sample(alph, size = k, replace = TRUE),
                                          collapse = "")))
  seqs <- c(seqs, paste0(setdiff(alph, substr(seqs[1], 1, 1)),
                         substr(seqs[1], 2, k)))
  
  # add sequences
  expect_equal(mutscan:::bk_add(seqs), n)
  expect_equal(mutscan:::bk_add(seqs[1]), n)
  expect_equal(mutscan:::bk_new(seqs), n)
  expect_equal(mutscan:::bk_size(), n)
  
  # check presence
  expect_true(mutscan:::bk_has(seqs[1]))
  expect_false(mutscan:::bk_has("non_existing"))

  # print tree
  expect_output(mutscan:::bk_print(),
                paste0("current size: ", n, ".+", n, ".+", seqs[1], ", 0"))
  # search  similar sequences
  res00 <- mutscan:::bk_search(seqs[1], tol = 0)
  res01 <- mutscan:::bk_search(seqs[1], tol = 1)
  res24 <- mutscan:::bk_search(seqs[1], tol = 24)
  expect_is(res00, "character")
  # ground truth e.g. from:
  # sum(stringdist::stringdist(seqs[1], seqs, method = "hamming") <= 24)
  expect_length(res00, 1L)
  expect_length(res01, 5L)
  expect_length(res24, 14L)

  # remove sequences
  expect_equal(mutscan:::bk_remove("non_existing"), n)
  r <- ceiling(n * 0.7)
  expect_equal(mutscan:::bk_remove(seqs[seq_len(r)]), n - r)
  expect_equal(mutscan:::bk_clear(), 0L)
})


test_that("collapseSeqs works as expected", {
  # create random k-mers
  set.seed(41)
  k <- 6L
  n <- 1000L
  alph <- c("A","C","G","T")
  seqs <- unlist(lapply(seq_len(n),
                        function(i) paste(sample(alph, size = k, replace = TRUE),
                                          collapse = "")))

  # pre-fligh checks
  expect_error(collapseSeqs(seqs = 1))
  expect_error(collapseSeqs(seqs = character(0)))
  expect_error(collapseSeqs(seqs = c("A","AA")))
  expect_error(collapseSeqs(seqs = c("AA","AC"), counts = c("AA","AC")))
  expect_error(collapseSeqs(seqs = c("AA","AC"), counts = 1:3))
  expect_error(collapseSeqs(seqs = c("AA","AC"), method = "non_existing"))
  expect_error(collapseSeqs(seqs = c("AA","AC"), tol = -1))
  expect_error(collapseSeqs(seqs = c("AA","AC"), tol = 1:2))
  expect_error(collapseSeqs(seqs = c("AA","AC"), verbose = "no"))
  expect_error(collapseSeqs(seqs = c("AA","AC"), verbose = c(TRUE,FALSE)))
  
  # add sequences
  res1 <- collapseSeqs(seqs, tol = 1)
  res2 <- collapseSeqs(rev(seqs), seq_along(seqs), tol = 1)
  res3 <- collapseSeqs(seqs, tol = k)
  expect_is(res1, "list")
  expect_length(unlist(res1), length(unique(seqs)))
  expect_identical(res1, res2)
  expect_length(res3, 1L)
})
