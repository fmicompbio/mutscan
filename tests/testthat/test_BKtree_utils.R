context("BKtree_utils")

test_that("levenshtein_distance() works", {
  # ground truth e.g. from:
  # stringdist::stringdist(s1, s2, method = "lv")

  s1 <- "AAAAAAACCC"
  s2 <- "AAAAAACCCC"
  s3 <- "ACGTACGTACGT"
  s4 <- "CGTACGTACGTA"
  
  expect_error(levenshtein_distance(c(s1, s1), s2))
  expect_error(levenshtein_distance(s1, c(s2, s2)))

  expect_identical(levenshtein_distance(s1, s1), 0L)
  expect_identical(levenshtein_distance(s2, s2), 0L)
  expect_identical(levenshtein_distance(s3, s3), 0L)
  expect_identical(levenshtein_distance(s4, s4), 0L)
  expect_identical(levenshtein_distance(s1, s2), 1L)
  expect_identical(levenshtein_distance(s1, s3), 8L)
  expect_identical(levenshtein_distance(s1, s4), 9L)
  expect_identical(levenshtein_distance(s2, s3), 9L)
  expect_identical(levenshtein_distance(s2, s4), 9L)
  expect_identical(levenshtein_distance(s3, s4), 2L)
})


test_that("hamming_distance() works", {
  # ground truth e.g. from:
  # stringdist::stringdist(s1, s2, method = "hamming")
  
  s1 <- "AAAAAAACCC"
  s2 <- "AAAAAACCCC"
  s3 <- "ACGTACGTACGT"
  s4 <- "CGTACGTACGTA"
  
  expect_error(hamming_distance(c(s1, s1), s2))
  expect_error(hamming_distance(s1, c(s2, s2)))
  
  expect_identical(hamming_distance(s1, s1),  0L)
  expect_identical(hamming_distance(s2, s2),  0L)
  expect_identical(hamming_distance(s3, s3),  0L)
  expect_identical(hamming_distance(s4, s4),  0L)
  expect_identical(hamming_distance(s1, s2),  1L)
  expect_identical(hamming_distance(s2, s1),  1L)
  expect_identical(hamming_distance(s3, s4), 12L)
  expect_identical(hamming_distance(s4, s3), 12L)
})


# using Rcpp-module that exposes the BKtree class to R
test_that("low-level BKtree wrapper functions work as expected", {
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
  
  # load module and create tree instances
  bk <- Rcpp::Module("mod_BKtree", PACKAGE = "mutscan")
  BKtree <- bk$BKtree
  expect_error(new(BKtree, "error"))
  tree <- new(BKtree)
  tree2 <- new(BKtree, "levenshtein")
  expect_is(bk, "Module")
  expect_is(BKtree, "C++Class")
  expect_is(tree, "Rcpp_BKtree")
  expect_is(tree2, "Rcpp_BKtree")
  expect_identical(tree$distance_metric(), "hamming")
  expect_identical(tree2$distance_metric(), "levenshtein")
  
  # add sequences
  expect_equal(tree$size, 0)
  expect_equal(tree2$size, 0)
  for (s in seqs) {
    tree$insert(s)
    tree2$insert(s)
  }
  expect_equal(tree$size, n)
  expect_equal(tree2$size, n)
  tree$insert(seqs[1])
  tree2$insert(seqs[1])
  expect_equal(tree$size, n)
  expect_equal(tree2$size, n)
  
  # check presence
  expect_true(tree$has(seqs[1], 0))
  expect_true(tree2$has(seqs[1], 0))
  expect_true(tree$has(seqs[1], k))
  expect_true(tree2$has(seqs[1], k))
  expect_false(tree$has("non_existing", 0))
  expect_false(tree2$has("non_existing", 0))
  
  # get first element
  expect_identical(tree$first(), seqs[1])
  expect_identical(tree2$first(), seqs[1])
  
  # memory size
  m <- tree$capacity()
  expect_is(m, "integer")
  expect_true(m > 0 && m < n * k * 8)
  expect_identical(tree2$capacity(), m)
  
  # print tree
  expect_output(tree$print(), paste0("^", n, ".+0: ", seqs[1]))
  expect_output(tree2$print(), paste0("^", n, ".+0: ", seqs[1]))
  
  # search  similar sequences
  # ... hamming distances
  # ground truth from: sum(stringdist::stringdist(seqs[1], seqs, method = "hamming") <= 24)
  res <- tree$search(seqs[1], 0)
  expect_is(res, "character")
  expect_length(res, 1L)
  expect_length(tree$search(seqs[1],  1),  5L)
  expect_length(tree$search(seqs[1], 24), 14L)
  expect_length(tree$search(seqs[2],  2),  1L)
  expect_length(tree$search(seqs[2], 17),  1L)
  expect_length(tree$search(seqs[2], 23),  9L)
  expect_length(tree$search(seqs[2], 26), 16L)
  expect_length(tree$search(seqs[2], 27), 20L)
  expect_length(tree$search(seqs[2], 30), 20L)
  expect_true(all(tree$search(seqs[2], 23) %in% seqs))
  # ... levenshtein distances
  # ground truth from: sum(stringdist::stringdist(seqs[1], seqs, method = "lv") <= 24)
  res <- tree2$search(seqs[1], 0)
  expect_is(res, "character")
  expect_length(res, 1L)
  expect_length(tree2$search(seqs[1],  1),  5L)
  expect_length(tree2$search(seqs[1], 18),  8L)
  expect_length(tree2$search(seqs[2],  2),  1L)
  expect_length(tree2$search(seqs[2], 17),  3L)
  expect_length(tree2$search(seqs[2], 18),  6L)
  expect_length(tree2$search(seqs[2], 19), 13L)
  expect_length(tree2$search(seqs[2], 20), 18L)
  expect_length(tree2$search(seqs[2], 30), 20L)
  expect_true(all(tree2$search(seqs[2], 19) %in% seqs))
  
  # remove sequences
  expect_identical(tree$remove("non_existing"), NULL)
  expect_identical(tree2$remove("non_existing"), NULL)
  expect_equal(tree$size, n)
  expect_equal(tree2$size, n)
  r <- ceiling(n * 0.7)
  for (i in seq_len(r)) {
    tree$remove(seqs[i])
    tree2$remove(seqs[i])
  }
  expect_equal(tree$size, n - r)
  expect_equal(tree2$size, n - r)
  expect_identical(tree$remove_all(), NULL)
  expect_identical(tree2$remove_all(), NULL)
  expect_equal(tree$size, 0)
  expect_equal(tree2$size, 0)
  
  # check existing elements
  expect_identical(tree$remove_all(), NULL)
  expect_identical(tree2$remove_all(), NULL)
  for (s in seqs) {
    tree$insert(s)
    tree2$insert(s)
  }
  expect_is(res <- tree$get_all(), "character")
  expect_is(res2 <- tree2$get_all(), "character")
  expect_identical(seqs, res)
  expect_identical(seqs, res2)
  set.seed(43)
  i <- sample(length(seqs), round(length(seqs) / 3))
  for (ii in i) {
    tree$remove(seqs[ii])
    tree2$remove(seqs[ii])
  }
  expect_is(res <- tree$get_all(), "character")
  expect_is(res2 <- tree2$get_all(), "character")
  expect_identical(seqs[-i], res)
  expect_identical(seqs[-i], res2)
})
