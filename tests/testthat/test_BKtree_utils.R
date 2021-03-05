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




test_that("hamming_shift_distance() works", {
  # ground truth e.g. from:
  # s <- 1; stringdist::stringdist(substr(s1, 1, 15 - s), substr(s2, 1 + s, 15), "hamming") + 2 * s
  # s <- 1; stringdist::stringdist(substr(s2, 1, 15 - s), substr(s1, 1 + s, 15), "hamming") + 2 * s

  s1 <- "AAACAAAACGTTGCC"
  s2 <- "ACAAAACGTTGCCCC"
  s3 <- "AGTCATGCTTAGAAA"
  s4 <- "AAAAGTCATGCTTAG"
  
  expect_error(hamming_shift_distance(c(s1, s1), s2))
  expect_error(hamming_shift_distance(s1, c(s2, s2)))
  expect_error(hamming_shift_distance(s1, s2, "error"))
  
  expect_identical(hamming_shift_distance(s1, s1),  0L)
  expect_identical(hamming_shift_distance(s2, s2),  0L)
  expect_identical(hamming_shift_distance(s3, s3),  0L)
  expect_identical(hamming_shift_distance(s4, s4),  0L)
  expect_identical(hamming_shift_distance(s1, s2,  0),  9L)
  expect_identical(hamming_shift_distance(s1, s2,  1),  9L)
  expect_identical(hamming_shift_distance(s1, s2,  2),  4L)
  expect_identical(hamming_shift_distance(s1, s2,  3),  4L)
  expect_identical(hamming_shift_distance(s1, s2, -1),  4L)
  expect_identical(hamming_shift_distance(s3, s4,  0), 11L)
  expect_identical(hamming_shift_distance(s3, s4,  1), 11L)
  expect_identical(hamming_shift_distance(s3, s4,  2), 11L)
  expect_identical(hamming_shift_distance(s3, s4,  3),  6L)
  expect_identical(hamming_shift_distance(s3, s4, -1),  6L)
})


# using Rcpp-module that exposes the BKtree class to R
test_that("low-level BKtree wrapper functions work as expected", {
  # create random k-mers
  set.seed(42)
  k <- 30L
  n <- 20
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
  expect_error(new(BKtree, seqs, "error"))
  expect_error(new(BKtree, seqs, "hamming_shift", "error"))
  tree <- new(BKtree)
  tree2 <- new(BKtree, "levenshtein")
  tree3 <- new(BKtree, seqs[c(seq.int(n), seq.int(n))], "hamming", -1)
  tree4 <- new(BKtree, character(0), "hamming", -1)
  tree5 <- new(BKtree, "hamming_shift", 5L)
  tree6 <- new(BKtree, "hamming_shift", 0L)
  
  expect_is(bk, "Module")
  expect_is(BKtree, "C++Class")
  expect_is(tree, "Rcpp_BKtree")
  expect_is(tree2, "Rcpp_BKtree")
  expect_is(tree3, "Rcpp_BKtree")
  expect_is(tree4, "Rcpp_BKtree")
  expect_is(tree5, "Rcpp_BKtree")
  expect_is(tree6, "Rcpp_BKtree")
  expect_identical(tree$distance_metric(), "hamming")
  expect_identical(tree2$distance_metric(), "levenshtein")
  expect_identical(tree3$distance_metric(), "hamming")
  expect_identical(tree4$distance_metric(), "hamming")
  expect_identical(tree5$distance_metric(), "hamming_shift")
  expect_identical(tree6$distance_metric(), "hamming_shift")
  expect_identical(tree5$maximal_absolute_shift(), 5L)
  expect_identical(tree6$maximal_absolute_shift(), 0L)
  
  # add sequences
  expect_error(tree$insert(seqs))
  expect_identical(tree$size, 0)
  expect_identical(tree2$size, 0)
  for (s in seqs) {
    tree$insert(s)
    tree2$insert(s)
  }
  expect_identical(tree$size, n)
  expect_identical(tree2$size, n)
  expect_identical(tree3$size, n)
  tree$insert(seqs[1])
  tree2$insert(seqs[1])
  expect_identical(tree$size, n)
  expect_identical(tree2$size, n)
  
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
  expect_is(m <- tree$capacity(), "integer")
  expect_true(m > 0 && m < n * k * 8)
  expect_identical(tree2$capacity(), m)
  expect_identical(tree3$remove(seqs[n]), NULL)
  expect_is(m3 <- tree3$capacity(), "integer")
  expect_true(m3 > 0 && m3 < (n - 1) * k * 8)
  expect_identical(tree3$insert(seqs[n]), NULL)
  
  # print tree
  expect_output(tree$print(), paste0("^", n, ".+0: ", seqs[1]))
  expect_output(tree2$print(), paste0("^", n, ".+0: ", seqs[1]))
  expect_identical(tree3$remove(seqs[n]), NULL)
  expect_output(tree3$print(), paste0("^", n - 1, ".+0: ", seqs[1]))
  expect_identical(tree3$insert(seqs[n]), NULL)
  expect_output(tree4$print(), "^0$")
  
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
  # ... hamming+shift distances
  expect_identical(tree4$size, 0)
  expect_identical(tree5$size, 0)
  expect_identical(tree6$size, 0)
  for (s in c("CGATCGATGCAA", "AACGATCGATGC")) {
    tree4$insert(s)
    tree5$insert(s)
    tree6$insert(s)
  }
  expect_identical(tree4$size, 2)
  expect_identical(tree5$size, 2)
  expect_identical(tree6$size, 2)
  expect_length(tree4$search("CGATCGATGCAA",  0), 1L)
  expect_length(tree4$search("CGATCGATGCAA", 11), 1L)
  expect_length(tree4$search("CGATCGATGCAA", 12), 2L)
  expect_length(tree5$search("CGATCGATGCAA",  0), 1L)
  expect_length(tree5$search("CGATCGATGCAA",  3), 1L)
  expect_length(tree5$search("CGATCGATGCAA",  4), 2L)
  expect_length(tree6$search("CGATCGATGCAA",  0), 1L)
  expect_length(tree6$search("CGATCGATGCAA", 11), 1L)
  expect_length(tree6$search("CGATCGATGCAA", 12), 2L)
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
  expect_identical(tree$size, n)
  expect_identical(tree2$size, n)
  r <- ceiling(n * 0.7)
  for (i in seq_len(r)) {
    tree$remove(seqs[i])
    tree2$remove(seqs[i])
  }
  expect_identical(tree$size, n - r)
  expect_identical(tree2$size, n - r)
  expect_identical(tree$remove_all(), NULL)
  expect_identical(tree2$remove_all(), NULL)
  expect_identical(tree$size, 0)
  expect_identical(tree2$size, 0)
  expect_identical(tree3$size, n)
  expect_is(tree3$remove(seqs[1]), "NULL")
  expect_identical(tree3$size, n - 1)
  expect_is(tree3$insert(seqs[1]), "NULL")
  expect_identical(tree3$size, n)
  
  # check existing elements
  expect_identical(tree$remove_all(), NULL)
  expect_identical(tree2$remove_all(), NULL)
  expect_identical(tree$search(seqs[1], 0), character(0))
  expect_identical(tree2$has(seqs[1], 0), FALSE)
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
