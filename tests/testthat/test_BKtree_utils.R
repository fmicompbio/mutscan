context("BKtree_utils")

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
  
  # load module and create tree instance
  bk <- Rcpp::Module("mod_BKtree", PACKAGE = "mutscan")
  BKtree <- bk$BKtree
  tree <- new(BKtree)
  expect_is(bk, "Module")
  expect_is(BKtree, "C++Class")
  expect_is(tree, "Rcpp_BKtree")
  
  # add sequences
  expect_equal(tree$size, 0)
  for (s in seqs)
    tree$insert(s)
  expect_equal(tree$size, n)
  tree$insert(seqs[1])
  expect_equal(tree$size, n)

  # check presence
  expect_true(tree$has(seqs[1], 0))
  expect_true(tree$has(seqs[1], k))
  expect_false(tree$has("non_existing", 0))
  
  # get first element
  expect_identical(tree$first(), seqs[1])
  
  # memory size
  m <- tree$capacity()
  expect_is(m, "integer")
  expect_true(m > 0 && m < n * k * 8)
  
  # print tree
  expect_output(tree$print(), paste0("^", n, ".+0: ", seqs[1]))

  # search  similar sequences
  res00 <- tree$search(seqs[1], 0)
  res01 <- tree$search(seqs[1], 1)
  res24 <- tree$search(seqs[1], 24)
  expect_is(res00, "character")
  # ground truth e.g. from:
  # sum(stringdist::stringdist(seqs[1], seqs, method = "hamming") <= 24)
  expect_length(res00, 1L)
  expect_length(res01, 5L)
  expect_length(res24, 14L)

  expect_length(tree$search(seqs[2],  2),  1L)
  expect_length(tree$search(seqs[2], 17),  1L)
  expect_length(tree$search(seqs[2], 23),  9L)
  expect_length(tree$search(seqs[2], 26), 16L)
  expect_length(tree$search(seqs[2], 27), 20L)
  expect_length(tree$search(seqs[2], 30), 20L)

  # remove sequences
  tree$remove("non_existing")
  expect_equal(tree$size, n)
  r <- ceiling(n * 0.7)
  for (i in seq_len(r))
    tree$remove(seqs[i])
  expect_equal(tree$size, n - r)
  tree$remove_all()
  expect_equal(tree$size, 0)

  # check existing elements
  tree$remove_all()
  for (s in seqs)
    tree$insert(s)
  res <- tree$get_all()
  expect_identical(seqs, res)
  set.seed(43)
  i <- sample(length(seqs), round(length(seqs) / 3))
  for (ii in i)
    tree$remove(seqs[ii])
  res <- tree$get_all()
  expect_identical(seqs[-i], res)
})

# using Rcpp wrapper function and tree instance global_tree (global variable)
# test_that("low-level BKtree wrapper functions work as expected", {
#   # create random k-mers
#   set.seed(42)
#   k <- 30L
#   n <- 20L
#   alph <- c("A","C","G","T","N")
#   seqs <- unlist(lapply(seq_len(n - length(alph) + 1),
#                         function(i) paste(sample(alph, size = k, replace = TRUE),
#                                           collapse = "")))
#   seqs <- c(seqs, paste0(setdiff(alph, substr(seqs[1], 1, 1)),
#                          substr(seqs[1], 2, k)))
#   
#   # add sequences
#   expect_equal(mutscan:::bk_add(seqs), n)
#   expect_equal(mutscan:::bk_add(seqs[1]), n)
#   expect_equal(mutscan:::bk_new(seqs), n)
#   expect_equal(mutscan:::bk_size(), n)
#   
#   # check presence
#   expect_true(mutscan:::bk_has(seqs[1]))
#   expect_false(mutscan:::bk_has("non_existing"))
# 
#   # print tree
#   expect_output(mutscan:::bk_print(),
#                 paste0("current size: ", n, ".+", n, ".+", seqs[1], ", 0"))
#   # search  similar sequences
#   res00 <- mutscan:::bk_search(seqs[1], tol = 0)
#   res01 <- mutscan:::bk_search(seqs[1], tol = 1)
#   res24 <- mutscan:::bk_search(seqs[1], tol = 24)
#   expect_is(res00, "character")
#   # ground truth e.g. from:
#   # sum(stringdist::stringdist(seqs[1], seqs, method = "hamming") <= 24)
#   expect_length(res00, 1L)
#   expect_length(res01, 5L)
#   expect_length(res24, 14L)
# 
#   # remove sequences
#   expect_equal(mutscan:::bk_remove("non_existing"), n)
#   r <- ceiling(n * 0.7)
#   expect_equal(mutscan:::bk_remove(seqs[seq_len(r)]), n - r)
#   expect_equal(mutscan:::bk_clear(), 0L)
# })


# test_that("collapseSeqs works as expected", {
#   # create random k-mers
#   set.seed(41)
#   k <- 6L
#   n <- 1000L
#   alph <- c("A","C","G","T")
#   seqs <- unlist(lapply(seq_len(n),
#                         function(i) paste(sample(alph, size = k, replace = TRUE),
#                                           collapse = "")))
# 
#   # pre-fligh checks
#   expect_error(collapseSeqs(seqs = 1))
#   expect_error(collapseSeqs(seqs = character(0)))
#   expect_error(collapseSeqs(seqs = c("A","AA")))
#   expect_error(collapseSeqs(seqs = c("AA","AC"), counts = c("AA","AC")))
#   expect_error(collapseSeqs(seqs = c("AA","AC"), counts = 1:3))
#   expect_error(collapseSeqs(seqs = c("AA","AC"), method = "non_existing"))
#   expect_error(collapseSeqs(seqs = c("AA","AC"), tol = -1))
#   expect_error(collapseSeqs(seqs = c("AA","AC"), tol = 1:2))
#   expect_error(collapseSeqs(seqs = c("AA","AC"), verbose = "no"))
#   expect_error(collapseSeqs(seqs = c("AA","AC"), verbose = c(TRUE,FALSE)))
#   
#   # add sequences
#   res1 <- collapseSeqs(seqs, tol = 1)
#   res2 <- collapseSeqs(rev(seqs), seq_along(seqs), tol = 1)
#   res3 <- collapseSeqs(seqs, tol = k)
#   expect_is(res1, "list")
#   expect_length(unlist(res1), length(unique(seqs)))
#   expect_identical(res1, res2)
#   expect_length(res3, 1L)
# })
