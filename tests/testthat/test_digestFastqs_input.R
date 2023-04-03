test_that("digestFastqs fails with incorrect arguments", {
  ## example "trans" fastq files 
  fqt1 <- system.file("extdata/transInput_1.fastq.gz", package = "mutscan")
  fqt2 <- system.file("extdata/transInput_2.fastq.gz", package = "mutscan")
  
  ## default arguments
  Ldef <- list(
    fastqForward = fqt1, fastqReverse = fqt2, 
    mergeForwardReverse = FALSE, 
    minOverlap = 0, maxOverlap = 0, maxFracMismatchOverlap = 0, greedyOverlap = TRUE, 
    revComplForward = FALSE, revComplReverse = FALSE,
    elementsForward = "SUCV", elementsReverse = "SUCV",
    elementLengthsForward = c(1, 10, 18, 96),
    elementLengthsReverse = c(1, 8, 20, 96),
    adapterForward = "GGAAGAGCACACGTC", 
    adapterReverse = "GGAAGAGCGTCGTGT",
    primerForward = c(""),
    primerReverse = c(""),
    wildTypeForward = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
    wildTypeReverse = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT", 
    constantForward = "AACCGGAGGAGGGAGCTG", 
    constantReverse = "GAAAAAGGAAGCTGGAGAGA", 
    avePhredMinForward = 20.0, avePhredMinReverse = 20.0, 
    variableNMaxForward = 0, variableNMaxReverse = 0, 
    umiNMax = 0,
    nbrMutatedCodonsMaxForward = 1,
    nbrMutatedCodonsMaxReverse = 1,
    nbrMutatedBasesMaxForward = -1,
    nbrMutatedBasesMaxReverse = -1,
    forbiddenMutatedCodonsForward = "NNW",
    forbiddenMutatedCodonsReverse = "NNW",
    useTreeWTmatch = TRUE,
    collapseToWTForward = FALSE,
    collapseToWTReverse = FALSE,
    mutatedPhredMinForward = 0.0, 
    mutatedPhredMinReverse = 0.0,
    mutNameDelimiter = ".",
    constantMaxDistForward = -1,
    constantMaxDistReverse = -1,
    variableCollapseMaxDist = 0,
    variableCollapseMinReads = 0,
    variableCollapseMinRatio = 0,
    umiCollapseMaxDist = 0,
    filteredReadsFastqForward = "",
    filteredReadsFastqReverse = "",
    maxNReads = -1, verbose = FALSE,
    nThreads = 1, chunkSize = 1000, 
    maxReadLength = 1024
  )
  
  ## Incorrect specification of merging/rev complementing arguments
  for (var in c("mergeForwardReverse", "revComplForward", "revComplReverse")) {
    L <- Ldef; L[[var]] <- "str"
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- ""
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- c(TRUE, TRUE)
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- c(TRUE, FALSE)
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- 1
    expect_error(do.call(digestFastqs, L))
  }
  
  ## Nonexistent fastqForward/fastqReverse file
  L <- Ldef; L$fastqForward <- "nonexistent.fastq.gz"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqForward <- ""
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqForward <- NULL
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqReverse <- "nonexistent.fastq.gz"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqReverse <- ""
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqForward <- c(fqt1, fqt2)
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L$fastqReverse <- c(fqt1, fqt2)
  expect_error(do.call(digestFastqs, L))
  
  ## No reverse reads if mergeForwardReverse is TRUE
  L <- Ldef; L[["fastqReverse"]] <- NULL; L[["mergeForwardReverse"]] <- TRUE
  expect_error(do.call(digestFastqs, L))
  
  ## Wrong type of numeric argument
  for (var in c("avePhredMinForward", "avePhredMinReverse", "variableNMaxForward",
                "variableNMaxReverse", "umiNMax", 
                "nbrMutatedCodonsMaxForward", "nbrMutatedCodonsMaxReverse", 
                "nbrMutatedBasesMaxForward", "nbrMutatedBasesMaxReverse", 
                "mutatedPhredMinForward", "mutatedPhredMinReverse",
                "variableCollapseMaxDist", "umiCollapseMaxDist", "variableCollapseMinReads",
                "variableCollapseMinRatio",
                "constantMaxDistForward", "constantMaxDistReverse", "maxNReads",
                "nThreads", "chunkSize", "maxReadLength")) {
    L <- Ldef; L[[var]] <- "str"
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- c(1, 2)
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- -2
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- TRUE
    expect_error(do.call(digestFastqs, L))
    if (!(var %in% c("nbrMutatedCodonsMaxForward", "nbrMutatedCodonsMaxReverse", 
                     "nbrMutatedBasesMaxForward", "nbrMutatedBasesMaxReverse",
                     "constantMaxDistForward", "constantMaxDistReverse",
                     "maxNReads", "maxReadLength"))) {
      L <- Ldef; L[[var]] <- -1
      expect_error(do.call(digestFastqs, L))
    }
    if (var %in% c("nThreads", "chunkSize", "maxReadLength")) {
      L <- Ldef; L[[var]] <- 0
      expect_error(do.call(digestFastqs, L))
    }
  }
  
  ## Both or none of max number of codons and bases specified
  ## Note that this should only give an error if a wildtype sequence is specified (which it is here)
  L <- Ldef; L[["nbrMutatedCodonsMaxForward"]] <- L[["nbrMutatedBasesMaxForward"]] <- 1
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["nbrMutatedCodonsMaxReverse"]] <- L[["nbrMutatedBasesMaxReverse"]] <- 1
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["nbrMutatedCodonsMaxForward"]] <- L[["nbrMutatedBasesMaxForward"]] <- -1
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["nbrMutatedCodonsMaxReverse"]] <- L[["nbrMutatedBasesMaxReverse"]] <- -1
  expect_error(do.call(digestFastqs, L))
  
  
  for (var in c("elementLengthsForward", "elementLengthsReverse")) {
    L <- Ldef; L[[var]] <- as.character(L[[var]])
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- as.logical(L[[var]])
    expect_error(do.call(digestFastqs, L))
  }
  
  for (var in c("elementsForward", "elementsReverse")) {
    L <- Ldef; L[[var]] <- "ABC"
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- "UPPV"
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- c("UPV", "CUS")
    expect_error(do.call(digestFastqs, L))
  }
  
  for (var in c("Forward", "Reverse")) {
    L <- Ldef; L[[paste0("elements", var)]] <- "CUV"
    L[[paste0("elementLengths", var)]] <- c(-1, -1, 4)
    expect_error(do.call(digestFastqs, L))
    L[[paste0("elementLengths", var)]] <- c(1, 4)
    expect_error(do.call(digestFastqs, L))
    
    L <- Ldef; L[[paste0("elementLengths", var)]][2] <- -2
    expect_error(do.call(digestFastqs, L))
    
    L <- Ldef; L[[paste0("elements", var)]] <- ""
    L[[paste0("elementLengths", var)]] <- c()
    expect_error(do.call(digestFastqs, L))
    
    L <- Ldef; L[[paste0("elements", var)]] <- "CUVPUV"
    L[[paste0("elementLengths", var)]] <- c(-1, -1, 4, 4, -1, 2)
    expect_error(do.call(digestFastqs, L))
    L[[paste0("elementLengths", var)]] <- c(-1, 1, 4, 4, -1, -1)
    expect_error(do.call(digestFastqs, L))
    L[[paste0("elementLengths", var)]] <- c(1, 4)
    expect_error(do.call(digestFastqs, L))
  }
  
  ## Invalid sequences
  for (var in c("adapterForward", "adapterReverse", "wildTypeForward",
                "wildTypeReverse", "constantForward", "constantReverse")) {
    if (!(var %in% c("constantForward", "constantReverse"))) {
      L <- Ldef; L[[var]] <- c("ACGT", "ACGT")
      expect_error(do.call(digestFastqs, L))
    }
    L <- Ldef; L[[var]] <- 1
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- "EF"
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- paste0(L[[var]], " ")
    expect_error(do.call(digestFastqs, L))
  }
  
  for (var in c("primerForward", "primerReverse")) {
    L <- Ldef; L[[var]] <- 1
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- "EF"
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- c("ACGT", "EF")
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[var]] <- "ACGT "
    expect_error(do.call(digestFastqs, L))
  }
  
  ## Wild type sequence not in named vector (or unnamed string)
  L <- Ldef
  L$wildTypeForward <- c("ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA", 
                         "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT")
  expect_error(do.call(digestFastqs, L))
  
  ## Duplicated wild type sequences
  L <- Ldef
  L$wildTypeForward <- c(wt1 = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA",
                         wt2 = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA")
  expect_error(do.call(digestFastqs, L))
  
  ## Constant sequence of wrong length
  for (var in c("constantForward", "constantReverse")) {
    L <- Ldef; L[[var]] <- substr(L[[var]], 1, 10)
    expect_error(do.call(digestFastqs, L))
  }
  
  ## Invalid forbidden codons
  L <- Ldef; L[["forbiddenMutatedCodonsForward"]] <- c("EFI")
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["forbiddenMutatedCodonsReverse"]] <- c("EFI")
  expect_error(do.call(digestFastqs, L))
  
  ## Invalid mutNameDelimiter
  L <- Ldef; L[["mutNameDelimiter"]] <- ""
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["mutNameDelimiter"]] <- "__"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["mutNameDelimiter"]] <- "f"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["mutNameDelimiter"]] <- "r"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["mutNameDelimiter"]] <- 1
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["mutNameDelimiter"]] <- c("a", "B")
  expect_error(do.call(digestFastqs, L))
  
  ## Invalid value of useTreeWTmatch
  L <- Ldef; L[["useTreeWTmatch"]] <- 2
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["useTreeWTmatch"]] <- "TRUE"
  expect_error(do.call(digestFastqs, L))
  
  ## Invalid value of output FASTQ files
  L <- Ldef; L[["filteredReadsFastqForward"]] <- 1; L[["filteredReadsFastqReverse"]] <- 1
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["filteredReadsFastqForward"]] <- c("file1.gz", "file2.gz")
  L[["filteredReadsFastqReverse"]] <- "file3.gz"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["filteredReadsFastqForward"]] <- TRUE; L[["filteredReadsFastqReverse"]] <- FALSE
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["filteredReadsFastqForward"]] <- "file1.fastq"
  L[["filteredReadsFastqReverse"]] <- "file3.fastq"
  expect_error(do.call(digestFastqs, L))
  
  ## Invalid combination of input and output FASTQ
  L <- Ldef; L[["fastqReverse"]] <- NULL; L[["filteredReadsFastqForward"]] <- "file1.fastq.gz";
  L[["filteredReadsFastqReverse"]] <- "file2.fastq.gz"
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["filteredReadsFastqForward"]] <- "file1.fastq.gz";
  L[["filteredReadsFastqReverse"]] <- ""
  expect_error(do.call(digestFastqs, L))
  L <- Ldef; L[["filteredReadsFastqForward"]] <- "";
  L[["filteredReadsFastqReverse"]] <- "file2.fastq.gz"
  expect_error(do.call(digestFastqs, L))
  
  ## Invalid value of verbose
  for (v in c("verbose", "collapseToWTForward", "collapseToWTReverse")) {
    L <- Ldef; L[[v]] <- 2
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[v]] <- "TRUE"
    expect_error(do.call(digestFastqs, L))
    L <- Ldef; L[[v]] <- c(TRUE, FALSE)
    expect_error(do.call(digestFastqs, L))
  }
  
  ## Too long read
  expect_error(digestFastqs(
      fastqForward = system.file("extdata", "cisInput_1.fastq.gz",
                                 package = "mutscan"),
      elementsForward = "V", elementLengthsForward = -1, 
      maxReadLength = 100),
      "Encountered a read exceeding the maximal allowed length")
  ## Check border cases
  expect_error(digestFastqs(
      fastqForward = system.file("extdata", "cisInput_1.fastq.gz",
                                 package = "mutscan"),
      elementsForward = "V", elementLengthsForward = -1, 
      maxReadLength = 124),
      "Encountered a read exceeding the maximal allowed length")
  ## Works if the maxReadLength is specified accordingly
  out <- digestFastqs(
      fastqForward = system.file("extdata", "cisInput_1.fastq.gz",
                                 package = "mutscan"),
      elementsForward = "V", elementLengthsForward = -1, 
      maxReadLength = 125)
  expect_type(out, "list")
})
