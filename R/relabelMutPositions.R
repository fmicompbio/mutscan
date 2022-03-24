#' Relabel the positions of mutations in the designated ID
#' 
#' @param se SummarizedExperiment object, with row names of the form
#'     XX{.}AA{.}NNN, where XX is the name of the reference sequence, AA is the
#'     position of the mutated codon, and NNN is the mutated codon or amino 
#'     acid. {.} is the delimiter, to be specified in the 
#'     \code{mutNameDelimiter} argument. For rows corresponding to sequences 
#'     with multiple mutated codons, the row names contain multiple names of 
#'     the form above in a single string, separated by "_".
#' @param conversionTable \code{data.frame} with at least three columns:
#'     \itemize{
#'     \item seqname The reference sequence name (should match XX in the 
#'     mutation name)
#'     \item position The codon position (should match AA in the mutation name)
#'     \item name The new name for the codon (will replace AA in the mutation 
#'     name, if the reference sequence matches seqname)
#'     } 
#' @param mutNameDelimiter The delimiter used in the mutation name ({.} above).
#' 
#' @author Charlotte Soneson
#' 
#' @return A SummarizedExperiment object with modified row names.
#' 
#' @importFrom S4Vectors unstrsplit
#' @importFrom utils relist
#' @importFrom BiocGenerics rownames
#' 
#' @examples
#' x <- readRDS(system.file("extdata", "GSE102901_cis_se.rds",
#'                          package = "mutscan"))
#' conversionTable <- data.frame(seqname = "f", position = 0:32) 
#' conversionTable$name = paste0((conversionTable$position - 1) %/% 7 + 1, 
#'                               c("", rep(letters[1:7], 6))[1:33])
#' out <- relabelMutPositions(x, conversionTable)
#' 
relabelMutPositions <- function(se, conversionTable, mutNameDelimiter = ".") {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertVector(x = colnames(conversionTable), type = "character")
    .assertScalar(x = mutNameDelimiter, type = "character")
    stopifnot(all(c("position", "name", "seqname") %in% colnames(conversionTable)))
    
    conversionTable$position <- as.character(conversionTable$position)
    spl <- base::strsplit(rownames(se), "_")
    unl <- base::unlist(spl)
    unl <- lapply(base::strsplit(unl, mutNameDelimiter, fixed = TRUE), 
                  function(w) {
                      idx <- which(conversionTable$seqname == w[1] & 
                                       conversionTable$position == w[2])
                      if (length(idx) > 0) {
                          w[2] <- conversionTable$name[idx]
                      }
                      w
                  })
    unl <- S4Vectors::unstrsplit(unl, sep = mutNameDelimiter)
    spl <- utils::relist(unl, skeleton = spl)
    spl <- S4Vectors::unstrsplit(spl, "_")
    BiocGenerics::rownames(se) <- spl
    se
}
