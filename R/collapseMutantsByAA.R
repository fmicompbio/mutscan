## This is from Biostrings
GENETIC_CODE <- structure(c(TTT = "F", TTC = "F", TTA = "L", TTG = "L", TCT = "S", 
                            TCC = "S", TCA = "S", TCG = "S", TAT = "Y", TAC = "Y", TAA = "*", 
                            TAG = "*", TGT = "C", TGC = "C", TGA = "*", TGG = "W", CTT = "L", 
                            CTC = "L", CTA = "L", CTG = "L", CCT = "P", CCC = "P", CCA = "P", 
                            CCG = "P", CAT = "H", CAC = "H", CAA = "Q", CAG = "Q", CGT = "R", 
                            CGC = "R", CGA = "R", CGG = "R", ATT = "I", ATC = "I", ATA = "I", 
                            ATG = "M", ACT = "T", ACC = "T", ACA = "T", ACG = "T", AAT = "N", 
                            AAC = "N", AAA = "K", AAG = "K", AGT = "S", AGC = "S", AGA = "R", 
                            AGG = "R", GTT = "V", GTC = "V", GTA = "V", GTG = "V", GCT = "A", 
                            GCC = "A", GCA = "A", GCG = "A", GAT = "D", GAC = "D", GAA = "E", 
                            GAG = "E", GGT = "G", GGC = "G", GGA = "G", GGG = "G"), alt_init_codons = c("TTG", 
                                                                                                        "CTG"))

collapseMutantsByAA <- function(se) {
  tmp <- base::strsplit(rownames(se), split = "_", fixed = TRUE)
  unl <- unlist(tmp, use.names = FALSE)
  unl <- paste0(substr(unl, 1, nchar(unl) - 3), GENETIC_CODE[substr(unl, nchar(unl) - 2, nchar(unl))])
  tmp <- relist(unl, skeleton = tmp)
  tmp <- S4Vectors::unstrsplit(tmp, sep = "_")
  tmp[rownames(se) == "WT"] <- "WT"

  collapsedCounts <- Matrix.utils::aggregate.Matrix(x = assay(se, "counts"), 
                                                    groupings = factor(tmp), 
                                                    fun = "colSums") 
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = collapsedCounts),
    colData = colData(se), 
    metadata = metadata(se)
  )
}