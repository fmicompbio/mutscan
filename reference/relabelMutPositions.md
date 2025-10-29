# Relabel the positions of mutations in the designated ID

Relabel the positions of mutations in the designated ID

## Usage

``` r
relabelMutPositions(se, conversionTable, mutNameDelimiter = ".")
```

## Arguments

- se:

  SummarizedExperiment object, with row names of the form XX{.}AA{.}NNN,
  where XX is the name of the reference sequence, AA is the position of
  the mutated codon, and NNN is the mutated codon or amino acid. {.} is
  the delimiter, to be specified in the `mutNameDelimiter` argument. For
  rows corresponding to sequences with multiple mutated codons, the row
  names contain multiple names of the form above in a single string,
  separated by "\_".

- conversionTable:

  `data.frame` with at least three columns:

  - seqname The reference sequence name (should match XX in the mutation
    name)

  - position The codon position (should match AA in the mutation name)

  - name The new name for the codon (will replace AA in the mutation
    name, if the reference sequence matches seqname)

- mutNameDelimiter:

  The delimiter used in the mutation name ({.} above).

## Value

A SummarizedExperiment object with modified row names.

## Author

Charlotte Soneson

## Examples

``` r
x <- readRDS(system.file("extdata", "GSE102901_cis_se.rds",
                         package = "mutscan"))
conversionTable <- data.frame(seqname = "f", position = 0:32) 
conversionTable$name = paste0((conversionTable$position - 1) %/% 7 + 1, 
                              c("", rep(letters[1:7], 6))[1:33])
out <- relabelMutPositions(x, conversionTable)
```
