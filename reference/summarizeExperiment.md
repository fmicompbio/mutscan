# Summarize and collapse multiple mutational scanning experiments

Combine multiple sequence lists (as returned by
[`digestFastqs`](https://fmicompbio.github.io/mutscan/reference/digestFastqs.md)
into a
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html),
with observed variable sequences (sequence pairs) in rows and samples in
columns.

## Usage

``` r
summarizeExperiment(x, coldata, countType = "umis")
```

## Arguments

- x:

  A named list of objects returned by
  [`digestFastqs`](https://fmicompbio.github.io/mutscan/reference/digestFastqs.md).
  Names are used to link the objects to the metadata provided in
  `coldata`.

- coldata:

  A `data.frame` with at least one column "Name", which will be used to
  link to objects in `x`. A potentially subset and reordered version of
  `coldata` is stored in the `colData` of the returned
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).

- countType:

  Either "reads" or "umis". If "reads", the "count" assay of the
  returned object will contain the observed number of reads for each
  sequence (pair). If "umis", the "count" assay will contain the number
  of unique UMIs observed for each sequence (pair).

## Value

A
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
`x` with

- assays(x)\$counts:

  containing the observed number of sequences or sequence pairs (if
  `countType` = "reads"), or the observed number of unique UMIs for each
  sequence or sequence pair (if `countType` = "umis").

- rowData(x):

  containing the unique sequences or sequence pairs.

- colData(x):

  containing the metadata provided by `coldata`.

## Author

Michael Stadler, Charlotte Soneson

## Examples

``` r
## Input sample
inp <- digestFastqs(
    fastqForward = system.file("extdata", "cisInput_1.fastq.gz", 
                               package = "mutscan"), 
    elementsForward = "SUCV", elementLengthsForward = c(1, 10, 18, 96), 
    constantForward = "AACCGGAGGAGGGAGCTG", 
    wildTypeForward = c(FOS = paste0(
        "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTC", 
        "TGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA")),
    nbrMutatedCodonsMaxForward = 1
)
## Output sample
outp <- digestFastqs(
    fastqForward = system.file("extdata", "cisOutput_1.fastq.gz", 
                               package = "mutscan"), 
    elementsForward = "SUCV", elementLengthsForward = c(1, 10, 18, 96), 
    constantForward = "AACCGGAGGAGGGAGCTG", 
    wildTypeForward = c(FOS = paste0(
        "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTC", 
        "TGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA")),
    nbrMutatedCodonsMaxForward = 1
)
## Combine
se <- summarizeExperiment(
    x = list(r1inp = inp, r1outp = outp), 
    coldata = data.frame(Name = c("r1inp", "r1outp"), 
                         Condition = c("input", "output"), 
                         Replicate = c("rep1", "rep1")),
    countType = "umis"
)
se
#> class: SummarizedExperiment 
#> dim: 134 2 
#> metadata(4): parameters errorStatistics countType mutNameDelimiter
#> assays(1): counts
#> rownames(134): FOS.0.WT FOS.1.AAT ... FOS.9.GAA FOS.9.GAG
#> rowData names(19): mutantName sequence ... mutationTypes varLengths
#> colnames(2): r1inp r1outp
#> colData names(19): Name Condition ... f13_nbrTooManyBestConstantHits
#>   nbrRetained
```
