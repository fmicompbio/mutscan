# Process an experiment with multiple variable sequences

This function enables the processing of data sets with multiple variable
sequences, which should potentially be handled in different ways. For
example, a barcode association experiment with two variable sequences
(the barcode and the biological variant) that need to be processed
differently, e.g. in terms of matching to wildtype sequences or
collapsing of similar sequences. In contrast, while `digestFastqs` allow
the specification of multiple variable sequences (within each of the
forward and reverse reads), they will be concatenated and processed as a
single unit.

## Usage

``` r
linkMultipleVariants(combinedDigestParams = list(), ...)
```

## Arguments

- combinedDigestParams:

  A named list of arguments to `digestFastqs` for the combined ("naive")
  run.

- ...:

  Additional arguments providing arguments to `digestFastqs` for the
  separate runs (processing each variable sequence in turn). Each
  argument must be a named list of arguments to `digestFastqs`. In
  addition, arguments `collapseMaxDist`, `collapseMinScore` and
  `collapseMinRatio` can be specified, and will be passed on to
  `collapseMutantsBySimilarity`.

## Value

A list with the following elements:

- countAggregated - a `tibble` with columns corresponding to each of the
  variable sequences, and a column with the total observed read count
  for the combination.

- convSeparate - a list of conversion tables from the respective
  separate runs.

- outCombined - the `digestFastqs` output for the combined run.

## Details

linkMultipleVariants will process the input in the following way:

- First, run `digestFastqs` with the parameters provided in
  `combinedDigestParams`. Typically, this will be a "naive" counting
  run, where the frequencies of all observed variants are tabulated. The
  variable sequences within the forward and reverse reads, respectively,
  will be processed as a single sequence.

- Next, run `digestFastqs` with each of the additional parameter sets
  provided (`...`). Each of these should correspond to a single variable
  sequence from the combined run (i.e., if there are two Vs in the
  element specifications in the combined run, there should be two
  additional parameter sets provided, each corresponding to the
  processing of one variable sequence part). It is assumed that the
  order of the additional arguments correspond to the order of the
  variable sequences in the combined run, in such a way that if the
  variable sequences extracted in each of the separate runs are
  concatenated in the order that the parameter sets are provided to
  `linkMultipleVariants`, they will form the variable sequence extracted
  in the combined run.

- The result of each of the separate runs is a 'conversion table',
  containing the final set of identified sequence variants as well as
  all individual sequences corresponding to each of them. This is then
  combined with the count table from the combined, "naive" run in order
  to create an aggregated count table. More precisely, each sequence in
  the combined run is split into the constituent variable sequences, and
  each variable sequence is then matched to the output from the right
  separate run, from which the final feature ID (mutant name, or
  collapsed sequence) will be extracted and used to replace the original
  sequence in the combined count table. Once all the matches are done,
  rows with NAs (where no match could be found in the separate run) are
  removed and the counts are aggregated across all identical
  combinations of variable sequences.

In order to define the `elementsForward` and `elementsReverse` arguments
for the separate runs, a strategy that often works is to simply copy the
arguments from the combined run, and successively replace all but one of
the 'V's by 'S'. This will effectively process one variable sequence at
the time, while keeping all other elements of the reads consistent
(since this can affect e.g. filtering criteria). Note that to process
individual variable sequences in the reverse read, you also need to swap
the 'forward' and 'reverse' specifications (since `digestFastqs`
requires a forward read).

## Author

Charlotte Soneson, Michael Stadler

## Examples

``` r
fqFile <- system.file("extdata", "cisInput_1.fastq.gz", 
                      package = "mutscan")
out <- linkMultipleVariants(
    combinedDigestParams = list(fastqForward = fqFile, 
                                elementsForward = "SVCV", 
                                elementLengthsForward = c(1, 10, 18, 96)),
    # the first variable sequence is the UMI
    umi = list(fastqForward = fqFile, elementsForward = "SVCS",
               elementLengthsForward = c(1, 10, 18, 96)),
    # the second variable sequence is the amplicon variant
    var = list(fastqForward = fqFile, elementsForward = "SSCV",
               elementLengthsForward = c(1, 10, 18, 96), 
               collapseMaxDist = 3, collapseMinScore = 1)
)
# conversion tables
lapply(out$convSeparate, head)
#> $umi
#> # A tibble: 6 × 2
#>   mutantName sequence  
#>   <chr>      <chr>     
#> 1 AAAAATATTT AAAAATATTT
#> 2 AAAAGCCTCG AAAAGCCTCG
#> 3 AAACATCAAC AAACATCAAC
#> 4 AAACCGGAGG AAACCGGAGG
#> 5 AAACTACCGG AAACTACCGG
#> 6 AAATAAAGAC AAATAAAGAC
#> 
#> $var
#> # A tibble: 6 × 2
#>   mutantName                                                            sequence
#>   <chr>                                                                 <chr>   
#> 1 AAATATAACGTTGACGATGTAGCTTTAGGTGTCTGTAAAACAGGTGCCGAAGAAGCTGGAGTAACAGA… AAATATA…
#> 2 AAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAA… AAGAGCA…
#> 3 AAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAA… AAGAGCA…
#> 4 AATGATACACTCCAAACGGAGACAGACCAACTAGAATATGAGAAGTCTGCTTTGCAGACCGAGATTGC… AATGATA…
#> 5 ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGC… AATGATA…
#> 6 ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGC… AATGATA…
#> 
# aggregated count table
head(out$countAggregated)
#> # A tibble: 6 × 3
#>   umi        var                                                        nbrReads
#>   <chr>      <chr>                                                         <int>
#> 1 AAACATCAAC ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAG…        1
#> 2 AAACCGGAGG CTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATT…        1
#> 3 AAACTACCGG ACTGATACACTCCAAGCGTAGACAGACCAACTAGAAGTTGAGAAGTCTGCTCTGCAG…        1
#> 4 AAATCCAGAC ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAG…        1
#> 5 AAATGAACGA ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAG…        1
#> 6 AAATGATCCG ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAG…        1
```
