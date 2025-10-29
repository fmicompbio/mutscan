# Create a conversion table for collapsing similar sequences

Create a conversion table for collapsing similar sequences

## Usage

``` r
groupSimilarSequences(
  seqs,
  scores,
  collapseMaxDist = 0,
  collapseMinScore = 0,
  collapseMinRatio = 0,
  verbose = FALSE
)
```

## Arguments

- seqs:

  Character vector with nucleotide sequences (or pairs of sequences
  concatenated with "\_") to be collapsed. The sequences must all be of
  the same length.

- scores:

  Numeric vector of "scores" for the sequences. Typically the total
  read/UMI count. A higher score will be preferred when deciding which
  sequence to use as the representative for a group of collapsed
  sequences.

- collapseMaxDist:

  Numeric scalar defining the tolerance for collapsing similar
  sequences. If the value is in \[0, 1), it defines the maximal Hamming
  distance in terms of a fraction of sequence length:
  (`round(collapseMaxDist * nchar(sequence))`). A value greater or equal
  to 1 is rounded and directly used as the maximum allowed Hamming
  distance. Note that sequences can only be collapsed if they are all of
  the same length. The default value is 0.

- collapseMinScore:

  Numeric scalar, indicating the minimum score required for a sequence
  to be considered as a representative for a group of similar sequences
  (i.e., to allow other sequences to be collapsed into it). The default
  value is 0.

- collapseMinRatio:

  Numeric scalar. During collapsing of similar sequences, a
  low-frequency sequence will be collapsed with a higher-frequency
  sequence only if the ratio between the high-frequency and the
  low-frequency scores is at least this high. A value of 0 indicates
  that no such check is performed.

- verbose:

  Logical scalar, whether to print progress messages.

## Value

A data.frame with two columns, containing the input sequences and the
representatives for the groups resulting from grouping similar
sequences, respectively.

## Author

Michael Stadler, Charlotte Soneson

## Examples

``` r
seqs <- c("AACGTAGCA", "ACCGTAGCA", "AACGGAGCA", "ATCGGAGCA", "TGAGGCATA")
scores <- c(5, 1, 3, 1, 8)
groupSimilarSequences(seqs = seqs, scores = scores, 
                      collapseMaxDist = 1, collapseMinScore = 0, 
                      collapseMinRatio = 0, verbose = FALSE)
#>    sequence representative
#> 1 AACGTAGCA      AACGTAGCA
#> 2 ACCGTAGCA      AACGTAGCA
#> 3 AACGGAGCA      AACGTAGCA
#> 4 ATCGGAGCA      ATCGGAGCA
#> 5 TGAGGCATA      TGAGGCATA
                            
```
