# Read, filter and digest sequences from fastq file(s).

Read sequences for one or a pair of fastq files and digest them (extract
umis, constant and variable parts, filter, extract mismatch information
from constant and count the observed unique variable parts).
Alternatively, primer sequences could be specified, in which case the
sequence immediately following the primer will be considered the
variable sequence.

## Usage

``` r
digestFastqs(
  fastqForward,
  fastqReverse = NULL,
  mergeForwardReverse = FALSE,
  minOverlap = 0,
  maxOverlap = 0,
  minMergedLength = 0,
  maxMergedLength = 0,
  maxFracMismatchOverlap = 1,
  greedyOverlap = TRUE,
  revComplForward = FALSE,
  revComplReverse = FALSE,
  adapterForward = "",
  adapterReverse = "",
  elementsForward = "",
  elementLengthsForward = numeric(0),
  elementsReverse = "",
  elementLengthsReverse = numeric(0),
  primerForward = c(""),
  primerReverse = c(""),
  wildTypeForward = "",
  wildTypeReverse = "",
  constantForward = c(""),
  constantReverse = c(""),
  avePhredMinForward = 20,
  avePhredMinReverse = 20,
  variableNMaxForward = 0,
  variableNMaxReverse = 0,
  umiNMax = 0,
  nbrMutatedCodonsMaxForward = 1,
  nbrMutatedCodonsMaxReverse = 1,
  nbrMutatedBasesMaxForward = -1,
  nbrMutatedBasesMaxReverse = -1,
  forbiddenMutatedCodonsForward = "",
  forbiddenMutatedCodonsReverse = "",
  useTreeWTmatch = FALSE,
  collapseToWTForward = FALSE,
  collapseToWTReverse = FALSE,
  mutatedPhredMinForward = 0,
  mutatedPhredMinReverse = 0,
  mutNameDelimiter = ".",
  constantMaxDistForward = -1,
  constantMaxDistReverse = -1,
  umiCollapseMaxDist = 0,
  filteredReadsFastqForward = "",
  filteredReadsFastqReverse = "",
  maxNReads = -1,
  verbose = FALSE,
  nThreads = 1,
  chunkSize = 1e+05,
  maxReadLength = 1024
)
```

## Arguments

- fastqForward, fastqReverse:

  Character vectors, paths to gzipped FASTQ files corresponding to
  forward and reverse reads, respectively. If more than one
  forward/reverse sequence file is given, they need to be provided in
  the same order. Note that if multiple fastq files are provided, they
  are all assumed to correspond to the same sample, and will effectively
  be concatenated.

- mergeForwardReverse:

  Logical scalar, whether to fuse the forward and reverse variable
  sequences.

- minOverlap, maxOverlap:

  Numeric scalar, the minimal and maximal allowed overlap between the
  forward and reverse reads when merging. Only used if
  `mergeForwardReverse` is `TRUE`. If set to 0, only overlaps covering
  the full length of the shortest of the two reads will be considered.

- minMergedLength, maxMergedLength:

  Numeric scalar, the minimal and maximal allowed total length of the
  merged product (if `mergeForwardReverse` is `TRUE`). If set to 0, any
  length is allowed.

- maxFracMismatchOverlap:

  Numeric scalar, maximal mismatch rate in the overlap. Only used if
  `mergeForwardReverse` is `TRUE`.

- greedyOverlap:

  Logical scalar. If `TRUE`, the first overlap satisfying `minOverlap`,
  `maxOverlap`, `minMergedLength`, `maxMergedLength` and
  `maxFracMismatchOverlap` will be retained. If `FALSE`, all valid
  overlaps will be scored and the one with the highest score (largest
  number of matches) will be retained.

- revComplForward, revComplReverse:

  Logical scalar, whether to reverse complement the forward/reverse
  variable and constant sequences, respectively.

- adapterForward, adapterReverse:

  Character scalars, the adapter sequence for forward/reverse reads,
  respectively. If a forward/reverse read contains the corresponding
  adapter sequence, the sequence pair will be filtered out. If set to
  `NULL`, no adapter filtering is performed. The number of filtered read
  pairs are reported in the return value.

- elementsForward, elementsReverse:

  Character scalars representing the composition of the forward and
  reverse reads, respectively. The strings should consist only of the
  letters S (skip), C (constant), U (umi), P (primer), V (variable), and
  cover the full extent of the read. Most combinations are allowed (and
  a given letter can appear multiple times), but there can be at most
  one occurrence of P. If a given letter is included multiple times, the
  corresponding sequences will be concatenated in the output.

- elementLengthsForward, elementLengthsReverse:

  Numeric vectors containing the lengths of each read component from
  `elementsForward`/`elementsReverse`, respectively. If the length of
  one element is set to -1, it will be inferred from the other lengths
  (as the remainder of the read). At most one number (or one number on
  each side of the primer P) can be set to -1. The indicated length of
  the primer is not used (instead it's inferred from the provided primer
  sequence) and can also be set to -1.

- primerForward, primerReverse:

  Character vectors, representing the primer sequence(s) for
  forward/reverse reads, respectively. Only read pairs that contain
  perfect matches to both the forward and reverse primers (if given)
  will be retained. Multiple primers can be specified - they will be
  considered in order and the first match will be used.

- wildTypeForward, wildTypeReverse:

  Character scalars or named character vectors, the wild type sequence
  for the forward and reverse variable region. If given as a single
  string, the reference sequence will be named 'f' (for forward) or 'r'
  (for reverse).

- constantForward, constantReverse:

  Character vectors giving, the expected constant forward and reverse
  sequences. If more than one sequence is provided, they must all have
  the same length.

- avePhredMinForward, avePhredMinReverse:

  Numeric scalar, the minimum average Phred score in the variable region
  for a read to be retained. If a read pair contains both forward and
  reverse variable regions, the minimum average Phred score has to be
  achieved in both for a read pair to be retained.

- variableNMaxForward, variableNMaxReverse:

  Numeric scalar, the maximum number of Ns allowed in the variable
  region for a read to be retained.

- umiNMax:

  Numeric scalar, the maximum number of Ns allowed in the UMI for a read
  to be retained.

- nbrMutatedCodonsMaxForward, nbrMutatedCodonsMaxReverse:

  Numeric scalar, the maximum number of mutated codons that are allowed.
  Note that for the forward and reverse sequence, respectively, exactly
  one of `nbrMutatedCodonsMax` and `nbrMutatedBasesMax` must be -1, and
  the other must be a non-negative number. The one that is not -1 will
  be used to filter and name the identified mutants.

- nbrMutatedBasesMaxForward, nbrMutatedBasesMaxReverse:

  Numeric scalar, the maximum number of mutated bases that are allowed.
  Note that for the forward and reverse sequence, respectively, exactly
  one of `nbrMutatedCodonsMax` and `nbrMutatedBasesMax` must be -1, and
  the other must be a non-negative number. The one that is not -1 will
  be used to filter and name the identified mutants.

- forbiddenMutatedCodonsForward, forbiddenMutatedCodonsReverse:

  Character vector of codons (can contain ambiguous IUPAC characters,
  see
  [`IUPAC_CODE_MAP`](https://rdrr.io/pkg/Biostrings/man/IUPAC_CODE_MAP.html)).
  If a read pair contains a mutated codon matching this pattern, it will
  be filtered out.

- useTreeWTmatch:

  Logical scalar. Should a tree-based matching to wild type sequences be
  used if possible? If the number of allowed mismatches is small, and
  the number of wild type sequences is large, this is typically faster.

- collapseToWTForward, collapseToWTReverse:

  Logical scalar, indicating whether to just represent the observed
  variable sequence by the closest wildtype sequence rather than
  retaining the information about the mutations.

- mutatedPhredMinForward, mutatedPhredMinReverse:

  Numeric scalar, the minimum Phred score of a mutated base for the read
  to be retained. If any mutated base has a Phred score lower than
  `mutatedPhredMin`, the read (pair) will be discarded.

- mutNameDelimiter:

  Character scalar, the delimiter used in the naming of mutants.
  Generally, mutants will be named as XX{.}YY{.}NNN, where XX is the
  closest provided reference sequence, YY is the mutated base or codon
  number (depending on whether `nbrMutatedBases*` or `nbrMutatedCodons*`
  is specified), and NNN is the mutated base or codon. Here, {.} is the
  provided `mutNameDelimiter`. The delimiter must be a single character
  (not "\_"), and can not appear in any of the provided reference
  sequence names.

- constantMaxDistForward, constantMaxDistReverse:

  Numeric scalars, the maximum allowed Hamming distance between the
  extracted and expected constant sequence. If multiple constant
  sequences are provided, the most similar one is used. Reads with a
  larger distance to the expected constant sequence are discarded. If
  set to -1, no filtering is done.

- umiCollapseMaxDist:

  Numeric scalar defining the tolerances for collapsing similar UMI
  sequences. If the value is in \[0, 1), it defines the maximal Hamming
  distance in terms of a fraction of sequence length:
  (`round(umiCollapseMaxDist * nchar(umiSeq))`). A value greater or
  equal to 1 is rounded and directly used as the maximum allowed Hamming
  distance.

- filteredReadsFastqForward, filteredReadsFastqReverse:

  Character scalars, the names of a (pair of) FASTQ file(s) where
  filtered-out reads will be written. The name(s) should end in .gz (the
  output will always be compressed). If empty, filtered reads will not
  be written to a file.

- maxNReads:

  Integer scalar, the maximum number of reads to process. The first
  `maxNReads` read (pairs) in the FASTQ file(s) will be used. If set to
  -1, all reads in the FASTQ file(s) will be processed.

- verbose:

  Logical scalar, whether to print out progress messages.

- nThreads:

  Numeric scalar, the number of threads to use for parallel processing.

- chunkSize:

  Numeric scalar, the number of read (pairs) to keep in memory for
  parallel processing. Reduce from the default value if you run out of
  memory.

- maxReadLength:

  Numeric scalar, the maximum allowed read length. Longer read lengths
  lead to higher memory allocation, and may require the `chunkSize` to
  be decreased.

## Value

A list with four entries:

- summaryTable:

  A `data.frame` that contains, for each observed mutation combination,
  the corresponding variable region sequences (or pair of sequences),
  the number of observed such sequences, and the number of unique UMIs
  observed for the sequence. It also has additional columns:
  'maxNbrReads' contains the number of reads for the most frequent
  observed sequence represented by the feature (only relevant if similar
  variable regions are collapsed). 'nbrMutBases', 'nbrMutCodons' and
  'nbrMutAAs' give the number of mutated bases, codons or amino acids in
  each variant. Alternative variant names based on base, codon or amino
  acid sequence are provided in columns mutantNameBase',
  'mutantNameCodon', 'mutantNameAA'. In addition, mutantNameBaseHGVS'
  and 'mutantNameAAHGVS' give base- and amino acid-based names following
  the HGVS nomenclature (https://varnomen.hgvs.org/). Please note that
  the provided reference sequence names are used for the HGVS sequence
  identifiers. It is up to the user to use appropriately named reference
  sequences in order to obtain valid HGVS variant names.

- filterSummary:

  A `data.frame` that contains the number of input reads, the number of
  reads filtered out in the processing, and the number of retained
  reads. The filters are named according to the convention "fxx_filter",
  where "xx" indicates the order in which the filters were applied, and
  "filter" indicates the type of filter. Note that filters are applied
  successively, and the reads filtered out in one step are not
  considered for successive filtering steps.

- errorStatistics:

  A `data.frame` that contains, for each Phred quality score between 0
  and 99, the number of bases in the extracted constant sequences with
  that quality score that match/mismatch with the provided reference
  constant sequence.

- parameters:

  A `list` with all parameter settings that were used in the processing.
  Also contains the version of the package and the time of processing.

## Details

The processing of a read pair goes as follows:

1.  Search for perfect matches to forward/reverse adapter sequences,
    filter out the read pair if a match is found in either the forward
    or reverse read.

2.  If primer sequences are provided, search for perfect matches, and
    filter out the read pair if not all provided primer sequences can be
    found.

3.  Extract the UMI, constant and variable sequence from forward and
    reverse reads, based on the definition of the respective read
    composition.

4.  If requested, collapse forward and reverse variable regions by
    retaining, for each position, the base with the highest reported
    base quality.

5.  Filter out the read (pair) if the average quality in the variable
    region is below `avePhredMinForward`/`avePhredMinReverse`, in either
    the forward or reverse read (or the merged read).

6.  Filter out the read (pair) if the number of Ns in the variable
    region exceeds `variableNMaxForward`/`variableNMaxReverse`.

7.  Filter out the read (pair) if the number of Ns in the combined
    forward and reverse UMI sequence exceeds `umiNMax`

8.  If one or more wild type sequences (for the variable region) are
    provided, find the mismatches between the (forward/reverse) variable
    region and the provided wild type sequence (if more than one wild
    type sequence is provided, first find the one that is closest to the
    read).

9.  Filter out the read (pair) if any mutated base has a quality below
    `mutatedPhredMinForward`/`mutatedPhredMinReverse`.

10. Filter out the read (pair) if the number of mutated codons exceeds
    `nbrMutatedCodonsMaxForward`/`nbrMutatedCodonsMaxReverse`.

11. Filter out the read (pair) if any of the mutated codons match any of
    the codons encoded by
    `forbiddenMutatedCodonsForward`/`forbiddenMutatedCodonsReverse`.

12. Assign a 'mutation name' to the read (pair). This name is a
    combination of parts of the form XX{.}YY{.}NNN, where XX is the name
    of the most similar reference sequence, YY is the mutated codon
    number, and NNN is the mutated codon. {.} is a delimiter, specified
    via `mutNameDelimiter`. If no wildtype sequences are provided, the
    variable sequence will be used as the mutation name'.

Based on the retained reads following this filtering process, count the
number of reads, and the number of unique UMIs, for each variable
sequence (or pair of variable sequences).

## Examples

``` r
## See the vignette for complete worked-out examples for different types of 
## data sets

## ---------------------------------------------------------------------- ## 
## Process a single-end data set, assume that the full read represents  
## the variable region
out <- digestFastqs(
    fastqForward = system.file("extdata", "cisInput_1.fastq.gz", 
                               package = "mutscan"), 
    elementsForward = "V", elementLengthsForward = -1
)
## Table with read counts and mutant information
head(out$summaryTable)
#>                                                                                                                      mutantName
#> 1 AAAACTACCGGAACCGGAGGAGGGAGCTGACTGATACACTCCAAGCGTAGACAGACCAACTAGAAGTTGAGAAGTCTGCTCTGCAGACCGAGATTGCAAACCTGCAGAAGGAGAAGGAAAAACTA
#> 2 AAACACATGTCAACCGGAGGAGGGAGCTGAATGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 3 AAACCCGCGACAACCGGAGGAGGGAGCTGACTGATACACTGCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAAGCTGCTAAAGGAGAAGGAAAAACTA
#> 4 AAACCTTGGTGAACCGGAGGAGGGAGCTGACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 5 AAACTAGCACCAACCGGAGGAGGGAGCTGACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 6 AAATTCGCGTCAACCGGAGGAGGGAGCTGACTGATACACTCCAAGAGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#>                                                                                                                        sequence
#> 1 AAAACTACCGGAACCGGAGGAGGGAGCTGACTGATACACTCCAAGCGTAGACAGACCAACTAGAAGTTGAGAAGTCTGCTCTGCAGACCGAGATTGCAAACCTGCAGAAGGAGAAGGAAAAACTA
#> 2 AAACACATGTCAACCGGAGGAGGGAGCTGAATGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 3 AAACCCGCGACAACCGGAGGAGGGAGCTGACTGATACACTGCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAAGCTGCTAAAGGAGAAGGAAAAACTA
#> 4 AAACCTTGGTGAACCGGAGGAGGGAGCTGACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 5 AAACTAGCACCAACCGGAGGAGGGAGCTGACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 6 AAATTCGCGTCAACCGGAGGAGGGAGCTGACTGATACACTCCAAGAGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#>   nbrReads maxNbrReads nbrUmis nbrMutBases nbrMutCodons nbrMutAAs varLengths
#> 1        1           1       0           0            0         0        125
#> 2        1           1       0           0            0         0        125
#> 3        1           1       0           0            0         0        125
#> 4        1           1       0           0            0         0        125
#> 5        1           1       0           0            0         0        125
#> 6        1           1       0           0            0         0        125
#>                                                                                                                  mutantNameBase
#> 1 AAAACTACCGGAACCGGAGGAGGGAGCTGACTGATACACTCCAAGCGTAGACAGACCAACTAGAAGTTGAGAAGTCTGCTCTGCAGACCGAGATTGCAAACCTGCAGAAGGAGAAGGAAAAACTA
#> 2 AAACACATGTCAACCGGAGGAGGGAGCTGAATGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 3 AAACCCGCGACAACCGGAGGAGGGAGCTGACTGATACACTGCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAAGCTGCTAAAGGAGAAGGAAAAACTA
#> 4 AAACCTTGGTGAACCGGAGGAGGGAGCTGACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 5 AAACTAGCACCAACCGGAGGAGGGAGCTGACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 6 AAATTCGCGTCAACCGGAGGAGGGAGCTGACTGATACACTCCAAGAGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#>                                                                                                                 mutantNameCodon
#> 1 AAAACTACCGGAACCGGAGGAGGGAGCTGACTGATACACTCCAAGCGTAGACAGACCAACTAGAAGTTGAGAAGTCTGCTCTGCAGACCGAGATTGCAAACCTGCAGAAGGAGAAGGAAAAACTA
#> 2 AAACACATGTCAACCGGAGGAGGGAGCTGAATGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 3 AAACCCGCGACAACCGGAGGAGGGAGCTGACTGATACACTGCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAAGCTGCTAAAGGAGAAGGAAAAACTA
#> 4 AAACCTTGGTGAACCGGAGGAGGGAGCTGACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 5 AAACTAGCACCAACCGGAGGAGGGAGCTGACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 6 AAATTCGCGTCAACCGGAGGAGGGAGCTGACTGATACACTCCAAGAGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#>   mutantNameBaseHGVS                              mutantNameAA mutantNameAAHGVS
#> 1                    KTTGTGGGS*LIHSKRRQTN*KLRSLLCRPRLQTCRRRRKN                 
#> 2                    KHMSTGGGS*MIHSKRRQTN*KMRSLLCRPRLPTC*RRRKN                 
#> 3                    KPATTGGGS*LIHCKRRQTN*KMRSLLCRPRLPSC*RRRKN                 
#> 4                    KPW*TGGGS*LIHSKRRQTN*KMRSLLCRPRLPTC*RRRKN                 
#> 5                    KLAPTGGGS*LIHSKRRQTN*KMRSLLCRPRLPTC*RRRKN                 
#> 6                    KFASTGGGS*LIHSKRRQTN*KMRSLLCRPRLPTC*RRRKN                 
#>   mutationTypes                                sequenceAA
#> 1               KTTGTGGGS*LIHSKRRQTN*KLRSLLCRPRLQTCRRRRKN
#> 2               KHMSTGGGS*MIHSKRRQTN*KMRSLLCRPRLPTC*RRRKN
#> 3               KPATTGGGS*LIHCKRRQTN*KMRSLLCRPRLPSC*RRRKN
#> 4               KPW*TGGGS*LIHSKRRQTN*KMRSLLCRPRLPTC*RRRKN
#> 5               KLAPTGGGS*LIHSKRRQTN*KMRSLLCRPRLPTC*RRRKN
#> 6               KFASTGGGS*LIHSKRRQTN*KMRSLLCRPRLPTC*RRRKN
## Filter summary
out$filterSummary
#>   nbrTotal f1_nbrAdapter f2_nbrNoPrimer f3_nbrReadWrongLength
#> 1     1000             0              0                     0
#>   f4_nbrNoValidOverlap f5_nbrAvgVarQualTooLow f6_nbrTooManyNinVar
#> 1                    0                      1                 593
#>   f7_nbrTooManyNinUMI f8_nbrTooManyBestWTHits f9_nbrMutQualTooLow
#> 1                   0                       0                   0
#>   f10a_nbrTooManyMutCodons f10b_nbrTooManyMutBases f11_nbrForbiddenCodons
#> 1                        0                       0                      0
#>   f12_nbrTooManyMutConstant f13_nbrTooManyBestConstantHits nbrRetained
#> 1                         0                              0         406

## ---------------------------------------------------------------------- ## 
## Process a single-end data set, specify the read as a combination of 
## UMI, constant region and variable region (skip the first base)
out <- digestFastqs(
    fastqForward = system.file("extdata", "cisInput_1.fastq.gz", 
                               package = "mutscan"), 
    elementsForward = "SUCV", elementLengthsForward = c(1, 10, 18, 96), 
    constantForward = "AACCGGAGGAGGGAGCTG"
)
## Table with read counts and mutant information
head(out$summaryTable)
#>                                                                                         mutantName
#> 1 AAATATAACGTTGACGATGTAGCTTTAGGTGTCTGTAAAACAGGTGCCGAAGAAGCTGGAGTAACAGAAGTGAGAACCAGCTTATCAGAAAAAAAG
#> 2 AAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAATAAACAAAGTACTTTCTATTTTCTAT
#> 3 AAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAATAATACTATTTTTTTTGTTTATTTCCA
#> 4 AATGATACACTCCAAACGGAGACAGACCAACTAGAATATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCGGAAAGAGGAGGAAAAACTA
#> 5 AATGATACACTCCAAGCGGAGACAGAACAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 6 AATGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#>                                                                                           sequence
#> 1 AAATATAACGTTGACGATGTAGCTTTAGGTGTCTGTAAAACAGGTGCCGAAGAAGCTGGAGTAACAGAAGTGAGAACCAGCTTATCAGAAAAAAAG
#> 2 AAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAATAAACAAAGTACTTTCTATTTTCTAT
#> 3 AAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAATAATACTATTTTTTTTGTTTATTTCCA
#> 4 AATGATACACTCCAAACGGAGACAGACCAACTAGAATATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCGGAAAGAGGAGGAAAAACTA
#> 5 AATGATACACTCCAAGCGGAGACAGAACAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 6 AATGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#>   nbrReads maxNbrReads nbrUmis nbrMutBases nbrMutCodons nbrMutAAs varLengths
#> 1        1           1       1           0            0         0         96
#> 2        1           1       1           0            0         0         96
#> 3        1           1       1           0            0         0         96
#> 4        1           1       1           0            0         0         96
#> 5        1           1       1           0            0         0         96
#> 6        2           2       2           0            0         0         96
#>                                                                                     mutantNameBase
#> 1 AAATATAACGTTGACGATGTAGCTTTAGGTGTCTGTAAAACAGGTGCCGAAGAAGCTGGAGTAACAGAAGTGAGAACCAGCTTATCAGAAAAAAAG
#> 2 AAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAATAAACAAAGTACTTTCTATTTTCTAT
#> 3 AAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAATAATACTATTTTTTTTGTTTATTTCCA
#> 4 AATGATACACTCCAAACGGAGACAGACCAACTAGAATATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCGGAAAGAGGAGGAAAAACTA
#> 5 AATGATACACTCCAAGCGGAGACAGAACAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 6 AATGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#>                                                                                    mutantNameCodon
#> 1 AAATATAACGTTGACGATGTAGCTTTAGGTGTCTGTAAAACAGGTGCCGAAGAAGCTGGAGTAACAGAAGTGAGAACCAGCTTATCAGAAAAAAAG
#> 2 AAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAATAAACAAAGTACTTTCTATTTTCTAT
#> 3 AAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAATAATACTATTTTTTTTGTTTATTTCCA
#> 4 AATGATACACTCCAAACGGAGACAGACCAACTAGAATATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCGGAAAGAGGAGGAAAAACTA
#> 5 AATGATACACTCCAAGCGGAGACAGAACAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 6 AATGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#>   mutantNameBaseHGVS                     mutantNameAA mutantNameAAHGVS
#> 1                    KYNVDDVALGVCKTGAEEAGVTEVRTSLSEKK                 
#> 2                    KSTRLNSSHITISYAVFCLKKKKINKVLSIFY                 
#> 3                    KSTRLNSSHITISYAVFCLKKKK*YYFFCLFP                 
#> 4                    NDTLQTETDQLEYEKSALQTEIANLRKEEEKL                 
#> 5                    NDTLQAETEQLEDEKSALQTEIANLLKEKEKL                 
#> 6                    NDTLQAETDQLEDEKSALQTEIANLLKEKEKL                 
#>   mutationTypes                       sequenceAA
#> 1               KYNVDDVALGVCKTGAEEAGVTEVRTSLSEKK
#> 2               KSTRLNSSHITISYAVFCLKKKKINKVLSIFY
#> 3               KSTRLNSSHITISYAVFCLKKKK*YYFFCLFP
#> 4               NDTLQTETDQLEYEKSALQTEIANLRKEEEKL
#> 5               NDTLQAETEQLEDEKSALQTEIANLLKEKEKL
#> 6               NDTLQAETDQLEDEKSALQTEIANLLKEKEKL
## Filter summary
out$filterSummary
#>   nbrTotal f1_nbrAdapter f2_nbrNoPrimer f3_nbrReadWrongLength
#> 1     1000             0              0                     0
#>   f4_nbrNoValidOverlap f5_nbrAvgVarQualTooLow f6_nbrTooManyNinVar
#> 1                    0                      3                 591
#>   f7_nbrTooManyNinUMI f8_nbrTooManyBestWTHits f9_nbrMutQualTooLow
#> 1                   0                       0                   0
#>   f10a_nbrTooManyMutCodons f10b_nbrTooManyMutBases f11_nbrForbiddenCodons
#> 1                        0                       0                      0
#>   f12_nbrTooManyMutConstant f13_nbrTooManyBestConstantHits nbrRetained
#> 1                         0                              0         406
## Error statistics
out$errorStatistics
#>     PhredQuality nbrMatchForward nbrMismatchForward nbrMatchReverse
#> 1              0               0                  0               0
#> 2              1               0                  0               0
#> 3              2               0                  0               0
#> 4              3               0                  0               0
#> 5              4               0                  0               0
#> 6              5               0                  0               0
#> 7              6               0                  0               0
#> 8              7               0                  0               0
#> 9              8               0                  0               0
#> 10             9               0                  0               0
#> 11            10               0                  0               0
#> 12            11               0                  0               0
#> 13            12               0                  0               0
#> 14            13               0                  0               0
#> 15            14             121                 66               0
#> 16            15               0                  0               0
#> 17            16               0                  0               0
#> 18            17               0                  0               0
#> 19            18               0                  0               0
#> 20            19               0                  0               0
#> 21            20               0                  0               0
#> 22            21               0                  0               0
#> 23            22              30                  1               0
#> 24            23               0                  0               0
#> 25            24               0                  0               0
#> 26            25               0                  0               0
#> 27            26               0                  0               0
#> 28            27             208                 56               0
#> 29            28               0                  0               0
#> 30            29               0                  0               0
#> 31            30               0                  0               0
#> 32            31               0                  0               0
#> 33            32               0                  0               0
#> 34            33             393                 96               0
#> 35            34               0                  0               0
#> 36            35               0                  0               0
#> 37            36               0                  0               0
#> 38            37            5929                408               0
#> 39            38               0                  0               0
#> 40            39               0                  0               0
#> 41            40               0                  0               0
#> 42            41               0                  0               0
#> 43            42               0                  0               0
#> 44            43               0                  0               0
#> 45            44               0                  0               0
#> 46            45               0                  0               0
#> 47            46               0                  0               0
#> 48            47               0                  0               0
#> 49            48               0                  0               0
#> 50            49               0                  0               0
#> 51            50               0                  0               0
#> 52            51               0                  0               0
#> 53            52               0                  0               0
#> 54            53               0                  0               0
#> 55            54               0                  0               0
#> 56            55               0                  0               0
#> 57            56               0                  0               0
#> 58            57               0                  0               0
#> 59            58               0                  0               0
#> 60            59               0                  0               0
#> 61            60               0                  0               0
#> 62            61               0                  0               0
#> 63            62               0                  0               0
#> 64            63               0                  0               0
#> 65            64               0                  0               0
#> 66            65               0                  0               0
#> 67            66               0                  0               0
#> 68            67               0                  0               0
#> 69            68               0                  0               0
#> 70            69               0                  0               0
#> 71            70               0                  0               0
#> 72            71               0                  0               0
#> 73            72               0                  0               0
#> 74            73               0                  0               0
#> 75            74               0                  0               0
#> 76            75               0                  0               0
#> 77            76               0                  0               0
#> 78            77               0                  0               0
#> 79            78               0                  0               0
#> 80            79               0                  0               0
#> 81            80               0                  0               0
#> 82            81               0                  0               0
#> 83            82               0                  0               0
#> 84            83               0                  0               0
#> 85            84               0                  0               0
#> 86            85               0                  0               0
#> 87            86               0                  0               0
#> 88            87               0                  0               0
#> 89            88               0                  0               0
#> 90            89               0                  0               0
#> 91            90               0                  0               0
#> 92            91               0                  0               0
#> 93            92               0                  0               0
#> 94            93               0                  0               0
#> 95            94               0                  0               0
#> 96            95               0                  0               0
#> 97            96               0                  0               0
#> 98            97               0                  0               0
#> 99            98               0                  0               0
#> 100           99               0                  0               0
#>     nbrMismatchReverse
#> 1                    0
#> 2                    0
#> 3                    0
#> 4                    0
#> 5                    0
#> 6                    0
#> 7                    0
#> 8                    0
#> 9                    0
#> 10                   0
#> 11                   0
#> 12                   0
#> 13                   0
#> 14                   0
#> 15                   0
#> 16                   0
#> 17                   0
#> 18                   0
#> 19                   0
#> 20                   0
#> 21                   0
#> 22                   0
#> 23                   0
#> 24                   0
#> 25                   0
#> 26                   0
#> 27                   0
#> 28                   0
#> 29                   0
#> 30                   0
#> 31                   0
#> 32                   0
#> 33                   0
#> 34                   0
#> 35                   0
#> 36                   0
#> 37                   0
#> 38                   0
#> 39                   0
#> 40                   0
#> 41                   0
#> 42                   0
#> 43                   0
#> 44                   0
#> 45                   0
#> 46                   0
#> 47                   0
#> 48                   0
#> 49                   0
#> 50                   0
#> 51                   0
#> 52                   0
#> 53                   0
#> 54                   0
#> 55                   0
#> 56                   0
#> 57                   0
#> 58                   0
#> 59                   0
#> 60                   0
#> 61                   0
#> 62                   0
#> 63                   0
#> 64                   0
#> 65                   0
#> 66                   0
#> 67                   0
#> 68                   0
#> 69                   0
#> 70                   0
#> 71                   0
#> 72                   0
#> 73                   0
#> 74                   0
#> 75                   0
#> 76                   0
#> 77                   0
#> 78                   0
#> 79                   0
#> 80                   0
#> 81                   0
#> 82                   0
#> 83                   0
#> 84                   0
#> 85                   0
#> 86                   0
#> 87                   0
#> 88                   0
#> 89                   0
#> 90                   0
#> 91                   0
#> 92                   0
#> 93                   0
#> 94                   0
#> 95                   0
#> 96                   0
#> 97                   0
#> 98                   0
#> 99                   0
#> 100                  0

## ---------------------------------------------------------------------- ## 
## Process a single-end data set, specify the read as a combination of 
## UMI, constant region and variable region (skip the first base), provide 
## the wild type sequence to compare the variable region to and limit the 
## number of allowed mutated codons to 1
out <- digestFastqs(
    fastqForward = system.file("extdata", "cisInput_1.fastq.gz", 
                               package = "mutscan"), 
    elementsForward = "SUCV", elementLengthsForward = c(1, 10, 18, 96), 
    constantForward = "AACCGGAGGAGGGAGCTG", 
    wildTypeForward = c(FOS = paste0(
        "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTC", 
        "TGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA")),
    nbrMutatedCodonsMaxForward = 1
)
## Table with read counts and mutant information
head(out$summaryTable)
#>   mutantName
#> 1   FOS.0.WT
#> 2  FOS.1.AAT
#> 3  FOS.1.ACC
#> 4  FOS.1.ACG
#> 5  FOS.1.CCT
#> 6 FOS.10.CAT
#>                                                                                           sequence
#> 1 ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 2 AATGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 3 ACCGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 4 ACGGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 5 CCTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 6 ACTGATACACTCCAAGCGGAGACAGACCATCTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#>   nbrReads maxNbrReads nbrUmis nbrMutBases nbrMutCodons nbrMutAAs varLengths
#> 1       30          30      30           0            0         0         96
#> 2        2           2       2           1            1         1         96
#> 3        1           1       1           1            1         0         96
#> 4        1           1       1           1            1         0         96
#> 5        1           1       1           1            1         1         96
#> 6        1           1       1           1            1         1         96
#>   mutantNameBase mutantNameCodon mutantNameBaseHGVS mutantNameAA
#> 1       FOS.0.WT        FOS.0.WT              FOS:c     FOS.0.WT
#> 2        FOS.2.A       FOS.1.AAT         FOS:c.2C>A      FOS.1.N
#> 3        FOS.3.C       FOS.1.ACC         FOS:c.3T>C     FOS.0.WT
#> 4        FOS.3.G       FOS.1.ACG         FOS:c.3T>G     FOS.0.WT
#> 5        FOS.1.C       FOS.1.CCT         FOS:c.1A>C      FOS.1.P
#> 6       FOS.30.T      FOS.10.CAT        FOS:c.30A>T     FOS.10.H
#>   mutantNameAAHGVS mutationTypes                       sequenceAA
#> 1            FOS:p               TDTLQAETDQLEDEKSALQTEIANLLKEKEKL
#> 2  FOS:p.(Thr1Asn) nonsynonymous NDTLQAETDQLEDEKSALQTEIANLLKEKEKL
#> 3            FOS:p        silent TDTLQAETDQLEDEKSALQTEIANLLKEKEKL
#> 4            FOS:p        silent TDTLQAETDQLEDEKSALQTEIANLLKEKEKL
#> 5  FOS:p.(Thr1Pro) nonsynonymous PDTLQAETDQLEDEKSALQTEIANLLKEKEKL
#> 6 FOS:p.(Gln10His) nonsynonymous TDTLQAETDHLEDEKSALQTEIANLLKEKEKL
## Filter summary
out$filterSummary
#>   nbrTotal f1_nbrAdapter f2_nbrNoPrimer f3_nbrReadWrongLength
#> 1     1000             0              0                     0
#>   f4_nbrNoValidOverlap f5_nbrAvgVarQualTooLow f6_nbrTooManyNinVar
#> 1                    0                      3                 591
#>   f7_nbrTooManyNinUMI f8_nbrTooManyBestWTHits f9_nbrMutQualTooLow
#> 1                   0                       0                   0
#>   f10a_nbrTooManyMutCodons f10b_nbrTooManyMutBases f11_nbrForbiddenCodons
#> 1                      297                       0                      0
#>   f12_nbrTooManyMutConstant f13_nbrTooManyBestConstantHits nbrRetained
#> 1                         0                              0         109
## Error statistics
out$errorStatistics
#>     PhredQuality nbrMatchForward nbrMismatchForward nbrMatchReverse
#> 1              0               0                  0               0
#> 2              1               0                  0               0
#> 3              2               0                  0               0
#> 4              3               0                  0               0
#> 5              4               0                  0               0
#> 6              5               0                  0               0
#> 7              6               0                  0               0
#> 8              7               0                  0               0
#> 9              8               0                  0               0
#> 10             9               0                  0               0
#> 11            10               0                  0               0
#> 12            11               0                  0               0
#> 13            12               0                  0               0
#> 14            13               0                  0               0
#> 15            14              17                  0               0
#> 16            15               0                  0               0
#> 17            16               0                  0               0
#> 18            17               0                  0               0
#> 19            18               0                  0               0
#> 20            19               0                  0               0
#> 21            20               0                  0               0
#> 22            21               0                  0               0
#> 23            22               5                  0               0
#> 24            23               0                  0               0
#> 25            24               0                  0               0
#> 26            25               0                  0               0
#> 27            26               0                  0               0
#> 28            27              31                  0               0
#> 29            28               0                  0               0
#> 30            29               0                  0               0
#> 31            30               0                  0               0
#> 32            31               0                  0               0
#> 33            32               0                  0               0
#> 34            33              88                  1               0
#> 35            34               0                  0               0
#> 36            35               0                  0               0
#> 37            36               0                  0               0
#> 38            37            1819                  1               0
#> 39            38               0                  0               0
#> 40            39               0                  0               0
#> 41            40               0                  0               0
#> 42            41               0                  0               0
#> 43            42               0                  0               0
#> 44            43               0                  0               0
#> 45            44               0                  0               0
#> 46            45               0                  0               0
#> 47            46               0                  0               0
#> 48            47               0                  0               0
#> 49            48               0                  0               0
#> 50            49               0                  0               0
#> 51            50               0                  0               0
#> 52            51               0                  0               0
#> 53            52               0                  0               0
#> 54            53               0                  0               0
#> 55            54               0                  0               0
#> 56            55               0                  0               0
#> 57            56               0                  0               0
#> 58            57               0                  0               0
#> 59            58               0                  0               0
#> 60            59               0                  0               0
#> 61            60               0                  0               0
#> 62            61               0                  0               0
#> 63            62               0                  0               0
#> 64            63               0                  0               0
#> 65            64               0                  0               0
#> 66            65               0                  0               0
#> 67            66               0                  0               0
#> 68            67               0                  0               0
#> 69            68               0                  0               0
#> 70            69               0                  0               0
#> 71            70               0                  0               0
#> 72            71               0                  0               0
#> 73            72               0                  0               0
#> 74            73               0                  0               0
#> 75            74               0                  0               0
#> 76            75               0                  0               0
#> 77            76               0                  0               0
#> 78            77               0                  0               0
#> 79            78               0                  0               0
#> 80            79               0                  0               0
#> 81            80               0                  0               0
#> 82            81               0                  0               0
#> 83            82               0                  0               0
#> 84            83               0                  0               0
#> 85            84               0                  0               0
#> 86            85               0                  0               0
#> 87            86               0                  0               0
#> 88            87               0                  0               0
#> 89            88               0                  0               0
#> 90            89               0                  0               0
#> 91            90               0                  0               0
#> 92            91               0                  0               0
#> 93            92               0                  0               0
#> 94            93               0                  0               0
#> 95            94               0                  0               0
#> 96            95               0                  0               0
#> 97            96               0                  0               0
#> 98            97               0                  0               0
#> 99            98               0                  0               0
#> 100           99               0                  0               0
#>     nbrMismatchReverse
#> 1                    0
#> 2                    0
#> 3                    0
#> 4                    0
#> 5                    0
#> 6                    0
#> 7                    0
#> 8                    0
#> 9                    0
#> 10                   0
#> 11                   0
#> 12                   0
#> 13                   0
#> 14                   0
#> 15                   0
#> 16                   0
#> 17                   0
#> 18                   0
#> 19                   0
#> 20                   0
#> 21                   0
#> 22                   0
#> 23                   0
#> 24                   0
#> 25                   0
#> 26                   0
#> 27                   0
#> 28                   0
#> 29                   0
#> 30                   0
#> 31                   0
#> 32                   0
#> 33                   0
#> 34                   0
#> 35                   0
#> 36                   0
#> 37                   0
#> 38                   0
#> 39                   0
#> 40                   0
#> 41                   0
#> 42                   0
#> 43                   0
#> 44                   0
#> 45                   0
#> 46                   0
#> 47                   0
#> 48                   0
#> 49                   0
#> 50                   0
#> 51                   0
#> 52                   0
#> 53                   0
#> 54                   0
#> 55                   0
#> 56                   0
#> 57                   0
#> 58                   0
#> 59                   0
#> 60                   0
#> 61                   0
#> 62                   0
#> 63                   0
#> 64                   0
#> 65                   0
#> 66                   0
#> 67                   0
#> 68                   0
#> 69                   0
#> 70                   0
#> 71                   0
#> 72                   0
#> 73                   0
#> 74                   0
#> 75                   0
#> 76                   0
#> 77                   0
#> 78                   0
#> 79                   0
#> 80                   0
#> 81                   0
#> 82                   0
#> 83                   0
#> 84                   0
#> 85                   0
#> 86                   0
#> 87                   0
#> 88                   0
#> 89                   0
#> 90                   0
#> 91                   0
#> 92                   0
#> 93                   0
#> 94                   0
#> 95                   0
#> 96                   0
#> 97                   0
#> 98                   0
#> 99                   0
#> 100                  0

## ---------------------------------------------------------------------- ## 
## Process a paired-end data set where both the forward and reverse reads 
## contain the same variable region and thus should be merged to generate 
## the final variable sequence, specify the reads as a combination of 
## UMI, constant region and variable region (skip the first and/or last
## base), provide the wild type sequence to compare the variable region to 
## and limit the number of allowed mutated codons to 1
out <- digestFastqs(
    fastqForward = system.file("extdata", "cisInput_1.fastq.gz", 
                               package = "mutscan"),
    fastqReverse = system.file("extdata", "cisInput_2.fastq.gz",
                               package = "mutscan"), 
    mergeForwardReverse = TRUE, 
    revComplForward = FALSE, revComplReverse = TRUE, 
    elementsForward = "SUCV", elementLengthsForward = c(1, 10, 18, 96),
    elementsReverse = "SUCVS", elementLengthsReverse = c(1, 7, 17, 96, -1),
    constantForward = "AACCGGAGGAGGGAGCTG", 
    constantReverse = "GAGTTCATCCTGGCAGC",
    wildTypeForward = c(FOS = paste0(
        "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTC", 
        "TGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA")),
    nbrMutatedCodonsMaxForward = 1
)
## Table with read counts and mutant information
head(out$summaryTable)
#>   mutantName
#> 1   FOS.0.WT
#> 2  FOS.1.AAT
#> 3  FOS.1.ACC
#> 4  FOS.1.ACG
#> 5  FOS.1.CCT
#> 6  FOS.1.CGT
#>                                                                                           sequence
#> 1 ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 2 AATGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 3 ACCGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 4 ACGGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 5 CCTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#> 6 CGTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA
#>   nbrReads maxNbrReads nbrUmis nbrMutBases nbrMutCodons nbrMutAAs varLengths
#> 1       77          77      77           0            0         0         96
#> 2        2           2       2           1            1         1         96
#> 3        2           2       2           1            1         0         96
#> 4        1           1       1           1            1         0         96
#> 5        1           1       1           1            1         1         96
#> 6        1           1       1           2            1         1         96
#>    mutantNameBase mutantNameCodon mutantNameBaseHGVS mutantNameAA
#> 1        FOS.0.WT        FOS.0.WT              FOS:c     FOS.0.WT
#> 2         FOS.2.A       FOS.1.AAT         FOS:c.2C>A      FOS.1.N
#> 3         FOS.3.C       FOS.1.ACC         FOS:c.3T>C     FOS.0.WT
#> 4         FOS.3.G       FOS.1.ACG         FOS:c.3T>G     FOS.0.WT
#> 5         FOS.1.C       FOS.1.CCT         FOS:c.1A>C      FOS.1.P
#> 6 FOS.1.C_FOS.2.G       FOS.1.CGT  FOS:c.1_2delinsCG      FOS.1.R
#>   mutantNameAAHGVS mutationTypes                       sequenceAA
#> 1            FOS:p               TDTLQAETDQLEDEKSALQTEIANLLKEKEKL
#> 2  FOS:p.(Thr1Asn) nonsynonymous NDTLQAETDQLEDEKSALQTEIANLLKEKEKL
#> 3            FOS:p        silent TDTLQAETDQLEDEKSALQTEIANLLKEKEKL
#> 4            FOS:p        silent TDTLQAETDQLEDEKSALQTEIANLLKEKEKL
#> 5  FOS:p.(Thr1Pro) nonsynonymous PDTLQAETDQLEDEKSALQTEIANLLKEKEKL
#> 6  FOS:p.(Thr1Arg) nonsynonymous RDTLQAETDQLEDEKSALQTEIANLLKEKEKL
## Filter summary
out$filterSummary
#>   nbrTotal f1_nbrAdapter f2_nbrNoPrimer f3_nbrReadWrongLength
#> 1     1000             0              0                     0
#>   f4_nbrNoValidOverlap f5_nbrAvgVarQualTooLow f6_nbrTooManyNinVar
#> 1                    0                      0                  52
#>   f7_nbrTooManyNinUMI f8_nbrTooManyBestWTHits f9_nbrMutQualTooLow
#> 1                   0                       0                   0
#>   f10a_nbrTooManyMutCodons f10b_nbrTooManyMutBases f11_nbrForbiddenCodons
#> 1                      699                       0                      0
#>   f12_nbrTooManyMutConstant f13_nbrTooManyBestConstantHits nbrRetained
#> 1                         0                              0         249
## Error statistics
out$errorStatistics
#>     PhredQuality nbrMatchForward nbrMismatchForward nbrMatchReverse
#> 1              0               0                  0               0
#> 2              1               0                  0               0
#> 3              2               0                  0               0
#> 4              3               0                  0               0
#> 5              4               0                  0               0
#> 6              5               0                  0               0
#> 7              6               0                  0               0
#> 8              7               0                  0               0
#> 9              8               0                  0               0
#> 10             9               0                  0               0
#> 11            10               0                  0               0
#> 12            11               0                  0               0
#> 13            12               0                  0               0
#> 14            13               0                  0               0
#> 15            14              69                  2              53
#> 16            15               0                  0               0
#> 17            16               0                  0               0
#> 18            17               0                  0               0
#> 19            18               0                  0               0
#> 20            19               0                  0               0
#> 21            20               0                  0               0
#> 22            21               0                  0               0
#> 23            22              23                  0               0
#> 24            23               0                  0               0
#> 25            24               0                  0               0
#> 26            25               0                  0               0
#> 27            26               0                  0               0
#> 28            27             115                  1              74
#> 29            28               0                  0               0
#> 30            29               0                  0               0
#> 31            30               0                  0               0
#> 32            31               0                  0               0
#> 33            32               0                  0               0
#> 34            33             252                  1             143
#> 35            34               0                  0               0
#> 36            35               0                  0               0
#> 37            36               0                  0               0
#> 38            37            4017                  2            3930
#> 39            38               0                  0               0
#> 40            39               0                  0               0
#> 41            40               0                  0               0
#> 42            41               0                  0               0
#> 43            42               0                  0               0
#> 44            43               0                  0               0
#> 45            44               0                  0               0
#> 46            45               0                  0               0
#> 47            46               0                  0               0
#> 48            47               0                  0               0
#> 49            48               0                  0               0
#> 50            49               0                  0               0
#> 51            50               0                  0               0
#> 52            51               0                  0               0
#> 53            52               0                  0               0
#> 54            53               0                  0               0
#> 55            54               0                  0               0
#> 56            55               0                  0               0
#> 57            56               0                  0               0
#> 58            57               0                  0               0
#> 59            58               0                  0               0
#> 60            59               0                  0               0
#> 61            60               0                  0               0
#> 62            61               0                  0               0
#> 63            62               0                  0               0
#> 64            63               0                  0               0
#> 65            64               0                  0               0
#> 66            65               0                  0               0
#> 67            66               0                  0               0
#> 68            67               0                  0               0
#> 69            68               0                  0               0
#> 70            69               0                  0               0
#> 71            70               0                  0               0
#> 72            71               0                  0               0
#> 73            72               0                  0               0
#> 74            73               0                  0               0
#> 75            74               0                  0               0
#> 76            75               0                  0               0
#> 77            76               0                  0               0
#> 78            77               0                  0               0
#> 79            78               0                  0               0
#> 80            79               0                  0               0
#> 81            80               0                  0               0
#> 82            81               0                  0               0
#> 83            82               0                  0               0
#> 84            83               0                  0               0
#> 85            84               0                  0               0
#> 86            85               0                  0               0
#> 87            86               0                  0               0
#> 88            87               0                  0               0
#> 89            88               0                  0               0
#> 90            89               0                  0               0
#> 91            90               0                  0               0
#> 92            91               0                  0               0
#> 93            92               0                  0               0
#> 94            93               0                  0               0
#> 95            94               0                  0               0
#> 96            95               0                  0               0
#> 97            96               0                  0               0
#> 98            97               0                  0               0
#> 99            98               0                  0               0
#> 100           99               0                  0               0
#>     nbrMismatchReverse
#> 1                    0
#> 2                    0
#> 3                    6
#> 4                    0
#> 5                    0
#> 6                    0
#> 7                    0
#> 8                    0
#> 9                    0
#> 10                   0
#> 11                   0
#> 12                   0
#> 13                   0
#> 14                   0
#> 15                   9
#> 16                   0
#> 17                   0
#> 18                   0
#> 19                   0
#> 20                   0
#> 21                   0
#> 22                   0
#> 23                   0
#> 24                   0
#> 25                   0
#> 26                   0
#> 27                   0
#> 28                   0
#> 29                   0
#> 30                   0
#> 31                   0
#> 32                   0
#> 33                   0
#> 34                   3
#> 35                   0
#> 36                   0
#> 37                   0
#> 38                  15
#> 39                   0
#> 40                   0
#> 41                   0
#> 42                   0
#> 43                   0
#> 44                   0
#> 45                   0
#> 46                   0
#> 47                   0
#> 48                   0
#> 49                   0
#> 50                   0
#> 51                   0
#> 52                   0
#> 53                   0
#> 54                   0
#> 55                   0
#> 56                   0
#> 57                   0
#> 58                   0
#> 59                   0
#> 60                   0
#> 61                   0
#> 62                   0
#> 63                   0
#> 64                   0
#> 65                   0
#> 66                   0
#> 67                   0
#> 68                   0
#> 69                   0
#> 70                   0
#> 71                   0
#> 72                   0
#> 73                   0
#> 74                   0
#> 75                   0
#> 76                   0
#> 77                   0
#> 78                   0
#> 79                   0
#> 80                   0
#> 81                   0
#> 82                   0
#> 83                   0
#> 84                   0
#> 85                   0
#> 86                   0
#> 87                   0
#> 88                   0
#> 89                   0
#> 90                   0
#> 91                   0
#> 92                   0
#> 93                   0
#> 94                   0
#> 95                   0
#> 96                   0
#> 97                   0
#> 98                   0
#> 99                   0
#> 100                  0

## ---------------------------------------------------------------------- ## 
## Process a paired-end data set where the forward and reverse reads 
## contain variable regions corresponding to different proteins, and thus 
## should not be merged, specify the reads as a combination of 
## UMI, constant region and variable region (skip the first base), provide 
## the wild type sequence to compare the variable region to and limit the 
## number of allowed mutated codons to 1
out <- digestFastqs(
    fastqForward = system.file("extdata", "transInput_1.fastq.gz", 
                               package = "mutscan"),
    fastqReverse = system.file("extdata", "transInput_2.fastq.gz",
                               package = "mutscan"), 
    mergeForwardReverse = FALSE,  
    elementsForward = "SUCV", elementLengthsForward = c(1, 10, 18, 96),
    elementsReverse = "SUCV", elementLengthsReverse = c(1, 8, 20, 96),
    constantForward = "AACCGGAGGAGGGAGCTG", 
    constantReverse = "GAAAAAGGAAGCTGGAGAGA",
    wildTypeForward = c(FOS = paste0(
        "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTC", 
        "TGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA")),
    wildTypeReverse = c(JUN = paste0(
        "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTC", 
        "GGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT")), 
    nbrMutatedCodonsMaxForward = 1,
    nbrMutatedCodonsMaxReverse = 1
)
## Table with read counts and mutant information
head(out$summaryTable)
#>            mutantName
#> 1 FOS.0.WT_JUN.13.CCC
#> 2 FOS.0.WT_JUN.13.CTC
#> 3  FOS.0.WT_JUN.2.TCC
#> 4 FOS.0.WT_JUN.20.ACC
#> 5 FOS.0.WT_JUN.30.AGG
#> 6 FOS.0.WT_JUN.30.GGG
#>                                                                                                                                                                                            sequence
#> 1 ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA_ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAACCCCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT
#> 2 ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA_ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAACTCCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT
#> 3 ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA_ATCTCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT
#> 4 ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA_ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGACCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT
#> 5 ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA_ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGAGGCAGCTT
#> 6 ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA_ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGGGCAGCTT
#>   nbrReads maxNbrReads nbrUmis nbrMutBases nbrMutCodons nbrMutAAs varLengths
#> 1        1           1       1           2            1         1      96_96
#> 2        1           1       1           3            1         1      96_96
#> 3        1           1       1           1            1         1      96_96
#> 4        1           1       1           1            1         1      96_96
#> 5        1           1       1           3            1         1      96_96
#> 6        1           1       1           2            1         1      96_96
#>                        mutantNameBase     mutantNameCodon
#> 1          FOS.0.WT_JUN.37.C_JUN.39.C FOS.0.WT_JUN.13.CCC
#> 2 FOS.0.WT_JUN.37.C_JUN.38.T_JUN.39.C FOS.0.WT_JUN.13.CTC
#> 3                    FOS.0.WT_JUN.4.T  FOS.0.WT_JUN.2.TCC
#> 4                   FOS.0.WT_JUN.58.A FOS.0.WT_JUN.20.ACC
#> 5 FOS.0.WT_JUN.88.A_JUN.89.G_JUN.90.G FOS.0.WT_JUN.30.AGG
#> 6          FOS.0.WT_JUN.89.G_JUN.90.G FOS.0.WT_JUN.30.GGG
#>           mutantNameBaseHGVS      mutantNameAA       mutantNameAAHGVS
#> 1 FOS:c_JUN:c.37_39delinsCCC FOS.0.WT_JUN.13.P FOS:p_JUN:p.(Ala13Pro)
#> 2 FOS:c_JUN:c.37_39delinsCTC FOS.0.WT_JUN.13.L FOS:p_JUN:p.(Ala13Leu)
#> 3           FOS:c_JUN:c.4G>T  FOS.0.WT_JUN.2.S  FOS:p_JUN:p.(Ala2Ser)
#> 4          FOS:c_JUN:c.58T>A FOS.0.WT_JUN.20.T FOS:p_JUN:p.(Ser20Thr)
#> 5 FOS:c_JUN:c.88_90delinsAGG FOS.0.WT_JUN.30.R FOS:p_JUN:p.(Ala30Arg)
#> 6  FOS:c_JUN:c.89_90delinsGG FOS.0.WT_JUN.30.G FOS:p_JUN:p.(Ala30Gly)
#>   mutationTypes
#> 1 nonsynonymous
#> 2 nonsynonymous
#> 3 nonsynonymous
#> 4 nonsynonymous
#> 5 nonsynonymous
#> 6 nonsynonymous
#>                                                          sequenceAA
#> 1 TDTLQAETDQLEDEKSALQTEIANLLKEKEKL_IARLEEKVKTLKPQNSELASTANMLREQVAQL
#> 2 TDTLQAETDQLEDEKSALQTEIANLLKEKEKL_IARLEEKVKTLKLQNSELASTANMLREQVAQL
#> 3 TDTLQAETDQLEDEKSALQTEIANLLKEKEKL_ISRLEEKVKTLKAQNSELASTANMLREQVAQL
#> 4 TDTLQAETDQLEDEKSALQTEIANLLKEKEKL_IARLEEKVKTLKAQNSELATTANMLREQVAQL
#> 5 TDTLQAETDQLEDEKSALQTEIANLLKEKEKL_IARLEEKVKTLKAQNSELASTANMLREQVRQL
#> 6 TDTLQAETDQLEDEKSALQTEIANLLKEKEKL_IARLEEKVKTLKAQNSELASTANMLREQVGQL
## Filter summary
out$filterSummary
#>   nbrTotal f1_nbrAdapter f2_nbrNoPrimer f3_nbrReadWrongLength
#> 1     1000             0              0                     0
#>   f4_nbrNoValidOverlap f5_nbrAvgVarQualTooLow f6_nbrTooManyNinVar
#> 1                    0                     17                   0
#>   f7_nbrTooManyNinUMI f8_nbrTooManyBestWTHits f9_nbrMutQualTooLow
#> 1                   0                       0                   0
#>   f10a_nbrTooManyMutCodons f10b_nbrTooManyMutBases f11_nbrForbiddenCodons
#> 1                      698                       0                      0
#>   f12_nbrTooManyMutConstant f13_nbrTooManyBestConstantHits nbrRetained
#> 1                         0                              0         285
## Error statistics
out$errorStatistics
#>     PhredQuality nbrMatchForward nbrMismatchForward nbrMatchReverse
#> 1              0               0                  0               0
#> 2              1               0                  0               0
#> 3              2               0                  0               0
#> 4              3               0                  0               0
#> 5              4               0                  0               0
#> 6              5               0                  0               0
#> 7              6               0                  0               0
#> 8              7               0                  0               0
#> 9              8               0                  0               0
#> 10             9               0                  0               0
#> 11            10               0                  0               0
#> 12            11               0                  0               0
#> 13            12               0                  0               0
#> 14            13               0                  0               0
#> 15            14             160                 11             206
#> 16            15               0                  0               0
#> 17            16               0                  0               0
#> 18            17               0                  0               0
#> 19            18               0                  0               0
#> 20            19               0                  0               0
#> 21            20               0                  0               0
#> 22            21               0                  0               0
#> 23            22              53                  0              14
#> 24            23               0                  0               0
#> 25            24               0                  0               0
#> 26            25               0                  0               0
#> 27            26               0                  0               0
#> 28            27             307                  4             486
#> 29            28               0                  0               0
#> 30            29               0                  0               0
#> 31            30               0                  0               0
#> 32            31               0                  0               0
#> 33            32               0                  0               0
#> 34            33             486                  0             468
#> 35            34               0                  0               0
#> 36            35               0                  0               0
#> 37            36               0                  0               0
#> 38            37            4108                  1            4505
#> 39            38               0                  0               0
#> 40            39               0                  0               0
#> 41            40               0                  0               0
#> 42            41               0                  0               0
#> 43            42               0                  0               0
#> 44            43               0                  0               0
#> 45            44               0                  0               0
#> 46            45               0                  0               0
#> 47            46               0                  0               0
#> 48            47               0                  0               0
#> 49            48               0                  0               0
#> 50            49               0                  0               0
#> 51            50               0                  0               0
#> 52            51               0                  0               0
#> 53            52               0                  0               0
#> 54            53               0                  0               0
#> 55            54               0                  0               0
#> 56            55               0                  0               0
#> 57            56               0                  0               0
#> 58            57               0                  0               0
#> 59            58               0                  0               0
#> 60            59               0                  0               0
#> 61            60               0                  0               0
#> 62            61               0                  0               0
#> 63            62               0                  0               0
#> 64            63               0                  0               0
#> 65            64               0                  0               0
#> 66            65               0                  0               0
#> 67            66               0                  0               0
#> 68            67               0                  0               0
#> 69            68               0                  0               0
#> 70            69               0                  0               0
#> 71            70               0                  0               0
#> 72            71               0                  0               0
#> 73            72               0                  0               0
#> 74            73               0                  0               0
#> 75            74               0                  0               0
#> 76            75               0                  0               0
#> 77            76               0                  0               0
#> 78            77               0                  0               0
#> 79            78               0                  0               0
#> 80            79               0                  0               0
#> 81            80               0                  0               0
#> 82            81               0                  0               0
#> 83            82               0                  0               0
#> 84            83               0                  0               0
#> 85            84               0                  0               0
#> 86            85               0                  0               0
#> 87            86               0                  0               0
#> 88            87               0                  0               0
#> 89            88               0                  0               0
#> 90            89               0                  0               0
#> 91            90               0                  0               0
#> 92            91               0                  0               0
#> 93            92               0                  0               0
#> 94            93               0                  0               0
#> 95            94               0                  0               0
#> 96            95               0                  0               0
#> 97            96               0                  0               0
#> 98            97               0                  0               0
#> 99            98               0                  0               0
#> 100           99               0                  0               0
#>     nbrMismatchReverse
#> 1                    0
#> 2                    0
#> 3                    0
#> 4                    0
#> 5                    0
#> 6                    0
#> 7                    0
#> 8                    0
#> 9                    0
#> 10                   0
#> 11                   0
#> 12                   0
#> 13                   0
#> 14                   0
#> 15                  17
#> 16                   0
#> 17                   0
#> 18                   0
#> 19                   0
#> 20                   0
#> 21                   0
#> 22                   0
#> 23                   0
#> 24                   0
#> 25                   0
#> 26                   0
#> 27                   0
#> 28                   3
#> 29                   0
#> 30                   0
#> 31                   0
#> 32                   0
#> 33                   0
#> 34                   1
#> 35                   0
#> 36                   0
#> 37                   0
#> 38                   0
#> 39                   0
#> 40                   0
#> 41                   0
#> 42                   0
#> 43                   0
#> 44                   0
#> 45                   0
#> 46                   0
#> 47                   0
#> 48                   0
#> 49                   0
#> 50                   0
#> 51                   0
#> 52                   0
#> 53                   0
#> 54                   0
#> 55                   0
#> 56                   0
#> 57                   0
#> 58                   0
#> 59                   0
#> 60                   0
#> 61                   0
#> 62                   0
#> 63                   0
#> 64                   0
#> 65                   0
#> 66                   0
#> 67                   0
#> 68                   0
#> 69                   0
#> 70                   0
#> 71                   0
#> 72                   0
#> 73                   0
#> 74                   0
#> 75                   0
#> 76                   0
#> 77                   0
#> 78                   0
#> 79                   0
#> 80                   0
#> 81                   0
#> 82                   0
#> 83                   0
#> 84                   0
#> 85                   0
#> 86                   0
#> 87                   0
#> 88                   0
#> 89                   0
#> 90                   0
#> 91                   0
#> 92                   0
#> 93                   0
#> 94                   0
#> 95                   0
#> 96                   0
#> 97                   0
#> 98                   0
#> 99                   0
#> 100                  0
```
