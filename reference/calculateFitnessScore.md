# Calculate fitness scores.

Using sequence counts before and after selection, calculate fitness
scores as described by Diss and Lehner (2018).

## Usage

``` r
calculateFitnessScore(
  se,
  pairingCol,
  ODCols,
  comparison,
  WTrows,
  selAssay = "counts"
)
```

## Arguments

- se:

  SummarizedExperiment object as returned by
  [`summarizeExperiment`](https://fmicompbio.github.io/mutscan/reference/summarizeExperiment.md).

- pairingCol:

  Name of column in `colData(se)` with replicate/pairing information.
  Samples with the same value in this column will be paired.

- ODCols:

  Name(s) of column(s) in `colData(se)` with OD values (numeric), used
  to normalize for different numbers of cells.

- comparison:

  3-element character vector of the form
  `(column, numerator, denominator)`. `column` is the name of the column
  in `colData(se)` with experimental conditions. `numerator` and
  `denominator` define the comparison, e.g.
  `c("cond", "output", "input")` will look in the `"cond"` column and
  calculate fitness for the ratio of `"output"` over `"input"` counts.

- WTrows:

  Vector of row names that will be used as the reference when
  calculating fitness scores. If more than one value is provided, the
  average of the corresponding fitness scores is used as a reference. If
  NULL, no division by WT scores will be done.

- selAssay:

  Assay to select from `se` for the analysis.

## Value

A numeric vector with fitness scores.

## References

"The genetic landscape of a physical interaction." Diss G and Lehner B.
Elife. 2018;7:e32472. doi: 10.7554/eLife.32472.

## Author

Michael Stadler and Charlotte Soneson

## Examples

``` r
se <- readRDS(system.file("extdata", "GSE102901_cis_se.rds", 
                          package = "mutscan"))
## Check that the wildtype sequence is present in the data
stopifnot("f.0.WT" %in% rownames(se))
## Calculate PPI scores as defined in Diss & Lehner (2018)
ppis <- calculateFitnessScore(
    se = se, pairingCol = "Replicate", 
    ODCols = c("OD1", "OD2"),
    comparison = c("Condition", "cis_output", "cis_input"),
    WTrows = "f.0.WT")
## Matrix with PPI scores for each replicate
head(ppis)
#>         cis_output_vs_cis_input_repl1 cis_output_vs_cis_input_repl2
#> f.0.WT                      1.0000000                     1.0000000
#> f.1.AAC                     0.9763609                     0.9173416
#> f.1.AAG                     0.8887226                     0.8775120
#> f.1.ACC                     1.0105154                     1.0096229
#> f.1.ACG                     0.9881472                     0.9905645
#> f.1.AGC                     0.8233725                     0.8718729
#>         cis_output_vs_cis_input_repl3
#> f.0.WT                      1.0000000
#> f.1.AAC                     0.8861805
#> f.1.AAG                     0.9407093
#> f.1.ACC                     1.0188111
#> f.1.ACG                     0.9884163
#> f.1.AGC                     0.9669444
```
