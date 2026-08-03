# Calculate logFCs relative to WT using edgeR

Calculate logFCs and associated p-values for a given comparison, using
either limma or the Negative Binomial quasi-likelihood framework of
edgeR. The observed counts for the WT variants can be used as offsets in
the model.

## Usage

``` r
calculateRelativeFC(
  se,
  design,
  coef = NULL,
  contrast = NULL,
  WTrows = NULL,
  selAssay = "counts",
  pseudocount = 1,
  method = "edgeR",
  normMethod = ifelse(is.null(WTrows), "TMM", "sum")
)
```

## Arguments

- se:

  SummarizedExperiment object.

- design:

  Design matrix. The rows of the design matrix must be in the same order
  as the columns in `se`.

- coef:

  Coefficient(s) to test with edgeR or limma.

- contrast:

  Numeric contrast to test with edgeR or limma.

- WTrows:

  Vector of row names that will be used as the reference when
  calculating logFCs and statistics. If more than one value is provided,
  the sum of the corresponding counts is used to generate offsets. If
  `NULL`, offsets will be defined as the effective library sizes (using
  TMM normalization factors).

- selAssay:

  Assay to select from `se` for the analysis.

- pseudocount:

  Pseudocount to add when calculating log-fold changes.

- method:

  Either 'edgeR' or 'limma'. If set to 'limma', voom is used to
  transform the counts and estimate observation weights before applying
  limma. In this case, the results also contain the standard errors of
  the logFCs.

- normMethod:

  Character scalar indicating which normalization method should be used
  to calculate size factors. Should be either `"TMM"` or `"csaw"` when
  `WTrows` is `NULL`, and `"geomean"` or `"sum"` when `WTrows` is
  provided.

## Value

A `data.frame` with output from the statistical testing framework (edgeR
or limma).

## Author

Charlotte Soneson, Michael Stadler

## Examples

``` r
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following object is masked from ‘package:utils’:
#> 
#>     data
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, scale, sequence, table,
#>     tapply, transform, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
se <- readRDS(system.file("extdata", "GSE102901_cis_se.rds",
                          package = "mutscan"))[1:200, ]
design <- model.matrix(~ Replicate + Condition,
                       data = colData(se))
                              
## Calculate "absolute" log-fold changes with edgeR
res <- calculateRelativeFC(se, design, coef = "Conditioncis_output", 
                           method = "edgeR")
head(res)
#>             logFC    logCPM         F       PValue          FDR logFC_shrunk
#> f.0.WT  1.9721335 20.083627 59.456865 1.858697e-11 3.860020e-10    1.9721331
#> f.1.AAC 1.3121715  4.688782 16.841571 9.097826e-05 3.084009e-04    1.2983948
#> f.1.AAG 1.0388250  5.677955 12.906277 5.405369e-04 1.637991e-03    1.0342284
#> f.1.ACC 2.0874104 12.862431 62.641335 7.172371e-12 1.793093e-10    2.0873513
#> f.1.ACG 1.8698754 13.513672 51.069321 2.539192e-10 2.821324e-09    1.8698449
#> f.1.AGC 0.9302662  4.754181  8.738407 4.004634e-03 1.013831e-02    0.9228443
#>         df.total df.prior df.test
#> f.0.WT  87.43786 84.43757       1
#> f.1.AAC 87.43697 84.43757       1
#> f.1.AAG 87.43763 84.43757       1
#> f.1.ACC 87.43786 84.43757       1
#> f.1.ACG 87.43786 84.43757       1
#> f.1.AGC 87.43719 84.43757       1
## Calculate log-fold changes relative to the WT sequence with edgeR
stopifnot("f.0.WT" %in% rownames(se))
res <- calculateRelativeFC(se, design, coef = "Conditioncis_output", 
                           method = "edgeR", WTrows = "f.0.WT")
head(res)
#>                 logFC    logCPM         F       PValue          FDR
#> f.0.WT   4.209908e-16 19.757846  0.000000 1.000000e+00 1.000000e+00
#> f.1.AAC -6.419635e-01  4.553595  9.542087 1.063414e-02 1.321011e-02
#> f.1.AAG -8.818448e-01  5.655807 24.029119 5.109289e-04 9.461646e-04
#> f.1.ACC  1.207040e-01 12.506900 29.621010 2.244723e-04 4.628294e-04
#> f.1.ACG -1.016662e-01 13.217924 66.399725 6.572830e-06 2.590725e-05
#> f.1.AGC -1.014998e+00  4.738963 19.367513 1.136289e-03 1.832724e-03
#>          logFC_shrunk df.total df.prior df.test
#> f.0.WT   5.110403e-16 10.69190 7.691777       1
#> f.1.AAC -6.667144e-01 10.69040 7.691777       1
#> f.1.AAG -8.824063e-01 10.69151 7.691777       1
#> f.1.ACC  1.207002e-01 10.69190 7.691777       1
#> f.1.ACG -1.016641e-01 10.69190 7.691777       1
#> f.1.AGC -1.031223e+00 10.69092 7.691777       1

## Calculate "absolute" log-fold changes with limma
res <- calculateRelativeFC(se, design, coef = "Conditioncis_output", 
                           method = "limma")
head(res)
#>             logFC      CI.L     CI.R   AveExpr         t      P.Value
#> f.0.WT  1.9474812 1.5678409 2.327121 19.757630 10.074558 3.681165e-22
#> f.1.AAC 1.2888967 0.6658454 1.911948  4.417300  4.062743 5.494443e-05
#> f.1.AAG 1.0502295 0.5417248 1.558734  5.581250  4.056153 5.647426e-05
#> f.1.ACC 2.0651448 1.6200963 2.510193 12.504065  9.113142 1.188299e-18
#> f.1.ACG 1.8429080 1.4051989 2.280617 13.217517  8.268811 8.735619e-16
#> f.1.AGC 0.9204368 0.3160363 1.524837  4.627474  2.990847 2.896169e-03
#>            adj.P.Val         B  se.logFC df.total df.prior
#> f.0.WT  1.840583e-20 43.331488 0.1933069      600      Inf
#> f.1.AAC 1.978864e-04  1.551562 0.3172479      600      Inf
#> f.1.AAG 1.981553e-04  1.382359 0.2589226      600      Inf
#> f.1.ACC 3.395141e-17 34.216900 0.2266117      600      Inf
#> f.1.ACG 1.247946e-14 26.953437 0.2228746      600      Inf
#> f.1.AGC 7.723116e-03 -2.174097 0.3077512      600      Inf
## Calculate log-fold changes relative to the WT sequence with limma
stopifnot("f.0.WT" %in% rownames(se))
res <- calculateRelativeFC(se, design, coef = "Conditioncis_output", 
                           method = "limma", WTrows = "f.0.WT")
head(res)
#>                 logFC        CI.L        CI.R   AveExpr             t
#> f.0.WT   9.531162e-08 -0.03804257  0.03804276 19.757630  5.947923e-06
#> f.1.AAC -6.025567e-01 -1.20358142 -0.00153195  4.417300 -2.380104e+00
#> f.1.AAG -8.904935e-01 -1.41357231 -0.36741467  5.581250 -4.041607e+00
#> f.1.ACC  1.233362e-01  0.05177181  0.19490069 12.504065  4.091513e+00
#> f.1.ACG -1.001837e-01 -0.12510198 -0.07526547 13.217517 -9.544871e+00
#> f.1.AGC -9.808048e-01 -1.66423494 -0.29737457  4.627474 -3.407053e+00
#>              P.Value    adj.P.Val          B   se.logFC df.total df.prior
#> f.0.WT  9.999954e-01 0.9999954234 -9.7096960 0.01602435 6.865421 3.865421
#> f.1.AAC 4.956108e-02 0.0563194038 -4.8577614 0.25316397 6.865421 3.865421
#> f.1.AAG 5.126411e-03 0.0080100165 -2.7505519 0.22033155 6.865421 3.865421
#> f.1.ACC 4.817458e-03 0.0076467589 -4.9657926 0.03014441 6.865421 3.865421
#> f.1.ACG 3.289888e-05 0.0002741573  0.1447601 0.01049608 6.865421 3.865421
#> f.1.AGC 1.167439e-02 0.0159923207 -3.3884774 0.28787484 6.865421 3.865421
```
