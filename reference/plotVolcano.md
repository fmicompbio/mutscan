# Construct a volcano plot

Construct a volcano plot

## Usage

``` r
plotVolcano(
  res,
  logFCCol = NULL,
  pvalCol = NULL,
  padjCol = NULL,
  padjThreshold = 0.05,
  pointSize = "small",
  interactivePlot = FALSE,
  nTopToLabel = 0
)
```

## Arguments

- res:

  `data.frame` (typically output from
  [`calculateRelativeFC()`](https://fmicompbio.github.io/mutscan/reference/calculateRelativeFC.md))
  with columns corresponding to the log-fold change (`logFC`), p-value
  (`PValue` or `P.Value`) and significance (`FDR` or `adj.P.Val`).

- logFCCol, pvalCol, padjCol:

  Character scalars indicating the columns from `res` that will be used
  to represent the logFC (x-axis), p-value (y-axis) and adjusted p-value
  (used for coloring). If `NULL` (default), pre-specified values will be
  used depending on the available columns (`"logFC"`, `"PValue"` or
  `"P.Value"`, and `"FDR"` or `"adj.P.Val"`, respectively).

- padjThreshold:

  Numeric scalar indicating the adjusted p-value threshold to use for
  coloring the points. All features with adjusted p-value below the
  treshold will be shown in red.

- pointSize:

  Either `"small"` or `"large"`, indicating which of the two available
  plot styles that will be used.

- interactivePlot:

  Logical scalar, indicating whether an interactive plot should be
  returned, in which one can hover over the individual points and obtain
  further information.

- nTopToLabel:

  Numeric scalar, indicating the number of points that should be labeled
  in the plot. The points will be ranked by the `pvalCol` column, and
  the top `nTopToLabel` values will be labeled by the corresponding row
  names. Only used if `interactivePlot` is `FALSE`.

## Value

If `interactivePlot` is `TRUE`, a `plotly` object. If `interactivePlot`
is `FALSE`, a `ggplot2` object.

## Author

Charlotte Soneson

## Examples

``` r
library(SummarizedExperiment)
se <- readRDS(system.file("extdata", "GSE102901_cis_se.rds",
                          package = "mutscan"))[1:200, ]
design <- model.matrix(~ Replicate + Condition,
                       data = colData(se))
res <- calculateRelativeFC(se, design, coef = "Conditioncis_output")
plotVolcano(res, pointSize = "large", nTopToLabel = 3)

```
