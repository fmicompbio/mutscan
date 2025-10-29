# Generate QC report

Generate QC report

## Usage

``` r
generateQCReport(
  se,
  outFile,
  reportTitle = "mutscan QC report",
  forceOverwrite = FALSE,
  ...
)
```

## Arguments

- se:

  A `SummarizedExperiment` object, typically generated with
  [`summarizeExperiment()`](https://fmicompbio.github.io/mutscan/reference/summarizeExperiment.md).

- outFile:

  Character string providing the name of the output file. Should have
  the extension `.html`.

- reportTitle:

  Character string specifying the title of the QC report.

- forceOverwrite:

  Logical scalar, indicating whether an existing file with the same name
  as `outFile` should be overwritten.

- ...:

  Additional parameters to be forwarded to
  [`render`](https://pkgs.rstudio.com/rmarkdown/reference/render.html),
  for example `quiet = TRUE`.

## Value

Invisibly, the path to the generated html file.

## See also

[`render`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)
used to render the html output file.

## Author

Charlotte Soneson

## Examples

``` r
## Load SummarizedExperiment object
se <- readRDS(system.file("extdata", "GSE102901_cis_se.rds",
                          package = "mutscan"))
## Define output file
outfile <- tempfile(fileext = ".html")

## Generate QC report
generateQCReport(se, outfile)
```
