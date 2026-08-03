# Calculate distances to the nearest string

Given a character vector, calculate the distance for each element to the
nearest neighbor amongst all the other elements.

## Usage

``` r
calcNearestStringDist(x, metric = "hamming", nThreads = 1L)
```

## Arguments

- x:

  A character vector.

- metric:

  A character scalar defining the string distance metric. One of
  `"hamming"` (default), `"hamming_shift"` or `"levenshtein"`.

- nThreads:

  numeric(1), number of threads to use for parallel processing.

## Value

An integer vector of the same length as `x`.

## Examples

``` r
calcNearestStringDist(c("lazy", "hazy", "cozy"))
#> [1] 1 1 2
calcNearestStringDist(c("lazy", "hazy", "cozy"), metric = "hamming_shift")
#> [1] 1 1 2
calcNearestStringDist(c("lazy", "hazy", "crazy"), metric = "levenshtein")
#> [1] 1 1 2
```
