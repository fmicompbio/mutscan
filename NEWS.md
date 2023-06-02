# mutscan 0.2.37

* accommodate changes in edgeR version 3.42 (Bioconductor 3.17)
* mutscan manuscript published and available at https://doi.org/10.1186/s13059-023-02967-0

# mutscan 0.2.36

* Check that reads don't exceed the maximal allowed length
* Add parameter to specify maximal read length

# mutscan 0.2.35

* Add alternative names for variants (including HGVS identifiers)

# mutscan 0.2.34

* Expand examples in function documentation

# mutscan 0.2.33

* Replace Matrix.utils::aggregate.Matrix (removed from CRAN) by DelayedArray::rowsum

# mutscan 0.2.32

* Expand FASTQ file paths automatically
* Filter out reads where the variable region is longer than the best matching WT sequence

# mutscan 0.2.31

* Add option to include identity line in pairs plots
* Select points to label in result plots based on nominal p-value
* Add option to manually set the correlation range for coloring pairs plots

# mutscan 0.2.30

* Swap to MIT license

# mutscan 0.2.29

* Include degrees of freedom in output from calculateRelativeFC

# mutscan 0.2.28

* Rename calculatePPIScore to calculateFitnessScore
* Add amino acid-related information (AA sequence, mutant name, number of mutated amino acid) to digestFastqs output, propagate through summarizeExperiment and collapseMutantsByAA. 

# mutscan 0.2.27

* Add linkMultipleVariants function

# mutscan 0.2.26

* Use arithmetic instead of geometric mean for PPI calculations

# mutscan 0.2.25

* Only add names to WT sequences if they are unnamed

# mutscan 0.2.24

* Added function to calculate nearest string distances

# mutscan 0.2.23

* Changed default value of forbidden codons to ""
* Make result plots more flexible
* Propagate information about the number of mutated bases/codons in `summarizeExperiment()`

# mutscan 0.2.22

* Added argument to specify the title of the QC report.
* Internal refactoring of argument checking.

# mutscan 0.2.21

* Added a `NEWS.md` file to track changes to the package.
* Added the `generateQCReport()` function.
* Added functions to plot results from `calculateRelativeFC()`.
* Added functions to plot distributions of counts as well as total counts.
