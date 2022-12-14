
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pearsonCorrelateCells

<!-- badges: start -->
<!-- badges: end -->

`pearsonCorrelateCells` is an R package that works with [10X Visium
Spatial
Transcriptomic](https://www.10xgenomics.com/products/spatial-gene-expression)
data. The package is developed to allow the identification and
visualization of cells that exhibit gene expression patterns correlated
with that of a comparator. Correlations are calculated using Pearson
Product-Moment Correlation. The comparator is a dataframe that contains
gene names and corresponding expression values.

Gene expression data in the comparator can be taken from any source.
Here are some anticipated use cases:

1)  The comparator is averaged or differential expression of a cluster
    in RNA-sequencing data. The goal is to find cells in 10X Visium that
    are most similar to this cluster.

2)  The comparator is averaged gene expression of samples of interest in
    bulk RNA-sequencing data

The idea of this package came from a conversation with Professor Scott
Yuzwa (personal communication, September 23, 2022). According to
Professor Yuzwa, some researchers use pearson correlations to identify
neural stem cells in spatial transcriptomic data of brains of mice with
a stroke.

This functions in this package work with 10X Visium Data that is loaded
as a Seurat object. To learn how to load 10X Visium datasets as Seurat
objects, please visit the [Seurat package vignette for working with 10X
Visium](https://satijalab.org/seurat/articles/spatial_vignette.html).

## Installation

You can install the development version of pearsonCorrelateCells from
[GitHub](https://github.com/) with:

``` r
require("devtools")
devtools::install_github("jessielam415/pearsonCorrelateCells", build_vignettes = TRUE)
library("pearsonCorrelateCells")
```

The shiny app in this package demonstrates the package’s utility by
allowing users to upload csv files of gene expression data. The gene
expression data is correlated with each spot in a sample Visium spatial
transcriptomic dataset. The chosen dataset is the 10X Genomics Visium
Mouse Brain Dataset. To run the shinyApp:

``` r
pearsonCorrelateCells::runPearsonCorrelateCells()
```

## Overview

``` r
ls("package:pearsonCorrelateCells")
browseVignettes("pearsonCorrelateCells")
```

`pearsonCorrelateCells` has 2 functions:

-   `findCorrelatedCells`: Computes the pearson correlation of the gene
    expression of the top expressed genes in a comparator dataframe and
    each spot in 10X Visium spatial dataset. Saves the p-value and
    correlation coefficient in the Visium object metadata for each spot.
-   `visualizeCells`: Visualizes the significance of the pearson
    correlations of each spot obtained from the findCorrelatedCells
    function using the p-value \<= 0.05 cutoff. Outputs a plot of each
    visium spot coloured by significance overlaying the Visium image.

![](./inst/extdata/pearsonCorrelateCellsWorkflow.png)

## Contributions

The author of the package is Wing Chung Jessie Lam.

The *findCorrelatedCells* function makes use of the \[\] functions from
the \[\] r package to \[\]

The *visualizeCells* function makes use the \[\] functions from the \[\]
r package to \[\]

## References

Bache S, Wickham H (2022). magrittr: A Forward-Pipe Operator for R.
<https://magrittr.tidyverse.org>,
<https://github.com/tidyverse/magrittr>.

Hahsler M, Nagar A (2019). rBLAST: R Interface for the Basic Local
Alignment Search Tool. R package version 0.99.2.
<https://github.com/mhahsler/rBLAST>.

Silva, A. (2022). Anjalisilva/TestingPackage: A Simple R Package
Illustrating Components of an R Package: 2019-2022 BCB410H - Applied
Bioinformatics, University of Toronto, Canada. GitHub.
<https://github.com/anjalisilva/TestingPackage>

Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis.
Springer-Verlag New York. ISBN 978-3-319-24277-4.
<https://ggplot2.tidyverse.org>

Wickham H, François R, Henry L, Müller K (2022). dplyr: A Grammar of
Data Manipulation. <https://dplyr.tidyverse.org>,
<https://github.com/tidyverse/dplyr>.

## Acknowledgements

This package was developed as part of an assessment for 2022 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `ComparePseudogenes` welcomes issues, enhancement requests, and
other contributions. To submit an issue, use the GitHub issues.
