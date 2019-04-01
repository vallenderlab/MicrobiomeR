[![Build Status](https://travis-ci.com/vallenderlab/MicrobiomeR.svg?branch=master)](https://travis-ci.com/vallenderlab/MicrobiomeR)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Coverage status](https://codecov.io/gh/vallenderlab/MicrobiomeR/branch/master/graph/badge.svg)](https://codecov.io/github/vallenderlab/MicrobiomeR?branch=master)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.01299/status.svg)](https://doi.org/10.21105/joss.01299)

# MicrobiomeR <a href='https://microbiomer.vallenderlab.science/'><img src="man/figures/logo.png" align="right" height=300/></a>

An R package for microbiome analysis that incorporates [phyloseq](https://github.com/joey711/phyloseq), 
[metacoder](https://github.com/grunwaldlab/metacoder), [taxa](https://github.com/ropensci/taxa), and [microbiome](https://github.com/microbiome/microbiome/) in order to standardize and simplify common microbiome workflows.

## Installation

We are currently not on CRAN or Bioconductor:

```r
library(devtools) # Load the devtools package
install_github("vallenderlab/MicrobiomeR") # Install the package
```

Please visit https://microbiomer.vallenderlab.science/ for extensive documentation of the package.

## Workflow Features

- Standardization of data wrangling.
    - Phyloseq for data import.
    - Taxa for the primary data object (**Taxmap**).
    - Proprietary data formatting and validation.
- Phyloseq inspired filtering for `taxa::taxmap` objects.
    - Metacoder/taxa for mainstream filtering.
    - Proprietary basic filtering for samples, taxonomies, and OTUs.
    - Proprietary advanced filtering (phyloseq-style).
    - Other Proprietary filtering functions for observation data.
- Metacoder enabled statistical analysis functions.
- Various visualization options.
    - Output Directories
    - Color Palettes
    - Heat Trees
    - Correlation Plots
    - Stacked Bar Plot
    - Alpha Diversity Plot
    - Ordination Plot

## Contributing to `MicrobiomeR`

* Our code style is based on Google's R style developed by Hadley Wickham.
* Please note that the `MicrobiomeR` project is released with a [Contributor Code of Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.
* View our [contributing](.github/CONTRIBUTING.md) documentation for more details.

