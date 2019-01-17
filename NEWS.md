# MicrobiomeR 0.2.3

## Documentation Updates

* Updated documentation for `create_pub_table` function.

## Tests Added

* Added tests for utils.R

# MicrobiomeR 0.2.2

## Functional Changes

* Added support for multiple treatments for some plots
  * `get_heat_tree_plots` now produces a `metacoder::heat_tree_matrix`.
  * `correlation_plot` now produces multiple plots for data with more that 2 treatments.
  * Added the `get_treatment_matrix` function
  * Added the `get_correlation_data` function
* Updated `analysis` vignette to demonstrate more than 2 treatment groups.
* Updated internal data files.
  * Added new metadata file with 3 treatment groups (`nephele_metadata3.txt`).
  * Formatted Treatment Group metadata.
* Updated public datasets.
  * Removed Greengenes datasets.
  * Added Silva data with 3 treatment groups.
* Added the `color-palettes` vignette.

# MicrobiomeR 0.2.1

## Features Added

* Added permanova function for generating permanova stats.
* Added top_coefficients_barplot for generating a plot of the top coefficients output from the permanova function.

## Tests Added

* Added tests for permanova.R.

# MicrobiomeR 0.2.0

## Documentation Updates

* Added vignettes to package and pkgdown config.
  * About
  * Introduction
  * Data Wrangling
  * Filtering
  * Analysis
* Added reference sections to pkgdown config corresponding to @family tag.
  * Import
  * Formatting
  * Validation
  * Filtering (Basic, Advanced, Other)
  * Visualization
  * Color Palettes

## Functional Changes

* Renamed `sample_filter` to `sample_id_filter`.
* Renamed `get_phyloseq_obj` to `create_phyloseq`.
* Renamed `format_metacoder_object` to `as_custom_format`.
* Updated `vlookup` and `get_color_palette` parameters.
* Added the `root_phyloseq_tree` function.
* Fixed minor bugs in code and docs.

# MicrobiomeR 0.1.3
 
* Added custom css for pkgdown.
* Added more sections including authors and vignettes to _pkgdown.yml.
* Fixed `@importFrom` sections for crayon.

# MicrobiomeR 0.1.2
 
* Added `ordination.R` for ordination plots.
* Added tests for ordination plots.

# MicrobiomeR 0.1.1
 
* Added loggging and better warning messages.

# MicrobiomeR 0.1.0
 
* Changed to a semantic versioning scheme.
* Added `barplot.R` file of stacked barplot functions.
* Improved `get_alpha_diversity_measures` function.
* Added a `NEWS.md` file to track changes to the package.
* Added simple tests for `stacked_barplot` and `alpha_diversity_plot`.

# MicrobiomeR 0.0.9.2

* Fixes to heat tree plots

# MicrobiomeR 0.0.9.1

* Added `barplot.R` file of stacked barplot functions.
