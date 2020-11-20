# rain (development version)

# rain 0.5.1

## Minor text improvements

* In `R/pca_report.R`,
    + Remove unnecessary lower condition for outliers.
    + Update obsolete captions.
    + Perform independent association testing between technical variables and principal components.
    + Order association plot according to pvalues and principal components.
* In `R/compute_pca.R`,
    + Fix explained variance histogram.

# rain 0.5.0

## Breaking Changes

* In `R/pca_report.R`,
    + Code refactoring using `data.table`.
* In `R/compute_pca.R`,
    + Code refactoring using `data.table`.
    + Tweak ethnicity figure.
    + Add an ethnicity figure focussing on the studied cohort.

# rain 0.4.1

## Minor text improvements

* In `R/pca_report.R`,
    + Change text for outliers in the last part of the PCA report.

# rain 0.4.0

## Minor improvements and fixes

* In `DESCRIPTION`,
    + Use `patchwork` for ethnicity plot.
    + Remove `gcta64` from system requirements.
* In `R/compute_pca.R`,
    + Use `R.utils` to count number of SNPs from `*.bim` file.
    + Fix number of SNPs in ggplot legend.
    + Fix a missplaced `.data` pronoun.
    + Core code rewrite.
    + Change messages prefix.
* In `R/pca_report.R`, 
    + Add packages prefix for `tidyselect`.
    + Update to `ggplot2` `v3.3.0`.

# rain 0.3.0

## Minor improvements and fixes

- Remove `covr`.
- Remove `testthat`.

# rain 0.2.0

## New features

* `pca_report()` allows to compute an analysis report using principal
    component analysis from
    [flashpca](https://github.com/gabraham/flashpca) tool.  
    The function can be used in a chunk within a Rmarkdown
    document/script with `results="asis"` to render the report 
    (from [CARoT v0.4.0](https://github.com/omicsr/CARoT/tree/v0.4.0)).

# rain 0.1.0

## New features

* `estimate_ethnicity()` allows to format VCF files and compute the
    genomic components (and some figures) for ethnicity. ([CARoT v0.4.0](https://github.com/omicsr/CARoT/tree/v0.4.0))
