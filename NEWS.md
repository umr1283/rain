# rain 0.3.2

## Minor improvements and fixes

* Use `patchwork` for ethnicity plot.
* In `/R/compute_pca.R`, 
    - fix a missplaced `.data` pronoun.
    - core code rewrite.
* In `/R/pca_report.R`, add packages prefix for `tidyselect`.
* Remove `gcta64` from system requirements.
* Remove dependency to `R.utils`.

# rain 0.3.1

## Minor improvements and fixes

* Use `file.path()` instead of `paste()`.
* Change messages prefix.

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
