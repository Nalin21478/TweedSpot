# TweedSpot

TweedSpot is an R package for detecting spatially variable genes (SVGs) in
spatial omics data with spatial generalized additive models using a Tweedie
count family.

The package supports three analysis modes:

- agnostic SVG detection across all spatial locations
- cell-type-specific SVG detection using spot-level cell-type proportions
- cell-type-specific SVG detection using per-location cell-type labels

## Installation

TweedSpot is intended to be installed from GitHub.

```r
install.packages("remotes")
remotes::install_github("himelmallick/TweedSpot", dependencies = TRUE)
```

Because TweedSpot depends on Bioconductor infrastructure, install missing
Bioconductor packages first if needed:

```r
install.packages("BiocManager")
BiocManager::install(c(
  "SpatialExperiment",
  "SummarizedExperiment",
  "SingleCellExperiment",
  "S4Vectors",
  "BiocParallel",
  "scran",
  "scuttle",
  "STexampleData"
))
```

## Quick Start

```r
library(TweedSpot)
library(SpatialExperiment)
library(STexampleData)

spe <- Visium_humanDLPFC()
spe <- spe[, colData(spe)$in_tissue == 1]

spe_res <- tweedspot(
  input = spe,
  assay_name = "counts"
)

head(SummarizedExperiment::rowData(spe_res)[, c(
  "tweedspot_stat",
  "tweedspot_pval",
  "tweedspot_padj"
)])
```

## Example With Filtering

```r
library(TweedSpot)
library(SpatialExperiment)
library(STexampleData)
library(nnSVG)
library(scran)

spe <- Visium_humanDLPFC()
spe <- spe[, colData(spe)$in_tissue == 1]

spe <- filter_genes(
  spe,
  filter_genes_ncounts = 5,
  filter_genes_pcspots = 5
)

spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

spe_res <- tweedspot(
  input = spe,
  assay_name = "counts",
  two_part = FALSE,
  family = "tw",
  fit_method = "REML",
  use_bam = TRUE,
  bam_discrete = TRUE,
  smooth_k = 20
)
```

## Main Function

TweedSpot currently exposes a single user-facing function:

- `tweedspot()`: run agnostic or cell-type-specific SVG detection on a
  `SpatialExperiment`

Key arguments:

- `input`: a `SpatialExperiment`
- `assay_name`: assay to model, typically `"counts"`
- `W`: optional spot-by-cell-type proportion matrix
- `celltype`: optional cell-type labels
- `two_part`: optional two-part agnostic model
- `use_bam`: use `mgcv::bam()` for faster fitting on larger datasets
- `smooth_k`: optional smooth basis size to trade flexibility for speed

## Output

In agnostic mode, TweedSpot writes these per-gene fields to `rowData(input)`:

- `tweedspot_stat`
- `tweedspot_pval`
- `tweedspot_padj`

In cell-type-specific modes, TweedSpot stores genes-by-cell-type p-value and
adjusted p-value matrices in `rowData(input)`, and a long-format results table
in `metadata(input)$tweedspot`.

## Status

The package is under active development. A vignette and broader test coverage
can be added on top of this installable package skeleton.

## License

MIT
