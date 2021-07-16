
<!-- README.md is generated from README.Rmd. Please edit that file -->
# MCView

<!-- badges: start -->
<!-- badges: end -->
The goal of MCView is to provide an interactive API of metacell objects analysis.

## Installation

``` r
remotes::install_github("tanaylab/MCView")
```

## Usage

See the [vignette](https://tanaylab.github.io/MCView/articles/getting-started.html).

Tl;dr:

#### Initialize a new project:

``` r
MCView::create_project("PBMC")
```

#### Import a dataset:

``` r
MCView::import_dataset(project = "PBMC", 
               dataset = "PBMC163k", 
               anndata_file = "raw/metacells.h5ad")
)
```

#### Run the app:

You can run the app from R using:

``` r
MCView::run_app("PBMC", launch.browser = TRUE)
```

Or from shell:

``` bash
Rscript start_app.R PBMC
```