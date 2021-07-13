
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

See the [vignette](https://tanaylab.github.io/MCview/articles/getting-started.html).

Tl;dr:

1.  initialize a new project:

``` r
MCView::create_project("PBMC")
```

1.  Edit `PBMC/config/config.yaml`.

See the [config vignette](https://tanaylab.github.io/MCview/articles/config.html) for more details.

1.  Import the project:

``` r
MCView::import("PBMC")
```

1.  Run the app:

You can run the app from R using:

``` r
MCView::run_app(project = "PBMC", port = 5555, host = "127.0.0.1")
```

Or from shell:

``` bash
Rscript start_app.R PBMC
```
