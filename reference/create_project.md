# Create a configuration folder for a project

Create a project directory with the default configuration files and
directory structure. If `edit_config == TRUE`, a text editor would be
opened in order to edit the project `config.yaml` file.

## Usage

``` r
create_project(
  project,
  title = "MCView",
  tabs = NULL,
  selected_gene1 = NULL,
  selected_gene2 = NULL,
  selected_mc1 = NULL,
  selected_mc2 = NULL,
  datasets = NULL,
  other_params = NULL,
  edit_config = TRUE,
  atlas = FALSE
)
```

## Arguments

- project:

  path of the project

- title:

  The title of the app. This would be shown on the top left of the
  screen.

- tabs:

  Controls which tabs to show in the left sidebar and their order.
  Options are: "QC", "Projection-QC", "Manifold", "Genes", "Query",
  "Atlas", "Markers", "Gene modules", "Projected-fold", "Diff.
  Expression", "Cell types", "Flow", "Annotate", "About". When NULL -
  default tabs would be set. For projects with atlas projections, please
  set `atlas` to TRUE.

- selected_gene1, selected_gene2:

  The default genes that would be selected (in any screen with gene
  selection). If this parameter is missing, the 2 genes with highest
  max(expr)-min(expr) in the first dataset would be chosen.

- selected_mc1, selected_mc2:

  The default metacells that would be selected in the Diff. Expression
  tab.

- datasets:

  A named list with additional per-dataset parameters. Current
  parameters include default visualization properties of projection and
  scatter plots.

- other_params:

  Named list of additional parameters such as projection_point_size,
  projection_point_stroke, scatters_point_size and scatters_stroke_size

- edit_config:

  open file editor for config file editing

- atlas:

  use default configuration for atlas projections (relevant only when
  `tabs` is NULL)

## Examples

``` r
if (FALSE) { # \dontrun{
dir.create("raw")
download.file(
    "http://www.wisdom.weizmann.ac.il/~atanay/metac_data/PBMC_processed.tar.gz",
    "raw/PBMC_processed.tar.gz"
)
untar("raw/PBMC_processed.tar.gz", exdir = "raw")
create_project("PBMC")
} # }
```
