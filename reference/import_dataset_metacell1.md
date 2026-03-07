# Import a dataset to an MCView project from metacell R package

Read objects from `metacell` R package and import a metacell dataset to
MCView.

## Usage

``` r
import_dataset_metacell1(
  project,
  dataset,
  scdb,
  matrix,
  mc,
  mc2d,
  metacell_types_file,
  cell_type_colors_file,
  gene_modules_file = NULL,
  gene_modules_k = NULL,
  calc_gg_cor = TRUE,
  network = NULL,
  time_annotation_file = NULL,
  time_bin_field = NULL,
  metadata_fields = NULL,
  categorical = c(),
  ...
)
```

## Arguments

- project:

  path to the project

- dataset:

  name for the dataset, e.g. "PBMC". The name of the dataset can only
  contain alphanumeric characters, dots, dashes and underscores.

- scdb:

  path to R metacell single cell RNA database

- matrix:

  name of the umi matrix to use

- mc:

  name of the metacell object to use

- mc2d:

  name of the 2d projection object to use

- metacell_types_file:

  path to a tabular file (csv,tsv) with cell type assignement for each
  metacell. The file should have a column named "metacell" with the
  metacell ids and another column named "cell_type" or "cluster" with
  the cell type assignment. Metacell ids that do not exists in the data
  would be ignored. In addition, the file can have a column named "age"
  or "mc_age" with age metadata per metacell

- cell_type_colors_file:

  path to a tabular file (csv,tsv) with color assignement for each cell
  type. The file should have a column named "cell_type" or "cluster"
  with the cell types and another column named "color" with the color
  assignment. Cell types that do not exist in the metacell types would
  be ignored.

- gene_modules_file:

  path to a tabular file (csv,tsv) with assignment of genes to gene
  modules. Should have a field named "gene" with the gene name and a
  field named "module" with the name of the gene module.

- gene_modules_k:

  number of clusters for initial gene module calculation. If NULL - the
  number of clusters would be determined such that an gene module would
  contain 16 genes on average.

- calc_gg_cor:

  Calculate top 30 correlated and anti-correlated genes for each gene.
  This computation can be heavy for large datasets or weaker machines,
  so you can set `calc_gg_cor=FALSE` to skip it. Note that then this
  feature would be missing from the app.

- network:

  name of the network object to use (optional)

- time_annotation_file:

  file with names for time bins (optional, only relevant with
  networks/flows). Should have a field named "time_bin" with the time
  bin id and another field named "time_desc" which contains the
  description of the time bin

- time_bin_field:

  name of a field in `cell_metadata` which contains time bin per cell
  (optional)

- metadata_fields:

  names of fields `mat@cell_metadata` which contains metadata per cell
  to be summarized using `cell_metadata_to_metacell`.  
  The fields should can be either numeric or categorical.  
  You can use `cell_metadata_to_metacell` to convert from categorical to
  a numeric score (e.g. by using fraction of the category).

- categorical:

  metadata fields that should be treated as categorical (optional)

- ...:

  Arguments passed on to [`create_project`](create_project.md)

  `title`

  :   The title of the app. This would be shown on the top left of the
      screen.

  `tabs`

  :   Controls which tabs to show in the left sidebar and their order.
      Options are: "QC", "Projection-QC", "Manifold", "Genes", "Query",
      "Atlas", "Markers", "Gene modules", "Projected-fold", "Diff.
      Expression", "Cell types", "Flow", "Annotate", "About". When
      NULL - default tabs would be set. For projects with atlas
      projections, please set `atlas` to TRUE.

  `selected_gene1,selected_gene2`

  :   The default genes that would be selected (in any screen with gene
      selection). If this parameter is missing, the 2 genes with highest
      max(expr)-min(expr) in the first dataset would be chosen.

  `selected_mc1,selected_mc2`

  :   The default metacells that would be selected in the Diff.
      Expression tab.

  `datasets`

  :   A named list with additional per-dataset parameters. Current
      parameters include default visualization properties of projection
      and scatter plots.

  `other_params`

  :   Named list of additional parameters such as projection_point_size,
      projection_point_stroke, scatters_point_size and
      scatters_stroke_size

  `edit_config`

  :   open file editor for config file editing

  `atlas`

  :   use default configuration for atlas projections (relevant only
      when `tabs` is NULL)

## Details

The result would be a directory under `project/cache/dataset` which
would contain objects used by MCView shiny app (such as the metacell
matrix).

In addition, you can supply file with type assignment for each metacell
(`metacell_types_file`) and a file with color assignment for each
metacell type (`cell_type_colors_file`).

Make sure that you have the R `metacell` package installed in order to
use this function.

`network`, `time_annotation_file` and `time_bin_field` are only relevant
if you computed flows/networks for your dataset and therefore are
optional.

In order to add time annotation to your dataset you will have to:

- 1\. Add a column named "mc_age" or "age" to `metacell_types_file` with
  time per metacell

- 2\. Create a `time_annotation_file` with id for each time bin and
  description

## Examples

``` r
if (FALSE) { # \dontrun{
import_dataset_metacell1(
    "embflow",
    "153embs",
    scdb = "raw/scrna_db",
    matrix = "embs",
    mc = "embs",
    mc2d = "embs",
    metacell_types_file = "raw/metacell-types.csv",
    cell_type_colors_file = "raw/cell-type-colors.csv",
    network = "embs",
    time_annotation_file = "raw/time-annot.tsv",
    time_bin_field = "age_group"
)
} # }
```
