# Import a dataset to an MCView project

Read an `anndata` file which is the output of python `metacells`
package, and import the metacell dataset to MCView. Each project can
have multiple datasets which can be in the app using the right sidebar.

## Usage

``` r
import_dataset(
  project,
  dataset,
  anndata_file,
  cell_type_field = NULL,
  metacell_types_file = NULL,
  cell_type_colors_file = NULL,
  outliers_anndata_file = NULL,
  cluster_metacells = TRUE,
  cluster_k = NULL,
  metadata_fields = NULL,
  metadata = NULL,
  metadata_colors = NULL,
  cell_metadata = NULL,
  cell_to_metacell = NULL,
  gene_modules_file = NULL,
  gene_modules_k = NULL,
  calc_gg_cor = TRUE,
  gene_names = NULL,
  metacell_graphs = NULL,
  atlas_project = NULL,
  atlas_dataset = NULL,
  projection_weights_file = NULL,
  copy_atlas = TRUE,
  minimal_max_log_fraction = -13,
  minimal_relative_log_fraction = 2,
  umap_anchors = NULL,
  umap_config = NULL,
  min_umap_log_expr = -14,
  genes_per_anchor = 30,
  layout = NULL,
  default_graph = NULL,
  overwrite = TRUE,
  copy_source_file = FALSE,
  ...
)
```

## Arguments

- project:

  path to the project

- dataset:

  name for the dataset, e.g. "PBMC". The name of the dataset can only
  contain alphanumeric characters, dots, dashes and underscores.

- anndata_file:

  path to `h5ad` file which contains the output of metacell2 pipeline
  (metacells python package) or a loaded anndata object of the same
  format.

- cell_type_field:

  name of a field in the anndata `object$obs` which contains a cell type
  (optional). If the field doesn't exist and `metacell_types_file` is
  missing, MCView would first look for a field named 'projected_type',
  'type', 'cell_type' or 'cluster' at `object$obs` (in this order), and
  if it doesn't exists MCView would cluster the metacell matrix using
  kmeans++ algorithm (from the `tglkmeans` package).

- metacell_types_file:

  path to a tabular file (csv,tsv) with cell type assignement for each
  metacell. The file should have a column named "metacell" with the
  metacell ids and another column named "cell_type", or "cluster" with
  the cell type assignment. Metacell ids that do not exists in the data
  would be ignored.  
  If this parameter and `cell_type_field` are missing and
  `cluster_metacells=TRUE`, MCView would cluster the metacell matrix
  using kmeans++ algorithm (from the `tglkmeans` package).  
  If the file has a field named 'color' and
  `cell_type_colors_file=NULL`, the cell types colors would be used.

- cell_type_colors_file:

  path to a tabular file (csv,tsv) with color assignement for each cell
  type. The file should have a column named "cell_type" or "cluster"
  with the cell types and another column named "color" with the color
  assignment.  
  In case the metacell types (given by a file or from the
  `anndata_file`) contain types which are not present in the cell type
  colors file, MCView would use the `chameleon` package to assign a
  color for them in addition to the cell type colors.  
  If this is missing, and `metacell_types_file` did not have a 'color'
  field, MCView would use the `chameleon` package to assign a color for
  each cell type.  
  When an atlas is given (using `atlas_project` and `atlas_dataset`), if
  the cell types are the same as the atlas, the atlas colors would be
  used.

- outliers_anndata_file:

  path to anndata file with outliers (optional). This would enable, by
  default, the following tabs: \["Outliers", "Similar-fold",
  "Deviant-fold"\]. See the metacells python package for more details.

- cluster_metacells:

  When TRUE and no metacell type is given (via `metacell_types_file` or
  `cell_type_field` - implicit and explicit), MCView would cluster the
  metacell matrix using kmeans++ algorithm (from the `tglkmeans`
  package).

- cluster_k:

  number of clusters for initial metacell clustering. If NULL - the
  number of clusters would be determined such that a metacell would
  contain 16 cells on average.

- metadata_fields:

  names of fields in the anndata `object$obs` which contains metadata
  for each metacell.  
  The fields should can be either numeric or categorical.  
  You can use `cell_metadata_to_metacell` to convert from categorical to
  a numeric score (e.g. by using fraction of the category). You can use
  'all' in order to import all the fields of the anndata object.

- metadata:

  can be either a data frame with a column named "metacell" with the
  metacell id and other metadata columns or a name of a delimited file
  which contains such data frame. See `metadata_fields` for other
  details.

- metadata_colors:

  a named list with colors for each metadata column, or a name of a yaml
  file with such list. For numerical metadata columns, colors should be
  given as a list where the first element is a vector of colors and the
  second element is a vector of breaks.  
  Note that in case the breaks are out of the range of the metadata
  values, the colors would be used for the minimum and maximum values.  
  If only colors are given breaks would be implicitly determined from
  the minimum and maximum of the metadata field.  
  For categorical metadata columns, color can be given either as a named
  vector where names are the categories and the values are the colors,
  or as a named list where the first element named 'colors' holds the
  colors, and the second element called 'categories' holds the
  categories.

- cell_metadata:

  data frame with a column named "cell_id" with the cell id and other
  metadata columns, or a name of a delimited file which contains such
  data frame. For activating the "Samples" tab, the data frame should
  have an additional column named "samp_id" with a sample identifier per
  cell (e.g., batch id, patient etc.). Optionally, a column named
  "metacell" can be added to the data frame, which will be used instead
  of the `cell_to_metacell` parameter.

- cell_to_metacell:

  data frame with a column named "cell_id" with cell id and another
  column named "metacell" with the metacell the cell is part of, or a
  name of a delimited file which contains such data frame. If NULL, the
  metacell will be inferred from the 'metacell' column in
  `cell_metadata`.

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

- gene_names:

  use alternative gene names (optional). A data frame with a column
  called 'gene_name' with the original gene name (as it appears at the
  'h5ad' file) and another column called 'alt_name' with the gene name
  to use in MCView. Genes that do not appear at the table would not be
  changed.

- metacell_graphs:

  a named list of metacell graphs or files containing metacell graphs.
  Each graph should be a data frame columns named "from", "to" and
  "weight" with the ids of the metacells and the weight of the edge. If
  the list is not named, the names would be 'graph1', 'graph2' and so
  on. Note that the graph cannot be named "metacell" as this is reserved
  for the metacell graph.

- atlas_project:

  path to and `MCView` project which contains the atlas.

- atlas_dataset:

  name of the atlas dataset

- projection_weights_file:

  Path to a tabular file (csv,tsv) with the following fields "query",
  "atlas" and "weight". The file is an output of `metacells` projection
  algorithm.

- copy_atlas:

  copy atlas MCView to the current project. If FALSE - a symbolic link
  would be created instead.

- minimal_max_log_fraction:

  When choosing marker genes: take only genes with at least one value
  (in log fraction units - normalized egc) above this threshold

- minimal_relative_log_fraction:

  When choosing marker genes: take only genes with relative log fraction
  (mc_fp) above this this value

- umap_anchors:

  a vector of gene names to use for UMAP calculation. If NULL, the umap
  from the anndata object would be used.

- umap_config:

  a named list with UMAP configuration. See
  [`umap::umap`](https://rdrr.io/pkg/umap/man/umap.html) for more
  details. When NULL, the default configuration would be used, except
  for: min_dist=0.96, n_neighbors=10, n_epoch=500.

- min_umap_log_expr:

  minimal log2 expression for genes to use for UMAP calculation.

- genes_per_anchor:

  number of genes to use for each umap anchor.

- layout:

  a data frame with a column named "metacell" with the metacell id and
  other columns with the x and y coordinates of the metacell. If NULL,
  the layout would be taken from the anndata object.

- default_graph:

  a data frame with a column named "from", "to" and "weight" with the
  ids of the metacells and the weight of the edge. If NULL, the graph
  would be taken from the anndata object.

- overwrite:

  if a dataset with the same name already exists, overwrite it.
  Otherwise, an error would be thrown.

- copy_source_file:

  if TRUE, copy the source file to the project cache directory. If
  FALSE, create a symbolic link to the source file.

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

## Value

invisibly returns an `AnnDataR6` object of the read `anndata_file`

## Details

The function would create a directory under `project/cache/dataset`
which would contain objects used by MCView shiny app (such as the
metacell matrix).

In addition, you can supply file with type assignment for each metacell
(`metacell_types_file`) and a file with color assignment for each
metacell type (`cell_type_colors_file`).

## Examples

``` r
if (FALSE) { # \dontrun{
dir.create("raw")
download.file(
    "http://www.wisdom.weizmann.ac.il/~atanay/metac_data/PBMC_processed.tar.gz",
    "raw/PBMC_processed.tar.gz"
)
untar("raw/PBMC_processed.tar.gz", exdir = "raw")
import_dataset("PBMC", "PBMC163k", "raw/metacells.h5ad")
} # }
```
