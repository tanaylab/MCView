# Update color assignment for each cell type

Change the color assignments for each cell type to the ones listed at
`cell_type_colors_file`.

## Usage

``` r
update_cell_type_colors(project, dataset, cell_type_colors_file)
```

## Arguments

- project:

  path to the project directory

- dataset:

  name for the dataset, e.g. "PBMC"

- cell_type_colors_file:

  path to a tabular file (csv,tsv) with color assignement for each cell
  type. The file should have a column named "cell_type" or "cluster"
  with the cell types and another column named "color" with the color
  assignment. Cell types that do not exist in the metacell types would
  be ignored, so if you changed the names of cell types you would have
  to also update the metacell types (using `update_metacell_types`). The
  function also accepts output of the 'export' button from the
  application annotation page. If this parameter is missing, MCView
  would use the `chameleon` package to assign a color for each cell
  type.

## Details

This is usually done after a first iteration of annotation using the
"Annotate" tab in the MCView annotation, which can export a valid
`cell_type_colors_file`.

The file should have a column named "cell_type" or "cluster" with the
cell types and another column named "color" with the color assignment.
Note that the exported file from the \_\_MCView\_\_ app contains
additional fields which will be ignored in this function.

Under the hood - MCView updates a file named "cell_type_colors.tsv"
under `project/cache/dataset`, which can also be edited manually.

## Examples

``` r
if (FALSE) { # \dontrun{
update_metacell_types("PBMC", "PBMC163k", "raw/cluster-colors.csv")
} # }
```
