# Update cell type assignment for each metacell

Change the cell type assignments for each metacell to the ones listed at
`metacell_types_file`.

## Usage

``` r
update_metacell_types(project, dataset, metacell_types_file)
```

## Arguments

- project:

  path to the project directory

- dataset:

  name for the dataset, e.g. "PBMC"

- metacell_types_file:

  path to a tabular file (csv,tsv) with cell type assignement for each
  metacell. The file should have a column named "metacell" with the
  metacell ids and another column named "cell_type" or "cluster" with
  the cell type assignment. Metacell ids that do not exists in the data
  would be ignored.

## Details

This is usually done after a first iteration of annotation using the
"Annotate" tab in the MCView annotation, which can export a valid
`metacell_types_file`. The file should have a column named "metacell"
with the metacell ids and another column named "cell_type" or "cluster"
with the cell type assignment.

Note that the exported file from the \_\_MCView\_\_ app contains
additional fields which will be ignored in this function.

Under the hood - MCView updates a file named "metacell_types.tsv" under
`project/cache/dataset`, which can also be edited manually.

If the file contains an additional 'color' field, the cell type colors
would be updated as well.

## Examples

``` r
if (FALSE) { # \dontrun{
update_metacell_types("PBMC", "PBMC163k", "raw/metacell-clusters.csv")
} # }
```
