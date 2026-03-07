# Update default 2D projection for a dataset

Update default 2D projection for a dataset

## Usage

``` r
update_2d_projection(project, dataset, layout, graph)
```

## Arguments

- project:

  path to the project

- dataset:

  name for the dataset, e.g. "PBMC". The name of the dataset can only
  contain alphanumeric characters, dots, dashes and underscores.

- layout:

  a data frame with a column named "metacell" with the metacell id and
  other columns with the x and y coordinates of the metacell.

- graph:

  a data frame with a column named "from", "to" and "weight" with the
  ids of the metacells and the weight of the edge. If NULL, the graph
  would be taken from the anndata object.
