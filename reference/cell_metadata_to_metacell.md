# Convert cell metadata to metacell metadata

Summarise cell metadata for each metacell. Metadata fields can be either
numeric and then the summary function `func` is applied for the values
of each field, or categorical metadata fields which are expanded to
multiple metadata columns with the fraction of cells (in each metacell)
for every category. Variables that are either character, factor or are
explicitly set at `categorical` are treated as categorical.

`cell_metadata_to_metacell` converts cell metadata to metacell metadata
from data frames.  
`cell_metadata_to_metacell_from_h5ad` extracts metadata fields and
cell_to_metacell from cells h5ad file and then runs
`cell_metadata_to_metacell`.  
`cell_metadata_to_metacell_from_metacell1` extracts metadata fields and
cell_to_metacell from metacell1 scdb and then runs
`cell_metadata_to_metacell`.

## Usage

``` r
cell_metadata_to_metacell(
  cell_metadata,
  cell_to_metacell,
  func = mean,
  categorical = c()
)

cell_metadata_to_metacell_from_metacell1(
  scdb,
  matrix,
  mc,
  metadata_fields,
  func = mean,
  categorical = c()
)

cell_metadata_to_metacell_from_h5ad(
  anndata_file,
  metadata_fields,
  func = mean,
  categorical = c(),
  rm_outliers = TRUE
)
```

## Arguments

- cell_metadata:

  data frame with a column named "cell_id" with the cell id and other
  metadata columns, or a name of a delimited file which contains such
  data frame.

- cell_to_metacell:

  data frame with a column named "cell_id" with cell id and another
  column named "metacell" with the metacell the cell is part of, or a
  name of a delimited file which contains such data frame.

- func:

  summary function for the cell metadata for non categorical metadata
  columns (e.g. mean, median, sum)

- categorical:

  a vector with names of categorical variables. The returned data frame
  would have a column for each category where the values are the
  fraction of cells with the category in each metacell.

- scdb, matrix, mc:

  scdb, matrix and mc objects from metacell1. See
  `import_dataset_metacell1` for more information.

- metadata_fields:

  names of fields in the anndata `object$obs` which contains metadata
  for each cell.

- anndata_file:

  path to `h5ad` file which contains the output of metacell2 pipeline
  (metacells python package).

- rm_outliers:

  do not calculate statistics for cells that are marked as outliers
  (`outiler=TRUE` in `object$obs`) (only relevant when running
  `cell_metadata_to_metacell_from_h5ad`)

## Value

A data frame with a column named "metacell" and the metadata columns
from `cell_metadata` summarized for each metacell using `func` for
non-categorical variables, and a column for each category of the
categorical metadata variables in`cell_metadata`, where the values are
the fraction of cells with the category in each metacell.

## Functions

- `cell_metadata_to_metacell_from_metacell1()`:

- `cell_metadata_to_metacell_from_h5ad()`:

## Examples

``` r
set.seed(60427)
n_cells <- 5e6
cell_metadata <- tibble::tibble(
    cell_id = 1:n_cells,
    md1 = sample(1:5, size = n_cells, replace = TRUE),
    md2 = rnorm(n = n_cells),
    md_categorical1 = sample(paste0("batch", 1:5), size = n_cells, replace = TRUE),
    md_categorical2 = sample(1:5, size = n_cells, replace = TRUE)
)

cell_to_metacell <- tibble::tibble(
    cell_id = 1:n_cells,
    metacell = sample(0:1535, size = n_cells, replace = TRUE)
)
metadata <- cell_metadata_to_metacell(
    cell_metadata[, 1:3],
    cell_to_metacell
)
#> ℹ Numerical variables: md1, md2
head(metadata)
#> # A tibble: 6 × 3
#>   metacell   md1      md2
#>      <int> <dbl>    <dbl>
#> 1        0  2.98 -0.0219 
#> 2        1  2.96  0.0152 
#> 3        2  2.99  0.00563
#> 4        3  3.01 -0.0268 
#> 5        4  3.01  0.0116 
#> 6        5  2.96  0.00804

metadata1 <- cell_metadata_to_metacell(
    cell_metadata[, 1:3], cell_to_metacell,
    func = function(x) x * 2
)
#> ℹ Numerical variables: md1, md2
#> Error in summarise(.tbl, !!!funs): ℹ In argument: `md1 = (function (x) ...`.
#> ℹ In group 1: `metacell = 0`.
#> Caused by error:
#> ! `md1` must be size 1, not 3260.
#> ℹ To return more or less than 1 row per group, use `reframe()`.
head(metadata1)
#> Error: object 'metadata1' not found


metadata3 <- cell_metadata_to_metacell(
    cell_metadata,
    cell_to_metacell,
    categorical = c("md_categorical1", "md_categorical2")
)
#> ℹ Categorical variables: md_categorical1, md_categorical2
#> ℹ Numerical variables: md1, md2

head(metadata3)
#> # A tibble: 6 × 13
#>   metacell md_categorical1: batc…¹ md_categorical1: bat…² md_categorical1: bat…³
#>      <int>                   <dbl>                  <dbl>                  <dbl>
#> 1        0                   0.195                  0.204                  0.205
#> 2        1                   0.188                  0.202                  0.210
#> 3        2                   0.221                  0.199                  0.193
#> 4        3                   0.192                  0.203                  0.213
#> 5        4                   0.188                  0.205                  0.205
#> 6        5                   0.206                  0.193                  0.203
#> # ℹ abbreviated names: ¹​`md_categorical1: batch1`, ²​`md_categorical1: batch2`,
#> #   ³​`md_categorical1: batch3`
#> # ℹ 9 more variables: `md_categorical1: batch4` <dbl>,
#> #   `md_categorical1: batch5` <dbl>, `md_categorical2: 1` <dbl>,
#> #   `md_categorical2: 2` <dbl>, `md_categorical2: 3` <dbl>,
#> #   `md_categorical2: 4` <dbl>, `md_categorical2: 5` <dbl>, md1 <dbl>,
#> #   md2 <dbl>
if (FALSE) { # \dontrun{
cell_metadata_to_metacell_from_h5ad("cells.h5ad", c("pile", "age", "batch"), categorical = "batch")
} # }
```
