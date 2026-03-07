# Import cell metadata to an MCView dataset

Import metadata which is at the cell level to MCView. The metadata can
be summarised to the metacell level by setting `summarise_md` to TRUE,
in which case it could be shown at the "Genes" and "Markers" tabs. In
order to view data at the samples level, an additional sample identifier
should be given as a column named "samp_id" in the `cell_metadata` data
frame.

## Usage

``` r
import_cell_metadata(
  project,
  dataset,
  cell_metadata,
  cell_to_metacell = NULL,
  summarise_md = FALSE,
  add_samples_tab = TRUE,
  ...
)
```

## Arguments

- project:

  path to the project

- dataset:

  name for the dataset, e.g. "PBMC". The name of the dataset can only
  contain alphanumeric characters, dots, dashes and underscores.

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

- summarise_md:

  summarise cell metadata to the metacell level.

- add_samples_tab:

  add the 'Samples' tab to the config file if it doesn't exist

- ...:

  Arguments passed on to
  [`cell_metadata_to_metacell`](cell_metadata_to_metacell.md)

  `func`

  :   summary function for the cell metadata for non categorical
      metadata columns (e.g. mean, median, sum)

  `categorical`

  :   a vector with names of categorical variables. The returned data
      frame would have a column for each category where the values are
      the fraction of cells with the category in each metacell.
