# Update metadata for a dataset

Update metadata for a dataset

## Usage

``` r
update_metadata(
  project,
  dataset,
  metadata = NULL,
  metadata_fields = NULL,
  anndata_file = NULL,
  overwrite = FALSE
)
```

## Arguments

- project:

  path to the project

- dataset:

  name for the dataset, e.g. "PBMC". The name of the dataset can only
  contain alphanumeric characters, dots, dashes and underscores.

- metadata:

  can be either a data frame with a column named "metacell" with the
  metacell id and other metadata columns or a name of a delimited file
  which contains such data frame. See `metadata_fields` for other
  details.

- metadata_fields:

  names of fields in the anndata `object$obs` which contains metadata
  for each metacell.  
  The fields should can be either numeric or categorical.  
  You can use `cell_metadata_to_metacell` to convert from categorical to
  a numeric score (e.g. by using fraction of the category). You can use
  'all' in order to import all the fields of the anndata object.

- anndata_file:

  path to `h5ad` file which contains the output of metacell2 pipeline
  (metacells python package) or a loaded anndata object of the same
  format.

- overwrite:

  overwrite all existing metadata. If `FALSE` - would override only
  existing metadata fields.
