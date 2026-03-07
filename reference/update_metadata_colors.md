# Update metadata colors for a dataset

Update metadata colors for a dataset

## Usage

``` r
update_metadata_colors(project, dataset, metadata_colors, overwrite = FALSE)
```

## Arguments

- project:

  path to the project

- dataset:

  name for the dataset, e.g. "PBMC". The name of the dataset can only
  contain alphanumeric characters, dots, dashes and underscores.

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

- overwrite:

  overwrite all existing colors. If `FALSE` - would override only the
  colors of existing metadata fields.
