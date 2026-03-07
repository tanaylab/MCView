# Calculate correlations for multiple genes individually

Calculate correlations for multiple genes individually

## Usage

``` r
calc_individual_correlations(
  dataset,
  genes,
  n_top = 30,
  cell_type_filter = NULL,
  threshold = 0,
  atlas = FALSE
)
```

## Arguments

- dataset:

  Dataset name

- genes:

  Vector of gene names

- n_top:

  Number of top correlations to return per gene

- cell_type_filter:

  Optional vector of cell types to filter metacells

- threshold:

  Minimum correlation threshold to include

## Value

Data frame with correlation results
