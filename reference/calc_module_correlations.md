# Calculate correlations for genes as a module/anchor

Calculate correlations for genes as a module/anchor

## Usage

``` r
calc_module_correlations(
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

  Vector of gene names to treat as module

- n_top:

  Number of top correlations to return

- cell_type_filter:

  Optional vector of cell types to filter metacells

- threshold:

  Minimum correlation threshold to include

## Value

Data frame with correlation results
