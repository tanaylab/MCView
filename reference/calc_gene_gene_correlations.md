# Calculate gene-gene correlations (correlation between input genes)

Calculate gene-gene correlations (correlation between input genes)

## Usage

``` r
calc_gene_gene_correlations(
  dataset,
  genes,
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

- cell_type_filter:

  Optional vector of cell types to filter metacells

- threshold:

  Minimum correlation threshold to include

## Value

Data frame with gene-gene correlation results
