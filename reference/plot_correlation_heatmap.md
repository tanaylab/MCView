# Create correlation heatmap for multiple genes

Create correlation heatmap for multiple genes

## Usage

``` r
plot_correlation_heatmap(
  correlation_results,
  input_genes,
  dataset,
  cluster = TRUE,
  max_genes = 50,
  mask_low_correlations = FALSE,
  correlation_mode = "individual"
)
```

## Arguments

- correlation_results:

  Data frame from calc_individual_correlations

- input_genes:

  Vector of input gene names

- dataset:

  Dataset name

- cluster:

  Whether to cluster genes

- max_genes:

  Maximum number of genes to include in heatmap
