# Update gene modules

Change the gene modules of an MCView app.

## Usage

``` r
update_gene_modules(project, dataset, gene_modules_file)
```

## Arguments

- project:

  path to the project

- dataset:

  name for the dataset, e.g. "PBMC". The name of the dataset can only
  contain alphanumeric characters, dots, dashes and underscores.

- gene_modules_file:

  path to a tabular file (csv,tsv) with assignment of genes to gene
  modules. Should have a field named "gene" with the gene name and a
  field named "module" with the name of the gene module.

## Details

This is usually done after a first iteration of annotation using the
"Gene modules" tab in the MCView annotation, which can export a valid
`gene_modules_file`. The file should have a column named "gene" with the
gene names and another column named "module" with the id of the gene
module. Note that the exported file from the \_\_MCView\_\_ app might
contain additional fields which will be ignored in this function.

Under the hood - MCView updates a file named "gene_modules.tsv" under
`project/cache/dataset`, which can also be edited manually.

## Examples

``` r
if (FALSE) { # \dontrun{
update_gene_modules("PBMC", "PBMC163k", "raw/gene-modules.csv")
} # }
```
