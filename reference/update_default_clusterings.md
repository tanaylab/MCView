# Store the clusterings of the default heatmaps

Cluster the metacells of the heatmap that each heatmap tab (Markers,
Gene modules, Inner-fold, Stdev-fold, Projected-fold and Outliers)
starts with, and store the clusterings in the project, so that the app
doesn't cluster them the first time a tab is opened, which takes minutes
for large datasets. `import_dataset` already does this. Use this
function for projects imported with an older version of MCView, or after
changing the gene modules with `update_gene_modules`.

## Usage

``` r
update_default_clusterings(project, overwrite = FALSE)
```

## Arguments

- project:

  path to the project directory

- overwrite:

  recompute the clusterings of datasets that already have them

## Details

Under the hood - MCView writes a file named "default_clusterings.qs2"
under `project/cache/dataset`, and skips datasets whose directory is a
symlink to outside the project. When that directory isn't writable and
the environment variable `MCVIEW_CLUSTERING_CACHE_DIR` is set, the
clusterings are stored in that directory instead. Apps that set the
variable (e.g. all the apps of a shiny server) also store there every
clustering they compute and look for clusterings there, so each
clustering is computed once per server rather than once per app process.
(`MCVIEW_SHINY_CACHE_DIR` keeps the apps' plots and tables across
restarts in the same way, see [`run_app`](run_app.md).)

## Examples

``` r
if (FALSE) { # \dontrun{
update_default_clusterings("PBMC")
} # }
```
