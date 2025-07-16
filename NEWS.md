# MCView 0.2.48

* Added controls for plotly download format, width, height and scale.
* Fix: new cell types loaded from annotation file were not updated in all tabs.

# MCView 0.2.47

* Fix: diff. expr. crashed when comparing metacells.

# MCView 0.2.46

* Added "Move up" and "Move down" buttons to the 'Cell types' tab.
* Added "Random seed" input to the 'Manifold' tab.

# MCView 0.2.45

* Reorder boxplot in cell types tab by cell type colors table. 
* Added "Show correlations" button, and compute correlation based on selected cell types in the 'Cell types' tab.

# MCView 0.2.44

* Added "Use DE genes" button to the 'Markers' tab, that uses the significantly differentially expressed genes from the 'Diff. Expr.' tab.

# MCView 0.2.43

* Added "Highlight selected genes" option to the 'Markers' tab.
* Added help for the 'Markers' tab.
* Fix: changing `force_cell_type` did not update the heatmap when in zoom mode.
* Fix: show "Enrichment type" toggle when filtering by clipboard, and allow updating genes.
* Added profiling mode for debug purposes.

# MCView 0.2.42 

* Fix: apps did not start when `light_version` was TRUE.

# MCView 0.2.41

* Projects are now deployment-ready by default. An app.R file is automatically added during dataset import.
* Unified the project and bundle structures - no need for a separate bundle creation step for basic deployment.
* Added "Download data" button.
* Keep the creation command in the config directory.
* Added `copy_source_file` parameter to import_dataset to optionally copy the source h5ad file instead of creating a symlink.
* Improved project directory handling with streamlined overwrite logic.
* Added missing dataset name detection in import_dataset.

# MCView 0.2.40 

* Show metadata in markers heatmap tooltip.
* Added search to cell type selectors.
* Reset on apply by default in `Annotate` tab.

# MCView 0.2.39

* Added `log2` option to the scatter plots.
* Fix: 'Stroke' was not working in the 2D projection plots.

# MCView 0.2.38

* Show GeneCards buttons also when Show-Correlations is off.

# MCView 0.2.36

* Fix: Metadata/Gene correlations were not shown.

# MCView 0.2.35

* Changed default of 'Show correlations' to FALSE.

# MCView 0.2.34

* Added an option to view corrected gene expression in the 'Genes' tab (discussion #170).
* Added cell type selectors to the 'Annotate' tab.
* Turn off the correlations on the 'Genes' tab by default when selecting a subset of the cell types. 

# MCView 0.2.33

* Added cell type selector to the 'Genes' tab.
* Added correlation between the genes/metadata shown in the 'Genes' tab.
* Fix: rescale metadata colors on zoom in 'Markers' tab.
* Fix: Adding new cell type prevented the use of diff. expression tab.
* Disabled 'select all' in 'Markers' gene selector.

# MCView 0.2.32

* Search genes by `startsWith` instead of `contains`.
* Slightly larger default point size in gene/gene plots.
* Do not show cell types with no metacells in the 'Markers' tab (issue #161).

# MCView 0.2.31 

* Do not reorder tabs when they are set in the config file.

# MCView 0.2.30

* Added an option to include metadata when downloading marker heatmap.

# MCView 0.2.29 

* Added "Enrichment type" toggle to the 'Markers' tab.

# MCView 0.2.28

* Added an option to show only fitted genes in the projected fold tab.
* Moved "Order by" and "Force cell type" to left sidebar in the 'Markers' tab.
* Added "Order cell types by" option to the 'Markers' tab.
* Fix: legend did not appear in the 'Markers' tab when using `ggplot 3.5.0`.

# MCView 0.2.27

* Added an option to order by metadata when 'force cell type' is FALSE.
* Added lateral/noisy labels to gene selectors.

# MCView 0.2.26 

* Multiple fixes:
- Fix: loading of light_version bundles was broken.
- Fix: `projection_point_size` and `scatters_point_size` were not detected when added in `other_params`. 
- Fix: import warned when there were metacells withour an edge when providing `default_graph` parameter.

* `About` tab is now the first in light-version bundles. 

# MCView 0.2.25 

* Added an option to exclude lateral and noisy genes from marker heatmap.

# MCView 0.2.24

* Added 'Load genes' and 'Save genes' buttons to the 'Markers' tab.
* Added 'Top correlated' widgets to metadata and gene modules.
* Added some missing package imports.

# MCView 0.2.23

* Fix: plot limits were sometimes too low in differential expression plots.

# MCView 0.2.22

* Added an option to sort gene heatmaps by metadata within each cell type. 
* In 'Projection QC' tab, changed percentage of 'similar' to percentage of 'dissimilar' metacells.

# MCView 0.2.21

* Fixed gene name colors in the 'Gene modules' tab. 

# MCView 0.2.20 

* Fixed metadata colors at 2D projections (previously the breaks were simply ignored)

# MCView 0.2.19

* Added an option to use disk for shiny caching (`see shiny_cache_dir` and `shiny_cache_max_size` parameters in `create_bundle`)
* Fix: sidebar still showed spaces for excluded tabs.

# MCView 0.2.18

* Added 'overwrite' parameter to `import_dataset`.
* Fixed an edge case in which a cell type was not continuous in the markers heatmap order.

# MCView 0.2.17

* Fix: atlas tab was not working after the new UMAP feature.
* Use corrected expression in the differential expression in the Query tab. 
* Mark dissimilar metacells in metacell selectors.

# MCView 0.2.16

* Added gene colors legend to the 'Markers' tab.

# MCView 0.2.15

* Added 'any' cell type to gene metadata table at Query tab.
* Added a value box with % of fitted genes out of atlas markers.

# MCView 0.2.14

* Query tab improvements: added a correction factor line and an option to show only fitted genes.

# MCView 0.2.13 

* Added UMAP computation given anchor genes to the 'Manifold' tab.

# MCView 0.2.12

* Added `Sample types` plot to the `Samples` tab.
* Fixed a few bugs in the `Samples` tab.

# MCView 0.2.11

* Fix: fitted gene plot was not shown when 'corrected_factor' was missing from the input data.

# MCView 0.2.10

* Changed zero-fold per gene plot to show expression vs zero-fold. 
* Added two static plots to About tab.
* Added "top-correlated" panel to "Cell types" tab.
* Added "inner-std" plot to QC panel. 
* Added "Stdev-fold" tab.
* Added "Projection QC" tab
* Fix: differential expression between metacells did not work with most recent version of `qs`.
* Added Lasso selection in Annotation and Differential Expression tabs.
* Added "# of metacells with significant inner-fold" to the "QC" tab.
* Filter lateral/noisy genes from diff expression plots.
* Annotate genes with their lateral/noisy status + gene modules.

# MCView 0.2.9

* Fixed import according to newest metacells version.

# MCView 0.2.8

* Added the 'QC' tab. Please re-import your dataset in order to enable it.
* Moved some of 'gene/gene' and 'projection' controls to below the plot. 
* Added 'excluded_tabs' to the config file.
* Added a footer with MCView and metacells versions.
* Some modifications to the 'Genes' tab layout. 
* Added an option to color the 2D projection by an axis from the Gene/Gene plot. 

# MCView 0.2.7

* Added QC metrics to the manifold tab.
* Fix: import_cell_metadata was overwriting existing metadata
* Fix: import crashed when metacell_types had missing colors

# MCView 0.2.6

* Use new `metacells` version. Note that old `anndata` objects will not work with this version without conversion.
* Added an option to bundle a 'light' version without the option to change the genes on heatmaps.
* Fixed resolution of heatmap plots.
* Added a 'rename' modal to the 'Annotate' tab.
* Fixed order of negative correlations in the left panels. Fixes #127.
* Added an option to use 'gene_names' parameter with atlas projection. 
* Use the name of the dataset as title by default.
* Rasterized Type Prediction plot.

# MCView 0.2.5

* Added tab selector. 
* Runtime optimizations: load tabs only on click.
* Some css fixes.
* Changed default 'about' document.
* Added preloader spinner.
* Allow using `run_app` with the bundle folder.
* Added cell metadata support for metacell1 import. 
* Fix: clicking on a cell type on the manifold legend was sometimes not aligned with the gene/gene plot. 
* Fix: categorical metadata tooltip was wrong in projection plots.
* Changed metacell selectors to search with "contains" instead of "startWith" behaviour. 

# MCView 0.2.4

* Add an option to add metacells to clipboard by number 
* runtime optimizations: use virtualSelectInput 

# MCView 0.2.3

* Added support for gene modules ('Gene modules' tab). 
* Added 'Outliers' tab.
* Added Clipboard functionality.
* Added support for alternative gene names. 
* Added control over the number of clusters / whether to cluster at all in the import phase.
* Added new selection controls to Annotate tab ('Gears' symbol on the upper right of 'Metacell Annotation' box).
* Bug fix: heatmaps crashed when only a single gene was chosen.  
* Bug fix: automatically color cell types missing from cell type colors file.
* Bug fix: Diff expr. crashed when there was only a single cell type.
* Bug fix: Wrong tooltip in metadata/metadata scatter plot. 
* Bug fix: Annotation colors in markers heatmap were scaled relative to current zoom.
* Added control for marker selection parameters in `import_dataset`.
# MCView 0.2.2

* Added support for atlas projections ('Atlas', 'Query' and 'Projected-fold' tabs)
* Added tooltip and interactive clicks to heatmaps ('Markers' and 'Projected-fold' tabs)
* Added zoom for heatmaps. 
* Homogenize Annotate, Diff. Expr and Genes tabs. 
* Added 'Cell types' tab. 
* Added support for cell metadata using the 'Samples' tab. 
* Allow coloring by categorical metadata annotations. 
* Moved flow graphs to their own tab named 'Flow'. 
* `create_project` and `import_dataset` now accept the config file parameters. 
* Added option to force x and y limits to be the same in gene/gene plots.
* Added an option to show x=y line in gene/gene plots.

# MCView 0.2.1

* Added Metadata tab.
* Added Markers tab.
* Added inner-fold matrix to markers tab. 
* Added an option to compare cell types on Metacell tab. 
* Added a download button.
* Added a busy spinner. 
* `create_project` now only takes a path (instead of path + project name)
* implicitly create a project when importing
* import cell type annotation using a single file that contains both metacell types and cell type colors. 
* Bug fixes: issues #48, #51, #52, #60, #62 and more. 
* Changed gene selectors on the manifold tab to be on-demand (to reduce initialization time).
* Added (some) caching.

# MCView 0.2.1

* A few bug fixes. See issues #36, #37, #38, #40, #41 and #42. 

# MCView 0.2.0

* First stable version.


