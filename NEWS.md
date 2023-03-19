# MCView 0.2.9.9000

* Changed zero-fold per gene plot to show expression vs zero-fold. 
* Added two static plots to About tab.
* Added "top-correlated" panel to "Cell types" tab.
* Added "inner-std" plot to QC panel. 
* Added "Stdev-fold" tab.

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


