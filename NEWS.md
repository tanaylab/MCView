# MCView 0.2.2

* Added support for atlas projections ('Atlas', 'Query' and 'Projected-fold' tabs)
* Added tooltip and interactive clicks to heatmaps ('Markers' and 'Projected-fold' tabs)
* Added zoom for heatmaps. 
* Homogenzie Annotate, Diff. Expr and Genes tabs. 
* Added 'Cell types' tab. 
* Added support for cell metadata using the 'Samples' tab. 
* Allow coloring by categorical metadata annotations. 
* Moved flow graphs to their own tab named 'Flow'. 
* `create_project` and `import_dataset` now accept the config file parameters. 
* Added option to force x and y limits to be the same in gene/gene plots.

# MCView 0.2.1

* Added Metadata tab.
* Added Markers tab.
* Added inner-fold matrix to markers tab. 
* Added an option to compare cell types on Metacell tab. 
* Added a download button.
* Added a busy spinner. 
* `create_project` now only takes path (instead of path + project name)
* implicitly create a project when importing
* import cell type annotation using a single file which contains both metacell types and cell type colors. 
* Bug fixes: issues #48, #51, #52, #60, #62 and more. 
* Changed gene selectors on manifold tab to be on demand (to reduce initialization time).
* Added (some) caching.

# MCView 0.2.1

* A few bug fixes. See issues #36, #37, #38, #40, #41 and #42. 

# MCView 0.2.0

* First stable version.


