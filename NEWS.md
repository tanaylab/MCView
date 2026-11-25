# MCView 0.2.65

Fixes:

* Annotate: "Apply" and "Create New Cell Type" could change other metacells than those selected in the table; a rejected selection no longer creates an empty cell type.
* Metacells without a cell type were hidden from the projection and the scatter of the Genes, Annotate and Gene modules tabs.
* Heatmaps: stayed blank after a cell type was renamed, merged or deleted; without "Force cell type", metacells with few non-zero genes were dropped and the others shifted; the legend of the cell types was cut and the genes legend was not shown; the color bar now says what the values are.
* Plots came back blank after visiting several tabs (the browser ran out of WebGL contexts).
* Clipboard: "Remove selected" removed another metacell.
* Atlas: coloring the projection by a metadata field of the atlas failed when the query doesn't have that field.
* Diff. expression: an annotation edit reset the chosen metacells or cell types; lasso selections were not added to a non-empty group B; empty groups show a message.
* Gene modules: a gene could be added twice to a module; genes of a loaded file that are not in the dataset broke the module views; projections colored by a module follow its edits.
* Gene correlation: the cell type filter follows the session's annotation; a gene without counts got noise correlations on some platforms; more than 30 correlations per gene can be requested; the exported file keeps its header.
* Cell types: "Facet by: Cell type" failed with cell types on the x axis; the cell types are shown in their defined order, and cell types named by numbers are ordered 1, 2, ..., 10 (also in projects imported before).
* Annotate: the metacells of a deleted cell type are shown as "(Missing)" and can be annotated again; the tab is usable in a narrow window; editing gene modules no longer resets the gene, metadata field or module that other plots show.
* Flow: gene trajectories used a wrong time window when the flow didn't start at the first time bin.
* Query: "(dissimilar)" labels and the comparison with a previous projection are matched by metacell; "Show only fitted genes" of a group uses the group's cell types.
* Samples: diff. expression with a sample that lacks the chosen cell types shows a message instead of failing.
* `update_metadata_colors(overwrite = FALSE)` discarded the previous colors; `update_metadata()` failed when `metacell` was not the first column; `light_version` did not exclude the Gene correlation tab; the `update_*()` functions support old bundles with a `project` directory.
* QC: the value boxes, the percentile labels and the hover texts are readable at every window width; the ECDF axis label matches the curve. The control bar's tabs, the Gene correlation buttons, the Annotate color codes, the projection and scatter legends and the plot settings sidebars were fixed for narrow windows and many cell types. The Stdev-fold heatmap has a color scale that shows its values.

Performance:

* Heatmaps are drawn ~5 times faster, and are not drawn again when an unrelated input changes.
* Diff. expression of cell types and groups no longer copies the expression matrix; the comparison, its scatter and the projections of the Annotate, Diff. Expression and Query tabs are cached; selecting a row of the table highlights the gene without redrawing the scatter.
* Faster gene correlations, default markers, lookups of stored correlations and projections colored by cell type; sparse matrices stored as triplets are read as compressed columns.
* Default point sizes are rounded, so that users with different screens share the cached plots.
* Apps start without reading the expression matrix: the names of the genes and metacells are stored next to it (`mc_dimnames`; written at import, and at the first start of projects imported before, when their directory is writable). A 33K genes x 5.3K metacells app answers ~2 s sooner and uses 1.3 GB less memory until a tab needs expression data.
* Heatmaps: "Update genes" with local enrichment is 2 to 4 times faster and uses a fraction of the memory.
* The top correlated genes of a gene, the default colors of metadata fields and the gene labels are computed once per process.
* Diff. expression: the metacell selectors send ~3.5 times less data and no longer add an element per metacell to the page.

# MCView 0.2.64

* Apps open 2 to 3 s faster when they start.
* `MCVIEW_SHINY_CACHE_DIR` (or `shiny_cache_dir` in the config) keeps an app's plots across restarts, so it opens faster after one.

# MCView 0.2.63

* Projects imported before June 2023: the QC tab's max inner-fold and zero-fold plots work again, and projects imported with MCView 0.2.7 to 0.2.9 that have metadata no longer disconnect.

# MCView 0.2.62

Performance:

* Apps start about twice as fast and use 2 to 6 times less memory.
* Tabs load when first opened and send much less data to the browser.
* Heatmaps open without clustering the metacells (`update_default_clusterings()` for existing projects; `MCVIEW_CLUSTERING_CACHE_DIR` shares clusterings between a server's apps).
* Faster gene correlations, diff. expression, and cell type and sample comparisons.
* Recomputing the 2D projection no longer blocks the other users of the app.
* Restart a running app after re-importing its project or calling an `update_*()` function.

Fixes:

* Query tab: metacell mode tested the dataset's first metacell instead of the chosen one.
* Gene correlation tab: Module mode ignored the cell type filter, and "Individual genes" failed for a gene without counts.
* Samples tab: the sample types plot failed (since 0.2.48), diff. expression could fail or show stale results, and changing the projection coloring reset the chosen samples.
* Heatmaps: "Order by" was reset when the markers changed, zooming re-clustered the metacells, and the Markers heatmap failed when the project wasn't writable.
* Atlas tab: coloring by query cell type drew stacked duplicate points and didn't follow annotation edits.
* A failed 2D projection recompute, a zero-size plot, or a browser reconnecting after a restart could break the app.
* The clipboard-colored projection didn't refresh, and hovers on "Age" metadata showed the first metacell's age.
* Import stops with a clear message when Python anndata >= 0.13 or pandas >= 3 is installed.
* Multi-dataset projects: heatmaps, Annotate colors and the Genes metadata choices now follow a dataset switch.

# MCView 0.2.55

* Made `qs` an optional (`Suggests`) dependency instead of a hard requirement. New projects are serialized with `qs2`; `qs` is only needed to read legacy `.qs` caches from pre-0.2.54 projects, and a clear message now points to re-importing (or installing `qs`) when it is missing. This removes the archived-from-CRAN `qs` (and its `stringfish` build conflict) from the default install and the Docker image.
* Docker: fixed `start_app.R` file permissions so the container runs under an arbitrary non-root UID (Kubernetes / OpenShift). The launcher is now owned by `root:0` with group/world read+execute, avoiding `cannot open file './start_app.R': Permission denied`.

# MCView 0.2.54

* Switched serialization from `qs` to `qs2` (`qs` was removed from CRAN). New cache files are saved in `.qs2` format. Existing projects with `.qs` cache files are still supported via the legacy `qs` package (installed from GitHub).

# MCView 0.2.53

* Fix: stroke was not working in the atlas tab.
* Fix: "Gene", "Gene module" and "Metadata" selectors were not working in the atlas tab.

# MCView 0.2.52

* Fix: 's' was considered as a new line character in gene list input in the 'Gene correlation' tab.

# MCView 0.2.51

* Added categorical filtering to the Markers tab.
* Added sina/violin plots to the Cell types tab and added the option to set x axis to categorical metadata. 
* Added facets to the Cell types tab.

# MCView 0.2.50

* Added "Gene correlation" tab.

# MCView 0.2.49

* Added "Apply" button (and apply behaviour) to the 'Markers' tab. 

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


