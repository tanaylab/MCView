#' Create help modal content for the heatmap
#'
#' @param mode The mode of the heatmap (Markers, Inner, Stdev, etc.)
#' @return A modalDialog UI element
#' @noRd
create_heatmap_help_modal <- function(mode = "Markers") {
    modalDialog(
        title = paste0("Heatmap Help"),
        h4("About This Heatmap"),
        if (mode == "Markers") {
            tagList(
                p("Each row represents a gene, and each column represents a metacell. The colors indicate the log-fold change
          of gene expression relative to the median across all metacells."),
                p("Genes are sorted by their expression pattern to highlight cell type-specific marker genes.")
            )
        } else if (mode == "Inner") {
            tagList(
                p("This heatmap shows the inner-fold matrix."),
                p("For each metacell, for each gene, the inner_fold is the strongest (highest absolute value) deviant_fold of any of the cells contained in the metacell"),
                p("For each cell, for each gene, the deviant_fold holds the fold factor (log base 2) between the expression level of the gene in the cell and the metacell it belongs to (or the most similar metacell for outlier cells). This uses the same (strong) normalization factor we use when computing deviant (outlier) cells, so for outliers, you should see some (non-excluded, non-noisy) genes with a fold factor above 3 (8x), or some (non-excluded, noisy) genes with a fold factor above 5 (32x), which justify why we haven't merged that cell into a metacell; for cells grouped into metacells, you shouldn't see (many) such genes. If there is a large number of outlier cells and a few non-noisy genes have a high fold factor for many of them, you should consider marking these genes as noisy and recomputing the metacells. If they are already marked as noisy, you may want to completely exclude them.")
            )
        } else if (mode == "Stdev") {
            tagList(
                p("This heatmap shows the inner_stdev_log matrix."),
                p("For each metacell, for each gene, the standard deviation of the log (base 2) of the fraction of the gene across the cells of the metacell. Ideally, the standard deviation should be ~1/3rd of the deviants_min_gene_fold_factor (which is 3 by default), indicating that (all)most cells are within that maximal fold factor. In practice we may see higher values - the lower, the better. Both this and the inner_fold can be used for quality control over the consistency of the gene expression in the metacell.")
            )
        } else {
            p("")
        },
        h4("Metacell Ordering"),
        p("By default, metacells are ordered based on the value selected in 'Order by' in the sidebar (available under 'Order cell types by' dropdown)."),
        p("When 'Hierarchical clustering' is enabled, the dendrogram is sorted by the difference between the gene covering the largest number of metacells, and the gene that is the most anti-correlated to it."),
        p("When 'Force cell type' is enabled, metacells from the same cell type are grouped together. The order of cell types is determined by the 'Order cell types by' setting."),
        h4("Cell Type Filtering"),
        p("You can filter the heatmap to show only specific cell types by selecting them from the 'Cell types' dropdown in the sidebar."),
        h4("Color Scale"),
        p("The colors represent log-fold change values:"),
        tags$ul(
            tags$li(span(style = "color: blue; font-weight: bold", "Blue"), " - Under-expressed genes (negative values)"),
            tags$li(span(style = "color: white; font-weight: bold; background-color: #eee;", "White"), " - Neutral expression (values near zero)"),
            tags$li(span(style = "color: red; font-weight: bold", "Red"), " - Over-expressed genes (positive values)")
        ),
        h4("Gene Labels"),
        p("Gene labels are colored to indicate special categories:"),
        tags$ul(
            tags$li(span(style = "color: blue", "Blue"), " - Lateral genes (expressed in multiple cell types)"),
            tags$li(span(style = "color: red", "Red"), " - Noisy genes (high variability)"),
            tags$li(span(style = "color: purple", "Purple"), " - Both lateral and noisy"),
            tags$li(span(style = "color: darkgray", "Gray"), " - Disjoined genes"),
            tags$li(span(style = "color: black", "Black"), " - Standard marker genes"),
            tags$li(span(style = "color: #012901", "Dark green"), " - Module gene (when using gene modules)")
        ),
        p("Note that when the number of genes is large, some genes may appear only at one side of the heatmap."),
        h4("Gene Selection and Highlighting"),
        p("You can highlight specific genes by selecting them from the gene list in the sidebar and clicking 'Highlight selected genes'."),
        p("The 'Update genes' button recalculates the set of displayed genes based on the current heatmap view. Genes are selected based on:"),
        tags$ul(
            tags$li("Maximum log fraction threshold: Selects genes with at least one expression value above a threshold"),
            tags$li("Relative log fraction: Filters genes based on their expression variation relative to the median"),
            tags$li("Expression range: each metacell is choosing the top N genes with the highest expression range")
        ),
        p("The number of genes displayed can be adjusted using the 'Number of genes' on the left sidebar."),
        h4("Enrichment Display Options"),
        p("When viewing a sub-matrix (after zooming in or filtering), you can choose between:"),
        tags$ul(
            tags$li(strong("Global"), " - Shows enrichment values calculated from the complete dataset"),
            tags$li(strong("Local"), " - Recalculates enrichment values using only the visible metacells")
        ),
        h4("Metadata Annotations"),
        p("You can add metadata annotations to the heatmap by selecting metadata fields from the 'Metadata' dropdown in the sidebar."),
        h4("Interacting With The Heatmap"),
        tags$ul(
            tags$li(strong("Hover"), " - See details about genes and metacells"),
            tags$li(strong("Double-click"), " - Select a gene or metacell for further analysis"),
            tags$li(strong("Brush"), " - Select a range of metacells to zoom or analyze")
        ),
        h4("Cell Type Bars"),
        p("The colored bars above and below the heatmap indicate the cell type for each metacell column."),
        h4("Advanced Features"),
        p("Additional functionality available in the left sidebar:"),
        tags$ul(
            tags$li(strong("Highlight genes"), " - Select genes from the list and click to highlight them in the plot"),
            tags$li(strong("Brush action"), " - Choose between zooming in on a region or selecting metacells"),
            tags$li(strong("Include lateral / noisy"), " - Include lateral and noisy genes in the heatmap"),
            tags$li("Load and save gene sets for later use"),
            tags$li("Download the heatmap data as a CSV file")
        ),
        p("Additional options are available in the \"settings\" on the right:"),
        tags$ul(
            tags$li(strong("Filter by clipboard"), " - Show only metacells that are currently in the clipboard"),
            tags$li(strong("Fold change range"), " - Adjust the color scale range"),
            tags$li(strong("Other settings"), " - Additional options for the heatmap")
        ),
        footer = tagList(
            modalButton("Close")
        ),
        size = "l",
        easyClose = TRUE,
        fade = TRUE
    )
}
