load_outliers <- function(outliers_anndata_file, project, dataset, gene_names = NULL) {
    cache_dir <- project_cache_dir(project)
    cli_alert_info("Loading anndata file: {.file {outliers_anndata_file}}")
    outliers <- anndata::read_h5ad(outliers_anndata_file)
    cli_alert_info("Processing outliers")

    # outliers_mat <- t(outliers$X)
    # rownames(outliers_mat) <- modify_gene_names(rownames(outliers_mat), gene_names)
    # serialize_shiny_data(outliers_mat, "outliers_mat", dataset = dataset, cache_dir = cache_dir)

    outliers_metadata <- outliers$obs %>%
        rownames_to_column("cell_id") %>%
        as_tibble()
    serialize_shiny_data(outliers_metadata, "outliers_metadata", dataset = dataset, cache_dir = cache_dir, flat = TRUE)

    if (!has_name(outliers_metadata, "most_similar")) {
        cli_abort("File {.file {outliers_anndata_file}} does not contain the column {.field most_similar}.")
    }

    if (is.null(outliers$layers[["deviant_fold"]])) {
        cli_abort("File {.file {outliers_anndata_file}} does not contain the layer {.field deviant_fold}.")
    }
    deviant_fold_mat <- t(outliers$layers[["deviant_fold"]])
    rownames(deviant_fold_mat) <- modify_gene_names(rownames(deviant_fold_mat), gene_names)
    serialize_shiny_data(deviant_fold_mat, "deviant_fold_mat", dataset = dataset, cache_dir = cache_dir)

    cli_alert_info("Calculating top inner-fold genes")
    deviant_fold_genes <- rownames(deviant_fold_mat)[Matrix::rowSums(deviant_fold_mat) > 0]
    deviant_fold_gene_metacells <- matrixStats::rowSums2(as.matrix(deviant_fold_mat[deviant_fold_genes, , drop = FALSE]) > 0)

    # fp here is the number of non-zero entries per gene
    marker_genes_deviant_fold <- tibble(gene = deviant_fold_genes, fp = deviant_fold_gene_metacells) %>%
        arrange(desc(fp))
    serialize_shiny_data(marker_genes_deviant_fold, "marker_genes_deviant_fold", dataset = dataset, cache_dir = cache_dir)

    cli_alert_info("Add the {.field \"Outliers\"} tab to your config file to view the outliers.")
}
