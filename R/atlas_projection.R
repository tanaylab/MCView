import_atlas <- function(query, atlas_project, atlas_dataset, dataset, cache_dir, copy_atlas) {
    cli_alert_info("Reading dataset {.file {atlas_dataset}} at project: {.file {atlas_project}}")
    verify_app_cache(atlas_project, datasets = atlas_dataset)

    atlas_new_path <- fs::path(cache_dir, dataset, "atlas")
    if (copy_atlas) {
        fs::dir_copy(
            fs::path(project_cache_dir(atlas_project), atlas_dataset),
            atlas_new_path,
            overwrite = TRUE
        )
    } else {
        if (fs::file_exists(atlas_new_path)) {
            fs::file_remove(atlas_new_path)
        }

        fs::link_create(
            normalizePath(fs::path(project_cache_dir(atlas_project), atlas_dataset)),
            atlas_new_path
        )
    }

    load_all_data(cache_dir, datasets = dataset)

    required_fields <- c("type", "charted")
    query_md <- query$obs %>%
        rownames_to_column("metacell") %>%
        as_tibble()

    proj_types <- unique(query_md$type)
    required_fields <- c(required_fields, paste0("fraction_", proj_types))

    purrr::walk(required_fields, ~ {
        if (!has_name(query_md, .x)) {
            cli_abort("Query h5ad file does not have the required field: '{.file {.x}}'")
        }
    })

    query_atlas_cell_type_fracs <- query_md %>%
        select(metacell, starts_with("fraction_")) %>%
        pivot_longer(
            cols = starts_with("fraction_"),
            names_to = "type",
            names_prefix = "fraction_",
            values_to = "fraction",
        )

    serialize_shiny_data(query_atlas_cell_type_fracs, "query_atlas_cell_type_fracs", dataset = dataset, cache_dir = cache_dir, flat = TRUE)

    proj_metacell_types <- query_md %>%
        select(metacell, cell_type = type, charted)

    # Take metacell metadata from the query
    query_types <- get_mc_data(dataset, "metacell_types", atlas = FALSE) %>%
        select(metacell, top1_gene, top2_gene, top1_lfp, top2_lfp)
    proj_metacell_types <- proj_metacell_types %>% left_join(query_types, by = "metacell")

    # Take the colors from the atlas
    atlas_colors <- get_mc_data(dataset, "cell_type_colors", atlas = TRUE)
    proj_metacell_types <- proj_metacell_types %>%
        left_join(atlas_colors %>% distinct(cell_type, mc_col = color), by = "cell_type")

    serialize_shiny_data(proj_metacell_types, "projected_metacell_types", dataset = dataset, cache_dir = cache_dir, flat = TRUE)

    # Add charted to regular metadata
    prev_metadata <- get_mc_data(dataset, "metadata")
    if (is.null(prev_metadata)) {
        metadata <- proj_metacell_types %>%
            select(metacell, charted)
    } else {
        if (has_name(prev_metadata, "charted")) {
            if (!all(prev_metadata$charted == proj_metacell_types$charted)) {
                cli_abort("Metacell metadata includes a field named {.field charted}. Please rename it in order to run MCView in atlas mode.")
            }
            metadata <- prev_metadata
        } else {
            metadata <- prev_metadata %>% left_join(
                proj_metacell_types %>% select(metacell, charted),
                by = "metacell"
            )
        }
    }
    metadata <- metadata %>%
        mutate(charted = ifelse(charted, "charted", "un-charted"))
    serialize_shiny_data(metadata, "metadata", dataset = dataset, cache_dir = cache_dir, flat = TRUE)

    # set colors for charted
    charting_colors <- c("charted" = "darkgreen", "un-charted" = "darkred")
    metadata_colors <- get_mc_data(dataset, "metadata_colors")
    if (is.null(metadata_colors)) {
        metadata_colors <- list()
    }
    metadata_colors[["charted"]] <- charting_colors
    serialize_shiny_data(metadata_colors, "metadata_colors", dataset = dataset, cache_dir = cache_dir)

    if (is.null(query$layers[["projected"]])) {
        cli_abort("Query h5ad is missing the '{.file projected}' layer")
    }

    projected_mat <- Matrix::t(query$layers[["projected"]])
    serialize_shiny_data(projected_mat, "projected_mat", dataset = dataset, cache_dir = cache_dir)

    projected_mat_sum <- colSums(projected_mat)
    serialize_shiny_data(projected_mat_sum, "projected_mat_sum", dataset = dataset, cache_dir = cache_dir)

    if (is.null(query$layers[["projected_fold"]])) {
        cli_abort("Query h5ad is missing the '{.file projected_fold}' layer")
    }

    projected_fold <- Matrix::t(query$layers[["projected_fold"]])
    serialize_shiny_data(projected_fold, "projected_fold", dataset = dataset, cache_dir = cache_dir)

    cli_alert_info("Calculating top atlas-query fold genes")
    forbidden <- query$var$forbidden_gene
    marker_genes_projected <- select_top_fold_genes(projected_fold[!forbidden, ])
    serialize_shiny_data(marker_genes_projected, "marker_genes_projected", dataset = dataset, cache_dir = cache_dir)

    # disjoined genes
    atlas_mat <- get_mc_data(dataset, "mc_mat", atlas = TRUE)
    query_mat <- get_mc_data(dataset, "mc_mat")
    disjoined_genes <- setdiff(rownames(atlas_mat), rownames(query_mat))
    serialize_shiny_data(disjoined_genes, "disjoined_genes", dataset = dataset, cache_dir = cache_dir)

    # TODO: systematic genes

    cli_alert_success("succesfully imported atlas projections")
}
