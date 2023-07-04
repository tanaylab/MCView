import_atlas <- function(query, atlas_project, atlas_dataset, projection_weights_file, dataset, cache_dir, copy_atlas, gene_names = NULL) {
    cli_alert_info("Reading dataset {.file {atlas_dataset}} at project: {.file {atlas_project}}")
    if (!fs::dir_exists(fs::path(project_cache_dir(atlas_project), atlas_dataset))) {
        # check if project cache dir only has a single dataset
        if (length(fs::dir_ls(project_cache_dir(atlas_project))) == 1) {
            atlas_dataset_new <- basename(fs::dir_ls(project_cache_dir(atlas_project))[[1]])
            cli_alert_info("Project {.file {atlas_project}} only has a single dataset ({.file {atlas_dataset_new}}), using that one instead of {.file {atlas_dataset}}")
            atlas_dataset <- atlas_dataset_new
        } else {
            cli_abort("Atlas dataset {.file {atlas_dataset}} does not exist at project: {.file {atlas_project}}")
        }
    }
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

    required_fields <- c("projected_type", "similar")
    query_md <- query$obs %>%
        rownames_to_column("metacell") %>%
        as_tibble()

    purrr::walk(required_fields, ~ {
        if (!has_name(query_md, .x)) {
            cli_abort("Query h5ad file does not have the required field: '{.file {.x}}'")
        }
    })

    # disjoined genes
    atlas_mat <- get_mc_data(dataset, "mc_mat", atlas = TRUE)
    query_mat <- get_mc_data(dataset, "mc_mat")
    disjoined_genes_no_atlas <- setdiff(rownames(query_mat), rownames(atlas_mat))
    disjoined_genes_no_query <- setdiff(rownames(atlas_mat), rownames(query_mat))

    serialize_shiny_data(disjoined_genes_no_atlas, "disjoined_genes_no_atlas", dataset = dataset, cache_dir = cache_dir)
    serialize_shiny_data(disjoined_genes_no_query, "disjoined_genes_no_query", dataset = dataset, cache_dir = cache_dir)

    # change the UMI matrix to use the corrected UMI counts
    mc_sum <- query$obs$total_atlas_umis
    serialize_shiny_data(mc_sum, "mc_sum", dataset = dataset, cache_dir = cache_dir)

    mc_mat_corrected <- t(query$layers[["corrected_fraction"]] * mc_sum)
    if (!is.null(gene_names)) {
        rownames(mc_mat_corrected) <- modify_gene_names(rownames(mc_mat_corrected), gene_names)
    }

    # if corrected_fraction is a sparse matrix, convert it to dense
    if (methods::is(mc_mat_corrected, "dgCMatrix")) {
        mc_mat_corrected <- as.matrix(mc_mat_corrected)
    }
    serialize_shiny_data(mc_mat_corrected, "mc_mat_corrected", dataset = dataset, cache_dir = cache_dir)

    proj_weights <- tgutil::fread(projection_weights_file) %>%
        mutate(atlas = as.character(atlas), query = as.character(query)) %>%
        as_tibble()
    if (!all(rlang::has_name(proj_weights, c("query", "atlas", "weight")))) {
        cli_abort(".{file {projection_weights_file}} should have fields named 'query', 'atlas' and 'weight'")
    }
    serialize_shiny_data(proj_weights, "proj_weights", dataset = dataset, cache_dir = cache_dir)

    proj_types <- unique(query_md$projected_type)
    atlas_metacell_types <- get_mc_data(dataset, "metacell_types", atlas = TRUE)
    validate_proj_weights(proj_weights, colnames(query_mat), atlas_metacell_types$metacell)

    query_atlas_cell_type_fracs <- proj_weights %>%
        left_join(
            atlas_metacell_types %>% select(atlas = metacell, type = cell_type),
            by = "atlas"
        ) %>%
        group_by(query, type) %>%
        summarise(fraction = sum(weight), .groups = "drop") %>%
        rename(metacell = query)

    query_atlas_cell_type_fracs <- query_atlas_cell_type_fracs %>%
        tidyr::complete(metacell, type, fill = list(fraction = 0))

    serialize_shiny_data(query_atlas_cell_type_fracs, "query_atlas_cell_type_fracs", dataset = dataset, cache_dir = cache_dir, flat = TRUE)

    proj_metacell_types <- query_md %>%
        select(metacell, cell_type = projected_type, similar)

    # Take metacell metadata from the query
    query_types <- get_mc_data(dataset, "metacell_types", atlas = FALSE) %>%
        select(metacell, top1_gene, top2_gene, top1_lfp, top2_lfp)
    proj_metacell_types <- proj_metacell_types %>% left_join(query_types, by = "metacell")

    # Take the colors from the atlas
    atlas_colors <- get_mc_data(dataset, "cell_type_colors", atlas = TRUE)
    proj_metacell_types <- proj_metacell_types %>%
        left_join(atlas_colors %>% distinct(cell_type, mc_col = color), by = "cell_type")

    serialize_shiny_data(proj_metacell_types, "projected_metacell_types", dataset = dataset, cache_dir = cache_dir, flat = TRUE)

    # Add similar to regular metadata
    prev_metadata <- get_mc_data(dataset, "metadata", atlas = FALSE)

    if (is.null(prev_metadata)) {
        metadata <- proj_metacell_types %>%
            select(metacell, similar)
    } else {
        if (has_name(prev_metadata, "similar")) {
            if (!all(prev_metadata$similar == proj_metacell_types$similar)) {
                cli_warn("Metacell metadata includes a field named {.field similar}. Overwriting it with the one from the projection file.")
                prev_metadata$similar <- proj_metacell_types$similar
            }
            metadata <- prev_metadata
        } else {
            metadata <- prev_metadata %>%
                mutate(metacell = as.character(metacell)) %>%
                left_join(
                    proj_metacell_types %>% mutate(metacell = as.character(metacell)) %>% select(metacell, similar),
                    by = "metacell"
                )
        }
    }
    metadata <- metadata %>%
        mutate(similar = ifelse(similar, "similar", "dissimilar"))
    serialize_shiny_data(metadata, "metadata", dataset = dataset, cache_dir = cache_dir, flat = TRUE)

    # set colors for similar
    similarity_colors <- c("similar" = "darkgreen", "dissimilar" = "darkred")
    metadata_colors <- get_mc_data(dataset, "metadata_colors")
    if (is.null(metadata_colors)) {
        metadata_colors <- list()
    }
    metadata_colors[["similar"]] <- similarity_colors
    serialize_shiny_data(metadata_colors, "metadata_colors", dataset = dataset, cache_dir = cache_dir)


    if (is.null(query$layers[["projected_fraction"]])) {
        cli_abort("Query h5ad is missing the '{.file projected_fraction}' layer")
    }
    projected_mat_sum <- query$obs$total_atlas_umis
    serialize_shiny_data(projected_mat_sum, "projected_mat_sum", dataset = dataset, cache_dir = cache_dir)
    projected_mat <- t(query$layers[["projected_fraction"]] * projected_mat_sum)
    rownames(projected_mat) <- colnames(query$X)
    colnames(projected_mat) <- rownames(query$X)
    if (!is.null(gene_names)) {
        rownames(projected_mat) <- modify_gene_names(rownames(projected_mat), gene_names)
    }
    serialize_shiny_data(projected_mat, "projected_mat", dataset = dataset, cache_dir = cache_dir)

    if (is.null(query$layers[["projected_fold"]])) {
        cli_abort("Query h5ad is missing the '{.file projected_fold}' layer")
    }

    projected_fold <- t(query$layers[["projected_fold"]])
    rownames(projected_fold) <- colnames(query$X)
    colnames(projected_fold) <- rownames(query$X)
    if (!is.null(gene_names)) {
        rownames(projected_fold) <- modify_gene_names(rownames(projected_fold), gene_names)
    }
    serialize_shiny_data(projected_fold, "projected_fold", dataset = dataset, cache_dir = cache_dir)

    cli_alert_info("Calculating top atlas-query fold genes")
    lateral_field <- "lateral_gene"
    lateral <- query$var[, lateral_field]

    marker_genes_projected <- select_top_fold_genes_per_metacell(projected_fold[!lateral, ], minimal_relative_log_fraction = -Inf, use_abs = TRUE, genes_per_metacell = 50)
    serialize_shiny_data(marker_genes_projected, "marker_genes_projected", dataset = dataset, cache_dir = cache_dir)

    serialize_shiny_data(query$uns$project_max_projection_fold_factor, "project_max_projection_fold_factor", dataset = dataset, cache_dir = cache_dir)

    # Gene metadata
    gene_md <- query$var %>%
        rownames_to_column("gene") %>%
        as_tibble()

    cell_type_gene_md <- gene_md %>%
        select(
            gene,
            contains("_of_"),
            any_of("lateral_gene"),
            any_of("atlas_lateral_gene"),
            any_of("atlas_marker_gene"),
            any_of("correction_factor"),
            any_of("correlated_gene")
        ) %>%
        pivot_longer(
            contains("_of_"),
            names_to = c("dtype", "cell_type"),
            names_pattern = c("^(.+)_of_(.+)$")
        ) %>%
        spread(dtype, value) %>%
        select(
            gene, cell_type, everything()
        )

    if (!is.null(gene_names)) {
        cell_type_gene_md$gene <- modify_gene_names(cell_type_gene_md$gene, gene_names)
    }

    serialize_shiny_data(cell_type_gene_md, "gene_metadata", dataset = dataset, cache_dir = cache_dir, flat = TRUE)


    cli_alert_success("succesfully imported atlas projections")
}


validate_proj_weights <- function(proj_weights, query_mcs, atlas_mcs) {
    # validate that all proj_weights$query metacells are in query_mcs
    missing_mcs <- setdiff(proj_weights$query, query_mcs)
    if (length(missing_mcs) > 0) {
        cli_abort("The following query metacells are missing from the query h5ad: {.val {missing_mcs}}")
    }

    # validate that all query_mcs are in proj_weights$query
    missing_mcs <- setdiff(query_mcs, proj_weights$query)
    if (length(missing_mcs) > 0) {
        cli_abort("The following query metacells are missing from the projection file: {.val {missing_mcs}}")
    }

    # validate that all proj_weights$atlas metacells are in atlas_mcs
    missing_mcs <- setdiff(proj_weights$atlas, atlas_mcs)
    if (length(missing_mcs) > 0) {
        cli_abort("The following atlas metacells are missing from the atlas h5ad: {.val {missing_mcs}}")
    }
}
