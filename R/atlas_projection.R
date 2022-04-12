import_atlas <- function(query, atlas_project, atlas_dataset, projection_weights_file, dataset, cache_dir, copy_atlas) {
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

    required_fields <- c("projected_type", "similar")
    query_md <- query$obs %>%
        rownames_to_column("metacell") %>%
        as_tibble()

    purrr::walk(required_fields, ~ {
        if (!has_name(query_md, .x)) {
            cli_abort("Query h5ad file does not have the required field: '{.file {.x}}'")
        }
    })

    proj_weights <- tgutil::fread(projection_weights_file) %>%
        mutate(atlas = as.character(atlas), query = as.character(query)) %>%
        as_tibble()
    if (!all(rlang::has_name(proj_weights, c("query", "atlas", "weight")))) {
        cli_abort(".{file {projection_weights_file}} should have fields named 'query', 'atlas' and 'weight'")
    }
    serialize_shiny_data(proj_weights, "proj_weights", dataset = dataset, cache_dir = cache_dir)

    proj_types <- unique(query_md$projected_type)
    atlas_metacell_types <- get_mc_data(dataset, "metacell_types", atlas = TRUE)

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
    prev_metadata <- get_mc_data(dataset, "metadata")
    if (is.null(prev_metadata)) {
        metadata <- proj_metacell_types %>%
            select(metacell, similar)
    } else {
        if (has_name(prev_metadata, "similar")) {
            if (!all(prev_metadata$similar == proj_metacell_types$similar)) {
                cli_abort("Metacell metadata includes a field named {.field similar}. Please rename it in order to run MCView in atlas mode.")
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

    if (is.null(query$layers[["projected"]])) {
        cli_abort("Query h5ad is missing the '{.file projected}' layer")
    }
    projected_mat <- t(query$layers[["projected"]])
    rownames(projected_mat) <- colnames(query$X)
    colnames(projected_mat) <- rownames(query$X)
    serialize_shiny_data(projected_mat, "projected_mat", dataset = dataset, cache_dir = cache_dir)

    projected_mat_sum <- colSums(projected_mat)
    serialize_shiny_data(projected_mat_sum, "projected_mat_sum", dataset = dataset, cache_dir = cache_dir)

    if (is.null(query$layers[["projected_fold"]])) {
        cli_abort("Query h5ad is missing the '{.file projected_fold}' layer")
    }

    projected_fold <- t(query$layers[["projected_fold"]])
    rownames(projected_fold) <- colnames(query$X)
    colnames(projected_fold) <- rownames(query$X)
    serialize_shiny_data(projected_fold, "projected_fold", dataset = dataset, cache_dir = cache_dir)

    cli_alert_info("Calculating top atlas-query fold genes")
    forbidden <- query$var$forbidden_gene    
    marker_genes_projected <- select_top_fold_genes_per_metacell(projected_fold[!forbidden, ], minimal_relative_log_fraction = -Inf, use_abs = TRUE, genes_per_metacell = 50)
    serialize_shiny_data(marker_genes_projected, "marker_genes_projected", dataset = dataset, cache_dir = cache_dir)

    serialize_shiny_data(query$uns$project_max_projection_fold_factor, "project_max_projection_fold_factor", dataset = dataset, cache_dir = cache_dir)

    # disjoined genes
    atlas_mat <- get_mc_data(dataset, "mc_mat", atlas = TRUE)
    query_mat <- get_mc_data(dataset, "mc_mat")
    disjoined_genes_no_atlas <- setdiff(rownames(query_mat), rownames(atlas_mat))
    disjoined_genes_no_query <- setdiff(rownames(atlas_mat), rownames(query_mat))

    serialize_shiny_data(disjoined_genes_no_atlas, "disjoined_genes_no_atlas", dataset = dataset, cache_dir = cache_dir)
    serialize_shiny_data(disjoined_genes_no_query, "disjoined_genes_no_query", dataset = dataset, cache_dir = cache_dir)

    if (is.null(query$var$systematic_gene)) {
        systematic_genes <- c()
    } else {
        systematic_genes <- rownames(query$var)[query$var$systematic_gene]
    }

    serialize_shiny_data(systematic_genes, "systematic_genes", dataset = dataset, cache_dir = cache_dir)

    # Gene metadata
    gene_md <- query$var %>%
        rownames_to_column("gene") %>%
        as_tibble()

    cell_type_gene_md <- gene_md %>%
        select(
            gene,
            contains("_of_"),
            forbidden_gene,
            ignored_gene,
            any_of("atlas_significant_gene"),
            correction_factor,
            projected_correlation,
            correlated_gene,
            uncorrelated_gene,
            glob_biased_gene = biased_gene,
            glob_ignored_gene = ignored_gene,
            glob_systematic_gene = systematic_gene
        ) %>%
        pivot_longer(
            contains("_of_"),
            names_to = c("dtype", "cell_type"),
            names_pattern = c("^(.+)_of_(.+)$")
        ) %>%
        spread(dtype, value) %>%
        select(
            gene, cell_type, biased_gene, systematic_gene, ignored_gene, projected_correlation, correlated_gene, uncorrelated_gene, everything()
        ) %>%
        arrange(
            desc(biased_gene),
            desc(systematic_gene),
            desc(glob_biased_gene),
            desc(glob_systematic_gene),
            desc(ignored_gene),
            desc(glob_ignored_gene)
        )

    serialize_shiny_data(cell_type_gene_md, "gene_metadata", dataset = dataset, cache_dir = cache_dir, flat = TRUE)


    cli_alert_success("succesfully imported atlas projections")
}


plot_type_predictions_bar <- function(dataset) {
    plotly::renderPlotly({
        df_fracs <- get_mc_data(dataset(), "query_atlas_cell_type_fracs")
        req(!is.null(df_fracs))
        req(has_atlas(dataset()))
        atlas_colors <- get_mc_data(dataset(), "cell_type_colors", atlas = TRUE)
        atlas_colors <- atlas_colors %>%
            select(cell_type, color) %>%
            deframe()

        fracs_mat <- df_fracs %>%
            spread(type, fraction) %>%
            column_to_rownames("metacell") %>%
            as.matrix()

        ord <- slanter::slanted_orders(fracs_mat)

        df_fracs <- df_fracs %>%
            mutate(metacell = factor(metacell, levels = rownames(fracs_mat)[ord$rows])) %>%
            mutate(type = factor(type, levels = colnames(fracs_mat)[ord$cols]))

        p <- df_fracs %>%
            ggplot(aes(x = metacell, y = fraction, fill = type, customdata = metacell)) +
            geom_col() +
            scale_fill_manual(name = "", values = atlas_colors) +
            theme(
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
            ) +
            scale_y_continuous(expand = c(0, 0)) +
            ylab("Fraction") +
            xlab("Metacell") +
            guides(fill = "none")


        fig <- plotly::ggplotly(
            p,
            source = "type_prediction_bar"
        ) %>%
            sanitize_plotly_buttons() %>%
            plotly::hide_legend()

        return(fig)
    })
}
