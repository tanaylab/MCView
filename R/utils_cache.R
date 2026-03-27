# ==============================================================================
# DAF-Based Data Access Functions
# ==============================================================================

get_cell_type_data <- function(dataset, atlas = FALSE) {
    daf_obj <- get_daf_for_query(dataset, atlas)
    if (is.null(daf_obj)) {
        return(NULL)
    }

    cell_type_colors <- convert_daf_to_mcview(daf_obj, "cell_type_colors", atlas)

    if (is.null(cell_type_colors)) {
        return(NULL)
    }

    cell_type_colors <- cell_type_colors %>%
        filter(!is.na(cell_type)) %>%
        mutate(cell_type = as.character(cell_type))
    return(cell_type_colors)
}

ensure_metacell_types_fields <- function(metacell_types) {
    if (is.null(metacell_types)) {
        return(NULL)
    }

    metacell_types <- as_tibble(metacell_types)

    if (!has_name(metacell_types, "metacell") && !is.null(rownames(metacell_types))) {
        metacell_types <- metacell_types %>%
            tibble::rownames_to_column("metacell")
    }

    if (!has_name(metacell_types, "top1_gene")) {
        metacell_types$top1_gene <- NA_character_
    }
    if (!has_name(metacell_types, "top2_gene")) {
        metacell_types$top2_gene <- NA_character_
    }
    if (!has_name(metacell_types, "top1_lfp")) {
        metacell_types$top1_lfp <- NA_real_
    }
    if (!has_name(metacell_types, "top2_lfp")) {
        metacell_types$top2_lfp <- NA_real_
    }

    if (!has_name(metacell_types, "mc_col")) {
        if (has_name(metacell_types, "mc_col.x") && has_name(metacell_types, "mc_col.y")) {
            metacell_types$mc_col <- ifelse(
                is.na(metacell_types$mc_col.x),
                metacell_types$mc_col.y,
                metacell_types$mc_col.x
            )
        } else if (has_name(metacell_types, "mc_col.x")) {
            metacell_types$mc_col <- metacell_types$mc_col.x
        } else if (has_name(metacell_types, "mc_col.y")) {
            metacell_types$mc_col <- metacell_types$mc_col.y
        }
    }

    metacell_types <- metacell_types %>%
        select(-any_of(c("mc_col.x", "mc_col.y")))

    metacell_types
}

get_metacell_types_data <- function(dataset, atlas = FALSE) {
    # Return cached result if available (metacell_types is static per session)
    cache_key <- if (atlas) "metacell_types_atlas" else "metacell_types"
    mc_data <- mcv_get("mc_data")
    cached <- mc_data[[dataset]][[cache_key]]
    if (!is.null(cached)) {
        return(cached)
    }

    daf_obj <- get_daf_for_query(dataset, atlas)
    if (is.null(daf_obj)) {
        return(NULL)
    }

    metacell_types <- convert_daf_to_mcview(daf_obj, "metacell_types", atlas)

    if (is.null(metacell_types)) {
        return(NULL)
    }

    if (!is.factor(metacell_types$cell_type)) {
        metacell_types$cell_type <- factor(metacell_types$cell_type)
    }
    metacell_types <- metacell_types %>%
        mutate(cell_type = as.character(forcats::fct_na_value_to_level(cell_type, "(Missing)"))) %>%
        mutate(metacell = as.character(metacell))

    # convert_daf_metacell_types() already populates mc_col from the DAF
    # (either per-metacell mc_col vector or type-level color vector).
    # Only fetch cell_type_colors and join if mc_col is missing, avoiding
    # a redundant DAF round-trip.
    if (!has_name(metacell_types, "mc_col")) {
        cell_type_colors <- get_cell_type_data(dataset, atlas = atlas)
        metacell_types <- metacell_types %>%
            left_join(cell_type_colors %>% select(cell_type, mc_col = color), by = "cell_type")
    }

    metacell_types <- ensure_metacell_types_fields(metacell_types)

    # Cache the result for subsequent calls
    mc_data <- mcv_get("mc_data")
    mc_data[[dataset]][[cache_key]] <- metacell_types
    mcv_set("mc_data", mc_data)

    return(metacell_types)
}

get_mc_color_key <- function(dataset) {
    get_metacell_types_data(dataset) %>%
        distinct(cell_type, mc_col) %>%
        rename(group = cell_type, color = mc_col) %>%
        as.data.frame()
}

get_metadata <- function(dataset, atlas = FALSE) {
    daf_obj <- get_daf_for_query(dataset, atlas)
    if (is.null(daf_obj)) {
        return(NULL)
    }

    metadata <- convert_daf_to_mcview(daf_obj, "metadata", atlas)

    if (!is.null(metadata)) {
        mc_qc_metadata <- convert_daf_to_mcview(daf_obj, "mc_qc_metadata", atlas)
        if (!is.null(mc_qc_metadata) && ncol(mc_qc_metadata) > 1) {
            # Both tibbles share the same metacell axis in the same order,
            # so we can bind columns directly instead of left_join().
            # left_join forces materialization of jlview ALTREP columns;
            # direct column assignment preserves zero-copy views.
            existing_cols <- colnames(metadata)
            for (col in colnames(mc_qc_metadata)) {
                if (col != "metacell" && !(col %in% existing_cols)) {
                    metadata[[col]] <- mc_qc_metadata[[col]]
                }
            }
        }
    }

    return(metadata)
}

#' Get metacell data from DAF
#'
#' @param dataset Dataset name
#' @param var_name Variable name (e.g., "mc_mat", "mc_sum", "mc2d")
#' @param atlas Whether to get atlas data
#' @return Data in MCView format, or NULL if not available
#' @export
get_mc_data <- function(dataset, var_name, atlas = FALSE) {
    # Handle special cases that need post-processing
    if (var_name == "metacell_types") {
        return(get_metacell_types_data(dataset, atlas = atlas))
    } else if (var_name == "cell_type_colors") {
        return(get_cell_type_data(dataset, atlas = atlas))
    } else if (var_name == "metadata") {
        return(get_metadata(dataset, atlas = atlas))
    }

    # Memoize DERIVED data that requires non-trivial assembly beyond raw DAF
    # reads.  Variables that are thin wrappers around a single dafr::get_vector()
    # or dafr::get_matrix() call are NOT cached here because dafr's own R-side
    # cache (keyed by version counter) already avoids the Julia round-trip.
    #
    # Removed from cache (dafr handles these):
    #   mc_sum         - single get_vector + names()
    #   lateral_genes  - get_vector + flag filter
    #   noisy_genes    - get_vector + flag filter
    #   projected_fold  - single get_matrix
    static_vars <- c(
        "mc_mat", "mc2d",
        "marker_genes", "gene_modules",
        "marker_genes_projected",
        "mc_qc_metadata", "gene_qc", "gg_mc_top_cor",
        "metacell_graphs", "qc_stats", "cell_metadata",
        "mc_mat_corrected", "projected_mat",
        "proj_weights", "query_atlas_cell_type_fracs", "query_md"
    )

    cache_key <- if (atlas) paste0(var_name, "_atlas") else var_name
    if (var_name %in% static_vars) {
        mc_data <- mcv_get("mc_data")
        cached <- mc_data[[dataset]][[cache_key]]
        if (!is.null(cached)) {
            return(cached)
        }
    }

    daf_obj <- get_daf_for_query(dataset, atlas)
    if (is.null(daf_obj)) {
        return(NULL)
    }

    result <- convert_daf_to_mcview(daf_obj, var_name, atlas)

    # Store in mc_data cache for static variables
    if (!is.null(result) && var_name %in% static_vars) {
        mc_data <- mcv_get("mc_data")
        mc_data[[dataset]][[cache_key]] <- result
        mcv_set("mc_data", mc_data)
    }

    result
}


get_mc_config <- function(dataset, var_name) {
    config <- mcv_get("config")
    if (is.null(config$datasets)) {
        return(NULL)
    }
    if (is.null(config$datasets[[dataset]])) {
        return(NULL)
    }
    config$datasets[[dataset]][[var_name]]
}

has_network <- function(dataset) {
    # Check DAF axis instead of loading full graph data
    daf_obj <- get_dataset_daf(dataset)
    if (!is.null(daf_obj)) {
        return(dafr::has_axis(daf_obj, "metacell_graph"))
    }
    !is.null(get_mc_data(dataset, "mc_network"))
}

has_time <- function(dataset) {
    daf_obj <- get_dataset_daf(dataset)
    if (is.null(daf_obj)) return(FALSE)
    tryCatch({
        dafr::has_vector(daf_obj, "metacell", "time")
    }, error = function(e) FALSE)
}

has_metadata <- function(dataset, atlas = FALSE) {
    !is.null(get_mc_data(dataset, "metadata", atlas = atlas))
}

has_cell_metadata <- function(dataset) {
    # Check DAF axis instead of loading full cell metadata
    daf_obj <- get_dataset_daf(dataset)
    if (!is.null(daf_obj)) {
        return(dafr::has_axis(daf_obj, "cell"))
    }
    !is.null(get_mc_data(dataset, "cell_metadata"))
}

has_samples <- function(dataset) {
    # Check DAF cell vector existence: mcview_sample_property scalar names the
    # source vector, falling back to literal "samp_id".
    daf_obj <- get_dataset_daf(dataset)
    if (!is.null(daf_obj) && dafr::has_axis(daf_obj, "cell")) {
        # Check mcview_sample_property scalar first
        if (dafr::has_scalar(daf_obj, "mcview_sample_property")) {
            prop_name <- dafr::get_scalar(daf_obj, "mcview_sample_property")
            if (dafr::has_vector(daf_obj, "cell", prop_name)) {
                return(TRUE)
            }
        }
        if (dafr::has_vector(daf_obj, "cell", "samp_id")) {
            return(TRUE)
        }
    }
    # Also check for cell-level grouping fields from the cells DAF
    if (has_cell_gene_umis(dataset) && length(get_cell_grouping_fields(dataset)) > 0) {
        return(TRUE)
    }
    if (!has_cell_metadata(dataset)) {
        return(FALSE)
    }
    cell_md <- get_mc_data(dataset, "cell_metadata")
    return(rlang::has_name(cell_md, "samp_id"))
}

has_projection <- function(dataset) {
    # Check DAF vector existence instead of loading full query metadata
    daf_obj <- get_dataset_daf(dataset)
    if (!is.null(daf_obj)) {
        return(dafr::has_vector(daf_obj, "metacell", "projected_type"))
    }
    !is.null(get_mc_data(dataset, "query_md"))
}

has_corrected <- function(dataset) {
    # Check DAF metadata instead of loading the full corrected matrix
    daf_obj <- get_dataset_daf(dataset)
    if (!is.null(daf_obj)) {
        return(dafr::has_matrix(daf_obj, "metacell", "gene", "corrected_fraction"))
    }
    !is.null(get_mc_data(dataset, "mc_mat_corrected"))
}

has_atlas <- function(dataset) {
    # 'dataset' is unused but kept for API consistency with other has_* functions
    !is.null(get_atlas_daf())
}

any_has_atlas <- function() {
    has_atlas(dataset_names()[1])
}

calc_samp_mc_count <- function(dataset) {
    metadata <- get_mc_data(dataset, "cell_metadata")
    samp_mc_count <- metadata %>%
        filter(metacell != -1) %>%
        filter(!is.na(samp_id)) %>%
        count(samp_id, metacell) %>%
        spread(metacell, n, fill = 0) %>%
        column_to_rownames("samp_id") %>%
        as.matrix()
    mc_data <- mcv_get("mc_data")
    mc_data[[dataset]][["samp_mc_count"]] <- samp_mc_count
    mcv_set("mc_data", mc_data)
    return(samp_mc_count)
}

get_samp_mc_count <- function(dataset) {
    return(get_mc_data(dataset, "samp_mc_count") %||% calc_samp_mc_count(dataset))
}

calc_samp_mc_frac <- function(dataset) {
    samp_mc_count <- get_samp_mc_count(dataset)
    samp_mc_frac <- samp_mc_count / rowSums(samp_mc_count)
    mc_data <- mcv_get("mc_data")
    mc_data[[dataset]][["samp_mc_frac"]] <- samp_mc_frac
    mcv_set("mc_data", mc_data)
    return(samp_mc_frac)
}

get_samp_mc_frac <- function(dataset) {
    return(get_mc_data(dataset, "samp_mc_frac") %||% calc_samp_mc_frac(dataset))
}

calc_samp_metadata <- function(dataset) {
    metadata <- get_mc_data(dataset, "cell_metadata")
    if (is.null(metadata)) {
        return(NULL)
    }
    metadata <- metadata %>% select(-any_of(c("metacell", "cell_id", "outlier")))
    samp_columns <- metadata %>%
        group_by(samp_id) %>%
        summarise(across(everything(), n_distinct))
    samp_columns <- colnames(samp_columns)[purrr::map_lgl(colnames(samp_columns), ~ all(samp_columns[[.x]] == 1))]
    samp_md <- metadata %>%
        select(samp_id, samp_columns) %>%
        arrange(samp_id) %>%
        distinct(samp_id, .keep_all = TRUE)

    mc_data <- mcv_get("mc_data")
    mc_data[[dataset]][["samp_metadata"]] <- samp_md
    mcv_set("mc_data", mc_data)
    return(samp_md)
}

get_samp_metadata <- function(dataset) {
    if (!has_cell_metadata(dataset)) {
        return(NULL)
    }
    return(get_mc_data(dataset, "samp_metadata") %||% calc_samp_metadata(dataset))
}

calc_samples_list <- function(dataset) {
    if (!has_cell_metadata(dataset)) {
        return(NULL)
    }
    samp_md <- get_samp_metadata(dataset)
    samp_list <- sort(unique(samp_md$samp_id))
    mc_data <- mcv_get("mc_data")
    mc_data[[dataset]][["samp_list"]] <- samp_list
    mcv_set("mc_data", mc_data)
    return(samp_list)
}

get_samples_list <- function(dataset) {
    get_mc_data(dataset, "samp_list") %||% calc_samples_list(dataset)
}

# Remove standard identity fields from metadata field names
filter_metadata_field_names <- function(fields) {
    fields <- fields[!(fields %in% c("samp_id", "cell_id"))]
    fields <- fields[!grepl("samp_id: ", fields)]
    fields
}

dataset_metadata_fields <- function(dataset, atlas = FALSE) {
    cache_key <- if (atlas) "metadata_fields_atlas" else "metadata_fields"
    mc_data <- mcv_get("mc_data")
    if (!is.null(mc_data[[dataset]][[cache_key]])) {
        return(mc_data[[dataset]][[cache_key]])
    }

    daf_obj <- if (atlas) get_atlas_daf() else get_dataset_daf(dataset)
    if (!is.null(daf_obj)) {
        fields <- tryCatch(daf_obj["@ metacell : ?"], error = function(e) NULL)
        if (is.null(fields)) {
            fields <- tryCatch(dafr::vectors_set(daf_obj, "metacell"), error = function(e) character(0))
        }
        core_fields <- c("type", "x", "y", "u", "v", "umap_x", "umap_y", "total_UMIs", "n_cells", "n_cell")
        fields <- setdiff(fields, core_fields)
        fields <- filter_metadata_field_names(fields)
        fields <- fields[!grepl("^mcview_cache_", fields)]
        mc_data[[dataset]][[cache_key]] <- fields
        mcv_set("mc_data", mc_data)
        return(fields)
    }

    metadata <- get_mc_data(dataset, "metadata", atlas = atlas)
    if (is.null(metadata)) {
        return(c())
    }
    fields <- colnames(metadata)
    fields <- fields[fields != "metacell"]
    fields <- filter_metadata_field_names(fields)
    mc_data[[dataset]][[cache_key]] <- fields
    mcv_set("mc_data", mc_data)
    return(fields)
}

dataset_metadata_fields_numeric <- function(dataset, atlas = FALSE) {
    cache_key <- if (atlas) "metadata_fields_numeric_atlas" else "metadata_fields_numeric"
    mc_data <- mcv_get("mc_data")
    if (!is.null(mc_data[[dataset]][[cache_key]])) {
        return(mc_data[[dataset]][[cache_key]])
    }
    fields <- dataset_metadata_fields(dataset, atlas = atlas)
    df <- get_mc_data(dataset, "metadata", atlas = atlas)
    numeric_f <- purrr::map_lgl(fields, ~ is_numeric_field(df, .x))
    fields <- fields[numeric_f]
    mc_data[[dataset]][[cache_key]] <- fields
    mcv_set("mc_data", mc_data)
    return(fields)
}

dataset_metadata_fields_categorical <- function(dataset, atlas = FALSE) {
    cache_key <- if (atlas) "metadata_fields_categorical_atlas" else "metadata_fields_categorical"
    mc_data <- mcv_get("mc_data")
    if (!is.null(mc_data[[dataset]][[cache_key]])) {
        return(mc_data[[dataset]][[cache_key]])
    }

    fields <- dataset_metadata_fields(dataset, atlas = atlas)
    df <- get_mc_data(dataset, "metadata", atlas = atlas)

    if (is.null(df) || length(fields) == 0) {
        return(c())
    }

    # Filter to non-numeric fields first
    numeric_f <- purrr::map_lgl(fields, ~ is_numeric_field(df, .x))
    categorical_fields <- fields[!numeric_f]

    # Filter out fields with only one unique non-empty value (nothing to filter)
    categorical_fields <- purrr::keep(categorical_fields, function(field) {
        values <- df[[field]]
        clean_values <- values[!is.na(values) & !is.null(values) & values != ""]
        length(unique(clean_values)) > 1
    })

    mc_data[[dataset]][[cache_key]] <- categorical_fields
    mcv_set("mc_data", mc_data)
    return(categorical_fields)
}

dataset_cell_metadata_fields <- function(dataset, atlas = FALSE) {
    metadata <- get_mc_data(dataset, "cell_metadata", atlas = atlas)
    if (is.null(metadata)) {
        return(c())
    }
    fields <- colnames(metadata)
    return(filter_metadata_field_names(fields))
}

dataset_cell_metadata_fields_numeric <- function(dataset, atlas = FALSE) {
    fields <- dataset_cell_metadata_fields(dataset, atlas = atlas)
    df <- get_mc_data(dataset, "cell_metadata")
    numeric_f <- purrr::map_lgl(fields, ~ is_numeric_field(df, .x))
    return(fields[numeric_f])
}

dataset_cell_metadata_fields_categorical <- function(dataset, atlas = FALSE) {
    all_fields <- dataset_cell_metadata_fields(dataset, atlas = atlas)
    numeric_fields <- dataset_cell_metadata_fields_numeric(dataset, atlas = atlas)
    return(setdiff(all_fields, numeric_fields))
}

is_numeric_field <- function(df, field) {
    if (is.null(df)) {
        return(FALSE)
    }

    if (is.null(field)) {
        return(FALSE)
    }

    if (is.character(df[[field]])) {
        return(FALSE)
    }

    if (is.logical(df[[field]])) {
        return(FALSE)
    }

    if (is.numeric(df[[field]])) {
        return(TRUE)
    }

    return(FALSE)
}

get_metacell_ids <- function(dataset) {
    daf_obj <- get_dataset_daf(dataset)
    if (is.null(daf_obj)) {
        return(NULL)
    }
    dafr::axis_entries(daf_obj, "metacell")
}

has_gg_mc_top_cor <- function(dataset) {
    daf_obj <- get_dataset_daf(dataset)
    if (!is.null(daf_obj)) {
        return(dafr::has_axis(daf_obj, "gg_mc_top_cor") ||
            dafr::has_axis(daf_obj, "mcview_cache_gg_mc_top_cor"))
    }
    !is.null(get_mc_data(dataset, "gg_mc_top_cor"))
}

get_mc_sum <- function(dataset, atlas = FALSE) {
    mc_sum <- get_mc_data(dataset, "mc_sum", atlas = atlas)
    # convert_daf_mc_sum() already returns a named vector from dafr;
    # skip redundant names assignment to avoid COW materialization of jlview.
    if (is.null(names(mc_sum))) {
        daf_obj <- get_daf_for_query(dataset, atlas)
        if (!is.null(daf_obj)) {
            names(mc_sum) <- dafr::axis_entries(daf_obj, "metacell")
        }
    }
    return(mc_sum)
}

get_gene_qc <- function(dataset) {
    get_mc_data(dataset, "gene_qc")
}
