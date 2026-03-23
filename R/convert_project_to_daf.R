# convert_project_to_daf.R - Convert MCView project cache to DAF format

#' Convert MCView project to DAF format
#'
#' Converts all datasets in an MCView project cache directory to DAF format.
#'
#' @param project Path to MCView project
#' @param output Path to output directory for DAF data
#' @param format DAF format: "files" for directory-based, "h5" for HDF5
#' @param datasets Optional vector of dataset names to convert (default: all)
#' @param precompute_cache Whether to precompute cache after conversion
#' @param cache_daf_root Optional path to write auxiliary DAF caches (default: project/cache_daf)
#'
#' @export
convert_project_to_daf <- function(project,
                                   output,
                                   format = c("files", "h5"),
                                   datasets = NULL,
                                   precompute_cache = FALSE,
                                   cache_daf_root = NULL) {
    format <- match.arg(format)

    cache_dir <- project_cache_dir(project)
    if (!fs::dir_exists(cache_dir)) {
        cli_abort("Project cache directory not found: {.path {cache_dir}}")
    }

    if (is.null(datasets)) {
        datasets <- list.dirs(cache_dir, recursive = FALSE, full.names = FALSE)
        datasets <- datasets[datasets != "atlas"] # Skip atlas subdirectory
    }

    about_markdown <- read_about_markdown(project)

    fs::dir_create(output)
    cache_daf_root <- cache_daf_root %||% fs::path(project, "cache_daf")

    for (dataset in datasets) {
        cli_alert_info("Converting dataset: {dataset}")
        dataset_cache <- fs::path(cache_dir, dataset)

        if (!fs::dir_exists(dataset_cache)) {
            cli_warn("Dataset cache not found: {.path {dataset_cache}}")
            next
        }

        output_path <- if (format == "h5") {
            fs::path(output, paste0(dataset, ".h5"))
        } else {
            fs::path(output, dataset)
        }

        convert_dataset_to_daf(
            dataset_cache,
            output_path,
            format,
            about_markdown = about_markdown,
            precompute_cache = precompute_cache,
            cache_daf_root = cache_daf_root
        )

        # Convert atlas if present
        atlas_cache <- fs::path(cache_dir, dataset, "atlas")
        if (fs::dir_exists(atlas_cache)) {
            cli_alert_info("Converting atlas for dataset: {dataset}")
            atlas_output <- if (format == "h5") {
                fs::path(output, paste0(dataset, "_atlas.h5"))
            } else {
                fs::path(output, paste0(dataset, "_atlas"))
            }
            convert_dataset_to_daf(
                atlas_cache,
                atlas_output,
                format,
                about_markdown = about_markdown,
                precompute_cache = precompute_cache,
                cache_daf_root = cache_daf_root
            )
        }
    }

    cli_alert_success("Conversion complete. Output: {.path {output}}")
}

#' Convert a single dataset cache to DAF
#'
#' @param cache_dir Path to dataset cache directory
#' @param output_path Path to output DAF
#' @param format DAF format: "files" or "h5"
#'
#' @noRd
convert_dataset_to_daf <- function(cache_dir,
                                   output_path,
                                   format,
                                   about_markdown = NULL,
                                   precompute_cache = FALSE,
                                   cache_daf_root = NULL) {
    # Read core data from cache
    mc_mat <- read_cache_file(cache_dir, "mc_mat")
    mc_sum <- read_cache_file(cache_dir, "mc_sum")
    mc2d <- read_cache_file(cache_dir, "mc2d")
    metacell_types <- read_cache_file(cache_dir, "metacell_types")
    cell_type_colors <- read_cache_file(cache_dir, "cell_type_colors")

    if (is.null(mc_mat)) {
        cli_abort("Required file mc_mat not found in {.path {cache_dir}}")
    }

    # Get axes (ensure all are character vectors for DAF)
    gene_names <- as.character(rownames(mc_mat))
    metacell_names <- as.character(colnames(mc_mat))

    # Get type names from colors or metacell_types
    if (!is.null(cell_type_colors)) {
        type_names <- as.character(cell_type_colors$cell_type)
    } else if (!is.null(metacell_types)) {
        type_names <- as.character(unique(metacell_types$cell_type))
    } else {
        cli_abort("No cell type information found")
    }

    # Remove existing output if present
    if (file.exists(output_path) || dir.exists(output_path)) {
        cli_alert_warning("Removing existing DAF at {.path {output_path}}")
        unlink(output_path, recursive = TRUE, force = TRUE)
        # Wait for filesystem to sync (NFS can be slow)
        Sys.sleep(1)
        # Verify it's gone, retry if needed
        retries <- 0
        while ((file.exists(output_path) || dir.exists(output_path)) && retries < 5) {
            Sys.sleep(1)
            unlink(output_path, recursive = TRUE, force = TRUE)
            retries <- retries + 1
        }
        if (file.exists(output_path) || dir.exists(output_path)) {
            cli_abort("Could not remove existing DAF at {.path {output_path}}")
        }
    }

    # Create DAF object
    daf <- if (format == "h5") {
        dafr::h5df(output_path, mode = "w+")
    } else {
        fs::dir_create(output_path)
        dafr::files_daf(output_path, mode = "w+")
    }

    # Add axes
    daf <- dafr::add_axis(daf, "gene", gene_names)
    daf <- dafr::add_axis(daf, "metacell", metacell_names)
    daf <- dafr::add_axis(daf, "type", type_names)

    if (!is.null(cache_daf_root)) {
        daf <- dafr::set_scalar(daf, "mcview_cache_daf_root", cache_daf_root)
    }

    # Add core matrix (transpose to metacell x gene)
    mc_mat_transposed <- Matrix::t(mc_mat)
    daf <- dafr::set_matrix(daf, "metacell", "gene", "UMIs", mc_mat_transposed)

    # Add mc_sum vector
    if (!is.null(mc_sum)) {
        daf <- dafr::set_vector(daf, "metacell", "total_UMIs", as.numeric(mc_sum))
    } else {
        # Compute from matrix
        daf <- dafr::set_vector(daf, "metacell", "total_UMIs", Matrix::rowSums(mc_mat_transposed))
    }

    # Add 2D coordinates
    if (!is.null(mc2d)) {
        # Handle different mc2d formats
        if (is.data.frame(mc2d)) {
            # Tibble format with metacell, x, y columns
            mc2d_ordered <- mc2d[match(metacell_names, as.character(mc2d$metacell)), ]
            x_coords <- as.numeric(mc2d_ordered$x)
            y_coords <- as.numeric(mc2d_ordered$y)
        } else if (is.list(mc2d) && "mc_x" %in% names(mc2d)) {
            # List format with mc_id, mc_x, mc_y
            idx <- match(metacell_names, as.character(mc2d$mc_id))
            x_coords <- as.numeric(mc2d$mc_x[idx])
            y_coords <- as.numeric(mc2d$mc_y[idx])
        } else {
            cli_warn("Unknown mc2d format, skipping coordinates")
            x_coords <- NULL
            y_coords <- NULL
        }

        if (!is.null(x_coords) && !is.null(y_coords)) {
            daf <- dafr::set_vector(daf, "metacell", "x", x_coords)
            daf <- dafr::set_vector(daf, "metacell", "y", y_coords)
        }
    }

    # Add metacell types
    if (!is.null(metacell_types)) {
        mc_types_ordered <- metacell_types[match(metacell_names, as.character(metacell_types$metacell)), ]
        daf <- dafr::set_vector(daf, "metacell", "type", as.character(mc_types_ordered$cell_type))

        # Add n_cell if present
        if ("n_cell" %in% names(metacell_types)) {
            daf <- dafr::set_vector(daf, "metacell", "n_cell", as.numeric(mc_types_ordered$n_cell))
        }
        if ("top1_gene" %in% names(metacell_types)) {
            vec <- as.character(mc_types_ordered$top1_gene)
            if (!any(is.na(vec))) {
                daf <- dafr::set_vector(daf, "metacell", "top1_gene", vec)
            }
        }
        if ("top2_gene" %in% names(metacell_types)) {
            vec <- as.character(mc_types_ordered$top2_gene)
            if (!any(is.na(vec))) {
                daf <- dafr::set_vector(daf, "metacell", "top2_gene", vec)
            }
        }
        if ("top1_lfp" %in% names(metacell_types)) {
            vec <- as.numeric(mc_types_ordered$top1_lfp)
            if (!any(is.na(vec))) {
                daf <- dafr::set_vector(daf, "metacell", "top1_lfp", vec)
            }
        }
        if ("top2_lfp" %in% names(metacell_types)) {
            vec <- as.numeric(mc_types_ordered$top2_lfp)
            if (!any(is.na(vec))) {
                daf <- dafr::set_vector(daf, "metacell", "top2_lfp", vec)
            }
        }
        if ("mc_col" %in% names(metacell_types)) {
            vec <- as.character(mc_types_ordered$mc_col)
            if (!any(is.na(vec))) {
                daf <- dafr::set_vector(daf, "metacell", "mc_col", vec)
            }
        }
    }

    # Add cell type colors
    if (!is.null(cell_type_colors)) {
        colors_ordered <- cell_type_colors[match(type_names, as.character(cell_type_colors$cell_type)), ]
        daf <- dafr::set_vector(daf, "type", "color", as.character(colors_ordered$color))
    }

    # Add optional gene vectors
    add_gene_flags(daf, cache_dir, gene_names)

    # Add optional matrices
    add_optional_matrices(daf, cache_dir, gene_names, metacell_names)

    # Add metadata
    add_metadata_vectors(daf, cache_dir, metacell_names)

    # Add projection data (for atlas-projected datasets)
    add_projection_data(daf, cache_dir, metacell_names)

    # Add gene QC vectors
    add_gene_qc_vectors(daf, cache_dir, gene_names)

    # Add zero-fold gene statistics
    add_gene_zero_fold(daf, cache_dir)

    # Add gene correlations
    add_gene_correlations(daf, cache_dir)

    # Add metacell graphs
    add_metacell_graphs(daf, cache_dir)

    # Add QC stats
    add_qc_stats(daf, cache_dir)

    # Add about markdown
    add_about_markdown(daf, about_markdown)

    if (isTRUE(precompute_cache)) {
        precompute_daf_cache(daf)
    }

    cli_alert_success("Created DAF at {.path {output_path}}")
}

#' Read file from cache directory
#' @noRd
read_cache_file <- function(cache_dir, name) {
    # Try .qs format
    qs_path <- fs::path(cache_dir, paste0(name, ".qs"))
    if (fs::file_exists(qs_path)) {
        return(qs::qread(qs_path))
    }

    # Try .tsv format
    tsv_path <- fs::path(cache_dir, paste0(name, ".tsv"))
    if (fs::file_exists(tsv_path)) {
        return(data.table::fread(tsv_path, data.table = FALSE) %>% tibble::as_tibble())
    }

    # Try .csv format
    csv_path <- fs::path(cache_dir, paste0(name, ".csv"))
    if (fs::file_exists(csv_path)) {
        return(data.table::fread(csv_path, data.table = FALSE) %>% tibble::as_tibble())
    }

    return(NULL)
}

#' Add gene flag vectors to DAF
#' @noRd
add_gene_flags <- function(daf, cache_dir, gene_names) {
    # Lateral genes
    lateral <- read_cache_file(cache_dir, "lateral_genes")
    if (!is.null(lateral) && length(lateral) > 0) {
        is_lateral <- gene_names %in% as.character(lateral)
        daf <- dafr::set_vector(daf, "gene", "is_lateral", as.logical(is_lateral))
    }

    # Noisy genes
    noisy <- read_cache_file(cache_dir, "noisy_genes")
    if (!is.null(noisy) && length(noisy) > 0) {
        is_noisy <- gene_names %in% as.character(noisy)
        daf <- dafr::set_vector(daf, "gene", "is_noisy", as.logical(is_noisy))
    }

    # Marker genes
    marker <- read_cache_file(cache_dir, "marker_genes")
    if (!is.null(marker) && is.data.frame(marker) && nrow(marker) > 0) {
        is_marker <- gene_names %in% as.character(marker$gene)
        daf <- dafr::set_vector(daf, "gene", "is_marker", as.logical(is_marker))
    }

    # Gene modules - only add if all genes have modules (DAF doesn't allow NAs)
    gene_modules <- read_cache_file(cache_dir, "gene_modules")
    if (!is.null(gene_modules) && is.data.frame(gene_modules) && nrow(gene_modules) > 0) {
        # Check if we have modules for all genes
        idx <- match(gene_names, as.character(gene_modules$gene))
        if (!any(is.na(idx))) {
            module_vec <- as.character(gene_modules$module[idx])
            if (!any(is.na(module_vec))) {
                daf <- dafr::set_vector(daf, "gene", "module", module_vec)
            } else {
                cli_alert_info("Skipping gene modules (contains NA module values)")
            }
        } else {
            cli_alert_info("Skipping gene modules (not all genes have module assignments)")
        }
    }

    invisible(daf)
}

#' Add optional matrices to DAF
#' @noRd
add_optional_matrices <- function(daf, cache_dir, gene_names, metacell_names) {
    # Helper: reindex a matrix to match target row/col names, filling missing with 0
    reindex_matrix <- function(mat, target_rows, target_cols) {
        # Get available rows/cols
        avail_rows <- as.character(rownames(mat))
        avail_cols <- as.character(colnames(mat))

        # Metacells must all be present (columns)
        col_idx <- match(target_cols, avail_cols)
        if (any(is.na(col_idx))) {
            return(NULL)
        }

        # For genes (rows), allow partial overlap - fill missing with 0
        row_idx <- match(target_rows, avail_rows)
        has_all_rows <- !any(is.na(row_idx))

        if (has_all_rows) {
            return(as.matrix(mat[row_idx, col_idx, drop = FALSE]))
        }

        # Partial overlap: build dense matrix, fill available rows
        result <- matrix(0, nrow = length(target_rows), ncol = length(target_cols),
                         dimnames = list(target_rows, target_cols))
        present <- !is.na(row_idx)
        result[present, ] <- as.matrix(mat[row_idx[present], col_idx, drop = FALSE])
        result
    }

    # Inner fold matrix
    inner_fold <- read_cache_file(cache_dir, "inner_fold_mat")
    if (!is.null(inner_fold) && !is.null(dim(inner_fold))) {
        mat <- reindex_matrix(inner_fold, gene_names, metacell_names)
        if (!is.null(mat)) {
            daf <- dafr::set_matrix(daf, "gene", "metacell", "inner_fold", mat)
        } else {
            cli_alert_info("Skipping inner_fold matrix (metacell mismatch)")
        }
    }

    # Inner stdev matrix
    inner_stdev <- read_cache_file(cache_dir, "inner_stdev_mat")
    if (!is.null(inner_stdev) && !is.null(dim(inner_stdev))) {
        mat <- reindex_matrix(inner_stdev, gene_names, metacell_names)
        if (!is.null(mat)) {
            daf <- dafr::set_matrix(daf, "gene", "metacell", "inner_stdev_log", mat)
        } else {
            cli_alert_info("Skipping inner_stdev_log matrix (metacell mismatch)")
        }
    }

    # Projected fold matrix
    projected_fold <- read_cache_file(cache_dir, "projected_fold")
    if (!is.null(projected_fold) && !is.null(dim(projected_fold))) {
        mat <- reindex_matrix(projected_fold, gene_names, metacell_names)
        if (!is.null(mat)) {
            daf <- dafr::set_matrix(daf, "gene", "metacell", "projected_fold", mat)
        } else {
            cli_alert_info("Skipping projected_fold matrix (metacell mismatch)")
        }
    }

    # Corrected UMI matrix (mc_mat_corrected)
    mc_mat_corrected <- read_cache_file(cache_dir, "mc_mat_corrected")
    if (!is.null(mc_mat_corrected) && !is.null(dim(mc_mat_corrected))) {
        row_idx <- match(gene_names, as.character(rownames(mc_mat_corrected)))
        col_idx <- match(metacell_names, as.character(colnames(mc_mat_corrected)))
        # Skip if indices have NAs (genes/metacells don't match)
        if (!any(is.na(row_idx)) && !any(is.na(col_idx))) {
            mc_mat_corrected <- mc_mat_corrected[row_idx, col_idx, drop = FALSE]
            # Store as corrected_fraction (normalized by total UMIs)
            mc_sum <- read_cache_file(cache_dir, "mc_sum")
            if (!is.null(mc_sum)) {
                corrected_frac <- t(t(mc_mat_corrected) / mc_sum[colnames(mc_mat_corrected)])
                # Check for NAs - skip if present
                if (!any(is.na(corrected_frac))) {
                    daf <- dafr::set_matrix(daf, "metacell", "gene", "corrected_fraction", t(as.matrix(corrected_frac)))
                } else {
                    cli_alert_info("Skipping corrected_fraction matrix (contains NA values)")
                }
            }
        } else {
            cli_alert_info("Skipping corrected_fraction matrix (gene/metacell mismatch)")
        }
    }

    # Projected (expected) matrix
    projected_mat <- read_cache_file(cache_dir, "projected_mat")
    if (!is.null(projected_mat) && !is.null(dim(projected_mat))) {
        row_idx <- match(gene_names, as.character(rownames(projected_mat)))
        col_idx <- match(metacell_names, as.character(colnames(projected_mat)))
        # Skip if indices have NAs
        if (!any(is.na(row_idx)) && !any(is.na(col_idx))) {
            projected_mat <- projected_mat[row_idx, col_idx, drop = FALSE]
            # Store as projected_fraction (normalized by total UMIs)
            mc_sum <- mc_sum %||% read_cache_file(cache_dir, "mc_sum")
            if (!is.null(mc_sum)) {
                projected_frac <- t(t(projected_mat) / mc_sum[colnames(projected_mat)])
                # Check for NAs
                if (!any(is.na(projected_frac))) {
                    daf <- dafr::set_matrix(daf, "metacell", "gene", "projected_fraction", t(as.matrix(projected_frac)))
                } else {
                    cli_alert_info("Skipping projected_fraction matrix (contains NA values)")
                }
            }
        } else {
            cli_alert_info("Skipping projected_fraction matrix (gene/metacell mismatch)")
        }
    }

    invisible(daf)
}

#' Add projection weights to DAF
#' @noRd
add_projection_data <- function(daf, cache_dir, metacell_names) {
    # Projection weights
    proj_weights <- read_cache_file(cache_dir, "proj_weights")
    if (!is.null(proj_weights) && is.data.frame(proj_weights) && nrow(proj_weights) > 0) {
        # proj_weights has columns: query, atlas, weight
        # We need to add this as an edge matrix if both axes exist
        # For now, store as a scalar JSON or skip if complex
        cli_alert_info("Projection weights found - storing as serialized data")
        # Store as JSON scalar for now (DAF edge support is complex)
        proj_weights_json <- jsonlite::toJSON(proj_weights, auto_unbox = TRUE)
        daf <- dafr::set_scalar(daf, "projection_weights_json", as.character(proj_weights_json))
    }

    # Projected metacell types (skip if already set from metadata)
    projected_types <- read_cache_file(cache_dir, "projected_metacell_types")
    if (!is.null(projected_types) && is.data.frame(projected_types) && nrow(projected_types) > 0) {
        if (!dafr::has_vector(daf, "metacell", "projected_type")) {
            # Store projected_type as a metacell vector
            idx <- match(metacell_names, as.character(projected_types$metacell))
            if (!any(is.na(idx))) {
                proj_type_vec <- as.character(projected_types$cell_type[idx])
                if (!any(is.na(proj_type_vec))) {
                    daf <- dafr::set_vector(daf, "metacell", "projected_type", proj_type_vec)
                }
            }
        }
    }

    # Query atlas cell type fractions
    query_fracs <- read_cache_file(cache_dir, "query_atlas_cell_type_fracs")
    if (!is.null(query_fracs) && is.data.frame(query_fracs) && nrow(query_fracs) > 0) {
        # Store as JSON scalar (complex matrix with metacell x cell_type)
        cli_alert_info("Query atlas cell type fractions found - storing as serialized data")
        query_fracs_json <- jsonlite::toJSON(query_fracs, auto_unbox = TRUE)
        daf <- dafr::set_scalar(daf, "query_cell_type_fracs_json", as.character(query_fracs_json))
    }

    invisible(daf)
}

#' Sanitize column name for DAF
#' @noRd
sanitize_daf_name <- function(name) {
    # Replace / and \ with underscore (they're invalid in file paths)
    # Also replace other problematic characters
    gsub("[/\\\\:*?\"<>|]", "_", name)
}

#' Coerce a vector to a DAF-compatible type
#'
#' Converts factors to character, and ensures numeric/logical/character types
#' are explicit. Returns NULL for list columns.
#'
#' @param vec A vector from a data frame column
#' @return Coerced vector, or NULL if the column type is unsupported (e.g. list)
#' @noRd
coerce_vec_for_daf <- function(vec) {
    if (is.list(vec)) {
        return(NULL)
    }
    if (is.factor(vec)) {
        return(as.character(vec))
    }
    if (is.numeric(vec)) {
        return(as.numeric(vec))
    }
    if (is.logical(vec)) {
        return(as.logical(vec))
    }
    as.character(vec)
}

#' Set DAF vectors from a data frame
#'
#' Iterates over columns of a data frame, performing type coercion and writing
#' each as a vector on the given axis. Skips columns that already exist in the
#' DAF, columns containing NAs, and list columns.
#'
#' @param daf DAF object (writable)
#' @param axis Axis name (e.g. "metacell", "gene")
#' @param df Data frame whose columns will become vectors
#' @param id_col Name of the ID column to skip (e.g. "metacell", "gene")
#' @param skip_existing If TRUE, skip columns that already exist as vectors on the axis
#'
#' @return The DAF object (invisibly), potentially modified
#' @noRd
set_daf_vectors_from_df <- function(daf, axis, df, id_col, skip_existing = TRUE) {
    for (col in setdiff(names(df), id_col)) {
        safe_col <- sanitize_daf_name(col)

        if (skip_existing && dafr::has_vector(daf, axis, safe_col)) {
            next
        }

        vec <- coerce_vec_for_daf(df[[col]])
        if (is.null(vec)) {
            next
        }

        if (!any(is.na(vec))) {
            daf <- dafr::set_vector(daf, axis, safe_col, vec)
        }
    }

    invisible(daf)
}

#' Add metadata vectors to DAF
#' @noRd
add_metadata_vectors <- function(daf, cache_dir, metacell_names) {
    # Read metadata
    metadata <- read_cache_file(cache_dir, "metadata")
    if (!is.null(metadata) && is.data.frame(metadata) && nrow(metadata) > 0) {
        metadata <- metadata[match(metacell_names, as.character(metadata$metacell)), ]
        daf <- set_daf_vectors_from_df(daf, "metacell", metadata, "metacell")
    }

    # Read QC metadata
    mc_qc <- read_cache_file(cache_dir, "mc_qc_metadata")
    if (!is.null(mc_qc) && is.data.frame(mc_qc) && nrow(mc_qc) > 0) {
        mc_qc <- mc_qc[match(metacell_names, as.character(mc_qc$metacell)), ]
        daf <- set_daf_vectors_from_df(daf, "metacell", mc_qc, "metacell")
    }

    invisible(daf)
}

#' Add gene QC vectors to DAF
#' @noRd
add_gene_qc_vectors <- function(daf, cache_dir, gene_names) {
    gene_qc <- read_cache_file(cache_dir, "gene_qc")
    if (!is.null(gene_qc) && is.data.frame(gene_qc) && nrow(gene_qc) > 0) {
        gene_qc <- gene_qc[match(gene_names, as.character(gene_qc$gene)), ]
        daf <- set_daf_vectors_from_df(daf, "gene", gene_qc, "gene")
    }

    invisible(daf)
}

#' Add gene zero-fold statistics to DAF
#' @noRd
add_gene_zero_fold <- function(daf, cache_dir) {
    gene_zero_fold <- read_cache_file(cache_dir, "gene_zero_fold")
    if (!is.null(gene_zero_fold) && is.data.frame(gene_zero_fold) && nrow(gene_zero_fold) > 0) {
        required_cols <- c("gene", "metacell", "zero_fold", "avg", "obs", "exp", "type")
        if (all(required_cols %in% colnames(gene_zero_fold))) {
            gene_zero_fold <- gene_zero_fold %>%
                filter(!is.na(gene), !is.na(metacell), !is.na(zero_fold), !is.na(avg), !is.na(obs), !is.na(exp), !is.na(type))

            if (nrow(gene_zero_fold) > 0) {
                axis_name <- "gene_zero_fold"
                if (!dafr::has_axis(daf, axis_name)) {
                    daf <- dafr::add_axis(daf, axis_name, as.character(seq_len(nrow(gene_zero_fold))))
                }

                daf <- dafr::set_vector(daf, axis_name, "gene", as.character(gene_zero_fold$gene))
                daf <- dafr::set_vector(daf, axis_name, "metacell", as.character(gene_zero_fold$metacell))
                daf <- dafr::set_vector(daf, axis_name, "zero_fold", as.numeric(gene_zero_fold$zero_fold))
                daf <- dafr::set_vector(daf, axis_name, "avg", as.numeric(gene_zero_fold$avg))
                daf <- dafr::set_vector(daf, axis_name, "obs", as.numeric(gene_zero_fold$obs))
                daf <- dafr::set_vector(daf, axis_name, "exp", as.numeric(gene_zero_fold$exp))
                daf <- dafr::set_vector(daf, axis_name, "type", as.character(gene_zero_fold$type))
            }
        }
    }

    invisible(daf)
}

#' Add gene correlations to DAF
#' @noRd
add_gene_correlations <- function(daf, cache_dir) {
    gg_mc_top_cor <- read_cache_file(cache_dir, "gg_mc_top_cor")
    if (!is.null(gg_mc_top_cor) && is.data.frame(gg_mc_top_cor) && nrow(gg_mc_top_cor) > 0) {
        required_cols <- c("gene1", "gene2", "cor", "type")
        if (all(required_cols %in% colnames(gg_mc_top_cor))) {
            gg_mc_top_cor <- gg_mc_top_cor %>%
                filter(!is.na(gene1), !is.na(gene2), !is.na(cor), !is.na(type))

            if (nrow(gg_mc_top_cor) > 0) {
                axis_name <- "gg_mc_top_cor"
                if (!dafr::has_axis(daf, axis_name)) {
                    daf <- dafr::add_axis(daf, axis_name, as.character(seq_len(nrow(gg_mc_top_cor))))
                }

                daf <- dafr::set_vector(daf, axis_name, "gene1", as.character(gg_mc_top_cor$gene1))
                daf <- dafr::set_vector(daf, axis_name, "gene2", as.character(gg_mc_top_cor$gene2))
                daf <- dafr::set_vector(daf, axis_name, "cor", as.numeric(gg_mc_top_cor$cor))
                daf <- dafr::set_vector(daf, axis_name, "type", as.character(gg_mc_top_cor$type))
            }
        }
    }

    invisible(daf)
}

#' Add metacell graph(s) to DAF
#' @noRd
add_metacell_graphs <- function(daf, cache_dir) {
    graphs <- list()

    mc2d <- read_cache_file(cache_dir, "mc2d")
    if (is.list(mc2d) && !is.null(mc2d$graph) && is.data.frame(mc2d$graph) && nrow(mc2d$graph) > 0) {
        graphs[["metacell"]] <- mc2d$graph
    }

    metacell_graphs <- read_cache_file(cache_dir, "metacell_graphs")
    if (!is.null(metacell_graphs) && is.list(metacell_graphs) && length(metacell_graphs) > 0) {
        for (graph_name in names(metacell_graphs)) {
            if (is.null(graph_name) || graph_name == "" || graph_name == "metacell") {
                next
            }
            if (is.data.frame(metacell_graphs[[graph_name]]) && nrow(metacell_graphs[[graph_name]]) > 0) {
                graphs[[graph_name]] <- metacell_graphs[[graph_name]]
            }
        }
    }

    if (length(graphs) == 0) {
        return(invisible(daf))
    }

    graph_df <- purrr::imap_dfr(graphs, function(graph, graph_name) {
        if (all(c("mc1", "mc2", "weight") %in% colnames(graph))) {
            tibble(
                graph_name = graph_name,
                from = as.character(graph$mc1),
                to = as.character(graph$mc2),
                weight = as.numeric(graph$weight)
            )
        } else if (all(c("from", "to", "weight") %in% colnames(graph))) {
            tibble(
                graph_name = graph_name,
                from = as.character(graph$from),
                to = as.character(graph$to),
                weight = as.numeric(graph$weight)
            )
        } else {
            tibble(graph_name = character(), from = character(), to = character(), weight = numeric())
        }
    }) %>%
        filter(!is.na(graph_name), !is.na(from), !is.na(to), !is.na(weight))

    if (nrow(graph_df) == 0) {
        return(invisible(daf))
    }

    axis_name <- "metacell_graph"
    if (!dafr::has_axis(daf, axis_name)) {
        daf <- dafr::add_axis(daf, axis_name, as.character(seq_len(nrow(graph_df))))
    }

    daf <- dafr::set_vector(daf, axis_name, "graph_name", as.character(graph_df$graph_name))
    daf <- dafr::set_vector(daf, axis_name, "from", as.character(graph_df$from))
    daf <- dafr::set_vector(daf, axis_name, "to", as.character(graph_df$to))
    daf <- dafr::set_vector(daf, axis_name, "weight", as.numeric(graph_df$weight))

    invisible(daf)
}

#' Add QC stats scalars to DAF
#' @noRd
add_qc_stats <- function(daf, cache_dir) {
    qc_stats <- read_cache_file(cache_dir, "qc_stats")
    if (!is.null(qc_stats) && is.list(qc_stats) && length(qc_stats) > 0) {
        for (name in names(qc_stats)) {
            if (length(qc_stats[[name]]) == 1 && !is.na(qc_stats[[name]])) {
                scalar_name <- glue("qc_stats_{name}")
                daf <- dafr::set_scalar(daf, scalar_name, qc_stats[[name]])
            }
        }
    }

    invisible(daf)
}

#' Add about markdown to DAF
#' @noRd
add_about_markdown <- function(daf, about_markdown) {
    if (is.null(about_markdown) || !nzchar(about_markdown)) {
        return(invisible(daf))
    }

    daf <- dafr::set_scalar(daf, "mcview_about_markdown", about_markdown)
    invisible(daf)
}

#' Read about markdown from project config
#' @noRd
read_about_markdown <- function(project) {
    about_file <- project_about_file(project)
    if (!fs::file_exists(about_file)) {
        return(NULL)
    }

    raw_text <- paste(readLines(about_file, warn = FALSE), collapse = "\n")
    clean_about_markdown(raw_text)
}

#' Get project cache directory
#' @noRd
project_cache_dir <- function(project) {
    fs::path(project, "cache")
}
