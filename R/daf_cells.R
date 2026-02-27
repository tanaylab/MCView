# daf_cells.R - Cell-level data access for MCView
# Provides access to cell-level DAF data (cells_clean) and
# flexible group-based analysis (pseudobulk, composition, DE, QC).

# ==============================================================================
# Cells DAF Management
# ==============================================================================

#' Get the cells DAF object for a dataset
#'
#' Returns a chained DAF that combines the metacells DAF with the cells DAF,
#' providing unified access to cell-level metadata, cell-gene UMIs, and
#' cell-metacell mapping.
#'
#' @param dataset Dataset name
#' @return Chained Daf object or NULL if no cells DAF is configured
#' @export
get_cells_daf <- function(dataset) {
    mc_data <- mcv_get("mc_data")
    if (is.null(mc_data) || is.null(mc_data[[dataset]])) {
        return(NULL)
    }

    # Return cached chained DAF if available
    if (!is.null(mc_data[[dataset]][["chained_cells_daf"]])) {
        return(mc_data[[dataset]][["chained_cells_daf"]])
    }

    # Check if a cells DAF path was configured
    cells_daf_path <- mc_data[[dataset]][["cells_daf_path"]]
    if (is.null(cells_daf_path) || !fs::dir_exists(cells_daf_path)) {
        return(NULL)
    }

    # Load the cells DAF and chain it with the metacells DAF
    tryCatch(
        {
            cells_daf <- dafr::files_daf(cells_daf_path, mode = "r")
            metacells_daf <- mc_data[[dataset]][["daf_obj"]]
            if (is.null(metacells_daf)) {
                metacells_daf <- mc_data[[dataset]][["base_daf"]]
            }
            if (is.null(metacells_daf)) {
                cli_warn("No metacells DAF found for dataset '{dataset}'; returning cells DAF only")
                mc_data[[dataset]][["cells_daf"]] <- cells_daf
                mc_data[[dataset]][["chained_cells_daf"]] <- cells_daf
                mcv_set("mc_data", mc_data)
                return(cells_daf)
            }

            chained <- dafr::chain_reader(
                list(metacells_daf, cells_daf),
                name = paste0(dataset, ".chained")
            )

            # Cache both the raw cells DAF and the chained version
            mc_data[[dataset]][["cells_daf"]] <- cells_daf
            mc_data[[dataset]][["chained_cells_daf"]] <- chained
            mcv_set("mc_data", mc_data)

            return(chained)
        },
        error = function(e) {
            cli_warn("Failed to load cells DAF from '{cells_daf_path}': {conditionMessage(e)}")
            return(NULL)
        }
    )
}

#' Set the cells DAF path for a dataset
#'
#' Configures the path to a cells-level DAF (e.g., cells_clean) so that
#' get_cells_daf() can load and chain it with the metacells DAF.
#'
#' @param dataset Dataset name
#' @param cells_daf_path Path to the cells DAF directory
#' @return TRUE invisibly on success
#' @export
set_cells_daf <- function(dataset, cells_daf_path) {
    mc_data <- mcv_get("mc_data")
    if (is.null(mc_data) || is.null(mc_data[[dataset]])) {
        cli_abort("Dataset '{dataset}' not found in mc_data. Initialize the metacells DAF first.")
    }

    if (!fs::dir_exists(cells_daf_path)) {
        cli_abort("Cells DAF path does not exist: {cells_daf_path}")
    }

    # Clear any previously cached cells DAF
    mc_data[[dataset]][["cells_daf_path"]] <- cells_daf_path
    mc_data[[dataset]][["cells_daf"]] <- NULL
    mc_data[[dataset]][["chained_cells_daf"]] <- NULL
    mcv_set("mc_data", mc_data)

    cli_alert_success("Cells DAF path set for '{dataset}': {cells_daf_path}")
    invisible(TRUE)
}

#' Check if cell-level gene UMIs are available
#'
#' @param dataset Dataset name
#' @return TRUE if a cells DAF with cell x gene UMIs matrix is available
#' @export
has_cell_gene_umis <- function(dataset) {
    cells_daf <- get_cells_daf(dataset)
    if (is.null(cells_daf)) {
        return(FALSE)
    }
    tryCatch(
        dafr::has_matrix(cells_daf, "cell", "gene", "UMIs"),
        error = function(e) FALSE
    )
}

# ==============================================================================
# Cell Metadata Access
# ==============================================================================

#' Get cell-to-value mapping for a categorical field
#'
#' Returns a named character vector mapping cell IDs to their values for
#' the specified field.
#'
#' @param dataset Dataset name
#' @param field_name Name of the cell-level vector property
#' @return Named character vector (cell_id -> field_value), or NULL if unavailable
#' @export
get_cell_field_map <- function(dataset, field_name) {
    cells_daf <- get_cells_daf(dataset)
    if (is.null(cells_daf)) {
        return(NULL)
    }

    if (!dafr::has_vector(cells_daf, "cell", field_name)) {
        return(NULL)
    }

    vec <- dafr::get_vector(cells_daf, "cell", field_name)
    result <- as.character(vec)
    if (is.null(names(result))) {
        names(result) <- dafr::axis_entries(cells_daf, "cell")
    }
    result
}

#' List available categorical fields for grouping cells
#'
#' Queries the cells DAF for all cell-level vector properties and filters
#' to categorical (string) fields with reasonable cardinality (< 1000 unique values).
#'
#' @param dataset Dataset name
#' @param max_cardinality Maximum number of unique values to consider a field
#'   suitable for grouping (default: 1000)
#' @return Character vector of field names, or character(0) if none available
#' @export
get_cell_grouping_fields <- function(dataset, max_cardinality = 1000) {
    cells_daf <- get_cells_daf(dataset)
    if (is.null(cells_daf)) {
        return(character(0))
    }

    # Get all vector properties on the cell axis
    all_vectors <- tryCatch(
        dafr::vectors_set(cells_daf, "cell"),
        error = function(e) character(0)
    )

    if (length(all_vectors) == 0) {
        return(character(0))
    }

    # Exclude internal/identity fields
    exclude_fields <- c("metacell", "name")
    candidates <- setdiff(all_vectors, exclude_fields)

    # Filter to string/categorical fields with reasonable cardinality
    grouping_fields <- character(0)
    for (field in candidates) {
        tryCatch(
            {
                vec <- dafr::get_vector(cells_daf, "cell", field)
                if (is.character(vec)) {
                    n_unique <- length(unique(vec))
                    if (n_unique > 1 && n_unique <= max_cardinality) {
                        grouping_fields <- c(grouping_fields, field)
                    }
                }
            },
            error = function(e) NULL
        )
    }

    return(grouping_fields)
}

#' Get cell-to-metacell mapping
#'
#' Returns a named character vector mapping cell IDs to their assigned metacell.
#'
#' @param dataset Dataset name
#' @return Named character vector (cell_id -> metacell_id), or NULL if unavailable
#' @export
get_cell_metacell_map <- function(dataset) {
    cells_daf <- get_cells_daf(dataset)
    if (is.null(cells_daf)) {
        # Fallback: try the metacells DAF directly
        daf_obj <- get_dataset_daf(dataset)
        if (!is.null(daf_obj) && dafr::has_axis(daf_obj, "cell") &&
            dafr::has_vector(daf_obj, "cell", "metacell")) {
            vec <- dafr::get_vector(daf_obj, "cell", "metacell")
            result <- as.character(vec)
            if (is.null(names(result))) {
                names(result) <- dafr::axis_entries(daf_obj, "cell")
            }
            return(result)
        }
        return(NULL)
    }

    if (!dafr::has_vector(cells_daf, "cell", "metacell")) {
        return(NULL)
    }

    vec <- dafr::get_vector(cells_daf, "cell", "metacell")
    result <- as.character(vec)
    if (is.null(names(result))) {
        names(result) <- dafr::axis_entries(cells_daf, "cell")
    }
    result
}

# ==============================================================================
# Cell Type View Helper
# ==============================================================================

#' Create a DAF view filtered to specific cell types
#'
#' Uses DAF's viewer() with a mask on the cell axis to restrict to cells
#' whose metacell's type matches the requested cell types. The filter uses
#' a regex match on the chained lookup `metacell : type`.
#'
#' @param cells_daf Chained cells DAF (cells + metacells)
#' @param cell_types Character vector of cell types to include
#' @param metacell_types Optional metacell_types tibble (unused, for API compat)
#' @return A filtered DAF view
#' @noRd
make_cell_type_view <- function(cells_daf, cell_types, metacell_types = NULL) {
    # Escape regex special characters in type names and build OR pattern
    escaped_types <- gsub("([.+*?^${}()|\\[\\]\\\\])", "\\\\\\1", cell_types)
    types_pattern <- paste(escaped_types, collapse = "|")
    filter_query <- glue::glue("/ cell & metacell : type ~ ^({types_pattern})$")
    dafr::viewer(cells_daf, axes = list("cell" = filter_query))
}

# ==============================================================================
# Group Composition
# ==============================================================================

#' Get cell-count-based composition grouped by a categorical field
#'
#' Groups cells by the specified field and by their metacell-derived cell type,
#' then counts cells in each group x cell_type combination.
#'
#' @param dataset Dataset name
#' @param group_field Cell-level categorical field to group by
#'   (e.g., "batch_set_id", "embryo", "age_group")
#' @param metacell_types Current metacell_types tibble (for dynamic cell type
#'   assignment). If NULL, fetched from the dataset.
#' @return Tibble with columns: group_id, cell_type, n_cells, fraction
#' @export
get_group_cell_type_composition <- function(dataset, group_field, metacell_types = NULL) {
    # Get cell -> group mapping
    group_vec <- get_cell_field_map(dataset, group_field)
    if (is.null(group_vec)) {
        cli_abort("Field '{group_field}' not found on cell axis for dataset '{dataset}'")
    }

    # Get cell -> metacell mapping
    mc_vec <- get_cell_metacell_map(dataset)
    if (is.null(mc_vec)) {
        cli_abort("Cell-to-metacell mapping not found for dataset '{dataset}'")
    }

    # Get metacell -> cell_type mapping
    if (is.null(metacell_types)) {
        metacell_types <- get_metacell_types_data(dataset)
    }
    if (is.null(metacell_types)) {
        cli_abort("metacell_types not available for dataset '{dataset}'")
    }

    mc_to_type <- stats::setNames(
        as.character(metacell_types$cell_type),
        as.character(metacell_types$metacell)
    )

    # Use only cells that have both a group and a valid metacell
    cell_names <- names(group_vec)
    if (is.null(cell_names)) {
        cell_names <- names(mc_vec)
    }
    common_cells <- intersect(names(group_vec), names(mc_vec))
    if (length(common_cells) == 0) {
        # Try positional matching if names are not set
        n <- min(length(group_vec), length(mc_vec))
        group_sub <- group_vec[seq_len(n)]
        mc_sub <- mc_vec[seq_len(n)]
    } else {
        group_sub <- group_vec[common_cells]
        mc_sub <- mc_vec[common_cells]
    }

    # Map metacell to cell type
    cell_types <- mc_to_type[mc_sub]

    # Filter out cells with unknown metacells or types
    valid <- !is.na(cell_types) & !is.na(group_sub)
    group_sub <- group_sub[valid]
    cell_types <- cell_types[valid]

    # Count
    counts_df <- tibble::tibble(
        group_id = group_sub,
        cell_type = cell_types
    ) %>%
        dplyr::count(group_id, cell_type, name = "n_cells") %>%
        dplyr::group_by(group_id) %>%
        dplyr::mutate(
            group_total = sum(n_cells),
            fraction = n_cells / group_total
        ) %>%
        dplyr::ungroup()

    # Add Wilson score confidence intervals for each fraction
    ci <- calc_wilson_ci(counts_df$n_cells, counts_df$group_total)
    counts_df$fraction_lower <- ci$lower
    counts_df$fraction_upper <- ci$upper
    counts_df$group_total <- NULL

    return(counts_df)
}

# ==============================================================================
# Confidence Intervals
# ==============================================================================

#' Wilson score confidence interval for a proportion
#'
#' Computes the Wilson score interval, which is well-behaved even for
#' small sample sizes and proportions near 0 or 1.
#'
#' @param k Number of successes (e.g., cells of a given type)
#' @param n Total count (e.g., total cells in a group)
#' @param conf_level Confidence level (default: 0.95)
#' @return List with elements: lower, upper, estimate (the raw proportion)
#' @export
calc_wilson_ci <- function(k, n, conf_level = 0.95) {
    z <- qnorm(1 - (1 - conf_level) / 2)
    p_hat <- k / n
    denom <- 1 + z^2 / n
    center <- (p_hat + z^2 / (2 * n)) / denom
    margin <- z * sqrt((p_hat * (1 - p_hat) + z^2 / (4 * n)) / n) / denom
    list(lower = pmax(0, center - margin), upper = pmin(1, center + margin), estimate = p_hat)
}

# ==============================================================================
# Per-Gene Pseudobulk (On-Demand)
# ==============================================================================

#' Get per-group expression for a single gene
#'
#' Queries single-gene UMIs across all cells from the cells DAF, then
#' aggregates (sums) by the specified grouping field.
#'
#' @param dataset Dataset name
#' @param gene Gene name
#' @param group_field Cell-level field to group by
#' @param cell_types Optional cell types to filter by (via metacell -> type)
#' @param metacell_types Current metacell_types tibble (for cell type filtering)
#' @return Named numeric vector: group_value -> sum of UMIs for this gene
#' @export
get_group_gene_expression <- function(dataset, gene, group_field,
                                      cell_types = NULL, metacell_types = NULL) {
    cells_daf <- get_cells_daf(dataset)
    if (is.null(cells_daf)) {
        cli_abort("Cells DAF not available for dataset '{dataset}'")
    }

    # Use DAF GroupBy query, with optional cell type view for filtering
    query_daf <- cells_daf
    if (!is.null(cell_types)) {
        query_daf <- make_cell_type_view(cells_daf, cell_types, metacell_types)
    }

    gene_escaped <- escape_daf_value(gene)
    query_str <- glue::glue("/ cell / gene = {gene_escaped} : UMIs @ {group_field} %> Sum")
    result <- tryCatch(
        query_daf[query_str],
        error = function(e) {
            cli_abort("Failed to query gene '{gene}' UMIs grouped by '{group_field}': {conditionMessage(e)}")
        }
    )
    return(stats::setNames(as.numeric(result), names(result)))
}

# ==============================================================================
# Multi-Gene Pseudobulk Matrix
# ==============================================================================

#' Get pseudobulk matrix for multiple genes
#'
#' Queries per-gene UMIs across cells and aggregates by group, returning
#' a genes x groups matrix suitable for DE analysis.
#'
#' @param dataset Dataset name
#' @param genes Gene names. If NULL, uses all genes (can be expensive).
#' @param group_field Cell-level field to group by
#' @param cell_types Optional cell type filter
#' @param metacell_types Current metacell_types tibble
#' @return Matrix: genes (rows) x groups (columns)
#' @export
get_group_pseudobulk_mat <- function(dataset, genes, group_field,
                                     cell_types = NULL, metacell_types = NULL) {
    cells_daf <- get_cells_daf(dataset)
    if (is.null(cells_daf)) {
        cli_abort("Cells DAF not available for dataset '{dataset}'")
    }

    # If genes is NULL, get all genes from the gene axis
    if (is.null(genes)) {
        genes <- dafr::axis_entries(cells_daf, "gene")
        cli_alert_warning("Querying all {length(genes)} genes - this may be slow")
    }

    # Use DAF's native GroupBy query for the aggregation, with an optional
    # viewer to filter cells by type. All computation happens in Julia.
    query_daf <- cells_daf
    if (!is.null(cell_types)) {
        query_daf <- make_cell_type_view(cells_daf, cell_types, metacell_types)
    }

    query_str <- glue::glue("/ cell / gene : UMIs @ {group_field} %> Sum")
    daf_mat <- tryCatch(
        {
            # DAF returns group_values x gene; transpose to gene x group_values
            t(query_daf[query_str])
        },
        error = function(e) {
            cli_abort("Failed to compute pseudobulk grouped by '{group_field}': {conditionMessage(e)}")
        }
    )

    # Subset to requested genes if needed
    if (!is.null(genes) && length(genes) < nrow(daf_mat)) {
        gene_idx <- match(genes, rownames(daf_mat))
        gene_idx <- gene_idx[!is.na(gene_idx)]
        daf_mat <- daf_mat[gene_idx, , drop = FALSE]
    }

    return(daf_mat)
}

# ==============================================================================
# Group Differential Expression
# ==============================================================================

#' Differential expression between two groups defined by a categorical field
#'
#' Computes pseudobulk DE between two sets of values of a categorical field.
#' Output format matches calc_diff_expr() used elsewhere in MCView.
#'
#' @param dataset Dataset name
#' @param group_field Cell-level categorical field
#' @param group1_values Values of group_field for group 1
#' @param group2_values Values of group_field for group 2
#' @param cell_types Optional cell type filter
#' @param metacell_types Current metacell_types tibble
#' @param diff_thresh Log2 fold-change threshold for significance (default: 1.5)
#' @param pval_thresh p-value threshold for significance (default: 0.01)
#' @return Tibble matching calc_diff_expr output format:
#'   gene, group1 (egc), group2 (egc), diff, pval, col
#' @export
calc_group_diff_expr <- function(dataset, group_field,
                                 group1_values, group2_values,
                                 cell_types = NULL, metacell_types = NULL,
                                 diff_thresh = 1.5, pval_thresh = 0.01) {
    cells_daf <- get_cells_daf(dataset)
    if (is.null(cells_daf)) {
        cli_abort("Cells DAF not available for dataset '{dataset}'")
    }

    # Validate group field exists before running the expensive query
    if (!dafr::has_vector(cells_daf, "cell", group_field)) {
        cli_abort("Field '{group_field}' not found on cell axis")
    }

    # Use DAF's native GroupBy query with an optional cell type view.
    # Single query computes grouped sums for all genes at once in Julia.
    query_daf <- cells_daf
    if (!is.null(cell_types)) {
        query_daf <- make_cell_type_view(cells_daf, cell_types, metacell_types)
    }

    group1_name <- paste(group1_values, collapse = "+")
    group2_name <- paste(group2_values, collapse = "+")

    query_str <- glue::glue("/ cell / gene : UMIs @ {group_field} %> Sum")
    daf_result <- tryCatch(
        {
            # Returns group_values x gene; transpose to gene x group_values
            t(query_daf[query_str])
        },
        error = function(e) {
            cli_abort("Failed to compute DE pseudobulk grouped by '{group_field}': {conditionMessage(e)}")
        }
    )

    # Sum across group values within each group
    g1_cols <- intersect(group1_values, colnames(daf_result))
    g2_cols <- intersect(group2_values, colnames(daf_result))
    if (length(g1_cols) == 0 || length(g2_cols) == 0) {
        cli_abort("One or both groups have zero cells after filtering")
    }
    umis1 <- rowSums(daf_result[, g1_cols, drop = FALSE])
    umis2 <- rowSums(daf_result[, g2_cols, drop = FALSE])
    all_genes <- rownames(daf_result)

    # Compute enriched gene counts (EGC) -- fraction of total
    egc_epsilon <- mcv_get("egc_epsilon") %||% 1e-5
    total1 <- sum(umis1)
    total2 <- sum(umis2)

    egc1 <- umis1 / total1 + egc_epsilon
    egc2 <- umis2 / total2 + egc_epsilon

    # Build DE dataframe
    mat <- cbind(umis1, umis2)
    colnames(mat) <- c(group1_name, group2_name)
    rownames(mat) <- all_genes

    egc <- cbind(egc1, egc2)
    colnames(egc) <- c(group1_name, group2_name)
    rownames(egc) <- all_genes

    df <- calc_diff_expr(mat, egc, c(group1_name, group2_name), diff_thresh, pval_thresh)

    return(df)
}

# ==============================================================================
# QC Stats
# ==============================================================================

#' Per-group QC metrics
#'
#' Computes basic QC statistics per group: number of cells, total UMIs,
#' and median UMIs per cell.
#'
#' @param dataset Dataset name
#' @param group_field Cell-level field to group by
#' @return Tibble with columns: group_id, n_cells, total_umis, median_umis_per_cell
#' @export
get_group_qc_stats <- function(dataset, group_field) {
    cells_daf <- get_cells_daf(dataset)
    if (is.null(cells_daf)) {
        cli_abort("Cells DAF not available for dataset '{dataset}'")
    }

    # Use DAF's GroupBy queries for efficient per-group aggregation
    has_total_umis <- dafr::has_vector(cells_daf, "cell", "total_UMIs")

    # Count cells per group
    count_query <- glue::glue("/ cell : {group_field} @ {group_field} %> Count")
    counts <- tryCatch(cells_daf[count_query], error = function(e) NULL)
    if (is.null(counts)) {
        cli_abort("Failed to compute cell counts by '{group_field}'")
    }

    if (!has_total_umis) {
        return(tibble::tibble(
            group_id = names(counts),
            n_cells = as.integer(counts),
            total_umis = NA_real_,
            median_umis_per_cell = NA_real_
        ) %>% dplyr::arrange(group_id))
    }

    # Sum and median of total_UMIs per group
    sum_query <- glue::glue("/ cell : total_UMIs @ {group_field} %> Sum")
    med_query <- glue::glue("/ cell : total_UMIs @ {group_field} %> Median")
    total_sums <- cells_daf[sum_query]
    medians <- cells_daf[med_query]

    # Align by group name
    groups <- names(counts)
    tibble::tibble(
        group_id = groups,
        n_cells = as.integer(counts[groups]),
        total_umis = as.numeric(total_sums[groups]),
        median_umis_per_cell = as.numeric(medians[groups])
    ) %>% dplyr::arrange(group_id)
}

# ==============================================================================
# Backward Compatibility: Default Sample Field
# ==============================================================================

#' Get the default grouping field for "samples" behavior
#'
#' Checks for common sample-like fields in order of preference:
#' samp_id, batch_set_id, then falls back to the first available
#' categorical field.
#'
#' @param dataset Dataset name
#' @return Field name (character), or NULL if none available
#' @export
get_default_sample_field <- function(dataset) {
    # First check the metacells DAF for samp_id (original pattern)
    daf_obj <- get_dataset_daf(dataset)
    if (!is.null(daf_obj) && dafr::has_axis(daf_obj, "cell") &&
        dafr::has_vector(daf_obj, "cell", "samp_id")) {
        return("samp_id")
    }

    # Check the cells DAF for common sample-like fields
    cells_daf <- get_cells_daf(dataset)
    if (is.null(cells_daf)) {
        return(NULL)
    }

    preferred_fields <- c("samp_id", "batch_set_id", "batch", "sample")
    for (field in preferred_fields) {
        if (dafr::has_vector(cells_daf, "cell", field)) {
            return(field)
        }
    }

    # Fall back to first available categorical grouping field
    grouping_fields <- get_cell_grouping_fields(dataset)
    if (length(grouping_fields) > 0) {
        return(grouping_fields[1])
    }

    return(NULL)
}
