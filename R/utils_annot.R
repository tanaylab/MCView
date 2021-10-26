sanitize_metacell_types <- function(metacell_types, cell_type_colors, dataset) {
    # Look for cell types that are missing in colors
    if (!all(metacell_types$cell_type %in% cell_type_colors$cell_type)) {
        bad_inds <- !is.na(metacell_types$cell_type) & (metacell_types$cell_type != "(Missing)") & !(metacell_types$cell_type %in% cell_type_colors$cell_type)
        if (sum(bad_inds) > 0) {
            bad_cell_types <- unique(metacell_types$cell_type[bad_inds]) %>% paste(collapse = ", ")
            showNotification(glue("The following cell types are invalid: {bad_cell_types}"), type = "warning", duration = 10)
            metacell_types$cell_type[bad_inds] <- NA
            metacell_types$cell_type[bad_inds] <- NA
        }
    }

    return(metacell_types)
}
