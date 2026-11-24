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

# All the cell types that metacells can have: those with a color and "(Missing)", the type of the metacells
# without annotation. The default of the cell type filters, so that unannotated metacells are shown.
all_cell_types <- function(cell_type_colors) {
    unique(c(as.character(cell_type_colors$cell_type), "(Missing)"))
}

# Cell types as shown to the user: metacells without a cell type (NA, e.g. after their type was deleted) are "(Missing)"
cell_type_or_missing <- function(cell_type) {
    cell_type <- as.character(cell_type)
    cell_type[is.na(cell_type)] <- "(Missing)"
    cell_type
}
