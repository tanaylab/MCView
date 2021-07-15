sanitize_metacell_types <- function(metacell_types, cell_type_colors, dataset) {
    # Look for cell_type_ids that are missing
    if (!all(metacell_types$cell_type_id %in% cell_type_colors$cell_type_id)) {
        bad_inds <- !is.na(metacell_types$cell_type_id) & !(metacell_types$cell_type_id %in% cell_type_colors$cell_type_id)
        if (length(bad_inds) > 0) {
            bad_cell_types <- unique(metacell_types$cell_type_id[bad_inds]) %>% paste(collapse = ", ")
            showNotification(glue("The following cell types ids are invalid: {bad_cell_types}"), type = "warning", duration = 10)
            metacell_types$cell_type_id[bad_inds] <- NA
            metacell_types$cell_type[bad_inds] <- NA
        }
    }

    # Look for cell_type_ids with the wrong name
    colliding_cell_types <- cell_type_colors %>%
        distinct(cell_type_id, cell_type) %>%
        left_join(metacell_types %>%
            select(cell_type_id, new_cell_type = cell_type) %>%
            distinct(cell_type_id, new_cell_type),
        by = "cell_type_id"
        ) %>%
        filter(cell_type != new_cell_type)

    if (nrow(colliding_cell_types) > 0) {
        cell_types_str <- purrr::map2_chr(colliding_cell_types$cell_type, colliding_cell_types$new_cell_type, ~ glue("{.x}={.y}")) %>% paste(collapse = ",")
        showNotification(glue("The names of the following cell type ids was changed to the one already loaded in cell type annotation: {cell_types_str}"), type = "warning", duration = 10)
    }

    metacell_types <- metacell_types %>%
        select(-cell_type) %>%
        left_join(cell_type_colors %>%
            distinct(cell_type_id, cell_type), by = "cell_type_id")

    return(metacell_types)
}

cell_type_to_cell_type_id <- function(cell_type, cell_type_colors) {
    cell_type_colors$cell_type_id[cell_type_colors$cell_type == cell_type]
}
