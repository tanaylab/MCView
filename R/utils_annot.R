sanitize_metacell_type<- function(mc_annot, cell_type_annot, dataset) {
    # Look for cell_type_ids that are missing
    if (!all(mc_annot$cell_type_id %in% cell_type_annot$cell_type_id)) {
        bad_inds <- !is.na(mc_annot$cell_type_id) & !(mc_annot$cell_type_id %in% cell_type_annot$cell_type_id)
        if (length(bad_inds) > 0) {
            bad_cell_types <- unique(mc_annot$cell_type_id[bad_inds]) %>% paste(collapse = ", ")
            showNotification(glue("The following cell types ids are invalid: {bad_cell_types}"), type = "warning", duration = 10)
            mc_annot$cell_type_id[bad_inds] <- NA
            mc_annot$cell_type[bad_inds] <- NA
        }
    }

    # Look for cell_type_ids with the wrong name
    colliding_cell_types <- cell_type_color%>%
        distinct(cell_type_id, cell_type) %>%
        left_join(metacell_type%>%
            select(cell_type_id, new_cell_type = cell_type) %>%
            distinct(cell_type_id, new_cell_type),
        by = "cell_type_id"
        ) %>%
        filter(cell_type != new_cell_type)

    if (nrow(colliding_cell_types) > 0) {
        cell_types_str <- purrr::map2_chr(colliding_cell_types$cell_type, colliding_cell_types$new_cell_type, ~ glue("{.x}={.y}")) %>% paste(collapse = ",")
        showNotification(glue("The names of the following cell type ids was changed to the one already loaded in cell type annotation: {cell_types_str}"), type = "warning", duration = 10)
    }

    metacell_type<- metacell_type%>%
        select(-cell_type) %>%
        left_join(cell_type_color%>%
            distinct(cell_type_id, cell_type), by = "cell_type_id")

    return(mc_annot)
}

cell_type_to_cell_type_id <- function(cell_type, cell_type_annot) {
    cell_type_annot$cell_type_id[cell_type_annot$cell_type == cell_type]
}
