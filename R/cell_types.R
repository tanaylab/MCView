# ==============================================================================
# Cell Type Parsing Utilities
# ==============================================================================

#' Parse cell type colors from file or data frame
#'
#' @param cell_type_colors A file path or data frame with cell type colors
#' @return A tibble with cell_type, color, and order columns
#' @noRd
parse_cell_type_colors <- function(cell_type_colors) {
    if (is.character(cell_type_colors)) {
        file <- cell_type_colors
        cell_type_colors <- fread(cell_type_colors) %>% as_tibble()
    } else {
        file <- "cell_type_colors"
    }


    if (!has_name(cell_type_colors, "cell_type") && !has_name(cell_type_colors, "cluster")) {
        cli_abort("{.field {file}} should have a column named {.field cell_type} or {.field cluster}")
    }

    if (!has_name(cell_type_colors, "color")) {
        cli_abort("{.field {file}} should have a column named {.field color}")
    }

    if (rlang::has_name(cell_type_colors, "cluster")) {
        cell_type_colors <- cell_type_colors %>% rename(cell_type = cluster)
    }

    cell_type_colors <- cell_type_colors %>%
        filter(!is.na(cell_type), !is.na(color)) %>%
        filter(cell_type != "(Missing)") %>%
        distinct(cell_type, .keep_all = TRUE)

    if (!has_name(cell_type_colors, "order")) {
        cell_type_colors <- cell_type_colors %>% mutate(order = 1:n())
    }

    cell_type_colors <- cell_type_colors %>%
        distinct(cell_type, color, order) %>%
        select(cell_type, color, order)

    n_colors <- cell_type_colors %>%
        count(cell_type) %>%
        pull(n)

    if (any(n_colors > 1)) {
        cli_abort("Some cell types appear more than once with different colors.")
    }

    return(cell_type_colors)
}


parse_metacell_types <- function(metacell_types, metacells = NULL) {
    if (is.character(metacell_types)) {
        metacell_types <- fread(metacell_types) %>% as_tibble()
    }

    if (!has_name(metacell_types, "metacell")) {
        cli_abort("{.field {file}} should have a column named {.field metacell}")
    }

    if (!has_name(metacell_types, "cell_type") && !has_name(metacell_types, "cluster")) {
        cli_abort("{.field {file}} should have a column named {.field cell_type} or {.field cluster}")
    }

    if (rlang::has_name(metacell_types, "cluster")) {
        metacell_types <- metacell_types %>% rename(cell_type = cluster)
    }

    metacell_types <- metacell_types %>%
        select(any_of(c("metacell", "cell_type", "age", "mc_age")))


    if ("age" %in% colnames(metacell_types)) {
        metacell_types <- metacell_types %>%
            rename(mc_age = age)
    }

    metacell_types <- metacell_types %>%
        mutate(metacell = as.character(metacell))

    if (!is.null(metacells)) {
        unknown_metacells <- metacell_types$metacell[!(metacell_types$metacell %in% metacells)]
        if (length(unknown_metacells) > 0) {
            mcs <- paste(unknown_metacells, collapse = ", ")
            cli_abort("Metacell types contains metacells that are missing from the data: {.field {mcs}}")
        }

        missing_metacells <- metacells[!(metacells %in% metacell_types$metacell)]
        if (length(missing_metacells) > 0) {
            mcs <- paste(missing_metacells, collapse = ", ")
            cli_warn("Some metacells are missing from metacell types: {.field {mcs}}")
        }
    }

    return(metacell_types)
}
