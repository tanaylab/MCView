#' Update cell type assignment for each metacell
#'
#' Change the cell type assignments for each metacell to the ones listed at \code{metacell_types_file}.
#'
#' This is usually done after a first iteration of annotation using the "Annotate" tab in the MCView annotation, which can export a valid \code{metacell_types_file}.
#' The file should have a column named "metacell" with the metacell ids and another
#' column named "cell_type" or "cluster" with the cell type assignment.
#'
#' Note that the exported file from the __MCView__ app contains additional fields
#' which will be ignored in this function.
#'
#' Under the hood - MCView updates a file named "metacell_types.tsv" under \code{project/cache/dataset}, which can also be edited manually.
#'
#' If the file contains an additional 'color' field, the cell type colors would be updated as well.
#'
#' @param project path to the project directory
#' @param dataset name for the dataset, e.g. "PBMC"
#' @param metacell_types_file path to a tabular file (csv,tsv) with cell type assignement for
#' each metacell. The file should have a column named "metacell" with the metacell ids and another
#' column named "cell_type" or "cluster" with the cell type assignment. Metacell ids that do
#' not exists in the data would be ignored.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' update_metacell_types("PBMC", "PBMC163k", "raw/metacell-clusters.csv")
#' }
#'
#' @export
update_metacell_types <- function(project, dataset, metacell_types_file) {
    # an old bundle has the project in a "project" subdirectory
    project <- init_project_dir(project)
    verify_app_cache(project)

    prev_metacell_types <- load_shiny_data("metacell_types", dataset, project_cache_dir(project))
    prev_metacell_types <- prev_metacell_types %>%
        mutate(metacell = as.character(metacell))

    raw_df <- fread(metacell_types_file)
    metacell_types <- parse_metacell_types(raw_df)

    metacell_types <- prev_metacell_types %>%
        select(-cell_type) %>%
        left_join(metacell_types %>% select(metacell, cell_type), by = "metacell")

    serialize_shiny_data(metacell_types, "metacell_types", dataset = dataset, cache_dir = project_cache_dir(project), flat = TRUE)

    cli_alert_success("Succesfully changed metacell cell type assignments")

    if (has_name(raw_df, "color")) {
        cli_alert_info("File has a field named 'color', updating also cell type colors.")
        update_cell_type_colors(project, dataset, raw_df)
    }
}


#' Update color assignment for each cell type
#'
#' Change the color assignments for each cell type to the ones listed at \code{cell_type_colors_file}.
#'
#' This is usually done after a first iteration of annotation using the "Annotate" tab in the MCView annotation, which can
#' export a valid \code{cell_type_colors_file}.
#'
#' The file should have a column named "cell_type" or "cluster" with the cell types and another column named "color" with the color assignment.
#' Note that the exported file from the __MCView__ app contains additional fields which will be
#' ignored in this function.
#'
#' Under the hood - MCView updates a file named "cell_type_colors.tsv" under \code{project/cache/dataset}, which can also be edited manually.
#'
#' @param project path to the project directory
#' @param dataset name for the dataset, e.g. "PBMC"
#' @param cell_type_colors_file path to a tabular file (csv,tsv) with color assignement for
#' each cell type. The file should have a column named "cell_type" or "cluster" with the
#' cell types and another column named "color" with the color assignment. Cell types that do not
#' exist in the metacell types would be ignored, so if you changed the names of cell types you would have to also
#' update the metacell types (using \code{update_metacell_types}).
#' The function also accepts output of the 'export' button from the application annotation page.
#' If this parameter is missing, MCView would use the \code{chameleon} package to assign a color for each cell type.
#'
#'
#' @examples
#' \dontrun{
#' update_metacell_types("PBMC", "PBMC163k", "raw/cluster-colors.csv")
#' }
#'
#' @export
update_cell_type_colors <- function(project, dataset, cell_type_colors_file) {
    # an old bundle has the project in a "project" subdirectory
    project <- init_project_dir(project)
    verify_app_cache(project)

    cell_type_colors <- parse_cell_type_colors(cell_type_colors_file)

    serialize_shiny_data(cell_type_colors, "cell_type_colors", dataset = dataset, cache_dir = project_cache_dir(project), flat = TRUE)

    cli_alert_success("Succesfully changed cell type color assignments")
}

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

# The natural order of cell type names: the numbers in the names are compared as numbers, so that cell types named
# 1, 2, 10 (or C1, C2, C10) are ordered so and not 1, 10, 2 as sort() orders them. Names without numbers are ordered
# as sort() orders them, and missing names last.
natural_order_cell_types <- function(cell_types) {
    cell_types <- as.character(cell_types)
    # pad every number with zeros to the same width, so that comparing the names as text compares the numbers
    padded <- cell_types
    named <- !is.na(cell_types)
    padded_named <- cell_types[named]
    digits <- gregexpr("[0-9]+", padded_named)
    regmatches(padded_named, digits) <- lapply(regmatches(padded_named, digits), function(d) {
        d <- sub("^0+(?=[0-9])", "", d, perl = TRUE)
        paste0(strrep("0", pmax(0, 30 - nchar(d))), d)
    })
    padded[named] <- padded_named
    return(order(padded, cell_types))
}

natural_sort_cell_types <- function(cell_types) {
    return(cell_types[natural_order_cell_types(cell_types)])
}

# The order in which to show cell types (e.g. the rows of a cell type colors table): their own, unless it is the
# order of sort(), which natural_sort_cell_types() replaces (cell types named by numbers were ordered 1, 10, 2 by
# default in projects imported before the default order was natural). "(Missing)" keeps its place.
cell_types_display_order <- function(cell_types) {
    cell_types <- as.character(cell_types)
    named <- which(is.na(cell_types) | cell_types != "(Missing)")
    if (length(named) < 2 || is.unsorted(order(cell_types[named]))) {
        return(seq_along(cell_types))
    }
    ord <- seq_along(cell_types)
    ord[named] <- named[natural_order_cell_types(cell_types[named])]
    return(ord)
}
