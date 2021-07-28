#' Update metadata for a dataset
#'
#'
#' @param overwrite overwrite all existing metadata. If \code{FALSE} - would override only existing metadata fields.
#'
#'
#' @inheritParams import_dataset
#'
#' @export
update_metadata <- function(project,
                            dataset,
                            metadata = NULL,
                            metadata_fields = NULL,
                            anndata_file = NULL,
                            overwrite = FALSE) {
    cache_dir <- project_cache_dir(project)

    if (!is.null(metadata_fields)) {
        if (is.null(anndata_file)) {
            cli_abort("You have to provide the {.code anndata_file} parameter with {.code metadata_fields}")
        } else {
            if (!fs::file_exists(anndata_file)) {
                cli_abort("{anndata_file} doesn't exist. Maybe there is a typo?")
            }

            cli_alert_info("Reading {.file {anndata_file}}")
            library(anndata)
            adata <- anndata::read_h5ad(anndata_file)
        }
    }

    cli_alert_info("Processing metadata")
    metadata <- load_metadata(metadata, metadata_fields, adata)

    prev_metadata_file <- fs::path(cache_dir, dataset, "metadata.tsv")
    if (fs::file_exists(prev_metadata_file)) {
        cli_alert_info("Merging with previous metadata")
        prev_metadata <- fread(prev_metadata_file) %>% as_tibble()
        metadata <- metadata %>%
            select(-any_of(colnames(prev_metadata)[-1])) %>%
            left_join(prev_metadata, by = "metacell")
    }


    serialize_shiny_data(
        metadata %>% select(metacell, everything()),
        "metadata",
        dataset = dataset,
        cache_dir = cache_dir,
        flat = TRUE
    )

    cli_alert_success("Successfully updated metadata of dataset {.field {dataset}} at {.path {project}} project")
}

#' Update metadata for a dataset
#'
#'
#' @param overwrite overwrite all existing colors. If \code{FALSE} - would override only
#' the colors of existing metadata fields.
#'
#'
#' @inheritParams import_dataset
#'
#' @export
update_metadata_colors <- function(project,
                                   dataset,
                                   metadata_colors,
                                   overwrite = FALSE) {
    cache_dir <- project_cache_dir(project)

    cli_alert_info("Processing metadata colors")
    if (is.character(metadata_colors)) {
        metadata_colors <- yaml::read_yaml(metadata_colors) %>% as_tibble()
    }

    new_metadata_colors <- parse_metadata_colors(metadata_colors, metadata)

    prev_colors_file <- fs::path(cache_dir, dataset, "metadata_colors.qs")
    if (fs::file_exists(prev_metadata_file)) {
        metadata_colors <- qs::qread(prev_colors_file)
        for (f in names(metadata_colors)) {
            metadata_colors[[f]] <- NULL
        }
        metadata_colors <- c(metadata_colors, new_metadata_colors)
    } else {
        metadata_colors <- new_metadata_colors
    }

    serialize_shiny_data(metadata_colors, "metadata_colors", dataset = dataset, cache_dir = cache_dir)
}



#' Convert cell metadata to metacell metadata
#'
#' @export
cell_metadata_to_metacell <- function() {

}



load_metadata <- function(metadata, metadata_fields, metacells, adata) {
    metadata_df <- NULL
    if (!is.null(metadata)) {
        if (is.character(metadata)) {
            metadata <- tgutil::fread(metadata) %>% as_tibble()
        }
        metadata_df <- parse_metadata(metadata, metacells)
    }

    if (!is.null(metadata_fields)) {
        purrr::walk(metadata_fields, ~ {
            if (!(.x %in% colnames(adata$obs))) {
                cli_abort("{.field {.x}} is not present in the h5ad object.")
            }
        })

        metadata_df_obs <- adata$obs %>%
            select(one_of(metadata_fields)) %>%
            rownames_to_column("metacell")
        metadata_df_obs <- parse_metadata(metadata_df_obs, metacells)

        if (!is.null(metadata_df)) {
            duplicated_metadata <- intersect(colnames(metadata), colnames(metadata_df))
            duplicated_metadata <- duplicated_metadata[duplicated_metadata != "metacell"]
            if (length(duplicated_metadata) > 0) {
                cli_abort("The metadata fields {.field {fields}} appear both in {.code metadata} and {.code metadata_fields}", fields = paste(duplicated_metadata, collapse = ", "))
            }

            metadata <- metadata_df_obs %>% full_join(metadata_df, by = "metacell")
        } else {
            metadata <- metadata_df_obs
        }
    } else {
        metadata <- metadata_df
    }

    return(metadata)
}

parse_metadata <- function(metadata, metacells) {
    metadata <- metadata %>% as_tibble()

    if (!has_name(metadata_colors, "metacell")) {
        cli_abort("metadata doesn't have a field named {.field metacell}")
    }

    unknown_metacells <- metadata$metacell[!(metadata$metacell %in% metacells)]
    if (length(unknown_metacells) > 0) {
        cli_abort("Metadata contains metacells that are missing from the data: {.field {mcs}}", mcs = paste(unknown_metacells, collapse = ", "))
    }

    missing_metacells <- metacells[!(metacells %in% metadata$metacell)]
    if (length(missing_metacells) > 0) {
        cli_warning("Some metacells are missing from metadata: {.field {mcs}}", mcs = paste(missing_metacells, collapse = ", "))
    }

    fields <- metadata %>%
        select(-metacell) %>%
        colnames()

    metadata <- metadata %>%
        mutate_at(one_of(fields), as.numeric)

    purrr::walk(fields, ~ {
        if (all(is.na(metadata[[.x]]))) {
            cli_abort("The metadata variable {.field {.x}} is all NA. Is it numeric? Convert categorical variables to numeric using {.code cell_metadata_to_metacell} with {.code categorical=TRUE}")
        }
    })

    return(metadata)
}

parse_metadata_colors <- function(metadata_colors, metadata) {
    if (is.null(metadata)) {
        cli_abort("{.code metadata_colors} is set but {.code metadata} is {.code NULL}. Did you forget to set the {.code metadata} parameter?")
    }

    metadata_colors <- purrr::imap(metadata_colors, ~ {
        if (!has_name(metadata, .y)) {
            cli_abort("{.code metadata} doesen't have a field named {.field {.y}} but it exists in {.code metadata_colors}")
        }

        if (length(.x) == 0 || length(.x[[1]]) == 0) {
            cli_abort("{.field {.y}} in {.code metadata_colors} does not contain any colors.")
        }

        colors <- .x[[1]]

        if (length(.x) > 1) {
            breaks <- .x[[2]]
            if (length(breaks) != length(colors)) {
                cli_abort("In metadata colors field {.field {.y}}: length of {.code breaks} shuold be equal to {.code colors}")
            }
            return(list(
                colors = colors,
                breaks = breaks
            ))
        } else {
            return(list(
                colors = colors
            ))
        }
    })

    return(metadata_colors)
}
