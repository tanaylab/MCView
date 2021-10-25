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
    } else {
        adata <- NULL
    }

    cli_alert_info("Processing metadata")
    metacells <- get_metacell_ids(project, dataset)
    metadata <- load_metadata(metadata, metadata_fields, metacells, adata)

    prev_metadata_file <- fs::path(cache_dir, dataset, "metadata.tsv")
    if (fs::file_exists(prev_metadata_file) && !overwrite) {
        cli_alert_info("Merging with previous metadata")
        prev_metadata <- fread(prev_metadata_file) %>% as_tibble()
        metadata <- prev_metadata %>%
            select(-any_of(colnames(metadata)[-1])) %>%
            left_join(metadata, by = "metacell")
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

#' Update metadata colors for a dataset
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
    prev_metadata_file <- fs::path(cache_dir, dataset, "metadata.tsv")
    if (!fs::file_exists(prev_metadata_file)) {
        cli_abort("No metadata found for project {.field {project}}. Please call {.code update_metadata} to add it and then run {.code update_metadata_colors} again.")
    }
    metadata <- tgutil::fread(prev_metadata_file) %>% as_tibble()

    cli_alert_info("Processing metadata colors")
    if (is.character(metadata_colors)) {
        metadata_colors <- yaml::read_yaml(metadata_colors) %>% as_tibble()
    }

    new_metadata_colors <- parse_metadata_colors(metadata_colors, metadata)
    
    prev_colors_file <- fs::path(cache_dir, dataset, "metadata_colors.qs")
    if (fs::file_exists(prev_colors_file) && !overwrite) {
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
#' @param cell_metadata data frame with a column named "cell_id" with
#' the cell id and other metadata columns, or a name of a delimited file which
#' contains such data frame. if \code{categorical = TRUE} the data frame can have only
#' a single metadata categorical column.
#' @param cell_to_metacell data frame with a column named "cell_id" with cell id and
#' another column named "metacell" with the metacell the cell is part of, or a
#' name of a delimited file which contains such data frame.
#' @param func summary function for the cell metadata for non categorical metadata columns
#' (e.g. mean, median, sum)
#' @param categorical is the metadata categorical. if \code{TRUE} - \code{cell_metadata}
#' can have only a single metadata column, and the returned data frame would have
#' a column for each category where the values are the fraction of cells with the
#' category in each metacell.
#' @param anndata_file path to \code{h5ad} file which contains the output of metacell2 pipeline (metacells python package).
#' @param metadata_fields names of fields in the anndata \code{object$obs} which contains metadata for each cell.
#' @param rm_outliers do not calculate statistics for cells that are marked as outliers (\code{outiler=TRUE} in \code{object$obs})
#'
#'
#' @return if \code{categorical=FALSE} - a data frame with a column named "metacell" and
#' the metadata columns from \code{cell_metadata} summarized for each metacell using
#' \code{func}. If \code{categorical=TRUE} - a data frame with a column named "metacell"
#' and a column for each category of the *single* categorical metadata variable in
#' \code{cell_metadata}, where the values are the fraction of cells with the category
#' in each metacell.
#'
#' @description
#' Summarise cell metadata for each metacell. Metadata fields can be either numeric and then the summary function \code{func} is
#' applied for the values of each field, or a single categorical metadata field which is expanded to multiple metadata columns with
#' the fraction of cells (in each metacell) for every category.
#' \code{cell_metadata_to_metacell} converts cell metadata to metacell metadata from data frames.
#' \code{cell_metadata_to_metacell_from_h5ad} extracts metadata fields and cell_to_metacell from cells h5ad file and
#' then runs \code{cell_metadata_to_metacell}.
#'
#' @examples
#' set.seed(60427)
#' n_cells <- 5e6
#' cell_metadata <- tibble(
#'     cell_id = 1:n_cells,
#'     md1 = sample(1:5, size = n_cells, replace = TRUE),
#'     md2 = rnorm(n = n_cells)
#' )
#' cell_To_metacell <- tibble(
#'     cell_id = 1:n_cells,
#'     metacell = sample(0:1535, size = n_cells, replace = TRUE)
#' )
#' metadata <- cell_metadata_to_metacell(cell_metadata, cell_To_metacell)
#' head(metadata)
#'
#' metadata1 <- cell_metadata_to_metacell(cell_metadata, cell_To_metacell, func = function(x) x * 2)
#' head(metadata1)
#'
#' cell_metadata_categorical <- tibble(
#'     cell_id = 1:n_cells,
#'     md1 = sample(paste0("batch", 1:5), size = n_cells, replace = TRUE)
#' )
#'
#' metadata3 <- cell_metadata_to_metacell(cell_metadata_categorical, cell_To_metacell, categorical = TRUE)
#' head(metadata3)
#' \dontrun{
#' cell_metadata_to_metacell_from_h5ad("cells.h5ad", c("pile", "age"))
#' cell_metadata_to_metacell_from_h5ad("cells.h5ad", "batch", categorical = TRUE)
#' }
#'
#' @export
cell_metadata_to_metacell <- function(cell_metadata, cell_to_metacell, func = mean, categorical = FALSE) {
    if (is.character(cell_metadata)) {
        cell_metadata <- tgutil::fread(cell_metadata) %>% as_tibble()
    }

    if (colnames(cell_metadata)[1] != "cell_id") {
        cli_abort("First column of {.code cell_metadata} is not named {.field cell_id} (it is named {.field {colnames(cell_metadata)[1]}})")
    }

    if (colnames(cell_to_metacell)[1] != "cell_id") {
        cli_abort("First column of {.code cell_to_metacell} is not named {.field cell_id} (it is named {.field {colnames(cell_to_metacell)[1]}})")
    }

    if (colnames(cell_to_metacell)[2] != "metacell") {
        cli_abort("Second column of {.code cell_to_metacell} is not named {.field metacell} (it is named {.field {colnames(cell_to_metacell)[2]}})")
    }

    if (categorical) {
        if (ncol(cell_metadata) != 2) {
            cli_abort("When {.code categorical=TRUE}, {.code cell_metadata} should have only two columns, {.field cell_id} and the categorical metadata.")
        }

        cell_metadata <- cell_metadata %>%
            left_join(cell_to_metacell, by = "cell_id")

        colnames(cell_metadata)[2] <- "cat_var"

        metadata <- cell_metadata %>%
            select(metacell, cat_var) %>%
            group_by(metacell, cat_var) %>%
            summarise(cat_var_n = n(), .groups = "drop") %>%
            group_by(metacell) %>%
            mutate(cat_varp = cat_var_n / sum(cat_var_n)) %>%
            select(-cat_var_n) %>%
            spread(cat_var, cat_varp, fill = 0) %>%
            ungroup()
    } else {
        fields <- cell_metadata %>%
            select(-cell_id) %>%
            colnames()

        cell_metadata <- cell_metadata %>%
            mutate_at(fields, as.numeric)

        purrr::walk(fields, ~ {
            if (all(is.na(cell_metadata[[.x]]))) {
                cli_abort("The cell_metadata variable {.field {.x}} is all NA. Is it numeric? Convert categorical variables to numeric by setting {.code categorical=TRUE}")
            }
        })

        cell_metadata <- cell_metadata %>%
            left_join(cell_to_metacell, by = "cell_id")

        metadata <- cell_metadata %>%
            select(-cell_id) %>%
            group_by(metacell) %>%
            summarise_at(vars(fields), func) %>%
            ungroup()
    }

    return(metadata)
}

#'
#' @describeIn cell_metadata_to_metacell
#'
#' @export
cell_metadata_to_metacell_from_h5ad <- function(anndata_file, metadata_fields, func = mean, categorical = FALSE, rm_outliers = TRUE) {
    if (!fs::file_exists(anndata_file)) {
        cli_abort("{anndata_file} doesn't exist. Maybe there is a typo?")
    }

    cli_alert_info("Reading {.file {anndata_file}}")
    library(anndata)
    adata <- anndata::read_h5ad(anndata_file)

    if (!("metacell" %in% colnames(adata$obs))) {
        cli_abort("h5ad object doesn't have a 'metacell' field.")
    }

    purrr::walk(metadata_fields, ~ {
        if (!(.x %in% colnames(adata$obs))) {
            cli_abort("{.field {.x}} is not present in the h5ad object.")
        }
    })

    if (categorical && length(metadata_fields) > 1) {
        cli_abort("{.code metadata_fields} should have a single field when {.code categorical=TRUE}")
    }

    df <- adata$obs %>%
        rownames_to_column("cell_id")

    if (rm_outliers && has_name(adata$obs, "outlier")) {
        df <- df %>% filter(!outlier)
    }

    cell_metadata <- df %>%
        select(cell_id, one_of(metadata_fields))

    cell_to_metacell <- df %>%
        select(cell_id, metacell)


    cell_metadata_to_metacell(cell_metadata = cell_metadata, cell_to_metacell = cell_to_metacell, func = func, categorical = categorical)
}


load_metadata <- function(metadata, metadata_fields, metacells, adata = NULL) {
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
                fields <- paste(duplicated_metadata, collapse = ", ")
                cli_abort("The metadata fields {.field {fields}} appear both in {.code metadata} and {.code metadata_fields}")
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

    if (!has_name(metadata, "metacell")) {
        cli_abort("metadata doesn't have a field named {.field metacell}")
    }

    unknown_metacells <- metadata$metacell[!(metadata$metacell %in% metacells)]
    if (length(unknown_metacells) > 0) {
        mcs <- paste(unknown_metacells, collapse = ", ")
        cli_abort("Metadata contains metacells that are missing from the data: {.field {mcs}}")
    }

    missing_metacells <- metacells[!(metacells %in% metadata$metacell)]
    if (length(missing_metacells) > 0) {
        mcs <- paste(missing_metacells, collapse = ", ")
        cli_warn("Some metacells are missing from metadata: {.field {mcs}}")
    }

    fields <- metadata %>%
        select(-metacell) %>%
        colnames()

    metadata <- metadata %>%
        mutate_at(fields, as.numeric)

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
