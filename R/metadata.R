# ==============================================================================
# Metadata Conversion Utilities
# ==============================================================================

#' Convert cell metadata to metacell metadata
#'
#' @param cell_metadata data frame with a column named "cell_id" with
#' the cell id and other metadata columns, or a name of a delimited file which
#' contains such data frame.
#' @param cell_to_metacell data frame with a column named "cell_id" with cell id and
#' another column named "metacell" with the metacell the cell is part of, or a
#' name of a delimited file which contains such data frame.
#' @param func summary function for the cell metadata for non categorical metadata columns
#' (e.g. mean, median, sum)
#' @param categorical a vector with names of categorical variables. The returned data frame would have
#' a column for each category where the values are the fraction of cells with the
#' category in each metacell.
#' @param anndata_file path to \code{h5ad} file which contains the output of metacell2 pipeline (metacells python package).
#' @param metadata_fields names of fields in the anndata \code{object$obs} which contains metadata for each cell.
#' @param rm_outliers do not calculate statistics for cells that are marked as outliers (\code{outiler=TRUE} in \code{object$obs}) (only relevant when running \code{cell_metadata_to_metacell_from_h5ad})
#'
#'
#' @return A data frame with a column named "metacell" and
#' the metadata columns from \code{cell_metadata} summarized for each metacell using
#' \code{func} for non-categorical variables, and a column for each category of the categorical metadata variables
#' in\code{cell_metadata}, where the values are the fraction of cells with the category in each metacell.
#'
#' @description
#' Summarise cell metadata for each metacell. Metadata fields can be either numeric and then the summary function \code{func} is applied for the values of each field, or categorical metadata fields which are expanded to multiple
#' metadata columns with the fraction of cells (in each metacell) for every category. Variables that are either character,
#' factor or are explicitly set at \code{categorical} are treated as categorical.
#'
#'
#' \code{cell_metadata_to_metacell} converts cell metadata to metacell metadata from data frames. \cr
#' \code{cell_metadata_to_metacell_from_h5ad} extracts metadata fields and cell_to_metacell from cells h5ad file and
#' then runs \code{cell_metadata_to_metacell}.
#'
#' @examples
#' set.seed(60427)
#' n_cells <- 5e6
#' cell_metadata <- tibble::tibble(
#'     cell_id = 1:n_cells,
#'     md1 = sample(1:5, size = n_cells, replace = TRUE),
#'     md2 = rnorm(n = n_cells),
#'     md_categorical1 = sample(paste0("batch", 1:5), size = n_cells, replace = TRUE),
#'     md_categorical2 = sample(1:5, size = n_cells, replace = TRUE)
#' )
#'
#' cell_to_metacell <- tibble::tibble(
#'     cell_id = 1:n_cells,
#'     metacell = sample(0:1535, size = n_cells, replace = TRUE)
#' )
#' metadata <- cell_metadata_to_metacell(
#'     cell_metadata[, 1:3],
#'     cell_to_metacell
#' )
#' head(metadata)
#'
#' metadata1 <- cell_metadata_to_metacell(
#'     cell_metadata[, 1:3], cell_to_metacell,
#'     func = function(x) x * 2
#' )
#' head(metadata1)
#'
#'
#' metadata3 <- cell_metadata_to_metacell(
#'     cell_metadata,
#'     cell_to_metacell,
#'     categorical = c("md_categorical1", "md_categorical2")
#' )
#'
#' head(metadata3)
#' \dontrun{
#' cell_metadata_to_metacell_from_h5ad("cells.h5ad", c("pile", "age", "batch"), categorical = "batch")
#' }
#'
#' @export
cell_metadata_to_metacell <- function(cell_metadata, cell_to_metacell, func = mean, categorical = c()) {
    if (is.character(cell_metadata)) {
        cell_metadata <- tgutil::fread(cell_metadata) %>% as_tibble()
    }

    if (is.character(cell_to_metacell)) {
        cell_to_metacell <- tgutil::fread(cell_to_metacell) %>% as_tibble()
    }

    if (colnames(cell_metadata)[1] != "cell_id") {
        cli_abort("First column of {.code cell_metadata} is not named {.field cell_id} (it is named {.field {colnames(cell_metadata)[1]}})")
    }

    if (colnames(cell_to_metacell)[2] != "metacell") {
        cli_abort("Second column of {.code cell_to_metacell} is not named {.field metacell} (it is named {.field {colnames(cell_to_metacell)[2]}})")
    }

    if (length(categorical) > 0) {
        purrr::walk(categorical, ~ {
            if (!rlang::has_name(cell_metadata, .x)) {
                cli_abort("The field {.field {.x}} doesn't exist in {.field cell_metadata}")
            }
        })
    }

    md_vars <- colnames(cell_metadata)[-1]

    categorical_f <- purrr::map_lgl(cell_metadata[, -1], ~ is.factor(.x) || is.character(.x))
    categorical_vars <- sort(unique(c(categorical, md_vars[categorical_f])))
    numerical_vars <- sort(setdiff(md_vars, categorical_vars))

    cell_metadata <- cell_metadata %>%
        left_join(cell_to_metacell, by = "cell_id")

    metadata <- cell_to_metacell %>%
        distinct(metacell) %>%
        arrange(metacell)

    if (length(categorical_vars) > 0) {
        cli_alert_info(glue("Categorical variables: {vars}", vars = paste(categorical_vars, collapse = ", ")))
        for (v in categorical_vars) {
            metadata_var <- cell_metadata %>%
                select(metacell, !!sym(v)) %>%
                group_by(metacell, !!sym(v)) %>%
                summarise(cat_var_n = n(), .groups = "drop") %>%
                group_by(metacell) %>%
                mutate(cat_varp = cat_var_n / sum(cat_var_n)) %>%
                select(-cat_var_n) %>%
                spread(!!sym(v), cat_varp, fill = 0) %>%
                ungroup()
            colnames(metadata_var) <- c("metacell", paste0(v, ": ", colnames(metadata_var)[-1]))
            metadata <- metadata %>% left_join(metadata_var, by = "metacell")
        }
    }

    if (length(numerical_vars) > 0) {
        cli_alert_info(glue("Numerical variables: {vars}", vars = paste(numerical_vars, collapse = ", ")))

        cell_metadata <- cell_metadata %>%
            mutate(across(all_of(numerical_vars), as.numeric))
        purrr::walk(numerical_vars, ~ {
            if (all(is.na(cell_metadata[[.x]]))) {
                cli_abort("The cell_metadata variable {.field {.x}} is all NA. Is it numeric? Convert categorical variables to numeric by adding it to the {.code categorical} parameter")
            }
        })

        metadata_numeric <- cell_metadata %>%
            select(-cell_id) %>%
            group_by(metacell) %>%
            summarise(across(all_of(numerical_vars), func)) %>%
            ungroup()

        metadata <- metadata %>%
            left_join(metadata_numeric, by = "metacell")
    }

    metadata <- metadata %>%
        filter(metacell != "Outliers")

    return(metadata)
}

#'
#' @describeIn cell_metadata_to_metacell
#'
#' @export
cell_metadata_to_metacell_from_h5ad <- function(anndata_file, metadata_fields, func = mean, categorical = c(), rm_outliers = TRUE) {
    if (!fs::file_exists(anndata_file)) {
        cli_abort("{anndata_file} doesn't exist. Maybe there is a typo?")
    }

    cli_alert_info("Reading {.file {anndata_file}}")
    adata <- anndata::read_h5ad(anndata_file)

    if (!("metacell_name" %in% colnames(adata$obs))) {
        cli_abort("h5ad object doesn't have a 'metacell_name' field.")
    }

    if (length(metadata_fields) == 1 && metadata_fields == "all") {
        metadata_fields <- colnames(adata$obs)
        forbidden_fields <- c("metacell", "outlier", "metacell_name")
        metadata_fields <- metadata_fields[!(metadata_fields %in% forbidden_fields)]
    } else {
        purrr::walk(metadata_fields, ~ {
            if (!(.x %in% colnames(adata$obs))) {
                cli_abort("{.field {.x}} is not present in the h5ad object.")
            }
        })
    }


    df <- adata$obs %>%
        rownames_to_column("cell_id")

    if (rm_outliers && has_name(adata$obs, "outlier")) {
        df <- df %>% filter(!outlier)
    }

    cell_metadata <- df %>%
        select(cell_id, one_of(metadata_fields))

    cell_to_metacell <- df %>%
        select(cell_id, metacell = metacell_name)


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
        if (length(metadata_fields) == 1 && metadata_fields == "all") {
            metadata_fields <- colnames(adata$obs)
            forbidden_fields <- c("metacell")
            metadata_fields <- metadata_fields[!(metadata_fields %in% forbidden_fields)]
        } else {
            purrr::walk(metadata_fields, ~ {
                if (!(.x %in% colnames(adata$obs))) {
                    cli_abort("{.field {.x}} is not present in the h5ad object.")
                }
            })
        }

        metadata_df_obs <- adata$obs %>%
            select(one_of(metadata_fields)) %>%
            rownames_to_column("metacell")
        metadata_df_obs <- parse_metadata(metadata_df_obs, metacells)

        if (!is.null(metadata_df)) {
            duplicated_metadata <- intersect(colnames(metadata_df), colnames(metadata_df_obs))
            duplicated_metadata <- duplicated_metadata[duplicated_metadata != "metacell"]
            if (length(duplicated_metadata) > 0) {
                fields <- paste(duplicated_metadata, collapse = ", ")
                cli_abort("The metadata fields {.field {fields}} appear both in {.code metadata} (file or dataframe) and {.code metadata_fields} (h5ad file)")
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

    metadata <- metadata %>% filter(!is.na(metacell))

    metadata$metacell <- as.character(metadata$metacell)
    metacells <- as.character(metacells)

    unknown_metacells <- metadata$metacell[!(metadata$metacell %in% metacells)]
    unknown_metacells <- unknown_metacells[unknown_metacells != -1]
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

        if (length(.x) == 0) {
            cli_abort("{.field {.y}} in {.code metadata_colors} does not contain any colors.")
        }

        if (is_numeric_field(metadata, .y)) {
            if (length(.x[[1]]) == 0) {
                cli_abort("{.field {.y}} in {.code metadata_colors} does not contain any colors.")
            }

            colors <- .x[[1]]

            if (length(.x) > 1) {
                breaks <- .x[[2]]
                if (length(breaks) != length(colors)) {
                    cli_abort("In metadata colors field {.field {.y}}: length of {.code breaks} should be equal to {.code colors}")
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
        } else {
            if (is.list(.x) && length(names(.x)) == 2 && all(names(.x) == c("colors", "categories"))) {
                res <- .x$colors
                names(res) <- .x$categories
            } else {
                if (is.null(names(.x))) {
                    cli_abort("In metadata colors field {.field {.y}}, color vector doesn't have any names.")
                }

                res <- .x
            }
            return(res)
        }
    })

    return(metadata_colors)
}
