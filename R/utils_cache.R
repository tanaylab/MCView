serialize_shiny_data <- function(object, name, dataset, cache_dir, df2mat = FALSE, preset = "fast", flat = FALSE, ...) {
    dataset_dir <- fs::path(cache_dir, dataset)

    if (!fs::dir_exists(dataset_dir)) {
        fs::dir_create(dataset_dir)
    }

    if (df2mat) {
        object <- as.data.frame(object) %>% rownames_to_column("__rowname__")
    }

    if (flat) {
        fwrite(object, fs::path(dataset_dir, glue("{name}.tsv")), sep = "\t")
    } else {
        qs::qsave(object, fs::path(dataset_dir, glue("{name}.qs")), preset = preset, ...)
    }

    cli_alert_success_verbose("saved {.field {name}}")
}

load_shiny_data <- function(name, dataset, cache_dir, atlas = FALSE) {
    if (atlas) {
        cache_dir <- fs::path(cache_dir, dataset, "atlas")
    } else {
        cache_dir <- fs::path(cache_dir, dataset)
    }

    flat_file <- fs::path(cache_dir, glue("{name}.tsv"))
    if (fs::file_exists(flat_file)) {
        object <- fread(flat_file) %>% as_tibble()
        if (has_name(object, "metacell")) {
            object$metacell <- as.character(object$metacell)
        }
    } else {
        object <- qs::qread(fs::path(cache_dir, glue("{name}.qs")))
    }

    if (is.data.frame(object) && rlang::has_name(object, "__rowname__")) {
        object <- object %>%
            remove_rownames() %>%
            column_to_rownames("__rowname__") %>%
            as.matrix()
    }
    return(object)
}

load_all_mc_data_atlas <- function(dataset, cache_dir) {
    atlas_dir <- fs::path(cache_dir, dataset, "atlas")
    if (fs::dir_exists(atlas_dir)) {
        files <- list.files(atlas_dir, pattern = "*\\.(qs|tsv|csv)")

        if (is.null(mc_data[[dataset]])) {
            mc_data[[dataset]] <<- list()
        }

        if (is.null(mc_data[[dataset]]$atlas)) {
            mc_data[[dataset]]$atlas <<- list()
        }

        for (fn in files) {
            var_name <- basename(fn) %>%
                sub("\\.qs$", "", .) %>%
                sub("\\.tsv$", "", .)
            obj <- load_shiny_data(var_name, dataset, cache_dir, atlas = TRUE)

            mc_data[[dataset]]$atlas[[var_name]] <<- obj
        }
    }
}

load_all_mc_data <- function(dataset, cache_dir) {
    atlas_dir <- fs::path(cache_dir, dataset, "atlas")
    if (fs::dir_exists(atlas_dir)) {
        load_all_mc_data_atlas(dataset, cache_dir)
    }

    files <- list.files(fs::path(cache_dir, dataset), pattern = "*\\.(qs|tsv|csv)")

    if (is.null(mc_data[[dataset]])) {
        mc_data[[dataset]] <<- list()
    }

    for (fn in files) {
        var_name <- basename(fn) %>%
            sub("\\.qs$", "", .) %>%
            sub("\\.tsv$", "", .)
        obj <- load_shiny_data(var_name, dataset, cache_dir)

        mc_data[[dataset]][[var_name]] <<- obj
    }
}

verify_app_cache <- function(project, required_files = c("mc_mat.qs", "mc_sum.qs", "mc2d.qs", "metacell_types.tsv", "cell_type_colors.tsv"), datasets = NULL) {
    cache_dir <- project_cache_dir(project)
    if (is.null(datasets)) {
        datasets <- dataset_ls(project)
    }

    for (dataset in datasets) {
        dataset_dir <- fs::path(cache_dir, dataset)
        for (file in required_files) {
            if (!fs::file_exists(fs::path(dataset_dir, file))) {
                cli_abort("The file {.file {file}} is missing in {.file {dataset_dir}}. Did you forget to import?")
            }
        }
    }
}

load_all_data <- function(cache_dir, datasets = NULL) {
    if (is.null(datasets)) {
        datasets <- dataset_ls(project)
    }

    mc_data <<- list()

    purrr::walk(datasets, ~ load_all_mc_data(dataset = .x, cache_dir = cache_dir))
}

get_cell_type_data <- function(dataset, atlas = FALSE) {
    if (atlas) {
        cell_type_colors <- mc_data[[dataset]]$atlas[["cell_type_colors"]]
    } else {
        cell_type_colors <- mc_data[[dataset]][["cell_type_colors"]]
    }

    cell_type_colors <- cell_type_colors %>%
        filter(!is.na(cell_type)) %>%
        mutate(cell_type = as.character(cell_type))
    return(cell_type_colors)
}

get_metacell_types_data <- function(dataset, atlas = FALSE) {
    if (atlas) {
        metacell_types <- mc_data[[dataset]]$atlas[["metacell_types"]]
    } else {
        metacell_types <- mc_data[[dataset]][["metacell_types"]]
    }

    if (!is.factor(metacell_types$cell_type)) {
        metacell_types$cell_type <- factor(metacell_types$cell_type)
    }
    metacell_types <- metacell_types %>%
        mutate(cell_type = as.character(forcats::fct_na_value_to_level(cell_type, "(Missing)"))) %>%
        mutate(metacell = as.character(metacell))

    cell_type_colors <- get_cell_type_data(dataset, atlas = atlas)

    metacell_types <- metacell_types %>%
        left_join(cell_type_colors %>% select(cell_type, mc_col = color), by = "cell_type")

    return(metacell_types)
}

get_mc_color_key <- function(dataset) {
    get_metacell_types_data(dataset) %>%
        distinct(cell_type, mc_col) %>%
        rename(group = cell_type, color = mc_col) %>%
        as.data.frame()
}

get_metadata <- function(dataset, atlas = FALSE) {
    if (atlas) {
        return(mc_data[[dataset]]$atlas[["metadata"]])
    } else {
        metadata <- mc_data[[dataset]][["metadata"]]
    }

    if (!is.null(metadata)) {
        mc_qc_metadata <- mc_data[[dataset]][["mc_qc_metadata"]]
        if (!is.null(mc_qc_metadata)) {
            same_colnames <- intersect(colnames(mc_qc_metadata), colnames(metadata))
            same_colnames <- same_colnames[same_colnames != "metacell"]

            mc_qc_metadata <- mc_qc_metadata %>%
                select(-any_of(same_colnames))

            if (!is.null(mc_qc_metadata)) {
                metadata <- metadata %>%
                    left_join(mc_qc_metadata, by = "metacell")
            }
        }
    }


    return(metadata)
}

get_mc_data <- function(dataset, var_name, atlas = FALSE) {
    if (var_name == "metacell_types") {
        return(get_metacell_types_data(dataset, atlas = atlas))
    } else if (var_name == "cell_type_colors") {
        return(get_cell_type_data(dataset, atlas = atlas))
    } else if (var_name == "metadata") {
        return(get_metadata(dataset, atlas = atlas))
    }

    if (atlas) {
        return(mc_data[[dataset]]$atlas[[var_name]])
    } else {
        return(mc_data[[dataset]][[var_name]])
    }
}

get_mc_config <- function(dataset, var_name) {
    if (is.null(config$datasets)) {
        return(NULL)
    }
    if (is.null(config$datasets[[dataset]])) {
        return(NULL)
    }
    config$datasets[[dataset]][[var_name]]
}

has_network <- function(dataset) {
    !is.null(get_mc_data(dataset, "mc_network"))
}

has_time <- function(dataset) {
    !is.null(get_mc_data(dataset, "time_annot"))
}

has_metadata <- function(dataset, atlas = FALSE) {
    !is.null(get_mc_data(dataset, "metadata", atlas = atlas))
}

has_cell_metadata <- function(dataset) {
    !is.null(get_mc_data(dataset, "cell_metadata"))
}

has_samples <- function(dataset) {
    if (!has_cell_metadata(dataset)) {
        return(FALSE)
    }
    cell_md <- get_mc_data(dataset, "cell_metadata")
    return(rlang::has_name(cell_md, "samp_id"))
}

has_projection <- function(dataset) {
    !is.null(get_mc_data(dataset, "query_md"))
}

has_corrected <- function(dataset) {
    !is.null(get_mc_data(dataset, "mc_mat_corrected"))
}

has_atlas <- function(dataset) {
    !is.null(mc_data[[dataset]]$atlas)
}

any_has_atlas <- function(project) {
    any(purrr::map_lgl(dataset_ls(project), has_atlas))
}

calc_samp_mc_count <- function(dataset) {
    metadata <- get_mc_data(dataset, "cell_metadata")
    samp_mc_count <- metadata %>%
        filter(metacell != -1) %>%
        filter(!is.na(samp_id)) %>%
        count(samp_id, metacell) %>%
        spread(metacell, n, fill = 0) %>%
        column_to_rownames("samp_id") %>%
        as.matrix()
    mc_data[[dataset]][["samp_mc_count"]] <<- samp_mc_count
    return(samp_mc_count)
}

get_samp_mc_count <- function(dataset) {
    return(get_mc_data(dataset, "samp_mc_count") %||% calc_samp_mc_count(dataset))
}

calc_samp_mc_frac <- function(dataset) {
    metadata <- get_mc_data(dataset, "cell_metadata")
    samp_mc_count <- get_samp_mc_count(dataset)
    samp_mc_frac <- samp_mc_count / rowSums(samp_mc_count)
    mc_data[[dataset]][["samp_mc_frac"]] <<- samp_mc_frac
    return(samp_mc_frac)
}

get_samp_mc_frac <- function(dataset) {
    return(get_mc_data(dataset, "samp_mc_frac") %||% calc_samp_mc_frac(dataset))
}

calc_samp_metadata <- function(dataset) {
    metadata <- get_mc_data(dataset, "cell_metadata") %>% select(-any_of(c("metacell", "cell_id", "outlier")))
    samp_columns <- metadata %>%
        group_by(samp_id) %>%
        summarise_all(n_distinct)
    samp_columns <- colnames(samp_columns)[purrr::map_lgl(colnames(samp_columns), ~ all(samp_columns[[.x]] == 1))]
    samp_md <- metadata %>%
        select(samp_id, samp_columns) %>%
        arrange(samp_id) %>%
        distinct(samp_id, .keep_all = TRUE)

    mc_data[[dataset]][["samp_metadata"]] <<- samp_md
    return(samp_md)
}

get_samp_metadata <- function(dataset) {
    if (!has_cell_metadata(dataset)) {
        return(NULL)
    }
    return(get_mc_data(dataset, "samp_metadata") %||% calc_samp_metadata(dataset))
}

calc_samples_list <- function(dataset) {
    if (!has_cell_metadata(dataset)) {
        return(NULL)
    }
    samp_md <- get_samp_metadata(dataset)
    samp_list <- sort(unique(samp_md$samp_id))
    mc_data[[dataset]][["samp_list"]] <<- samp_list
    return(samp_list)
}

get_samples_list <- function(dataset) {
    get_mc_data(dataset, "samp_list") %||% calc_samples_list(dataset)
}

dataset_metadata_fields <- function(dataset, atlas = FALSE) {
    metadata <- get_mc_data(dataset, "metadata", atlas = atlas)
    if (is.null(metadata)) {
        return(c())
    }
    fields <- colnames(metadata)
    fields <- fields[fields != "metacell"]
    fields <- fields[!(fields %in% c("samp_id", "cell_id"))]
    fields <- fields[!grepl("samp_id: ", fields)]
    return(fields)
}

dataset_metadata_fields_numeric <- function(dataset, atlas = FALSE) {
    fields <- dataset_metadata_fields(dataset, atlas = atlas)
    df <- get_mc_data(dataset, "metadata", atlas = atlas)
    numeric_f <- purrr::map_lgl(fields, ~ is_numeric_field(df, .x))
    return(fields[numeric_f])
}

dataset_cell_metadata_fields <- function(dataset, atlas = FALSE) {
    metadata <- get_mc_data(dataset, "cell_metadata", atlas = atlas)
    if (is.null(metadata)) {
        return(c())
    }
    fields <- colnames(metadata)
    fields <- fields[!(fields %in% c("samp_id", "cell_id"))]
    fields <- fields[!grepl("samp_id: ", fields)]

    return(fields)
}

dataset_cell_metadata_fields_numeric <- function(dataset, atlas = FALSE) {
    fields <- dataset_cell_metadata_fields(dataset, atlas = atlas)
    df <- get_mc_data(dataset, "cell_metadata")
    numeric_f <- purrr::map_lgl(fields, ~ is_numeric_field(df, .x))
    return(fields[numeric_f])
}

dataset_cell_metadata_fields_categorical <- function(dataset, atlas = FALSE) {
    fields <- dataset_cell_metadata_fields(dataset, atlas = atlas)
    df <- get_mc_data(dataset, "cell_metadata")
    numeric_f <- purrr::map_lgl(fields, ~ is_numeric_field(df, .x))
    return(fields[!numeric_f])
}

is_numeric_field <- function(df, field) {
    if (is.null(df)) {
        return(FALSE)
    }

    if (is.null(field)) {
        return(FALSE)
    }

    if (is.character(df[[field]])) {
        return(FALSE)
    }

    if (is.logical(df[[field]])) {
        return(FALSE)
    }

    if (is.numeric(df[[field]])) {
        return(TRUE)
    }

    return(FALSE)
}

get_metacell_ids <- function(project, dataset) {
    colnames(load_shiny_data("mc_mat", dataset, project_cache_dir(project)))
}

has_gg_mc_top_cor <- function(project, dataset) {
    !is.null(get_mc_data(dataset, "gg_mc_top_cor"))
}

get_mc_sum <- function(dataset, atlas = FALSE) {
    mc_sum <- get_mc_data(dataset, "mc_sum", atlas = atlas)
    mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
    names(mc_sum) <- colnames(mc_mat)

    return(mc_sum)
}

get_gene_qc <- function(dataset) {
    gene_qc <- get_mc_data(dataset, "gene_qc")
    if (is.null(gene_qc)) {
        gene_qc <- get_mc_data(dataset, "gene_inner_fold")
    }

    return(gene_qc)
}
