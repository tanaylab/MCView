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

load_shiny_data <- function(name, dataset, cache_dir) {
    cache_dir <- fs::path(cache_dir, dataset)

    flat_file <- fs::path(cache_dir, glue("{name}.tsv"))
    if (fs::file_exists(flat_file)) {
        object <- fread(flat_file) %>% as_tibble()
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

load_all_mc_data <- function(dataset, cache_dir) {
    files <- list.files(fs::path(cache_dir, dataset), pattern = "*\\.(qs|tsv|csv)")

    mc_data[[dataset]] <<- list()

    for (fn in files) {
        var_name <- basename(fn) %>%
            sub("\\.qs$", "", .) %>%
            sub("\\.tsv$", "", .)
        obj <- load_shiny_data(var_name, dataset, cache_dir)

        mc_data[[dataset]][[var_name]] <<- obj
    }    
}

verify_app_cache <- function(project, required_files = c("mc_mat.qs", "mc_sum.qs", "mc2d.qs", "gg_mc_top_cor.qs", "metacell_types.tsv", "cell_type_colors.tsv") ){
    cache_dir <- project_cache_dir(project)
    datasets <- dataset_ls(project)     
    
    for (dataset in datasets){
        dataset_dir <- fs::path(cache_dir, dataset)
        for (file in required_files){
            if (!fs::file_exists(fs::path(dataset_dir, file))){
                cli_abort("The file {.file {file}} is missing in {.file {dataset_dir}}. Did you forget to import?")
            }
        }
    }
}

load_all_data <- function(cache_dir) {
    datasets <- dataset_ls(project)

    mc_data <<- list()

    purrr::walk(datasets, ~ load_all_mc_data(dataset = .x, cache_dir = cache_dir))
}

get_cell_type_data <- function(dataset) {
    cell_type_colors <- mc_data[[dataset]][["cell_type_colors"]]
    cell_type_colors <- cell_type_colors %>%
        mutate(cell_type = as.character(cell_type)) %>%
        mutate(cell_type_id = as.character(1:n()))
    return(cell_type_colors)
}

get_metacell_types_data <- function(dataset) {
    metacell_types <- mc_data[[dataset]][["metacell_types"]]
    if (!is.factor(metacell_types$cell_type)) {
        metacell_types$cell_type <- factor(metacell_types$cell_type)
    }
    metacell_types <- metacell_types %>%
        mutate(cell_type = as.character(forcats::fct_explicit_na(cell_type))) %>%
        mutate(metacell = as.character(metacell))

    cell_type_colors <- get_cell_type_data(dataset)

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

get_mc_data <- function(dataset, var_name) {
    if (var_name == "metacell_types") {
        return(get_metacell_types_data(dataset))
    } else if (var_name == "cell_type_colors") {
        return(get_cell_type_data(dataset))
    }

    mc_data[[dataset]][[var_name]]
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
    !is.null(get_mc_config(dataset, "network"))
}

has_time <- function(dataset) {
    !is.null(get_mc_config(dataset, "time_bin_field")) && has_name(get_mc_config(dataset, "metacell_types"), "mc_age")
}
