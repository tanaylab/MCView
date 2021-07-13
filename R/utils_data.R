serialize_shiny_data <- function(object, name, dataset, data_dir, df2mat = FALSE, preset = "fast", verbose = getOption("MCView.verbose"), ...) {
    if (df2mat) {
        object <- as.data.frame(object) %>% rownames_to_column("__rowname__")
    }

    qs::qsave(object, fs::path(data_dir, dataset, glue("{name}.qs")), preset = preset, ...)

    if (!is.null(verbose) && verbose) {
        cli_alert_info(name)
    }
}

load_shiny_data <- function(name, dataset, data_dir) {
    data_dir <- fs::path(data_dir, dataset)

    object <- qs::qread(fs::path(data_dir, glue("{name}.qs")))
    if (is.data.frame(object) && rlang::has_name(object, "__rowname__")) {
        object <- object %>%
            remove_rownames() %>%
            column_to_rownames("__rowname__") %>%
            as.matrix()
    }
    return(object)
}

load_all_mc_data <- function(dataset, data_dir) {
    files <- list.files(fs::path(data_dir, dataset), pattern = "*.qs")

    mc_data[[dataset]] <<- list()

    for (fn in files) {
        var_name <- basename(fn) %>% sub("\\.qs$", "", .)
        obj <- load_shiny_data(var_name, dataset, data_dir)

        mc_data[[dataset]][[var_name]] <<- obj
    }

    # initialize top correlated genes for each dataset
    mc_data[[dataset]]$top_cor_genes <- list()
}

load_all_data <- function(data_dir) {
    metacells_data_dirs <- list.files(data_dir, full.names = FALSE)
    metacells <- names(config$metacells)

    for (mc in metacells) {
        if (!(mc %in% metacells_data_dirs)) {
            cli_abort("{mc} dataset doesn't have any data. Did you forget to import it?")
        }
    }

    mc_data <<- list()

    purrr::walk(metacells, ~ load_all_mc_data(dataset = .x, data_dir = data_dir))
}

get_cell_type_data <- function(dataset) {
    cell_type_annot <- mc_data[[dataset]][["cell_type_annot"]]
    cell_type_annot <- cell_type_annot %>% mutate(cell_type_id = as.character(1:n()))
    return(cell_type_annot)
}

get_mc_annot_data <- function(dataset) {
    mc_annot <- mc_data[[dataset]][["mc_annot"]]
    if (!is.factor(mc_annot$cell_type)) {
        mc_annot$cell_type <- factor(mc_annot$cell_type)
    }
    mc_annot <- mc_annot %>% mutate(cell_type = as.character(forcats::fct_explicit_na(cell_type)))
    return(mc_annot)
}

get_mc_color_key <- function(dataset) {
    get_mc_annot_data(dataset) %>%
        distinct(cell_type, mc_col) %>%
        rename(group = cell_type, color = mc_col) %>%
        as.data.frame()
}

get_mc_data <- function(dataset, var_name) {
    if (var_name == "mc_annot") {
        return(get_mc_annot_data(dataset))
    } else if (var_name == "cell_type_annot") {
        return(get_cell_type_data(dataset))
    }

    mc_data[[dataset]][[var_name]]
}

get_mc_config <- function(dataset, var_name) {
    config$metacells[[dataset]][[var_name]]
}

dataset_ls <- function() {
    names(config$metacells)
}

has_network <- function(dataset) {
    !is.null(get_mc_config(dataset, "network"))
}

has_time <- function(dataset) {
    !is.null(get_mc_config(dataset, "time_bin_field")) && has_name(get_mc_config(dataset, "mc_annot"), "mc_age")
}
