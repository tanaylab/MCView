
verify_config_file <- function(config) {
    required_fields <- c(
        "title",
        "tabs",
        "help",
        "selected_mc1",
        "selected_mc2"
    )

    for (field in required_fields) {
        if (!(field %in% names(config))) {
            cli_abort("The field {.field {field}} does not exist in the config file.")
        }
    }

    # TODO: verify datasets config

    invisible(TRUE)
}


#' Initialize application config file
#'
#' @param project path of the project to initialize
#'
#' @noRd
init_config <- function(project) {
    config_file <- project_config_file(project)
    help_file <- project_help_file(project)
    project <<- project
    about_file <<- fs::path_abs(project_about_file(project))
    cache_dir <<- project_cache_dir(project)

    config <<- yaml::read_yaml(config_file)
    verify_config_file(config)

    verify_app_cache(project)

    help_config <<- yaml::read_yaml(help_file)
}

init_defs <- function() {
    egc_epsilon <<- 1e-5
    options(spinner.type = 6)
    gene_names <<- get_gene_names()

    theme_set(theme_classic())

    init_tab_defs()

    init_selected_genes()
}

init_selected_genes <- function() {
    # if selected genes are not set - choose them from the first dataset
    if (is.null(config$selected_gene1) || is.null(config$selected_gene2)) {
        mc_egc <- get_mc_egc(dataset_ls(project)[1])
        minmax <- matrixStats::rowMaxs(mc_egc, na.rm = TRUE) - matrixStats::rowMins(mc_egc, na.rm = TRUE)
        names(minmax) <- rownames(mc_egc)
        genes <- names(head(sort(minmax, decreasing = TRUE), n = 2))
    }

    default_gene1 <<- config$selected_gene1 %||% genes[1]
    default_gene2 <<- config$selected_gene2 %||% genes[2]
}


init_tab_defs <- function() {
    tab_defs <<- list(
        "About" = list(
            title = "About",
            module_name = "about",
            icon = "info"
        ),
        "Manifold" = list(
            title = "Manifold",
            module_name = "manifold",
            icon = "project-diagram"
        ),
        "Genes" = list(
            title = "Genes",
            module_name = "gene_mc",
            icon = "wind"
        ),
        "Metacells" = list(
            title = "Metacells",
            module_name = "mc_mc",
            icon = "bar-chart"
        ),
        "Annotate" = list(
            title = "Annotate",
            module_name = "annotate",
            icon = "sticky-note"
        )
    )

    default_tabs <- c("About", "Genes", "Metacells")

    if (!is.null(config$tabs)) {
        purrr::walk(config$tabs, ~ {
            if (!(.x %in% names(tab_defs))) {
                cli_abort("{.x} is not a valid tab name. Update `tabs` in your configuration file.")
            }
        })
        tab_defs <<- tab_defs[config$tabs]
    } else {
        tab_defs <<- tab_defs[default_tabs]
    }
}

get_gene_names <- function() {
    gene_names <- rownames(get_mc_data(dataset_ls(project)[1], "mc_mat"))
    # We remove gene names that are too long in order to reduce pickerInput search bar width
    gene_names <- gene_names[stringr::str_length(gene_names) <= 30]
    return(gene_names)
}


#' Access files in the current app
#'
#' @param ... Character vector specifying directory and or file to
#'     point to inside the current package.
#'
#' @noRd
app_sys <- function(...) {
    system.file(..., package = "MCView")
}


#' Read App Config
#'
#' @param value Value to retrieve from the config file.
#' @param config R_CONFIG_ACTIVE value.
#' @param use_parent Logical, scan the parent directory for config file.
#'
#' @importFrom config get
#'
#' @noRd
get_golem_config <- function(value,
                             config = Sys.getenv("R_CONFIG_ACTIVE", "default"),
                             use_parent = TRUE) {
    config::get(
        value = value,
        config = config,
        # Modify this if your config file is somewhere else:
        file = app_sys("golem-config.yml"),
        use_parent = use_parent
    )
}