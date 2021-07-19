#' Import a dataset to an MCView project
#'
#' Read an \code{anndata} file which is the output of python \code{metacells} package, 
#' and import the metacell dataset to MCView. Each project can have multiple datasets
#' which can be in the app using the right sidebar. 
#' 
#' The function would create a directory under \code{project/cache/dataset} which
#' would contain objects used by MCView shiny app (such as the metacell matrix).
#' 
#' In addition, you can supply file with type assignment for each metacell
#' (\code{metacell_types_file}) and a file with color assignment for each metacell type
#' (\code{cell_type_colors_file}).
#' 
#'
#'
#' @param project path to the project
#' @param dataset name for the dataset, e.g. "PBMC"
#' @param anndata_file path to \code{h5ad} file which contains the output of metacell2 pipeline (metacells python package).
#' @param cell_type_field name of a field in the anndata object$obs which contains a cell type (optional).
#' If the field doesn't exist and \code{metacell_types_file} are missing, MCView would cluster the
#' metacell matrix using kmeans++ algorithm (from the \code{tglkmeans} package).
#' If \code{metacell_types} parameter is set this field is ignored.
#' @param metacell_types_file path to a tabular file (csv,tsv) with cell type assignement for
#' each metacell. The file should have a column named "metacell" with the metacell ids and another
#' column named "cell_type" or "cluster" with the cell type assignment. Metacell ids that do
#' not exists in the data would be ignored.
#' If this parameter and \code{cell_type_field} are missing, MCView would cluster the
#' metacell matrix using kmeans++ algorithm (from the \code{tglkmeans} package).
#' @param cell_type_colors_file path to a tabular file (csv,tsv) with color assignement for
#' each cell type. The file should have a column named "cell_type" or "cluster" with the
#' cell types and another column named "color" with the color assignment. Cell types that do not
#' exist in the metacell types would be ignored.
#' If this is missing, MCView would use the \code{chameleon} package to assign a color for each cell type.
#'
#'
#' @examples
#' \dontrun{
#' dir.create("raw")
#' download.file(
#'      "http://www.wisdom.weizmann.ac.il/~atanay/metac_data/PBMC_processed.tar.gz", 
#'      "raw/PBMC_processed.tar.gz"
#' )
#' untar("raw/PBMC_processed.tar.gz", exdir = "raw")
#' create_project("PBMC")
#' import_dataset("PBMC", "PBMC163k", "raw/metacells.h5ad")
#' }
#'
#' @export
import_dataset <- function(project, dataset, anndata_file, cell_type_field = "cluster", metacell_types_file = NULL, cell_type_colors_file = NULL) {
    verbose <- !is.null(getOption("MCView.verbose")) && getOption("MCView.verbose")
    verify_project_dir(project)

    cli_alert_info("Importing {.field {dataset}}")

    cache_dir <- project_cache_dir(project)

    if (!fs::file_exists(anndata_file)) {
        cli_abort("{anndata_file} doesn't exist. Maybe there is a typo?")
    }

    cli_alert_info("Reading {.file {anndata_file}}")
    library(anndata)
    adata <- anndata::read_h5ad(anndata_file)

    cli_alert_info("Processing metacell matrix")
    mc_mat <- t(adata$X)
    serialize_shiny_data(mc_mat, "mc_mat", dataset = dataset, cache_dir = cache_dir)

    mc_sum <- colSums(mc_mat)
    serialize_shiny_data(mc_sum, "mc_sum", dataset = dataset, cache_dir = cache_dir)

    cli_alert_info("Processing 2d projection")
    graph <- Matrix::summary(adata$obsp$obs_outgoing_weights) %>%
        as.data.frame()

    graph <- graph %>% mutate(i = rownames(adata$obs)[i], j = rownames(adata$obs)[j])

    mc2d_list <- list(
        graph = tibble(mc1 = graph[, 1], mc2 = graph[, 2], weight = graph[, 3]),
        mc_id = rownames(adata$obs),
        mc_x = adata$obs %>% select(umap_x) %>% tibble::rownames_to_column("mc") %>% tibble::deframe(),
        mc_y = adata$obs %>% select(umap_y) %>% tibble::rownames_to_column("mc") %>% tibble::deframe()
    )
    serialize_shiny_data(mc2d_list, "mc2d", dataset = dataset, cache_dir = cache_dir)

    cli_alert_info("Calculating top genes per metacell")
    mc_egc <- t(t(mc_mat) / mc_sum)

    mc_egc_norm <- mc_egc + 1e-5
    mc_fp <- mc_egc_norm / apply(mc_egc_norm, 1, median, na.rm = TRUE)

    mc_genes_top2 <- apply(mc_fp, 2, function(fp) {
        top_ind <- order(-fp)[1:2]
        return(rownames(mc_fp)[top_ind])
    })

    mc_genes_top2 <- mc_genes_top2 %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("metacell") %>%
        rlang::set_names(c("metacell", "top1_gene", "top2_gene")) %>%
        tibble::remove_rownames() %>%
        distinct(metacell, .keep_all = TRUE) %>%
        mutate(metacell = as.character(metacell))


    if (!is.null(metacell_types_file)) {
        if (!is.null(cell_type_field)) {
            cli_alert_warning("{.field cell_type_field} was ignored since {.field metacell_types_file} was set.")
        }
        cli_alert_info("Loading metacell type annotations from {.file {metacell_types_file}}")
        metacell_types <- parse_metacell_types(metacell_types_file)
    } else {
        if (!is.null(cell_type_field) && !is.null(adata$obs[[cell_type_field]])) {
            cli_alert_info("Taking cell type annotations from {.field {cell_type_field}} field in the anndata object")
            metacell_types <- adata$obs %>%
                select(cell_type = !!cell_type_field) %>%
                rownames_to_column("metacell") %>%
                as_tibble()
        } else {
            cli_alert_info("Clustering in order to get initial annotation.")
            # we generate clustering as initial annotation
            feat_mat <- mc_egc[adata$var$top_feature_gene, ]
            km <- cluster_egc(feat_mat, verbose = verbose)
            metacell_types <- km$clusters %>%
                rename(cell_type = cluster) %>%
                mutate(cell_type = as.character(cell_type))
        }
    }

    if (!is.null(cell_type_colors_file)) {
        cli_alert_info("Loading cell type color annotations from {.file {cell_type_colors_file}}")
        cell_type_colors <- parse_cell_type_colors(cell_type_colors_file)
    } else {
        cli_alert_info("Generating cell type colors using {.pkg chameleon} package.")
        color_of_clusters <- chameleon::data_colors(t(mc_egc[adata$var$top_feature_gene, ]), groups = metacell_types$cell_type)

        cell_type_colors <- enframe(color_of_clusters, name = "cell_type", value = "color") %>%
            mutate(order = 1:n())
    }

    cell_type_colors <- cell_type_colors %>%
        arrange(as.numeric(order)) %>%
        mutate(cell_type = factor(cell_type), cell_type = forcats::fct_inorder(cell_type))
    serialize_shiny_data(cell_type_colors, "cell_type_colors", dataset = dataset, cache_dir = cache_dir, flat = TRUE)

    metacell_types <- metacell_types %>%
        arrange(as.numeric(metacell)) %>%
        left_join(mc_genes_top2, by = "metacell")

    # filter metacell_types to contain only metacells that are included in mc_egc
    metacell_types <- metacell_types %>%
        as.data.frame() %>%
        column_to_rownames("metacell")
    metacell_types <- metacell_types[colnames(mc_egc), ]
    metacell_types <- metacell_types %>%
        rownames_to_column("metacell") %>%
        as_tibble()

    # add the expression (log2) of top genes per metacell
    metacell_types$top1_lfp <- purrr::map2_dbl(metacell_types$metacell, metacell_types$top1_gene, ~ log2(mc_fp[.y, .x]))
    metacell_types$top2_lfp <- purrr::map2_dbl(metacell_types$metacell, metacell_types$top2_gene, ~ log2(mc_fp[.y, .x]))

    serialize_shiny_data(metacell_types, "metacell_types", dataset = dataset, cache_dir = cache_dir, flat = TRUE)


    cli_alert_info("Calculating top 30 correlated and anti-correlated genes for each gene")
    gg_mc_top_cor <- calc_gg_mc_top_cor(mc_egc, k = 30)

    serialize_shiny_data(gg_mc_top_cor, "gg_mc_top_cor", dataset = dataset, cache_dir = cache_dir)

    cli_alert_success("{.field {dataset}} dataset imported succesfully to {.path {project}} project")
}

#' Remove a dataset from a project
#'
#' @param project path to the project directory
#' @param dataset name of the dataset to remove
#'
#' @examples
#' \dontrun{
#' dataset_rm("PBMC", "PBMC163k")
#' }
#'
#' @export
dataset_rm <- function(project, dataset) {
    fs::dir_delete(fs::path(project_cache_dir(project), dataset))
}

#' List datasets in a project
#'
#' @param project path to the project directory
#'
#' @return names of the existing datasets in the project
#'
#' @examples
#' \dontrun{
#' dataset_ls("PBMC")
#' }
#'
#' @export
dataset_ls <- function(project) {
    basename(fs::dir_ls(project_cache_dir(project), type = c("directory", "symlink")))
}

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
    verify_project_dir(project)
    verify_app_cache(project)

    prev_metacell_types <- load_shiny_data("metacell_types", dataset, project_cache_dir(project))
    prev_metacell_types <- prev_metacell_types %>%
        mutate(metacell = as.character(metacell))

    metacell_types <- parse_metacell_types(metacell_types_file)

    metacell_types <- prev_metacell_types %>%
        select(-cell_type) %>%
        left_join(metacell_types %>% select(metacell, cell_type), by = "metacell")

    serialize_shiny_data(metacell_types, "metacell_types", dataset = dataset, cache_dir = project_cache_dir(project), flat = TRUE)

    cli_alert_success("Succesfully changed metacell cell type assignments")
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
#' update the metacell types (using \code{update_metacell_types})
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
    verify_project_dir(project)
    verify_app_cache(project)

    cell_type_colors <- parse_cell_type_colors(cell_type_colors_file)

    serialize_shiny_data(cell_type_colors, "cell_type_colors", dataset = dataset, cache_dir = project_cache_dir(project), flat = TRUE)

    cli_alert_success("Succesfully changed cell type color assignments")
}


parse_cell_type_colors <- function(file) {
    cell_type_colors <- fread(file) %>% as_tibble()

    if (!has_name(cell_type_colors, "cell_type") && !has_name(cell_type_colors, "cluster")) {
        cli_abort("{.field {file}} should have a column named {.field cell_type} or {.field cluster}")
    }

    if (!has_name(cell_type_colors, "color")) {
        cli_abort("{.field {file}} should have a column named {.field color}")
    }

    if (rlang::has_name(cell_type_colors, "cluster")) {
        cell_type_colors <- cell_type_colors %>% rename(cell_type = cluster)
    }

    if (!has_name(cell_type_colors, "order")) {
        cell_type_colors <- cell_type_colors %>% mutate(order = 1:n())
    }

    cell_type_colors <- cell_type_colors %>%
        select(cell_type, color, order)

    return(cell_type_colors)
}


parse_metacell_types <- function(file) {
    metacell_types <- fread(file) %>% as_tibble()

    if (!has_name(metacell_types, "metacell")) {
        cli_abort("{.field {file}} should have a column named {.field metacell}")
    }

    if (!has_name(metacell_types, "cell_type") && !has_name(cell_type_colors, "cluster")) {
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

    return(metacell_types)
}

#' calculate the k top correlated and anti correlated genes for each gene
#'
#' @param egc egc matrix (normalized metacell counts per gene)
#'
#' @noRd
calc_gg_mc_top_cor <- function(egc, k = 30, egc_epsilon = 1e-5) {
    lfp <- log2(egc + egc_epsilon)

    cm <- tgs_cor(t(lfp), y = t(lfp), pairwise.complete.obs = TRUE, spearman = FALSE)

    top_cor <- tgs_knn(cm, knn = k, diag = FALSE) %>%
        rename(gene1 = col1, gene2 = col2, cor = val) %>%
        mutate(type = "pos")
    anti_cor <- tgs_knn(-cm, knn = k, diag = FALSE) %>%
        rename(gene1 = col1, gene2 = col2, cor = val) %>%
        mutate(cor = -cor) %>%
        mutate(type = "neg")

    gg_mc_top_cor <- bind_rows(top_cor, anti_cor)

    gg_mc_top_cor <- gg_mc_top_cor %>%
        filter(!is.na(cor))

    return(gg_mc_top_cor)
}

cli_alert_verbose <- function(...) {
    verbose <- !is.null(getOption("MCView.verbose")) && getOption("MCView.verbose")
    if (verbose) {
        cli_alert_info(...)
    }
}

cli_alert_success_verbose <- function(...) {
    verbose <- !is.null(getOption("MCView.verbose")) && getOption("MCView.verbose")
    if (verbose) {
        cli_alert_success(...)
    }
}
