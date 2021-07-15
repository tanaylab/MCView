#' Import a dataset to an MCView project
#'
#' This would read would read an \code{anndata} file which is the output of
#' python \code{metacells} package, and import the metacell dataset to MCView.
#' The result would be a directory under \code{project/cache/dataset} which
#' would contain objects used by MCView shiny app (such as the metacell matrix).
#' In addition, you can supply file with type assignment for each metacell
#' (\code{metacell_types_file}) and a file with color assignment for each metacell type
#' (\code{cell_type_colors_file}).
#'
#'
#' @param project path to the project
#' @param dataset name for the dataset, e.g. "PBMC"
#' @param anndata_file path to \code{h5ad} file which contains the output of metacell2 pipeline (metacells python package).
#' @param cell_type_field name of a field in the anndata object$obs which contains a cell type (optional).
#' If this parameter and \code{metacell_types_file} are missing, MCView would cluster the
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
#' download.file("http://www.wisdom.weizmann.ac.il/~atanay/metac_data/PBMC_processed.tar.gz", "raw/#' PBMC_processed.tar.gz")
#' untar("raw/PBMC_processed.tar.gz", exdir = "raw")
#' create_project("PBMC")
#' import_dataset("PBMC", "PBMC163k", "raw/metacells.h5ad")
#' }
#'
#' @export
import_dataset <- function(project, dataset, anndata_file, cell_type_field = NULL, metacell_types_file = NULL, cell_type_colors_file = NULL) {
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
        if (!is.null(cell_type_field)) {
            cli_alert_info("Taking cell type annotations from {.field {cell_type_field}} field in the anndata object")
            metacell_types <- adata$obs %>%
                select(cell_type = !!cell_type_field) %>%
                rownames_to_column("metacell") %>%
                as_tibble()
        } else {
            cli_alert_info("Clustering in order to get initial annotation.")
            # we generate clustering as initial annotation
            feat_mat <- mc_egc[adata$var$top_feature_gene, ]
            km <- cluster_egc(feat_mat, verbose = getOption("MCView.verbose"))
            metacell_types <- km$clusters %>%
                rename(cell_type = cluster) %>%
                mutate(cell_type = as.character(cell_type))
        }
    }

    if (!is.null(cell_type_colors_file)) {
        cli_alert_info("Loading metacell type annotations from {.file {cell_type_colors_file}}")
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
#' 
#' \dontrun{
#' dataset_rm("PBMC", "PBMC163k")
#' } 
#' 
#' @export
dataset_rm <- function(project, dataset){
    fs::dir_delete(fs::path(project_cache_dir(project), dataset))
}

#' List datasets in a project
#' 
#' @project path to the project directory
#' 
#' @return names of the existing datasets in the project
#' 
#' @examples 
#' 
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
#' This would change the cell type assignments for each metacell to the ones listed at \code{metacell_types_file}. 
#' This is usually done after a first iteration of annotation using the "Annotate" tab in the MCView annotation, which can 
#' export a valid \code{metacell_types_file}. 
#' The file should have a column named "metacell" with the metacell ids and another
#' column named "cell_type" or "cluster" with the cell type assignment.
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
#' 
#' \dontrun{
#'  update_metacell_types("PBMC", "PBMC163k", "raw/metacell-clusters.csv")
#' }
#' 
#' @export
update_metacell_types <- function(project, dataset, metacell_types_file){
    verify_project_dir(project)
    verify_app_cache(project)

    prev_metacell_types <- load_shiny_data("metacell_types", dataset, project_cache_dir(project))

    metacell_types <- parse_metacell_types(metacell_types_file)

    metacell_types <- metacell_types %>% 
        select(metacell, cell_type) %>% 
        left_join(prev_metacell_types %>% select(-cell_type), by = "metacell")

    serialize_shiny_data(metacell_types, "metacell_types", dataset = dataset, cache_dir = project_cache_dir(project), flat = TRUE)
    
    cli_alert_success("Succesfully changed metacell cell type assignments")
}


#' Update color assignment for each cell type
#' 
#' This would change the color assignments for each cell type to the ones listed at \code{cell_type_colors_file}. 
#' This is usually done after a first iteration of annotation using the "Annotate" tab in the MCView annotation, which can 
#' export a valid \code{cell_type_colors_file}. 
#' The file should have a column named "cell_type" or "cluster" with the cell types and another column named "color" with the color assignment. 
#' Under the hood - MCView updates a file named "cell_type_colors.tsv" under \code{project/cache/dataset}, which can also be edited manually.
#' 
#' 
#' @param cell_type_colors_file path to a tabular file (csv,tsv) with color assignement for
#' each cell type. The file should have a column named "cell_type" or "cluster" with the
#' cell types and another column named "color" with the color assignment. Cell types that do not
#' exist in the metacell types would be ignored, so if you changed the names of cell types you would have to also 
#' update the metacell types (using \code{update_metacell_types})
#' If this parameter is missing, MCView would use the \code{chameleon} package to assign a color for each cell type.
#' 
#' 
#' @examples
#' 
#' \dontrun{
#'  update_metacell_types("PBMC", "PBMC163k", "raw/cluster-colors.csv")
#' }
#' 
#'
#' @export
update_cell_type_colors <- function(project, dataset, cell_type_colors_file){
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
        select(any_of(c("metacell", "cell_type", "age")))

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


# #' Initialize from metacell R package
# #'
# #' @param config list with config parameters
# #' @param dataset name of the dataset
# #' @param cache_dir path of the data directory
# #'
# #' @noRd
# init_metacell <- function(config, dataset, cache_dir) {
#     fs::dir_create(cache_dir, dataset, recurse = TRUE)
#     if (!is.null(config$anndata)) {
#         return(init_metacell2(config, dataset, cache_dir))
#     }

#     library(metacell)

#     init_temp_scdb(config, dataset)
#     mc <- scdb_mc(config$mc)

#     mc_egc <- mc@e_gc

#     mc_fp <- mc@mc_fp

#     mat <- scdb_mat(config$matrix)

#     mc_mat <- tgs_matrix_tapply(mat@mat[, names(mc@mc)], mc@mc, sum, na.rm = TRUE) %>% t()
#     serialize_shiny_data(mc_mat, "mc_mat", dataset = dataset, cache_dir = cache_dir)

#     mc_sum <- colSums(mc_mat)
#     serialize_shiny_data(mc_sum, "mc_sum", dataset = dataset, cache_dir = cache_dir)

#     mc2d <- scdb_mc2d(config$mc2d)

#     if (is.null(names(mc2d@mc_x))) {
#         names(mc2d@mc_x) <- 1:length(mc2d@mc_x)
#         warning(glue("mc2d@mc_x doesn't have names. Setting to 1:{length(mc2d@mc_x)})"))
#     }

#     if (is.null(names(mc2d@mc_y))) {
#         names(mc2d@mc_y) <- 1:length(mc2d@mc_y)
#         warning(glue("mc2d@mc_y doesn't have names. Setting to 1:{length(mc2d@mc_y)})"))
#     }


#     mc2d_list <- list(
#         graph = mc2d@graph,
#         mc_id = mc2d@mc_id,
#         mc_x = mc2d@mc_x,
#         mc_y = mc2d@mc_y
#     )
#     serialize_shiny_data(mc2d_list, "mc2d", dataset = dataset, cache_dir = cache_dir)


#     mc_genes_top2 <- apply(mc@mc_fp, 2, function(fp) {
#         top_ind <- order(-fp)[1:2]
#         return(rownames(mc@mc_fp)[top_ind])
#     })

#     mc_genes_top2 <- mc_genes_top2 %>%
#         t() %>%
#         as.data.frame() %>%
#         rownames_to_column("metacell") %>%
#         rlang::set_names(c("metacell", "top1_gene", "top2_gene")) %>%
#         tibble::remove_rownames() %>%
#         distinct(metacell, .keep_all = TRUE) %>%
#         mutate(metacell = as.character(metacell))

#     metacell_types<- parse_metacell_types(metacell_types_file)

#     metacell_types<- metacell_types%>%
#         arrange(as.numeric(metacell)) %>%
#         left_join(mc_genes_top2, by = "metacell")

#     # filter metacell_typeto contain only metacells that are within mc_egc
#     metacell_types<- metacell_types%>%
#         as.data.frame() %>%
#         column_to_rownames("metacell")
#     metacell_types<- metacell_types[colnames(mc_egc), ]
#     metacell_types<- metacell_types%>%
#         rownames_to_column("metacell") %>%
#         as_tibble()

#     if (!is.null(cell_type_colors_file)) {
#         cell_type_colors<- parse_cell_type_colors(cell_type_colors_file)
#     } else {
#         cell_type_colors<- metacell_types%>%
#             distinct(cell_type, mc_col) %>%
#             rename(color = mc_col) %>%
#             arrange(cell_type) %>%
#             mutate(order = 1:n())
#     }

#     cell_type_colors<- cell_type_colors%>%
#         arrange(as.numeric(order)) %>%
#         mutate(cell_type = factor(cell_type), cell_type = forcats::fct_inorder(cell_type))
#     serialize_shiny_data(cell_type_colors, "cell_type_colors", dataset = dataset, cache_dir = cache_dir)

#     metacell_types$top1_lfp <- purrr::map2_dbl(metacell_types$metacell, metacell_types$top1_gene, ~ log2(mc_fp[.y, .x]))
#     metacell_types$top2_lfp <- purrr::map2_dbl(metacell_types$metacell, metacell_types$top2_gene, ~ log2(mc_fp[.y, .x]))
#     serialize_shiny_data(metacell_types, "metacell_types", dataset = dataset, cache_dir = cache_dir)

#     # Top 30 correlated and anti-correlated genes for each gene
#     gg_mc_top_cor <- calc_gg_mc_top_cor(mc@e_gc, k = 30)

#     serialize_shiny_data(gg_mc_top_cor, "gg_mc_top_cor", dataset = dataset, cache_dir = cache_dir)

#     if (!is.null(config$time_bin_field)) {
#         mc_ag <- table(mc@mc, mat@cell_metadata[names(mc@mc), config$time_bin_field])
#         mc_ag_n <- t(t(mc_ag) / colSums(mc_ag))
#         serialize_shiny_data(mc_ag, "mc_ag", dataset = dataset, cache_dir = cache_dir)
#         serialize_shiny_data(mc_ag_n, "mc_ag_n", dataset = dataset, cache_dir = cache_dir)
#     }

#     if (!is.null(config$time_annot) && !is.null(config$time_bin_field)) {
#         time_annot <- fread(config$time_annot$fn) %>% as_tibble()
#         serialize_shiny_data(time_annot, "time_annot", dataset = dataset, cache_dir = cache_dir)

#         cell_md <- mat@cell_metadata[names(mc@mc), ] %>%
#             rownames_to_column("cell_id") %>%
#             mutate(metacell = as.character(mc@mc)) %>%
#             left_join(metacell_types)

#         if (config$time_bin_field != "time_bin") {
#             if (rlang::has_name(cell_md, "time_bin")) {
#                 cell_md <- cell_md %>%
#                     select(-time_bin)
#             }
#             cell_md <- cell_md %>%
#                 rename(time_bin = !!config$time_bin_field) %>%
#                 as_tibble()
#         }

#         min_time_bin <- min(cell_md$time_bin, na.rm = TRUE)
#         max_time_bin <- max(cell_md$time_bin, na.rm = TRUE)
#         obs_time_bins <- sort(unique(as.numeric(cell_md$time_bin)))

#         # in case time bins are not from 1:max_t
#         if (min_time_bin != 1 || length(obs_time_bins) != length(1:max_time_bin) || !all(1:max_time_bin == obs_time_bins)) {
#             cell_md <- cell_md %>% mutate(time_bin = as.numeric(factor(time_bin, levels = obs_time_bins)))
#         }

#         cell_md <- cell_md %>%
#             left_join(time_annot, by = "time_bin")

#         type_ag <- cell_md %>%
#             count(cell_type, time_bin) %>%
#             filter(!is.na(time_bin)) %>%
#             spread(time_bin, n, fill = 0) %>%
#             as.data.frame() %>%
#             column_to_rownames("cell_type") %>%
#             as.matrix()
#         serialize_shiny_data(type_ag, "type_ag", dataset = dataset, df2mat = TRUE, cache_dir = cache_dir)

#         mc_ag <- cell_md %>%
#             count(metacell, time_bin) %>%
#             filter(!is.na(time_bin)) %>%
#             spread(time_bin, n, fill = 0) %>%
#             as.data.frame() %>%
#             column_to_rownames("metacell") %>%
#             as.matrix()
#         serialize_shiny_data(mc_ag, "mc_ag", dataset = dataset, df2mat = TRUE, cache_dir = cache_dir)
#     }

#     if (!is.null(config$network)) {
#         mc_network <- scdb_mctnetwork(config$network)
#         serialize_shiny_data(mc_network@network, "mc_network", dataset = dataset, cache_dir = cache_dir)

#         # calculate flows of every metacell
#         metacells <- as.character(sort(as.numeric(unique(mc@mc))))
#         future::plan(future::multicore, workers = min(20, future::availableCores(constraints = "multicore")))
#         options(future.globals.maxSize = 1e11)
#         mct_probs_trans <- furrr::future_map(1:length(metacells), ~ {
#             cli_alert_info("prop_flow: ({.x}/{length(metacells)}")
#             mct_propagate_flow_through_metacell(mct = mc_network, m = .x)
#         })
#         names(mct_probs_trans) <- metacells
#         serialize_shiny_data(mct_probs_trans, "mct_probs_trans", dataset = dataset, cache_dir = cache_dir)

#         # calculate order of metacells in the flow chart
#         # this order is static and will always be the same
#         # mc_rank <- mctnetwork_mc_rank_from_color_ord(config$network, mc@color_key$color)
#         mc_rank <- mctnetwork_mc_rank_from_color_ord(config$network, metacell_types$mc_col, cell_type_colors$color)
#         mc_rank["-2"] <- 0
#         mc_rank["-1"] <- length(mc_rank) / 2
#         serialize_shiny_data(mc_rank, "mc_rank", dataset = dataset, cache_dir = cache_dir)



#         type_flow <- mctnetwork_get_type_flows(mc_network, 1, ncol(type_ag))
#         serialize_shiny_data(type_flow, "type_flow", dataset = dataset, cache_dir = cache_dir)
#     }
# }
