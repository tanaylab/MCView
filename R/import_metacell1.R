#' Import a dataset to an MCView project from metacell R package
#'
#' Read objects from \code{metacell} R package and import a
#' metacell dataset to MCView.
#'
#' The result would be a directory under \code{project/cache/dataset} which
#' would contain objects used by MCView shiny app (such as the metacell matrix).
#'
#' In addition, you can supply file with type assignment for each metacell
#' (\code{metacell_types_file}) and a file with color assignment for each metacell type
#' (\code{cell_type_colors_file}).
#'
#' Make sure that you have the R \code{metacell} package installed in order to use
#' this function.
#'
#' \code{network}, \code{time_annotation_file} and \code{time_bin_field} are only relevant
#' if you computed flows/networks for your dataset and therefore are optional.
#'
#' In order to add time annotation to your dataset you will have to:
#' \itemize{
#'  \item{1. }{Add a column named "mc_age" or "age" to \code{metacell_types_file} with time per metacell}
#'  \item{2. }{Create a \code{time_annotation_file} with id for each time bin and description}
#' }
#'
#'
#' @param scdb path to R metacell single cell RNA database
#' @param matrix  name of the umi matrix to use
#' @param mc  name of the metacell object to use
#' @param mc2d  name of the 2d projection object to use
#' @param metacell_types_file path to a tabular file (csv,tsv) with cell type assignement for
#' each metacell. The file should have a column named "metacell" with the metacell ids and another
#' column named "cell_type" or "cluster" with the cell type assignment. Metacell ids that do
#' not exists in the data would be ignored. In addition, the file can have a column named
#' "age" or "mc_age" with age metadata per metacell
#' @param cell_type_colors_file path to a tabular file (csv,tsv) with color assignement for
#' each cell type. The file should have a column named "cell_type" or "cluster" with the
#' cell types and another column named "color" with the color assignment. Cell types that do not
#' exist in the metacell types would be ignored.
#' @param network  name of the network object to use (optional)
#' @param time_annotation_file file with names for time bins (optional, only relevant with networks/flows). Should have a field named "time_bin" with the time bin id and another field named "time_desc" which contains the description of the time bin
#' @param time_bin_field name of a field in \code{cell_metadata} which contains time bin per cell (optional)
#' @param metadata_fields names of fields \code{mat@cell_metadata} which contains metadata per cell to be summarized using \code{cell_metadata_to_metacell}. \cr
#' The fields should can be either numeric or categorical. \cr
#'  You can use \code{cell_metadata_to_metacell} to convert from categorical to a numeric score (e.g. by using fraction of the category).
#' @param categorical metadata fields that should be treated as categorical (optional)
#'
#' @examples
#' \dontrun{
#' import_dataset_metacell1(
#'     "embflow",
#'     "153embs",
#'     scdb = "raw/scrna_db",
#'     matrix = "embs",
#'     mc = "embs",
#'     mc2d = "embs",
#'     metacell_types_file = "raw/metacell-types.csv",
#'     cell_type_colors_file = "raw/cell-type-colors.csv",
#'     network = "embs",
#'     time_annotation_file = "raw/time-annot.tsv",
#'     time_bin_field = "age_group"
#' )
#' }
#'
#' @inheritParams import_dataset
#' @inheritDotParams create_project
#'
#' @export
import_dataset_metacell1 <- function(project,
                                     dataset,
                                     scdb,
                                     matrix,
                                     mc,
                                     mc2d,
                                     metacell_types_file,
                                     cell_type_colors_file,
                                     gene_modules_file = NULL,
                                     gene_modules_k = NULL,
                                     calc_gg_cor = TRUE,
                                     network = NULL,
                                     time_annotation_file = NULL,
                                     time_bin_field = NULL,
                                     metadata_fields = NULL,
                                     categorical = c(),
                                     ...) {
    if (!requireNamespace("metacell", quietly = TRUE)) {
        stop("Please install metacell R package in order to use this function")
    }

    verbose <- !is.null(getOption("MCView.verbose")) && getOption("MCView.verbose")
    init_project_dir(project, create = TRUE, ...)

    cli_alert_info("Importing {.field {dataset}}")

    cache_dir <- project_cache_dir(project)

    library(metacell)

    mc_name <- mc
    init_temp_scdb(scdb, matrix, mc_name, mc2d, network, dataset)
    mc <- scdb_mc(mc_name)

    mc_egc <- mc@e_gc

    mc_fp <- mc@mc_fp

    mat <- scdb_mat(matrix)

    mc_mat <- tgs_matrix_tapply(mat@mat[, names(mc@mc)], mc@mc, sum, na.rm = TRUE) %>% t()
    serialize_shiny_data(mc_mat, "mc_mat", dataset = dataset, cache_dir = cache_dir)

    metacells <- colnames(mc_mat)

    mc_sum <- colSums(mc_mat)
    serialize_shiny_data(mc_sum, "mc_sum", dataset = dataset, cache_dir = cache_dir)

    mc2d <- scdb_mc2d(mc2d)

    if (is.null(names(mc2d@mc_x))) {
        names(mc2d@mc_x) <- 1:length(mc2d@mc_x)
        warning(glue("mc2d@mc_x doesn't have names. Setting to 1:{length(mc2d@mc_x)})"))
    }

    if (is.null(names(mc2d@mc_y))) {
        names(mc2d@mc_y) <- 1:length(mc2d@mc_y)
        warning(glue("mc2d@mc_y doesn't have names. Setting to 1:{length(mc2d@mc_y)})"))
    }

    mc2d_list <- list(
        graph = mc2d@graph,
        mc_id = mc2d@mc_id,
        mc_x = mc2d@mc_x,
        mc_y = mc2d@mc_y
    )
    serialize_shiny_data(mc2d_list, "mc2d", dataset = dataset, cache_dir = cache_dir)

    cli_alert_info("Calculating top genes per metacell (marker genes)")
    marker_genes <- calc_marker_genes(mc_egc, 20)
    serialize_shiny_data(marker_genes, "marker_genes", dataset = dataset, cache_dir = cache_dir)

    mc_genes_top2 <- marker_genes %>%
        group_by(metacell) %>%
        slice(1:2) %>%
        ungroup() %>%
        pivot_wider(names_from = "rank", values_from = c("gene", "fp"))

    colnames(mc_genes_top2) <- c("metacell", "top1_gene", "top2_gene", "top1_lfp", "top2_lfp")


    cli_alert_info("Loading metacell type annotations from {.file {metacell_types_file}}")
    metacell_types <- parse_metacell_types(metacell_types_file, metacells)


    cli_alert_info("Loading cell type color annotations from {.file {cell_type_colors_file}}")
    cell_type_colors <- parse_cell_type_colors(cell_type_colors_file)

    cell_type_colors <- cell_type_colors %>%
        arrange(as.numeric(order)) %>%
        mutate(cell_type = factor(cell_type), cell_type = forcats::fct_inorder(cell_type))
    serialize_shiny_data(cell_type_colors, "cell_type_colors", dataset = dataset, cache_dir = cache_dir, flat = TRUE)

    metacell_types <- metacell_types %>%
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

    if (!is.null(metadata_fields)) {
        metadata <- cell_metadata_to_metacell_from_metacell1(
            scdb = scdb,
            matrix = matrix,
            mc = mc_name,
            metadata_fields = metadata_fields,
            categorical = categorical
        )
        update_metadata(project, dataset, metadata)
    }

    if (calc_gg_cor) {
        cli_alert_info("Calculating top 30 correlated and anti-correlated genes for each gene")
        # Top 30 correlated and anti-correlated genes for each gene
        gg_mc_top_cor <- calc_gg_mc_top_cor(mc@e_gc, k = 30)
        serialize_shiny_data(gg_mc_top_cor, "gg_mc_top_cor", dataset = dataset, cache_dir = cache_dir)
    } else {
        cli_alert_info("Skipping calculation of top 30 correlated and anti-correlated genes for each gene. Some features in the app would not be available")
    }

    if (!is.null(gene_modules_file)) {
        gene_modules <- parse_gene_modules_file(gene_modules_file)
    } else {
        gene_modules <- calc_gene_modules(mc_mat, k = gene_modules_k)
    }
    serialize_shiny_data(gene_modules, "gene_modules", dataset = dataset, cache_dir = cache_dir, flat = TRUE)

    if (!is.null(time_bin_field)) {
        mc_ag <- table(mc@mc, mat@cell_metadata[names(mc@mc), time_bin_field])
        mc_ag_n <- t(t(mc_ag) / colSums(mc_ag))
        serialize_shiny_data(mc_ag, "mc_ag", dataset = dataset, cache_dir = cache_dir)
        serialize_shiny_data(mc_ag_n, "mc_ag_n", dataset = dataset, cache_dir = cache_dir)
    }

    if (!is.null(time_annotation_file)) {
        time_annot <- fread(time_annotation_file) %>% as_tibble()
        serialize_shiny_data(time_annot, "time_annot", dataset = dataset, cache_dir = cache_dir)

        if (!is.null(time_bin_field)) {
            cell_md <- mat@cell_metadata[names(mc@mc), ] %>%
                rownames_to_column("cell_id") %>%
                mutate(metacell = as.character(mc@mc)) %>%
                left_join(metacell_types)

            if (time_bin_field != "time_bin") {
                if (rlang::has_name(cell_md, "time_bin")) {
                    cell_md <- cell_md %>%
                        select(-time_bin)
                }
                cell_md <- cell_md %>%
                    rename(time_bin = !!time_bin_field) %>%
                    as_tibble()
            }

            min_time_bin <- min(cell_md$time_bin, na.rm = TRUE)
            max_time_bin <- max(cell_md$time_bin, na.rm = TRUE)
            obs_time_bins <- sort(unique(as.numeric(cell_md$time_bin)))

            # in case time bins are not from 1:max_t
            if (min_time_bin != 1 || length(obs_time_bins) != length(1:max_time_bin) || !all(1:max_time_bin == obs_time_bins)) {
                cell_md <- cell_md %>% mutate(time_bin = as.numeric(factor(time_bin, levels = obs_time_bins)))
            }

            cell_md <- cell_md %>%
                left_join(time_annot, by = "time_bin")

            type_ag <- cell_md %>%
                count(cell_type, time_bin) %>%
                filter(!is.na(time_bin)) %>%
                spread(time_bin, n, fill = 0) %>%
                as.data.frame() %>%
                column_to_rownames("cell_type") %>%
                as.matrix()
            serialize_shiny_data(type_ag, "type_ag", dataset = dataset, df2mat = TRUE, cache_dir = cache_dir)

            mc_ag <- cell_md %>%
                count(metacell, time_bin) %>%
                filter(!is.na(time_bin)) %>%
                spread(time_bin, n, fill = 0) %>%
                as.data.frame() %>%
                column_to_rownames("metacell") %>%
                as.matrix()
            serialize_shiny_data(mc_ag, "mc_ag", dataset = dataset, df2mat = TRUE, cache_dir = cache_dir)
        }
    }



    if (!is.null(network)) {
        mc_network <- scdb_mctnetwork(network)
        serialize_shiny_data(mc_network@network, "mc_network", dataset = dataset, cache_dir = cache_dir)

        # calculate flows of every metacell
        metacells <- as.character(sort(as.numeric(unique(mc@mc))))
        future::plan(future::multicore, workers = min(20, future::availableCores(constraints = "multicore")))
        options(future.globals.maxSize = 1e11)
        mct_probs_trans <- furrr::future_map(1:length(metacells), ~ {
            cli_alert_info("prop_flow: ({.x}/{length(metacells)}")
            mct_propagate_flow_through_metacell(mct = mc_network, m = .x)
        })
        names(mct_probs_trans) <- metacells
        serialize_shiny_data(mct_probs_trans, "mct_probs_trans", dataset = dataset, cache_dir = cache_dir)

        # calculate order of metacells in the flow chart
        # this order is static and will always be the same
        mc_rank <- mctnetwork_mc_rank_from_color_ord(
            network,
            metacell_types %>% left_join(cell_type_colors, by = "cell_type") %>% pull(color),
            cell_type_colors$color
        )
        mc_rank["-2"] <- 0
        mc_rank["-1"] <- length(mc_rank) / 2
        serialize_shiny_data(mc_rank, "mc_rank", dataset = dataset, cache_dir = cache_dir)



        type_flow <- mctnetwork_get_type_flows(mc_network, 1, ncol(type_ag))
        serialize_shiny_data(type_flow, "type_flow", dataset = dataset, cache_dir = cache_dir)
    }
}
