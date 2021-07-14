#' Import a project to MCView
#'
#' This would read the project \code{config.yaml} file and would import the metacell datasets listed
#' within it to MCView. The result would be a directory for each dataset under \code{project/cache} which
#' would contain objects used by MCView shiny app.
#'
#'
#' @param project path to the project
#' 
#' @examples
#'\dontrun{
#' dir.create("raw")
#' download.file("http://www.wisdom.weizmann.ac.il/~atanay/metac_data/PBMC_processed.tar.gz", "raw/#' PBMC_processed.tar.gz")
#' untar("raw/PBMC_processed.tar.gz", exdir = "raw")
#' create_project("PBMC")
#' import("PBMC")
#' }
#'
#' @export
import <- function(project) {
    verify_project_dir(project)

    config <- yaml::read_yaml(project_config_file(project))

    verify_config_file(config)

    for (dataset in names(config$metacells)) {
        cli_alert_info("Importing {dataset}")
        mc_config <- config$metacells[[dataset]]

        cache_dir <- project_cache_dir(project)
        init_metacell(mc_config, dataset, cache_dir)
    }
}


#' Initialize from metacell R package
#'
#' @param config list with config parameters
#' @param dataset name of the dataset
#' @param cache_dir path of the data directory
#'
#' @noRd
init_metacell <- function(config, dataset, cache_dir) {
    fs::dir_create(cache_dir, dataset, recurse = TRUE)
    if (!is.null(config$anndata)) {
        return(init_metacell2(config, dataset, cache_dir))
    }

    library(metacell)

    init_temp_scdb(config, dataset)
    mc <- scdb_mc(config$mc)

    mc_egc <- mc@e_gc

    mc_fp <- mc@mc_fp

    mat <- scdb_mat(config$matrix)

    mc_mat <- tgs_matrix_tapply(mat@mat[, names(mc@mc)], mc@mc, sum, na.rm = TRUE) %>% t()
    serialize_shiny_data(mc_mat, "mc_mat", dataset = dataset, cache_dir = cache_dir)

    mc_sum <- colSums(mc_mat)
    serialize_shiny_data(mc_sum, "mc_sum", dataset = dataset, cache_dir = cache_dir)

    mc2d <- scdb_mc2d(config$mc2d)

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


    mc_genes_top2 <- apply(mc@mc_fp, 2, function(fp) {
        top_ind <- order(-fp)[1:2]
        return(rownames(mc@mc_fp)[top_ind])
    })

    mc_genes_top2 <- mc_genes_top2 %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("metacell") %>%
        rlang::set_names(c("metacell", "top1_gene", "top2_gene")) %>%
        tibble::remove_rownames() %>%
        distinct(metacell, .keep_all = TRUE) %>%
        mutate(metacell = as.character(metacell))

    mc_annot <- parse_mc_annot(config$mc_annot)

    mc_annot <- mc_annot %>%
        arrange(as.numeric(metacell)) %>%
        left_join(mc_genes_top2, by = "metacell")

    # filter mc_annot to contain only metacells that are within mc_egc
    mc_annot <- mc_annot %>%
        as.data.frame() %>%
        column_to_rownames("metacell")
    mc_annot <- mc_annot[colnames(mc_egc), ]
    mc_annot <- mc_annot %>%
        rownames_to_column("metacell") %>%
        as_tibble()

    if (!is.null(config$cell_type_annot)) {
        cell_type_annot <- parse_cell_type_annot(config$cell_type_annot)
    } else {
        cell_type_annot <- mc_annot %>%
            distinct(cell_type, mc_col) %>%
            rename(color = mc_col) %>%
            arrange(cell_type) %>%
            mutate(order = 1:n())
    }

    cell_type_annot <- cell_type_annot %>%
        arrange(as.numeric(order)) %>%
        mutate(cell_type = factor(cell_type), cell_type = forcats::fct_inorder(cell_type))
    serialize_shiny_data(cell_type_annot, "cell_type_annot", dataset = dataset, cache_dir = cache_dir)

    mc_annot$top1_lfp <- purrr::map2_dbl(mc_annot$metacell, mc_annot$top1_gene, ~ log2(mc_fp[.y, .x]))
    mc_annot$top2_lfp <- purrr::map2_dbl(mc_annot$metacell, mc_annot$top2_gene, ~ log2(mc_fp[.y, .x]))
    serialize_shiny_data(mc_annot, "mc_annot", dataset = dataset, cache_dir = cache_dir)

    # Top 30 correlated and anti-correlated genes for each gene
    gg_mc_top_cor <- calc_gg_mc_top_cor(mc@e_gc, k = 30)

    serialize_shiny_data(gg_mc_top_cor, "gg_mc_top_cor", dataset = dataset, cache_dir = cache_dir)

    if (!is.null(config$time_bin_field)) {
        mc_ag <- table(mc@mc, mat@cell_metadata[names(mc@mc), config$time_bin_field])
        mc_ag_n <- t(t(mc_ag) / colSums(mc_ag))
        serialize_shiny_data(mc_ag, "mc_ag", dataset = dataset, cache_dir = cache_dir)
        serialize_shiny_data(mc_ag_n, "mc_ag_n", dataset = dataset, cache_dir = cache_dir)
    }

    if (!is.null(config$time_annot) && !is.null(config$time_bin_field)) {
        time_annot <- fread(config$time_annot$fn) %>% as_tibble()
        serialize_shiny_data(time_annot, "time_annot", dataset = dataset, cache_dir = cache_dir)

        cell_md <- mat@cell_metadata[names(mc@mc), ] %>%
            rownames_to_column("cell_id") %>%
            mutate(metacell = as.character(mc@mc)) %>%
            left_join(mc_annot)

        if (config$time_bin_field != "time_bin") {
            if (rlang::has_name(cell_md, "time_bin")) {
                cell_md <- cell_md %>%
                    select(-time_bin)
            }
            cell_md <- cell_md %>%
                rename(time_bin = !!config$time_bin_field) %>%
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

    if (!is.null(config$network)) {
        mc_network <- scdb_mctnetwork(config$network)
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
        # mc_rank <- mctnetwork_mc_rank_from_color_ord(config$network, mc@color_key$color)
        mc_rank <- mctnetwork_mc_rank_from_color_ord(config$network, mc_annot$mc_col, cell_type_annot$color)
        mc_rank["-2"] <- 0
        mc_rank["-1"] <- length(mc_rank) / 2
        serialize_shiny_data(mc_rank, "mc_rank", dataset = dataset, cache_dir = cache_dir)



        type_flow <- mctnetwork_get_type_flows(mc_network, 1, ncol(type_ag))
        serialize_shiny_data(type_flow, "type_flow", dataset = dataset, cache_dir = cache_dir)
    }
}


#' Initialize from metacells python package objects
#'
#' @param config list with config parameters
#' @param dataset name of the dataset
#' @param cache_dir path of the data directory
#'
#' @noRd
init_metacell2 <- function(config, dataset, cache_dir) {
    library(anndata)
    adata <- anndata::read_h5ad(config$anndata)

    mc_mat <- t(adata$X)
    serialize_shiny_data(mc_mat, "mc_mat", dataset = dataset, cache_dir = cache_dir)

    mc_sum <- colSums(mc_mat)
    serialize_shiny_data(mc_sum, "mc_sum", dataset = dataset, cache_dir = cache_dir)

    mc_egc <- t(t(mc_mat) / mc_sum)

    mc_egc_norm <- mc_egc + 1e-5
    mc_fp <- mc_egc_norm / apply(mc_egc_norm, 1, median, na.rm = TRUE)

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


    if (!is.null(config$mc_annot)) {
        mc_annot <- parse_mc_annot(config$mc_annot)
    } else {
        if (!is.null(config$cell_type_field)) {
            mc_annot <- adata$obs %>%
                select(cell_type = !!config$cell_type_field) %>%
                rownames_to_column("metacell") %>%
                as_tibble()
        } else {
            # we generate clustering as initial annotation
            feat_mat <- mc_egc[adata$var$top_feature_gene, ]
            km <- cluster_egc(feat_mat, k = config$cluster_k)
            mc_annot <- km$clusters %>% 
                rename(cell_type = cluster) %>% 
                mutate(cell_type = as.character(cell_type))
        }
    }

    if (!is.null(config$cell_type_annot)) {
        cell_type_annot <- parse_cell_type_annot(config$cell_type_annot)
    } else {
        color_of_clusters <- chameleon::data_colors(t(mc_egc[adata$var$top_feature_gene, ]), groups = mc_annot$cell_type)

        cell_type_annot <- enframe(color_of_clusters, name = "cell_type", value = "color") %>%
            mutate(order = 1:n())
    }

    mc_annot <- mc_annot %>% left_join(cell_type_annot %>% select(cell_type, mc_col = color), by = "cell_type")

    cell_type_annot <- cell_type_annot %>%
        arrange(as.numeric(order)) %>%
        mutate(cell_type = factor(cell_type), cell_type = forcats::fct_inorder(cell_type))
    serialize_shiny_data(cell_type_annot, "cell_type_annot", dataset = dataset, cache_dir = cache_dir)

    mc_annot <- mc_annot %>%
        arrange(as.numeric(metacell)) %>%
        left_join(mc_genes_top2, by = "metacell")

    # filter mc_annot to contain only metacells that are within mc_egc
    mc_annot <- mc_annot %>%
        as.data.frame() %>%
        column_to_rownames("metacell")
    mc_annot <- mc_annot[colnames(mc_egc), ]
    mc_annot <- mc_annot %>%
        rownames_to_column("metacell") %>%
        as_tibble()

    mc_annot$top1_lfp <- purrr::map2_dbl(mc_annot$metacell, mc_annot$top1_gene, ~ log2(mc_fp[.y, .x]))
    mc_annot$top2_lfp <- purrr::map2_dbl(mc_annot$metacell, mc_annot$top2_gene, ~ log2(mc_fp[.y, .x]))

    serialize_shiny_data(mc_annot, "mc_annot", dataset = dataset, cache_dir = cache_dir)


    cli_alert_info("Calculating top 30 correlated and anti-correlated genes for each gene")
    gg_mc_top_cor <- calc_gg_mc_top_cor(mc_egc, k = 30)

    serialize_shiny_data(gg_mc_top_cor, "gg_mc_top_cor", dataset = dataset, cache_dir = cache_dir)
}

parse_cell_type_annot <- function(config) {
    cell_type_annot_raw <- fread(config$fn) %>% as_tibble()
    cell_type_annot <- tibble(
        cell_type = as.character(cell_type_annot_raw[[config$cell_type]]),
        color = cell_type_annot_raw[[config$color]]
    )
    if (!is.null(config$order)) {
        cell_type_annot <- cell_type_annot %>% mutate(order = cell_type_annot_raw[[config$order]])
    } else {
        cell_type_annot <- cell_type_annot %>% mutate(order = 1:n())
    }
    return(cell_type_annot)
}


parse_mc_annot <- function(config) {
    mc_annot_raw <- fread(config$fn) %>% as_tibble()
    mc_annot <- tibble(
        metacell = as.character(mc_annot_raw[[config$metacell]]),
        cell_type = as.character(mc_annot_raw[[config$cell_type]])
    )

    if (!is.null(config$cell_type_color)) {
        mc_annot <- mc_annot %>%
            mutate(mc_col = mc_annot_raw[[config$cell_type_color]])
    }


    if (!is.null(config$mc_age)) {
        mc_annot <- mc_annot %>% mutate(mc_age = mc_annot_raw[[config$mc_age]])
    }

    return(mc_annot)
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
