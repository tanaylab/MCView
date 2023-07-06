#' Compute a umap 2D projection based on gene anchors
#'
#' @noRd
compute_umap <- function(mc_egc, anchors, min_dist = 0.96, n_neighbors = 10, n_epoch = 500, min_log_expr = -14, genes_per_anchor = 30, config = NULL) {
    if (!all(anchors %in% rownames(mc_egc))) {
        cli::cli_abort("Umap gene{?s} {.val {anchors[!(anchors %in% rownames(mc_egc))]}} not found in metacell gene expression data")
    }

    legc <- log2(mc_egc + 1e-5)

    f_gcov <- matrixStats::rowMaxs(legc) >= min_log_expr
    legc <- legc[f_gcov, ]

    if (!all(anchors %in% rownames(legc))) {
        cli::cli_warn("Umap genes: {.val {anchors[!(anchors %in% rownames(legc))]}} do not have a minimum expression of {.val {min_log_expr}} in the metacell gene expression data. Please consider using a lower value for {.val {min_log_expr}}")
    }
    anchors <- anchors[anchors %in% rownames(legc)]

    if (length(anchors) < 2) {
        cli::cli_warn("At least 2 anchors are required to compute a umap projection")
        return(NULL)
    }

    knn <- tgs_cor_knn(t(legc[anchors, ]), t(legc), knn = genes_per_anchor)
    gmods <- knn$col1
    names(gmods) <- knn$col2
    legc_anchors <- legc[names(gmods), ]

    feats <- t(tgs_matrix_tapply(t(legc_anchors), gmods, mean, na.rm = TRUE))

    if (is.null(config)) {
        config <- umap::umap.defaults
        config$min_dist <- min_dist
        config$n_neighbors <- n_neighbors
        config$n_epoch <- n_epoch
    }


    uf <- umap::umap(feats, config = config)

    # graph
    idx_df <- uf$knn$indexes %>%
        as.data.frame() %>%
        rownames_to_column("mc1") %>%
        gather("k", "idx", -mc1) %>%
        as_tibble()

    distances_df <- uf$knn$distances %>%
        as.data.frame() %>%
        rownames_to_column("mc1") %>%
        gather("k", "weight", -mc1) %>%
        as_tibble()

    graph <- idx_df %>%
        left_join(distances_df, by = join_by(mc1, k)) %>%
        mutate(mc2 = rownames(feats)[idx]) %>%
        filter(mc1 != mc2) %>%
        select(mc1, mc2, weight)

    layout_df <- as.data.frame(uf$layout) %>%
        rownames_to_column("mc") %>%
        as_tibble()

    colnames(layout_df) <- c("mc", "umap_x", "umap_y")

    mc2d_list <- list(
        graph = graph,
        mc_id = layout_df$mc,
        mc_x = tibble::deframe(layout_df[, c("mc", "umap_x")]),
        mc_y = tibble::deframe(layout_df[, c("mc", "umap_y")])
    )

    return(mc2d_list)
}

load_default_2d_projection <- function(project, dataset, adata, mc_egc, umap_anchors, min_umap_log_expr, umap_config, genes_per_anchor) {
    cache_dir <- project_cache_dir(project)

    mc2d_list <- NULL
    if (!is.null(umap_anchors)) {
        cli_alert_info("Computing umap 2D projection based on gene anchors")
        mc2d_list <- compute_umap(mc_egc, umap_anchors, min_log_expr = min_umap_log_expr, config = umap_config, genes_per_anchor = genes_per_anchor)
        if (!is.null(mc2d_list)) {
            serialize_shiny_data(umap_anchors, "umap_anchors", dataset = dataset, cache_dir = cache_dir)
        }
    }

    if (is.null(mc2d_list)) {
        if (is.null(adata$obsp$obs_outgoing_weights)) {
            cli_abort_compute_for_mcview("$obsp$obs_outgoing_weights")
        }
        cli_alert_info("Using 2D projection from anndata object")
        graph <- Matrix::summary(adata$obsp$obs_outgoing_weights) %>%
            as.data.frame()

        graph <- graph %>% mutate(i = rownames(adata$obs)[i], j = rownames(adata$obs)[j])

        purrr::walk(c("x", "y"), ~ {
            if (is.null(adata$obs[[.x]])) {
                cli_abort_compute_for_mcview(glue("$obs${.x}"))
            }
        })

        mc2d_list <- list(
            graph = tibble(mc1 = graph[, 1], mc2 = graph[, 2], weight = graph[, 3]),
            mc_id = rownames(adata$obs),
            mc_x = adata$obs %>% select(umap_x = x) %>% tibble::rownames_to_column("mc") %>% tibble::deframe(),
            mc_y = adata$obs %>% select(umap_y = y) %>% tibble::rownames_to_column("mc") %>% tibble::deframe()
        )
    }
    serialize_shiny_data(mc2d_list, "mc2d", dataset = dataset, cache_dir = cache_dir)
}


#' Update default 2D projection for a dataset
#'
#' @param layout a data frame with a column named "metacell" with the metacell id and other columns with the x and y coordinates of the metacell.
#' @param graph a data frame with a column named "from", "to" and "weight" with the ids of the metacells and the weight of the edge. If NULL, the graph would be taken from the anndata object.
#'
#' @inheritParams import_dataset
#'
#' @export
update_2d_projection <- function(project, dataset, layout, graph) {
    cache_dir <- project_cache_dir(project)

    cli_alert_info("Loading 2D projection for dataset {.val {dataset}}")

    metacells <- get_metacell_ids(project, dataset)

    if (is.character(layout)) {
        cli_alert_info("Loading layout from {.val {layout}}")
        layout <- fread(layout)
    }

    purrr::walk(c("metacell", "x", "y"), ~ {
        if (!.x %in% colnames(layout)) {
            cli_abort("Column {.val {.x}} not found in {.field layout}")
        }
    })

    if (is.character(graph)) {
        cli_alert_info("Loading graph from {.val {graph}}")
        graph <- fread(graph)
    }

    purrr::walk(c("from", "to", "weight"), ~ {
        if (!.x %in% colnames(graph)) {
            cli_abort("Column {.val {.x}} not found in {.field graph}")
        }
    })


    mc2d_list <- list(
        graph = graph %>% rename(mc1 = from, mc2 = to),
        mc_id = layout$metacell,
        mc_x = tibble::deframe(layout[, c("metacell", "x")]),
        mc_y = tibble::deframe(layout[, c("metacell", "y")])
    )

    serialize_shiny_data(mc2d_list, "mc2d", dataset = dataset, cache_dir = cache_dir)
}
