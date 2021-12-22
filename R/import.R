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
#' @param cell_type_field name of a field in the anndata \code{object$obs} which contains a cell type (optional).
#' If the field doesn't exist and \code{metacell_types_file} are missing, MCView would first look
#' for a field named 'type', 'cell_type' or 'cluster' at \code{object$obs}, and if it doesn't exists
#' MCView would cluster the metacell matrix using kmeans++ algorithm (from the \code{tglkmeans} package).
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
#' @param metadata_fields names of fields in the anndata \code{object$obs} which contains metadata for each metacell.
#' The fields should *always* be numeric - if you have cell categorical annotations use
#' \code{cell_metadata_to_metacell} with \code{categorical=TRUE} to convert them to a
#' numeric score (e.g. using fraction of the category). You can use 'all' in order to import all the numeric fields
#'  of the anndata object.
#' @param metadata can be either a data frame with a column named "metacell" with the metacell id and other metadata columns
#' or a name of a delimited file which contains such data frame.
#' @param metadata_colors a named list with colors for each metadata column, or a name of a yaml file with such list.
#' Colors should be given as a list where the first element is a vector of colors and the second element is a vector of breaks.
#' If only colors are given breaks would be implicitly determined from the minimum and maximum of the metadata field.
#' @param calc_gg_cor calculate top 30 correlated and anti-correlated genes for each gene. This computation can be heavy
#' for large datasets or weaker machines, so you can set \code{calc_gg_cor=FALSE} to skip it. Note that then this feature
#' would be missing from the app.
#' @param atlas_project path to and \code{MCView} project which contains the atlas.
#' @param atlas_dataset name of the atlas dataset
#' @param copy_atlas copy atlas MCView to the current project. If FALSE - a symbolic link would be created instaed.
#'
#' @return invisibly returns an \code{AnnDataR6} object of the read \code{anndata_file}
#'
#'
#' @examples
#' \dontrun{
#' dir.create("raw")
#' download.file(
#'     "http://www.wisdom.weizmann.ac.il/~atanay/metac_data/PBMC_processed.tar.gz",
#'     "raw/PBMC_processed.tar.gz"
#' )
#' untar("raw/PBMC_processed.tar.gz", exdir = "raw")
#' create_project("PBMC")
#' import_dataset("PBMC", "PBMC163k", "raw/metacells.h5ad")
#' }
#'
#' @export
import_dataset <- function(project,
                           dataset,
                           anndata_file,
                           cell_type_field = NULL,
                           metacell_types_file = NULL,
                           cell_type_colors_file = NULL,
                           metadata_fields = NULL,
                           metadata = NULL,
                           metadata_colors = NULL,
                           calc_gg_cor = TRUE,
                           atlas_project = NULL,
                           atlas_dataset = NULL,
                           copy_atlas = TRUE) {
    verbose <- !is.null(getOption("MCView.verbose")) && getOption("MCView.verbose")
    verify_project_dir(project, create = TRUE)

    cli_alert_info("Importing {.field {dataset}}")

    cache_dir <- project_cache_dir(project)

    if (!fs::file_exists(anndata_file)) {
        cli_abort("{anndata_file} doesn't exist. Maybe there is a typo?")
    }

    cli_alert_info("Reading {.file {anndata_file}}")
    adata <- anndata::read_h5ad(anndata_file)

    cli_alert_info("Processing metacell matrix")
    mc_mat <- t(adata$X)
    serialize_shiny_data(mc_mat, "mc_mat", dataset = dataset, cache_dir = cache_dir)

    metacells <- colnames(mc_mat)

    mc_sum <- colSums(mc_mat)
    serialize_shiny_data(mc_sum, "mc_sum", dataset = dataset, cache_dir = cache_dir)

    metacells <- rownames(adata$obs)

    metadata <- load_metadata(metadata, metadata_fields, metacells, adata)
    if (!is.null(metadata)) {
        serialize_shiny_data(
            metadata %>% select(metacell, everything()),
            "metadata",
            dataset = dataset,
            cache_dir = cache_dir,
            flat = TRUE
        )
    }

    if (!is.null(metadata_colors)) {
        cli_alert_info("Processing metadata colors")
        if (is.character(metadata_colors)) {
            metadata_colors <- yaml::read_yaml(metadata_colors) %>% as_tibble()
        }
        metadata_colors <- parse_metadata_colors(metadata_colors, metadata)
        serialize_shiny_data(metadata_colors, "metadata_colors", dataset = dataset, cache_dir = cache_dir)
    }

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

    mc_egc <- t(t(mc_mat) / mc_sum)

    cli_alert_info("Calculating top genes per metacell (marker genes)")
    forbidden <- adata$var$forbidden_gene %||% rep(TRUE, nrow(mc_egc))
    marker_genes <- calc_marker_genes(mc_egc[!forbidden, ], 20)
    serialize_shiny_data(marker_genes, "marker_genes", dataset = dataset, cache_dir = cache_dir)

    # serialize the inner fold matrix (if exists)
    if (!is.null(adata$layers[["inner_fold"]])) {
        cli_alert_info("Processing inner-folds matrix")
        inner_fold_mat <- Matrix::t(adata$layers[["inner_fold"]])
        serialize_shiny_data(inner_fold_mat, "inner_fold_mat", dataset = dataset, cache_dir = cache_dir)

        cli_alert_info("Calculating top inner-fold genes")
        marker_genes_inner_fold <- select_top_fold_genes(inner_fold_mat[!forbidden, ])
        serialize_shiny_data(marker_genes_inner_fold, "marker_genes_inner_fold", dataset = dataset, cache_dir = cache_dir)
    }

    mc_genes_top2 <- marker_genes %>%
        group_by(metacell) %>%
        slice(1:2) %>%
        ungroup() %>%
        pivot_wider(names_from = "rank", values_from = c("gene", "fp"))

    colnames(mc_genes_top2) <- c("metacell", "top1_gene", "top2_gene", "top1_lfp", "top2_lfp")

    if (!is.null(metacell_types_file)) {
        if (!is.null(cell_type_field)) {
            cli_alert_warning("{.field cell_type_field} was ignored since {.field metacell_types_file} was set.")
        }
        cli_alert_info("Loading metacell type annotations from {.file {metacell_types_file}}")
        metacell_types <- parse_metacell_types(metacell_types_file, metacells)
    } else {
        if (!is.null(cell_type_field) && !is.null(adata$obs[[cell_type_field]])) {
            cli_alert_info("Taking cell type annotations from {.field {cell_type_field}} field in the anndata object")
            metacell_types <- adata$obs %>%
                select(cell_type = !!cell_type_field) %>%
                rownames_to_column("metacell") %>%
                as_tibble()
        } else if (any(c("type", "cell_type", "cluster") %in% colnames(adata$obs))) {
            cell_type_field <- colnames(adata$obs)[colnames(adata$obs) %in% c("type", "cell_type", "cluster")]
            cell_type_field <- cell_type_field[1]
            cli_alert_info("Taking cell type annotations from {.field {cell_type_field}} field in the anndata object")
            metacell_types <- adata$obs %>%
                select(cell_type = !!cell_type_field) %>%
                rownames_to_column("metacell") %>%
                as_tibble()
        } else {
            cli_alert_info("Clustering in order to get initial annotation.")
            # we generate clustering as initial annotation
            if (rlang::has_name(adata$var, "top_feature_gene")) {
                feat_mat <- mc_egc[adata$var$top_feature_gene, ]
            } else if (rlang::has_name(adata$var, "feature_gene")) {
                feat_mat <- mc_egc[adata$var$feature_gene, ]
            } else {
                cli_abort("{anndata_file} object doesn't have a 'var' field named 'top_feature_gene' or 'feature_gene'")
            }

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
        if (rlang::has_name(adata$var, "top_feature_gene")) {
            feat_mat <- mc_egc[adata$var$top_feature_gene, ]
        } else if (rlang::has_name(adata$var, "feature_gene")) {
            feat_mat <- mc_egc[adata$var$feature_gene, ]
        } else {
            cli_abort("{anndata_file} object doesn't have a 'var' field named 'top_feature_gene' or 'feature_gene'")
        }

        if (all(paste0("umap_", c("x", "y", "u")) %in% colnames(adata$obs))) {
            cli_alert_info("Coloring using pre-calculated 3D umap")
            color_of_clusters <- chameleon::data_colors(adata$obs[, paste0("umap_", c("x", "y", "u"))], groups = metacell_types$cell_type, run_umap = FALSE)
        } else {
            cli_alert_info("Coloring using umap on feature matrix")
            color_of_clusters <- chameleon::data_colors(t(feat_mat), groups = metacell_types$cell_type)
        }

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

    # filter metacell_types to contain only metacells that are included in mc_mat
    metacell_types <- metacell_types %>%
        filter(metacell %in% colnames(mc_mat)) %>%
        as_tibble()

    serialize_shiny_data(metacell_types, "metacell_types", dataset = dataset, cache_dir = cache_dir, flat = TRUE)

    if (calc_gg_cor) {
        if (!is.null(adata$varp$var_similarity)) {
            gg_mc_top_cor <- Matrix::summary(adata$varp$var_similarity) %>%
                rlang::set_names(c("gene1", "gene2", "cor")) %>%
                mutate(
                    gene1 = adata$var_names[gene1],
                    gene2 = adata$var_names[gene2]
                ) %>%
                filter(gene1 != gene2) %>%
                as.data.frame()
            gg_mc_top_cor <- gg_mc_top_cor %>%
                arrange(gene1, desc(cor)) %>%
                group_by(gene1) %>%
                slice(1:30) %>%
                ungroup() %>%
                mutate(type = "pos")
            # TODO: add the anti-correlated when Oren would implement this
        } else {
            cli_alert_info("Calculating top 30 correlated and anti-correlated genes for each gene")
            gg_mc_top_cor <- calc_gg_mc_top_cor(mc_egc, k = 30)
        }
        serialize_shiny_data(gg_mc_top_cor, "gg_mc_top_cor", dataset = dataset, cache_dir = cache_dir)
    } else {
        cli_alert_info("Skipping calculation of top 30 correlated and anti-correlated genes for each gene. Some features in the app would not be available")
    }

    if (!is.null(atlas_project)) {
        if (is.null(atlas_dataset)) {
            cli_abort("Please provide {.code atlas_dataset}")
        }

        import_atlas(adata, atlas_project, atlas_dataset, dataset = dataset, cache_dir = cache_dir, copy_atlas)
    }

    cli_alert_success("{.field {dataset}} dataset imported succesfully to {.path {project}} project")
    invisible(adata)
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
