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
#' If the field doesn't exist and \code{metacell_types_file} is missing, MCView would first look
#' for a field named 'projected_type', 'type', 'cell_type' or 'cluster' at \code{object$obs} (in this order), and if it doesn't exists MCView would cluster the metacell matrix using kmeans++ algorithm (from the \code{tglkmeans} package).
#' @param metacell_types_file path to a tabular file (csv,tsv) with cell type assignement for
#' each metacell. The file should have a column named "metacell" with the metacell ids and another
#' column named "cell_type", or "cluster" with the cell type assignment. Metacell ids that do
#' not exists in the data would be ignored. \cr
#' If this parameter and \code{cell_type_field} are missing and \code{cluster_metacells=TRUE}, MCView would cluster the
#' metacell matrix using kmeans++ algorithm (from the \code{tglkmeans} package). \cr
#' If the file has a field named 'color' and \code{cell_type_colors_file=NULL}, the cell types colors would
#' be used.
#' @param cell_type_colors_file path to a tabular file (csv,tsv) with color assignement for
#' each cell type. The file should have a column named "cell_type" or "cluster" with the
#' cell types and another column named "color" with the color assignment. \cr
#' In case the metacell types (given by a file or from the \code{anndata_file}) contain types which are not present in the cell type colors file, MCView would use the \code{chameleon} package to assign a color for them in addition to the cell type colors. \cr
#' If this is missing, and \code{metacell_types_file} did not have a 'color' field, MCView would use the \code{chameleon} package to assign a color for each cell type. \cr
#' When an atlas is given (using \code{atlas_project} and \code{atlas_dataset}), if the cell types
#' are the same as the atlas, the atlas colors would be used.
#' @param outliers_anndata_file path to anndata file with outliers (optional). This would enable, by default,
#' the following tabs: ["Outliers", "Similar-fold", "Deviant-fold"]. See the metacells python package for more details.
#' @param cluster_metacells When TRUE and no metacell type is given (via \code{metacell_types_file} or \code{cell_type_field}), MCView would cluster the metacell matrix using kmeans++ algorithm (from the \code{tglkmeans} package).
#' @param cluster_k number of clusters for initial metacell clustering. If NULL - the number of clusters would be determined such that a metacell would contain 16 cells on average.
#' @param metadata_fields names of fields in the anndata \code{object$obs} which contains metadata for each metacell. \cr
#' The fields should can be either numeric or categorical. \cr
#'  You can use \code{cell_metadata_to_metacell} to convert from categorical to a numeric score (e.g. by using fraction of the category). You can use 'all' in order to import all the fields
#'  of the anndata object.
#' @param metadata can be either a data frame with a column named "metacell" with the metacell id and other metadata columns
#' or a name of a delimited file which contains such data frame. See \code{metadata_fields} for other details.
#' @param metadata_colors a named list with colors for each metadata column, or a name of a yaml file with such list.
#' For numerical metadata columns, colors should be given as a list where the first element is a vector of colors and the second element is a vector of breaks. \cr
#' If only colors are given breaks would be implicitly determined from the minimum and maximum of the metadata field. \cr
#' For categorical metadata columns, color can be given either as a named vector where names are the categories and the values are the colors, or as a named list where the first element named 'colors' holds the colors, and the second element
#' called 'categories' holds the categories.
#' @param gene_modules_file path to a tabular file (csv,tsv) with assignment of genes to gene modules. Should have a field named "gene" with the gene name and a field named "module" with the name of the gene module.
#' @param gene_modules_k number of clusters for initial gene module calculation. If NULL - the number of clusters would be determined such that an gene module would contain 16 genes on average.
#' @param calc_gg_cor Calculate top 30 correlated and anti-correlated genes for each gene. This computation can be heavy for large datasets or weaker machines, so you can set \code{calc_gg_cor=FALSE} to skip it. Note that then this feature would be missing from the app.
#' @param gene_names use alternative gene names (optional). A data frame with a column called 'gene_name' with the original gene name (as it appears at the 'h5ad' file) and another column called 'alt_name' with the gene name to use in MCView. Genes that do not appear at the table would not be changed.
#' @param metacell_graphs a named list of metacell graphs or files containing metacell graphs. Each graph should be a data frame columns named "from", "to" and "weight" with the ids of the metacells and the weight of the edge. If the list is not named, the names would be 'graph1', 'graph2' and so on. Note that the graph cannot be named "metacell" as this is reserved for the metacell graph.
#' @param atlas_project path to and \code{MCView} project which contains the atlas.
#' @param atlas_dataset name of the atlas dataset
#' @param projection_weights_file Path to a tabular file (csv,tsv) with the following fields "query", "atlas" and "weight". The file is an output of \code{metacells} projection algorithm.
#' @param copy_atlas copy atlas MCView to the current project. If FALSE - a symbolic link would be created instead.
#' @param minimal_max_log_fraction When choosing marker genes: take only genes with at least one value
#' (in log fraction units - normalized egc) above this threshold
#' @param minimal_relative_log_fraction When choosing marker genes: take only genes with relative
#' log fraction (mc_fp) above this this value
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
#' @inheritDotParams create_project
#' @inheritParams import_cell_metadata
#' @export
import_dataset <- function(project,
                           dataset,
                           anndata_file,
                           cell_type_field = NULL,
                           metacell_types_file = NULL,
                           cell_type_colors_file = NULL,
                           outliers_anndata_file = NULL,
                           cluster_metacells = TRUE,
                           cluster_k = NULL,
                           metadata_fields = NULL,
                           metadata = NULL,
                           metadata_colors = NULL,
                           cell_metadata = NULL,
                           cell_to_metacell = NULL,
                           gene_modules_file = NULL,
                           gene_modules_k = NULL,
                           calc_gg_cor = TRUE,
                           gene_names = NULL,
                           metacell_graphs = NULL,
                           atlas_project = NULL,
                           atlas_dataset = NULL,
                           projection_weights_file = NULL,
                           copy_atlas = TRUE,
                           minimal_max_log_fraction = -13,
                           minimal_relative_log_fraction = 2,
                           ...) {
    verbose <- !is.null(getOption("MCView.verbose")) && getOption("MCView.verbose")
    verify_project_dir(project, create = TRUE, atlas = !is.null(atlas_project), ...)

    cli_alert_info("Importing {.field {dataset}}")

    cache_dir <- project_cache_dir(project)

    if (!fs::file_exists(anndata_file)) {
        cli_abort("{anndata_file} doesn't exist. Maybe there is a typo?")
    }

    cli_alert_info("Reading {.file {anndata_file}}")
    adata <- anndata::read_h5ad(anndata_file)

    if (rlang::has_name(adata$obs, "hidden")) {
        adata <- adata[!adata$obs$hidden, ]
    }

    # sort the gene names lexycographically
    adata <- adata[, sort(colnames(adata$X))]

    cli_alert_info("Processing metacell matrix")
    mc_mat <- t(adata$X)
    rownames(mc_mat) <- modify_gene_names(rownames(mc_mat), gene_names)

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
    if (is.null(adata$obsp$obs_outgoing_weights)) {
        cli_abort_compute_for_mcview("$obsp$obs_outgoing_weights")
    }
    graph <- Matrix::summary(adata$obsp$obs_outgoing_weights) %>%
        as.data.frame()

    graph <- graph %>% mutate(i = rownames(adata$obs)[i], j = rownames(adata$obs)[j])

    purrr::walk(c("umap_x", "umap_y"), ~ {
        if (is.null(adata$obs[[.x]])) {
            cli_abort_compute_for_mcview(glue("$obs${.x}"))
        }
    })

    mc2d_list <- list(
        graph = tibble(mc1 = graph[, 1], mc2 = graph[, 2], weight = graph[, 3]),
        mc_id = rownames(adata$obs),
        mc_x = adata$obs %>% select(umap_x) %>% tibble::rownames_to_column("mc") %>% tibble::deframe(),
        mc_y = adata$obs %>% select(umap_y) %>% tibble::rownames_to_column("mc") %>% tibble::deframe()
    )
    serialize_shiny_data(mc2d_list, "mc2d", dataset = dataset, cache_dir = cache_dir)

    if (!is.null(metacell_graphs)) {
        cli_alert_info("Processing metacell graphs")
        metacell_graphs <- read_metacell_graphs(metacell_graphs, metacells)
        serialize_shiny_data(metacell_graphs, "metacell_graphs", dataset = dataset, cache_dir = cache_dir)
    }

    mc_egc <- t(t(mc_mat) / mc_sum)

    cli_alert_info("Calculating top genes per metacell (marker genes)")
    forbidden_field <- "lateral_gene"
    if (!("lateral_gene" %in% colnames(adata$var)) && "forbidden_gene" %in% colnames(adata$var)) {
        forbidden_field <- "forbidden_gene"
    }
    if (is.null(adata$var[, forbidden_field])) {
        forbidden_genes <- c()
        forbidden <- rep(FALSE, nrow(mc_egc))
    } else {
        forbidden <- adata$var[, forbidden_field]
        forbidden_genes <- rownames(adata$var)[adata$var[, forbidden_field]]
        forbidden_genes <- modify_gene_names(forbidden_genes, gene_names)
    }

    serialize_shiny_data(forbidden_genes, "forbidden_genes", dataset = dataset, cache_dir = cache_dir)
    marker_genes <- calc_marker_genes(mc_egc[!forbidden, ], 20, minimal_max_log_fraction = minimal_max_log_fraction, minimal_relative_log_fraction = minimal_relative_log_fraction)
    serialize_shiny_data(marker_genes, "marker_genes", dataset = dataset, cache_dir = cache_dir)

    # serialize the inner fold matrix (if exists)
    if (!is.null(adata$layers[["inner_fold"]])) {
        cli_alert_info("Processing inner-folds matrix")
        inner_fold_mat <- t(adata$layers[["inner_fold"]])
        rownames(inner_fold_mat) <- rownames(mc_mat)
        colnames(inner_fold_mat) <- colnames(mc_mat)
        serialize_shiny_data(inner_fold_mat, "inner_fold_mat", dataset = dataset, cache_dir = cache_dir)

        cli_alert_info("Calculating top inner-fold genes")
        inner_fold_genes <- rownames(inner_fold_mat)[Matrix::rowSums(inner_fold_mat) > 0]
        inner_fold_gene_metacells <- matrixStats::rowSums2(as.matrix(inner_fold_mat[inner_fold_genes, , drop = FALSE]) > 0)
        # fp here is the number of non-zero entries per gene
        marker_genes_inner_fold <- tibble(gene = inner_fold_genes, fp = inner_fold_gene_metacells) %>%
            arrange(desc(fp))
        serialize_shiny_data(marker_genes_inner_fold, "marker_genes_inner_fold", dataset = dataset, cache_dir = cache_dir)
        cli_alert_info("Add the {.field \"Inner-fold\"} tab to your config file to view the inner-fold matrix")
    }

    mc_genes_top2 <- marker_genes %>%
        group_by(metacell) %>%
        slice(1:2) %>%
        ungroup() %>%
        pivot_wider(names_from = "rank", values_from = c("gene", "fp"))

    colnames(mc_genes_top2) <- c("metacell", "top1_gene", "top2_gene", "top1_lfp", "top2_lfp")

    cell_type_colors <- NULL
    if (!is.null(metacell_types_file)) {
        if (!is.null(cell_type_field)) {
            cli_alert_warning("{.field cell_type_field} was ignored since {.field metacell_types_file} was set.")
        }
        cli_alert_info("Loading metacell type annotations from {.file {metacell_types_file}}")
        metacell_types <- fread(metacell_types_file)
        if (has_name(metacell_types, "color")) {
            if (is.null(cell_type_colors_file)) {
                cli_alert_info("Loading cell type colors from the 'color' field at {.file {metacell_types_file}}")
                cell_type_colors <- parse_cell_type_colors(metacell_types)
            }
        }
        metacell_types <- parse_metacell_types(metacell_types, metacells)
    } else {
        if (!is.null(cell_type_field) && !is.null(adata$obs[[cell_type_field]])) {
            cli_alert_info("Taking cell type annotations from {.field {cell_type_field}} field in the anndata object")
            metacell_types <- adata$obs %>%
                select(cell_type = !!cell_type_field) %>%
                rownames_to_column("metacell") %>%
                as_tibble()
        } else if (any(c("projected_type", "type", "cell_type", "cluster") %in% colnames(adata$obs))) {
            cell_type_field <- colnames(adata$obs)[colnames(adata$obs) %in% c("projected_type", "type", "cell_type", "cluster")]
            cell_type_field <- cell_type_field[1]
            cli_alert_info("Taking cell type annotations from {.field {cell_type_field}} field in the anndata object")
            metacell_types <- adata$obs %>%
                select(cell_type = !!cell_type_field) %>%
                rownames_to_column("metacell") %>%
                as_tibble()
        } else if (cluster_metacells) {
            cli_alert_info("Clustering in order to get initial annotation.")
            # we generate clustering as initial annotation
            if (rlang::has_name(adata$var, "top_feature_gene")) {
                feat_mat <- mc_egc[adata$var$top_feature_gene, ]
            } else if (rlang::has_name(adata$var, "feature_gene")) {
                feat_mat <- mc_egc[adata$var$feature_gene, ]
            } else {
                cli_abort("{anndata_file} object doesn't have a 'var' field named 'top_feature_gene' or 'feature_gene'")
            }

            km <- cluster_egc(feat_mat, verbose = verbose, k = cluster_k)
            metacell_types <- km$clusters %>%
                rename(cell_type = cluster) %>%
                mutate(cell_type = as.character(cell_type))
        } else {
            metacell_types <- tibble(metacell = metacells, cell_type = NA)
        }
    }

    metacell_types <- metacell_types %>% mutate(cell_type = forcats::fct_explicit_na(factor(cell_type)))

    if (!is.null(cell_type_colors_file)) {
        cli_alert_info("Loading cell type color annotations from {.file {cell_type_colors_file}}")
        cell_type_colors <- parse_cell_type_colors(cell_type_colors_file)

        # Add cell types that are missing
        missing_cell_types <- setdiff(unique(metacell_types$cell_type), cell_type_colors$cell_type)
        missing_cell_types <- missing_cell_types[missing_cell_types != "(Missing)"]

        if (length(missing_cell_types) > 0) {
            cli_alert_warning("The following cell types are missing from the color annotations: {.field {missing_cell_types}}. Adding them to the color annotations with random colors.")

            new_cell_type_colors <- color_cell_types(adata, mc_egc, metacell_types) %>%
                filter(cell_type %in% missing_cell_types)

            cell_type_colors <- bind_rows(
                cell_type_colors,
                new_cell_type_colors
            ) %>%
                distinct(cell_type, color) %>%
                mutate(order = 1:n())
        }
    } else if (is.null(cell_type_colors)) {
        if (!is.null(atlas_dataset) && !is.null(atlas_project)) { # use atlas colors
            atlas_colors <- fread(fs::path(project_cache_dir(atlas_project), atlas_dataset, "cell_type_colors.tsv")) %>% as_tibble()
            atlas_colors <- atlas_colors %>%
                bind_rows(
                    tibble(cell_type = c("Dissimilar", "Mixture", "Doublet"), color = c("gray", "darkgray", "black"), order = max(atlas_colors$order) + 1:3)
                ) %>%
                distinct(cell_type, color, .keep_all = TRUE)
            # if we are using the atlas cell types - use their colors
            if (all(metacell_types$cell_type %in% atlas_colors$cell_type)) {
                cli_alert_info("Loading cell type color annotations from {.field atlas}")
                cell_type_colors <- atlas_colors
            } else {
                types_without_colors <- setdiff(unique(metacell_types$cell_type), atlas_colors$cell_type)
                types_without_colors <- paste(types_without_colors, collapse = ", ")
                cli_abort("The following cell types do not have colors at the atlas: {.field {types_without_colors}}. To fix it either provide {.code cell_type_colors_file} or add the cell type(s) to the atlas colors.")
            }
        } else {
            cell_type_colors <- color_cell_types(adata, mc_egc, metacell_types)
        }
    }

    cell_type_colors <- cell_type_colors %>%
        arrange(as.numeric(order)) %>%
        mutate(cell_type = factor(cell_type), cell_type = forcats::fct_inorder(cell_type))
    serialize_shiny_data(cell_type_colors, "cell_type_colors", dataset = dataset, cache_dir = cache_dir, flat = TRUE)

    metacell_types <- metacell_types %>%
        left_join(mc_genes_top2, by = "metacell")

    # filter metacell_types to contain only metacells that are included in mc_mat
    metacell_types <- metacell_types %>%
        filter(metacell %in% colnames(mc_mat)) %>%
        as_tibble()

    serialize_shiny_data(metacell_types, "metacell_types", dataset = dataset, cache_dir = cache_dir, flat = TRUE)

    if (!is.null(gene_modules_file)) {
        gene_modules <- parse_gene_modules_file(gene_modules_file)
    } else {
        gene_modules <- calc_gene_modules(mc_mat[!forbidden, ], k = gene_modules_k)
    }
    serialize_shiny_data(gene_modules, "gene_modules", dataset = dataset, cache_dir = cache_dir, flat = TRUE)


    if (!is.null(adata$varp$var_similarity)) {
        if (!methods::is(adata$varp$var_similarity, "sparseMatrix")) {
            cli_abort("{.field var_similarity} matrix is not a sparse matrix. This probably means that you are running an old version of the {.field metacells} python moudle. Please update the module, rerun {.field compute_for_mcview} and try again.")
        }
        cli_alert_info("Loading previously calculated 30 correlated and anti-correlated genes for each gene")

        gg_mc_top_cor <- Matrix::summary(adata$varp$var_similarity) %>%
            rlang::set_names(c("gene1", "gene2", "cor")) %>%
            mutate(
                gene1 = adata$var_names[gene1],
                gene2 = adata$var_names[gene2]
            ) %>%
            filter(gene1 != gene2) %>%
            as.data.frame()

        gg_mc_top_cor_pos <- gg_mc_top_cor %>%
            arrange(gene1, desc(cor)) %>%
            group_by(gene1) %>%
            slice(1:15) %>%
            ungroup() %>%
            mutate(type = "pos")
        gg_mc_top_cor_neg <- gg_mc_top_cor %>%
            arrange(gene1, cor) %>%
            group_by(gene1) %>%
            slice(1:15) %>%
            ungroup() %>%
            mutate(type = "neg")
        gg_mc_top_cor <- bind_rows(
            gg_mc_top_cor_pos,
            gg_mc_top_cor_neg
        ) %>%
            arrange(gene1, desc(cor))
        gg_mc_top_cor <- gg_mc_top_cor %>%
            mutate(gene1 = modify_gene_names(gene1, gene_names), gene2 = modify_gene_names(gene2, gene_names))
        serialize_shiny_data(gg_mc_top_cor, "gg_mc_top_cor", dataset = dataset, cache_dir = cache_dir)
    } else {
        if (calc_gg_cor) {
            cli_alert_info("Calculating top 30 correlated and anti-correlated genes for each gene")
            gg_mc_top_cor <- calc_gg_mc_top_cor(mc_egc, k = 30)
            serialize_shiny_data(gg_mc_top_cor, "gg_mc_top_cor", dataset = dataset, cache_dir = cache_dir)
        } else {
            cli_alert_info("Skipping calculation of top 30 correlated and anti-correlated genes for each gene. Some features in the app would not be available")
        }
    }


    if (!is.null(cell_metadata)) {
        if (is.null(cell_to_metacell)) {
            cli_abort("Please provide also {.field cell_to_metacell} in order to import {.field cell_metadata}")
        }
        import_cell_metadata(project, dataset, cell_metadata, cell_to_metacell)
    }

    if (!is.null(outliers_anndata_file)) {
        load_outliers(outliers_anndata_file, project, dataset, gene_names = gene_names)
    }


    if (!is.null(atlas_project)) {
        if (is.null(atlas_dataset)) {
            cli_abort("Please provide {.code atlas_dataset}")
        }

        if (is.null(projection_weights_file)) {
            cli_abort("Please provide {.code projection_weights_file}")
        }

        if (!is.null(gene_names)) {
            cli_warn("Use {.field gene_names} with atlas only if you are sure that the atlas and the dataset use the same gene names. Otherwise, the atlas will not be imported correctly.")
        }

        import_atlas(adata, atlas_project, atlas_dataset, projection_weights_file, dataset = dataset, cache_dir = cache_dir, copy_atlas = copy_atlas, gene_names = gene_names)
    }

    cli_alert_success("{.field {dataset}} dataset imported succesfully to {.path {project}} project")
    cli::cli_ul("You can now run the app using: {.field run_app(\"{project}\")}")
    cli::cli_ul("or create a bundle using: {.field create_bundle(\"{project}\", name = \"name_of_bundle\")}")
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

color_cell_types <- function(adata, mc_egc, metacell_types) {
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

    return(cell_type_colors)
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

cli_abort_compute_for_mcview <- function(field) {
    cli_abort("{.field {field}} is missing from the h5ad file. Did you remember to run {.code compute_for_mcview} using the metacells python package?")
}
