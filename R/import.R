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
#' @param dataset name for the dataset, e.g. "PBMC". The name of the dataset can only contain alphanumeric characters, dots, dashes and underscores.
#' @param anndata_file path to \code{h5ad} file which contains the output of metacell2 pipeline (metacells python package) or a loaded anndata object of the same format.
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
#' @param cluster_metacells When TRUE and no metacell type is given (via \code{metacell_types_file} or \code{cell_type_field} - implicit and explicit), MCView would cluster the metacell matrix using kmeans++ algorithm (from the \code{tglkmeans} package).
#' @param cluster_k number of clusters for initial metacell clustering. If NULL - the number of clusters would be determined such that a metacell would contain 16 cells on average.
#' @param metadata_fields names of fields in the anndata \code{object$obs} which contains metadata for each metacell. \cr
#' The fields should can be either numeric or categorical. \cr
#'  You can use \code{cell_metadata_to_metacell} to convert from categorical to a numeric score (e.g. by using fraction of the category). You can use 'all' in order to import all the fields
#'  of the anndata object.
#' @param metadata can be either a data frame with a column named "metacell" with the metacell id and other metadata columns
#' or a name of a delimited file which contains such data frame. See \code{metadata_fields} for other details.
#' @param metadata_colors a named list with colors for each metadata column, or a name of a yaml file with such list.
#' For numerical metadata columns, colors should be given as a list where the first element is a vector of colors and the second element is a vector of breaks. \cr
#' Note that in case the breaks are out of the range of the metadata values, the colors would be used for the minimum and maximum values. \cr
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
#' @param umap_anchors a vector of gene names to use for UMAP calculation. If NULL, the umap from the anndata object would be used.
#' @param umap_config a named list with UMAP configuration. See \code{umap::umap} for more details. When NULL, the default configuration would be used, except for: min_dist=0.96, n_neighbors=10, n_epoch=500.
#' @param min_umap_log_expr minimal log2 expression for genes to use for UMAP calculation.
#' @param genes_per_anchor number of genes to use for each umap anchor.
#' @param layout a data frame with a column named "metacell" with the metacell id and other columns with the x and y coordinates of the metacell. If NULL, the layout would be taken from the anndata object.
#' @param default_graph a data frame with a column named "from", "to" and "weight" with the ids of the metacells and the weight of the edge. If NULL, the graph would be taken from the anndata object.
#' @param overwrite if a dataset with the same name already exists, overwrite it. Otherwise, an error would be thrown.
#' @param copy_source_file if TRUE, copy the source file to the project cache directory. If FALSE, create a symbolic link to the source file.
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
                           umap_anchors = NULL,
                           umap_config = NULL,
                           min_umap_log_expr = -14,
                           genes_per_anchor = 30,
                           layout = NULL,
                           default_graph = NULL,
                           overwrite = TRUE,
                           copy_source_file = FALSE,
                           ...) {
    if (missing(project)) {
        cli::cli_abort("Please provide a {.field project} path")
    }

    if (missing(anndata_file)) {
        cli::cli_abort("Please provide a path to an {.field anndata_file} (output of the metacells python package)")
    }

    verbose <- !is.null(getOption("MCView.verbose")) && getOption("MCView.verbose")
    project <- init_project_dir(project, create = TRUE, overwrite = overwrite, atlas = !is.null(atlas_project), ...)

    if (missing(dataset)) {
        dataset <- basename(project)
        cli::cli_alert_info("Setting dataset name to {.field {dataset}}. Set {.field dataset} to a different name if you want to change it.")
    }

    if (!grepl("^[A-Za-z0-9_.-]+$", dataset)) {
        cli_abort("Dataset name can only contain letters, numbers, '.', '-' and '_'")
    }

    cli_alert_info("Importing {.field {dataset}}")

    # if dataset directory already exists and overwrite is FALSE - throw an error
    if (fs::dir_exists(fs::path(project_cache_dir(project), dataset))) {
        if (!overwrite) {
            cli_abort("Dataset {.field {dataset}} already exists in project {.field {project}}. If you want to overwrite it, set {.field overwrite=TRUE}.")
        } else {
            cli_alert_warning("Dataset {.field {dataset}} already exists in project {.field {project}}. Overwriting it.")
        }
        # delete the dataset directory
        fs::dir_delete(fs::path(project_cache_dir(project), dataset))
    }

    cache_dir <- project_cache_dir(project)

    # if anndata is an anndata object - just use it
    if (inherits(anndata_file, "AnnDataR6")) {
        cli_alert_info("Using the provided anndata object")
        adata <- anndata_file
    } else {
        if (!fs::file_exists(anndata_file)) {
            cli_abort("{.file {anndata_file}} doesn't exist. Maybe there is a typo?")
        }

        cli_alert_info("Reading {.file {anndata_file}}")
        adata <- anndata::read_h5ad(anndata_file)
    }

    if (is.null(adata$uns$mcview_format)) {
        cli_abort("The anndata file {.file {anndata_file}} is missing the {.field mcview_format} field. Did you run the {.code compute_for_mcview()} function in the {.pkg metacells} package?")
    }

    if (adata$uns$mcview_format != "1.0") {
        cli_abort("The file {.file {anndata_file}} was created by an old version of {.pkg metacells}. Please convert it to the new format and try again. The conversion script can be found at: {.url https://github.com/tanaylab/metacells/blob/master/bin/convert_0.8_to_0.9.py}")
    }

    if (rlang::has_name(adata$obs, "hidden")) {
        adata <- adata[!adata$obs$hidden, ]
    }

    # sort the gene names lexycographically
    adata <- adata[, sort(colnames(adata$X))]

    cli_alert_info("Processing metacell matrix")
    if (is.null(adata$obs$total_umis)) {
        cli_abort("The file {.file {anndata_file}} doesn't contain the metacell sum. Please convert it to the new format and try again.")
    }
    mc_sum <- adata$obs$total_umis
    names(mc_sum) <- rownames(adata$obs)
    serialize_shiny_data(mc_sum, "mc_sum", dataset = dataset, cache_dir = cache_dir)

    mc_mat <- t(adata$X * mc_sum)
    mc_mat <- as.matrix(mc_mat)
    rownames(mc_mat) <- modify_gene_names(rownames(mc_mat), gene_names)

    serialize_shiny_data(mc_mat, "mc_mat", dataset = dataset, cache_dir = cache_dir)

    metacells <- colnames(mc_mat)

    mc_egc <- t(t(mc_mat) / mc_sum)

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
    if (!is.null(layout) && !is.null(default_graph)) {
        update_2d_projection(project, dataset, layout, default_graph)
    } else {
        load_default_2d_projection(project, dataset, adata, mc_egc, umap_anchors, min_umap_log_expr, umap_config, genes_per_anchor)
    }

    if (!is.null(metacell_graphs)) {
        cli_alert_info("Processing metacell graphs")
        metacell_graphs <- read_metacell_graphs(metacell_graphs, metacells)
        serialize_shiny_data(metacell_graphs, "metacell_graphs", dataset = dataset, cache_dir = cache_dir)
    }

    cli_alert_info("Calculating top genes per metacell (marker genes)")
    lateral_field <- "lateral_gene"

    if (is.null(adata$var[, lateral_field])) {
        lateral_genes <- c()
        lateral <- rep(FALSE, nrow(mc_egc))
    } else {
        lateral <- adata$var[, lateral_field]
        lateral_genes <- rownames(adata$var)[adata$var[, lateral_field]]
        lateral_genes <- modify_gene_names(lateral_genes, gene_names)
    }

    serialize_shiny_data(lateral_genes, "lateral_genes", dataset = dataset, cache_dir = cache_dir)

    if (has_name(adata$var, "noisy_gene")) {
        noisy <- adata$var[, "noisy_gene"]
        noisy_genes <- rownames(adata$var)[adata$var[, "noisy_gene"]]
        noisy_genes <- modify_gene_names(noisy_genes, gene_names)
    } else {
        noisy <- rep(FALSE, nrow(mc_egc))
        noisy_genes <- c()
    }

    serialize_shiny_data(noisy_genes, "noisy_genes", dataset = dataset, cache_dir = cache_dir)


    if (is.null(adata$var[, "marker_gene"])) {
        marker_genes <- calc_marker_genes(mc_egc[!lateral, ], 20, minimal_max_log_fraction = minimal_max_log_fraction, minimal_relative_log_fraction = minimal_relative_log_fraction)
    } else {
        marker_genes <- calc_marker_genes(mc_egc[adata$var[, "marker_gene"], ], 20, minimal_max_log_fraction = minimal_max_log_fraction, minimal_relative_log_fraction = minimal_relative_log_fraction)
    }
    serialize_shiny_data(marker_genes, "marker_genes", dataset = dataset, cache_dir = cache_dir)

    # cache metacell correlations of default marker genes
    cli_alert_info("Calculating metacell correlations of default marker genes")
    default_markers <- choose_markers(marker_genes, 100)
    m_norm <- mc_egc[default_markers, ] + 1e-5
    mc_fp_markers <- m_norm / apply(m_norm, 1, median, na.rm = TRUE)
    zero_mcs <- colSums(abs(mc_fp_markers) > 0) < 2
    mc_fp_markers <- mc_fp_markers[, !zero_mcs, drop = FALSE]
    default_markers_dist <- tgs_dist(tgs_cor(mc_fp_markers, pairwise.complete.obs = TRUE))
    serialize_shiny_data(default_markers_dist, "default_markers_dist", dataset = dataset, cache_dir = cache_dir)
    serialize_shiny_data(default_markers, "default_markers", dataset = dataset, cache_dir = cache_dir)

    # serialize the inner fold matrix (if exists)
    if (!is.null(adata$layers[["inner_fold"]])) {
        cli_alert_info("Processing inner-folds matrix")
        inner_fold_mat <- t(adata$layers[["inner_fold"]])
        rownames(inner_fold_mat) <- rownames(mc_mat)
        colnames(inner_fold_mat) <- colnames(mc_mat)
        inner_fold_mat <- inner_fold_mat[!noisy, , drop = FALSE]

        serialize_shiny_data(inner_fold_mat, "inner_fold_mat", dataset = dataset, cache_dir = cache_dir)

        cli_alert_info("Calculating top inner-fold genes")
        inner_fold_genes <- rownames(inner_fold_mat)[Matrix::rowSums(inner_fold_mat) > 0]
        inner_fold_gene_metacells <- matrixStats::rowSums2(as.matrix(inner_fold_mat[inner_fold_genes, , drop = FALSE]) > 0)
        # fp here is the number of non-zero entries per gene
        marker_genes_inner_fold <- tibble(gene = inner_fold_genes, fp = inner_fold_gene_metacells) %>%
            arrange(desc(fp))
        serialize_shiny_data(marker_genes_inner_fold, "marker_genes_inner_fold", dataset = dataset, cache_dir = cache_dir)
        add_tab("Inner-fold", project)
    }

    if (!is.null(adata$layers[["inner_stdev_log"]])) {
        cli_alert_info("Processing inner-stdev matrix")
        inner_stdev_mat <- t(adata$layers[["inner_stdev_log"]])
        rownames(inner_stdev_mat) <- rownames(mc_mat)
        colnames(inner_stdev_mat) <- colnames(mc_mat)
        inner_stdev_mat <- inner_stdev_mat[!noisy, , drop = FALSE]

        serialize_shiny_data(inner_stdev_mat, "inner_stdev_mat", dataset = dataset, cache_dir = cache_dir)

        cli_alert_info("Calculating top inner-stdev genes")
        inner_stdev_genes <- rownames(inner_stdev_mat)[Matrix::rowSums(inner_stdev_mat) > 0]
        inner_stdev_gene_metacells <- matrixStats::rowSums2(as.matrix(inner_stdev_mat[inner_stdev_genes, , drop = FALSE]) > 0)
        # fp here is the number of non-zero entries per gene
        marker_genes_inner_stdev <- tibble(gene = inner_stdev_genes, fp = inner_stdev_gene_metacells) %>%
            arrange(desc(fp))
        serialize_shiny_data(marker_genes_inner_stdev, "marker_genes_inner_stdev", dataset = dataset, cache_dir = cache_dir)
        add_tab("Stdev-fold", project)
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
        if (!is.null(cell_type_field)) {
            if (!has_name(adata$obs, cell_type_field)) {
                cli_abort("{.field {cell_type_field}} is not a field in the anndata object")
            }

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
            feat_mat <- mc_egc[choose_markers(marker_genes, 1e3), ]

            km <- cluster_egc(feat_mat, verbose = verbose, k = cluster_k)
            metacell_types <- km$clusters %>%
                rename(cell_type = cluster) %>%
                mutate(cell_type = as.character(cell_type))
        } else {
            metacell_types <- tibble(metacell = metacells, cell_type = NA)
        }
    }

    metacell_types <- metacell_types %>% mutate(cell_type = forcats::fct_na_value_to_level(factor(cell_type), "(Missing)"))

    if (!is.null(cell_type_colors_file)) {
        cli_alert_info("Loading cell type color annotations from {.file {cell_type_colors_file}}")
        cell_type_colors <- parse_cell_type_colors(cell_type_colors_file)

        # Add cell types that are missing
        missing_cell_types <- setdiff(unique(metacell_types$cell_type), cell_type_colors$cell_type)
        missing_cell_types <- missing_cell_types[missing_cell_types != "(Missing)"]

        if (length(missing_cell_types) > 0) {
            cli_alert_warning("The following cell types are missing from the color annotations: {.field {missing_cell_types}}. Adding them to the color annotations with random colors.")
            new_cell_type_colors <- color_cell_types(adata, mc_egc, metacell_types, choose_markers(marker_genes, 1e3)) %>%
                filter(cell_type %in% missing_cell_types) %>%
                mutate(cell_type = as.character(cell_type))

            cell_type_colors <- bind_rows(
                cell_type_colors %>% mutate(cell_type = as.character(cell_type)),
                new_cell_type_colors
            ) %>%
                distinct(cell_type, color) %>%
                mutate(order = 1:n())
        }
    } else if (is.null(cell_type_colors)) {
        if (!is.null(atlas_dataset) && !is.null(atlas_project)) { # use atlas colors
            atlas_colors <- fread(fs::path(project_cache_dir(atlas_project), atlas_dataset, "cell_type_colors.tsv"), colClasses = c("cell_type" = "character", "color" = "character")) %>% as_tibble()
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
            cell_type_colors <- color_cell_types(adata, mc_egc, metacell_types, choose_markers(marker_genes, 1e3))
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
        gene_modules <- calc_gene_modules(mc_mat[!lateral, ], k = gene_modules_k)
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

    # cell metadata
    if (!is.null(cell_metadata)) {
        import_cell_metadata(project, dataset, cell_metadata, cell_to_metacell)
    }

    if (!is.null(outliers_anndata_file)) {
        load_outliers(outliers_anndata_file, project, dataset, gene_names = gene_names)
    }

    mc_qc_metadata <- adata$obs %>%
        rownames_to_column("metacell") %>%
        select(metacell, umis = total_umis, cells = grouped)

    if (!is.null(adata$layers[["inner_fold"]])) {
        max_inner_fold <- apply(inner_fold_mat[!(rownames(inner_fold_mat) %in% noisy_genes), , drop = FALSE], 2, max, na.rm = TRUE) %>%
            tibble::enframe(name = "metacell", value = "max_inner_fold")
        max_inner_fold_no_lateral <- apply(inner_fold_mat[!(rownames(inner_fold_mat) %in% noisy_genes) & !(rownames(inner_fold_mat) %in% lateral_genes), , drop = FALSE], 2, max, na.rm = TRUE) %>%
            tibble::enframe(name = "metacell", value = "max_inner_fold_no_lateral")
        mc_qc_metadata <- mc_qc_metadata %>%
            left_join(max_inner_fold, by = "metacell") %>%
            left_join(max_inner_fold_no_lateral, by = "metacell")
    }

    if (!is.null(adata$layers[["inner_stdev_log"]])) {
        max_inner_stdev_log <- apply(inner_stdev_mat, 2, max, na.rm = TRUE) %>%
            tibble::enframe(name = "metacell", value = "max_inner_stdev_log")
        mc_qc_metadata <- mc_qc_metadata %>%
            left_join(max_inner_stdev_log, by = "metacell")
    }

    if (!is.null(adata$layers[["zeros"]]) && has_name(adata$obs, "__zeros_downsample_umis")) {
        obs_zeros <- adata$layers[["zeros"]]

        # expected number of zeros assuming a poisson distribution: e^(-T * lambda)*N where T is the number of UMIs (downsampled), lambda is the average number of UMIs per cell and N is the number of cells
        exp_zeros <- exp(-as.matrix(adata$X) * adata$obs$`__zeros_downsample_umis`) * adata$obs$grouped
        obs_zeros <- obs_zeros[, !noisy, drop = FALSE]
        exp_zeros <- exp_zeros[, !noisy, drop = FALSE]
        zero_fold <- log(obs_zeros + 1) - log(exp_zeros + 1)
        gene_max_folds <- matrixStats::colMaxs(zero_fold)
        names(gene_max_folds) <- colnames(zero_fold)
        gene_avgs <- log2(rowMeans(mc_egc) + 1e-5)
        idxs <- apply(zero_fold, 2, which.max)

        gene_zero_fold_df <- tibble::enframe(gene_max_folds, name = "gene", value = "zero_fold") %>%
            mutate(
                gene = modify_gene_names(gene, gene_names),
                avg = gene_avgs[gene],
                obs = purrr::map2_dbl(gene, idxs, ~ obs_zeros[.y, .x]),
                exp = purrr::map2_dbl(gene, idxs, ~ exp_zeros[.y, .x]),
                metacell = rownames(obs_zeros)[idxs],
                type = case_when(
                    (gene %in% lateral_genes) & (gene %in% noisy_genes) ~ "lateral, noisy",
                    gene %in% noisy_genes ~ "noisy",
                    gene %in% lateral_genes ~ "lateral",
                    TRUE ~ "other"
                )
            ) %>%
            arrange(desc(zero_fold))

        serialize_shiny_data(gene_zero_fold_df, "gene_zero_fold", dataset = dataset, cache_dir = cache_dir)

        mc_max_folds <- matrixStats::rowMaxs(zero_fold)
        names(mc_max_folds) <- rownames(zero_fold)

        mc_qc_metadata <- mc_qc_metadata %>%
            mutate(zero_fold = mc_max_folds[metacell])
    }

    if (has_name(adata$obs, "projected_correlation")) {
        mc_qc_metadata <- mc_qc_metadata %>%
            mutate(projected_correlation = adata$obs$projected_correlation)
    }

    serialize_shiny_data(mc_qc_metadata, "mc_qc_metadata", dataset = dataset, cache_dir = cache_dir)

    qc_stats <- list(
        n_outliers = adata$uns$outliers,
        n_cells = sum(mc_qc_metadata$cells),
        n_umis = sum(mc_qc_metadata$umis),
        median_umis_per_metacell = median(mc_qc_metadata$umis, na.rm = TRUE),
        median_cells_per_metacell = median(mc_qc_metadata$cells, na.rm = TRUE)
    )
    serialize_shiny_data(qc_stats, "qc_stats", dataset = dataset, cache_dir = cache_dir)

    if (has_name(adata$var, "gene")) {
        cli::cli_abort("A column named {.field 'gene'} already exists in the var slot of the anndata object. Please rename it to avoid conflicts.")
    }

    gene_qc <- adata$var %>%
        rownames_to_column("gene") %>%
        select(gene) %>%
        mutate(max_expr = matrixStats::rowMaxs(mc_egc)) %>%
        mutate(type = case_when(
            (gene %in% lateral_genes) & (gene %in% noisy_genes) ~ "lateral, noisy",
            gene %in% noisy_genes ~ "noisy",
            gene %in% lateral_genes ~ "lateral",
            TRUE ~ "other"
        ))

    # genes QC
    if (has_name(adata$var, "significant_inner_folds_count")) {
        gene_qc <- gene_qc %>%
            mutate(significant_inner_folds_count = adata$var$significant_inner_folds_count) %>%
            arrange(desc(significant_inner_folds_count)) %>%
            as_tibble()
    }

    if (has_name(adata$var, "correction_factor")) {
        if (has_name(adata$var, "gene")) {
            cli::cli_abort("A column named {.field gene} already exists in the var slot of the anndata object. Please rename it to avoid conflicts.")
        }
        gene_qc <- gene_qc %>%
            left_join(adata$var %>% rownames_to_column("gene") %>% select(gene, correction_factor), by = "gene")
    }

    if (any(grepl("^fitted_gene", colnames(adata$var)))) {
        if (has_name(adata$var, "gene")) {
            cli::cli_abort("A column named {.field gene} already exists in the var slot of the anndata object. Please rename it to avoid conflicts.")
        }
        gene_qc <- gene_qc %>%
            left_join(adata$var %>% rownames_to_column("gene") %>% select(gene, starts_with("fitted_gene")), by = "gene")
    }

    serialize_shiny_data(gene_qc, "gene_qc", dataset = dataset, cache_dir = cache_dir)

    if (!is.null(atlas_dataset) && is.null(atlas_project)) {
        cli_abort("Please provide {.code atlas_project} if you provide {.code atlas_dataset}")
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

    # write the version of the package
    writeLines(
        as.character(utils::packageVersion("MCView")),
        project_version_file(project)
    )
    cli::cli_alert_info("MCView version: {.field {utils::packageVersion('MCView')}}")

    if (!is.null(adata$uns$metacells_algorithm)) {
        writeLines(
            as.character(adata$uns$metacells_algorithm),
            project_metacells_algorithm_file(project)
        )
        cli::cli_alert_info("Metacells algorithm version: {.field {adata$uns$metacells_algorithm}}")
    }

    # Add app.R file to make the project ready for direct deployment
    add_app_file(project)
    
    if (inherits(anndata_file, "AnnDataR6")){
        if (copy_source_file){
            anndata::write_h5ad(adata, source_metacells_file_path(project))
            cli::cli_alert_info("Copying the source anndata file to the project directory: {.file {source_metacells_file_path(project)}}")
        }
    } else {
        # Add metacells file symlink
        add_metacells_file(project, anndata_file, copy_source_file)
    }

    save_function_call(command_file_path(project), add_details = TRUE, project = project)
    cli::cli_alert_info("Saving the command to {.file {command_file_path(project)}}")

    cli_alert_success("{.field {dataset}} dataset imported succesfully to {.path {project}} project")
    cli::cli_ul("You can now run the app using: {.field run_app(\"{project}\")}")
    cli::cli_ul("The project is deployment-ready. Just copy the entire project directory to a Shiny server.")
    cli::cli_ul("For more advanced deployment options, you can create a bundle using: {.field create_bundle(\"{project}\", name = \"name_of_bundle\")}")
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

color_cell_types <- function(adata, mc_egc, metacell_types, marker_genes) {
    cli_alert_info("Generating cell type colors using {.pkg chameleon} package.")
    feat_mat <- mc_egc[marker_genes, ]

    if (all(c("u", "v", "w") %in% colnames(adata$obs))) {
        cli_alert_info("Coloring using pre-calculated 3D umap")
        umap_data <- adata$obs %>% select(any_of(c("u", "v", "w")))
        color_of_clusters <- chameleon::data_colors(umap_data, groups = metacell_types$cell_type, run_umap = FALSE)
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
