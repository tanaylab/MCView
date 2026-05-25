get_cell_type_colors <- function(dataset, cell_type_colors = NULL, na_color = "gray", atlas = FALSE) {
    if (is.null(cell_type_colors)) {
        cell_type_colors <- get_mc_data(dataset, "cell_type_colors", atlas = atlas)
    }

    if (!("(Missing)" %in% cell_type_colors$cell_type)) {
        res <- bind_rows(
            tibble(cell_type = "(Missing)", color = na_color),
            cell_type_colors %>%
                distinct(cell_type, color) %>%
                select(cell_type, color)
        ) %>%
            deframe()
    } else {
        res <- cell_type_colors %>%
            distinct(cell_type, color) %>%
            select(cell_type, color) %>%
            deframe()
    }

    return(res)
}

get_cached_cor_daf <- function(daf_obj, gene, type, exclude = NULL) {
    if (is.null(daf_obj)) {
        return(NULL)
    }

    get_cached_axis <- function(axis_name) {
        if (!dafr::has_axis(daf_obj, axis_name)) {
            return(NULL)
        }

        base_q <- dafr::Axis(axis_name) %>%
            dafr::BeginMask("gene1") %>%
            dafr::IsEqual(gene) %>%
            dafr::EndMask()

        df <- tryCatch(
            dafr::get_dataframe(daf_obj, base_q, columns = c("gene2", "cor", "type"), cache = TRUE),
            error = function(e) NULL
        )
        if (is.null(df) || nrow(df) == 0 || !all(c("gene2", "cor", "type") %in% names(df))) {
            return(NULL)
        }

        df <- tibble::as_tibble(df) %>%
            mutate(
                gene2 = as.character(gene2),
                cor = as.numeric(cor),
                type = as.character(type)
            )

        if (type %in% c("pos", "neg")) {
            df <- df %>%
                filter(type == !!type)
        }

        df <- df %>%
            filter(gene2 != gene)

        if (!is.null(exclude)) {
            df <- df %>%
                filter(!(gene2 %in% exclude))
        }

        df
    }

    axis_candidates <- c("mcview_cache_gg_mc_top_cor", "gg_mc_top_cor")
    for (axis_name in axis_candidates) {
        df <- get_cached_axis(axis_name)
        if (!is.null(df) && nrow(df) > 0) {
            return(df)
        }
    }

    NULL
}

get_cached_cor <- function(gg_mc_top_cor, gene, type, exclude = NULL) {
    df <- gg_mc_top_cor %>%
        filter(gene1 == gene)

    if (type %in% c("pos", "neg")) {
        df <- df %>%
            filter(type == !!type)
    }

    df <- df %>%
        filter(gene1 != gene2)

    if (!is.null(exclude)) {
        df <- df %>%
            filter(!(gene2 %in% exclude))
    }
    return(df)
}

calc_top_cors <- function(dataset, gene, type, data_vec, metacell_filter, exclude, atlas = FALSE) {
    # BLAS-accelerated correlation.
    # metacell_filter restricts the metacell set the correlation runs over,
    # on BOTH the gene path and the data_vec (module / metadata) path. When a
    # filter is set, x_mat is matched against the filtered colnames(lfp)
    # below, so data_vec is automatically subset to the same metacells - the
    # filtered lfp and the filtered data_vec stay aligned by name.
    # (Historically the data_vec path silently dropped metacell_filter, which
    # made the Gene Correlation tab's cell-type filter a no-op in module mode
    # while single-gene mode honoured it - an inconsistency, now removed.)
    effective_filter <- metacell_filter
    # The cached lfp_full and its paired transpose apply whenever no metacell
    # subset applies and we're not in atlas mode - i.e. the common no-filter
    # case (gene path or data_vec path). A filter forces a fresh filtered lfp.
    use_cached_lfp <- is.null(effective_filter) && !atlas
    if (use_cached_lfp) {
        lfp <- get_mc_lfp(dataset)
    } else {
        mc_egc <- get_mc_egc(dataset, atlas = atlas, metacells = effective_filter)
        lfp <- dafr::fast_log(mc_egc, eps = mcv_get("egc_epsilon"), base = 2)
    }
    # exclude only affects WHICH genes appear in the cor output - it doesn't
    # change x_mat (gene-row extract or data_vec by metacell name) or the
    # underlying lfp values, so we defer the actual row-removal until we
    # know which materialisation we need. For the cached-transpose path we
    # can column-subset full_t directly; for the t(lfp) fallback path we
    # row-subset lfp at that point. This avoids the ~250 ms row-subset of
    # the 28K-row lfp on every Gene Modules click (where exclude is the
    # module's own 16-ish genes - the subset cost is dominated by the
    # SIZE of the result, not the size of exclude).
    exclude_idx <- if (!is.null(exclude)) which(rownames(lfp) %in% exclude) else integer(0)
    excluded <- length(exclude_idx) > 0

    if (!is.null(data_vec)) {
        req(!atlas)
        x_mat <- as.matrix(data_vec[colnames(lfp)])
    } else {
        x_mat <- t(lfp[gene, , drop = FALSE])
    }
    # y_mat is metacell-rows x gene-cols. With the cached lfp_full_t we
    # can serve both the unfiltered case (whole-matrix reuse) and the
    # exclude case (column-subset) without ever materialising t(lfp) on
    # a row-subsetted lfp - column-subsetting the cached transpose is
    # ~4x faster than subset-then-transpose on OBK (~250 vs ~1000 ms).
    y_mat <- if (use_cached_lfp) {
        full_t <- get_mc_lfp_t(dataset)
        if (excluded) full_t[, -exclude_idx, drop = FALSE] else full_t
    } else if (excluded) {
        t(lfp[-exclude_idx, , drop = FALSE])
    } else {
        t(lfp)
    }

    cor_mat <- tryCatch(
        .Call(
            "tgs_cross_cor_blas",
            x_mat,
            y_mat,
            TRUE,
            FALSE,
            FALSE,
            0,
            new.env(parent = parent.frame()),
            PACKAGE = "tgstat"
        ),
        error = function(e) {
            prev_blas <- getOption("tgs_use.blas")
            if (!isTRUE(prev_blas)) {
                options(tgs_use.blas = TRUE)
                on.exit(options(tgs_use.blas = prev_blas), add = TRUE)
            }
            tgs_cor(x_mat, y = y_mat, pairwise.complete.obs = TRUE, spearman = FALSE)
        }
    )
    cors <- cor_mat[1, ]

    if (type == "pos") {
        pos <- head(sort(cors, decreasing = TRUE), n = 30)
        df <- enframe(pos, "gene2", "cor") %>%
            mutate(type = "pos")
    } else if (type == "neg") {
        neg <- head(sort(cors), n = 30)
        df <- enframe(neg, "gene2", "cor") %>%
            mutate(type = "neg")
    } else {
        pos <- head(sort(cors, decreasing = TRUE), n = 30)
        neg <- head(sort(cors), n = 30)
        df <- bind_rows(
            enframe(pos, "gene2", "cor") %>% mutate(type = "pos"),
            enframe(neg, "gene2", "cor") %>% mutate(type = "neg")
        )
    }

    return(df)
}

get_top_cor_gene <- function(dataset, gene, type = "pos", atlas = FALSE, data_vec = NULL, exclude = NULL, metacell_filter = NULL, markers_only = FALSE) {
    df <- NULL
    df_from_calc <- FALSE
    cache_used <- FALSE
    cache_key <- paste0(as.character(gene), "::", as.character(type))

    if (is.null(data_vec) && is.null(metacell_filter)) {
        mc_data <- mcv_get("mc_data")
        if (!is.null(mc_data[[dataset]][["top_cor_genes"]]) && has_name(mc_data[[dataset]][["top_cor_genes"]], cache_key)) {
            df <- mc_data[[dataset]][["top_cor_genes"]][[cache_key]]
            cache_used <- TRUE
        }
    }

    if (is.null(data_vec) && is.null(metacell_filter)) {
        daf_obj <- if (atlas) get_atlas_daf() else get_dataset_daf(dataset)
        if (is.null(df)) {
            df <- tryCatch(
                get_cached_cor_daf(daf_obj, gene, type = type, exclude = NULL),
                error = function(e) NULL
            )
            if (!is.null(df)) {
                mc_data <- mcv_get("mc_data")
                mc_data[[dataset]][["top_cor_genes"]][[cache_key]] <- df
                mcv_set("mc_data", mc_data)
            }
        }

        if (is.null(df)) {
            gg_mc_top_cor <- get_mc_data(dataset, "gg_mc_top_cor", atlas = atlas)
            if (!is.null(gg_mc_top_cor) && gene %in% gg_mc_top_cor$gene1) {
                df <- get_cached_cor(gg_mc_top_cor, gene, type = type, exclude = NULL)
                mc_data <- mcv_get("mc_data")
                mc_data[[dataset]][["top_cor_genes"]][[cache_key]] <- df
                mcv_set("mc_data", mc_data)
            }
        }
    }

    if (is.null(df)) {
        df <- calc_top_cors(dataset, gene, type, data_vec, metacell_filter, exclude, atlas = atlas)
        df_from_calc <- TRUE
    }

    if (df_from_calc && is.null(data_vec) && is.null(metacell_filter)) {
        daf_obj <- if (atlas) get_atlas_daf() else get_dataset_daf(dataset)
        cache_correlations_daf(daf_obj, gene, df)
        if (!cache_used) {
            mc_data <- mcv_get("mc_data")
            mc_data[[dataset]][["top_cor_genes"]][[cache_key]] <- df
            mcv_set("mc_data", mc_data)
        }
    }

    if (type %in% c("pos", "neg")) {
        df <- df %>%
            filter(type == !!type)
    }

    if (!is.null(exclude)) {
        df <- df %>%
            filter(!(gene2 %in% exclude))
    }

    if (isTRUE(markers_only)) {
        marker_genes <- tryCatch({
            get_mc_data(dataset, "marker_genes")$gene
        }, error = function(e) NULL)
        if (!is.null(marker_genes) && length(marker_genes) > 0) {
            df <- df %>%
                filter(gene2 %in% marker_genes)
        }
    }

    res <- df %>%
        mutate(
            gene2 = as.character(gene2),
            gene_name = gene_label(gene2, dataset),
            label = glue("{gene_name} ({cr})", cr = round(cor, digits = 2))
        ) %>%
        select(label, gene2) %>%
        tibble::deframe()

    return(res)
}

# ==============================================================================
# Correlation helpers for gene_correlation module
# ==============================================================================

#' Calculate correlations for multiple genes individually
#'
#' For each input gene, finds the top positively and negatively correlated
#' genes across all metacells (optionally filtered by cell type).
#'
#' @param dataset Dataset name
#' @param genes Vector of gene names
#' @param n_top Number of top correlations to return per gene
#' @param cell_type_filter Optional vector of cell types to filter metacells
#' @param threshold Minimum absolute correlation to include
#' @param atlas Logical; use atlas data
#' @return Data frame with columns: gene2, cor, input_gene, rank, correlation_type
calc_individual_correlations <- function(dataset, genes, n_top = 30, cell_type_filter = NULL, threshold = 0, atlas = FALSE) {
    available_genes <- gene_names(dataset, atlas = atlas)
    valid_genes <- genes[genes %in% available_genes]

    if (length(valid_genes) == 0) {
        stop("No valid genes found in dataset")
    }

    # Calculate metacell filter if cell types specified
    metacell_filter <- NULL
    if (!is.null(cell_type_filter)) {
        metacell_types <- get_mc_data(dataset, "metacell_types")
        if (!is.null(metacell_types)) {
            metacell_filter <- metacell_types %>%
                filter(cell_type %in% cell_type_filter) %>%
                pull(metacell)
        }
    }

    # Calculate correlations for each gene
    results <- purrr::map_dfr(valid_genes, function(gene) {
        cors <- calc_top_cors(dataset, gene, "both", NULL, metacell_filter, NULL, atlas = atlas)

        # Filter by threshold and take top correlations
        cors <- cors %>%
            filter(abs(cor) >= threshold)

        pos_cors <- cors %>%
            filter(cor > 0) %>%
            arrange(desc(cor)) %>%
            head(n_top)
        neg_cors <- cors %>%
            filter(cor < 0) %>%
            arrange(cor) %>%
            head(n_top)

        combined_cors <- bind_rows(pos_cors, neg_cors) %>%
            arrange(desc(abs(cor))) %>%
            head(n_top * 2) %>%
            mutate(
                input_gene = gene,
                rank = row_number(),
                correlation_type = ifelse(cor > 0, "positive", "negative")
            )

        return(combined_cors)
    })

    return(results)
}

#' Calculate correlations for genes as a module/anchor
#'
#' Treats a set of genes as a single module by averaging their expression,
#' then finds the top correlated genes with that combined pattern.
#'
#' @param dataset Dataset name
#' @param genes Vector of gene names to treat as module
#' @param n_top Number of top correlations to return
#' @param cell_type_filter Optional vector of cell types to filter metacells
#' @param threshold Minimum absolute correlation to include
#' @param atlas Logical; use atlas data
#' @return Data frame with columns: gene2, cor, input_gene, rank, correlation_type
calc_module_correlations <- function(dataset, genes, n_top = 30, cell_type_filter = NULL, threshold = 0, atlas = FALSE) {
    available_genes <- gene_names(dataset, atlas = atlas)
    valid_genes <- genes[genes %in% available_genes]

    if (length(valid_genes) == 0) {
        stop("No valid genes found in dataset")
    }

    if (length(valid_genes) == 1) {
        # Fall back to individual correlation for single gene
        return(calc_individual_correlations(dataset, valid_genes, n_top, cell_type_filter, threshold, atlas))
    }

    # Calculate metacell filter if cell types specified
    metacell_filter <- NULL
    if (!is.null(cell_type_filter)) {
        metacell_types <- get_mc_data(dataset, "metacell_types")
        if (!is.null(metacell_types)) {
            metacell_filter <- metacell_types %>%
                filter(cell_type %in% cell_type_filter) %>%
                pull(metacell)
        }
    }

    # Module expression = mean fraction over the gene set per metacell.
    # Pull only the requested gene rows (mask query) instead of the full
    # gene x metacell matrix.
    mc_egc <- get_mc_egc(dataset, genes = valid_genes, atlas = atlas)
    module_expr <- colMeans(mc_egc, na.rm = TRUE)

    # Use calc_top_cors with the module expression vector
    cors <- calc_top_cors(dataset, "module", "both", module_expr, metacell_filter, valid_genes, atlas = atlas)

    # Format results with truncated gene list for display
    gene_display <- if (length(valid_genes) <= 3) {
        paste(valid_genes, collapse = ",")
    } else {
        paste0(paste(head(valid_genes, 3), collapse = ","), ",...")
    }

    results <- cors %>%
        filter(abs(cor) >= threshold) %>%
        arrange(desc(abs(cor))) %>%
        head(n_top * 2) %>%
        mutate(
            input_gene = paste0("Module(", gene_display, ")"),
            rank = row_number(),
            correlation_type = ifelse(cor > 0, "positive", "negative")
        )

    return(results)
}

#' Calculate gene-gene correlations (correlation between input genes)
#'
#' Computes all pairwise correlations among the supplied genes.
#' Useful for exploring the internal structure of a gene set.
#'
#' @param dataset Dataset name
#' @param genes Vector of gene names
#' @param cell_type_filter Optional vector of cell types to filter metacells
#' @param threshold Minimum absolute correlation to include
#' @param atlas Logical; use atlas data
#' @return Data frame with columns: input_gene, gene2, cor, correlation_type, rank
calc_gene_gene_correlations <- function(dataset, genes, cell_type_filter = NULL, threshold = 0, atlas = FALSE) {
    available_genes <- gene_names(dataset, atlas = atlas)
    valid_genes <- genes[genes %in% available_genes]

    if (length(valid_genes) < 2) {
        stop("Need at least 2 valid genes for gene-gene correlation")
    }

    # Calculate metacell filter if cell types specified
    metacell_filter <- NULL
    if (!is.null(cell_type_filter)) {
        metacell_types <- get_mc_data(dataset, "metacell_types")
        if (!is.null(metacell_types)) {
            metacell_filter <- metacell_types %>%
                filter(cell_type %in% cell_type_filter) %>%
                pull(metacell)
        }
    }

    # Pull only the requested gene rows; mask query keeps the full gene x
    # metacell matrix off the heap.
    expr_data <- get_mc_egc(dataset, genes = valid_genes, atlas = atlas)
    if (!is.null(metacell_filter)) {
        metacell_filter <- intersect(metacell_filter, colnames(expr_data))
        expr_data <- expr_data[, metacell_filter, drop = FALSE]
    }

    # Calculate all pairwise correlations between genes (BLAS-accelerated)
    t_expr <- t(expr_data)
    cor_matrix <- tryCatch(
        .Call("tgs_cross_cor_blas", t_expr, t_expr, TRUE, FALSE, FALSE, 0, new.env(parent = parent.frame()), PACKAGE = "tgstat"),
        error = function(e) tgs_cor(t_expr, pairwise.complete.obs = TRUE, spearman = FALSE)
    )

    # Convert to long format, excluding self-correlations
    results <- expand.grid(
        input_gene = valid_genes,
        gene2 = valid_genes,
        stringsAsFactors = FALSE
    ) %>%
        filter(input_gene != gene2) %>%
        mutate(
            cor = mapply(function(g1, g2) cor_matrix[g1, g2], input_gene, gene2),
            correlation_type = ifelse(cor > 0, "positive", "negative")
        ) %>%
        filter(abs(cor) >= threshold) %>%
        mutate(
            rank = ave(abs(cor), input_gene, FUN = function(x) rank(-x))
        ) %>%
        arrange(input_gene, desc(abs(cor)))

    return(results)
}

#' Create correlation heatmap for multiple genes
#'
#' Renders a heatmap of the correlation matrix for the top correlated
#' genes.
#'
#' Tries ComplexHeatmap first, falls back to pheatmap, then to base R
#' \code{heatmap()}.  All three packages are loaded conditionally via
#' \code{requireNamespace()}.
#'
#' @param correlation_results Data frame from correlation calculations
#' @param input_genes Vector of input gene names
#' @param dataset Dataset name
#' @param cluster Whether to cluster genes
#' @param max_genes Maximum number of genes in heatmap
#' @param mask_low_correlations Mask weak correlations
#' @param correlation_mode "individual" or "module"
plot_correlation_heatmap <- function(correlation_results, input_genes, dataset,
                                     cluster = TRUE, max_genes = 50,
                                     mask_low_correlations = FALSE,
                                     correlation_mode = "individual") {
    # Get top correlated genes across all input genes
    top_genes <- correlation_results %>%
        arrange(desc(abs(cor))) %>%
        head(max_genes) %>%
        pull(gene2)

    # Combine input genes and top correlated genes
    all_genes <- unique(c(input_genes, top_genes))

    # Get expression matrix for the small target gene set only.
    genes_in_data <- intersect(all_genes, gene_names(dataset))

    if (length(genes_in_data) < 2) {
        plot.new()
        text(0.5, 0.5, "Not enough genes for heatmap\n(need at least 2 genes)",
            cex = 1.5, adj = c(0.5, 0.5)
        )
        return(invisible(NULL))
    }

    # Calculate correlation matrix (BLAS-accelerated)
    expr_subset <- get_mc_egc(dataset, genes = genes_in_data)
    t_expr <- t(expr_subset)
    cor_matrix <- tryCatch(
        .Call("tgs_cross_cor_blas", t_expr, t_expr, TRUE, FALSE, FALSE, 0, new.env(parent = parent.frame()), PACKAGE = "tgstat"),
        error = function(e) tgs_cor(t_expr, pairwise.complete.obs = TRUE, spearman = FALSE)
    )

    cor_matrix_display <- cor_matrix

    # Mask low correlations if requested
    if (mask_low_correlations) {
        cor_matrix_display[abs(cor_matrix_display) < 0.3] <- NA
    }

    # Create annotation for input genes
    gene_type_label <- if (correlation_mode == "module") "Module" else "Input"
    gene_annotation <- data.frame(
        Gene_Type = ifelse(genes_in_data %in% input_genes, gene_type_label, "Correlated"),
        row.names = genes_in_data
    )

    # Safety: disable clustering if too many NAs
    if (cluster && mask_low_correlations) {
        na_proportion <- sum(is.na(cor_matrix_display)) / length(cor_matrix_display)
        if (na_proportion > 0.5) {
            cluster <- FALSE
            warning("Disabled clustering due to too many masked correlations (>50% NAs)")
        }
    }
    if (cluster) {
        min_complete_pairs <- apply(cor_matrix_display, 1, function(x) sum(!is.na(x)))
        if (any(min_complete_pairs < 2)) {
            cluster <- FALSE
            warning("Disabled clustering due to insufficient complete correlation pairs")
        }
    }

    # Colour palette
    col_fun <- circlize::colorRamp2(
        breaks = c(-1, -0.5, 0, 0.5, 1),
        colors = c("#313695", "#74ADD1", "white", "#F46D43", "#A50026")
    )

    # Annotation colours
    if (correlation_mode == "module") {
        ann_colors <- list(Gene_Type = c("Module" = "#FF7F00", "Correlated" = "#1F78B4"))
    } else {
        ann_colors <- list(Gene_Type = c("Input" = "#E31A1C", "Correlated" = "#1F78B4"))
    }

    # Render heatmap -- prefer ComplexHeatmap > pheatmap > base
    if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        ComplexHeatmap::Heatmap(
            cor_matrix_display,
            name = "Correlation",
            col = col_fun,
            cluster_rows = cluster,
            cluster_columns = cluster,
            show_row_names = TRUE,
            show_column_names = TRUE,
            row_names_gp = grid::gpar(fontsize = 10),
            column_names_gp = grid::gpar(fontsize = 10),
            left_annotation = ComplexHeatmap::rowAnnotation(
                df = gene_annotation,
                col = ann_colors,
                width = grid::unit(0.5, "cm")
            ),
            top_annotation = ComplexHeatmap::columnAnnotation(
                df = gene_annotation,
                col = ann_colors,
                height = grid::unit(0.5, "cm")
            ),
            heatmap_legend_param = list(
                title = "Correlation",
                legend_direction = "vertical"
            )
        )
    } else if (requireNamespace("pheatmap", quietly = TRUE)) {
        pheatmap::pheatmap(
            cor_matrix_display,
            color = colorRampPalette(c("#313695", "#74ADD1", "white", "#F46D43", "#A50026"))(100),
            breaks = seq(-1, 1, length.out = 101),
            na_col = "grey90",
            cluster_rows = cluster,
            cluster_cols = cluster,
            annotation_row = gene_annotation,
            annotation_col = gene_annotation,
            annotation_colors = ann_colors,
            main = "Gene Correlations",
            fontsize = 10,
            fontsize_row = 8,
            fontsize_col = 8
        )
    } else {
        # Minimal base R fallback
        heatmap(cor_matrix_display,
            main = "Gene Correlations",
            col = colorRampPalette(c("blue", "white", "red"))(50)
        )
    }
}

#' Create correlation barplot for single gene
#'
#' Produces a horizontal bar chart of top correlations for a single gene
#' (or a module anchor) via ggplot2 + plotly.
#'
#' @param correlation_results Data frame from correlation calculations
#' @param gene Gene name to plot
#' @return A plotly htmlwidget
plot_correlation_barplot <- function(correlation_results, gene) {
    plot_data <- correlation_results %>%
        filter(input_gene == gene | grepl("Module", input_gene)) %>%
        arrange(desc(abs(cor))) %>%
        head(30) %>%
        mutate(
            gene2 = factor(gene2, levels = gene2[order(cor)]),
            cor_direction = ifelse(cor > 0, "Positive", "Negative")
        )

    if (nrow(plot_data) == 0) {
        p <- ggplot() +
            annotate("text",
                x = 0.5, y = 0.5,
                label = "No correlations found",
                size = 5
            ) +
            xlim(0, 1) +
            ylim(0, 1) +
            theme_void()
        return(plotly::ggplotly(p))
    }

    p <- ggplot(plot_data, aes(
        x = gene2, y = cor, fill = cor_direction,
        text = paste(
            "Gene:", gene2, "<br>",
            "Correlation:", round(cor, 3), "<br>",
            "Rank:", rank
        )
    )) +
        geom_col(alpha = 0.8) +
        coord_flip() +
        scale_fill_manual(
            values = c("Positive" = "#E31A1C", "Negative" = "#1F78B4"),
            name = "Direction"
        ) +
        labs(
            title = paste("Top Correlations with", gene),
            x = "Gene",
            y = "Correlation Coefficient"
        ) +
        theme_minimal() +
        theme(
            axis.text.y = element_text(size = 10),
            plot.title = element_text(size = 14, hjust = 0.5),
            legend.position = "bottom"
        )

    plotly::ggplotly(p, tooltip = "text") %>%
        plotly::config(displayModeBar = TRUE, displaylogo = FALSE)
}
