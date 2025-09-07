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
    mc_egc <- get_mc_egc(dataset, atlas = atlas)

    lfp <- log2(mc_egc + egc_epsilon)
    if (!is.null(exclude)) {
        lfp <- lfp[-which(rownames(lfp) %in% exclude), , drop = FALSE]
    }

    if (!is.null(data_vec)) {
        req(!atlas)
        cors <- tgs_cor(as.matrix(data_vec[colnames(lfp)]), t(lfp), pairwise.complete.obs = TRUE, spearman = FALSE)[1, ]
    } else {
        if (!is.null(metacell_filter)) {
            req(metacell_filter %in% colnames(lfp))
            lfp <- lfp[, metacell_filter, drop = FALSE]
        }
        cors <- tgs_cor(t(as.matrix(lfp[gene, , drop = FALSE])), t(lfp), pairwise.complete.obs = TRUE, spearman = FALSE)[1, ]
    }

    if (type == "pos") {
        cors <- head(sort(cors, decreasing = TRUE), n = 30)
    } else if (type == "neg") {
        cors <- head(sort(cors), n = 30)
    } else {
        cors <- c(head(sort(cors, decreasing = TRUE), n = 30), head(sort(cors), n = 30))
    }
    df <- enframe(cors, "gene2", "cor")

    return(df)
}

get_top_cor_gene <- function(dataset, gene, type = "pos", atlas = FALSE, data_vec = NULL, exclude = NULL, metacell_filter = NULL) {
    gg_mc_top_cor <- get_mc_data(dataset, "gg_mc_top_cor", atlas = atlas)
    if (is.null(data_vec) && is.null(metacell_filter) && gene %in% gg_mc_top_cor$gene1) {
        df <- get_cached_cor(gg_mc_top_cor, gene, type, exclude = exclude)
    } else {
        df <- calc_top_cors(dataset, gene, type, data_vec, metacell_filter, exclude, atlas = atlas)
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

# Correlation export utility functions

#' Calculate correlations for multiple genes individually
#'
#' @param dataset Dataset name
#' @param genes Vector of gene names
#' @param n_top Number of top correlations to return per gene
#' @param cell_type_filter Optional vector of cell types to filter metacells
#' @return Data frame with correlation results
calc_individual_correlations <- function(dataset, genes, n_top = 30, cell_type_filter = NULL, atlas = FALSE) {
    # Get available genes to validate input
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

        # Take top positive and negative correlations
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
#' @param dataset Dataset name
#' @param genes Vector of gene names to treat as module
#' @param n_top Number of top correlations to return
#' @param cell_type_filter Optional vector of cell types to filter metacells
#' @return Data frame with correlation results
calc_module_correlations <- function(dataset, genes, n_top = 30, cell_type_filter = NULL, atlas = FALSE) {
    # Get expression matrix
    mc_egc <- get_mc_egc(dataset, atlas = atlas)
    available_genes <- rownames(mc_egc)
    valid_genes <- genes[genes %in% available_genes]

    if (length(valid_genes) == 0) {
        stop("No valid genes found in dataset")
    }

    if (length(valid_genes) == 1) {
        # Fall back to individual correlation for single gene
        return(calc_individual_correlations(dataset, valid_genes, n_top, cell_type_filter, atlas))
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

    # Calculate module expression as mean of constituent genes
    module_expr <- colMeans(mc_egc[valid_genes, , drop = FALSE], na.rm = TRUE)

    # Use calc_top_cors with the module expression vector
    cors <- calc_top_cors(dataset, "module", "both", module_expr, metacell_filter, valid_genes, atlas = atlas)

    # Format results with truncated gene list for display
    gene_display <- if (length(valid_genes) <= 3) {
        paste(valid_genes, collapse = ",")
    } else {
        paste0(paste(head(valid_genes, 3), collapse = ","), ",...")
    }

    results <- cors %>%
        arrange(desc(abs(cor))) %>%
        head(n_top * 2) %>%
        mutate(
            input_gene = paste0("Module(", gene_display, ")"),
            rank = row_number(),
            correlation_type = ifelse(cor > 0, "positive", "negative")
        )

    return(results)
}

#' Create correlation heatmap for multiple genes
#'
#' @param correlation_results Data frame from calc_individual_correlations
#' @param input_genes Vector of input gene names
#' @param dataset Dataset name
#' @param threshold Correlation threshold for filtering
#' @param cluster Whether to cluster genes
#' @param max_genes Maximum number of genes to include in heatmap
plot_correlation_heatmap <- function(correlation_results, input_genes, dataset, threshold = 0.3,
                                     cluster = TRUE, max_genes = 50, mask_low_correlations = FALSE,
                                     correlation_mode = "individual") {
    # Get top correlated genes across all input genes
    top_genes <- correlation_results %>%
        filter(abs(cor) >= threshold) %>%
        arrange(desc(abs(cor))) %>%
        head(max_genes) %>%
        pull(gene2)

    # Combine input genes and top correlated genes
    # In module mode, always include all the individual module genes in the heatmap
    if (correlation_mode == "module") {
        all_genes <- unique(c(input_genes, top_genes))
    } else {
        all_genes <- unique(c(input_genes, top_genes))
    }

    # Get expression matrix for these genes
    mc_egc <- get_mc_egc(dataset)
    genes_in_data <- intersect(all_genes, rownames(mc_egc))

    if (length(genes_in_data) < 2) {
        # Create simple message plot
        plot.new()
        text(0.5, 0.5, "Not enough genes for heatmap\n(need at least 2 genes)",
            cex = 1.5, adj = c(0.5, 0.5)
        )
        return()
    }

    # Calculate correlation matrix
    expr_subset <- mc_egc[genes_in_data, , drop = FALSE]
    cor_matrix <- cor(t(expr_subset), use = "pairwise.complete.obs")

    # Remove genes that have no correlations above threshold (when threshold > 0)
    if (threshold > 0) {
        # Find genes that have at least one correlation above threshold (excluding self-correlation)
        genes_with_correlations <- apply(cor_matrix, 1, function(row) {
            any(abs(row) >= threshold & abs(row) < 1, na.rm = TRUE) # exclude diagonal (self-correlation = 1)
        })

        # Also keep input genes even if they don't meet threshold
        genes_to_keep <- genes_with_correlations | (rownames(cor_matrix) %in% input_genes)

        if (sum(genes_to_keep) < 2) {
            plot.new()
            text(0.5, 0.5, "Not enough genes above threshold\nfor heatmap",
                cex = 1.5, adj = c(0.5, 0.5)
            )
            return()
        }

        cor_matrix <- cor_matrix[genes_to_keep, genes_to_keep, drop = FALSE]
        genes_in_data <- rownames(cor_matrix)
    }

    # Optionally mask low correlations (make them NA for display)
    cor_matrix_display <- cor_matrix
    if (mask_low_correlations) {
        cor_matrix_display[abs(cor_matrix_display) < threshold] <- NA
    }

    # Create annotation for input genes
    if (correlation_mode == "module") {
        gene_annotation <- data.frame(
            Gene_Type = ifelse(genes_in_data %in% input_genes, "Module", "Correlated"),
            row.names = genes_in_data
        )
    } else {
        gene_annotation <- data.frame(
            Gene_Type = ifelse(genes_in_data %in% input_genes, "Input", "Correlated"),
            row.names = genes_in_data
        )
    }

    # Safety check for clustering - disable clustering if too many NAs
    if (cluster && mask_low_correlations) {
        na_proportion <- sum(is.na(cor_matrix_display)) / length(cor_matrix_display)
        if (na_proportion > 0.5) {
            cluster <- FALSE
            warning("Disabled clustering due to too many masked correlations (>50% NAs)")
        }
    }

    # Additional safety check - ensure matrix is suitable for clustering
    if (cluster) {
        # Check if we have enough non-NA values for clustering
        min_complete_pairs <- apply(cor_matrix_display, 1, function(x) sum(!is.na(x)))
        if (any(min_complete_pairs < 2)) {
            cluster <- FALSE
            warning("Disabled clustering due to insufficient complete correlation pairs")
        }
    }

    # Color palette
    col_fun <- circlize::colorRamp2(
        breaks = c(-1, -0.5, 0, 0.5, 1),
        colors = c("#313695", "#74ADD1", "white", "#F46D43", "#A50026")
    )

    # Annotation colors
    if (correlation_mode == "module") {
        ann_colors <- list(
            Gene_Type = c("Module" = "#FF7F00", "Correlated" = "#1F78B4")
        )
    } else {
        ann_colors <- list(
            Gene_Type = c("Input" = "#E31A1C", "Correlated" = "#1F78B4")
        )
    }

    # Create heatmap
    if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        # Use ComplexHeatmap if available
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
                width = unit(0.5, "cm")
            ),
            top_annotation = ComplexHeatmap::columnAnnotation(
                df = gene_annotation,
                col = ann_colors,
                height = unit(0.5, "cm")
            ),
            heatmap_legend_param = list(
                title = "Correlation",
                legend_direction = "vertical"
            )
        )
    } else {
        # Fallback to pheatmap
        if (requireNamespace("pheatmap", quietly = TRUE)) {
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
                main = paste("Gene Correlations (threshold >=", threshold, ")"),
                fontsize = 10,
                fontsize_row = 8,
                fontsize_col = 8
            )
        } else {
            # Simple base R heatmap
            heatmap(cor_matrix_display,
                main = paste("Gene Correlations (threshold >=", threshold, ")"),
                col = colorRampPalette(c("blue", "white", "red"))(50)
            )
        }
    }
}

#' Create correlation barplot for single gene
#'
#' @param correlation_results Data frame from correlation calculations
#' @param gene Gene name to plot
#' @param threshold Correlation threshold for filtering
plot_correlation_barplot <- function(correlation_results, gene, threshold = 0.3) {
    # Filter and prepare data
    plot_data <- correlation_results %>%
        filter(input_gene == gene | grepl("Module", input_gene)) %>%
        filter(abs(cor) >= threshold) %>%
        arrange(desc(abs(cor))) %>%
        head(30) %>%
        mutate(
            gene2 = factor(gene2, levels = gene2[order(cor)]),
            cor_direction = ifelse(cor > 0, "Positive", "Negative")
        )

    if (nrow(plot_data) == 0) {
        # Create empty plot with message
        p <- ggplot() +
            annotate("text",
                x = 0.5, y = 0.5,
                label = paste("No correlations above threshold", threshold),
                size = 5
            ) +
            xlim(0, 1) +
            ylim(0, 1) +
            theme_void()
        return(plotly::ggplotly(p))
    }

    # Create ggplot
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

    # Convert to plotly with custom tooltip
    plotly::ggplotly(p, tooltip = "text") %>%
        plotly::config(displayModeBar = TRUE, displaylogo = FALSE)
}
