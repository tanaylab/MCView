cell_type_gene_boxplot <- function(gene,
                                   dataset,
                                   cell_types = NULL,
                                   metacell_types = get_mc_data(dataset, "metacell_types"),
                                   cell_type_colors = get_mc_data(dataset, "cell_type_colors"),
                                   egc_gene = NULL,
                                   plot_type = "boxplot",
                                   custom_ylim = NULL,
                                   x_axis_var = "cell_type",
                                   x_axis_categories = NULL,
                                   facet_var = NULL,
                                   coord_flip = FALSE,
                                   log_scale = FALSE) {
    egc_gene <- egc_gene %||% get_gene_egc(gene, dataset) + egc_epsilon

    # Prepare base data
    df <- metacell_types %>%
        mutate(
            !!gene := egc_gene[metacell]
        )

    # Handle x-axis variable
    if (x_axis_var == "cell_type") {
        if (!is.null(cell_types)) {
            df <- df %>% filter(cell_type %in% cell_types)
        }

        if (nrow(df) == 0) {
            return(NULL)
        }

        df <- df %>%
            mutate(cell_type = factor(cell_type, levels = sort(as.character(cell_type_colors$cell_type)))) %>%
            mutate(cell_type = forcats::fct_na_value_to_level(cell_type, "(Missing)")) %>%
            rename(`X axis` = cell_type)

        x_label <- "Cell type"
    } else {
        # Handle metadata x-axis
        metadata <- get_mc_data(dataset, "metadata")
        req(metadata)

        df <- df %>%
            mutate(metacell = as.character(metacell)) %>%
            left_join(metadata %>% select(metacell, !!x_axis_var), by = "metacell")

        # Filter by selected categories
        if (!is.null(x_axis_categories)) {
            df <- df %>% filter(!!sym(x_axis_var) %in% x_axis_categories)
        }

        if (nrow(df) == 0) {
            return(NULL)
        }

        df <- df %>%
            mutate(!!sym(x_axis_var) := factor(!!sym(x_axis_var))) %>%
            mutate(!!sym(x_axis_var) := forcats::fct_na_value_to_level(!!sym(x_axis_var), "(Missing)")) %>%
            rename(`X axis` = !!x_axis_var)

        x_label <- x_axis_var
    }

    # Handle faceting
    if (!is.null(facet_var)) {
        if (facet_var == "cell_type") {
            # Facet by cell type (already have cell_type column)
            df <- df %>%
                mutate(`Facet variable` = factor(cell_type))
            facet_label <- "Cell type"
        } else {
            # Facet by metadata variable
            metadata <- get_mc_data(dataset, "metadata")
            if (!is.null(metadata) && facet_var %in% colnames(metadata)) {
                df <- df %>%
                    mutate(metacell = as.character(metacell)) %>%
                    left_join(metadata %>% select(metacell, !!facet_var), by = "metacell") %>%
                    mutate(`Facet variable` = factor(!!sym(facet_var)))
                facet_label <- facet_var
            } else {
                facet_var <- NULL # Invalid facet variable
            }
        }
    }

    # Handle colors
    if (x_axis_var == "cell_type") {
        col_to_ct <- get_cell_type_colors(dataset, cell_type_colors)
        fill_colors <- col_to_ct
        df <- df %>%
            mutate(`X axis` = factor(`X axis`, levels = names(col_to_ct)))
    } else {
        # For metadata x-axis, use default colors
        categories <- unique(df$`X axis`)
        fill_colors <- rainbow(length(categories))
        names(fill_colors) <- categories
    }

    # Handle y-axis scaling
    if (!is.null(custom_ylim)) {
        # Custom limits override everything
        if (log_scale) {
            y_scale <- scale_y_log10(limits = custom_ylim)
        } else {
            y_scale <- scale_y_continuous(limits = custom_ylim)
        }
    } else if (log_scale) {
        # Log scale without custom limits
        y_scale <- scale_y_log10()
    } else {
        # Default gene expression scaling (log2 with expression breaks)
        ylims <- expr_breaks
        ymax <- min(c(1:length(ylims))[ylims >= max(egc_gene)])
        ymin <- max(c(1:length(ylims))[ylims <= min(egc_gene)])
        y_scale <- scale_y_continuous(limits = c(ylims[ymin], ylims[ymax]), trans = "log2", breaks = ylims[ymin:ymax], labels = scales::scientific(ylims[ymin:ymax]))
    }

    p <- df %>%
        ggplot(aes(x = `X axis`, y = !!sym(gene), fill = `X axis`))

    if (plot_type == "boxplot") {
        p <- p + geom_boxplot()
    } else if (plot_type == "violin") {
        p <- p + geom_violin()
    } else if (plot_type == "sina") {
        p <- p + geom_boxplot() + ggforce::geom_sina()
    }

    # Create y-axis label
    y_label <- if (log_scale) {
        glue("{gene} Expression (log10)")
    } else {
        glue("{gene} Expression")
    }

    p <- p +
        scale_fill_manual(values = fill_colors) +
        y_scale +
        xlab(x_label) +
        ylab(y_label) +
        guides(fill = "none") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # Add faceting if specified
    if (!is.null(facet_var)) {
        p <- p + facet_wrap(~`Facet variable`, scales = "free_x")
    }

    # Add coordinate flipping if specified
    if (coord_flip) {
        p <- p + coord_flip()
    }

    return(p)
}

cell_type_metadata_boxplot <- function(var,
                                       dataset,
                                       cell_types = NULL,
                                       metadata = NULL,
                                       metacell_types = get_mc_data(dataset, "metacell_types"),
                                       cell_type_colors = get_mc_data(dataset, "cell_type_colors"),
                                       plot_type = "boxplot",
                                       custom_ylim = NULL,
                                       x_axis_var = "cell_type",
                                       x_axis_categories = NULL,
                                       facet_var = NULL,
                                       coord_flip = FALSE,
                                       log_scale = FALSE) {
    metadata <- metadata %||% get_mc_data(dataset, "metadata")

    metadata <- metadata %>%
        mutate(metacell = as.character(metacell))

    df <- metacell_types %>%
        select(-any_of(var)) %>%
        left_join(metadata %>% select(metacell, !!var), by = "metacell")

    # Handle x-axis variable
    if (x_axis_var == "cell_type") {
        if (!is.null(cell_types)) {
            df <- df %>% filter(cell_type %in% cell_types)
        }

        if (nrow(df) == 0) {
            return(NULL)
        }

        df <- df %>%
            mutate(cell_type = factor(cell_type, levels = sort(as.character(cell_type_colors$cell_type)))) %>%
            mutate(cell_type = forcats::fct_na_value_to_level(cell_type, "(Missing)")) %>%
            rename(`X axis` = cell_type)

        x_label <- "Cell type"
    } else {
        # Handle metadata x-axis - need to join with x-axis metadata if different from y-axis
        if (x_axis_var != var) {
            df <- df %>%
                left_join(metadata %>% select(metacell, !!x_axis_var), by = "metacell")
        }

        # Filter by selected categories
        if (!is.null(x_axis_categories)) {
            df <- df %>% filter(!!sym(x_axis_var) %in% x_axis_categories)
        }

        if (nrow(df) == 0) {
            return(NULL)
        }

        df <- df %>%
            mutate(!!sym(x_axis_var) := factor(!!sym(x_axis_var))) %>%
            mutate(!!sym(x_axis_var) := forcats::fct_na_value_to_level(!!sym(x_axis_var), "(Missing)")) %>%
            rename(`X axis` = !!x_axis_var)

        x_label <- x_axis_var
    }

    # Handle faceting
    if (!is.null(facet_var)) {
        if (facet_var == "cell_type") {
            # Facet by cell type (already have cell_type column)
            df <- df %>%
                mutate(`Facet variable` = factor(cell_type))
            facet_label <- "Cell type"
        } else {
            # Facet by metadata variable
            if (!is.null(metadata) && facet_var %in% colnames(metadata)) {
                # Need to join facet metadata if different from x-axis or y-axis
                if (facet_var != var && (x_axis_var == "cell_type" || facet_var != x_axis_var)) {
                    df <- df %>%
                        left_join(metadata %>% select(metacell, !!facet_var), by = "metacell")
                }
                df <- df %>%
                    mutate(`Facet variable` = factor(!!sym(facet_var)))
                facet_label <- facet_var
            } else {
                facet_var <- NULL # Invalid facet variable
            }
        }
    }

    # Handle colors
    if (x_axis_var == "cell_type") {
        col_to_ct <- get_cell_type_colors(dataset, cell_type_colors)
        fill_colors <- col_to_ct
        df <- df %>%
            mutate(`X axis` = factor(`X axis`, levels = names(col_to_ct)))
    } else {
        # For metadata x-axis, use default colors
        categories <- unique(df$`X axis`)
        fill_colors <- rainbow(length(categories))
        names(fill_colors) <- categories
    }

    p <- df %>%
        ggplot(aes(x = `X axis`, y = !!sym(var), fill = `X axis`))

    if (plot_type == "boxplot") {
        p <- p + geom_boxplot()
    } else if (plot_type == "violin") {
        p <- p + geom_violin()
    } else if (plot_type == "sina") {
        p <- p + geom_boxplot() + ggforce::geom_sina()
    }

    # Handle y-axis scaling
    if (!is.null(custom_ylim)) {
        # Custom limits override everything
        if (log_scale) {
            y_scale <- scale_y_log10(limits = custom_ylim)
        } else {
            y_scale <- scale_y_continuous(limits = custom_ylim)
        }
    } else if (log_scale) {
        # Log scale without custom limits
        y_scale <- scale_y_log10()
    } else {
        # Default scaling
        y_scale <- NULL
    }

    # Create y-axis label
    y_label <- if (log_scale) {
        glue("{var} (log10)")
    } else {
        var
    }

    p <- p +
        scale_fill_manual(values = fill_colors) +
        xlab(x_label) +
        ylab(y_label) +
        guides(fill = "none") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    if (!is.null(y_scale)) {
        p <- p + y_scale
    }

    # Add faceting if specified
    if (!is.null(facet_var)) {
        p <- p + facet_wrap(~`Facet variable`, scales = "free_x")
    }

    # Add coordinate flipping if specified
    if (coord_flip) {
        p <- p + coord_flip()
    }

    return(p)
}

cell_type_metadata_confusion <- function(var,
                                         dataset,
                                         cell_types = NULL,
                                         color_by = "Cell type",
                                         metadata = NULL,
                                         metacell_types = get_mc_data(dataset, "metacell_types"),
                                         cell_type_colors = get_mc_data(dataset, "cell_type_colors")) {
    metadata <- metadata %||% get_mc_data(dataset, "metadata")

    metadata <- metadata %>%
        mutate(metacell = as.character(metacell))

    df <- metacell_types %>%
        select(-any_of(var)) %>%
        left_join(metadata %>% select(metacell, !!var), by = "metacell")

    if (!is.null(cell_types)) {
        df <- df %>% filter(cell_type %in% cell_types)
    }

    if (nrow(df) == 0) {
        return(NULL)
    }

    df <- df %>%
        mutate(cell_type = factor(cell_type, levels = sort(as.character(cell_type_colors$cell_type)))) %>%
        mutate(cell_type = forcats::fct_na_value_to_level(cell_type, "(Missing)")) %>%
        mutate(var = factor(var), !!var := forcats::fct_na_value_to_level(!!sym(var), "(Missing)"))

    df <- df %>%
        count(cell_type, !!sym(var), .drop = FALSE) %>%
        group_by(cell_type) %>%
        mutate(n_tot = sum(n), p_cell_type = n / n_tot) %>%
        ungroup() %>%
        group_by(!!sym(var)) %>%
        mutate(n_tot_md = sum(n), p_md = n / n_tot_md) %>%
        ungroup() %>%
        tidyr::replace_na(replace = list(p_cell_type = 0, p_var = 0)) %>%
        mutate(
            `# of metacells` = n,
            `total # of cell type metacells` = n_tot,
            !!sym(glue("total # of {var} metacells")) := n_tot_md,
            `% of cell type metacells` = glue("{scales::percent(p_cell_type)} ({n}/{n_tot})"),
            !!sym(glue("% of {var} metacells")) := glue("{scales::percent(p_md)} ({n}/{n_tot_md})")
        ) %>%
        rename(`Cell type` = cell_type)

    if (color_by == "X axis") {
        df <- df %>% mutate(color = p_cell_type)
        label <- "% Cell type"
    } else {
        df <- df %>% mutate(color = p_md)
        label <- glue("% {var}")
    }

    fracs_mat <- df %>%
        distinct(`Cell type`, !!sym(var), color) %>%
        spread(`Cell type`, color) %>%
        column_to_rownames(var) %>%
        as.matrix()

    ord <- slanter::slanted_orders(fracs_mat)

    df <- df %>%
        mutate(`Cell type` = factor(`Cell type`, levels = colnames(fracs_mat)[ord$cols])) %>%
        mutate(!!sym(var) := factor(!!sym(var), levels = rownames(fracs_mat)[ord$rows]))

    p <- df %>%
        ggplot(aes(
            x = `Cell type`,
            y = !!sym(var),
            lab1 = `# of metacells`,
            lab2 = `total # of cell type metacells`,
            lab3 = !!sym(glue("total # of {var} metacells")),
            lab4 = `% of cell type metacells`,
            lab5 = !!sym(glue("% of {var} metacells")),
            fill = color
        )) +
        geom_tile() +
        scale_fill_gradientn(
            colors = c("white", "#F4A582", "#D6604D", "#B2182B", "#67001F", "black"),
            name = label,
            limits = c(0, 1)
        ) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    return(p)
}
