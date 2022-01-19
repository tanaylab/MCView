cell_type_gene_boxplot <- function(gene,
                                   dataset,
                                   cell_types = NULL,
                                   metacell_types = get_mc_data(dataset, "metacell_types"),
                                   cell_type_colors = get_mc_data(dataset, "cell_type_colors")) {
    egc_gene <- get_gene_egc(gene, dataset) + egc_epsilon

    df <- metacell_types %>%
        mutate(
            !!gene := egc_gene[metacell]
        )

    if (!is.null(cell_types)) {
        df <- df %>% filter(cell_type %in% cell_types)
    }

    if (nrow(df) == 0) {
        return(NULL)
    }

    df <- df %>%
        mutate(cell_type = factor(cell_type, levels = sort(as.character(cell_type_colors$cell_type)))) %>%
        mutate(cell_type = forcats::fct_explicit_na(cell_type)) %>%
        rename(
            `Cell type` = cell_type
        )

    col_to_ct <- get_cell_type_colors(dataset, cell_type_colors)

    ylims <- expr_breaks
    ymax <- min(c(1:length(ylims))[ylims >= max(egc_gene)])
    ymin <- max(c(1:length(ylims))[ylims <= min(egc_gene)])

    p <- df %>%
        ggplot(aes(x = `Cell type`, y = !!sym(gene), fill = `Cell type`)) +
        geom_boxplot() +
        scale_fill_manual(values = col_to_ct) +
        scale_y_continuous(limits = c(ylims[ymin], ylims[ymax]), trans = "log2", breaks = ylims[ymin:ymax], labels = scales::scientific(ylims[ymin:ymax])) +
        xlab("") +
        ylab(glue("{gene} Expression")) +
        guides(fill = "none") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    return(p)
}

cell_type_metadata_boxplot <- function(var,
                                       dataset,
                                       cell_types = NULL,
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
        mutate(cell_type = forcats::fct_explicit_na(cell_type)) %>%
        rename(
            `Cell type` = cell_type
        )

    col_to_ct <- get_cell_type_colors(dataset, cell_type_colors)

    p <- df %>%
        ggplot(aes(x = `Cell type`, y = !!sym(var), fill = `Cell type`)) +
        geom_boxplot() +
        scale_fill_manual(values = col_to_ct) +
        xlab("") +
        ylab(var) +
        guides(fill = "none") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

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
        mutate(cell_type = forcats::fct_explicit_na(cell_type)) %>%
        mutate(var = factor(var), !!var := forcats::fct_explicit_na(!!sym(var)))

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
            cell_type_stats = glue("{scales::percent(p_cell_type)} ({n}/{n_tot})"),
            md_stats = glue("{scales::percent(p_md)} ({n}/{n_tot_md})")
        ) %>%
        rename(`Cell type` = cell_type)

    if (color_by == "Cell type") {
        color_var <- "p_cell_type"
        label <- "% Type"
    } else {
        color_var <- "p_md"
        label <- glue("% {var}")
    }

    p <- df %>%
        ggplot(aes(x = `Cell type`, y = !!sym(var), cell_type_stats = cell_type_stats, md_stats = md_stats, fill = !!sym(color_var))) +
        geom_tile() +
        scale_fill_gradientn(
            colors = c("white", "#F4A582", "#D6604D", "#B2182B", "#67001F", "black"),
            name = label,
            limits = c(0, 1)
        ) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    return(p)
}
