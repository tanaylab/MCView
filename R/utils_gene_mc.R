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
    req(gene %in% rownames(mc_egc))
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
