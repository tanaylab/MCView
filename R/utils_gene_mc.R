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
            dafr::And("gene1") %>%
            dafr::IsEqual(gene)

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
    # Push metacell filter down to DAF query to avoid loading full matrix
    mc_egc <- get_mc_egc(dataset, atlas = atlas,
        metacells = if (is.null(data_vec)) metacell_filter else NULL)

    # Densify for correlation computation (tgs_cor requires plain numeric matrix)
    lfp <- as.matrix(log2(mc_egc + mcv_get("egc_epsilon")))
    if (!is.null(exclude)) {
        exclude_idx <- which(rownames(lfp) %in% exclude)
        if (length(exclude_idx) > 0) {
            lfp <- lfp[-exclude_idx, , drop = FALSE]
        }
    }

    if (!is.null(data_vec)) {
        req(!atlas)
        x_mat <- as.matrix(data_vec[colnames(lfp)])
    } else {
        x_mat <- t(lfp[gene, , drop = FALSE])
    }
    y_mat <- t(lfp)

    cor_mat <- tryCatch(
        .Call(
            "tgs_cross_cor_blas",
            x_mat,
            y_mat,
            TRUE,
            FALSE,
            FALSE,
            0,
            new.env(parent = parent.frame())
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

get_top_cor_gene <- function(dataset, gene, type = "pos", atlas = FALSE, data_vec = NULL, exclude = NULL, metacell_filter = NULL) {
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
