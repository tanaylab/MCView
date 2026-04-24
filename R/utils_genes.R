annotate_genes <- function(genes, dataset) {
    lateral_genes <- get_mc_data(dataset, "lateral_genes")
    noisy_genes <- get_mc_data(dataset, "noisy_genes")

    case_when(
        (genes %in% lateral_genes) & (genes %in% noisy_genes) ~ "lateral, noisy",
        genes %in% lateral_genes ~ "lateral",
        genes %in% noisy_genes ~ "noisy",
        TRUE ~ "other"
    )
}

gene_label <- function(genes, dataset, gene_modules = NULL) {
    type <- annotate_genes(genes, dataset)
    if (!is.null(gene_modules)) {
        modules <- tibble(gene = genes) %>%
            left_join(gene_modules, by = join_by(gene)) %>%
            pull(module) %>%
            as.character()
    }
    new_names <- ifelse(type == "other", genes, glue("{genes} ({type})"))

    if (!is.null(gene_modules)) {
        new_names <- ifelse(!is.na(modules), glue("{new_names}<br />Gene module: {modules}"), new_names)
    }

    return(new_names)
}

add_gene_modules <- function(genes, dataset, gene_modules) {
    modules <- tibble(gene = genes) %>%
        left_join(gene_modules, by = join_by(gene)) %>%
        pull(module) %>%
        as.character()

    return(ifelse(!is.na(modules), glue("{genes} ({modules})"), genes))
}

#' Find top-2 values per column of a matrix
#'
#' For each column, return the row indices and values of the two largest
#' finite elements. Columns with fewer than two finite entries get NA for
#' the missing positions. Ties are broken by row order (lowest index first).
#'
#' @param m Numeric matrix (genes x metacells in MCView's EGC convention).
#' @return A list with integer vectors `top1_idx`, `top2_idx` and numeric
#'   vectors `top1_val`, `top2_val`, each of length `ncol(m)`.
#' @noRd
top2_per_col <- function(m) {
    n_cols <- ncol(m)
    top1_idx <- integer(n_cols)
    top1_val <- numeric(n_cols)
    top2_idx <- integer(n_cols)
    top2_val <- numeric(n_cols)
    for (j in seq_len(n_cols)) {
        col <- m[, j]
        finite <- which(is.finite(col))
        if (length(finite) == 0L) {
            top1_idx[j] <- NA_integer_
            top1_val[j] <- NA_real_
            top2_idx[j] <- NA_integer_
            top2_val[j] <- NA_real_
        } else if (length(finite) == 1L) {
            top1_idx[j] <- finite
            top1_val[j] <- col[finite]
            top2_idx[j] <- NA_integer_
            top2_val[j] <- NA_real_
        } else {
            ord <- finite[order(col[finite], decreasing = TRUE)]
            top1_idx[j] <- ord[1]
            top1_val[j] <- col[ord[1]]
            top2_idx[j] <- ord[2]
            top2_val[j] <- col[ord[2]]
        }
    }
    list(
        top1_idx = top1_idx,
        top1_val = top1_val,
        top2_idx = top2_idx,
        top2_val = top2_val
    )
}

#' Find top-2 values per row of a matrix
#'
#' Row-oriented variant of [top2_per_col()].
#' @noRd
top2_per_row <- function(m) {
    top2_per_col(t(m))
}
