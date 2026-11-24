# The legend drawn beside the heatmaps: the colors of the gene labels and of the cell types, laid out
# (font size, number of columns, truncated labels) to fit the size of the plot.

#' The colors of the gene labels of a heatmap, as get_gene_colors() and the heatmap set them
#'
#' @param mode the mode of the heatmap. Only the gene modules heatmap has module rows.
#' @param highlight_color the color of the highlighted genes
#' @noRd
heatmap_gene_label_colors <- function(mode, highlight_color = "#fba236") {
    colors <- c(
        "lateral" = "blue",
        "noisy" = "red",
        "lateral+noisy" = "purple",
        "disjoined" = "darkgray",
        "module" = "#012901",
        "highlighted" = highlight_color %||% "#fba236"
    )
    if (mode != "Gene modules") {
        colors <- colors[names(colors) != "module"]
    }
    colors
}

#' The layout of the heatmap legend
#'
#' The genes legend is above the cell types legend. The cell types are in as few columns as the height
#' allows. The layout is the one with the largest font size (between \code{min_fontsize} and
#' \code{max_fontsize}) at which both legends fit in \code{width} x \code{height} inches, when at most
#' \code{max_truncated} cell type labels are truncated, each to at least \code{min_visible} of its width.
#' When there is none, it is the one with the smallest font size, and the labels are truncated to the width
#' of their columns.
#'
#' @param ct_widths,gene_widths the widths (in inches) of the labels of the cell types and of the genes at a font size of 1 point
#' @param width,height the size of the plot in inches
#' @param title_widths the widths of the titles of the cell types and genes legends at a font size of 1 point (bold)
#'
#' @return a list with the font sizes and the sizes (in inches) of the elements of the legend, and for each
#' of the genes and cell types legends: the number of rows, the x offset of each column, the column of each
#' label and the width that a label of each column may take (\code{label_width}, the cell types only).
#' @noRd
heatmap_legend_layout <- function(ct_widths, gene_widths, width, height, title_widths = c(0, 0),
                                  max_fontsize = 10, min_fontsize = 6.5, pad = 5 / 96,
                                  max_truncated = max(1, floor(0.1 * length(ct_widths))), min_visible = 0.55) {
    n_ct <- length(ct_widths)
    n_genes <- length(gene_widths)
    W <- max(width - 2 * pad, 0)
    H <- max(height - 2 * pad, 0)

    sizes <- function(fs, spacing) {
        pt <- fs / 72
        list(
            fontsize = fs,
            title_fontsize = fs * 1.1,
            line = spacing * pt, # the height of a row
            title = 1.7 * pt, # the height of a title
            key = 0.85 * pt, # the diameter of the cell type circles
            key_gap = 0.45 * pt, # between a circle and its label
            col_gap = 1.1 * pt, # between columns
            section_gap = 0.9 * pt, # between the genes and the cell types legends
            pad = pad
        )
    }

    # The widths of the columns, when they want `want` and have `total`: the widest ones are cut to the same
    # width, so that the labels lose as little as possible
    fit_widths <- function(want, total) {
        if (sum(want) <= total) {
            return(want)
        }
        w <- sort(want)
        m <- length(w)
        for (k in seq_len(m)) {
            cap <- (total - sum(w[seq_len(k - 1)])) / (m - k + 1)
            if (cap <= w[k]) break
        }
        pmin(want, max(cap, 0))
    }

    offsets <- function(col_w, col_gap) cumsum(c(0, col_w + col_gap))[seq_along(col_w)]

    # The genes legend, row-major: as few rows as fit in the width
    genes_layout <- function(s) {
        if (n_genes == 0) {
            return(list(nrow = 0, x = numeric(0), col = integer(0), height = 0, fits = TRUE))
        }
        w <- gene_widths * s$fontsize
        for (ncol in n_genes:1) {
            col <- (seq_len(n_genes) - 1) %% ncol + 1
            col_w <- as.vector(tapply(w, col, max))
            total <- sum(col_w) + s$col_gap * (ncol - 1)
            if (total <= W || ncol == 1) {
                return(list(
                    nrow = ceiling(n_genes / ncol), x = offsets(col_w, s$col_gap), col = col,
                    height = s$title + ceiling(n_genes / ncol) * s$line + if (n_ct > 0) s$section_gap else 0,
                    fits = total <= W && title_widths[2] * s$title_fontsize <= W
                ))
            }
        }
    }

    # The cell types legend, column-major: as few columns as the height allows
    ct_layout <- function(s, genes_height, nrow = NULL) {
        if (n_ct == 0) {
            return(list(nrow = 0, x = numeric(0), col = integer(0), label_width = numeric(0), fits = TRUE))
        }
        if (is.null(nrow)) {
            max_rows <- floor((H - genes_height - s$title) / s$line + 1e-9)
            if (max_rows < 1) {
                return(NULL)
            }
            nrow <- ceiling(n_ct / ceiling(n_ct / max_rows))
        }
        col <- (seq_len(n_ct) - 1) %/% nrow + 1
        key <- s$key + s$key_gap
        w <- ct_widths * s$fontsize
        col_w <- fit_widths(key + as.vector(tapply(w, col, max)), W - s$col_gap * (max(col) - 1))
        label_width <- pmax(col_w - key, 0)
        visible <- pmin(1, label_width[col] / pmax(w, 1e-9))
        list(
            nrow = nrow, x = offsets(col_w, s$col_gap), col = col, label_width = label_width,
            fits = sum(visible < 1 - 1e-6) <= max_truncated && min(visible) >= min_visible && title_widths[1] * s$title_fontsize <= W
        )
    }

    layout <- function(fs, spacing) {
        s <- sizes(fs, spacing)
        genes <- genes_layout(s)
        ct <- ct_layout(s, genes$height)
        c(s, list(genes = genes, cell_types = ct, fits = genes$fits && isTRUE(ct$fits)))
    }

    for (fs in seq(max_fontsize, min_fontsize, by = -0.25)) {
        for (spacing in c(1.4, 1.2)) {
            l <- layout(fs, spacing)
            if (l$fits) {
                return(l)
            }
        }
    }

    # Nothing fits: the smallest font and spacing, with the labels truncated to their columns. When not even
    # one row fits in the height, a single column, clipped at the bottom.
    l <- layout(min_fontsize, 1.2)
    if (is.null(l$cell_types)) {
        l$cell_types <- ct_layout(l, 0, nrow = n_ct)
    }
    l
}

# The width of each of the labels in inches, at a font size of 1 point, on the current device
label_widths <- function(labels, fontface = "plain") {
    vapply(labels, function(x) {
        grid::convertWidth(grid::grobWidth(grid::textGrob(x, gp = grid::gpar(fontsize = 100, fontface = fontface))), "in", valueOnly = TRUE) / 100
    }, numeric(1), USE.NAMES = FALSE)
}

# The label, cut with an ellipsis to at most max_width inches at the font size
truncate_label <- function(label, max_width, fontsize) {
    if (label_widths(label) * fontsize <= max_width) {
        return(label)
    }
    chars <- strsplit(label, "")[[1]]
    cut <- function(n) paste0(trimws(paste(chars[seq_len(n)], collapse = ""), "right"), "\u2026")
    lo <- 0
    hi <- length(chars) - 1
    while (lo < hi) { # the longest prefix that fits with the ellipsis
        mid <- ceiling((lo + hi) / 2)
        if (label_widths(cut(mid)) * fontsize <= max_width) {
            lo <- mid
        } else {
            hi <- mid - 1
        }
    }
    cut(lo)
}

#' The legend of the heatmaps
#'
#' The genes legend (the colors of the gene labels) above the cell types legend, from the top left, laid out
#' by heatmap_legend_layout() to fit in width x height inches. The text widths are measured on the current
#' device, so call it with the device of the plot open.
#'
#' @param cell_types,cell_type_colors the cell types and their colors (NULL for no cell types legend)
#' @param gene_colors a named vector of the colors of the gene label categories (NULL for no genes legend)
#' @param width,height the size of the plot in inches
#' @noRd
heatmap_legend_grob <- function(cell_types, cell_type_colors, gene_colors, width, height) {
    cell_types <- as.character(cell_types)
    gene_types <- names(gene_colors) %||% character(0)
    titles <- c("Cell types", "Genes")
    ct_widths <- label_widths(cell_types)
    layout <- heatmap_legend_layout(
        ct_widths,
        label_widths(gene_types),
        width = width,
        height = height,
        title_widths = label_widths(titles, fontface = "bold")
    )

    grobs <- list()
    x0 <- layout$pad
    top <- layout$pad # the distance from the top of the next row
    from_top <- function(y) grid::unit(1, "npc") - grid::unit(y, "in")
    title <- function(label) {
        g <- grid::textGrob(label,
            x = grid::unit(x0, "in"), y = from_top(top + layout$title / 2), just = c("left", "centre"),
            gp = grid::gpar(fontsize = layout$title_fontsize, fontface = "bold", col = "black")
        )
        top <<- top + layout$title
        g
    }

    if (length(gene_types) > 0) {
        l <- layout$genes
        grobs <- c(grobs, list(title(titles[2])))
        row <- (seq_along(gene_types) - 1) %/% max(l$col)
        grobs <- c(grobs, list(grid::textGrob(gene_types,
            x = grid::unit(x0 + l$x[l$col], "in"), y = from_top(top + (row + 0.5) * layout$line), just = c("left", "centre"),
            gp = grid::gpar(fontsize = layout$fontsize, col = unname(gene_colors))
        )))
        top <- top + l$nrow * layout$line + if (length(cell_types) > 0) layout$section_gap else 0
    }

    if (length(cell_types) > 0) {
        l <- layout$cell_types
        labels <- cell_types
        truncated <- ct_widths * layout$fontsize > l$label_width[l$col] + 1e-6
        labels[truncated] <- vapply(which(truncated), function(i) truncate_label(labels[i], l$label_width[l$col[i]], layout$fontsize), character(1))
        grobs <- c(grobs, list(title(titles[1])))
        x <- x0 + l$x[l$col]
        y <- from_top(top + ((seq_along(cell_types) - 1) %% l$nrow + 0.5) * layout$line)
        grobs <- c(grobs, list(
            grid::circleGrob(
                x = grid::unit(x + layout$key / 2, "in"), y = y, r = grid::unit(layout$key / 2, "in"),
                gp = grid::gpar(fill = cell_type_colors, col = "gray20", lwd = 0.5)
            ),
            grid::textGrob(labels,
                x = grid::unit(x + layout$key + layout$key_gap, "in"), y = y, just = c("left", "centre"),
                gp = grid::gpar(fontsize = layout$fontsize, col = "black")
            )
        ))
    }

    grid::gTree(children = do.call(grid::gList, grobs), name = "heatmap_legend")
}
