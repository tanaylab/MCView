# The bottom of the legend and the right end of its widest row, in inches
legend_extent <- function(l, n_ct, n_genes, ct_widths) {
    bottom <- l$pad
    if (n_genes > 0) {
        bottom <- bottom + l$title + l$genes$nrow * l$line + if (n_ct > 0) l$section_gap else 0
    }
    right <- 0
    if (n_ct > 0) {
        bottom <- bottom + l$title + l$cell_types$nrow * l$line
        label <- pmin(ct_widths * l$fontsize, l$cell_types$label_width[l$cell_types$col])
        right <- max(l$pad + l$cell_types$x[l$cell_types$col] + l$key + l$key_gap + label)
    }
    c(bottom = bottom, right = right)
}

test_that("the heatmap legend fits in the plot", {
    set.seed(60427)
    gene_widths <- c(0.33, 0.3, 0.6, 0.45, 0.6) / 10 # at 1 point, ~10pt labels of 30-60 px
    for (n in c(1, 5, 20, 41, 80)) {
        ct_widths <- runif(n, 0.2, 2.2) / 10
        for (size in list(c(193, 800), c(150, 576), c(120, 560), c(600, 400))) {
            width <- size[1] / 96
            height <- size[2] / 96
            l <- heatmap_legend_layout(ct_widths, gene_widths, width, height)
            expect_gte(l$fontsize, 6.5)
            expect_lte(l$fontsize, 10)
            expect_equal(sort(unique(l$cell_types$col)), seq_len(max(l$cell_types$col)))
            expect_true(all(table(l$cell_types$col) <= l$cell_types$nrow))
            ext <- legend_extent(l, n, length(gene_widths), ct_widths)
            expect_lte(ext[["bottom"]], height - l$pad + 1e-9)
            expect_lte(ext[["right"]], width - l$pad + 1e-9)
        }
    }
})

test_that("the heatmap legend uses one column and the largest font when there is room", {
    l <- heatmap_legend_layout(rep(0.1, 5), rep(0.05, 5), width = 3, height = 8)
    expect_equal(l$fontsize, 10)
    expect_equal(max(l$cell_types$col), 1)
    expect_true(all(l$cell_types$label_width >= 0.1 * 10))

    # many cell types in a short plot: more columns, and the labels are truncated to them
    l <- heatmap_legend_layout(rep(0.15, 80), rep(0.05, 5), width = 150 / 96, height = 400 / 96)
    expect_gt(max(l$cell_types$col), 1)
    expect_true(all(l$cell_types$label_width < 0.15 * l$fontsize))
})

test_that("the heatmap legend is drawn with or without each of its parts", {
    pdf(NULL)
    withr::defer(dev.off())
    gene_colors <- heatmap_gene_label_colors("Gene modules", "orange")
    expect_true("module" %in% names(gene_colors))
    expect_false("module" %in% names(heatmap_gene_label_colors("Markers")))
    cell_types <- c("Erythroid", "B cells", "Lateral & intermediate mesoderm")
    for (args in list(
        list(cell_types, c("red", "green", "blue"), gene_colors),
        list(NULL, NULL, gene_colors),
        list(cell_types, c("red", "green", "blue"), NULL),
        list(NULL, NULL, NULL)
    )) {
        g <- heatmap_legend_grob(args[[1]], args[[2]], args[[3]], width = 1, height = 1.5)
        expect_s3_class(g, "gTree")
        grid::grid.newpage()
        grid::grid.draw(g)
    }

    expect_equal(truncate_label("T cells", 10, 10), "T cells")
    short <- truncate_label("Lateral & intermediate mesoderm", 1, 10)
    expect_match(short, "^Lateral .*…$")
    expect_lte(label_widths(short) * 10, 1)
})
