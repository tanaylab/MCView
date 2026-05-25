# Regression test: cell_type_gene_boxplot() must use a supplied egc_gene
# verbatim. Line 14 was `egc_gene <- egc_gene %||% get_gene_egc(...) + eps`,
# which parses as `(egc_gene %||% get_gene_egc(...)) + eps` because %||% binds
# tighter than +. Callers on the gene-module path (mod_cell_type.R) already add
# egc_epsilon before passing egc_gene, so it got added a second time.

test_that("cell_type_gene_boxplot uses supplied egc_gene without re-adding epsilon", {
    prev_breaks <- mcv_get("expr_breaks")
    prev_eps <- mcv_get("egc_epsilon")
    mcv_set("expr_breaks", c(1e-5, 2e-5, 4e-5, 1e-4, 2e-4, 4e-4, 1e-3, 2e-3, 4e-3, 1e-2, 2e-2, 4e-2, 1e-1, 2e-1, 4e-1, 1))
    mcv_set("egc_epsilon", 1e-5)
    on.exit({
        mcv_set("expr_breaks", prev_breaks)
        mcv_set("egc_epsilon", prev_eps)
    }, add = TRUE)

    mts <- tibble::tibble(
        metacell = c("m1", "m2", "m3", "m4"),
        cell_type = c("A", "A", "B", "B")
    )
    ctc <- tibble::tibble(
        cell_type = c("A", "B"),
        color = c("#ff0000", "#00ff00")
    )
    egc_gene <- c(m1 = 0.10, m2 = 0.20, m3 = 0.30, m4 = 0.40)

    p <- cell_type_gene_boxplot(
        gene = "MyGene",
        dataset = "dummy",
        metacell_types = mts,
        cell_type_colors = ctc,
        egc_gene = egc_gene
    )

    expect_s3_class(p, "ggplot")
    # The plot data must carry the exact egc_gene values, not egc_gene + eps.
    vals <- p$data[["MyGene"]]
    names(vals) <- p$data$metacell
    expect_equal(vals[names(egc_gene)], egc_gene, tolerance = 1e-12)
})
