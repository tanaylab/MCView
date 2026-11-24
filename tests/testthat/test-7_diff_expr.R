# A toy dataset: 3 cell types of 10 metacells, with genes that differ between them
toy_diff_expr_data <- function() {
    set.seed(60427)
    n_genes <- 300
    metacells <- paste0("m", 1:30)
    types <- rep(c("A", "B", "C"), each = 10)
    means <- matrix(rexp(n_genes * 3, 1 / 5), n_genes, 3, dimnames = list(NULL, c("A", "B", "C")))
    means[1:30, "A"] <- means[1:30, "A"] * 20
    means[31:60, "B"] <- means[31:60, "B"] * 20
    mat <- matrix(rpois(n_genes * 30, means[, types]), n_genes, 30, dimnames = list(paste0("g", 1:n_genes), metacells))
    list(
        mat = mat,
        metacell_types = tibble(metacell = metacells, cell_type = types),
        cell_type_colors = tibble(cell_type = c("A", "B", "C"), color = c("red", "green", "blue"))
    )
}

with_toy_mc_data <- function(env, mc_data_env, envir = parent.frame()) {
    prev <- mget(c("mc_data", "egc_epsilon"), envir = globalenv(), ifnotfound = list(NULL, NULL))
    withr::defer(purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv())), envir = envir)
    assign("mc_data", list(toy = list2env(mc_data_env)), envir = globalenv()) # an environment per dataset, as in the app
    assign("egc_epsilon", 1e-5, envir = globalenv())
}

test_that("the cached diff. expression follows every input of the comparison", {
    toy <- toy_diff_expr_data()
    with_toy_mc_data(environment(), list(mc_mat = toy$mat, mc_sum = colSums(toy$mat)))

    metacell_types <- reactiveVal(toy$metacell_types)
    cell_type_colors <- reactiveVal(toy$cell_type_colors)
    globals <- reactiveValues(clipboard = character(0))
    groupA <- reactiveVal()
    groupB <- reactiveVal()

    server <- function(id) {
        moduleServer(id, function(input, output, session) {
            scatter_df <- mc_mc_gene_scatter_df_reactive(reactive("toy"), input, output, session, metacell_types, cell_type_colors, globals, groupA, groupB)
        })
    }

    # the comparison computed directly from the current inputs
    expected <- function(input) {
        metacell_filter <- if (isTRUE(input$filter_by_clipboard) && length(globals$clipboard) > 0) globals$clipboard
        if (input$mode == "MCs") {
            return(calc_mc_mc_gene_df("toy", input$metacell1, input$metacell2))
        }
        if (input$mode == "Types") {
            types_df <- metacell_types()
            if (!is.null(metacell_filter)) {
                types_df <- types_df %>% filter(metacell %in% metacell_filter)
            }
            return(calc_ct_ct_gene_df("toy", input$metacell1, input$metacell2, types_df))
        }
        types_df <- bind_rows(tibble(metacell = groupA(), cell_type = "Group A"), tibble(metacell = groupB(), cell_type = "Group B"))
        if (!is.null(metacell_filter)) {
            types_df <- types_df %>% filter(metacell %in% metacell_filter)
        }
        calc_ct_ct_gene_df("toy", "Group A", "Group B", types_df)
    }

    testServer(server, {
        check <- function() {
            session$flushReact()
            expect_equal(scatter_df(), expected(input))
        }
        session$setInputs(mode = "Types", metacell1 = "A", metacell2 = "B", filter_by_clipboard = FALSE)
        check()
        ab <- scatter_df()
        session$setInputs(metacell2 = "C")
        check()
        session$setInputs(metacell1 = "B")
        check()
        session$setInputs(metacell1 = "A", metacell2 = "B")
        check()
        expect_equal(scatter_df(), ab)

        # an annotation edit
        metacell_types(metacell_types() %>% mutate(cell_type = ifelse(metacell %in% c("m1", "m2", "m3"), "B", cell_type)))
        check()
        expect_false(isTRUE(all.equal(scatter_df(), ab)))
        metacell_types(toy$metacell_types)
        check()

        # the clipboard filter
        session$setInputs(filter_by_clipboard = TRUE)
        check()
        globals$clipboard <- c("m1", "m2", "m4", "m11", "m12", "m13")
        check()
        globals$clipboard <- c("m1", "m5", "m11", "m19")
        check()
        session$setInputs(filter_by_clipboard = FALSE)
        check()
        expect_equal(scatter_df(), ab)

        # metacells
        session$setInputs(mode = "MCs", metacell1 = "m1", metacell2 = "m15")
        check()
        session$setInputs(metacell2 = "m25")
        check()
        session$setInputs(metacell1 = "m2")
        check()

        # groups
        groupA(c("m1", "m2", "m11"))
        groupB(c("m21", "m22"))
        session$setInputs(mode = "Groups")
        check()
        groupA(c("m1", "m2"))
        check()
        groupB(c("m21", "m23", "m24"))
        check()
        globals$clipboard <- c("m1", "m21", "m24")
        session$setInputs(filter_by_clipboard = TRUE)
        check()
        session$setInputs(filter_by_clipboard = FALSE)
        check()

        session$setInputs(mode = "Types", metacell1 = "A", metacell2 = "B")
        check()
        expect_equal(scatter_df(), ab)
    })
})

test_that("metacell pickers mark the dissimilar metacells by metacell, not by position", {
    # metadata in another order, without m4
    md <- tibble(metacell = c("m3", "m1", "m2"), similar = c("dissimilar", "similar", "similar"))
    choices <- metacell_choices_with_similarity(c("m1", "m2", "m3", "m4"), md)
    expect_identical(unname(choices), c("m1", "m2", "m3", "m4"))
    expect_identical(names(choices), c("m1", "m2", "m3 (dissimilar)", "m4"))
    expect_identical(metacell_choices_with_similarity(c("m1", "m2"), tibble(metacell = c("m1", "m2"), x = 1)), c("m1", "m2"))
    expect_identical(metacell_choices_with_similarity(c("m1", "m2"), NULL), c("m1", "m2"))
})

test_that("the projection's similarity is added to the metadata by metacell", {
    proj <- tibble(metacell = c("m1", "m2", "m3", "m4"), cell_type = "A", similar = c(TRUE, FALSE, TRUE, FALSE))

    expect_identical(add_similar_to_metadata(NULL, proj), tibble(metacell = proj$metacell, similar = proj$similar))

    # another order, a metacell the projection doesn't have, numeric ids
    prev <- tibble(metacell = c(4, 2, 1, 5), x = 1:4)
    res <- add_similar_to_metadata(prev, proj %>% mutate(metacell = sub("m", "", metacell)))
    expect_identical(res$metacell, c("4", "2", "1", "5"))
    expect_identical(res$x, 1:4)
    expect_identical(res$similar, c(FALSE, FALSE, TRUE, NA))

    # a 'similar' field from a previous import (as strings), the same: no warning
    prev <- tibble(metacell = c("m3", "m1", "m2"), similar = c("similar", "similar", "dissimilar"))
    expect_no_warning(res <- add_similar_to_metadata(prev, proj))
    expect_identical(res$similar, c(TRUE, TRUE, FALSE))
    # a different one is overwritten, with a warning
    prev$similar <- c(FALSE, TRUE, FALSE)
    expect_warning(res <- add_similar_to_metadata(prev, proj), "similar")
    expect_identical(res$similar, c(TRUE, TRUE, FALSE))
})

test_that("'Show only fitted genes' of a group uses the cell types of the group's metacells", {
    gene_metadata <- tibble(
        gene = c("g1", "g2", "g3", "g4"),
        cell_type = c("A", "B", "C", "C"),
        fitted_gene = c(TRUE, TRUE, TRUE, FALSE)
    )
    with_toy_mc_data(environment(), list(gene_metadata = gene_metadata))
    df <- tibble(gene = paste0("g", 1:5), diff = 0)
    metacell_types <- function() tibble(metacell = c("m1", "m2", "m3"), cell_type = c("A", "B", "C"))
    # input$metacell1 is still the metacell chosen in the MC mode
    input <- list(show_only_fitted = TRUE, metacell1 = "m1")

    expect_identical(diff_expr_plot_df(df, input, function() "toy", "Group", metacell_types, group = function() c("m2", "m3"))$gene, c("g2", "g3"))
    expect_identical(diff_expr_plot_df(df, input, function() "toy", "MC", metacell_types)$gene, "g1")
    expect_identical(diff_expr_plot_df(df, list(show_only_fitted = TRUE, metacell1 = "C"), function() "toy", "Type", metacell_types)$gene, "g3")
    expect_identical(diff_expr_plot_df(df, list(show_only_fitted = FALSE), function() "toy", "Group", metacell_types)$gene, df$gene)
})
