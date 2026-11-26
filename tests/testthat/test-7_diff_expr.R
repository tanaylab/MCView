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

# mc_data of a dataset named "toy" (an environment per dataset, as in the app), until the calling test ends
with_toy_mc_data <- function(mc_data_env, envir = parent.frame()) {
    prev <- mget(c("mc_data", "egc_epsilon"), envir = globalenv(), ifnotfound = list(NULL, NULL))
    withr::defer(purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv())), envir = envir)
    assign("mc_data", list(toy = list2env(mc_data_env)), envir = globalenv())
    assign("egc_epsilon", 1e-5, envir = globalenv())
}

test_that("the cached diff. expression follows every input of the comparison", {
    toy <- toy_diff_expr_data()
    with_toy_mc_data(list(mc_mat = toy$mat, mc_sum = colSums(toy$mat)))

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

        # metacells that the dataset doesn't have (e.g. pasted from another dataset's clipboard) are left out
        groupA(c("m1", "x99", "m2"))
        session$flushReact()
        expect_equal(scatter_df(), calc_ct_ct_gene_df("toy", "Group A", "Group B", tibble(metacell = c("m1", "m2", groupB()), cell_type = rep(c("Group A", "Group B"), c(2, 3)))))
        # an empty group: a message instead of a blank plot
        groupA(NULL)
        session$flushReact()
        expect_error(scatter_df(), "Please add metacells to group A and to group B", class = "validation")
        groupA("x99")
        session$flushReact()
        expect_error(scatter_df(), "Please add metacells to group A and to group B", class = "validation")
        groupA(c("m1", "m2"))
        globals$clipboard <- c("m1", "m2")
        session$setInputs(filter_by_clipboard = TRUE)
        expect_error(scatter_df(), "group A or group B has no metacells in the clipboard", class = "validation")
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

test_that("the metacell select colors each metacell by a class per color, and sends labels only when they differ", {
    config_of <- function(tag) {
        script <- htmltools::tagQuery(tag)$find("script")$selectedTags()[[1]]
        jsonlite::fromJSON(as.character(script$children[[1]]))
    }
    tag <- metacell_select_input("mod-metacell1", "Metacell A", c("m1", "m2", "m3"), c("red", "#0000FF", "red"), selected = "m2")
    html <- as.character(tag)
    expect_match(html, ".mod-metacell1-dropbox .c1 {color: #FF0000;}", fixed = TRUE)
    expect_match(html, ".mod-metacell1-dropbox .c2 {color: #0000FF;}", fixed = TRUE)
    data <- config_of(tag)
    options <- data$options
    expect_identical(options$choices$value, c("m1", "m2", "m3"))
    expect_identical(options$choices$classNames, c("c1", "c2", "c1"))
    expect_null(options$choices$label)
    expect_identical(data$config$labelKey, "value")
    expect_identical(data$config$selectedValue, "m2")
    expect_identical(data$config$additionalDropboxClasses, "mod-metacell1-dropbox")

    # the "(dissimilar)" labels
    choices <- c(m1 = "m1", `m2 (dissimilar)` = "m2")
    options <- config_of(metacell_select_input("x", "A", choices, c("red", "blue"), selected = "m1"))$options
    expect_identical(options$choices$label, c("m1", "m2 (dissimilar)"))
    expect_identical(options$choices$value, c("m1", "m2"))
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
    with_toy_mc_data(list(gene_metadata = gene_metadata))
    df <- tibble(gene = paste0("g", 1:5), diff = 0)
    metacell_types <- function() tibble(metacell = c("m1", "m2", "m3"), cell_type = c("A", "B", "C"))
    # input$metacell1 is still the metacell chosen in the MC mode
    input <- list(show_only_fitted = TRUE, metacell1 = "m1")

    expect_identical(diff_expr_plot_df(df, input, function() "toy", "Group", metacell_types, group = function() c("m2", "m3"))$gene, c("g2", "g3"))
    expect_identical(diff_expr_plot_df(df, input, function() "toy", "MC", metacell_types)$gene, "g1")
    expect_identical(diff_expr_plot_df(df, list(show_only_fitted = TRUE, metacell1 = "C"), function() "toy", "Type", metacell_types)$gene, "g3")
    expect_identical(diff_expr_plot_df(df, list(show_only_fitted = FALSE), function() "toy", "Group", metacell_types)$gene, df$gene)
})

# The implementations before the samples were weighted without copying the matrix (reference)
old_get_samples_mat <- function(cell_types, metacell_types, dataset, samples = NULL) {
    mc_mat <- get_mc_data(dataset, "mc_mat")
    mc_types <- metacell_types %>%
        filter(cell_type %in% cell_types) %>%
        select(metacell, cell_type) %>%
        deframe()
    mc_mat <- mc_mat[, names(mc_types), drop = FALSE]
    samp_mc_frac <- get_samp_mc_frac(dataset)[samples %||% TRUE, names(mc_types), drop = FALSE]
    samp_sums <- rowSums(samp_mc_frac)
    samp_mc_frac <- samp_mc_frac[samp_sums > 0, , drop = FALSE]
    samp_mc_frac <- samp_mc_frac / rowSums(samp_mc_frac)
    mc_mat %*% t(samp_mc_frac)
}

old_get_samples_egc <- function(cell_types, metacell_types, dataset, samples = NULL) {
    mc_types <- metacell_types %>%
        filter(cell_type %in% cell_types) %>%
        select(metacell, cell_type) %>%
        deframe()
    mc_mat <- get_mc_data(dataset, "mc_mat")[, names(mc_types), drop = FALSE]
    mc_sum <- get_mc_sum(dataset)[names(mc_types)]
    samp_mc_frac <- get_samp_mc_frac(dataset)[samples %||% TRUE, names(mc_types), drop = FALSE]
    samp_mc_frac <- samp_mc_frac / rowSums(samp_mc_frac)
    mc_mat %*% (t(samp_mc_frac) / mc_sum)
}

# Cells of 4 samples in the toy metacells: S4 has cells only in metacells of type C
toy_cell_metadata <- function(metacell_types, without = character(0)) {
    set.seed(19)
    purrr::map_dfr(seq_len(nrow(metacell_types)), function(i) {
        samples <- if (metacell_types$cell_type[i] == "C") c("S1", "S2", "S3", "S4") else c("S1", "S2", "S3")
        n <- sample(5:20, 1)
        tibble(metacell = metacell_types$metacell[i], samp_id = sample(samples, n, replace = TRUE))
    }) %>%
        filter(!(metacell %in% without)) %>%
        mutate(cell_id = paste0("c", row_number()))
}

test_that("the samples' expression is computed without copying the matrix, as before", {
    toy <- toy_diff_expr_data()
    with_toy_mc_data(list(mc_mat = toy$mat, mc_sum = colSums(toy$mat) * 1.1, cell_metadata = toy_cell_metadata(toy$metacell_types)))
    mt <- toy$metacell_types

    for (cell_types in list(c("A", "B"), "B", c("A", "B", "C"))) {
        for (samples in list(NULL, c("S2", "S1"), "S3")) {
            new_mat <- get_samples_mat(cell_types, mt, "toy", samples = samples)
            expect_equal(new_mat, old_get_samples_mat(cell_types, mt, "toy", samples = samples))
            if (!("S4" %in% colnames(old_get_samples_mat(cell_types, mt, "toy", samples = NULL))) && is.null(samples)) {
                # a sample without cells of the cell types: left out (the old egc had NaN for it)
                expect_false("S4" %in% colnames(get_samples_egc(cell_types, mt, "toy", samples = samples)))
                expect_true(all(is.nan(old_get_samples_egc(cell_types, mt, "toy", samples = samples)[, "S4"])))
                expect_equal(get_samples_egc(cell_types, mt, "toy", samples = samples), old_get_samples_egc(cell_types, mt, "toy", samples = samples)[, c("S1", "S2", "S3")])
            } else {
                expect_equal(get_samples_egc(cell_types, mt, "toy", samples = samples), old_get_samples_egc(cell_types, mt, "toy", samples = samples))
            }
        }
    }
    # metacells in another order in the annotation than in the matrix
    shuffled <- mt[sample(nrow(mt)), ]
    expect_equal(get_samples_mat(c("A", "C"), shuffled, "toy"), old_get_samples_mat(c("A", "C"), mt, "toy"))
    expect_equal(get_samples_egc(c("A", "C"), shuffled, "toy"), old_get_samples_egc(c("A", "C"), mt, "toy"))
})

test_that("samples are compared when a metacell has no cells, and a sample without the cell types is reported", {
    toy <- toy_diff_expr_data()
    mt <- toy$metacell_types
    with_toy_mc_data(list(mc_mat = toy$mat, mc_sum = colSums(toy$mat), cell_metadata = toy_cell_metadata(mt, without = "m3")))
    expect_false("m3" %in% colnames(get_samp_mc_frac("toy")))

    # "subscript out of bounds" before
    expect_error(old_get_samples_mat("A", mt, "toy"))
    without_m3 <- mt %>% filter(metacell != "m3")
    expect_equal(get_samples_mat("A", mt, "toy"), old_get_samples_mat("A", without_m3, "toy"))
    expect_equal(get_samples_egc("A", mt, "toy"), old_get_samples_egc("A", without_m3, "toy")[, c("S1", "S2", "S3")])
    weights <- get_samples_mc_weights("A", mt, "toy")
    expect_identical(unname(weights["m3", ]), c(0, 0, 0))
    expect_equal(unname(colSums(weights)), c(1, 1, 1))

    df <- calc_samp_samp_gene_df("toy", "S1", "S2", mt, cell_types = c("A", "B"))
    expect_setequal(df$gene, rownames(toy$mat))
    expect_equal(df$S1[match(rownames(toy$mat), df$gene)], unname(old_get_samples_egc(c("A", "B"), without_m3, "toy")[, "S1"]) + 1e-5)

    # S4 has no metacells of type A: it is left out of the weights, which the Samples tab reports
    expect_identical(colnames(get_samples_mc_weights("A", mt, "toy", samples = c("S1", "S4"))), "S1")
    expect_identical(colnames(get_samples_mc_weights(c("A", "C"), mt, "toy", samples = c("S1", "S4"))), c("S1", "S4"))
})

test_that("a gene's expression in the samples is matched by metacell", {
    toy <- toy_diff_expr_data()
    mt <- toy$metacell_types
    with_toy_mc_data(list(mc_mat = toy$mat, mc_sum = colSums(toy$mat), cell_metadata = toy_cell_metadata(mt, without = "m3")))
    frac <- get_samp_mc_frac("toy")
    egc <- toy$mat["g1", ] / colSums(toy$mat)

    mcs <- colnames(frac)
    expect_equal(get_samples_gene_egc("g1", "toy"), (frac %*% egc[mcs])[, 1])
    # the chosen metacells, one without cells
    chosen <- c("m10", "m2", "m3", "m25")
    f <- frac[, c("m10", "m2", "m25")]
    f <- f / rowSums(f)
    expect_equal(get_samples_gene_egc("g1", "toy", chosen), (f %*% egc[c("m10", "m2", "m25")])[, 1])
})

test_that("the update functions write into the project of an old bundle (a 'project' subdirectory)", {
    bundle <- withr::local_tempdir()
    project <- file.path(bundle, "project")
    dataset_dir <- file.path(project, "cache", "ds")
    dir.create(file.path(project, "config"), recursive = TRUE)
    dir.create(dataset_dir, recursive = TRUE)
    writeLines("title: test", file.path(project, "config", "config.yaml"))
    for (f in c("mc_mat.qs2", "mc_sum.qs2", "mc2d.qs2")) {
        file.create(file.path(dataset_dir, f))
    }
    data.table::fwrite(tibble(metacell = c("m1", "m2"), cell_type = c("A", "B")), file.path(dataset_dir, "metacell_types.tsv"), sep = "\t")
    data.table::fwrite(tibble(cell_type = c("A", "B"), color = c("red", "blue")), file.path(dataset_dir, "cell_type_colors.tsv"), sep = "\t")
    expect_identical(as.character(init_project_dir(bundle)), project)

    modules_file <- file.path(bundle, "modules.csv")
    data.table::fwrite(tibble(gene = c("g1", "g2"), module = c("1", "2")), modules_file)
    suppressMessages(update_gene_modules(bundle, "ds", modules_file))
    expect_identical(data.table::fread(file.path(dataset_dir, "gene_modules.tsv"))$gene, c("g1", "g2"))

    types_file <- file.path(bundle, "types.csv")
    data.table::fwrite(tibble(metacell = c("m1", "m2"), cell_type = c("B", "C"), color = c("blue", "green")), types_file)
    suppressMessages(update_metacell_types(bundle, "ds", types_file))
    expect_identical(data.table::fread(file.path(dataset_dir, "metacell_types.tsv"))$cell_type, c("B", "C"))
    expect_identical(data.table::fread(file.path(dataset_dir, "cell_type_colors.tsv"))$color, c("blue", "green"))

    expect_false(dir.exists(file.path(bundle, "cache")))
})
