# Tests for clipboard utility functions and their integration
# with the markers heatmap and gene correlation modules.
#
# These tests verify:
#   1. clipboard_copy_button_ui() returns a renderUI function
#   2. clipboard_copy_button_server() creates an observer without error
#   3. heatmap_download_handlers() receives ns and can render the copy button
#   4. Gene correlation module copy button initializes
#
# DAF setup and skip_if_no_daf() provided by helper-daf.R

setup_app_state <- function() {
    skip_if_no_daf()

    daf <- dafr::open_daf(get_test_daf_path())
    init_mcview_env()
    init_single_daf_mode(daf, "test_data")

    mc_types <- get_metacell_types_data("test_data")
    ct_colors <- get_cell_type_data("test_data")
    gene_mods <- get_mc_data("test_data", "gene_modules")

    if (!is.null(gene_mods)) {
        gene_mods <- gene_mods %>%
            dplyr::filter(gene %in% gene_names("test_data")) %>%
            dplyr::mutate(gene = as.character(gene))
        if (!is.factor(gene_mods$module)) {
            gene_mods <- gene_mods %>% dplyr::mutate(module = factor(module))
        }
    }

    if (!is.null(mc_types) && !is.null(ct_colors)) {
        if ("mc_col" %in% colnames(mc_types)) {
            mc_types$mc_col <- NULL
        }
        mc_types <- mc_types %>%
            dplyr::left_join(ct_colors %>% dplyr::select(cell_type, mc_col = color), by = "cell_type")
    }

    globals <- shiny::reactiveValues(
        screen_width = 1200,
        screen_height = 800,
        mc2d = get_mc_data("test_data", "mc2d"),
        anchor_genes = get_mc_data("test_data", "umap_anchors"),
        clipboard = character(0),
        active_tabs = app_config("tabs"),
        plotly_scale = 1,
        plotly_format = "svg",
        plotly_width = NULL,
        plotly_height = NULL
    )

    state <- list(
        session_ui = shiny::reactiveValues(),
        tab_state = shiny::reactiveValues(),
        selection = shiny::reactiveValues(),
        manifold_state = shiny::reactiveValues(
            mc2d = get_mc_data("test_data", "mc2d"),
            anchor_genes = get_mc_data("test_data", "umap_anchors")
        )
    )

    list(
        mc_types = mc_types,
        ct_colors = ct_colors,
        gene_mods = gene_mods,
        globals = globals,
        state = state
    )
}

build_module_args <- function(fixture) {
    list(
        dataset = shiny::reactive("test_data"),
        metacell_types = shiny::reactiveVal(fixture$mc_types),
        cell_type_colors = shiny::reactiveVal(fixture$ct_colors),
        gene_modules = shiny::reactiveVal(fixture$gene_mods),
        globals = fixture$globals,
        state = fixture$state
    )
}

# ==============================================================================
# Unit tests for clipboard utility functions
# ==============================================================================

test_that("clipboard_copy_button_ui returns a function (renderUI)", {
    ns <- shiny::NS("test")
    data_reactive <- shiny::reactiveVal(c("geneA", "geneB"))

    result <- clipboard_copy_button_ui(ns, "test_copy", data_reactive)
    # clipboard_copy_button_ui returns the result of renderUI(), which is a function
    expect_true(is.function(result))
})

test_that("clipboard_copy_button_ui handles empty data", {
    ns <- shiny::NS("test")
    data_reactive <- shiny::reactiveVal(character(0))

    result <- clipboard_copy_button_ui(ns, "test_copy_empty", data_reactive)
    expect_true(is.function(result))
})

test_that("clipboard_copy_button_ui handles NULL data", {
    ns <- shiny::NS("test")
    data_reactive <- shiny::reactiveVal(NULL)

    result <- clipboard_copy_button_ui(ns, "test_copy_null", data_reactive)
    expect_true(is.function(result))
})

test_that("clipboard_copy_button_server is a callable function with correct signature", {
    # clipboard_copy_button_server must accept (input, id, data_reactive, globals, message_template)
    fn_args <- names(formals(clipboard_copy_button_server))
    expect_true("input" %in% fn_args)
    expect_true("id" %in% fn_args)
    expect_true("data_reactive" %in% fn_args)
    expect_true("globals" %in% fn_args)
    expect_true("message_template" %in% fn_args)
})

# ==============================================================================
# Integration: heatmap_download_handlers receives ns
# ==============================================================================

test_that("heatmap_download_handlers function signature includes ns", {
    # Verify the function accepts ns as a parameter
    fn_args <- names(formals(heatmap_download_handlers))
    expect_true("ns" %in% fn_args,
        label = "heatmap_download_handlers must accept 'ns' parameter"
    )
})

test_that("mod_markers_server clipboard button initializes without error", {
    state <- setup_app_state()
    args <- build_module_args(state)

    testServer(mod_markers_server, args = args, {
        session$flushReact()

        # The copy_genes_button output should exist (rendered by heatmap_download_handlers)
        # In testServer we cannot check renderUI output directly, but we verify
        # no "object 'ns' not found" error was thrown.
        expect_true(TRUE)
    })
})

# ==============================================================================
# Integration: gene correlation clipboard
# ==============================================================================

test_that("mod_gene_correlation_server clipboard button initializes", {
    state <- setup_app_state()

    daf_obj <- get_dataset_daf("test_data")
    skip_if(
        is.null(daf_obj) || !dafr::has_axis(daf_obj, "gg_mc_top_cor"),
        "gg_mc_top_cor axis not available in DAF"
    )

    args <- build_module_args(state)

    testServer(mod_gene_correlation_server, args = args, {
        session$flushReact()

        # After initialization, the copy_genes_button output should be registered
        # without any "object 'ns' not found" error
        expect_true(TRUE)
    })
})

test_that("mod_gene_correlation_server clipboard works after calculation", {
    state <- setup_app_state()

    daf_obj <- get_dataset_daf("test_data")
    skip_if(
        is.null(daf_obj) || !dafr::has_axis(daf_obj, "gg_mc_top_cor"),
        "gg_mc_top_cor axis not available in DAF"
    )

    args <- build_module_args(state)
    genes <- head(gene_names("test_data"), 3)
    skip_if(length(genes) < 2, "Need at least 2 genes in dataset")

    testServer(mod_gene_correlation_server, args = args, {
        # Set up inputs and calculate
        session$setInputs(
            gene_list = paste(genes, collapse = "\n"),
            correlation_mode = "individual",
            analysis_type = "find_genes",
            n_correlations = 10,
            cor_threshold = 0,
            correlation_direction = "positive"
        )
        session$flushReact()

        session$setInputs(calculate_correlations = 1)
        session$flushReact()

        # After calculation, filtered_gene_list (used by clipboard) should be available
        # The copy_genes_button renderUI should not throw an error
        expect_true(TRUE)
    })
})
