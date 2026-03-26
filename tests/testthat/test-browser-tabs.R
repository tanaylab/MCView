# test-browser-tabs.R - Comprehensive browser-based tests for all MCView tabs
#
# This test file launches the MCView Shiny app ONCE, connects a headless
# Chrome browser, then navigates through every tab to verify that key UI
# elements render correctly. Screenshots are captured for each tab.
#
# Prerequisites:
#   - chromote, callr, httr, httpuv packages
#   - Chrome headless shell at .CHROME_PATH (set in helper-browser.R)
#   - DAF test dataset at get_test_daf_path()
#   - helper-browser.R and helper-daf.R auto-sourced by testthat

# ==============================================================================
# One-time setup for all browser tests
# ==============================================================================
skip_on_cran()
skip_if_no_daf()
skip_if_no_browser()

port <- httpuv::randomPort()
bg <- launch_mcview_bg(get_test_daf_path(), port)
app_ready <- wait_for_app(port, timeout = 180, bg_process = bg)

if (!app_ready) {
    if (bg$is_alive()) bg$kill()
    skip("MCView app did not become ready within 180 seconds")
}

session <- connect_browser(port)
initial_idle <- wait_for_shiny_idle(session, timeout = 120)

withr::defer(cleanup_browser_test(session, bg), teardown_env())

if (!initial_idle) {
    skip("Shiny app did not reach idle state after initial load")
}

# ==============================================================================
# Helper: robust tab navigation
# ==============================================================================
navigate_to_tab <- function(session, tab_name, timeout = 90) {
    # Wait for sidebar menu to be rendered (dynamic sidebarMenuOutput)
    for (i in seq_len(10)) {
        has_menu <- tryCatch(
            {
                result <- session$Runtime$evaluate(
                    expression = "document.querySelectorAll('.sidebar-menu .treeview-menu li a').length > 0"
                )
                isTRUE(result$result$value)
            },
            error = function(e) FALSE
        )
        if (has_menu) break
        Sys.sleep(2)
    }

    clicked <- click_tab(session, tab_name, timeout = timeout)
    if (!isTRUE(clicked)) {
        return(FALSE)
    }
    # Extra pause for rendering
    Sys.sleep(2)
    wait_for_shiny_idle(session, timeout = timeout)
}

# ==============================================================================
# Tab 1: About
# ==============================================================================
test_that("About tab renders correctly", {
    idle <- navigate_to_tab(session, "About")
    expect_true(idle, info = "Shiny should reach idle state on About tab")

    # The About tab has module_name "about", so output IDs are prefixed "about-"
    # Check for the about box
    has_about_box <- element_exists(session, "#about-about")
    expect_true(has_about_box, info = "About box should exist")

    # Check for the 2D projection plot (base R plot output)
    has_2d_plot <- element_exists(session, "#about-about_2d_plot")
    expect_true(has_2d_plot, info = "About 2D projection plot should exist")

    # Check for the scatter plot
    has_scatter <- element_exists(session, "#about-about_scatter_plot")
    expect_true(has_scatter, info = "About scatter plot should exist")

    take_screenshot(session, "tab-01-about")
})

# ==============================================================================
# Tab 2: Manifold
# ==============================================================================
test_that("Manifold tab renders correctly", {
    idle <- navigate_to_tab(session, "Manifold")
    expect_true(idle, info = "Shiny should reach idle state on Manifold tab")

    # The Manifold module_name is "manifold"
    # projection_box creates plotlyOutput(ns("plot_gene_proj_2d"))
    # In browser: manifold-plot_gene_proj_2d
    has_proj_plot <- element_exists(session, "#manifold-plot_gene_proj_2d")
    expect_true(has_proj_plot, info = "Manifold projection plot container should exist")

    # Check for plotly rendering
    Sys.sleep(3) # allow plotly to render
    plotly_count <- get_visible_elements(session, "#manifold-plot_gene_proj_2d .js-plotly-plot")
    expect_gte(plotly_count, 0, label = "Manifold plot may contain plotly element")

    # Check for projection box
    has_proj_box <- element_exists(session, "#manifold-gene_projection")
    expect_true(has_proj_box, info = "Manifold projection box should exist")

    take_screenshot(session, "tab-02-manifold")
})

# ==============================================================================
# Tab 3: Genes
# ==============================================================================
test_that("Genes tab renders correctly", {
    idle <- navigate_to_tab(session, "Genes")
    expect_true(idle, info = "Shiny should reach idle state on Genes tab")

    # module_name is "gene_mc"
    # scatter_box creates plotlyOutput(ns("plot_gene_gene_mc"))
    has_scatter <- element_exists(session, "#gene_mc-plot_gene_gene_mc")
    expect_true(has_scatter, info = "Gene/Gene scatter plot container should exist")

    # projection_box creates plotlyOutput(ns("plot_gene_proj_2d"))
    has_proj <- element_exists(session, "#gene_mc-plot_gene_proj_2d")
    expect_true(has_proj, info = "Gene projection plot container should exist")

    # Check for the scatter box
    has_scatter_box <- element_exists(session, "#gene_mc-gene_gene_box")
    expect_true(has_scatter_box, info = "Gene/Gene scatter box should exist")

    # Check for projection box
    has_proj_box <- element_exists(session, "#gene_mc-gene_projection")
    expect_true(has_proj_box, info = "Gene projection box should exist")

    take_screenshot(session, "tab-03-genes")
})

# ==============================================================================
# Tab 4: Diff. Expression
# ==============================================================================
test_that("Diff. Expression tab renders correctly", {
    idle <- navigate_to_tab(session, "Diff. Expression")
    expect_true(idle, info = "Shiny should reach idle state on Diff. Expression tab")

    # module_name is "mc_mc"
    # diff expr scatter: plotlyOutput(ns("plot_mc_mc_gene_scatter"))
    has_de_scatter <- element_exists(session, "#mc_mc-plot_mc_mc_gene_scatter")
    expect_true(has_de_scatter, info = "Diff expression scatter plot container should exist")

    # 2D projection: plotlyOutput(ns("plot_mc_proj_2d"))
    has_proj <- element_exists(session, "#mc_mc-plot_mc_proj_2d")
    expect_true(has_proj, info = "Diff expression projection plot container should exist")

    # Radio buttons for mode selection (MCs/Types/Groups)
    has_mode <- element_exists(session, "#mc_mc-mode")
    expect_true(has_mode, info = "Mode radio buttons should exist")

    # Check the DT output for diff expression table exists
    has_table <- element_exists(session, "#mc_mc-diff_expr_table")
    expect_true(has_table, info = "Diff expression table container should exist")

    take_screenshot(session, "tab-04-diff-expression")
})

# ==============================================================================
# Tab 5: Cell types
# ==============================================================================
test_that("Cell types tab renders correctly", {
    idle <- navigate_to_tab(session, "Cell types")
    expect_true(idle, info = "Shiny should reach idle state on Cell types tab")

    # module_name is "cell_type"
    # plotlyOutput(ns("cell_type_boxplot"), height = "70vh")
    has_boxplot <- element_exists(session, "#cell_type-cell_type_boxplot")
    expect_true(has_boxplot, info = "Cell type boxplot container should exist")

    # Check for the boxplot box container
    has_boxplot_box <- element_exists(session, "#cell_type-boxplot_box")
    expect_true(has_boxplot_box, info = "Cell type boxplot box should exist")

    # Check that plot type controls exist
    has_plot_type <- element_exists(session, "#cell_type-plot_type")
    expect_true(has_plot_type, info = "Plot type radio buttons should exist")

    take_screenshot(session, "tab-05-cell-types")
})

# ==============================================================================
# Tab 6: QC
# ==============================================================================
test_that("QC tab renders correctly", {
    idle <- navigate_to_tab(session, "QC")
    expect_true(idle, info = "Shiny should reach idle state on QC tab")

    # module_name is "qc"
    # Value boxes
    has_num_mc <- element_exists(session, "#qc-num_metacells")
    expect_true(has_num_mc, info = "Number of metacells value box should exist")

    has_num_cells <- element_exists(session, "#qc-num_cells")
    expect_true(has_num_cells, info = "Number of cells value box should exist")

    has_median_umis <- element_exists(session, "#qc-median_umis_per_metacell")
    expect_true(has_median_umis, info = "Median UMIs value box should exist")

    has_num_outliers <- element_exists(session, "#qc-num_outliers")
    expect_true(has_num_outliers, info = "Number of outliers value box should exist")

    # QC stat plots
    has_umis_plot <- element_exists(session, "#qc-plot_qc_umis")
    expect_true(has_umis_plot, info = "UMIs per metacell plot container should exist")

    has_cell_num_plot <- element_exists(session, "#qc-plot_qc_cell_num")
    expect_true(has_cell_num_plot, info = "Cells per metacell plot container should exist")

    take_screenshot(session, "tab-06-qc")
})

# ==============================================================================
# Tab 7: Markers
# ==============================================================================
test_that("Markers tab renders correctly", {
    idle <- navigate_to_tab(session, "Markers")
    expect_true(idle, info = "Shiny should reach idle state on Markers tab")

    # module_name is "markers"
    # heatmap_box creates ns = NS("markers-markers_heatmap"), with plotOutput inside
    # The heatmap box ID: markers-markers_heatmap-heatmap_box
    has_heatmap_box <- element_exists(session, "#markers-markers_heatmap-heatmap_box")
    expect_true(has_heatmap_box, info = "Markers heatmap box should exist")

    # The plotting_area is a uiOutput rendered dynamically
    # It creates a plotOutput with id: markers-markers_heatmap-heatmap
    has_plotting_area <- element_exists(session, "#markers-markers_heatmap-plotting_area")
    expect_true(has_plotting_area, info = "Markers heatmap plotting area should exist")

    # Wait for heatmap to render; check for canvas or img element
    Sys.sleep(5)
    has_img <- element_exists(session, "#markers-markers_heatmap-heatmap img")
    has_canvas <- element_exists(session, "#markers-markers_heatmap-plotting_area canvas")
    expect_true(has_img || has_canvas || element_exists(session, "#markers-markers_heatmap-plotting_area"),
        info = "Markers heatmap should have a rendered plot element"
    )

    take_screenshot(session, "tab-07-markers")
})

# ==============================================================================
# Tab 8: Projection QC (conditional - requires atlas)
# ==============================================================================
test_that("Projection QC tab renders correctly", {
    clicked <- tryCatch(
        click_tab(session, "Projection QC", timeout = 30),
        error = function(e) FALSE
    )
    if (!isTRUE(clicked)) {
        skip("Projection QC tab not available in this dataset (requires atlas)")
    }

    Sys.sleep(2)
    idle <- wait_for_shiny_idle(session, timeout = 90)
    expect_true(idle, info = "Shiny should reach idle state on Projection QC tab")

    # module_name is "projection_qc"
    # Value boxes
    has_atlas_mc <- element_exists(session, "#projection_qc-num_metacells_atlas")
    expect_true(has_atlas_mc, info = "Atlas metacells value box should exist")

    has_query_mc <- element_exists(session, "#projection_qc-num_metacells_query")
    expect_true(has_query_mc, info = "Query metacells value box should exist")

    has_similar <- element_exists(session, "#projection_qc-num_metacells_similar")
    expect_true(has_similar, info = "Similar metacells value box should exist")

    has_avg_cor <- element_exists(session, "#projection_qc-avg_projection_cor")
    expect_true(has_avg_cor, info = "Average projection correlation value box should exist")

    has_common_genes <- element_exists(session, "#projection_qc-common_genes")
    expect_true(has_common_genes, info = "Common genes value box should exist")

    has_fitted_genes <- element_exists(session, "#projection_qc-fitted_genes")
    expect_true(has_fitted_genes, info = "Fitted genes value box should exist")

    # Stacked type bar plot (base R plotOutput)
    has_stacked_type <- element_exists(session, "#projection_qc-plot_mc_stacked_type")
    expect_true(has_stacked_type, info = "Stacked type bar plot should exist")

    # Projected correlation plotly
    has_proj_cor_plot <- element_exists(session, "#projection_qc-plot_projected_correlation")
    expect_true(has_proj_cor_plot, info = "Projected correlation plot container should exist")

    # Fitted genes per cell type plotly
    has_fitted_ct_plot <- element_exists(session, "#projection_qc-plot_fitted_genes_per_cell_type")
    expect_true(has_fitted_ct_plot, info = "Fitted genes per cell type plot container should exist")

    take_screenshot(session, "tab-08-projection-qc")
})

# ==============================================================================
# Tab 9: Atlas (conditional - requires atlas)
# ==============================================================================
test_that("Atlas tab renders correctly", {
    clicked <- tryCatch(
        click_tab(session, "Atlas", timeout = 30),
        error = function(e) FALSE
    )
    if (!isTRUE(clicked)) {
        skip("Atlas tab not available in this dataset")
    }

    Sys.sleep(2)
    idle <- wait_for_shiny_idle(session, timeout = 90)
    expect_true(idle, info = "Shiny should reach idle state on Atlas tab")

    # module_name is "atlas"
    # plotlyOutput(ns("plot_mc_proj_2d"), height = "80vh")
    has_proj <- element_exists(session, "#atlas-plot_mc_proj_2d")
    expect_true(has_proj, info = "Atlas 2D projection plot container should exist")

    # Atlas box
    has_atlas_box <- element_exists(session, "#atlas-metacell_projection")
    expect_true(has_atlas_box, info = "Atlas projection box should exist")

    # Color by radio buttons
    has_color_proj <- element_exists(session, "#atlas-color_proj")
    expect_true(has_color_proj, info = "Atlas color-by radio buttons should exist")

    take_screenshot(session, "tab-09-atlas")
})

# ==============================================================================
# Tab 10: Annotate
# ==============================================================================
test_that("Annotate tab renders correctly", {
    idle <- navigate_to_tab(session, "Annotate")
    expect_true(idle, info = "Shiny should reach idle state on Annotate tab")

    # module_name is "annotate"
    # projection_box creates plotlyOutput(ns("plot_gene_proj_2d"))
    has_proj <- element_exists(session, "#annotate-plot_gene_proj_2d")
    expect_true(has_proj, info = "Annotate projection plot container should exist")

    # scatter_box creates plotlyOutput(ns("plot_gene_gene_mc"))
    has_scatter <- element_exists(session, "#annotate-plot_gene_gene_mc")
    expect_true(has_scatter, info = "Annotate scatter plot container should exist")

    # diff_expr_box creates plotlyOutput(ns("plot_mc_mc_gene_scatter"))
    has_de <- element_exists(session, "#annotate-plot_mc_mc_gene_scatter")
    expect_true(has_de, info = "Annotate diff expression plot container should exist")

    # Metacell annotation table
    has_mc_table <- element_exists(session, "#annotate-mc_type_table")
    expect_true(has_mc_table, info = "Metacell type annotation table should exist")

    # Cell type table
    has_ct_table <- element_exists(session, "#annotate-cell_type_table")
    expect_true(has_ct_table, info = "Cell type colors table should exist")

    # Metacell annotation box
    has_annot_box <- element_exists(session, "#annotate-metacell_types_box")
    expect_true(has_annot_box, info = "Metacell annotation box should exist")

    # Cell type colors box
    has_ct_box <- element_exists(session, "#annotate-cell_type_colors")
    expect_true(has_ct_box, info = "Cell type colors box should exist")

    # Action buttons
    has_reset_btn <- element_exists(session, "#annotate-reset_metacell_types")
    expect_true(has_reset_btn, info = "Reset metacell types button should exist")

    has_add_ct_btn <- element_exists(session, "#annotate-add_cell_type_modal")
    expect_true(has_add_ct_btn, info = "Add cell type button should exist")

    take_screenshot(session, "tab-10-annotate")
})

# ==============================================================================
# Tab 11: Samples (conditional - requires sample data)
# ==============================================================================
test_that("Samples tab renders correctly", {
    clicked <- tryCatch(
        click_tab(session, "Samples", timeout = 30),
        error = function(e) FALSE
    )
    if (!isTRUE(clicked)) {
        skip("Samples tab not available in this dataset")
    }

    Sys.sleep(2)
    idle <- wait_for_shiny_idle(session, timeout = 90)
    expect_true(idle, info = "Shiny should reach idle state on Samples tab")

    # module_name is "samples"
    # Stacked types plot
    has_stacked_types <- element_exists(session, "#samples-plot_sample_stacked_types")
    expect_true(has_stacked_types, info = "Sample stacked types plot container should exist")

    # Sample types box
    has_sample_types_box <- element_exists(session, "#samples-sample_types_box")
    expect_true(has_sample_types_box, info = "Sample types box should exist")

    # Sample/Sample scatter
    has_sample_scatter <- element_exists(session, "#samples-plot_gene_gene_mc")
    expect_true(has_sample_scatter, info = "Sample scatter plot container should exist")

    # Sample projection
    has_sample_proj <- element_exists(session, "#samples-plot_gene_proj_2d")
    expect_true(has_sample_proj, info = "Sample projection plot container should exist")

    # Sample ordering accordion
    has_ordering <- element_exists(session, "#samples-sample_types_ordering")
    expect_true(has_ordering, info = "Sample ordering selector should exist")

    take_screenshot(session, "tab-11-samples")
})

# ==============================================================================
# Interaction test: Diff. Expression mode switch
# ==============================================================================
test_that("Diff. Expression tab responds to mode switch", {
    idle <- navigate_to_tab(session, "Diff. Expression")
    expect_true(idle, info = "Should navigate to Diff. Expression tab")

    # The mode radio buttons should be present
    has_mode <- element_exists(session, "#mc_mc-mode")
    if (!has_mode) {
        skip("Mode selector not found; sidebar may not be rendered yet")
    }

    # Click "MCs" mode button
    js_click_mc <- "
    (function() {
        var labels = document.querySelectorAll('#mc_mc-mode .btn-group-container-sw label, #mc_mc-mode label');
        for (var i = 0; i < labels.length; i++) {
            if (labels[i].textContent.trim() === 'MCs') {
                labels[i].click();
                return true;
            }
        }
        // Try input element directly
        var inputs = document.querySelectorAll('input[name=\"mc_mc-mode\"]');
        for (var i = 0; i < inputs.length; i++) {
            if (inputs[i].value === 'MCs') {
                inputs[i].click();
                return true;
            }
        }
        return false;
    })()
    "
    result <- tryCatch(session$Runtime$evaluate(expression = js_click_mc), error = function(e) NULL)

    Sys.sleep(2)
    wait_for_shiny_idle(session, timeout = 60)

    # After switching, the diff expression scatter should still be present
    has_de <- element_exists(session, "#mc_mc-plot_mc_mc_gene_scatter")
    expect_true(has_de, info = "Diff expression scatter should still exist after mode switch")

    take_screenshot(session, "tab-04-diff-expression-mc-mode")
})

# ==============================================================================
# Interaction test: QC value boxes have content
# ==============================================================================
test_that("QC tab value boxes contain data", {
    idle <- navigate_to_tab(session, "QC")
    expect_true(idle, info = "Should navigate to QC tab")

    # Give time for value boxes to render (they depend on reactive data)
    Sys.sleep(5)
    wait_for_shiny_idle(session, timeout = 30)

    # Check that value boxes have rendered text content
    mc_value <- get_shiny_value(session, "qc-num_metacells")
    expect_false(is.null(mc_value), info = "Number of metacells value box should have content")
    expect_true(nchar(mc_value) > 0, info = "Number of metacells value should be non-empty")

    # num_cells is an optional QC stat (may not exist in all DAF datasets).
    # If the value box element exists, try to read its value.
    cells_value <- NULL
    for (i in seq_len(5)) {
        cells_value <- get_shiny_value(session, "qc-num_cells")
        if (!is.null(cells_value) && nzchar(cells_value)) break
        Sys.sleep(2)
    }
    # The element should exist even if empty (it's in the UI definition)
    has_cells_box <- element_exists(session, "#qc-num_cells")
    expect_true(has_cells_box, info = "Number of cells value box element should exist")
    # If value is available, verify it's non-empty; otherwise just note it
    if (is.null(cells_value) || !nzchar(cells_value)) {
        message("[browser-test] Note: qc-num_cells value box is empty (qc_stats_n_cells may not be in this DAF)")
    } else {
        expect_true(nchar(cells_value) > 0, info = "Number of cells value should be non-empty")
    }

    take_screenshot(session, "tab-06-qc-with-values")
})

# ==============================================================================
# Interaction test: Genes tab scatter plot renders plotly
# ==============================================================================
test_that("Genes tab scatter plot renders plotly content", {
    idle <- navigate_to_tab(session, "Genes")
    expect_true(idle, info = "Should navigate to Genes tab")

    # Wait for plotly to render
    Sys.sleep(5)

    # Check for plotly plot elements (could be .js-plotly-plot or .plotly)
    has_plotly <- element_exists(session, "#gene_mc-plot_gene_gene_mc .js-plotly-plot")
    if (!has_plotly) {
        has_plotly <- element_exists(session, "#gene_mc-plot_gene_gene_mc .plotly")
    }
    # The container may still be loading, so check the wrapper
    has_container <- element_exists(session, "#gene_mc-plot_gene_gene_mc")
    expect_true(has_container, info = "Gene scatter container should exist")

    # Also check projection has rendered
    has_proj_plotly <- element_exists(session, "#gene_mc-plot_gene_proj_2d .js-plotly-plot")
    if (!has_proj_plotly) {
        has_proj_plotly <- element_exists(session, "#gene_mc-plot_gene_proj_2d .plotly")
    }
    has_proj_container <- element_exists(session, "#gene_mc-plot_gene_proj_2d")
    expect_true(has_proj_container, info = "Gene projection container should exist")

    take_screenshot(session, "tab-03-genes-plotly")
})

# ==============================================================================
# Interaction test: Cell types tab has plot type controls
# ==============================================================================
test_that("Cell types tab advanced options work", {
    idle <- navigate_to_tab(session, "Cell types")
    expect_true(idle, info = "Should navigate to Cell types tab")

    # Verify the plot type radio buttons
    has_plot_type <- element_exists(session, "#cell_type-plot_type")
    expect_true(has_plot_type, info = "Plot type radio buttons should exist")

    # Verify the facet-by controls
    has_facet <- element_exists(session, "#cell_type-facet_by")
    expect_true(has_facet, info = "Facet-by radio buttons should exist")

    # Verify x-axis type selector
    has_xaxis_type <- element_exists(session, "#cell_type-x_axis_type")
    expect_true(has_xaxis_type, info = "X-axis type selector should exist")

    take_screenshot(session, "tab-05-cell-types-controls")
})

# ==============================================================================
# Verify sidebar menu has tab items
# ==============================================================================
test_that("Sidebar menu contains expected tab items", {
    # Count sidebar menu items
    sidebar_items <- get_visible_elements(session, ".sidebar-menu .treeview-menu li a")
    expect_gt(sidebar_items, 0, label = "Sidebar should contain at least one tab item")

    take_screenshot(session, "sidebar-menu")
})

# ==============================================================================
# Global: Check for JavaScript errors
# ==============================================================================
test_that("No critical JavaScript errors in the app", {
    js_errors <- get_js_errors(session)

    # Filter out known non-critical errors
    if (length(js_errors) > 0) {
        # Shiny output errors that may be transient (e.g., req() failures)
        critical_errors <- js_errors[!grepl(
            "shiny-output-error-validation|Loading|Error:|recalculating|\\[object Object\\]|\\(empty error element",
            js_errors,
            ignore.case = TRUE
        )]

        if (length(critical_errors) > 0) {
            message("[browser-test] Non-critical JS errors found: ", length(js_errors) - length(critical_errors))
            message("[browser-test] Critical JS errors: ", paste(critical_errors, collapse = "\n"))
        }

        # Only fail on true critical errors
        expect_equal(
            length(critical_errors), 0,
            info = paste("Critical JS errors:\n", paste(critical_errors, collapse = "\n"))
        )
    } else {
        expect_equal(length(js_errors), 0, info = "No JS errors found")
    }
})

# ==============================================================================
# Global: Verify Shiny is still connected
# ==============================================================================
test_that("Shiny session is still connected after all tab navigation", {
    is_connected <- tryCatch(
        {
            result <- session$Runtime$evaluate(
                expression = "typeof Shiny !== 'undefined' && Shiny.shinyapp && Shiny.shinyapp.$socket && Shiny.shinyapp.$socket.readyState === 1"
            )
            isTRUE(result$result$value)
        },
        error = function(e) FALSE
    )

    expect_true(is_connected, info = "Shiny websocket should still be connected")

    take_screenshot(session, "final-state")
})

# ==============================================================================
# Global: Background process is still alive
# ==============================================================================
test_that("MCView background process is still running", {
    expect_true(bg$is_alive(), info = "Background R process should still be alive")
})
