# Smoke tests for MCView module server functions using shiny::testServer()
#
# These tests verify that each module's server logic can initialize without
# throwing errors when given a real DAF dataset. They do NOT test UI rendering
# or complex user interactions -- they confirm the reactive graph can be set up.
#
# Each test:
#   1. Opens the OBK DAF and initializes the MCView environment
#   2. Retrieves the data needed for module args (metacell_types, cell_type_colors, etc.)
#   3. Calls shiny::testServer() with the module's server function
#   4. Flushes reactive dependencies
#   5. Passes if no error is thrown

# ==============================================================================
# Setup
# ==============================================================================

test_daf_path <- "/home/obk/data/mcview/metacells_clean"

# Track whether DAF setup has been completed in this R session
.daf_setup_done_mod <- new.env(parent = emptyenv())
.daf_setup_done_mod$ok <- FALSE

skip_if_no_daf <- function() {
    skip_if(!dir.exists(test_daf_path), "OBK DAF not available")
    skip_if(!requireNamespace("dafr", quietly = TRUE), "dafr not installed")
    if (!.daf_setup_done_mod$ok) {
        tryCatch(
            {
                dafr::setup_daf(pkg_check = FALSE, julia_environment = "custom")
                .daf_setup_done_mod$ok <- TRUE
            },
            error = function(e) {
                skip(paste("dafr::setup_daf() failed:", conditionMessage(e)))
            }
        )
    }
}

# Helper: set up the full app state and return data needed by module args.
# Every test that calls this should have already called skip_if_no_daf().
setup_app_state <- function() {
    skip_if_no_daf()

    daf <- dafr::open_daf(test_daf_path)
    init_mcview_env()
    init_single_daf_mode(daf, "test_data")

    # Fetch data the same way app_server.R does
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

    # Add cell type color to metacell_types (mirrors app_server.R logic)
    if (!is.null(mc_types) && !is.null(ct_colors)) {
        if ("mc_col" %in% colnames(mc_types)) {
            mc_types$mc_col <- NULL
        }
        mc_types <- mc_types %>%
            dplyr::left_join(ct_colors %>% dplyr::select(cell_type, mc_col = color), by = "cell_type")
    }

    # Create mock globals matching app_server.R initialization
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

    list(
        mc_types = mc_types,
        ct_colors = ct_colors,
        gene_mods = gene_mods,
        globals = globals
    )
}

# Helper: build the standard args list for a module server.
# `state` is the output of setup_app_state().
build_module_args <- function(state) {
    list(
        dataset = shiny::reactive("test_data"),
        metacell_types = shiny::reactiveVal(state$mc_types),
        cell_type_colors = shiny::reactiveVal(state$ct_colors),
        gene_modules = shiny::reactiveVal(state$gene_mods),
        globals = state$globals
    )
}

# ==============================================================================
# Module tests
# ==============================================================================

# --- mod_about_server ---------------------------------------------------------

test_that("mod_about_server initializes without error", {
    state <- setup_app_state()
    args <- build_module_args(state)

    testServer(mod_about_server, args = args, {
        session$flushReact()
        expect_true(TRUE) # reached here without error
    })
})

# --- mod_qc_server -----------------------------------------------------------

test_that("mod_qc_server initializes without error", {
    state <- setup_app_state()
    args <- build_module_args(state)

    testServer(mod_qc_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_manifold_server ------------------------------------------------------

test_that("mod_manifold_server initializes without error", {
    state <- setup_app_state()
    args <- build_module_args(state)

    testServer(mod_manifold_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_gene_mc_server -------------------------------------------------------

test_that("mod_gene_mc_server initializes without error", {
    state <- setup_app_state()
    args <- build_module_args(state)

    testServer(mod_gene_mc_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_mc_mc_server (Diff. Expression) --------------------------------------

test_that("mod_mc_mc_server initializes without error", {
    state <- setup_app_state()
    args <- build_module_args(state)

    testServer(mod_mc_mc_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_markers_server -------------------------------------------------------

test_that("mod_markers_server initializes without error", {
    state <- setup_app_state()
    args <- build_module_args(state)

    testServer(mod_markers_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_cell_type_server -----------------------------------------------------

test_that("mod_cell_type_server initializes without error", {
    state <- setup_app_state()
    args <- build_module_args(state)

    testServer(mod_cell_type_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_inner_fold_server ----------------------------------------------------

test_that("mod_inner_fold_server initializes without error", {
    state <- setup_app_state()

    # Inner-fold requires gene x metacell inner_fold matrix in the DAF
    daf_obj <- get_dataset_daf("test_data")
    skip_if(
        is.null(daf_obj) || !dafr::has_matrix(daf_obj, "gene", "metacell", "inner_fold"),
        "inner_fold matrix not available in DAF"
    )

    args <- build_module_args(state)

    testServer(mod_inner_fold_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_stdev_fold_server ----------------------------------------------------

test_that("mod_stdev_fold_server initializes without error", {
    state <- setup_app_state()

    # Stdev-fold requires gene x metacell inner_stdev_log matrix in the DAF
    daf_obj <- get_dataset_daf("test_data")
    skip_if(
        is.null(daf_obj) || !dafr::has_matrix(daf_obj, "gene", "metacell", "inner_stdev_log"),
        "inner_stdev_log matrix not available in DAF"
    )

    args <- build_module_args(state)

    testServer(mod_stdev_fold_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_annotate_server ------------------------------------------------------

test_that("mod_annotate_server initializes without error", {
    state <- setup_app_state()
    args <- build_module_args(state)

    testServer(mod_annotate_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_gene_modules_server -------------------------------------------------

test_that("mod_gene_modules_server initializes without error", {
    state <- setup_app_state()

    # Gene modules tab requires gene_modules data
    skip_if(is.null(state$gene_mods), "gene_modules data not available")

    args <- build_module_args(state)

    testServer(mod_gene_modules_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_flow_server ----------------------------------------------------------

test_that("mod_flow_server initializes without error", {
    state <- setup_app_state()

    # Flow tab requires type_flow data
    type_flow <- tryCatch(
        get_mc_data("test_data", "type_flow"),
        error = function(e) NULL
    )
    skip_if(is.null(type_flow), "type_flow data not available for Flow tab")

    args <- build_module_args(state)

    testServer(mod_flow_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_outliers_server ------------------------------------------------------

test_that("mod_outliers_server initializes without error", {
    state <- setup_app_state()

    # Outliers tab needs outliers data
    outliers_metadata <- tryCatch(
        get_mc_data("test_data", "outliers_metadata"),
        error = function(e) NULL
    )
    skip_if(is.null(outliers_metadata), "outliers_metadata not available")

    args <- build_module_args(state)

    testServer(mod_outliers_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_samples_server -------------------------------------------------------

test_that("mod_samples_server initializes without error", {
    state <- setup_app_state()

    # Samples tab requires sample data in the DAF
    has_samples_data <- tryCatch(
        has_samples("test_data"),
        error = function(e) FALSE
    )
    skip_if(!has_samples_data, "Samples data not available in DAF")

    args <- build_module_args(state)

    testServer(mod_samples_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_atlas_server ---------------------------------------------------------

test_that("mod_atlas_server initializes without error", {
    state <- setup_app_state()

    # Atlas tab requires atlas data
    has_atlas_data <- tryCatch(
        has_atlas("test_data"),
        error = function(e) FALSE
    )
    skip_if(!has_atlas_data, "Atlas data not available in DAF")

    args <- build_module_args(state)

    testServer(mod_atlas_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_projection_qc_server ------------------------------------------------

test_that("mod_projection_qc_server initializes without error", {
    state <- setup_app_state()

    # Projection QC requires atlas data
    has_atlas_data <- tryCatch(
        has_atlas("test_data"),
        error = function(e) FALSE
    )
    skip_if(!has_atlas_data, "Atlas data not available for Projection QC tab")

    args <- build_module_args(state)

    testServer(mod_projection_qc_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_proj_fold_server -----------------------------------------------------

test_that("mod_proj_fold_server initializes without error", {
    state <- setup_app_state()

    # Projected-fold requires atlas data
    has_atlas_data <- tryCatch(
        has_atlas("test_data"),
        error = function(e) FALSE
    )
    skip_if(!has_atlas_data, "Atlas data not available for Projected-fold tab")

    args <- build_module_args(state)

    testServer(mod_proj_fold_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# --- mod_query_server ---------------------------------------------------------

test_that("mod_query_server initializes without error", {
    state <- setup_app_state()

    # Query tab requires atlas data
    has_atlas_data <- tryCatch(
        has_atlas("test_data"),
        error = function(e) FALSE
    )
    skip_if(!has_atlas_data, "Atlas data not available for Query tab")

    args <- build_module_args(state)

    testServer(mod_query_server, args = args, {
        session$flushReact()
        expect_true(TRUE)
    })
})

# ==============================================================================
# Deeper reactive tests (for modules where we can easily poke inputs)
# ==============================================================================

test_that("mod_qc_server renders value boxes when reactives resolve", {
    state <- setup_app_state()
    args <- build_module_args(state)

    testServer(mod_qc_server, args = args, {
        # Value boxes are renderValueBox outputs -- they should resolve without error
        # after flushing reactives. We cannot check their HTML in testServer, but
        # we can verify no error is thrown during the render cycle.
        session$flushReact()

        # Note: shinydashboard::renderValueBox outputs (e.g. num_metacells,
        # num_cells) are not exposed in names(output) within testServer because
        # they are not standard shiny render functions. Successful flush without
        # error is the key check here.
        expect_true(TRUE)
    })
})

test_that("mod_mc_mc_server handles mode input", {
    state <- setup_app_state()
    args <- build_module_args(state)

    testServer(mod_mc_mc_server, args = args, {
        # Set the compare mode to "Types" (the default in the UI)
        session$setInputs(mode = "Types")
        session$flushReact()
        expect_true(TRUE)
    })
})

test_that("mod_markers_server initializes markers reactiveVal", {
    state <- setup_app_state()
    args <- build_module_args(state)

    testServer(mod_markers_server, args = args, {
        session$flushReact()
        # The markers reactiveVal is created inside the module; we just verify
        # the module initializes cleanly.
        expect_true(TRUE)
    })
})

test_that("mod_annotate_server initializes annotation state", {
    state <- setup_app_state()
    args <- build_module_args(state)

    testServer(mod_annotate_server, args = args, {
        session$flushReact()
        # The annotation module creates several internal reactiveVals
        # (selected_metacell_types, to_show, last_chosen_cell_type).
        # Successful initialization without error is the key check.
        expect_true(TRUE)
    })
})
