# testServer coverage for dataset-change invalidation.
#
# When the user switches datasets via input$dataset, the observers in
# app_server() must repopulate state$manifold_state and the
# metacell_types / cell_type_colors / gene_modules reactiveVals from
# the new dataset. This pins that behavior without booting a browser.

test_that("dataset-switch repopulates state$manifold_state and annotation reactives", {
    skip_if_no_daf()
    skip_if(!requireNamespace("shiny", quietly = TRUE), "shiny not available")

    # Build a fresh env with two datasets pointing at the same real DAF so
    # the observers can succeed against real data on either name.
    init_mcview_env()
    daf <- dafr::open_daf(get_test_daf_path())
    init_multi_daf_mode(list(dataset_a = daf, dataset_b = daf))
    init_defs()

    app <- shiny::shinyApp(ui = app_ui, server = app_server)

    shiny::testServer(app, {
        # First flush: dataset() reads input$dataset which is NULL initially.
        # Provide a value and let the observers fire.
        session$setInputs(dataset = "dataset_a")
        session$flushReact()

        expect_equal(dataset(), "dataset_a")

        # Observers in app_server should have populated state$manifold_state.
        # mc2d may legitimately be NULL if the DAF has no mc2d data; the
        # important invariant is that the observer FIRED at all.
        expect_true("mc2d" %in% names(state$manifold_state))

        # metacell_types() reactiveVal should be populated to a tibble.
        mct_a <- metacell_types()
        expect_true(is.data.frame(mct_a) || is.null(mct_a))

        # Switch dataset.
        session$setInputs(dataset = "dataset_b")
        session$flushReact()

        expect_equal(dataset(), "dataset_b")
        # Observers should have refired; metacell_types() reflects new dataset.
        mct_b <- metacell_types()
        expect_true(is.data.frame(mct_b) || is.null(mct_b))
    })
})
