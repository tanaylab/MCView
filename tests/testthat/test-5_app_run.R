library(golem)

test_that("app ui", {
    init_config(project = project_dir)
    load_all_data(cache_dir = project_cache_dir(project_dir))
    init_defs()
    ui <- app_ui()
    expect_shinytaglist(ui)
})

test_that("default genes are the two with the largest range of expression", {
    init_config(project = project_dir)
    load_all_data(cache_dir = project_cache_dir(project_dir))
    config$selected_gene1 <<- NULL
    init_selected_genes()

    dataset <- dataset_ls(project_dir)[1]
    egc <- get_mc_egc(dataset)
    egc <- egc[!(rownames(egc) %in% c(get_mc_data(dataset, "lateral_genes"), get_mc_data(dataset, "noisy_genes"))), ]
    if (has_atlas(dataset)) {
        egc <- egc[rownames(egc) %in% gene_names(dataset, atlas = TRUE), ]
    }
    expr_range <- matrixStats::rowMaxs(egc) - matrixStats::rowMins(egc)
    names(expr_range) <- rownames(egc)
    expect_equal(c(default_gene1, default_gene2), names(head(sort(expr_range, decreasing = TRUE), n = 2)))
})

test_that("app server", {
    init_config(project = project_dir)
    load_all_data(cache_dir = project_cache_dir(project_dir))
    init_defs()
    server <- app_server
    expect_is(server, "function")
})


test_that(
    "app launches",
    {
        skip_on_ci()
        skip_on_cran()
        skip_on_travis()
        skip_on_appveyor()
        x <- processx::process$new(
            "R",
            c(
                "-e",
                glue("pkgload::load_all(here::here());run_app(project = '{project_dir}')")
            )
        )
        Sys.sleep(5)
        expect_true(x$is_alive())
        x$kill()
    }
)
