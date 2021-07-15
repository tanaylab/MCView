library(golem)

test_that("app ui", {
    init_config(project = project_dir)
    load_all_data(cache_dir = project_cache_dir(project_dir))
    init_defs()
    ui <- app_ui()
    expect_shinytaglist(ui)
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
