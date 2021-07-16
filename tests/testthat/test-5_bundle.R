test_that("create_bundle works", {
    MCView::create_bundle(project = project_dir, path = bundle_dir, name = "PBMC")


    project_bundle <- fs::path(bundle_dir, "PBMC")
    expect_true(fs::dir_exists(project_bundle))
    expect_true(fs::dir_exists(fs::path(project_bundle, "project")))
    expect_true(fs::dir_exists(fs::path(project_bundle, "project", "cache")))
    expect_true(fs::dir_exists(fs::path(project_bundle, "project", "cache", "PBMC163k")))
    expect_true(fs::dir_exists(fs::path(project_bundle, "project", "config")))
    expect_true(fs::file_exists(fs::path(project_bundle, "app.R")))
    expect_true(fs::file_exists(fs::path(project_bundle, "project", "config", "config.yaml")))
    expect_true(fs::file_exists(fs::path(project_bundle, "project", "config", "about.Rmd")))
    expect_true(fs::file_exists(fs::path(project_bundle, "project", "config", "help.yaml")))

    test_dataset(fs::path(project_bundle, "project"), "PBMC163k")

    config_file <- project_config_file(project_dir)
    config <- yaml::read_yaml(config_file)
    config_bundle <- yaml::read_yaml(fs::path(project_bundle, "project", "config", "config.yaml"))

    expect_equal(config, config_bundle)
})

test_that(
    "bundle app launches",
    {
        cur_wd <- getwd()
        withr::defer(setwd(cur_wd))
        setwd(fs::path(bundle_dir, "PBMC"))

        skip_on_cran()
        skip_on_travis()
        skip_on_appveyor()
        x <- processx::process$new(
            "Rscript",
            "app.R"
        )
        Sys.sleep(5)
        expect_true(x$is_alive())
        x$kill()
    }
)
